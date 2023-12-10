#include "PICTimeIntegrator_EM_Explicit.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_Explicit::initialize( const std::string&  a_restart_file_name )
{
  m_fields->copyEToVec( m_E );
  m_fields->copyBToVec( m_B );
  m_fields->copyEoldToVec( m_Eold );
  m_fields->copyBoldToVec( m_Bold );
}

int PICTimeIntegrator_EM_Explicit::prepForCheckpoint() const
{
  m_fields->copyEoldFromVec( m_Eold );
  m_fields->copyBoldFromVec( m_Bold );
  return 1;
}

void PICTimeIntegrator_EM_Explicit::preTimeStep(  const Real a_time,
                                                  const Real a_dt,
                                                  const int )
{  
  CH_TIME("PICTimeIntegrator_EM_Explicit::preTimeStep()");
     
  const Real halfDt = 0.5*a_dt;
  const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
  const int num_species = pic_species_ptr_vect.size();

  // complete advance of xp from t_{n} to t_{n+1/2}
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    if (a_time==0.0) { // initial advance of particles positions by 1/2 time step
      species->createInflowParticles( a_time, halfDt );
      species->advancePositionsExplicit( halfDt, false ); // false is correct here
    } 
    else {
      species->advancePositions_2ndHalf();
    }
    species->applyBCs( false );
    species->removeOutflowParticles(); // needed here for method used to achieve charge
                                       // conservation with outflow particles. Outflow parts
                                       // created here are those that have xpbar from last 
                                       // timeStep inside the physical domain and thus they 
                                       // are retained in the main particle list and their current
                                       // is deposited correctly. If we do not remove them here,
                                       // then their current will contribute during the 
                                       // next timeStep() via the deposit from the outflow list,
                                       // which should only contain particles with xpbar outside
                                       // the physical domain created during the half advance.
  }
  
  // update old particle values and create inflow particles  
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    species->updateOldParticlePositions();
    species->updateOldParticleVelocities();
    species->createInflowParticles( a_time, a_dt );
  }
  
  // complete advance of Eg from t_{n} to t_{n+1/2}
  if (a_time == 0.0) {
    m_fields->computeRHSElectricField( a_time, halfDt );
    m_fields->copyERHSToVec( m_FE );
    m_E = m_Eold + m_FE;
  } else {
    m_E = 2.0*m_E - m_Eold;
  }
  m_fields->copyEFromVec( m_E );
  if( m_fields->usePoisson() && a_time == 0.0 ) {
    m_pic_species->setChargeDensityOnNodes( m_fields->useFiltering() );
    const LevelData<NodeFArrayBox>& pic_rho = m_pic_species->getChargeDensityOnNodes();
    m_fields->enforceGaussLaw( pic_rho, a_time );
  }
  m_fields->applyBCs_electricField( a_time );
  m_fields->copyEToVec( m_E );

  // update old field values 
  m_Eold = m_E;
  m_Bold = m_B;
  
  return;
}

void PICTimeIntegrator_EM_Explicit::timeStep( const Real a_time,
                                              const Real a_dt,
                                              const int )
{  
  CH_TIME("PICTimeIntegrator_EM_Explicit::timeStep()");
  
  // explicit leap-frog time advance of particles and fields
  // B and vp are defined a whole time steps while E and xp are at half
   
  const Real halfDt = 0.5*a_dt;
  const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
  const int num_species = pic_species_ptr_vect.size();
 
  if (m_fields) {
 
    // Step 1: advance B from t_{n} to t_{n+1/2} using E_{n+1/2}
    m_fields->computeRHSMagneticField( a_time, halfDt );
    m_fields->copyBRHSToVec( m_FB );
    m_B = m_Bold + m_FB;
    m_fields->updatePhysicalState( m_B, a_time, b_only );
    m_fields->updateBoundaryProbes( a_dt );

    // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2}
    //         advance vp from t_{n} to t_{n+1}
    for (int sp=0; sp<num_species; sp++) {
      auto species(pic_species_ptr_vect[sp]);
      species->interpolateFieldsToParticles( *m_fields );
      species->addExternalFieldsToParticles( *m_fields ); 
      species->advanceVelocities( a_dt, false );
    }

  }

  // advance positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
  for (int sp=0; sp<num_species; sp++) {
     auto species(pic_species_ptr_vect[sp]);
    species->advancePositionsExplicit( a_dt, true );
    species->applyInertialForces( a_dt, false, true, true );
    species->applyBCs( true );
  }

  if (m_fields) {
 
    // complete advance of B from t_{n+1/2} to t_{n+1} and compute the curl
    m_B = 2.0*m_B - m_Bold;
    m_fields->updatePhysicalState( m_B, a_time, b_only );

    // compute current density at t_{n+1} and 
    // advance E from t_{n+1/2} to t_{n+1} using B_{n+1} and J_{n+1}
    if(m_fields->advanceE()) {
      m_pic_species->setCurrentDensity( true );
      if(m_fields->useFiltering()) { m_pic_species->filterJ( *m_fields, a_time ); }
      const LevelData<EdgeDataBox>& pic_J = m_pic_species->getCurrentDensity();
      const LevelData<NodeFArrayBox>& pic_Jv = m_pic_species->getVirtualCurrentDensity();
      m_fields->setCurrentDensity( pic_J, pic_Jv );
    }
    m_fields->computeRHSElectricField( a_time, halfDt );
    m_fields->copyERHSToVec( m_FE );
    m_E = m_Eold + m_FE;

    m_fields->copyEFromVec( m_E );
    if(m_fields->usePoisson()) {
      m_pic_species->setChargeDensityOnNodes( m_fields->useFiltering() );
      const LevelData<NodeFArrayBox>& pic_rho = m_pic_species->getChargeDensityOnNodes();
      m_fields->enforceGaussLaw( pic_rho, a_time );
    }
    m_fields->applyBCs_electricField( a_time );
    m_fields->copyEToVec( m_E );

  }

  return;
}

#include "NamespaceFooter.H"

