#include "PICTimeIntegrator_EM_Explicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_Explicit::initialize( const std::string&  a_restart_file_name )
{
  m_system->copyEToVec( m_E );
  m_system->copyBToVec( m_B );
  m_system->copyEoldToVec( m_Eold );
  m_system->copyBoldToVec( m_Bold );
}

int PICTimeIntegrator_EM_Explicit::prepForCheckpoint() const
{
  m_system->copyEoldFromVec( m_Eold );
  m_system->copyBoldFromVec( m_Bold );
  return 1;
}

void PICTimeIntegrator_EM_Explicit::preTimeStep(  const Real a_time,
                                                  const Real a_dt,
                                                  const int )
{  
  CH_TIME("PICTimeIntegrator_EM_Explicit::preTimeStep()");
     
  const Real halfDt = 0.5*a_dt;

  if (a_time == 0.0) {
    if(m_fields->usePoisson()) m_system->setChargeDensity();
    m_system->computeRHS( m_FE, m_B, a_time, halfDt, e_only );
    m_E = m_Eold + m_FE;
  } else {
    m_E = 2.0*m_E - m_Eold;
  }
  m_system->updatePhysicalState( m_E, a_time, e_only );
  
  // complete advance of xp from t_{n} to t_{n+1/2}
  for (int s=0; s<m_particles.size(); s++) {
    auto this_picSpecies(m_particles[s]);
    if (a_time==0.0) { // initial advance of particles positions by 1/2 time step
      this_picSpecies->createInflowParticles( a_time, halfDt );
      this_picSpecies->advancePositions( halfDt, false ); // false is correct here
    } 
    else {
      this_picSpecies->advancePositions_2ndHalf();
    }
    this_picSpecies->applyBCs( false );
    this_picSpecies->removeOutflowParticles(); // needed here for method used to achieve charge
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

  m_Eold = m_E;
  m_Bold = m_B;
  
  // update old particle values and create inflow particles  
  for (int s=0; s<m_particles.size(); s++) {
    auto this_picSpecies(m_particles[s]);
    this_picSpecies->updateOldParticlePositions();
    this_picSpecies->updateOldParticleVelocities();
    this_picSpecies->createInflowParticles( a_time, a_dt );
    //this_picSpecies->injectInflowParticles();
  }

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
 
  if (m_fields) {
 
    // Step 1: advance B from t_{n} to t_{n+1/2} using E_{n+1/2}
    m_system->computeRHS( m_FB, m_E, a_time, halfDt, b_only );
    m_B = m_Bold + m_FB;
    m_system->updatePhysicalState( m_B, a_time, b_only );

    // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2}
    //         advance vp from t_{n} to t_{n+1}
    for (int s=0; s<m_particles.size(); s++) {
      auto this_species(m_particles[s]);
      this_species->interpolateFieldsToParticles( *m_fields );
      this_species->addExternalFieldsToParticles( *m_fields ); 
      this_species->advanceVelocities( a_dt, false );
    }

  }

  if (m_old_ordering) {

    // scatter the particles: vp_{n+1} ==> vp'_{n+1}
    m_system->scatterParticles( a_dt );

    // advance positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
    for (int s=0; s<m_particles.size(); s++) {
      auto this_species(m_particles[s]);
      this_species->advancePositions( a_dt, true );
      // not compatible with Boris method for inertia
      this_species->applyBCs( true );
    }

  }
  else {

    // advance positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
    for (int s=0; s<m_particles.size(); s++) {
      auto this_species(m_particles[s]);
      this_species->advancePositions( a_dt, true );
      this_species->applyInertialForces( a_dt, false, true, true );
      this_species->applyBCs( true );
    }

    // scatter the particles: vp_{n+1} ==> vp'_{n+1}
    // doing here works well for numerical energy conservation test,
    // but fails badly for dynamicPinch and steadyStateShock sims
    // need to apply scattering after computing J
    //m_system->scatterParticles( a_dt );

  }

  if (m_fields) {
 
    // complete advance of B from t_{n+1/2} to t_{n+1} and compute the curl
    m_B = 2.0*m_B - m_Bold;
    m_system->updatePhysicalState( m_B, a_time, b_only );

    // compute current density at t_{n+1} and 
    // advance E from t_{n+1/2} to t_{n+1} using B_{n+1} and J_{n+1}
    if(m_fields->usePoisson()) m_system->setChargeDensity();
    if(m_fields->advanceE()) m_system->setCurrentDensity( true );
    m_system->computeRHS( m_FE, m_B, a_time, halfDt, e_only );
    m_E = m_Eold + m_FE;
    m_system->updatePhysicalState( m_E, a_time, e_only );

  }

  // scatter the particles: vp_{n+1} ==> vp'_{n+1}
  if(!m_old_ordering) m_system->scatterParticles( a_dt );
  
  return;
}

#include "NamespaceFooter.H"

