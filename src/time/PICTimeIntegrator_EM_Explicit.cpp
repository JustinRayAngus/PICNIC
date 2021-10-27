#include "PICTimeIntegrator_EM_Explicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_Explicit::initialize()
{
  m_system->copyEToVec( m_E );
  m_system->copyBToVec( m_B );
  m_Eold = m_E;
  m_Bold = m_B;
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
      this_picSpecies->advancePositions(halfDt,false); // false is correct here
    } else {
      this_picSpecies->advancePositions_2ndHalf();
    }
  }

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
    
  if (m_fields) {
 
    // Step 1: advance B from t_{n} to t_{n+1/2} using E_{n+1/2}
    m_system->computeRHS( m_FB, m_E, a_time, halfDt, b_only );
    m_B = m_Bold + m_FB;
    m_system->updatePhysicalState( m_B, a_time, b_only );

    // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance vp from t_{n} to t_{n+1}
    for (int s=0; s<m_particles.size(); s++) {
      auto this_species(m_particles[s]);
      this_species->interpolateFieldsToParticles( *m_fields );
      this_species->advanceVelocities( a_dt, false );
    }

  }
  
  // scatter the particles: vp_{n+1} ==> vp'_{n+1}
  m_system->scatterParticles( a_dt );

  // advance particle positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
  for (int s=0; s<m_particles.size(); s++) {
    auto this_species(m_particles[s]);
    this_species->advancePositions(halfDt,true);
  }
  
  if (m_fields) {
 
    // complete advance of B from t_{n+1/2} to t_{n+1} and compute the curl
    m_B = 2.0*m_B - m_Bold;
    m_system->updatePhysicalState( m_B, a_time, b_only );

    // compute current density at t_{n+1} and advance E from t_{n+1/2} to t_{n+1} 
    // using B_{n+1} and J_{n+1}
    if(m_fields->usePoisson()) m_system->setChargeDensity();
    if(m_fields->advanceE()) m_system->setCurrentDensity();
    m_system->computeRHS( m_FE, m_B, a_time, halfDt, e_only );
    m_E = m_Eold + m_FE;
    m_system->updatePhysicalState( m_E, a_time, e_only );
  }
  
  return;
}

#include "NamespaceFooter.H"

