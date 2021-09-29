#include "PICTimeIntegrator_EM_Explicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_Explicit::preTimeStep(  const Real a_time,
                                                  const Real a_dt,
                                                  const int )
{  
  CH_TIME("PICTimeIntegrator_EM_Explicit::preTimeStep()");
     
  Real cnormDt = a_dt*m_units->CvacNorm();
  Real cnormHalfDt = 0.5*cnormDt;

  if (m_fields) {
    // initial advance of electric field by 1/2 time step
    if (a_time==0.0) { 
      if(m_fields->usePoisson()) m_system->setChargeDensity();
      m_fields->advanceElectricField(a_time, cnormHalfDt);
    } else {
      m_fields->advanceElectricField_2ndHalf(a_time);
    }
  }

  // complete advance of xp from t_{n} to t_{n+1/2}
  for (int s=0; s<m_particles.size(); s++) {
    auto this_picSpecies(m_particles[s]);
    if(this_picSpecies->motion()) {
      // initial advance of particle positions by 1/2 time step
      if (a_time==0.0) { 
        this_picSpecies->createInflowParticles( a_time, cnormHalfDt );
        this_picSpecies->advancePositions(cnormHalfDt,false); // false is correct here
      } else {
        this_picSpecies->advancePositions_2ndHalf();
      }
    }
  }
  
  return;
}

void PICTimeIntegrator_EM_Explicit::timeStep( const Real a_time,
                                              const Real a_dt,
                                              const int )
{  
  CH_TIME("PICTimeIntegrator_EM_Explicit::advance()");
  
  // explicit leap-frog time advance of particles and fields
  // B and vp are defined a whole time steps while E and xp are at half
   
  Real cnormDt = a_dt*m_units->CvacNorm();
  Real cnormHalfDt = 0.5*cnormDt;
    
  Real cur_time = 0.0; // dummy. Add to passed arguements latter

  if (m_fields) {
 
    // Step 1: advance B from t_{n} to t_{n+1/2} using E_{n+1/2}
    m_fields->advanceMagneticField(a_time, cnormHalfDt);

    // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance vp from t_{n} to t_{n+1}
    for (int s=0; s<m_particles.size(); s++) {
      auto this_species(m_particles[s]);
      if(this_species->forces()) {
        this_species->interpolateFieldsToParticles( *m_fields );
        this_species->advanceVelocities( cnormDt, false );
      }
    }

  }
  
  // scatter the particles: vp_{n+1} ==> vp'_{n+1}
  m_system->scatterParticles( a_dt );

  // advance particle positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
  for (int s=0; s<m_particles.size(); s++) {
    auto this_species(m_particles[s]);
    if(this_species->motion()) this_species->advancePositions(cnormHalfDt,true);
  }
  
  if (m_fields) {
 
    // complete advance of B from t_{n+1/2} to t_{n+1} and compute the curl
    m_fields->advanceMagneticField_2ndHalf(a_time);

    // compute current density at t_{n+1} and advance E from t_{n+1/2} to t_{n+1} 
    // using B_{n+1} and J_{n+1}
    if(m_fields->advanceE()) m_system->setCurrentDensity();
    if(m_fields->usePoisson()) m_system->setChargeDensity();
    m_fields->advanceElectricField(a_time, cnormHalfDt);

  }
  
  return;
}

#include "NamespaceFooter.H"

