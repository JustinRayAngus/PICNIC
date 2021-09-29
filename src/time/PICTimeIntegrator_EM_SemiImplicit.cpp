#include "PICTimeIntegrator_EM_SemiImplicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_SemiImplicit::preTimeStep(  const Real a_time,
                                                      const Real a_dt,
                                                      const int )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::preTimeStep()");
     
  Real cnormDt = a_dt*m_units->CvacNorm();
  Real cnormHalfDt = 0.5*cnormDt;

  // complete advance of B from t_{n} to t_{n+1/2}
  if (m_fields) {
    if (a_time == 0.0) { // initial advance of magnetic field by 1/2 time step
      m_fields->advanceMagneticField(a_time, cnormHalfDt);
    } else {
      m_fields->advanceMagneticField_2ndHalf(a_time);
    }
  }

  return;
}

void PICTimeIntegrator_EM_SemiImplicit::timeStep( const Real a_time,
                                                  const Real a_dt,
                                                  const int )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::advance()");

  // semi-implicit semi-energy-conservative scheme by Chen 2020
  // B is defined at half time steps while all other quantities
  // (E, xp, and vp) are defined at whole time steps
   
  Real cnormDt = a_dt*m_units->CvacNorm();
  Real cnormHalfDt = 0.5*cnormDt;
 
  int iter = 0;
  Real norm_efield = 0, norm_efield0 = 0;

  while( 1 ) {

    // Step 1: advance particle positions from t_{n} to t_{n+1/2} using vp_{n+1/2}
    for (int s=0; s<m_particles.size(); s++) {
      auto this_picSpecies(m_particles[s]);
      if(m_fields && iter>0) {
        if(this_picSpecies->forces()) { // reposition outcasts from previous iteration and update velocity
          this_picSpecies->repositionOutcastsAndApplyForces( *m_fields, cnormDt, true );
        }
      }
      if(this_picSpecies->motion()) {
        this_picSpecies->advancePositions(cnormHalfDt,true);
      }
    }
    
    // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance vp from t_{n} to t_{n+1/2}
    if (m_fields) {
  
      for (int s=0; s<m_particles.size(); s++) {
        auto this_picSpecies(m_particles[s]);
        if(this_picSpecies->forces()) {
          this_picSpecies->interpolateFieldsToParticles( *m_fields );
          this_picSpecies->advanceVelocities( cnormDt, true );
        }
      }
   
      // Step 3: compute current density at t_{n+1/2} and advance E from t_n to t_{n+1/2} 
      if(m_fields->advanceE()) m_system->setCurrentDensity();
      m_fields->saveElectricField();
      m_fields->advanceElectricField(a_time, cnormHalfDt);
      norm_efield = m_fields->diffElectricField();
      if (iter == 0) norm_efield0 = norm_efield;

      if (!procID()) {
        printf("  iter = %3d,", iter);
        printf(" norm_efield = %1.4e (abs.), %1.4e (rel.)\n",
               norm_efield, norm_efield/norm_efield0 );
      }
      if (norm_efield < m_atol) {
        if (!procID()) {
          printf("  exiting: satisfied absolute tolerance (%1.3e).\n", 
                 m_atol);
        }
        break;
      }
      if (norm_efield/norm_efield0 < m_rtol) {
        if (!procID()) {
          printf("  exiting: satisfied relative tolerance (%1.3e).\n", 
                 m_rtol);
        }
        break;
      }

      if ( iter >= m_iter_max ) {
        if (!procID()) {
          printf("  exiting: iterations exceed max iterations (%d).\n", 
                 m_iter_max);
        }
        break;
      }
      
    }
    else {
       break;
    }
    
    iter = iter + 1;

  } // end iteration loop

  // Step 4: 2nd half-advance of xp and vp from t_{n+1/2} to t_{n+1}
  for (int s=0; s<m_particles.size(); s++) {
     auto this_picSpecies(m_particles[s]);
     if(this_picSpecies->motion()) this_picSpecies->advancePositions_2ndHalf();
     if(this_picSpecies->forces()) this_picSpecies->advanceVelocities_2ndHalf();
  }

  // Step 5: 2nd half-advance of E from t_{n+1/2} to t_{n+1} and 
  //         1st half-advance of B from t_{n+1/2} to t_{n+1} using E_{n+1}
  if (m_fields) {
     m_fields->advanceElectricField_2ndHalf(a_time);
     m_fields->advanceMagneticField(a_time, cnormHalfDt);
  }

  // Step 6: scatter the particles: vp_{n+1} ==> vp'_{n+1}
  m_system->scatterParticles( a_dt );

}

#include "NamespaceFooter.H"

