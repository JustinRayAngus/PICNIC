#include "PICTimeIntegrator_EM_ThetaImplicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_ThetaImplicit::timeStep( const Real a_time,
                                              const Real a_dt,
                                              const int )
{  
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::advance()");

  // fully-implicit energy-conservative scheme by Markidis and
  // Lapenta 2011. All quantities (E, B, xp, and vp) are defined 
  // at whole time steps
   
  Real cnormDt = a_dt*m_units->CvacNorm();
  Real cnormHalfDt = 0.5*cnormDt;
  Real cnormThetaDt = m_theta*cnormDt;
 
  int iter = 0;
  Real  norm_efield = 0, norm_efield0 = 0, 
        norm_bfield = 0, norm_bfield0 = 0;

  while(1) {

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
    
    if (m_fields) {
   
      // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2} and advance 
      // vp from t_{n} to t_{n+1/2}
      for (int s=0; s<m_particles.size(); s++) {
        auto this_picSpecies(m_particles[s]);
        if(this_picSpecies->forces()) {
          this_picSpecies->interpolateFieldsToParticles( *m_fields );
          this_picSpecies->advanceVelocities( cnormDt, true );
        }
      }
   
      // Step 3: compute current density at t_{n+1/2} and advance B and E 
      // from t_n to t_{n+1/2} 
      if (m_fields->advance()) {
        
        // advance B first
        m_fields->saveMagneticField();
        m_fields->advanceMagneticField(a_time, cnormThetaDt);
        norm_bfield = m_fields->diffMagneticField();

        // then advance E
        if(m_fields->advanceE()) m_system->setCurrentDensity();
        m_fields->saveElectricField();
        m_fields->advanceElectricField(a_time, cnormThetaDt);
        norm_efield = m_fields->diffElectricField();

        if (iter == 0) norm_bfield0 = norm_bfield;
        if (iter == 0) norm_efield0 = norm_efield;

        if (!procID()) {
          printf("  iter = %3d,", iter);
          printf(" norm_efield = %1.4e (abs.), %1.4e (rel.),",
                 norm_efield, norm_efield/norm_efield0 );
          printf(" norm_bfield = %1.4e (abs.), %1.4e (rel.)\n",
                 norm_bfield, norm_bfield/norm_bfield0 );
        }
        if ((norm_efield < m_atol) && (norm_bfield < m_atol)) {
          if (!procID()) {
            printf("  exiting: satisfied absolute tolerance (%1.3e).\n", 
                   m_atol);
          }
          break;
        }
        if (     (norm_efield/norm_efield0 < m_rtol)
             &&  (norm_bfield/norm_bfield0 < m_rtol) ) {
          if (!procID()) {
            printf("  exiting: satisfied relative tolerance (%1.3e).\n", 
                   m_rtol);
          }
          break;
        }

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

  // Step 4: 2nd half-advance of all quantities (E, B, xp, and vp)
  // from t_{n+1/2} to t_{n+1}
  for (int s=0; s<m_particles.size(); s++) {
    auto this_picSpecies(m_particles[s]);
    if(this_picSpecies->motion()) this_picSpecies->advancePositions_2ndHalf();
    if(this_picSpecies->forces()) this_picSpecies->advanceVelocities_2ndHalf();
  }

  if (m_fields) {
    m_fields->advanceElectricField_2ndHalf(a_time, m_theta);
    m_fields->advanceMagneticField_2ndHalf(a_time, m_theta);
  }

  // Step 6: scatter the particles: vp_{n+1} ==> vp'_{n+1}
  m_system->scatterParticles( a_dt );

}

#include "NamespaceFooter.H"

