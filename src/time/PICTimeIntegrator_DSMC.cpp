#include "PICTimeIntegrator_DSMC.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_DSMC::timeStep( const Real a_time,
                                       const Real a_dt,
                                       const int )
{  
  CH_TIME("System::advance_DSMC()");
  
  // explicit advance of particle positions 
  // + inertial forces
  // + velocity scatter
  // + special operators
   
  // advance particle positions from t_{n} to t_{n+1} using vp_{n}
  for (int s=0; s<m_particles.size(); s++) {
    auto this_picSpecies(m_particles[s]);
    this_picSpecies->advancePositions(a_dt);
    this_picSpecies->applyInertialForces(a_dt, false, true);
    this_picSpecies->applyBCs(false);
  }
  
  // scatter the particles: vp_{n+1} ==> vp'_{n+1}
  m_system->scatterParticles( a_dt );
  
  return;
}

#include "NamespaceFooter.H"

