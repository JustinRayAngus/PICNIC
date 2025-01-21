#include "PICTimeIntegrator_DSMC.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_DSMC::preTimeStep( const Real a_time,
                                          const Real a_dt,
                                          const int )
{  
  CH_TIME("PICTimeIntegrator_DSMC::preTimeStep()");
  
  // update old particle values and create inflow particles  
  const PicChargedSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getChargedPtrVect();
  const int num_species = pic_species_ptr_vect.size();
  for (int sp=0; sp<num_species; sp++) {
      auto species = pic_species_ptr_vect[sp];
      species->updateOldParticlePositions();
      species->updateOldParticleVelocities();
      species->createInflowParticles( a_time, a_dt );
  }
     
  return;
}

int PICTimeIntegrator_DSMC::timeStep( const Real  a_time,
                                      const Real  a_dt,
                                      const int )
{  
  CH_TIME("PICTimeIntegrator_DSMC::timeStep()");
  
  // explicit advance of particle positions 
  // + inertial forces
   
  // advance particle positions from t_{n} to t_{n+1} using vp_{n}
  const PicChargedSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getChargedPtrVect();
  const int num_species = pic_species_ptr_vect.size();
  for (int sp=0; sp<num_species; sp++) {
      auto species = pic_species_ptr_vect[sp];
      species->advancePositionsExplicit(a_dt);
      species->applyInertialForces(a_dt, false, true);
      species->applyBCs( false, a_time );
  }
  
  return 0;
}

#include "NamespaceFooter.H"

