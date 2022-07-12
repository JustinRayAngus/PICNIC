#include "PICTimeIntegrator_EM_SemiImplicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_SemiImplicit::define( System* const             a_sys,
                                                const PicSpeciesPtrVect&  a_particles,
                                                ElectroMagneticFields* const )
{
  CH_assert(!isDefined());

  m_system = a_sys;
  m_particles = a_particles;

  m_E.define(*m_system);
  m_Eold.define(m_E);
  m_E.printLoadBalanceInfo();

  m_B.define(m_system->getVectorSize(b_only));
  m_Bold.define(m_B);
  m_FB.define(m_B);

  m_nlsolver = new PicardSolver<ODEVector<System>, System>;
  m_nlsolver->define(m_E, m_system, nullptr, 0.5);

  m_is_defined = true;
}

void PICTimeIntegrator_EM_SemiImplicit::initialize()
{
  m_system->copyEToVec( m_E );
  m_system->copyBToVec( m_B );
  m_system->copyEoldToVec( m_Eold );
  m_system->copyBoldToVec( m_Bold );
}

int PICTimeIntegrator_EM_SemiImplicit::prepForCheckpoint() const
{
  m_system->copyEoldFromVec( m_Eold );
  m_system->copyBoldFromVec( m_Bold );
  return 1;
}

void PICTimeIntegrator_EM_SemiImplicit::preTimeStep(  const Real a_time,
                                                      const Real a_dt,
                                                      const int )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::preTimeStep()");
     
  if (a_time == 0.0) {
    m_system->computeRHS( m_FB, m_E, 0.0, a_dt, b_only );
    m_B = m_Bold + m_FB;
  } else {
    m_B = 2.0*m_B - m_Bold;
  }
  m_system->updatePhysicalState( m_B, a_time, b_only );

  m_Eold = m_E;
  m_Bold = m_B;

  return;
}

void PICTimeIntegrator_EM_SemiImplicit::timeStep( const Real a_time,
                                                  const Real a_dt,
                                                  const int )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::timeStep()");
  m_nlsolver->solve( m_E, m_Eold, a_time, a_dt );
}

void PICTimeIntegrator_EM_SemiImplicit::postTimeStep( const Real a_time,
                                                      const Real a_dt )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::postTimeStep()");
  
  const Real halfDt = 0.5*a_dt;

  for (int s=0; s<m_particles.size(); s++) {
     auto this_picSpecies(m_particles[s]);
     this_picSpecies->advanceVelocities_2ndHalf();
     this_picSpecies->advancePositions_2ndHalf();
     this_picSpecies->applyInertialForces(a_dt,true,true);
     this_picSpecies->applyBCs(false);
  }

  m_E = 2.0*m_E - m_Eold;
  m_system->updatePhysicalState( m_E, a_time, e_only );

  m_system->computeRHS( m_FB, m_E, a_time, halfDt, b_only );
  m_B = m_Bold + m_FB;
  m_system->updatePhysicalState( m_B, a_time, b_only );

  m_system->scatterParticles( a_dt );
}

#include "NamespaceFooter.H"

