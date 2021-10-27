#include "PICTimeIntegrator_EM_ThetaImplicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_ThetaImplicit::define(  System* const             a_sys,
                                                  const PicSpeciesPtrVect&  a_particles,
                                                  ElectroMagneticFields* const )
{
  CH_assert(!isDefined());

  m_system = a_sys;
  m_particles = a_particles;

  m_U.define(*m_system);
  m_Uold.define(m_U);
  m_U.printLoadBalanceInfo();

  ParmParse pp("pic_em_theta_implicit");
  pp.query("solver_type", m_nlsolver_type);

  if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    m_nlsolver = new PicardSolver<ODEVector<System>, System>;
    dynamic_cast<PicardSolver<ODEVector<System>, System>*>
      (m_nlsolver)->numBlocks( m_system->numPicardBlocks() );

  } else if (m_nlsolver_type == _NLSOLVER_NEWTON_) {

    m_nlsolver = new NewtonSolver<ODEVector<System>, System>;
    m_func = new EMResidualFunction<ODEVector<System>, System>;
    m_func->define( m_U, m_system, m_theta );

  } else {
    MayDay::Error("Invalid choice for PICTimeIntegrator_EM_ThetaImplicit::m_nlsolver_type");
  }
  m_nlsolver->define(m_U, m_system, m_func, m_theta);

  m_is_defined = true;
}

void PICTimeIntegrator_EM_ThetaImplicit::initialize()
{
  m_system->copySolutionToVec( m_U );
  m_Uold = m_U;
}

void PICTimeIntegrator_EM_ThetaImplicit::timeStep(  const Real a_time,
                                                    const Real a_dt,
                                                    const int )
{  
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::timeStep()");

  if (m_func) {
    m_func->curTime( a_time );
    m_func->curTimeStep( a_dt );
  }
  m_U = m_Uold;
  m_nlsolver->solve( m_U, m_Uold, a_time, a_dt );

}

void PICTimeIntegrator_EM_ThetaImplicit::postTimeStep(  const Real a_time,
                                                        const Real a_dt )
{  
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::postTimeStep()");

  for (int s=0; s<m_particles.size(); s++) {
    auto this_picSpecies(m_particles[s]);
    this_picSpecies->advancePositions_2ndHalf();
    this_picSpecies->advanceVelocities_2ndHalf();
  }

  m_U = (1.0/m_theta)*m_U + ((m_theta-1.0)/m_theta)*m_Uold;
  m_system->updatePhysicalState( m_U, a_time );
  m_system->scatterParticles( a_dt );
}

#include "NamespaceFooter.H"

