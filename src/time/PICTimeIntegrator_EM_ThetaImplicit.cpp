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
  pp.query("theta_parameter", m_theta);
  pp.query("solver_type", m_nlsolver_type);
  pp.query("pc_update_freq", m_pc_update_freq);
  pp.query("pc_update_newton", m_pc_update_newton);

  CH_assert(m_theta>=0.5 && m_theta<=1.0);

  if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    m_nlsolver = new PicardSolver<ODEVector<System>, System>;
    dynamic_cast<PicardSolver<ODEVector<System>, System>*>
      (m_nlsolver)->numBlocks( m_system->numPicardBlocks() );
  
  } else if (m_nlsolver_type == _NLSOLVER_NEWTON_) {

    m_nlsolver = new NewtonSolver<ODEVector<System>, System>;
    m_func = new EMResidualFunction<ODEVector<System>, System>;
    m_func->define( m_U, m_system, m_theta, m_pc_update_newton );
    
  } else if (     (m_nlsolver_type == _NLSOLVER_PETSCSNES_) 
              ||  (m_nlsolver_type == "petsc") /* backward compatibility */) {

  if(m_nlsolver_type=="petsc") m_nlsolver_type = _NLSOLVER_PETSCSNES_;

#ifdef with_petsc
    m_nlsolver = new PetscSNESWrapper<ODEVector<System>, System>;
    m_func = new EMResidualFunction<ODEVector<System>, System>;
    m_func->define( m_U, m_system, m_theta, m_pc_update_newton );
#else
    MayDay::Error("PETSc SNES solver requires compilation with PETSc");
#endif
    
  } else {
    MayDay::Error("Invalid choice for PICTimeIntegrator_EM_ThetaImplicit::m_nlsolver_type");
  }
  m_nlsolver->define(m_U, m_system, m_func, m_theta);
  m_nlsolver->setOutputIndent("  ");
  
  Real rtol, atol;
  int maxits;
  m_nlsolver->getParams(rtol,atol,maxits);    
  if( (m_nlsolver_type == _NLSOLVER_PICARD_) && (maxits % 2 == 0) ) {
      cout << "WARNING!!! iter_max = " << maxits << endl;
      cout << "theta_implicit time advance with picard solver " << endl;
      cout << "may not work properly with even iter_max " << endl;
  } 
  
  // look for option to do energy-conserving semi-implicit
  // this requires using JFNK. Diff in maxits for petsc vs native newton
  // is due to how counting sequence goes for nonlinear iterations
  // (newton: 0,1,2... petsc: 0,0,1,2...)
  pp.query("ec_semi_implicit", m_ec_semi_implicit);
  if(m_ec_semi_implicit) {
    if(m_nlsolver_type==_NLSOLVER_PICARD_) m_ec_semi_implicit = false;
    if(m_nlsolver_type==_NLSOLVER_NEWTON_) {
       CH_assert(maxits==1);
       m_system->DefineECSemiImpFlag();
    }
    if(m_nlsolver_type==_NLSOLVER_PETSCSNES_) {
       CH_assert(maxits==2);
       m_system->DefineECSemiImpFlag();
    }
  }
 

  printParams();

  m_is_defined = true;
}

void PICTimeIntegrator_EM_ThetaImplicit::printParams() const
{
   if(procID()>0) return;
   cout << "================== Time Solver ==================" << endl;
   cout << "advance method = PIC_EM_THETA_IMPLICIT" << endl;
   cout << "theta_parameter = " << m_theta << endl;
   cout << "solver_type = " << m_nlsolver_type << endl;
   if(m_nlsolver_type != _NLSOLVER_PICARD_) { 
   cout << "ec_semi_implicit = " << m_ec_semi_implicit << endl;
   cout << "pc_update_freq = " << m_pc_update_freq << endl;
   cout << "pc_update_newton = " << (m_pc_update_newton?"true":"false") << endl;
   }
   m_nlsolver->printParams();
   cout << "=================================================" << endl;
   cout << endl;
}

void PICTimeIntegrator_EM_ThetaImplicit::initialize( const std::string&  a_restart_file_name )
{
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::initialize()");

  m_system->copySolutionToVec( m_U );
  m_Uold = m_U;

  if (m_func) { // this doesnt work on restart for some reason...
                // it goes here and calls updatePreCondMat(), but doest work....
                // It does work if I call it from timeStep()....
                // Note that updatePreCondMat() is also called at petsc define
                // Not sure why we need to also call it inside timeStep() call 
                // for it to work....
    //int dummy_step = 0;
    //dynamic_cast<EMJacobianFunction<ODEVector<System>, System>*>
    //    ( &(m_func->getJacobian()) )->updatePreCondMat( m_U );
  }
   
  if(!a_restart_file_name.empty()) {
  
    // read solver probes
#ifdef CH_USE_HDF5
    
    if(!procID()) cout << "Reading solver probe data from restart file..." << endl << endl;

    HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );
    HDF5HeaderData header;
   
    const std::string solverGroup = std::string("solver");
    handle.setGroup(solverGroup);
    header.readFromFile( handle );
   
    m_last_l_exit_status  = header.m_int["l_exit_status"];
    m_last_nl_exit_status = header.m_int["nl_exit_status"];
    m_last_l_total_iter = header.m_int["l_total_iter"];
    m_last_nl_total_iter = header.m_int["nl_total_iter"];
    m_last_l_last_iter = header.m_int["l_last_iter"];
    m_last_nl_iter = header.m_int["nl_iter"];
    m_last_nl_abs_res = header.m_real["nl_abs_res"];
    m_last_nl_rel_res = header.m_real["nl_rel_res"];
  
    handle.close();
#endif

  }

}

void PICTimeIntegrator_EM_ThetaImplicit::timeStep(  const Real a_time,
                                                    const Real a_dt,
                                                    const int  a_step )
{  
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::timeStep()");

  if (m_func) {
    m_func->curTime( a_time );
    m_func->curTimeStep( a_dt );
    if (a_step%m_pc_update_freq == 0 || m_restart_pc_update_flag) {
      dynamic_cast<EMJacobianFunction<ODEVector<System>, System>*>
          ( &(m_func->getJacobian()) )->updatePreCondMat( m_U );
      m_restart_pc_update_flag = false;
    }
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
    //this_picSpecies->mergeSubOrbitParticles();
    this_picSpecies->advanceVelocities_2ndHalf( a_dt );
    this_picSpecies->advancePositions_2ndHalf();
    this_picSpecies->mergeSubOrbitParticles();
    this_picSpecies->applyInertialForces(a_dt,true,true);
    this_picSpecies->applyBCs(false);
  }

  m_U = (1.0/m_theta)*m_U + ((m_theta-1.0)/m_theta)*m_Uold;
  m_system->updatePhysicalState( m_U, a_time );
  m_system->scatterParticles( a_dt );

}

#include "NamespaceFooter.H"

