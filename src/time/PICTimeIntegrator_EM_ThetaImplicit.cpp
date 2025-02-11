#include "PICTimeIntegrator_EM_ThetaImplicit.H"
#include "CH_HDF5.H"
#include <sys/time.h>

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_ThetaImplicit::define( System* const               a_system,
                                                 PicSpeciesInterface* const  a_pic_species,
                                                 EMFields* const             a_fields )
{
  CH_assert(!isDefined());

  m_pic_species = a_pic_species;
  m_fields = a_fields;

  // parse implicit solver parameters
  ParmParse pp_old("pic_em_theta_implicit");
  if (pp_old.contains("solver_type")) {
    if (!procID()) {
      std::cout << "EXIT FAILURE!!!" << std::endl;
      std::cout << "Prefix pic_em_theta_implicit no longer accepted." << std::endl;
      std::cout << "Use implicit_solver prefix instead." << std::endl << std::endl;
    }
    exit(EXIT_FAILURE);
  }
  pp_old.query("curl_curl_form", m_curl2);
  pp_old.query("ec_semi_implicit", m_ec_semi_implicit);
  pp_old.query("theta_parameter", m_theta);
  CH_assert(m_theta>=0.5 && m_theta<=1.0);

  ParmParse pp("implicit_solver");
  parseCommonImplicitSolverParams( pp );

  m_vec_type = (m_curl2 ? curl2 : e_and_b);

  m_fields->computePrecondMatrixNNZ(m_vec_type);
  m_fields->defineVectorsAndDOFs(m_vec_type);

  m_U.define(m_fields->getVectorSize(m_vec_type));
  m_Uold.define(m_U);
  m_U.printLoadBalanceInfo();

  //backward compatibility
  if (m_nlsolver_type=="petsc") { m_nlsolver_type = _NLSOLVER_PETSCSNES_; }

  if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    if (m_curl2) MayDay::Error("Picard solver not yet implemeted for curl-curl form");
    m_nlsolver = new PicardSolver<ODEVector<EMFields>, PICTimeIntegrator>;
    dynamic_cast<PicardSolver<ODEVector<EMFields>, PICTimeIntegrator>*>
      (m_nlsolver)->numBlocks( m_fields->numPicardBlocks() );

  } else if (m_nlsolver_type == _NLSOLVER_NEWTON_) {

    m_nlsolver = new NewtonSolver<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func = new EMResidualFunction<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func->define( m_U, this, m_theta, m_pc_update_newton );

  } else if (m_nlsolver_type == _NLSOLVER_PETSCSNES_) {

#ifdef with_petsc
    m_nlsolver = new PetscSNESWrapper<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func = new EMResidualFunction<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func->define( m_U, this, m_theta, m_pc_update_newton );
#else
    MayDay::Error("PETSc SNES solver requires compilation with PETSc");
#endif

  } else {
    MayDay::Error("Invalid choice for PICTimeIntegrator_EM_ThetaImplicit::m_nlsolver_type");
  }
  m_nlsolver->define(m_U, this, m_func, m_theta);
  m_nlsolver->setOutputIndent("  ");

  Real rtol, atol;
  int maxits;
  m_nlsolver->getParams(rtol,atol,maxits);
  if( (m_nlsolver_type == _NLSOLVER_PICARD_) && (maxits % 2 == 0) ) {
      cout << "WARNING!!! iter_max = " << maxits << endl;
      cout << "theta_implicit time advance with picard solver " << endl;
      cout << "may not work properly with even iter_max " << endl;
  }

  // ec_semi_implicit is for predictor-corrector method that achieves
  // exact energy conservation, but only approximate charge conservation
  if(m_ec_semi_implicit) {
    if(m_nlsolver_type==_NLSOLVER_PETSCSNES_ || m_nlsolver_type==_NLSOLVER_PICARD_) {
      if(!procID()) {
        cout << "EXIT FAILURE!!!" << endl;
        cout << "m_nlsolver_type = petsc and picard are not a valid option" << endl;
        cout << "when using ec_semi_implicit = true. For petsc, this is " << endl;
	    cout << "owed to the ordering of calls in the SNES Newton loop." << endl;
        cout << "Use m_nlsolver_type = newton." << endl;
      }
      exit(EXIT_FAILURE);
    }
    if (m_curl2) {
      MayDay::Error("Curl-curl form not yet supported for ec_semi_implicit=true");
    }
  }

  printParams();
  m_is_defined = true;
}

void PICTimeIntegrator_EM_ThetaImplicit::printParams() const
{
   if (procID()>0) { return; }
   cout << "================== Time Solver ==================" << endl;
   cout << "advance method = PIC_EM_THETA_IMPLICIT" << endl;
   cout << "theta_parameter = " << m_theta << endl;
   cout << "curl-curl form = " << (m_curl2 ? "true" : "false") << "\n";
   cout << "solver_type = " << m_nlsolver_type << endl;
   cout << "require_convergnce = " << (m_require_convergence?"true":"false") << endl;
   if (m_nlsolver_type != _NLSOLVER_PICARD_) {
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

  if (m_curl2) {
    m_fields->copyEToVec(m_U);
  } else {
    int offset(0);
    m_fields->copyBToVec( m_U, offset );
    m_fields->copyEToVec( m_U, offset );
  }
  m_Uold = m_U;

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

void PICTimeIntegrator_EM_ThetaImplicit::resetVariables( const Real a_time,
                                                         const Real a_dt,
                                                         const int  a_step )
{
    CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::resetVariables()");

    // reset the field values to those at the start of the time step
    m_U = m_Uold;
    int offset = 0;
    m_fields->copyBFromVec( m_U, offset );
    m_fields->copyEFromVec( m_U, offset );

    // reset the nonlinear iteration to zero
    m_nl_iter = 0;

    // reset particle values to those at beginning of time step
    const PicChargedSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getChargedPtrVect();
    const int num_species = pic_species_ptr_vect.size();
    for (int sp=0; sp<num_species; sp++) {
        auto species(pic_species_ptr_vect[sp]);
        species->resetParticles();
    }

}

void PICTimeIntegrator_EM_ThetaImplicit::preTimeStep( const Real a_time,
                                                      const Real a_dt,
                                                      const int  a_step )
{
  m_Uold = m_U;
  m_fields->setEoldFromE(); // needed for absorbing BC
  m_fields->setBoldFromB();
  m_nl_iter = 0;

  // update solver verbosity
  bool this_verbosity = false;
  int cout_step_interval = 1;
  if(a_step<10) cout_step_interval = 1;
  else if(a_step<100) cout_step_interval = 10;
  else cout_step_interval = 100;

  if( ((a_step+1) % cout_step_interval)==0 ) {
     this_verbosity = true;
  }
  m_nlsolver->verbose(this_verbosity);

  // update old particle values and create inflow particles
  const PicChargedSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getChargedPtrVect();
  const int num_species = pic_species_ptr_vect.size();
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    species->updateOldParticlePositions();
    species->updateOldParticleVelocities();
    species->createInflowParticles( a_time, a_dt );
    species->injectInflowParticles();
  }

}


int PICTimeIntegrator_EM_ThetaImplicit::timeStep( const Real  a_old_time,
                                                  const Real  a_dt,
                                                  const int   a_step )
{
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::timeStep()");

  const Real half_time = a_old_time + a_dt/2.0;

  struct timeval pcassembly_start, pcassembly_end;
  struct timeval solve_start, solve_end;

  gettimeofday(&pcassembly_start,NULL);
  if (m_func) {
    m_func->curTime( half_time );
    m_func->curTimeStep( a_dt );
    if (a_step%m_pc_update_freq == 0 || m_restart_pc_update_flag) {
      dynamic_cast<EMJacobianFunction<ODEVector<EMFields>, PICTimeIntegrator>*>
          ( &(m_func->getJacobian()) )->updatePreCondMat( m_U );
      m_restart_pc_update_flag = false;
    }
  }
  gettimeofday(&pcassembly_end,NULL);
  {
    long long walltime;
    walltime = (  (pcassembly_end.tv_sec * 1000000   + pcassembly_end.tv_usec  )
                - (pcassembly_start.tv_sec * 1000000 + pcassembly_start.tv_usec));
    double pcassembly_runtime = (double) walltime / 1000000.0;
#ifdef CH_MPI
    MPI_Allreduce(&pcassembly_runtime,&m_pcassembly_walltime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
    m_pcassembly_walltime = pcassembly_runtime;
#endif
  }

  m_U = m_Uold;
  gettimeofday(&solve_start,NULL);
  // advance state (Eg, Bg, xp, vp) from old_time to half_time
  m_nlsolver->solve( m_U, m_Uold, half_time, a_dt );

  // check for solver failure
  int exit_status = m_nlsolver->getExitStatus();
  if (exit_status<0 && m_require_convergence) {
     return -1;
  }

  if (m_curl2) {
    m_fields->copyEFromVec( m_U );
    m_fields->applyBCs_electricField( half_time );
    m_fields->updateB(m_theta*a_dt, half_time);
  }
  if(m_ec_semi_implicit) {
    CH_assert(!m_curl2);
    int offset(0);
    m_fields->copyBFromVec( m_U, offset );
    m_fields->copyEFromVec( m_U, offset );
    m_fields->applyBCs_magneticField( half_time );
    m_fields->applyBCs_electricField( half_time );
    m_pic_species->postNewtonUpdate( *m_fields, half_time, a_dt);
  }
  gettimeofday(&solve_end,NULL);

  {
    long long walltime;
    walltime = (  (solve_end.tv_sec * 1000000   + solve_end.tv_usec  )
                - (solve_start.tv_sec * 1000000 + solve_start.tv_usec));
    double solver_runtime = (double) walltime / 1000000.0;
#ifdef CH_MPI
    MPI_Allreduce(&solver_runtime,&m_step_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
    m_step_wall_time = solver_runtime;
#endif
  }

  // update field boundary probes, and then advance fields to new time
  m_fields->updateBoundaryProbes( a_dt );
  const Real new_time = a_old_time + a_dt;
  if (m_curl2) {
    m_fields->extrapolateB(m_theta, new_time);
  }
  m_U = (1.0/m_theta)*m_U + ((m_theta-1.0)/m_theta)*m_Uold;
  m_fields->updatePhysicalState( m_U, new_time, (m_curl2?e_only:e_and_b) );

  // update particles from half_time to new time
  const PicChargedSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getChargedPtrVect();
  const int num_species = pic_species_ptr_vect.size();
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    species->advanceVelocities_2ndHalf( a_dt );
    species->advancePositions_2ndHalf();
    species->mergeSubOrbitParticles();
    species->applyBCs( false, new_time);
  }

  return 0;
}

void PICTimeIntegrator_EM_ThetaImplicit::preRHSOp( const ODEVector<EMFields>&  a_U,
                                                   const Real                  a_time,
                                                   const Real                  a_dt,
                                                   const int                   a_nl_iter,
                                                   const bool                  a_from_emjacobian )
{
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::preRHSOp()");

  // This function is called from nonlinear solvers prior to call to
  // computeRHS(). Here, the particles are updated using the current
  // state of the fields on the grid, which are then used to compute
  // the current density on the grid needed for computing the RHS.

  if (m_nlsolver_type != _NLSOLVER_PICARD_) {
    if (m_curl2) {
      m_fields->copyEFromVec( a_U );
      m_fields->applyBCs_electricField( a_time );
      m_fields->updateB(m_theta*a_dt, a_time);
    } else {
      int offset(0);
      m_fields->copyBFromVec( a_U, offset );
      m_fields->copyEFromVec( a_U, offset );
      m_fields->applyBCs_magneticField( a_time );
      m_fields->applyBCs_electricField( a_time );
    }
  }

  // update xp and vp and set the pic current density
  m_pic_species->preRHSOp( a_from_emjacobian, *m_fields, a_dt, m_nl_iter );
  const LevelData<EdgeDataBox>& pic_J = m_pic_species->getCurrentDensity();
  const LevelData<NodeFArrayBox>& pic_Jv = m_pic_species->getVirtualCurrentDensity();
  m_fields->setCurrentDensity( pic_J, pic_Jv );

  // a_nl_iter when using Petsc comes here with zero twice, and passing it
  // to m_pic_species->preRHSOp() doesn't work with the usage of this parameter,
  // which is to do a modified particle update on the iteration

  if(!a_from_emjacobian) m_nl_iter += 1;
  m_prerhsop_call = true;

}

void PICTimeIntegrator_EM_ThetaImplicit::computeRHS( ODEVector<EMFields>&  a_F,
                                               const ODEVector<EMFields>&,
                                               const Real                  a_time,
                                               const Real                  a_dt,
                                               const int                   a_block )
{
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::computeRHS()");
  CH_assert(m_prerhsop_call || (m_nlsolver_type == _NLSOLVER_PICARD_));

  // this function is called from the nonlinear solvers

#if CH_SPACEDIM==1
   const Real full_dt = a_dt/m_theta;
   m_fields->applyAbsorbingBCs( a_time, full_dt );
#endif

  if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    CH_assert(!m_curl2);
    if (a_block == _EM_PICARD_B_FIELD_BLOCK_) {
      m_fields->computeRHSMagneticField( a_dt );
      int offset = m_fields->vecOffsetBField();
      m_fields->copyBRHSToVec( a_F, offset );
    } else if (a_block == _EM_PICARD_E_FIELD_BLOCK_) {
      m_fields->computeRHSElectricField( a_dt );
      int offset = m_fields->vecOffsetEField();
      m_fields->copyERHSToVec( a_F, offset );
    }

  } else {

     m_fields->zeroRHS();
     if (m_curl2) {
       m_fields->computeRHSElectricField( a_dt );
       m_fields->copyERHSToVec( a_F );
     } else {
       m_fields->computeRHSMagneticField( a_dt );
       m_fields->computeRHSElectricField( a_dt );
       int offset(0);
       m_fields->copyBRHSToVec( a_F, offset );
       m_fields->copyERHSToVec( a_F, offset );
     }

  }

  m_prerhsop_call = false;
}

void PICTimeIntegrator_EM_ThetaImplicit::updatePhysicalState( ODEVector<EMFields>&  a_U,
                                                        const int                   a_block,
                                                        const Real                  a_time )
{
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::updatePhysicalState() from Picard solver");
  CH_assert(!m_curl2);

  if (a_block == _EM_PICARD_B_FIELD_BLOCK_) {
    int offset = m_fields->vecOffsetBField();
    m_fields->copyBFromVec( a_U, offset );
    m_fields->applyBCs_magneticField( a_time );
    offset = m_fields->vecOffsetBField();
    m_fields->copyBToVec( a_U, offset );
  } else if (a_block == _EM_PICARD_E_FIELD_BLOCK_) {
    int offset = m_fields->vecOffsetEField();
    m_fields->copyEFromVec( a_U, offset );
    m_fields->applyBCs_electricField( a_time );
    offset = m_fields->vecOffsetEField();
    m_fields->copyEToVec( a_U, offset );
  }

}

void PICTimeIntegrator_EM_ThetaImplicit::updatePrecondMat( BandedMatrix& a_Pmat,
                                                     const Real    a_time,
                                                     const Real    a_dt )
{
  CH_TIME("PICTimeIntegrator_EM_ThetaImplicit::updatePrecondMat()");

  a_Pmat.setToIdentityMatrix();
  if(m_pic_species->getSigmaxx().isDefined()) {
     m_pic_species->setMassMatricesForPC(*m_fields);
  }
  m_fields->assemblePrecondMatrix( a_Pmat,
                                   m_pic_species->getSigmaxx().isDefined(),
                                   a_dt,
                                   m_vec_type );
  a_Pmat.finalAssembly();
  //a_Pmat.writeToFile("pc_matrix.txt",1);
}

#include "NamespaceFooter.H"

