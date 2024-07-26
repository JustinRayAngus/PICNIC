#include "PICTimeIntegrator_EM_SemiImplicit.H"
#include "CH_HDF5.H"
#include <sys/time.h>

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_SemiImplicit::define( System* const               a_system,
                                                PicSpeciesInterface* const  a_pic_species,
                                                EMFields* const             a_fields )
{
  CH_assert(!isDefined());

  m_pic_species = a_pic_species;
  m_fields = a_fields;
  m_fields->computePrecondMatrixNNZ(e_only);

  m_E.define(m_fields->getVectorSize(e_only));
  m_Eold.define(m_E);
  m_E.printLoadBalanceInfo();

  m_B.define(m_fields->getVectorSize(b_only));
  m_Bold.define(m_B);
  m_FB.define(m_B);

  ParmParse pp("pic_em_semi_implicit");
  pp.query("solver_type", m_nlsolver_type);
  pp.query("pc_update_freq", m_pc_update_freq);
  pp.query("pc_update_newton", m_pc_update_newton);

  //backward compatibility
  if (m_nlsolver_type=="petsc") m_nlsolver_type = _NLSOLVER_PETSCSNES_;

  if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    m_nlsolver = new PicardSolver<ODEVector<EMFields>, PICTimeIntegrator>;
  
  } else if (m_nlsolver_type == _NLSOLVER_NEWTON_) {

    m_nlsolver = new NewtonSolver<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func = new EMResidualFunction<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func->define( m_E, this, 0.5, m_pc_update_newton );

  } else if (m_nlsolver_type == _NLSOLVER_PETSCSNES_) {

#ifdef with_petsc
    m_nlsolver = new PetscSNESWrapper<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func = new EMResidualFunction<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func->define( m_E, this, 0.5, m_pc_update_newton );
#else
    MayDay::Error("PETSc SNES solver requires compilation with PETSc");
#endif

  } else {
    MayDay::Error("Invalid choice for PICTimeIntegrator_EM_SemiImplicit::m_nlsolver_type");
  }
  m_nlsolver->define(m_E, this, m_func, 0.5);
  m_nlsolver->setOutputIndent("  ");

  printParams();
  m_is_defined = true;
}
    
void PICTimeIntegrator_EM_SemiImplicit::printParams() const
{
   if(procID()>0) return;
   cout << "================== Time Solver ==================" << endl;
   cout << "advance method = PIC_EM_SEMI_IMPLICIT" << endl;
   cout << "solver_type = " << m_nlsolver_type << endl;
   cout << "pc_update_freq = " << m_pc_update_freq << endl;
   cout << "pc_update_newton = " << (m_pc_update_newton?"true":"false") << endl;
   m_nlsolver->printParams();
   cout << "=================================================" << endl;
   cout << endl;
}

void PICTimeIntegrator_EM_SemiImplicit::initialize( const std::string&  a_restart_file_name )
{
  m_fields->copyEToVec( m_E );
  m_fields->copyBToVec( m_B );
  m_fields->copyEoldToVec( m_Eold );
  m_fields->copyBoldToVec( m_Bold );

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

int PICTimeIntegrator_EM_SemiImplicit::prepForCheckpoint() const
{
  m_fields->copyEoldFromVec( m_Eold );
  m_fields->copyBoldFromVec( m_Bold );
  return 1;
}

void PICTimeIntegrator_EM_SemiImplicit::preTimeStep(  const Real a_time,
                                                      const Real a_dt,
                                                      const int  a_step )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::preTimeStep()");
  
  const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
  const int num_species = pic_species_ptr_vect.size();
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
     
  // advance magnetic field to half time-step ahead
  if (a_step == 0) {
    const Real halfDt = 0.5*a_dt;
    m_fields->computeRHSMagneticField( halfDt );
    m_fields->copyBRHSToVec( m_FB );
    m_B = m_Bold + m_FB;
  } 

  m_Eold = m_E;
  m_Bold = m_B;
  m_fields->setBoldFromB(); // needed for energy diagnostic
  m_fields->setEoldFromE(); // needed for absorbing BC
   
  // update old particle values and create inflow particles  
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    species->updateOldParticlePositions();
    species->updateOldParticleVelocities();
    species->createInflowParticles( a_time, a_dt );
    species->injectInflowParticles();
  }

  return;
}

void PICTimeIntegrator_EM_SemiImplicit::timeStep( const Real a_old_time,
                                                  const Real a_dt,
                                                  const int  a_step )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::timeStep()");
  
  struct timeval pcassembly_start, pcassembly_end;
  struct timeval solve_start, solve_end;
  
  // advance system variables (Eg, Bg, xp, and vp) to 
  // new_time = old_time + dt. Note that Eg, xp, and vp 
  // live at whole time steps, while Bg lives at half
  
  const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
  const int num_species = pic_species_ptr_vect.size();
  const Real half_time = a_old_time + a_dt/2.0;
  
  gettimeofday(&pcassembly_start,NULL);
  if (m_func) {
    m_func->curTime( half_time );
    m_func->curTimeStep( a_dt );
    if (a_step%m_pc_update_freq == 0 || m_restart_pc_update_flag) {
      dynamic_cast<EMJacobianFunction<ODEVector<EMFields>, PICTimeIntegrator>*>
          ( &(m_func->getJacobian()) )->updatePreCondMat( m_E );
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

  // advance Eg, xp, and vp from old_time to half_time
  gettimeofday(&solve_start,NULL);
  m_nlsolver->solve( m_E, m_Eold, half_time, a_dt );
  gettimeofday(&solve_end,NULL);
  m_fields->updateBoundaryProbes( a_dt );
  
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

  // update Eg from half_time to new_time
  m_E = 2.0*m_E - m_Eold;
  const Real new_time = a_old_time + a_dt;
  m_fields->updatePhysicalState( m_E, new_time, e_only );
  
  // advance of Bg n+1/2 to n+3/2
  m_fields->computeRHSMagneticField( a_dt );
  m_fields->copyBRHSToVec( m_FB );
  m_B = m_Bold + m_FB;
  m_fields->updatePhysicalState( m_B, new_time + a_dt/2.0, b_only );

  // update particles from half_time to new_time
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    species->advanceVelocities_2ndHalf(a_dt);
    species->advancePositions_2ndHalf();
    species->mergeSubOrbitParticles();
    species->applyBCs( false, new_time );
  }

}

void PICTimeIntegrator_EM_SemiImplicit::preRHSOp( const ODEVector<EMFields>&  a_U,
                                                  const Real                  a_time,
                                                  const Real                  a_dt,
                                                  const int                   a_nl_iter,
                                                  const bool                  a_from_emjacobian )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::preRHSOp()");
  
  // This function is called from nonlinear solvers prior to call to 
  // computeRHS(). Here, the particles are updated using the current
  // state of the fields on the grid, which are then used to compute
  // the current density on the grid needed for computing the RHS.

  if (m_nlsolver_type != _NLSOLVER_PICARD_) {
    m_fields->copyEFromVec( a_U );
    m_fields->applyBCs_electricField( a_time );
  }

  // update xp and vp and set the pic current density
  m_pic_species->preRHSOp( a_from_emjacobian, *m_fields, a_dt, m_nl_iter );
  if(m_fields->useFiltering()) { m_pic_species->filterJ( *m_fields, a_time ); }
  const LevelData<EdgeDataBox>& pic_J = m_pic_species->getCurrentDensity();
  const LevelData<NodeFArrayBox>& pic_Jv = m_pic_species->getVirtualCurrentDensity();
  m_fields->setCurrentDensity( pic_J, pic_Jv );
  
  // a_nl_iter when using Petsc comes here with zero twice, and passing it
  // to m_pic_species->preRHSOp() doesn't work with the usage of this parameter,
  // which is to do a modified particle update on the iteration
  
  if(!a_from_emjacobian) m_nl_iter += 1;
     
}

void PICTimeIntegrator_EM_SemiImplicit::computeRHS( ODEVector<EMFields>&  a_F,
                                              const ODEVector<EMFields>&,
                                              const Real                  a_time,
                                              const Real                  a_dt,
					      const int )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::computeRHS()");

  // this function is called from the nonlinear solver

#if CH_SPACEDIM==1
  m_fields->applyAbsorbingBCs( a_time, 2.0*a_dt ); // a_dt here is actually a_dt/2.0
#endif

  m_fields->computeRHSElectricField( a_dt );
  m_fields->copyERHSToVec( a_F );

}

void PICTimeIntegrator_EM_SemiImplicit::updatePhysicalState( ODEVector<EMFields>&  a_U,
                                                       const int,
                                                       const Real                  a_time )
{
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::updatePhysicalState() from Picard solver");
  
  m_fields->updatePhysicalState( a_U, a_time, e_only );

}

void PICTimeIntegrator_EM_SemiImplicit::updatePrecondMat( BandedMatrix& a_Pmat,
                                                    const Real    a_time,
                                                    const Real    a_dt )
{
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::updatePrecondMat()");

  a_Pmat.setToIdentityMatrix();
  if(m_pic_species->getSigmaxx().isDefined()) {
     m_pic_species->setMassMatricesForPC(*m_fields);
  }
  m_fields->assemblePrecondMatrix( a_Pmat,
                                   m_pic_species->getSigmaxx().isDefined(),
                                   a_dt,
                                   e_only );
  a_Pmat.finalAssembly();
  //a_Pmat.writeToFile("pc_matrix.txt",1);
}

#include "NamespaceFooter.H"

