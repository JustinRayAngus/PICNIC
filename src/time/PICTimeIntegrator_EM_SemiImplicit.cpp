#include "PICTimeIntegrator_EM_SemiImplicit.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_SemiImplicit::define( PicSpeciesInterface* const  a_pic_species,
                                                EMFields* const             a_fields )
{
  CH_assert(!isDefined());

  m_pic_species = a_pic_species;
  m_fields = a_fields;

  m_E.define(m_fields->getVectorSize(e_only));
  m_Eold.define(m_E);
  m_E.printLoadBalanceInfo();

  m_B.define(m_fields->getVectorSize(b_only));
  m_Bold.define(m_B);
  m_FB.define(m_B);

  ParmParse pp("pic_em_semi_implicit");
  pp.query("solver_type", m_nlsolver_type);

  if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    m_nlsolver = new PicardSolver<ODEVector<EMFields>, PICTimeIntegrator>;
  
  } else if (m_nlsolver_type == _NLSOLVER_NEWTON_) {

    m_nlsolver = new NewtonSolver<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func = new EMResidualFunction<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func->define( m_E, this, 0.5, false );

  } else if (     (m_nlsolver_type == _NLSOLVER_PETSCSNES_)
              ||  (m_nlsolver_type == "petsc") /* backward compatibility */) {

  if(m_nlsolver_type=="petsc") m_nlsolver_type = _NLSOLVER_PETSCSNES_;

#ifdef with_petsc
    m_nlsolver = new PetscSNESWrapper<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func = new EMResidualFunction<ODEVector<EMFields>, PICTimeIntegrator>;
    m_func->define( m_E, this, 0.5, false );
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
    m_fields->computeRHSMagneticField( a_time, halfDt );
    m_fields->copyBRHSToVec( m_FB );
    m_B = m_Bold + m_FB;
  } else { m_B = 2.0*m_B - m_Bold; }
  const Real half_time = a_time + a_dt/2.0;
  m_fields->updatePhysicalState( m_B, half_time, b_only );

  m_Eold = m_E;
  m_Bold = m_B;
   
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
                                                  const int )
{  
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::timeStep()");
  
  // advance system variables (Eg, Bg, xp, and vp) to 
  // new_time = old_time + dt. Note that Eg, xp, and vp 
  // live at whole time steps, while Bg lives at half
  
  const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
  const int num_species = pic_species_ptr_vect.size();
  const Real half_time = a_old_time + a_dt/2.0;
  
  if (m_func) {
    m_func->curTime( half_time );
    m_func->curTimeStep( a_dt );
  }

  // advance Eg, xp, and vp from old_time to half_time
  m_nlsolver->solve( m_E, m_Eold, half_time, a_dt );
  m_fields->updateBoundaryProbes( a_dt );
  
  // update Eg from half_time to new_time
  m_E = 2.0*m_E - m_Eold;
  const Real new_time = a_old_time + a_dt;
  m_fields->updatePhysicalState( m_E, new_time, e_only );
  
  // perform halfDt advance of Bg to go from half_time to new_time
  const Real halfDt = 0.5*a_dt;
  m_fields->computeRHSMagneticField( new_time, halfDt );
  m_fields->copyBRHSToVec( m_FB );
  m_B = m_Bold + m_FB;
  m_fields->updatePhysicalState( m_B, new_time, b_only );

  // update particles from half_time to new_time
  for (int sp=0; sp<num_species; sp++) {
    auto species(pic_species_ptr_vect[sp]);
    species->advanceVelocities_2ndHalf(a_dt);
    species->advancePositions_2ndHalf();
    species->applyInertialForces(a_dt,true,true);
    species->mergeSubOrbitParticles();
    species->applyBCs(false);
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

  m_fields->computeRHSElectricField( a_time, a_dt );
  m_fields->copyERHSToVec( a_F );

}

void PICTimeIntegrator_EM_SemiImplicit::updatePhysicalState( ODEVector<EMFields>&  a_U,
                                                       const int,
                                                       const Real                  a_time )
{
  CH_TIME("PICTimeIntegrator_EM_SemiImplicit::updatePhysicalState() from Picard solver");
  
  m_fields->updatePhysicalState( a_U, a_time, e_only );

}

#include "NamespaceFooter.H"

