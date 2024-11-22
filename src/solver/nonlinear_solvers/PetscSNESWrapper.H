/* TODO
 - what is the purpose of m_iter? and function curIteration()?
*/

#ifndef _PetscSNESWrapper_H_
#define _PetscSNESWrapper_H_

#define _NLSOLVER_PETSCSNES_ "petsc_snes" 

#ifdef with_petsc

#include <iostream>
#include <string>

#include <petscsnes.h>
#include "NonlinearSolver.H"
#include "PetscCopiers.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class PetscSNESWrapper : public NonlinearSolver<T, Ops>
{
  public:

    PetscSNESWrapper<T,Ops>()
    {
      m_is_Defined = false;
      m_func = nullptr;
      m_op_indent = "";
      m_iter = -1;
      m_count = 0;
      m_count_lin = 0;

      /* default values */
      m_newton_maxits  = 100;
      m_newton_rtol    = 1.0e-6;
      m_newton_atol    = 0.0;
      m_newton_stol    = 1.0e-14;
      m_newton_maxstep = 1.0e14;
      m_newton_maxfail = 10;
      m_gmres_maxits   = 1000;
      m_gmres_rtol     = 1.0e-3;
      m_gmres_atol     = 0.0;

      m_newton_verbose = false;
      m_gmres_verbose = false;

      m_pc_type = "asm";
      if (numProc()==1) { m_pc_type = "ilu"; }
      m_pc_asm_overlap = 2; /*TODO: make this local size-dependent */
      m_pc_factor_levels = 2;
      m_sub_pc_type = "ilu";
      m_sub_pc_factor_levels = 2;

      m_linesearch_type = "basic"; // default for SNES in PETSc is backtracing (bt)

    }

    ~PetscSNESWrapper<T,Ops>()
    {
      SNESDestroy(&m_nlsolver);
      MatDestroy(&m_J);
      VecDestroy(&m_Y);
      return;
    }

    void define(  const T&, 
                  Ops* const, 
                  NonlinearFunction<T,Ops>* const,
                  const double );

    inline void solve(  T&        a_U, 
                        const T&  a_b, 
                        const double, 
                        const double ) const
    {
      CH_TIME("PetscSNESWrapper::solve()");
      CH_assert(isDefined());

      PetscCopiers::cpAppVecToPetsc<T>( m_Y, a_U );
      PetscCopiers::cpAppVecToPetsc<T>( m_R, a_b );
      SNESSolve( m_nlsolver, m_R, m_Y );
      PetscCopiers::cpAppVecFromPetsc<T>( m_Y, a_U );

      // update iteration counts 
      PetscInt nl_iter;
      SNESGetIterationNumber( m_nlsolver, &nl_iter );
      m_count += (int)nl_iter;
       
      PetscInt l_iter;
      SNESGetLinearSolveIterations( m_nlsolver, &l_iter );
      m_count_lin += (int)l_iter;

      SNESConvergedReason reason;
      SNESGetConvergedReason( m_nlsolver, &reason );
      this->m_exit_status = (int)reason;
      if(m_newton_maxits > 2) { // if maxits==2, using ec_semi_implicit
        if(!procID() && this->m_exit_status < 0) {
          cout << "SNES solver diverged: exit status = " << this->m_exit_status << endl;
          cout << "see https://petsc.org/release/docs/manualpages/SNES/SNESConvergedReason/" << endl;
        }
      }

      return;
    }


    bool isDefined() const
    {
      return m_is_Defined;
    }

    inline void getParams( double&  a_rtol,
                           double&  a_atol,
                           int&     a_maxits )
    {
      a_rtol = m_newton_rtol;
      a_atol = m_newton_atol;
      a_maxits = m_newton_maxits;
    }

    void printParams( ) const
    {
       if (procID()>0) { return; }

       PetscInt max_fail;
       SNESGetMaxNonlinearStepFailures( m_nlsolver, &max_fail );

       cout << " Petsc SNES solver parameters: " << endl;
       cout << "   solver_type = NEWTONLS" << endl;
       cout << "   linesearch_type = " << m_linesearch_type << endl;
       cout << "   newton_rtol     = " << m_newton_rtol << endl;
       cout << "   newton_atol     = " << m_newton_atol << endl;
       cout << "   newton_stol     = " << m_newton_stol << endl;
       cout << "   newton_maxits   = " << m_newton_maxits << endl;
       cout << "   newton_maxstep  = " << m_newton_maxstep << endl;
       cout << "   newton_maxfail  = " << (int)max_fail << endl;
       //cout << "   newton_maxfail = " << m_newton_maxfail << endl;
       cout << "   gmres_atol   = " << m_gmres_atol << endl;
       cout << "   gmres_rtol   = " << m_gmres_rtol << endl;
       cout << "   gmres_maxits = " << m_gmres_maxits << endl;
       cout << "   pc_type        = " << m_pc_type << endl;
       if (m_pc_type=="asm") {
           cout << "   pc_asm_overlap = " << m_pc_asm_overlap << endl;
           cout << "   sub_pc_type    = " << m_sub_pc_type << endl;
           if (m_sub_pc_type=="ilu") {
               cout << "   sub_pc_factor_levels = " << m_sub_pc_factor_levels << endl;
           }
       }
       else if (m_pc_type=="ilu") {
          cout << "   pc_factor_levels = " << m_pc_factor_levels << endl;
       }
    }

    virtual void verbose( const bool  a_verbose )
    { /*
       m_newton_verbose = a_verbose;
       m_gmres_verbose = a_verbose;
       if(m_gmres_verbose) {
         PetscOptionsSetValue( NULL, "-ksp_monitor" , NULL );
         PetscOptionsSetValue( NULL, "-ksp_monitor" , NULL );
         PetscOptionsSetValue( NULL, "-ksp_converged_reason" , NULL );
       }
       else {
         PetscOptionsClearValue( NULL, "-ksp_monitor" );
         PetscOptionsClearValue( NULL, "-ksp_converged_reason" );
       }
      */
    }

    virtual void setOutputIndent( const std::string& a_indent )
    { 
      m_op_indent = a_indent;
    }

    inline const std::string& opIndent() const
    {
      return m_op_indent;
    }

    inline NonlinearFunction<T,Ops>* getFunction()
    {
      return m_func;
    }

    inline int curIteration() const
    {
      return m_iter;
    }

    inline void curIteration(int a_iter)
    {
      m_iter = a_iter;
    }

    inline void norm0(double a_nl_norm0)
    {
      m_nl_norm0 = a_nl_norm0;
    }

    inline bool PCMatAvailable() const
    {
      return m_pcmat_available;
    }

    inline int updatePCMatrix(Mat a_P)
    {
      LinearFunction<T,Ops>& jac( m_func->getJacobian() );
      int retval = jac.copyPCMatrixToPetsc( a_P );
      return retval;
    }

    mutable T m_U, m_F;

    virtual void getConvergenceParams( int& a_l_exit_status, int& a_l_last_iter, int& a_l_total_iter,
                                       int& a_nl_exit_status, int& a_nl_iter, int& a_nl_total_iter,
                                    double& a_nl_abs_res, double& a_nl_rel_res ) const
    {
       // get nonlinear solver params (https://petsc.org/release/docs/manualpages/SNES/SNESConvergedReason/)
       SNESConvergedReason reason;
       SNESGetConvergedReason( m_nlsolver, &reason );
       a_nl_exit_status = (int)reason;

       PetscInt nl_iter;
       SNESGetIterationNumber( m_nlsolver, &nl_iter );
       a_nl_iter = (int)nl_iter;

       PetscReal nl_norm;
       SNESGetFunctionNorm( m_nlsolver, &nl_norm);

       // get linear solver params
       KSP lin_solver;
       SNESGetKSP( m_nlsolver, &lin_solver );

       KSPConvergedReason l_reason;
       KSPGetConvergedReason( lin_solver, &l_reason);
       a_l_exit_status = (int)l_reason;

       PetscInt l_iter;
       KSPGetIterationNumber( lin_solver, &l_iter );
       a_l_last_iter = (int)l_iter;

       a_l_total_iter = m_count_lin;
       a_nl_total_iter = m_count;
       a_nl_abs_res = (Real)nl_norm;
       a_nl_rel_res = (Real)nl_norm/m_nl_norm0;
    }

  protected:

    void parseParametersNewton( ParmParse& );
    void parseParametersGMRES( ParmParse& );
    void parseParametersPC( ParmParse& );

    void setOptions();

    bool m_is_Defined;

    int m_iter;

    int m_newton_maxits;
    double m_newton_rtol, m_newton_atol, m_newton_stol;
    double m_newton_maxstep;
    int m_newton_maxfail;
    bool m_newton_verbose;

    std::string m_pc_type;
    std::string m_sub_pc_type;
    int m_pc_asm_overlap;
    int m_pc_factor_levels;
    int m_sub_pc_factor_levels;

    int m_gmres_maxits;
    double m_gmres_rtol, m_gmres_atol;
    bool m_gmres_verbose;

    std::string m_op_indent;
    bool m_pcmat_available;

    NonlinearFunction<T,Ops>* m_func;

    mutable Vec m_Y, m_R;
    mutable Mat m_J, m_P;
    SNES m_nlsolver;
    std::string m_linesearch_type;

  private:
    
    mutable int m_count, m_count_lin;
    mutable double m_nl_norm0;
    
};

template <class T, class Ops>
static PetscErrorCode RHSFunction(  SNES  solver,
                                    Vec   a_Y,
                                    Vec   a_F,
                                    void* ctxt  )
{
  CH_TIME("PetscErrorCode::RHSFunction()");

  PetscSNESWrapper<T,Ops>* context
      = (PetscSNESWrapper<T,Ops>*) ctxt;
  NonlinearFunction<T,Ops>* func( context->getFunction() );

  int its;
  SNESGetIterationNumber( solver, &its );
  context->curIteration(its);

  if(its==0) {
     PetscReal nl_norm;
     SNESGetFunctionNorm( solver, &nl_norm);
     context->norm0(nl_norm);
  }

  PetscCopiers::cpAppVecFromPetsc<T>( a_Y, context->m_U );
  func->evalFunction( context->m_F, context->m_U, its );
  PetscCopiers::cpAppVecToPetsc<T>( a_F, context->m_F);

  return 0;
}

template <class T, class Ops> 
static PetscErrorCode JacobianFunction( SNES  a_nlsolver,
                                        Vec   a_Y,
                                        Mat   a_A,
                                        Mat   a_P,
                                        void* ctxt  )
{
  CH_TIME("PetscErrorCode::JacobianFunction()");

  PetscSNESWrapper<T,Ops> *context 
      = (PetscSNESWrapper<T,Ops>*) ctxt;

  KSP lin_solver;
  SNESGetKSP( a_nlsolver, &lin_solver );
  PC pc;
  KSPGetPC( lin_solver, &pc );
  PCType pctype;
  PCGetType( pc, &pctype );

  if (strcmp(pctype,PCNONE)) {
    PetscCopiers::cpAppVecFromPetsc<T>( a_Y, context->m_U );
    int err = context->updatePCMatrix(a_P);
    if (err) {
      if (!procID()) {
        printf( "%s(%s) Preconditioning matrix unavailable.\n",
                context->opIndent().c_str(), _NLSOLVER_PETSCSNES_ );
        printf( "%s(%s) Turning off preconditioning.\n",
                context->opIndent().c_str(), _NLSOLVER_PETSCSNES_ );
      }
      PCSetType( pc, PCNONE );
    }
  }

  return 0;
}

template <class T, class Ops> 
static PetscErrorCode applyJacobian( Mat a_A,
                                     Vec a_Y,
                                     Vec a_F )
{
  CH_TIME("PetscErrorCode::applyJacobian()");

  PetscSNESWrapper<T,Ops> *context;
  MatShellGetContext(a_A,&context);
  NonlinearFunction<T,Ops>* func( context->getFunction() );
  LinearFunction<T,Ops>* jac( &func->getJacobian() );

  PetscCopiers::cpAppVecFromPetsc<T>( a_Y, context->m_U );
  jac->applyOp( context->m_F, context->m_U, context->curIteration() );
  PetscCopiers::cpAppVecToPetsc<T>( a_F, context->m_F);

  return 0;
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::define( const T&                        a_U,
                                      Ops* const                      a_ops,
                                      NonlinearFunction<T,Ops>* const a_func,
                                      const double )
{
  CH_assert(!isDefined());

  ParmParse pp_newton ( "newton" );
  parseParametersNewton( pp_newton );

  ParmParse pp_gmres ( "gmres" );
  parseParametersGMRES( pp_gmres );

  ParmParse pp_pc ( "pc" );
  parseParametersPC( pp_pc );

  m_U.define(a_U);
  m_F.define(a_U);

  CH_assert( a_func->isDefined() );
  m_func = a_func;
  LinearFunction<T,Ops>& jac( m_func->getJacobian() );
  jac.updatePreCondMat( a_U );

  int m_local_size = a_U.getVectorSize();
  VecCreate( MPI_COMM_WORLD, &m_Y );
  VecSetSizes( m_Y, m_local_size, PETSC_DECIDE );
  VecSetUp( m_Y );
  VecDuplicate( m_Y, &m_R );

  SNESCreate( MPI_COMM_WORLD, &m_nlsolver );
  SNESSetFunction( m_nlsolver, NULL, RHSFunction<T,Ops>, this);

  MatCreateShell( MPI_COMM_WORLD,
                  m_local_size,
                  m_local_size,
                  PETSC_DETERMINE,
                  PETSC_DETERMINE,
                  this,
                  &m_J );
  MatShellSetOperation( m_J,
                        MATOP_MULT,
                        (void (*)(void))applyJacobian<T,Ops> );
  MatSetUp(m_J);

  MatCreate( MPI_COMM_WORLD, &m_P );
  MatSetSizes(  m_P, 
                m_local_size, m_local_size,
                PETSC_DETERMINE, PETSC_DETERMINE );
  MatSetType( m_P, MATAIJ );
  MatMPIAIJSetPreallocation( m_P, a_ops->numPCMatBands(), NULL, a_ops->numPCMatBands()-1, NULL);
  MatSeqAIJSetPreallocation( m_P, a_ops->numPCMatBands(), NULL);
  MatSetUp(m_P);
  m_pcmat_available = jac.copyPCMatrixToPetsc( m_P );

  SNESSetJacobian( m_nlsolver, m_J, m_P, JacobianFunction<T,Ops>, this);

  // specify nonlinear solver (https://petsc.org/release/docs/manualpages/SNES/SNESType/)
  SNESSetType( m_nlsolver, SNESNEWTONLS );

  // defaut line search method for SNESNEWTONLS is backtracking over L2 norm of function
  if(m_linesearch_type=="basic") { // basic: damping method that defaults to full Newton step
     SNESLineSearch linesearch;
     SNESGetLineSearch( m_nlsolver, &linesearch );
     SNESLineSearchSetType( linesearch, SNESLINESEARCHBASIC );
  }

  setOptions();
  SNESSetFromOptions( m_nlsolver );

  m_is_Defined = true;
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::setOptions()
{
  // the is_specified flag is to check if this parameter has already
  // been set from the .petscrc file, which has priority
  {
    char option[100] = "-snes_atol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_newton_atol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-snes_rtol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_newton_rtol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-snes_stol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_newton_stol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-snes_linesearch_maxstep";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_newton_maxstep );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-snes_max_fail";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%d", m_newton_maxfail );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  //SNESSetMaxNonlinearStepFailures( m_nlsolver, m_newton_maxfail );
  {
    char option[100] = "-snes_max_it";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%d", m_newton_maxits );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-snes_monitor";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified && m_newton_verbose) {
      //PetscOptionsSetValue( NULL, "-snes_view_pre" , NULL );
      PetscOptionsSetValue( NULL, option , NULL );
      PetscOptionsSetValue( NULL, "-snes_converged_reason" , NULL );
    }
  }

  {
    char option[100] = "-ksp_atol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_gmres_atol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-ksp_rtol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_gmres_rtol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-ksp_max_it";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%d", m_gmres_maxits );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-ksp_monitor";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      if (m_gmres_verbose) {
        PetscOptionsSetValue( NULL, option , NULL );
        PetscOptionsSetValue( NULL, "-ksp_converged_reason" , NULL );
      }
    }
  }

  // Set options for PETSc Preconditioner. Priorty is given to values set
  // in the .petscrc file. Some examples are below:
  // -pc_type lu
  //
  // -pc_type ilu
  // -pc_factor_levels 16
  //
  // -pc_type asm
  // -pc_asm_overlap 16
  // -sub_pc_type ilu
  // -sub_pc_factor_levels 8

  KSP lin_solver;
  SNESGetKSP( m_nlsolver, &lin_solver );
  PC precon;
  KSPGetPC( lin_solver, &precon );

  PetscBool pc_specified;
  PetscOptionsHasName( NULL, NULL, "-pc_type", &pc_specified );
  if (pc_specified) {
    char value[100];
    PetscOptionsGetString( NULL, NULL, "-pc_type", value,
                           sizeof(value), &pc_specified );
    std::string str_value(value);
    m_pc_type = str_value;
  }
  else {
    PetscOptionsSetValue( NULL, "-pc_type", m_pc_type.c_str());
  }

  if (!procID()) {
    printf( "%s(%s) Preconditioner: PETSc's %s.\n",
            m_op_indent.c_str(), _NLSOLVER_PETSCSNES_, m_pc_type.c_str() );
    cout << endl;
  }

  if (m_pc_type=="ilu") {
    char option[100] = "-pc_factor_levels";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (is_specified) {
      PetscOptionsGetInt( NULL, NULL, option, &m_pc_factor_levels, &is_specified );
    }
    else {
      char levels[100];
      sprintf(levels, "%i", m_pc_factor_levels );
      PetscOptionsSetValue( NULL, "-pc_factor_levels", levels);
    }
  }
  else if (m_pc_type=="asm") {
    char option[100] = "-pc_asm_overlap";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (is_specified) {
      PetscOptionsGetInt( NULL, NULL, option, &m_pc_asm_overlap, &is_specified );
    }
    else {
      char overlap[100];
      sprintf(overlap, "%i", m_pc_asm_overlap );
      PetscOptionsSetValue( NULL, "-pc_asm_overlap", overlap);
    }
    //
    PetscOptionsHasName( NULL, NULL, "-sub_pc_type", &is_specified );
    if (is_specified) {
      char value[100];
      PetscOptionsGetString( NULL, NULL, "-sub_pc_type", value,
                             sizeof(value), &pc_specified );
      std::string str_value(value);
      m_sub_pc_type = str_value;
    }
    else {
      PetscOptionsSetValue( NULL, "-sub_pc_type", m_sub_pc_type.c_str());
    }
    //
    if (m_sub_pc_type=="ilu") {
      char option[100] = "-sub_pc_factor_levels";
      PetscBool is_specified;
      PetscOptionsHasName( NULL, NULL, option, &is_specified );
      if (is_specified) {
        PetscOptionsGetInt( NULL, NULL, option, &m_sub_pc_factor_levels, &is_specified );
      }
      else {
        char levels[100];
        sprintf(levels, "%i", m_sub_pc_factor_levels );
        PetscOptionsSetValue( NULL, "-sub_pc_factor_levels", levels);
      }
    }
  }

  return;
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::parseParametersNewton( ParmParse& a_pp )
{
  a_pp.query("verbose", m_newton_verbose); /* verbosity */
  a_pp.query("line_search.type", m_linesearch_type); /* linesearch type */
  CH_assert(m_linesearch_type=="basic" || m_linesearch_type=="bt");
  a_pp.query("rtol",    m_newton_rtol);    /* relative tolerance */
  a_pp.query("atol",    m_newton_atol);    /* absolute tolerance */
  a_pp.query("stol",    m_newton_stol);    /* step tolerance */
  a_pp.query("maxits",  m_newton_maxits);  /* maximum iterations */
  a_pp.query("maxstep", m_newton_maxstep); /* max step for line search */
  a_pp.query("maxfail", m_newton_maxfail); /* maximum step failures */
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::parseParametersGMRES( ParmParse& a_pp )
{
  a_pp.query("verbose", m_gmres_verbose); /* verbosity */
  a_pp.query("rtol",    m_gmres_rtol);    /* relative tolerance */
  a_pp.query("atol",    m_gmres_atol);    /* absolute tolerance */
  a_pp.query("maxits",  m_gmres_maxits);  /* maximum iterations */
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::parseParametersPC( ParmParse& a_pp )
{
  a_pp.query("pc_type", m_pc_type);       /* preconditioner type (bjacobi, lu, ilu, asm, none) */
  if (m_pc_type=="asm") {
    a_pp.query("pc_asm_overlap", m_pc_asm_overlap);  /* overlap for asm pc */
    a_pp.query("sub_pc_type", m_sub_pc_type);  /* local linear solver for asm pc (ilu or lu)*/
    if (m_sub_pc_type=="ilu") {
      a_pp.query("sub_pc_factor_levels", m_sub_pc_factor_levels); /* factor levels for ilu */
    }
  }
  else if (m_pc_type=="ilu") {
    a_pp.query("pc_factor_levels", m_pc_factor_levels); /* factor levels for ilu */
  }
}

#include "NamespaceFooter.H"
#endif

#endif
