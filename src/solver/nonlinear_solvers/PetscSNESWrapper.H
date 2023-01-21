#ifndef _PetscSNESWrapper_H_
#define _PetscSNESWrapper_H_

#define _NLSOLVER_PETSC_ "petsc"

#ifdef with_petsc

#include <iostream>
#include <string>

#include <petscsnes.h>
#include "NonlinearSolver.H"

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
      m_newton_maxits = 10;
      m_newton_rtol   = 1e-6;
      m_newton_atol   = 1e-10;
      m_newton_stol   = 1e-14;
      m_gmres_maxits  = 100;
      m_gmres_rtol    = 1e-3;
      m_gmres_atol    = 1e-12;

      m_newton_verbose = true;
      m_gmres_verbose = true;

    }

    ~PetscSNESWrapper<T,Ops>()
    {
      SNESDestroy(&m_nlsolver);
      MatDestroy(&m_J);
      VecDestroy(&m_Y);
      return;
    }

    void cpAppVecToPetsc( Vec&, const T& ) const;
    void cpAppVecFromPetsc( const Vec&, T& ) const;

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
    
      cpAppVecToPetsc( m_Y, a_U );
      cpAppVecToPetsc( m_R, a_b );
      SNESSolve( m_nlsolver, m_R, m_Y );
      cpAppVecFromPetsc( m_Y, a_U );
      
      // update iteration counts 
      PetscInt nl_iter;
      SNESGetIterationNumber( m_nlsolver, &nl_iter );
      m_count += (int)nl_iter;
       
      PetscInt l_iter;
      SNESGetLinearSolveIterations( m_nlsolver, &l_iter );
      m_count_lin += (int)l_iter;

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
       if(procID()>0) return;
       cout << " Petsc SNES solver parameters: " << endl;
       cout << "   newton_atol   = " << m_newton_atol << endl;
       cout << "   newton_rtol   = " << m_newton_rtol << endl;
       cout << "   newton_maxits = " << m_newton_maxits << endl;
       cout << "   gmres_atol   = " << m_gmres_atol << endl;
       cout << "   gmres_rtol   = " << m_gmres_rtol << endl;
       cout << "   gmres_maxits = " << m_gmres_maxits << endl;
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
       // get nonlinear solver params (https://petsc.org/main/docs/manualpages/SNES/SNESConvergedReason/)
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

    void setOptions();

    bool m_is_Defined;

    int m_iter;

    int m_newton_maxits;
    double m_newton_rtol, m_newton_atol, m_newton_stol;
    bool m_newton_verbose;

    int m_gmres_maxits;
    double m_gmres_rtol, m_gmres_atol;
    bool m_gmres_verbose;

    std::string m_op_indent;
    bool m_pcmat_available;

    NonlinearFunction<T,Ops>* m_func;  

    mutable Vec m_Y, m_R;
    mutable Mat m_J, m_P;
    SNES m_nlsolver;

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

  context->cpAppVecFromPetsc( a_Y, context->m_U );
  func->evalFunction( context->m_F, context->m_U, its );
  context->cpAppVecToPetsc( a_F, context->m_F);

  return 0;
}

template <class T, class Ops> 
static PetscErrorCode JacobianFunction( SNES  a_solver,
                                        Vec   a_Y, 
                                        Mat   a_A, 
                                        Mat   a_P, 
                                        void* ctxt  )
{
  CH_TIME("PetscErrorCode::JacobianFunction()");
  
  PetscSNESWrapper<T,Ops> *context 
      = (PetscSNESWrapper<T,Ops>*) ctxt;

  KSP lin_solver;
  SNESGetKSP( a_solver, &lin_solver );
  PC pc;
  KSPGetPC( lin_solver, &pc );
  PCType pctype;
  PCGetType( pc, &pctype );

  if (strcmp(pctype,PCNONE)) {
    context->cpAppVecFromPetsc( a_Y, context->m_U );
    int err = context->updatePCMatrix(a_P);
    if (err) {
      if (!procID()) {
        printf( "%s(%s) Preconditioning matrix unavailable.\n",
                context->opIndent().c_str(), _NLSOLVER_PETSC_ );
        printf( "%s(%s) Turning off preconditioning.\n",
                context->opIndent().c_str(), _NLSOLVER_PETSC_ );
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

  context->cpAppVecFromPetsc( a_Y, context->m_U );
  jac->applyOp( context->m_F, context->m_U, context->curIteration() );
  context->cpAppVecToPetsc( a_F, context->m_F);

  return 0;
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::cpAppVecToPetsc( Vec&     a_Y, 
                                               const T& a_U ) const
{
  PetscScalar *Yarr;
  VecGetArray(a_Y,&Yarr);
  a_U.copyTo( (Real* const) Yarr );
  VecRestoreArray(a_Y,&Yarr);
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::cpAppVecFromPetsc( const Vec&  a_Y, 
                                                 T&          a_U ) const
{
  const PetscScalar *Yarr;
  VecGetArrayRead(a_Y,&Yarr);
  a_U.copyFrom( (const Real* const) Yarr );
  VecRestoreArrayRead(a_Y,&Yarr);
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::define( const T&                        a_U,
                                      Ops* const                      a_ops, 
                                      NonlinearFunction<T,Ops>* const a_func,
                                      const double )
{
  CH_assert(!isDefined());

  ParmParse ppnewton ( "newton" );
  parseParametersNewton( ppnewton );

  ParmParse ppgmres ( "gmres" );
  parseParametersGMRES( ppgmres );

  m_U.define(a_U);
  m_F.define(a_U);

  CH_assert( a_func->isDefined() );
  m_func = a_func;
  LinearFunction<T,Ops>& jac( m_func->getJacobian() );
  jac.updatePreCondMat( a_U, 0 );

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

  SNESSetType( m_nlsolver, SNESNEWTONLS );

  setOptions();
  SNESSetFromOptions( m_nlsolver );

  m_is_Defined = true;
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::setOptions()
{
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
    if (!is_specified) {
      if (m_newton_verbose) {
        PetscOptionsSetValue( NULL, option , NULL );
        PetscOptionsSetValue( NULL, "-snes_converged_reason" , NULL );
      }
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

  PetscBool pc_specified;
  PetscOptionsHasName( NULL, NULL, "-pc_type", &pc_specified );

  if (!pc_specified) {

    KSP lin_solver;
    SNESGetKSP( m_nlsolver, &lin_solver );
    PC precon;
    KSPGetPC( lin_solver, &precon );

    if (numProc() == 1) {
      if (SpaceDim == 1) {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: PETSc's native ILU (serial).\n",
                  m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "ilu");
        PetscOptionsSetValue( NULL, "-pc_factor_levels", "2");
      } else {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: PETSc's native LU (serial).\n",
                  m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "lu");
      }
    } else {
#if defined(PETSC_HAVE_HYPRE)
      if (SpaceDim == 1) {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: hypre's parallel ILU (euclid).\n",
                  m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "hypre");
        PetscOptionsSetValue( NULL, "-pc_hypre_type", "euclid");
        PetscOptionsSetValue( NULL, "-pc_hypre_euclid_level", "2");
      } else {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: None.\n",
                  m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "none");
      }
#elif defined(PETSC_HAVE_SUPERLU_DIST)
      if (SpaceDim == 1) {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: SuperLU_dist's parallel LU.\n",
                  m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "lu");
      } else {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: None.\n",
                  m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "none");
      }
#else
      if (!procID()) {
        printf( "%s(%s)Preconditioner: None available.\n",
                m_op_indent.c_str(), _NLSOLVER_PETSC_ );
        printf( "%s(%s)Preconditioner: Please configure PETSc with an appropriate package.\n",
                m_op_indent.c_str(), _NLSOLVER_PETSC_ );
      }
      PetscOptionsSetValue( NULL, "-pc_type", "none");
#endif
    }

  }

  return;
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::parseParametersNewton( ParmParse& a_pp )
{
  a_pp.query("verbose", m_newton_verbose); /* verbosity */
  a_pp.query("atol",    m_newton_atol);    /* absolute tolerance */
  a_pp.query("rtol",    m_newton_rtol);    /* relative tolerance */
  a_pp.query("maxits",  m_newton_maxits);  /* maximum iterations */
}

template <class T, class Ops>
void PetscSNESWrapper<T,Ops>::parseParametersGMRES( ParmParse& a_pp )
{
  a_pp.query("verbose", m_gmres_verbose); /* verbosity */
  a_pp.query("atol",    m_gmres_atol);    /* absolute tolerance */
  a_pp.query("rtol",    m_gmres_rtol);    /* relative tolerance */
  a_pp.query("maxits",  m_gmres_maxits);  /* maximum iterations */
}

#include "NamespaceFooter.H"
#endif

#endif
