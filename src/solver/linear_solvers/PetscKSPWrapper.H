#ifndef _PetscKSPWrapper_H_
#define _PetscKSPWrapper_H_

#define _LINSOLVER_PETSCKSP_ "petsc_ksp"

#ifdef with_petsc

#include <iostream>
#include <string>

#include <petscksp.h>
#include "LinearSolver.H"
#include "PetscCopiers.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class PetscKSPWrapper : public LinearFunctionSolver<T, Ops>
{
  public:

    PetscKSPWrapper<T,Ops>()
    {
      m_is_Defined = false;
      m_op = nullptr;
      m_iter = -1;
      m_count = 0;
      m_op_indent = "";

      /* default values */
      m_maxits  = 100;
      m_rtol    = 1e-3;
      m_atol    = 1e-12;
      m_verbose = false;
    }

    ~PetscKSPWrapper<T,Ops>()
    {
      MatDestroy(&m_A);
      VecDestroy(&m_Y);
      return;
    }

    void define( Ops* const, LinearFunction<T,Ops>* const, const T& );

    inline void solve(  T&        a_U,
                        const T&  a_b,
                        int )
    {
      CH_TIME("PetscKSPWrapper::solve()");
      CH_assert(isDefined());

      /* copy initial guess and RHS to PETSc's vectors */
      PetscCopiers::cpAppVecToPetsc<T>( m_Y, a_U );
      PetscCopiers::cpAppVecToPetsc<T>( m_R, a_b );

      /* update preconditioner */
      {
        PC pc;
        KSPGetPC(m_solver, &pc);
        PCType pctype;
        PCGetType(pc, &pctype);
        if (strcmp(pctype, PCNONE)) {
          Mat Pmat;
          KSPGetOperators(m_solver, NULL, &Pmat);
          m_op->copyPCMatrixToPetsc( Pmat );
        }
      }

      /* solve */
      KSPSolve( m_solver, m_R, m_Y );

      /* copy solution back */
      PetscCopiers::cpAppVecFromPetsc<T>( m_Y, a_U );
      
      // update iteration counts 
      KSPGetIterationNumber( m_solver, &m_iter );
      m_count += (int)m_iter;
       
      KSPConvergedReason reason;
      KSPGetConvergedReason( m_solver, &reason );
      m_exitStatus = (int)reason;
      if(!procID() && m_exitStatus < 0) {
         cout << "KSP solver diverged: exit status = " << m_exitStatus << endl;
         cout << "see https://petsc.org/release/docs/manualpages/KSP/KSPConvergedReason/" << endl;
      }

      return;
    }


    bool isDefined() const 
    { 
      return m_is_Defined; 
    }

    virtual void verbose( const bool  a_verbose ) 
    { 
//      m_verbose = a_verbose;
//      if(m_verbose) {
//        PetscOptionsSetValue( NULL, "-ksp_monitor" , NULL );
//        PetscOptionsSetValue( NULL, "-ksp_converged_reason" , NULL );
//      }
//      else {
//        PetscOptionsClearValue( NULL, "-ksp_monitor" );
//        PetscOptionsClearValue( NULL, "-ksp_converged_reason" );
//      }
    }

    virtual void setOutputIndent( const std::string& a_indent ) 
    { 
      m_op_indent = a_indent;
    }

    virtual int getConvergenceParams( int& a_last_iter, int& a_total_iter ) const
    { 
      a_last_iter = m_iter;
      a_total_iter = m_count;
      return m_exitStatus;
    }

    virtual void applyOp( T& a_F, const T& a_U)
    {
      m_op->applyOp(a_F, a_U, 0);
    }

    mutable T m_F, m_U;

  protected:

    void parseParameters( ParmParse& );

    void setOptions(const bool);

    bool m_is_Defined;
    int m_iter;
    int m_maxits;
    double m_rtol, m_atol;
    bool m_verbose;

    std::string m_op_indent;
    LinearFunction<T,Ops>* m_op;  

    mutable Vec m_Y, m_R;
    mutable Mat m_A, m_P;
    KSP m_solver;
    std::string m_linesearch_type;

    mutable int m_exitStatus;

  private:
    
    mutable int m_count;
    
};

template <class T, class Ops> 
static PetscErrorCode applyOperator( Mat a_A, 
                                     Vec a_Y, 
                                     Vec a_F )
{
  CH_TIME("PetscErrorCode::applyJacobian()");

  PetscKSPWrapper<T,Ops> *context;
  MatShellGetContext(a_A,&context);

  PetscCopiers::cpAppVecFromPetsc<T>( a_Y, context->m_U );
  context->applyOp( context->m_F, context->m_U );
  PetscCopiers::cpAppVecToPetsc<T>( a_F, context->m_F);

  return 0;
}

template <class T, class Ops>
void PetscKSPWrapper<T,Ops>::define(  Ops* const                   a_ops,
                                      LinearFunction<T,Ops>* const a_op,
                                      const T&                     a_U )
{
  CH_assert(!isDefined());

  ParmParse ppgmres ( "gmres" );
  parseParameters( ppgmres );

  m_U.define(a_U);
  m_F.define(a_U);

  m_op = a_op;
  m_op->updatePreCondMat( a_U );

  int m_local_size = a_U.getVectorSize();
  VecCreate( MPI_COMM_WORLD, &m_Y );
  VecSetSizes( m_Y, m_local_size, PETSC_DECIDE );
  VecSetUp( m_Y );
  VecDuplicate( m_Y, &m_R );

  /* Set up the matrix-free op */
  MatCreateShell( MPI_COMM_WORLD,
                  m_local_size,
                  m_local_size,
                  PETSC_DETERMINE,
                  PETSC_DETERMINE,
                  this,
                  &m_A );
  MatShellSetOperation( m_A,
                        MATOP_MULT,
                        (void (*)(void))applyOperator<T,Ops> );
  MatSetUp(m_A);

  /* set up the PC matrix */
  MatCreate( MPI_COMM_WORLD, &m_P );
  MatSetSizes(  m_P, 
                m_local_size, m_local_size,
                PETSC_DETERMINE, PETSC_DETERMINE );
  MatSetType( m_P, MATAIJ );
  MatMPIAIJSetPreallocation( m_P, a_ops->numPCMatBands(), NULL, a_ops->numPCMatBands()-1, NULL);
  MatSeqAIJSetPreallocation( m_P, a_ops->numPCMatBands(), NULL);
  MatSetOption(m_P, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
  MatSetOption(m_P, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetUp(m_P);
  bool flag_pcmat = m_op->copyPCMatrixToPetsc( m_P );


  /* Set up KSP */
  KSPCreate( MPI_COMM_WORLD, &m_solver );
  KSPSetOperators(m_solver, m_A, m_P);
  setOptions(!flag_pcmat);
  KSPSetFromOptions( m_solver );

  m_is_Defined = true;
}

template <class T, class Ops>
void PetscKSPWrapper<T,Ops>::setOptions(const bool a_flag_pcmat)
{
  {
    char option[100] = "-ksp_atol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_atol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-ksp_rtol";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%e", m_rtol );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-ksp_max_it";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      char val[100];
      sprintf(val, "%d", m_maxits );
      PetscOptionsSetValue( NULL, option , val );
    }
  }
  {
    char option[100] = "-ksp_monitor";
    PetscBool is_specified;
    PetscOptionsHasName( NULL, NULL, option, &is_specified );
    if (!is_specified) {
      if (m_verbose) {
        PetscOptionsSetValue( NULL, option , NULL );
        PetscOptionsSetValue( NULL, "-ksp_converged_reason" , NULL );
      }
    }
  }

  PetscBool pc_specified;
  PetscOptionsHasName( NULL, NULL, "-pc_type", &pc_specified );

  if (!pc_specified) {

    PC precon;
    KSPGetPC( m_solver, &precon );

    if (a_flag_pcmat) {

      if (numProc() == 1) {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: PETSc's native ILU (serial).\n",
                  m_op_indent.c_str(), _LINSOLVER_PETSCKSP_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "ilu");
        PetscOptionsSetValue( NULL, "-pc_factor_levels", "2");
      } else {
        if (!procID()) {
          printf( "%s(%s)Preconditioner: Additive Schwarz method with ILU on each block.\n",
                  m_op_indent.c_str(), _LINSOLVER_PETSCKSP_ );
        }
        PetscOptionsSetValue( NULL, "-pc_type", "asm");
        PetscOptionsSetValue( NULL, "-pc_asm_overlap", "2"); /*TODO: make this local size-dependent */
        PetscOptionsSetValue( NULL, "-sub_pc_type", "ilu");
        PetscOptionsSetValue( NULL, "-sub_pc_factor_levels", "2");
      }

    } else {
      if (!procID()) {
        printf( "%s(%s)Preconditioner: Matrix not available; not using preconditioning.\n",
                m_op_indent.c_str(), _LINSOLVER_PETSCKSP_ );
      }
      PetscOptionsSetValue( NULL, "-pc_type", "none");
    }

  }

  return;
}

template <class T, class Ops>
void PetscKSPWrapper<T,Ops>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("verbose", m_verbose); /* verbosity */
  a_pp.query("rtol",    m_rtol);    /* relative tolerance */
  a_pp.query("atol",    m_atol);    /* absolute tolerance */
  a_pp.query("maxits",  m_maxits);  /* maximum iterations */
}

#include "NamespaceFooter.H"
#endif

#endif
