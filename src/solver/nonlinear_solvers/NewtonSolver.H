#ifndef _NewtonSolver_H_
#define _NewtonSolver_H_

#include <iostream>
#include <string>

#include "NonlinearSolver.H"
#include "LinearFunctionSolverLibrary.H"
#include "EMResidualFunction.H"

#define _NLSOLVER_NEWTON_ "newton"

#include "NamespaceHeader.H"

template <class Vec, class Ops>
class NewtonSolver : public NonlinearSolver<Vec, Ops>
{
  public:

    NewtonSolver<Vec,Ops>()
    {
      m_is_Defined = false;
      m_func = nullptr;
      m_linear_solver = nullptr;
      m_op_indent = "";
      m_linear_solver_type = _LINSOLVER_GMRES_;
    }

    ~NewtonSolver<Vec,Ops>()
    {
      delete m_linear_solver;
      return;
    }


    void define(  const Vec&, 
                  Ops* const, 
                  NonlinearFunction<Vec,Ops>* const,
                  const double );

    void solve( Vec&, 
                const Vec&, 
                const double, 
                const double ) const;

    bool isDefined() const 
    { 
      return m_is_Defined; 
    }

    inline void getParams( double&  a_rtol,
                           double&  a_atol,
                           int&     a_maxits )
    {
      a_rtol = m_rtol;
      a_atol = m_atol;
      a_maxits = m_maxits;
    }
    
    void printParams( ) const
    {
       if(procID()>0) return;
       cout << " Newton solver parameters: " << endl;
       cout << "   atol   = " << m_atol << endl;
       cout << "   rtol   = " << m_rtol << endl;
       cout << "   stol   = " << m_stol << endl;
       cout << "   maxits = " << m_maxits << endl;
       cout << "   linear_solver = " << m_linear_solver_type << endl;
    }
    
    virtual void verbose( const bool  a_verbose ) 
    { 
       m_verbose = a_verbose && m_verbose_input;
       m_linear_solver->verbose(a_verbose);
    }

    virtual void setOutputIndent( const std::string& a_indent ) 
    { 
      m_op_indent = a_indent;
      if (m_linear_solver) {
        m_linear_solver->setOutputIndent(m_op_indent+"  ");
      }
    }
     
    virtual void getConvergenceParams( int& a_l_exit_status, int& a_l_last_iter, int& a_l_total_iter, 
                                       int& a_nl_exit_status, int& a_nl_iter, int& a_nl_total_iter,
                                    double& a_nl_abs_res, double& a_nl_rel_res ) const
    {
       a_l_exit_status = m_linear_solver->getConvergenceParams( a_l_last_iter, a_l_total_iter );
       a_nl_exit_status = m_exitStatus;
       a_nl_iter = m_its;
       a_nl_total_iter = m_count;
       a_nl_abs_res = m_norm;
       a_nl_rel_res = m_norm/m_norm0;
    }
    
  protected:

  private:

    void parseParameters( ParmParse& );

    bool m_is_Defined;
    int m_maxits;
    double m_rtol, m_atol, m_stol;
    bool m_verbose, m_verbose_input, m_isLinear;
    std::string m_linear_solver_type;
    std::string m_op_indent;

    NonlinearFunction<Vec,Ops>* m_func;  
    LinearFunctionSolver<Vec,Ops>* m_linear_solver;

    mutable double m_norm, m_norm0;
    mutable Vec m_dU, m_F;
    mutable std::string m_convergedReason;
    mutable int m_its, m_exitStatus, m_count;

};

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::define( const Vec&                        a_U,
                                    Ops* const                        a_ops, 
                                    NonlinearFunction<Vec,Ops>* const a_func,
                                    const double )
{
  CH_assert(!isDefined());

  /* default values */
  m_maxits  = 10;
  m_rtol    = 1e-6;
  m_atol    = 1e-12;
  m_stol    = 1e-14;
  m_verbose_input = true;

  ParmParse ppnewton ( "newton" );
  parseParameters( ppnewton );
  m_verbose = true;

  m_count   = 0;
  m_isLinear = false;

  m_F.define(a_U); 
  m_dU.define(a_U);

  CH_assert( a_func->isDefined() );
  m_func = a_func;

  if ( m_linear_solver_type == _LINSOLVER_GMRES_) {
     m_linear_solver = new GMRESSolver<Vec,Ops>;
  } else if ( m_linear_solver_type == _LINSOLVER_PETSCKSP_) {
#ifdef with_petsc
     m_linear_solver = new PetscKSPWrapper<Vec,Ops>;
#else
     MayDay::Error("PETSc KSP solver requires compilation with PETSc");
#endif
  } else {
     MayDay::Error("NewtonSolver<Vec,Ops>::define(): Unrecognized linear solver");
  }
  m_linear_solver->define(a_ops, &(m_func->getJacobian()), a_U);
  
  m_is_Defined = true;
}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("verbose",       m_verbose_input);       /* verbosity */
  a_pp.query("atol",          m_atol);                /* absolute tolerance */
  a_pp.query("rtol",          m_rtol);                /* relative tolerance */
  a_pp.query("stol",          m_stol);                /* step tolerance */
  a_pp.query("maxits",        m_maxits);              /* maximum iterations */
  a_pp.query("linear_solver", m_linear_solver_type);  /* linear solver type */
}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::solve(  Vec&          a_U, 
                                    const Vec&    a_b,
                                    const double,
                                    const double ) const
{
  CH_TIME("NewtonSolver::solve()");
  CH_assert(isDefined());

  for (m_its = 0; m_its < m_maxits; m_its++) {

    m_func->setRHS( a_b );
    m_func->evalFunction(m_F, a_U, m_its);

    /* compute norm of m_F */
    m_norm = m_F.norm();
    if (m_its == 0) {
      if (m_norm > 0) m_norm0 = m_norm;
      else            m_norm0 = 1.0;
    }

    if (m_verbose && (!procID())) {
      printf("%s(%s) iter = %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
              m_op_indent.c_str(), _NLSOLVER_NEWTON_, m_its, m_norm, m_norm/m_norm0 );
    }

    /* check for exit conditions */
    if (m_norm < m_atol) {
      if (m_verbose && (!procID())) {
        printf("%s(%s)exiting: satisfied absolute tolerance (%1.3e).\n", 
               m_op_indent.c_str(), _NLSOLVER_NEWTON_, m_atol);
      }
      m_convergedReason = "norm less than absolute tolerance";
      m_exitStatus = 2;
      break;
    }
    if (m_norm/m_norm0 < m_rtol) {
      if (m_verbose && (!procID())) {
        printf("%s(%s)exiting: satisfied relative tolerance (%1.3e).\n", 
               m_op_indent.c_str(), _NLSOLVER_NEWTON_, m_rtol);
      }
      m_convergedReason = "norm less than relative tolerance";
      m_exitStatus = 3;
      break;
    }
    if (m_norm > 100*m_norm0) {
      if (m_verbose && (!procID())) {
        printf( "%s(%s)exiting: newton solver diverged.\n",
                m_op_indent.c_str(), _NLSOLVER_NEWTON_);
      }
      m_convergedReason = "Newton solve diverged";
      m_exitStatus = -1;
      break;
    }

    /* Compute step size [Jac]dY = F */
    m_dU.zero();
    m_linear_solver->solve(m_dU, m_F);
    double step_norm = m_dU.norm();

    if (step_norm < m_stol) {
      if (m_verbose && (!procID())) {
        printf( "%s(%s)exiting: step norm less than tolerance (%1.3e).\n",
                m_op_indent.c_str(), _NLSOLVER_NEWTON_, m_stol );
      }
      m_convergedReason = "step size norm less than tolerance";
      m_exitStatus = 4;
      break;
    }

    /* Update solution Y = Y - dY */
    a_U -= m_dU;
  }
    
  if (m_its == m_maxits) {
    if (m_verbose && (!procID())) {
      printf("%s(%s)exiting: iterations exceed max iterations (%d).\n", 
             m_op_indent.c_str(), _NLSOLVER_NEWTON_, m_maxits);
    }
    m_convergedReason = "number of iterations exceeded maximum";
    m_exitStatus = 5;
  }
      
  if((!procID()) && (m_exitStatus<0)) {
    cout << "Warning: Newton solver diverged" << endl;
    cout << "iter = " << m_its << endl;
    cout << "rel_norm = " << m_norm/m_norm0 << endl;
  }

  if((!procID()) && (m_exitStatus==5)) {
    if(m_rtol>0.0) {
      cout << "Warning: Newton solver max iterations reached" << endl;
      cout << "iter = " << m_its << endl;
      cout << "rel_norm = " << m_norm/m_norm0 << endl;
    }
  }

  m_count += m_its;
  return;
}

#include "NamespaceFooter.H"
#endif
