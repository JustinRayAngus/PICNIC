#ifndef _NewtonSolver_H_
#define _NewtonSolver_H_

#include <iostream>
#include <string>

#include "NonlinearSolver.H"
#include "GMRESSolver.H"
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

    void setParams( const double, 
                    const double, 
                    const int )
    {
    }

  protected:

  private:

    void parseParameters( ParmParse& );

    bool m_is_Defined;
    int m_maxits;
    double m_rtol, m_atol, m_stol;
    bool m_verbose, m_isLinear;
    std::string m_linear_solver_type;
    double m_dtfac;

    NonlinearFunction<Vec,Ops>* m_func;  
    GMRESSolver<Vec,Ops>* m_linear_solver;

    mutable double m_norm, m_norm0;
    mutable Vec m_dU, m_F;
    mutable std::string m_convergedReason;
    mutable int m_its, m_exitStatus, m_count;

};

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::define( const Vec&                        a_U,
                                    Ops* const, 
                                    NonlinearFunction<Vec,Ops>* const a_func,
                                    const double                      a_dtfac )
{
  CH_assert(!isDefined());

  /* default values */
  m_maxits  = 10;
  m_rtol    = 1e-6;
  m_atol    = 1e-12;
  m_stol    = 1e-14;
  m_verbose = true;

  ParmParse ppnewton ( "newton" );
  parseParameters( ppnewton );

  m_dtfac = a_dtfac;

  m_count   = 0;
  m_isLinear = false;
  m_linear_solver_type = _LINSOLVER_GMRES_;

  m_F.define(a_U); 
  m_dU.define(a_U);

  CH_assert( a_func->isDefined() );
  m_func = a_func;

  if ( m_linear_solver_type == _LINSOLVER_GMRES_) {
     m_linear_solver = new GMRESSolver<Vec,Ops>;
  } else {
     MayDay::Error("NewtonSolver<Vec,Ops>::define(): Unrecognized linear solver");
  }
  m_linear_solver->define(&(m_func->getJacobian()), a_U);
  
  m_is_Defined = true;
}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("verbose", m_verbose); /* verbosity */
  a_pp.query("atol",    m_atol);    /* absolute tolerance */
  a_pp.query("rtol",    m_rtol);    /* relative tolerance */
  a_pp.query("maxits",  m_maxits);  /* maximum iterations */
}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::solve(  Vec&          a_U, 
                                    const Vec&    a_b,
                                    const double,
                                    const double ) const
{
  CH_TIME("NewtonSolver::solve()");
  CH_assert(isDefined());

  m_its = 0;
  while(1) {

    m_func->evalFunction(m_F, a_U, a_b, m_its);

    /* compute norm of m_F */
    m_norm = m_F.norm();
    if (m_its == 0) {
      if (m_norm > 0) m_norm0 = m_norm;
      else            m_norm0 = 1.0;
    }

    if (m_verbose && (!procID())) {
      printf("  (%s) iter = %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
              _NLSOLVER_NEWTON_, m_its, m_norm, m_norm/m_norm0 );
    }

    /* check for exit conditions */
    if (m_norm < m_atol) {
      if (!procID()) {
        printf("  exiting: satisfied absolute tolerance (%1.3e).\n", 
               m_atol);
      }
      m_convergedReason = "norm less than absolute tolerance";
      m_exitStatus = 0;
      break;
    }
    if (m_norm/m_norm0 < m_rtol) {
      if (!procID()) {
        printf("  exiting: satisfied relative tolerance (%1.3e).\n", 
               m_rtol);
      }
      m_convergedReason = "norm less than relative tolerance";
      m_exitStatus = 0;
      break;
    }
    if (m_its > m_maxits) {
      if (!procID()) {
        printf("  exiting: iterations exceed max iterations (%d).\n", 
               m_maxits);
      }
      m_convergedReason = "number of iterations exceeded maximum";
      m_exitStatus = 1;
      break;
    }
    if (m_norm > 100*m_norm0) {
      if (!procID()) {
        printf("  exiting: newton solver diverged.\n");
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
      if (!procID()) {
        printf( "  exiting: step norm less than tolerance (%1.3e).\n",
                m_stol );
      }
      m_convergedReason = "step size norm less than tolerance";
      m_exitStatus = 0;
      break;
    }

    /* Update solution Y = Y - dY */
    a_U -= m_dU;
    m_its++;
  }

  m_count += m_its;
  return;
}

#include "NamespaceFooter.H"
#endif
