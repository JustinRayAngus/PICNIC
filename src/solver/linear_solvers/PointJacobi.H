#ifndef _PointJacobiSolver_H_
#define _PointJacobiSolver_H_

#include <string>
#include "LinearFunctionSolver.H"
#include "ParmParse.H"
#include "parstream.H"

#define _LINSOLVER_JACOBI_ "jacobi"

#include "NamespaceHeader.H"

template <class T, class Ops>
class PointJacobiSolver : public LinearFunctionSolver<T, Ops>
{
  public:

    PointJacobiSolver<T,Ops>()
    { 
      m_is_defined = false; 
      m_op = nullptr; 
      m_verbose = true;
      m_omega = 0.9;
      m_its = 200;
      m_atol = 1e-10;
      m_rtol = 1e-10;
    }

    ~PointJacobiSolver<T,Ops>()  {}

    void define (Ops* const, LinearFunction<T,Ops>* const, const T&);
    void solve(T&, const T&, int a_its = -1);

    inline void verbose(const bool a_v) 
    {
      m_verbose = a_v;
    }

    virtual void setOutputIndent( const std::string& a_indent ) 
    { 
      m_op_indent = a_indent;
      if (m_op) {
        m_op->setOutputIndent(m_op_indent);
      }
    }

  private:
    
    bool m_is_defined;
    bool m_verbose;
    std::string m_op_indent;
    int m_its; /*!< Number of Jacobi iterations */

    double m_omega;

    T  m_y;     /*!< Work vector */

    double m_atol, m_rtol;

    LinearFunction<T,Ops>* m_op;

    void parseParameters(ParmParse& a_pp)
    {
      a_pp.query("its", m_its);
      a_pp.query("omega", m_omega);
      a_pp.query("verbose", m_verbose);
      a_pp.query("atol", m_atol);
      a_pp.query("rtol", m_rtol);
      return;
    }
};

template <class T, class Ops>
void PointJacobiSolver<T,Ops>::define(  Ops* const,
                                        LinearFunction<T,Ops>* const  a_op, 
                                        const T&                      a_x )
{
  m_op = a_op;
  CH_assert(m_op != NULL);

  ParmParse pp_Jac( _LINSOLVER_JACOBI_ );
  parseParameters( pp_Jac );

  m_y.define(a_x);
  m_is_defined = true;
}

template <class T, class Ops>
void PointJacobiSolver<T,Ops>::solve(T& a_x, const T& a_b, int a_its)
{
  CH_assert(m_is_defined);
  if (a_its < 0) a_its = m_its;
  double norm, norm0;

  std::vector<double> diag;
  m_op->getDiagonal(diag);

  int n(0);
  for (n = 0; n < a_its; n++) {

    m_op->residual(m_y, a_x, a_b);
    norm = m_y.norm();

    if (n == 0) norm0 = norm;
    if (m_verbose && (!procID())) {
      printf("%s(%s) iter %4d, residual = %1.6e (abs.), %1.6e (rel.)\n",
              m_op_indent.c_str(), _LINSOLVER_JACOBI_, 
              n, norm, norm/norm0 );
    }

    if ((norm < m_atol) || (norm/norm0 < m_rtol)) break;
    if (norm/norm0 > 1) {
      if (m_verbose && (!procID())) {
        printf( "%s(%s) solver diverged.\n", 
                m_op_indent.c_str(), _LINSOLVER_JACOBI_ );
      }
      break;
    }

    for (int i(0); i < diag.size(); i++) {
      a_x[i] -= m_omega * m_y[i] / diag[i];
    }

  }

  return;
}

#include "NamespaceFooter.H"
#endif
