#ifndef _GaussSeidelSolver_H_
#define _GaussSeidelSolver_H_

#include <string>
#include "ParmParse.H"
#include "parstream.H"
#include "LinearFunctionSolver.H"
#include "BandedMatrixGaussSeidelSolver.H"

#define _LINSOLVER_GAUSSSEIDEL_ "gauss_seidel"

#include "NamespaceHeader.H"

template <class T, class Ops>
class GaussSeidelSolver : public LinearFunctionSolver<T, Ops>
{
  public:

    GaussSeidelSolver<T,Ops>() 
    {
      m_is_defined = false; 
      m_op = nullptr; 
      m_solver = nullptr;
    }

    ~GaussSeidelSolver<T,Ops>()  
    {
      delete m_solver;
    }

    void define( Ops* const, LinearFunction<T,Ops>* const, const T& );
    void solve( T&, const T&, int a_its = -1 );

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
      if (m_solver) {
        m_solver->setOutputIndent(m_op_indent);
      }
    }

  private:
    
    bool m_is_defined;
    bool m_verbose;
    std::string m_op_indent;

    LinearFunction<T,Ops>* m_op;
    BandedMatrixSolver<T>* m_solver;
};

template <class T, class Ops>
void GaussSeidelSolver<T,Ops>::define(Ops* const,
                                      LinearFunction<T,Ops>* const  a_op, 
                                      const T&                      a_x )
{
  /* default values */
  m_op = a_op;
  CH_assert(m_op != NULL);

  m_solver = new BandedMatrixGaussSeidelSolver<T>;
  m_solver->define(a_x);

  m_is_defined = true;
}

template <class T, class Ops>
void GaussSeidelSolver<T,Ops>::solve(T& a_x, const T& a_b, int a_its)
{
  CH_assert(m_is_defined);
  m_solver->solve(  a_x,
                    m_op->getBandedMatrixForm(),
                    a_b,
                    a_its );
  return;
}

#include "NamespaceFooter.H"
#endif
