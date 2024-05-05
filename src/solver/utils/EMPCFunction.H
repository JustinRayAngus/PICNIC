#ifndef _EMPCFUNCTION_H_
#define _EMPCFUNCTION_H_

#include <string>

#include "PCFunction.H"
#include "BandedMatrix.H"
#include "BandedMatrixSolver.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class EMPCFunction: public PCFunction<T,Ops>
{
  public:

    EMPCFunction<T,Ops>()
    {
      m_is_Defined = false;
      m_is_linear = false;
      m_matrix_free = false;
      m_ops = nullptr;
      m_epsJFNK = 1e-6;
      m_op_indent = "";
    }

    ~EMPCFunction<T,Ops>() {}

    inline void applyOp(T& a_F, const T& a_U, int a_its = -1)
    {
      if (m_matrix_free) {
        applyRHSOp( a_F, a_U, a_its );
      } else {
        applyMatOp( a_F, a_U, a_its );
      }
    }

    inline void residual(T& a_R, const T& a_U, const T& a_b)
    {
      applyOp(a_R, a_U);
      a_R -= a_b;
    }

    inline void axby(T& a_Z, const T& a_X, const T& a_Y, double a_a, double a_b)
    {
      a_Z = a_a*a_X + a_b*a_Y;
    }

    inline void preCond(T& a_Y, const T& a_X) { a_Y.copy(a_X); }

    inline void create(T& a_Z, const T& a_U)
    {
      a_Z.define(a_U);
    }

    inline void assign(T& a_Z, const T& a_U) { a_Z = a_U; }

    inline void incr(T& a_Z, const T& a_U, double a_scale)
    {
      a_Z.increment(a_U,a_scale);
    }

    inline void scale(T& a_U, const double& a_scale)
    {
      a_U.scale(a_scale);
    }

    inline void setToZero (T& a_U)
    {
      a_U.zero();
    }

    inline double dotProduct(T& a_X, const T& a_U)
    {
      return(a_X.dotProduct(a_U));
    }

    inline double norm(T& a_U, int a_ord)
    {
      CH_assert(a_ord == 2);
      return(a_U.norm());
    }

    inline bool isDefined() const
    {
      return (m_is_Defined && m_Pmat.isDefined());
    }

    inline void setBaseSolution (const T& a_U)
    {
      m_Y0.copy(a_U);
      m_normY0 = m_Y0.norm();
    }

    inline void setBaseRHS (const T& a_F)
    {
      m_F0.copy(a_F);
    }

    inline void setJFNKEps (double a_eps)
    {
      m_epsJFNK = a_eps;
    }

    inline void setIsLinear (bool a_isLinear)
    {
      m_is_linear = a_isLinear;
    }

    inline void curTime( double a_time )
    {
      m_cur_time = a_time;
    }

    inline void curTimeStep( double a_dt )
    {
      m_dt = a_dt;
    }

    inline void getDiagonal(std::vector<double>& a_D) const
    {
      return m_Pmat.getDiagonal(a_D);
    }

    void define(const T&, Ops* const, const double);

    inline const BandedMatrix& getBandedMatrixForm() const { return m_Pmat; }
    inline BandedMatrix& getBandedMatrixForm() { return m_Pmat; }

    inline void update( const T& a_U )
    {
      CH_TIME("EMPCFunction::update()");
      m_ops->updatePrecondMat( m_Pmat, m_cur_time, m_dt*m_dtfac );
    }

    inline void setOutputIndent( const std::string& a_indent )
    {
      m_op_indent = a_indent;
      m_Pmat.setOutputIndent( a_indent );
    }

    inline void setMatrixFree( const bool a_bool )
    {
      m_matrix_free = true;
    }

  private:

    bool m_is_Defined, m_is_linear;
    bool m_matrix_free;
    double m_cur_time, m_dt, m_dtfac;
    double m_epsJFNK, m_normY0;

    std::string m_op_indent;

    T m_Z, m_Y0, m_F0, m_F;

    Ops* m_ops;
    BandedMatrix m_Pmat;
    BandedMatrixSolver<T> m_dummy_solver;

    void applyRHSOp( T&, const T&, int );
    void applyMatOp( T&, const T&, int );

};

template <class T, class Ops>
void EMPCFunction<T,Ops>::define( const T&      a_U,
                                  Ops* const    a_ops,
                                  const double  a_dtfac )
{
  m_Z.define  (a_U);
  m_Y0.define (a_U);
  m_F0.define (a_U);
  m_F.define  (a_U);

  m_ops = a_ops;
  m_dtfac = a_dtfac;

  m_Pmat.define( a_U.getVectorSize(), m_ops->numPCMatBands(), a_U.mpiOffset() );
  m_Pmat.setToIdentityMatrix();
  m_Pmat.finalAssembly();

  CH_assert(m_ops != NULL);
  m_is_Defined = true;
  return;
}

template <class T, class Ops>
void EMPCFunction<T,Ops>::applyMatOp(T& a_F, const T& a_U, int)
{
  CH_TIME("EMPCFunction::applyMatOp()");
  m_dummy_solver.applyOp( a_F, a_U, m_Pmat );
}

template <class T, class Ops>
void EMPCFunction<T,Ops>::applyRHSOp(T& a_F, const T& a_U, int)
{
  CH_TIME("EMPCFunction::applyRHSOp()");
  double normY = a_U.norm();

  if (normY < 1e-15) {

    a_F = 0;

  } else {

    double eps;
    if (m_is_linear) {
      eps = 1.0;
    } else {
      /* eps = error_rel * sqrt(1 + ||U||) / ||a||
       * M. Pernice and H. F. Walker, "NITSOL: A Newton Iterative Solver for
       * Nonlinear Systems", SIAM J. Sci. Stat. Comput.", 1998, vol 19,
       * pp. 302--318. */
      eps = m_epsJFNK * sqrt(1 + m_normY0) / normY;
    }

    m_Z = m_Y0 + eps*a_U;
    m_ops->preRHSOp( m_Z, m_cur_time, m_dt, -1, true );
    m_ops->computeRHS(m_F, m_Z, m_cur_time, m_dtfac*m_dt );

    a_F = a_U - (1.0/eps)*(m_F-m_F0);
  }

}

#include "NamespaceFooter.H"
#endif

