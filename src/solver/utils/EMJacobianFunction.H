#ifndef _EMJacobianFunction_H_
#define _EMJacobianFunction_H_

#include "PreconditionerLibrary.H"
#include "LinearFunction.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class EMJacobianFunction: public LinearFunction<T,Ops>
{
  public:

    EMJacobianFunction<T,Ops>()
    {
      m_is_Defined = false; 
      m_is_linear = false;
      m_epsJFNK = 1e-6;
      m_usePreCond = true;
      m_pc_type = _EM_MATRIX_PC_;
    }

    ~EMJacobianFunction<T,Ops>() 
    {
      delete m_preCond;
      return;
    }

    void applyOp(T& a_F, const T& a_U, int a_its = -1);

    inline void residual(T& a_R, const T& a_U, const T& a_b)
    {
      applyOp(a_R, a_U);
      a_R -= a_b;
    }

    inline void axby(T& a_Z, const T& a_X, const T& a_Y, double a_a, double a_b)
    {
      a_Z = a_a*a_X + a_b*a_Y;
    }

    inline void preCond(T& a_U, const T& a_X)     
    { 
      if (m_usePreCond) {
        a_U.zero();
        m_preCond->apply(a_U, a_X);
      } else {
        a_U.copy(a_X);
      }
      return;
    }

    inline void updatePreCondMat( const T&  a_X )
    {
      CH_TIME("EMJacobianFunction::updatePreCondMat()");
      if (m_usePreCond) { m_preCond->update( a_X ); }
      return;
    }

#ifdef with_petsc
    inline int copyPCMatrixToPetsc( Mat a_P )
    {
      if (m_usePreCond) {
        m_preCond->copyMatrixToPetsc( a_P );
        return 0;
      } else {
        return 1;
      }
    }
#endif

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

    inline bool isDefined() const { return m_is_Defined;  }

    inline void setBaseSolution (const T& a_U)
    { 
      m_Y0.copy(a_U); 
      m_normY0 = m_Y0.norm(); 
      if (m_usePreCond) m_preCond->setBaseSolution(a_U);
    }

    inline void setBaseRHS (const T& a_F)
    { 
      m_F0.copy(a_F); 
      if (m_usePreCond) m_preCond->setBaseRHS(a_F);
    }

    inline void setJFNKEps (double a_eps)
    { 
      m_epsJFNK = a_eps; 
      if (m_usePreCond) m_preCond->setJFNKEps(a_eps);
    }

    inline void setIsLinear (bool a_isLinear) 
    { 
      m_is_linear = a_isLinear; 
      if (m_usePreCond) m_preCond->setIsLinear(true);
    }

    inline void curTime( double a_time )
    {
      m_cur_time = a_time;
      if (m_usePreCond) m_preCond->curTime( a_time );
    }

    inline void curTimeStep( double a_dt )
    {
      m_dt = a_dt;
      if (m_usePreCond) m_preCond->curTimeStep( a_dt );
    }

    void define(const T&, Ops* const, const double);

    void setOutputIndent( const std::string& a_indent )
    {
      if (m_preCond) {
        m_preCond->setOutputIndent( a_indent+"  " );
      }
    }

  private:
  
    bool m_is_Defined, m_is_linear, m_usePreCond;
    double m_epsJFNK, m_normY0;
    double m_cur_time, m_dt, m_dtfac;
    std::string m_pc_type;

    T m_Z, m_Y0, m_F0, m_F;
    Ops* m_ops;
    Preconditioner<T,Ops>* m_preCond;
    
    void parseParameters( ParmParse& a_pp )
    {
      a_pp.query("jfnk_eps", m_epsJFNK);
      a_pp.query("with_pc",  m_usePreCond);
      if (m_usePreCond) {
        a_pp.query("pc_type",  m_pc_type);
      }
    }
};

template <class T, class Ops>
void EMJacobianFunction<T,Ops>::define( const T&     a_U, 
                                        Ops* const   a_ops,
                                        const double a_dtfac )
{
  m_Z.define  (a_U);
  m_Y0.define (a_U);
  m_F0.define (a_U);
  m_F.define  (a_U);

  m_ops = a_ops;
  m_dtfac = a_dtfac;

  ParmParse ppIJac( "jacobian" );
  parseParameters( ppIJac );

  m_preCond = nullptr;
  if (m_usePreCond) {
    if (m_pc_type == _EM_MATRIX_PC_) {
      m_preCond = new EMMatrixPreconditioner<T,Ops>;
    } else {
      MayDay::Error("Invalid choice of preconditioner type");
    }
    if (!procID()) {
      printf("EMJacobianFunction: preconditioner type is %s.\n", m_pc_type.c_str());
    }
    m_preCond->define(a_U, a_ops, a_dtfac);
    m_preCond->setJFNKEps( m_epsJFNK );
  }

  m_is_Defined = true;
  return;
}

template <class T, class Ops>
void EMJacobianFunction<T,Ops>::applyOp(T& a_F, const T& a_U, int a_its)
{
  CH_TIME("EMJacobianFunction::applyOp()");
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
