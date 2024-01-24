#ifndef _EMResidualFunction_H_
#define _EMResidualFunction_H_

#include "NonlinearFunction.H"
#include "EMJacobianFunction.H"

#include "NamespaceHeader.H"

template <class Vec, class Ops>
class EMResidualFunction: public NonlinearFunction<Vec,Ops>
{
  public:

    EMResidualFunction<Vec,Ops>() 
    {
      m_cur_time = 0.0;
      m_ops = nullptr;
      m_update_pc = true;

      // flag to update pc after call to preRHS() 
      // at t = 0 and on restart
      m_update_pc_init = true;

      m_is_Defined = false; 
    }

    ~EMResidualFunction<Vec,Ops>() {}

    void evalFunction(  Vec&, 
                        const Vec&, 
                        const int );
    
    void define(const Vec&, 
                Ops* const,
                const double,
                bool );

    inline bool isDefined() const
    { 
      return m_is_Defined;  
    }

    inline void curTime( const double a_time )
    {
      m_cur_time = a_time;
      m_jac.curTime( a_time );
    }

    inline void curTimeStep( const double a_dt )
    {
      m_dt = a_dt;
      m_jac.curTimeStep( a_dt );
    }

    inline EMJacobianFunction<Vec,Ops>& getJacobian()
    {
      return m_jac;
    }

    inline void setRHS( const Vec& a_R )
    {
      m_R = a_R;
    }

  private:

    bool m_is_Defined;
    double m_cur_time, m_dt, m_dtfac;
    bool m_update_pc;
    bool m_update_pc_init;

    Vec m_Z, m_F, m_R;
    Ops* m_ops;

    EMJacobianFunction<Vec,Ops> m_jac;
};

template <class Vec, class Ops>
void EMResidualFunction<Vec,Ops>::define( const Vec&    a_U, 
                                          Ops* const    a_ops,
                                          const double  a_dtfac,
                                          const bool    a_update_pc )
{
  m_Z.define(a_U);
  m_F.define(a_U);
  m_R.define(a_U);

  m_R = 0.0;

  m_dtfac = a_dtfac;

  m_ops = a_ops;
  m_jac.define( a_U, a_ops, a_dtfac );

  m_update_pc = a_update_pc;
  m_is_Defined = true;
}

template <class Vec, class Ops>
void EMResidualFunction<Vec,Ops>::evalFunction( Vec&        a_F, 
                                                const Vec&  a_U, 
                                                const int   a_its )
{
  CH_TIME("EMResidualFunction::evalFunction()");
  CH_assert(isDefined());

  m_ops->preRHSOp( a_U, m_cur_time, m_dt, a_its, false );
  m_ops->computeRHS( m_F, a_U, m_cur_time, m_dtfac*m_dt );

  /* set base Y and G(Y) for Jacobian */
  m_jac.setBaseSolution(a_U);
  m_jac.setBaseRHS(m_F);
  /* update preconditioner */
  if (m_update_pc || m_update_pc_init) m_jac.updatePreCondMat(a_U);
  m_update_pc_init = false;

  /* Compute ([Y - G(Y)] - b) */
  a_F = a_U - (m_F + m_R);

  return;
}

#include "NamespaceFooter.H"
#endif
