#ifndef _PIC_TIME_INTEGRATOR_EM_THETAIMPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_THETAIMPLICIT_H_

#include "PICTimeIntegrator.H"
#include "NonlinearSolverLibrary.H"
#include "ODEVector.H"
#include "EMResidualFunction.H"

#define PIC_EM_THETA_IMPLICIT "PIC_EM_THETA_IMPLICIT"
#define PICMC_FULLY_IMPLICIT "PICMC_FULLY_IMPLICIT"

#include "NamespaceHeader.H"

class PICTimeIntegrator_EM_ThetaImplicit : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_EM_ThetaImplicit() 
    {
      m_system = nullptr;
      m_particles.clear();

      m_nlsolver_type = _NLSOLVER_PICARD_;
      m_nlsolver = nullptr;
      m_func = nullptr;

      m_theta = 0.5;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_ThetaImplicit() 
    {
      delete m_nlsolver;
      if (m_func) delete m_func;
    }

    virtual 
    void define( System* const,
                 const PicSpeciesPtrVect&,
                 ElectroMagneticFields* const );

    virtual void initialize();

    inline virtual 
    bool isDefined() const { return m_is_defined; }

    inline virtual 
    PICTIType type() const { return PIC_EM_THETA_IMPLICIT; }

    inline void preTimeStep( const Real, const Real, const int )
    {
      m_Uold = m_U;
    }

    virtual void timeStep( const Real, const Real, const int );
    virtual void postTimeStep( const Real, const Real );

    virtual inline void theta( const Real a_theta ) { m_theta = a_theta; }

    virtual void setSolverParams( const Real  a_rtol, 
                                  const Real  a_atol, 
                                  const int   a_iter_max ) 
    { 
      m_nlsolver->setParams( a_rtol, a_atol, a_iter_max );
    }

  protected:

    bool m_is_defined;
    std::string m_nlsolver_type;

    System*                 m_system;
    PicSpeciesPtrVect       m_particles;

    Real m_theta;

    ODEVector<System> m_U, m_Uold;
    NonlinearFunction<ODEVector<System>, System>* m_func;
    NonlinearSolver<ODEVector<System>, System>* m_nlsolver;

  private:

};

#include "NamespaceFooter.H"
#endif

