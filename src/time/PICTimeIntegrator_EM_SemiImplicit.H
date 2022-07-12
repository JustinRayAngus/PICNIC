#ifndef _PIC_TIME_INTEGRATOR_EM_SEMI_IMPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_SEMI_IMPLICIT_H_

#include "PICTimeIntegrator.H"
#include "NonlinearSolverLibrary.H"
#include "ODEVector.H"

#define PIC_EM_SEMI_IMPLICIT "PIC_EM_SEMI_IMPLICIT"
#define PICMC_SEMI_IMPLICIT "PICMC_SEMI_IMPLICIT"

#include "NamespaceHeader.H"

class PICTimeIntegrator_EM_SemiImplicit : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_EM_SemiImplicit() 
    {
      m_system = nullptr;
      m_particles.clear();
      m_nlsolver = nullptr;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_SemiImplicit() 
    {
      delete m_nlsolver;
    }

    virtual 
    void define( System* const,
                 const PicSpeciesPtrVect&,
                 ElectroMagneticFields* const );

    virtual void initialize();
    
    virtual int prepForCheckpoint() const;

    inline virtual 
    bool isDefined() const { return m_is_defined; }

    inline virtual 
    PICTIType type() const { return PIC_EM_SEMI_IMPLICIT; }

    virtual void preTimeStep( const Real, const Real, const int );
    virtual void timeStep( const Real, const Real, const int );
    virtual void postTimeStep( const Real, const Real );

    virtual void setSolverParams( const Real  a_rtol, 
                                  const Real  a_atol, 
                                  const int   a_iter_max ) 
    { 
      m_nlsolver->setParams( a_rtol, a_atol, a_iter_max );
    }

  protected:

    bool m_is_defined;

    System*                 m_system;
    PicSpeciesPtrVect       m_particles;

    ODEVector<System> m_E, m_Eold;
    ODEVector<System> m_B, m_Bold, m_FB;
    NonlinearSolver<ODEVector<System>, System>* m_nlsolver;

  private:

};

#include "NamespaceFooter.H"
#endif

