#ifndef _PIC_TIME_INTEGRATOR_EM_SEMI_IMPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_SEMI_IMPLICIT_H_

#include "PICTimeIntegrator.H"
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
      m_fields = nullptr;
      m_units = nullptr;

      m_rtol = 1e-6;
      m_atol = 1e-12;
      m_iter_max = 0;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_SemiImplicit() {}

    virtual 
    void define( System* const                 a_sys,
                 const PicSpeciesPtrVect&      a_particles,
                 ElectroMagneticFields* const  a_fields,
                 const CodeUnits* const        a_units )
    {
      CH_assert(!isDefined());

      m_system = a_sys;
      m_particles = a_particles;
      m_fields = a_fields;
      m_units = a_units;

      m_is_defined = true;
    }

    inline virtual 
    bool isDefined() const { return m_is_defined; }

    inline virtual 
    PICTIType type() const { return PIC_EM_SEMI_IMPLICIT; }

    virtual void preTimeStep( const Real, const Real, const int );
    virtual void timeStep( const Real, const Real, const int );

    virtual void setSolverParams( const Real  a_rtol, 
                                  const Real  a_atol, 
                                  const int   a_iter_max ) 
    { 
      m_rtol = a_rtol;
      m_atol = a_atol;
      m_iter_max = a_iter_max;
    }

  protected:

    bool m_is_defined;

    System*                 m_system;
    PicSpeciesPtrVect       m_particles;
    ElectroMagneticFields*  m_fields;
    const CodeUnits*        m_units;

    Real m_rtol, m_atol;
    int m_iter_max;

  private:

};

#include "NamespaceFooter.H"
#endif

