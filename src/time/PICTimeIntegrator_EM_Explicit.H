#ifndef _PIC_TIME_INTEGRATOR_EM_EXPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_EXPLICIT_H_

#include "PICTimeIntegrator.H"
#define PIC_EM_EXPLICIT "PIC_EM_EXPLICIT"
#define PICMC_EXPLICIT "PICMC_EXPLICIT"

#include "NamespaceHeader.H"

class PICTimeIntegrator_EM_Explicit : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_EM_Explicit() 
    {
      m_system = nullptr;
      m_particles.clear();
      m_fields = nullptr;
      m_units = nullptr;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_Explicit() {}

    virtual 
    void define(  System* const                 a_sys,
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
    PICTIType type() const { return PIC_EM_EXPLICIT; }

    virtual void preTimeStep( const Real, const Real, const int );
    virtual void timeStep( const Real, const Real, const int );

  protected:

    bool m_is_defined;

    System*                 m_system;
    PicSpeciesPtrVect       m_particles;
    ElectroMagneticFields*  m_fields;
    const CodeUnits*        m_units;

  private:

};

#include "NamespaceFooter.H"
#endif

