#ifndef _PIC_TIME_INTEGRATOR_DSMC_H_
#define _PIC_TIME_INTEGRATOR_DSMC_H_

#include "PICTimeIntegrator.H"

#define PIC_DSMC "PIC_DSMC"

#include "NamespaceHeader.H"

class PICTimeIntegrator_DSMC : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_DSMC() 
    {
      m_system = nullptr;
      m_particles.clear();
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_DSMC() {}

    virtual 
    void define(  System* const                 a_sys,
                  const PicSpeciesPtrVect&      a_particles,
                  ElectroMagneticFields* const  a_fields )
    {
      CH_assert(!isDefined());

      m_system = a_sys;
      m_particles = a_particles;

      printParams();

      m_is_defined = true;
    }

    inline virtual 
    bool isDefined() const { return m_is_defined; }
    
    virtual void printParams() const 
    {
       if(procID()>0) return;
       cout << "================== Time Solver ==================" << endl;
       cout << "advance method = PIC_DSMC" << endl;
       cout << "=================================================" << endl;
       cout << endl;
    }

    inline virtual 
    PICTIType type() const { return PIC_DSMC; }

    virtual void preTimeStep( const Real, const Real, const int );
    virtual void timeStep( const Real, const Real, const int );

  protected:

    bool m_is_defined;

    System*                 m_system;
    PicSpeciesPtrVect       m_particles;

  private:

};

#include "NamespaceFooter.H"
#endif

