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
      m_pic_species = nullptr;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_DSMC() {}

    virtual 
    void define( PicSpeciesInterface* const    a_pic_species,
                 EMFields* const  a_fields )
    {
      CH_assert(!isDefined());

      m_pic_species = a_pic_species;

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

    PicSpeciesInterface*    m_pic_species;

  private:

};

#include "NamespaceFooter.H"
#endif

