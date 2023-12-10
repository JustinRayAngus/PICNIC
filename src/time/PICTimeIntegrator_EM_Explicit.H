#ifndef _PIC_TIME_INTEGRATOR_EM_EXPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_EXPLICIT_H_

#include "PICTimeIntegrator.H"
#include "ODEVector.H"

#define PIC_EM_EXPLICIT "PIC_EM_EXPLICIT"
#define PICMC_EXPLICIT "PICMC_EXPLICIT"

#include "NamespaceHeader.H"

class PICTimeIntegrator_EM_Explicit : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_EM_Explicit() 
    {
      m_pic_species = nullptr;
      m_fields = nullptr;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_Explicit() {}

    virtual 
    void define( PicSpeciesInterface* const  a_pic_species,
                 EMFields* const             a_fields )
    {
      CH_assert(!isDefined());

      m_pic_species = a_pic_species;
      m_fields = a_fields;

      m_E.define(m_fields->getVectorSize(e_only));
      m_Eold.define(m_E);
      m_FE.define(m_E);

      m_B.define(m_fields->getVectorSize(b_only));
      m_Bold.define(m_B);
      m_FB.define(m_B);
    
      // parse options for this integrator
      ParmParse pp("pic_em_explicit");

      m_is_defined = true;
      printParams();
      
    }

    inline
    virtual void printParams() const
    {
      if(procID()>0) return;
      cout << "================== Time Solver ==================" << endl;
      cout << "advance method = PIC_EM_EXPLICIT" << endl;
      cout << "=================================================" << endl;
      cout << endl;
    }

    virtual void initialize( const std::string&  a_restart_file_name );
    
    virtual int prepForCheckpoint() const;

    inline virtual 
    bool isDefined() const { return m_is_defined; }
    
    inline virtual 
    PICTIType type() const { return PIC_EM_EXPLICIT; }

    virtual void preTimeStep( const Real, const Real, const int );
    virtual void timeStep( const Real, const Real, const int );
    
  protected:

    bool m_is_defined;

    PicSpeciesInterface*    m_pic_species;
    EMFields*  m_fields;

    ODEVector<EMFields> m_E, m_Eold, m_FE;
    ODEVector<EMFields> m_B, m_Bold, m_FB;

  private:

};

#include "NamespaceFooter.H"
#endif

