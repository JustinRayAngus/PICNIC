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

    PICTimeIntegrator_EM_Explicit() {}

    virtual ~PICTimeIntegrator_EM_Explicit() {}

    virtual 
    void define( System* const  a_system,
                 PicSpeciesInterface* const  a_pic_species,
                 EMFields* const             a_fields )
    {
      CH_assert(!isDefined());

      m_system = a_system;
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
      pp.query("strang_splitting", m_strang_splitting);
      if(!m_strang_splitting) { pp.query("average_v_deposit", m_average_v_deposit); }
      if(m_strang_splitting || m_average_v_deposit) { m_post_push_scatter = true; }
      else { pp.query("post_push_scatter", m_post_push_scatter); }
      pp.query("init_half_advanceX", m_init_half_advanceX);
      pp.query("init_half_advanceE", m_init_half_advanceE);

      if (m_init_half_advanceX && procID()==0) {
         cout << "Warning!!! explicit em pic with init_half_advanceX = true can cause " << endl;
         cout << "           violation in charge conservation when using charge-conserving " << endl;
         cout << "           current deposition schemes." << endl;
      }

      m_is_defined = true;
      printParams();
      
    }

    inline
    virtual void printParams() const
    {
      if(procID()>0) return;
      cout << "================== Time Solver ==================" << endl;
      cout << "advance method = PIC_EM_EXPLICIT" << endl;
      cout << "init_half_advanceX = " << (m_init_half_advanceX?"true":"false") << endl;
      cout << "init_half_advanceE = " << (m_init_half_advanceE?"true":"false") << endl;
      cout << "strang_splitting   = " << (m_strang_splitting?"true":"false") << endl;
      cout << "average_v_deposit  = " << (m_average_v_deposit?"true":"false") << endl;
      cout << "post_push_scatter  = " << (m_post_push_scatter?"true":"false") << endl;
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
   // virtual void postTimeStep( const Real, const Real );
    
  protected:

    bool m_is_defined = false;
    bool m_init_half_advanceX = false;
    bool m_init_half_advanceE = false;
    bool m_strang_splitting = false;
    bool m_average_v_deposit = false;
    bool m_post_push_scatter = false;

    System* m_system = nullptr;
    PicSpeciesInterface* m_pic_species = nullptr;
    EMFields* m_fields = nullptr;

    ODEVector<EMFields> m_E, m_Eold, m_FE;
    ODEVector<EMFields> m_B, m_Bold, m_FB;

  private:

};

#include "NamespaceFooter.H"
#endif

