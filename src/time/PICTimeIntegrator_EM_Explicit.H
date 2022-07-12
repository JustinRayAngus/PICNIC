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
      m_system = nullptr;
      m_particles.clear();
      m_fields = nullptr;
      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_Explicit() {}

    virtual 
    void define(  System* const                 a_sys,
                  const PicSpeciesPtrVect&      a_particles,
                  ElectroMagneticFields* const  a_fields )
    {
      CH_assert(!isDefined());

      m_system = a_sys;
      m_particles = a_particles;
      m_fields = a_fields;

      m_E.define(m_fields->getVectorSize(e_only));
      m_Eold.define(m_E);
      m_FE.define(m_E);

      m_B.define(m_fields->getVectorSize(b_only));
      m_Bold.define(m_B);
      m_FB.define(m_B);
    
      // parse options for this integrator
      ParmParse pp("pic_em_explicit");

      // old ordering is not compatible with Boris method for inertia
      m_old_ordering = true;
      pp.query("old_ordering", m_old_ordering);

      if(!procID()) {
         cout << "advance method = PIC_EM_EXPLICIT" << endl;
         cout << "old ordering = " << m_old_ordering << endl;
         cout << endl;
      }

      m_is_defined = true;
      
      if(m_old_ordering) {
        for (int s=0; s<m_particles.size(); s++) {
          auto this_picSpecies(m_particles[s]);
          bool boris_inertia = this_picSpecies->borisInertia();
          if(boris_inertia) { 
            cout << "old ordering is incompatible with Boris method for inertia !!!"  << endl;
          }
          CH_assert(!boris_inertia);
        }
      }
    
    }

    virtual void initialize();
    
    virtual int prepForCheckpoint() const;

    inline virtual 
    bool isDefined() const { return m_is_defined; }

    inline virtual 
    PICTIType type() const { return PIC_EM_EXPLICIT; }

    virtual void preTimeStep( const Real, const Real, const int );
    virtual void timeStep( const Real, const Real, const int );

  protected:

    bool m_is_defined;
    bool m_old_ordering;

    System*                 m_system;
    PicSpeciesPtrVect       m_particles;
    ElectroMagneticFields*  m_fields;

    ODEVector<System> m_E, m_Eold, m_FE;
    ODEVector<System> m_B, m_Bold, m_FB;

  private:

};

#include "NamespaceFooter.H"
#endif

