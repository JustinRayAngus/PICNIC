#ifndef _SCATTERING_H_
#define _SCATTERING_H_

#include "REAL.H"
#include "PicSpeciesInterface.H"
#include "PicSpecies.H"
#include "LevelData.H"

#include "Vector.H"
#include "RefCountedPtr.H"

#include "DomainGrid.H"
#include "MathUtils.H"

#include "NamespaceHeader.H"

/**
 * Generic scattering interface.
*/

typedef enum { ELASTIC, 
               EXCITATION, 
               IONIZATION, 
               MONTE_CARLO_NULL, 
               COULOMB, 
               FUSION,
               CHARGE_EXCHANGE } ScatteringType;

typedef enum { PROBABILISTIC,
               CONSERVATIVE } WeightMethod;
      
class Scattering
{
   public:

      //Scattering( const PicSpecies&  a_picSpecies, 
      Scattering( const int  a_verbosity );
   
      virtual ~Scattering() {;}
  
      virtual ScatteringType getScatteringType() const = 0;
     
      virtual int species1() const = 0;
      
      virtual int species2() const = 0;
      
      virtual void setMeanFreeTime( const PicSpeciesInterface&  a_pic_species_intf ) const = 0;
      
      virtual void initialize( const PicSpeciesInterface&  a_pic_species_intf,
                               const DomainGrid&           a_mesh ) = 0;
      
      virtual Real scatterDt() const { return DBL_MAX; }
      
      virtual Real getDeltaEizn() const { return 0.0; }
      virtual void zeroDeltaEizn() const {;} 
      
      virtual Real getDeltaEexc() const { return 0.0; }
      virtual void zeroDeltaEexc() const {;} 
      
      virtual Real getDeltaEfusion() const { return 0.0; }
      virtual void zeroDeltaEfusion() const {;} 
      
      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const = 0;
      
      virtual void printParameters() const = 0;
   
   protected:
   
      bool m_verbosity;
  
};

typedef RefCountedPtr<Scattering> ScatteringPtr;
typedef Vector<ScatteringPtr> ScatteringPtrVect;

#include "NamespaceFooter.H"

#endif
