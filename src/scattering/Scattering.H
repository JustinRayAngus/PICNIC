#ifndef _SCATTERING_H_
#define _SCATTERING_H_

#include "REAL.H"
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
      
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity,
                                    const LevelData<FArrayBox>&  a_energyDensity ) const
      {};
      
      virtual void setMeanFreeTime( const LevelData<FArrayBox>&  a_numberDensity1,
                                    const LevelData<FArrayBox>&  a_energyDensity1,
                                    const LevelData<FArrayBox>&  a_numberDensity2,
                                    const LevelData<FArrayBox>&  a_energyDensity2 ) const
      {};
      
      virtual void initialize( const PicSpeciesPtrVect&  a_picSpeciesPtrVect,
                               const DomainGrid&         a_mesh )
      {};
      
      virtual Real scatterDt() const { return DBL_MAX; }
      
      virtual Real getDeltaEizn() const { return 0.0; }
      virtual void zeroDeltaEizn() const {;} 
      
      virtual Real getDeltaEexc() const { return 0.0; }
      virtual void zeroDeltaEexc() const {;} 
      
      virtual void applyScattering( PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                              const DomainGrid&         a_mesh,
                              const Real                a_dt_sec ) const = 0;
      
      virtual void printParameters() const = 0;
   
   protected:
   
      bool m_verbosity;
      bool m_volume_scale;
  
};

typedef RefCountedPtr<Scattering> ScatteringPtr;
typedef Vector<ScatteringPtr> ScatteringPtrVect;

#include "NamespaceFooter.H"

#endif
