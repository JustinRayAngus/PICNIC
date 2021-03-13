#ifndef _SCATTERING_H_
#define _SCATTERING_H_

#include "REAL.H"
#include "PicSpecies.H"
#include "LevelData.H"

#include "Vector.H"
#include "RefCountedPtr.H"

#include "DomainGrid.H"

#include "NamespaceHeader.H"

/**
 * Generic scattering interface.
*/
class Scattering
{
   public:

      // Constructor
      //
      //Scattering( const PicSpecies&  a_picSpecies, 
      Scattering( const int  a_verbosity );
   
      // Destructor
      //
      virtual ~Scattering() {;}
  
      virtual int species1() const = 0;
     
      virtual int species2() const = 0;
      
      // set mean free time
      //
      virtual void setMeanFreeTime( const DomainGrid&            a_mesh,
                                    const LevelData<FArrayBox>&  a_numberDensity,
                                    const LevelData<FArrayBox>&  a_energyDensity )
      {};
      
      // initialization
      //
      virtual void initialize( const DomainGrid&         a_mesh,
                               const PicSpeciesPtrVect&  a_picSpeciesPtrVect )
      {};
      
      virtual Real scatterDt() const = 0;

      // apply self scattering to passed pic species
      //
      virtual void applySelfScattering( PicSpecies&            a_picSpecies,
                                  const DomainGrid&            a_mesh,
                                  const LevelData<FArrayBox>&  a_numberDensity,
                                  const LevelData<FArrayBox>&  a_energyDensity,
                                  const Real                   a_dt ) const = 0;

      // Print object parameters.
      //
      virtual void printParameters() const = 0;
   
   protected:
   
      bool m_verbosity;
  
};

typedef RefCountedPtr<Scattering> ScatteringPtr;
typedef Vector<ScatteringPtr> ScatteringPtrVect;

#include "NamespaceFooter.H"

#endif
