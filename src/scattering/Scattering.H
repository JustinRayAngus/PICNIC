#ifndef _SCATTERING_H_
#define _SCATTERING_H_

#include "REAL.H"
#include "PicSpeciesInterface.H"
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
               BREMSSTRAHLUNG,
               IBREMSSTRAHLUNG,
               CHARGE_EXCHANGE } ScatteringType;

typedef enum { PROBABILISTIC,
               CONSERVATIVE } WeightMethod;
      
class Scattering
{
   public:

      Scattering( const int  a_verbosity );
   
      virtual ~Scattering() {;}
  
      virtual ScatteringType getScatteringType() const = 0;

      virtual std::string getScatteringSubTypeName() const { return "none"; }
     
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
      
      virtual Real getFusionQ() const { return 0.0; }
      virtual Real getFusionQb() const { return 0.0; }
      virtual Real getDeltaEfusion() const { return 0.0; }
      virtual void zeroDeltaEfusion() const {;} 
      
      virtual Real getDeltaEbremsstrahlung() const { return 0.0; }
      virtual void zeroDeltaEbremsstrahlung() const {;} 

      virtual Real getDeltaEIBremsstrahlung() const { return 0.0; }
      virtual void zeroDeltaEIBremsstrahlung() const {;}

      LevelData<FArrayBox>& getFusionProducts() { return m_fusionProducts; }
      const LevelData<FArrayBox>& getFusionProducts() const { return m_fusionProducts; }
      
      virtual void applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                              const DomainGrid&           a_mesh,
                              const Real                  a_dt_sec ) const = 0;
      
      virtual void printParameters() const = 0;

   protected:

    bool m_verbosity;

    // container for storing fusion product densities (diagnostic purposes only)
    // Actual products are produced probabilistically, but we also compute the
    // product density as would be obtained in a continuum method. This is to have an
    // in-situ fusion diagnostic without creating many low weight particle in
    // situations where the density of products is much less than that of the plasma.
    mutable LevelData<FArrayBox> m_fusionProducts;

};

typedef RefCountedPtr<Scattering> ScatteringPtr;
typedef Vector<ScatteringPtr> ScatteringPtrVect;

#include "NamespaceFooter.H"

#endif
