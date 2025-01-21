/***
 *
 * Base species class
 *
***/

#ifndef Species_h
#define Species_h

#include "ParmParse.H"
#include "DomainGrid.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "RefCountedPtr.H"

#include "CodeUnits.H"

#include "NamespaceHeader.H"

typedef enum { PIC_CHARGED, 
               PIC_PHOTON } SpeciesType;

class Species
{

public:

   Species( ParmParse&   a_ppsp,
      const int          a_species,
      const string       a_name,
      const DomainGrid&  a_mesh );

   virtual ~Species() {;}

   virtual SpeciesType getSpeciesType() const = 0;

   virtual void initialize( const CodeUnits&    a_units,
                            const Real          a_time,
                            const std::string&  a_restart_file_name ) = 0;

   const std::string& name() const { return m_name; }
   
   const int species_num() const { return m_species; }

   const Real mass() const { return m_mass; }

   const int charge() const { return m_charge; }

   bool motion() const { return m_motion; }

   bool forces() const { return m_forces; }

   bool scatter() const { return m_scatter; }

   bool isSpecies( const string&  name ) const;  

protected:

   int  m_species;     // species number identifier 
   std::string m_name; // species name

   Real m_mass = 0.0;  // species mass / electron mass
   int  m_charge = 0;  // species charge / electron charge

   bool m_motion = true;
   bool m_forces = true;
   bool m_scatter = false;

   const DomainGrid& m_mesh;

};

typedef RefCountedPtr<Species> SpeciesPtr;
typedef Vector<SpeciesPtr> SpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
