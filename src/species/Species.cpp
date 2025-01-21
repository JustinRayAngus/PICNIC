#include "Species.H"

#include "NamespaceHeader.H"

Species::Species( ParmParse&   a_ppspc,
            const int          a_species,
            const string       a_name,
            const DomainGrid&  a_mesh )
   : m_species(a_species),
     m_name(a_name),
     m_mesh(a_mesh)
{
   a_ppspc.get( "mass", m_mass );
   a_ppspc.get( "charge", m_charge ); 
   a_ppspc.query( "motion", m_motion );
   a_ppspc.query( "forces", m_forces );
   a_ppspc.query( "scatter", m_scatter );
}

bool Species::isSpecies( const string&  a_name ) const
{
   if (name()==a_name) { return true; }
   else { return false; }
}

#include "NamespaceFooter.H"

