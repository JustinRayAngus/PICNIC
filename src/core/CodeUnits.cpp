#include "CodeUnits.H"
#include "Constants.H"
#include <assert.h>

#include "NamespaceHeader.H"

inline void getPosDefUnit(
   Real& a_val, const std::string &a_name, ParmParse& a_pp )
{
   Real val = 0.0;
   a_pp.get( a_name.c_str(), val ); // a_name must exist or abort
   assert( val > 0.0 );
   a_val = val;
}

CodeUnits::CodeUnits( ParmParse& a_parm_parse )
{
   // Fundamental Characteristic Scales
   ParmParse ppunits( "units" );
   getPosDefUnit( m_scale[TEMPERATURE],    "temperature",    ppunits );
   getPosDefUnit( m_scale[LENGTH],         "length",         ppunits );
   getPosDefUnit( m_scale[TIME],           "time",           ppunits );

   // Universal Constants
   //Real pi = Constants::PI;
   //m_scale[CHARGE]       = Constants::ELEMENTARY_CHARGE;
   //m_scale[BOLTZMANN]    = Constants::BOLTZMANN_CONSTANT;
   //m_scale[PERMITTIVITY] = Constants::VACUUM_PERMITTIVITY;
   //m_scale[PERMEABILITY] = pi * 4.0e-7;

}

void CodeUnits::printParameters( const int a_procID) const
{
   if(!a_procID) {
      cout << "====================== Fundamental Code Units ======================" << endl;
      cout << "  TEMPERATURE  [eV]: " << m_scale[TEMPERATURE] << " (NOT USED YET)" << endl;
      cout << "  LENGTH        [m]: " << m_scale[LENGTH]      << " (NOT USED YET)" << endl;
      cout << "  TIME          [s]: " << m_scale[TIME]        << " (NOT USED YET)" << endl;
      cout << "====================================================================" << endl;
      cout << endl;
   }
}

#include "NamespaceFooter.H"
