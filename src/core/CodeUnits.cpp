#include "CodeUnits.H"
#include "ParmParse.H"
#include "PicnicConstants.H"
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

CodeUnits::CodeUnits()
{
   // Fundamental Characteristic Scales
   ParmParse ppunits( "units" );
   getPosDefUnit( m_scale[NUMBER_DENSITY], "number_density", ppunits );
   getPosDefUnit( m_scale[TEMPERATURE],    "temperature",    ppunits );
   getPosDefUnit( m_scale[LENGTH],         "length",         ppunits );
   getPosDefUnit( m_scale[TIME],           "time",           ppunits );

   CH_assert(m_scale[NUMBER_DENSITY]==1.0);
   CH_assert(m_scale[TEMPERATURE]==1.0);
   
   m_CVAC_NORM = Constants::CVAC*m_scale[TIME]/m_scale[LENGTH];
   m_WP_NORM = Constants::QE*sqrt(m_scale[NUMBER_DENSITY]/Constants::EP0/Constants::ME)*m_scale[TIME];

   Real Escale = 1.0e5; // kV/cm
   m_scale[ELECTRIC_FIELD] = Escale; // 1 kV/cm = 1.0e5 V/m
   m_scale[MAGNETIC_FIELD] = m_scale[ELECTRIC_FIELD]/Constants::CVAC; // 1 kV/cm/CVAC = 3.3356e-4 T

}

void CodeUnits::printParameters() const
{
   if(procID()>0) return;
   cout << "====================== Fundamental Code Units ======================" << endl;
   cout << "  NUMBER_DENISTY [1/m^3]: " << m_scale[NUMBER_DENSITY] << endl;
   cout << "  TEMPERATURE       [eV]: " << m_scale[TEMPERATURE] << endl;
   cout << "  LENGTH             [m]: " << m_scale[LENGTH]      << endl;
   cout << "  TIME               [s]: " << m_scale[TIME]        << endl;
   cout << "  ELECTRIC_FIELD   [V/m]: " << m_scale[ELECTRIC_FIELD] << endl;
   cout << "  MAGNETIC_FIELD     [T]: " << m_scale[MAGNETIC_FIELD] << endl;
   cout << "====================================================================" << endl;
   cout << endl;

}

#include "NamespaceFooter.H"
