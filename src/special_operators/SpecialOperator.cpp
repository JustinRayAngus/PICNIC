#include "SpecialOperator.H"

#include "NamespaceHeader.H"

SpecialOperator::SpecialOperator( const DomainGrid&  a_mesh,
                                  const CodeUnits&   a_units,
                                  const int&         a_verbosity )
   : m_mesh(a_mesh),
     m_units(a_units),
     m_verbosity(a_verbosity)
{
}

#include "NamespaceFooter.H"
