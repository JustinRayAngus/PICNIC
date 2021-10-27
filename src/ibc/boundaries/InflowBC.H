#ifndef _INFLOWBC_H_
#define _INFLOWBC_H_

#include "DomainGrid.H"
#include "CodeUnits.H"
#include "JustinsParticle.H"

#include "NamespaceHeader.H"

class InflowBC
{
public:

   InflowBC( const int&             a_bdry_layout_index,
             const std::string&     a_bdry_name,
             const int&             a_bdry_dir,
             const Side::LoHiSide&  a_bdry_side,
             const std::string&     a_species_name,
             const Real&            a_species_mass,
             const DomainGrid&      a_mesh,
             const CodeUnits&       a_units );

   ~InflowBC() {;}

   const int& getBoundaryIndex() const {return m_bdry_layout_index;}
   
   void apply( List<JustinsParticle>&  a_pList,
         const FArrayBox&              a_Xec,
         const Box&                    a_box,
         const Real&                   a_cnormDt );

   void printParameters();
  
private:

   void parseParameters();
 
   int m_bdry_layout_index;   

   std::string m_bdry_name;
   int m_bdry_dir;
   int m_bdry_side;
   RealVect m_dX;
   Real m_X_bdry;

   std::string m_species_name;

   Real m_density;
   std::array<Real,3> m_velocity;
   std::array<Real,3> m_temperature;
   std::array<Real,3> m_beta_VT;
   std::array<Real,3> m_beta_U;
   
   int m_discrete_samples;
   Real m_weight;
   
};

typedef RefCountedPtr<InflowBC> InflowBCPtr;
typedef Vector<InflowBCPtr> InflowBCPtrVect;

#include "NamespaceFooter.H"

#endif
