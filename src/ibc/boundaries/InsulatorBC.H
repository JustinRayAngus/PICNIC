#ifndef _INSULATORBC_H_
#define _INSULATORBC_H_

#include "TimeFunction.H"
#include "LevelData.H"
#include "DomainGrid.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

class InsulatorBC
{
public:

   InsulatorBC( const int&          a_bdry_layout_index,
                const std::string&  a_bdry_name,
                const DomainGrid&   a_mesh );

   ~InsulatorBC() {;}

   const int& getBvComp() const {return m_Bv_comp;}
   
   const RealVect& getXmin() const {return m_Xmin_insulator;}
  
   const RealVect& getXmax() const {return m_Xmax_insulator;}
  
   const bool isAbsorbing() const { return m_absorbing; }

   const LevelData<FArrayBox>& getICbinary_cc() const {return m_ICbinary_cc;}
  
   const LevelData<FluxBox>& getICbinary_fc() const {return m_ICbinary_fc;}
  
   const LevelData<EdgeDataBox>& getICbinary_ec() const {return m_ICbinary_ec;}
  
   const LevelData<NodeFArrayBox>& getICbinary_nc() const {return m_ICbinary_nc;}
   
   Real getTimeFunctionValue( const Real a_time ) const 
   { 
      Real this_ft;
      m_timeFunction->getValue(this_ft,a_time);
      return this_ft;
   }
   
private:

   void parseParameters(const std::string&  a_bdry_name);
   void defineInsulatorBinaries( const DomainGrid&  a_mesh );
 
   int m_Bv_comp = 0;
   int m_bdry_layout_index;   

   bool m_absorbing = false;

   RealVect m_Xmin_insulator; // should create a ProblemDomain for insulator...
   RealVect m_Xmax_insulator;
   
   RefCountedPtr<TimeFunction> m_timeFunction;

   LevelData<FArrayBox> m_ICbinary_cc;     // insulator conductor binary on cells
   LevelData<FluxBox> m_ICbinary_fc;       // insulator conductor binary on faces
   LevelData<EdgeDataBox> m_ICbinary_ec;   // insulator conductor binary on edges
   LevelData<NodeFArrayBox> m_ICbinary_nc; // insulator conductor binary on nodes
      
};

#include "NamespaceFooter.H"

#endif
