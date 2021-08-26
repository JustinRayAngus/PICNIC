#ifndef _FIELD_BC_H_
#define _FIELD_BC_H_

#include "DomainGrid.H"
#include "InsulatorBC.H"

#include "NamespaceHeader.H"

class FieldBC
{
   public:
   
      FieldBC( const DomainGrid&   a_mesh,
               const int           a_verbosity );

      ~FieldBC();
   
      void applyCellBC( LevelData<FArrayBox>&  a_dst,
                  const Real                   a_time );
       
      void applyFluxBC( LevelData<FluxBox>&  a_dst,
                  const Real                 a_time );
      
      void applyEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                  const Real                     a_time );
      
      void applyNodeBC( LevelData<NodeFArrayBox>&  a_dst,
                  const Real                       a_time );

      void setFluxBC( LevelData<FluxBox>&  a_dst,
                const LevelData<FluxBox>&  a_src,
                const Real                 a_time );
      
      void setEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                const LevelData<EdgeDataBox>&  a_src,
                const Real                     a_time );
      
   protected:
   
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
   
      RefCountedPtr<InsulatorBC> m_InsulatorBC;
   
   private:
  
   const DomainGrid& m_mesh;
   
   int m_verbosity;

   void parseParameters( ParmParse&  a_pp );
   
   void printParameters() const;

/*
   enum {INVALID=-1,
  	DIR0_LOWER, DIR0_UPPER,
#if CH_SPACEDIM>=2
	DIR1_LOWER, DIR1_UPPER,
#endif
#if CH_SPACEDIM==3
	DIR2_LOWER, DIR2_UPPER,
#endif
	NUM_BOUNDARIES};
*/

};


#include "NamespaceFooter.H"

#endif
