#ifndef _FIELD_BC_H_
#define _FIELD_BC_H_

#include "DomainGrid.H"
#include "ParmParse.H"
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
                  const Real                     a_time ) const;
      
      void applyNodeBC( LevelData<NodeFArrayBox>&  a_dst,
                  const Real                       a_time ) const;

      void setFluxBC( LevelData<FluxBox>&  a_dst,
                const LevelData<FluxBox>&  a_src,
                const Real                 a_time );
      
      void setEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                const LevelData<EdgeDataBox>&  a_src,
                const Real                     a_time );
      
      void applyCellPCMask( LevelData<FArrayBox>&  a_dst,
                      const Real                   a_time );
       
      void applyFluxPCMask( LevelData<FluxBox>&  a_dst,
                      const Real                 a_time );
      
      void applyEdgePCMask( LevelData<EdgeDataBox>&  a_dst,
                      const Real                     a_time ) const;
      
      void applyNodePCMask( LevelData<NodeFArrayBox>&  a_dst,
                      const Real                       a_time ) const;

      void applyOnAxisCurlBC( LevelData<EdgeDataBox>&  a_curlB,
                        const LevelData<FArrayBox>&    a_B );
      
      void applyOnAxisCurlBC( LevelData<NodeFArrayBox>&  a_curlB,
                        const LevelData<FArrayBox>&      a_B );

      void applyOnAxisDivBC( LevelData<NodeFArrayBox>&  a_divE,
                       const LevelData<EdgeDataBox>&    a_E );

      void computeIntSdA( RealVect&                  a_intSdA_lo, 
                          RealVect&                  a_intSdA_hi,
                    const LevelData<EdgeDataBox>&    a_E, 
                    const LevelData<FluxBox>&        a_B,
                    const LevelData<NodeFArrayBox>&  a_Ev,
                    const LevelData<FArrayBox>&      a_Bv );
    
      void applyToJ( LevelData<EdgeDataBox>&    a_J_inPlane,
                     LevelData<NodeFArrayBox>&  a_J_virtual ) const;
  
   protected:
   
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      bool m_conservative_wall; // if true, E for insulator/conductor BCs will be set
                                // such that energy is conserved  
 
      Vector<RefCountedPtr<InsulatorBC>> m_InsulatorBC;
   
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
