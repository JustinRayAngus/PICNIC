#ifndef _FLUID_VAR_BC_H_
#define _FLUID_VAR_BC_H_

#include "EdgeDataBox.H"
#include "DomainGrid.H"
//#include "InsulatorConductorBC.H"

#include "NamespaceHeader.H"

class FieldBC
{
   public:
   
      FieldBC( const std::string&  variable_name,
               const DomainGrid&   a_mesh,
               const int           verbosity );

      ~FieldBC();
   
      void applyCellBC( LevelData<FArrayBox>&   a_dst,
                  const Real                    a_time );
       
      void applyFluxBC( LevelData<FluxBox>&     a_dst,
                  const Real                    a_time );
      
      void applyEdgeBC( LevelData<EdgeDataBox>&  a_dst,
                  const Real                     a_time );

      /*
      void setFluxBC( const ElectroMagneticFields&  species_phys,
                              LevelData<FluxBox>&  dst,
                        const LevelData<FluxBox>&  src,
                        const Real                 time );
      
      void setEdgeBC( const ElectroMagneticFields&     species_phys,
                              LevelData<EdgeDataBox>& dst,
                        const LevelData<EdgeDataBox>& src,
                        const Real                    time );
      
      void setOnAxisCurlBC( LevelData<EdgeDataBox>&  a_curl_covar,
                              const LevelData<FArrayBox>&    a_By_phys,
                              const ElectroMagneticFields&            a_species_comp );
      
      void setInsulatorBC( const ElectroMagneticFields&    a_species_comp,
                                   LevelData<FArrayBox>&  a_dst,
                             const LevelData<FArrayBox>&  a_src,
                             const Real                   a_time );
   
      */

      void setDomainTypesAndNames( const BoundaryBoxLayoutPtrVect&  a_bdry_layout );
      
      //bool isInsulatorConductorBC() const {return m_insulator_conductor_bc;}

      bool isForVariable( const std::string& a_variable_name ) const
      {
         return (m_variable_name==a_variable_name);
      }
   
   protected:
   
      std::string getBCType( const int dir,
                             const int side );
      
      std::string getBdryName( const int dir,
                               const int side );
   
      void setNames();

      const std::string m_variable_name;
   
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      Vector<std::string> m_domain_bdry_name;
      Vector<std::string> m_domain_bc_type;
   
      //bool m_insulator_conductor_bc;
   
   private:
  
   const DomainGrid& m_mesh;
   
   int m_verbosity;

   void parseParameters( ParmParse&  a_pp );
   
   void printParameters() const;

   enum {INVALID=-1,
  	DIR0_LOWER, DIR0_UPPER,
#if CH_SPACEDIM>=2
	DIR1_LOWER, DIR1_UPPER,
#endif
#if CH_SPACEDIM==3
	DIR2_LOWER, DIR2_UPPER,
#endif
	NUM_BOUNDARIES};

};


#include "NamespaceFooter.H"

#endif
