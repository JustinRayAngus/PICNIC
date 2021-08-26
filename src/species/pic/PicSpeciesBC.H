#ifndef _PIC_SPECIES_BC_H_
#define _PIC_SPECIES_BC_H_

#include "DomainGrid.H"
#include "List.H"
#include "ParticleData.H"
#include "JustinsParticle.H"
#include "ParticleBC.H"
#include "InflowBC.H"

#include "NamespaceHeader.H"

class PicSpeciesBC
{
   public:
   
   PicSpeciesBC( const std::string&  a_species_name,
                 const Real&         a_species_mass,
                 const DomainGrid&   a_mesh,
                 const CodeUnits&    a_units,
                 const int           a_verbosity );

   ~PicSpeciesBC();
   
   void apply( List<JustinsParticle>&  a_outcast_list,
         const bool&                   a_intermediate_advance,
         const Real&                   a_time );
   
   void removeOutflowParticles();
   
   void repositionInflowParticles( ParticleData<JustinsParticle>&  a_data );
   
   void repositionOutflowParticles( List<JustinsParticle>&  a_outcast_list );
   
   void createInflowParticles( const Real&  a_time,
                               const Real&  a_dt );
   
   protected:
  
      void enforcePeriodic( List<JustinsParticle>&  a_list,
                      const int&      a_dir,
                      const Real&     a_leftEdge,
                      const Real&     a_rightEdge );
      
      void applyBC_Lo( List<JustinsParticle>&  a_list,
                       List<JustinsParticle>&  a_inflow_list,
                       List<JustinsParticle>&  a_outflow_list,
                 const std::string&            a_bc_type,
                 const int&                    a_dir,
                 const Real&                   a_leftEdge,
                 const bool&                   a_intermediate_advance );
      
      void applyBC_Hi( List<JustinsParticle>&  a_list,
                       List<JustinsParticle>&  a_inflow_list,
                       List<JustinsParticle>&  a_outflow_list,
                 const std::string&            a_bc_type,
                 const int&                    a_dir,
                 const Real&                   a_rightEdge,
                 const bool&                   a_intermediate_advance );
      
      
      void symmetry_Lo( List<JustinsParticle>&  a_list,
                  const int&                    a_dir,
                  const Real&                   a_leftEdge );
      
      void symmetry_Hi( List<JustinsParticle>&  a_list,
                  const int&                    a_dir,
                  const Real&                   a_rightEdge );
      
      void outflow_Lo( List<JustinsParticle>&  a_list,
                       List<JustinsParticle>&  a_outflow_list,
                 const int&                    a_dir,
                 const Real&                   a_leftEdge );
      
      void outflow_Hi( List<JustinsParticle>&  a_list,
                       List<JustinsParticle>&  a_outflow_list,
                 const int&                    a_dir,
                 const Real&                   a_rightEdge );
      
      void inflow_Lo( List<JustinsParticle>&  a_list,
                      List<JustinsParticle>&  a_inflow_list,
                const int&                    a_dir,
                const Real&                   a_leftEdge,
                const bool&                   a_intermediate_advance );
      
      void inflow_Hi( List<JustinsParticle>&  a_list,
                      List<JustinsParticle>&  a_inflow_list,
                const int&                    a_dir,
                const Real&                   a_rightEdge,
                const bool&                   a_intermediate_advance );
      
      void remove_Lo( List<JustinsParticle>&  a_list,
                const int&                    a_dir,
                const Real&                   a_leftEdge );
      
      void remove_Hi( List<JustinsParticle>&  a_list,
                const int&                    a_dir,
                const Real&                   a_rightEdge );
      
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      Vector<int> m_bc_binary;
   
   private:
      
      const string m_species_name;
      const Real m_species_mass;  

      const DomainGrid& m_mesh;
   
      int m_verbosity;

      Vector<List<JustinsParticle>>  m_outflow_list_vector;
      Vector<List<JustinsParticle>>  m_inflow_list_vector;

      InflowBCPtrVect m_inflow_bc_ptr_vect;

      void parseParameters( ParmParse&  a_pp,
                      const CodeUnits&  a_units );
   
      void printParameters() const;

};


#include "NamespaceFooter.H"

#endif
