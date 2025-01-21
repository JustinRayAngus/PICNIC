#ifndef _PIC_CHARGED_SPECIES_BC_H_
#define _PIC_CHARGED_SPECIES_BC_H_

#include "DomainGrid.H"
#include "ParmParse.H"
#include "List.H"
#include "ParticleData.H"
#include "JustinsParticle.H"
#include "ParticleBC.H"
#include "InflowBC.H"

#include "MeshInterp.H"

#include "NamespaceHeader.H"

class PicChargedSpeciesBC
{
   public:
   
   PicChargedSpeciesBC( const std::string&  a_species_name,
                 const Real          a_species_mass,
                 const int           a_species_charge,
                 const InterpType&   a_interpRhoToGrid,
                 const DomainGrid&   a_mesh,
                 const CodeUnits&    a_units,
                 const int           a_verbosity );

   ~PicChargedSpeciesBC();
   
   void apply( List<JustinsParticle>&     a_outcast_list,
               LevelData<NodeFArrayBox>&  a_surfaceCharge,
         const bool&                      a_intermediate_advance,
         const Real&                      a_time );
   
   void applyToRho( LevelData<FArrayBox>&  a_Rho );
   void applyToRho( LevelData<FluxBox>&  a_Rho );
   void applyToRho( LevelData<NodeFArrayBox>&  a_Rho );
   void applyToRhoInGhosts( LevelData<NodeFArrayBox>&  a_Rho );
   
   void removeOutflowParticles( LevelData<NodeFArrayBox>&  a_surfaceCharge );

   void zeroDeltas();
  
   const RealVect& getDeltaMassOut_lo() const { return m_delta_MassOut_lo; }
   const RealVect& getDeltaMassOut_hi() const { return m_delta_MassOut_hi; }
   const RealVect& getDeltaMomXOut_lo() const { return m_delta_MomXOut_lo; }
   const RealVect& getDeltaMomXOut_hi() const { return m_delta_MomXOut_hi; }
   const RealVect& getDeltaMomYOut_lo() const { return m_delta_MomYOut_lo; }
   const RealVect& getDeltaMomYOut_hi() const { return m_delta_MomYOut_hi; }
   const RealVect& getDeltaMomZOut_lo() const { return m_delta_MomZOut_lo; }
   const RealVect& getDeltaMomZOut_hi() const { return m_delta_MomZOut_hi; }
   const RealVect& getDeltaEnergyOut_lo() const { return m_delta_EnergyOut_lo; }
   const RealVect& getDeltaEnergyOut_hi() const { return m_delta_EnergyOut_hi; }
   
   const RealVect& getDeltaMassIn_lo() const { return m_delta_MassIn_lo; }
   const RealVect& getDeltaMassIn_hi() const { return m_delta_MassIn_hi; }
   const RealVect& getDeltaMomXIn_lo() const { return m_delta_MomXIn_lo; }
   const RealVect& getDeltaMomXIn_hi() const { return m_delta_MomXIn_hi; }
   const RealVect& getDeltaMomYIn_lo() const { return m_delta_MomYIn_lo; }
   const RealVect& getDeltaMomYIn_hi() const { return m_delta_MomYIn_hi; }
   const RealVect& getDeltaMomZIn_lo() const { return m_delta_MomZIn_lo; }
   const RealVect& getDeltaMomZIn_hi() const { return m_delta_MomZIn_hi; }
   const RealVect& getDeltaEnergyIn_lo() const { return m_delta_EnergyIn_lo; }
   const RealVect& getDeltaEnergyIn_hi() const { return m_delta_EnergyIn_hi; }
  
   const IntVect& getIsInflowBC_lo() const { return m_isInflowBC_lo; }
   const IntVect& getIsInflowBC_hi() const { return m_isInflowBC_hi; }
   
   const IntVect& getIsOutflowBC_lo() const { return m_isOutflowBC_lo; }
   const IntVect& getIsOutflowBC_hi() const { return m_isOutflowBC_hi; }
   
   Vector<List<JustinsParticle>>& getInflowListVect() {return m_inflow_list_vector;}
 
   void depositInflowOutflowJ( LevelData<EdgeDataBox>&    a_J_inPlane,
                               LevelData<NodeFArrayBox>&  a_J_virtual, 
                         const MeshInterp&                a_meshInterp, 
                         const InterpType                 a_interpJToGrid,
                         const Real                       a_cnormDt,
                         const bool                       a_from_explicit_solver );
   
   void createInflowParticles( const Real&  a_time,
                               const Real&  a_cnormDt,
                               const ParticleData<JustinsParticle>&  a_data );
   
   void injectInflowParticles( ParticleData<JustinsParticle>&  a_data,
                               LevelData<NodeFArrayBox>&       a_surfaceCharge );
      
   protected:
  
      void enforcePeriodic( List<JustinsParticle>&  a_list,
                      const int                     a_dir,
                      const Real                    a_leftEdge,
                      const Real                    a_rightEdge );
      
      void axis( List<JustinsParticle>&  a_list,
           const int                     a_dir );
      
      void symmetry_Lo( List<JustinsParticle>&  a_list,
                  const int                     a_dir,
                  const Real                    a_leftEdge );
      
      void symmetry_Hi( List<JustinsParticle>&  a_list,
                  const int                     a_dir,
                  const Real                    a_rightEdge );
      
      void outflow_Lo( List<JustinsParticle>&  a_list,
                       List<JustinsParticle>&  a_outflow_list,
                 const int                     a_dir,
                 const Real                    a_leftEdge );
      
      void outflow_Hi( List<JustinsParticle>&  a_list,
                       List<JustinsParticle>&  a_outflow_list,
                 const int                     a_dir,
                 const Real                    a_rightEdge );
      
      void inflow_Lo( List<JustinsParticle>&  a_list,
                      List<JustinsParticle>&  a_inflow_list,
                const int                     a_dir,
                const Real                    a_leftEdge,
                const bool                    a_intermediate_advance );
      
      void inflow_Hi( List<JustinsParticle>&  a_list,
                      List<JustinsParticle>&  a_inflow_list,
                const int                     a_dir,
                const Real                    a_rightEdge,
                const bool                    a_intermediate_advance );

      void depositSurfaceChargeLo( FArrayBox&              a_sigma,
                             const Box&                    a_bdry_box,
		             const List<JustinsParticle>&  a_list,
                             const int                     a_bdry_dir,
                             const Real                    a_leftEdge );
      
      void depositSurfaceChargeHi( FArrayBox&              a_sigma,
                             const Box&                    a_bdry_box,
		             const List<JustinsParticle>&  a_list,
                             const int                     a_bdry_dir,
                             const Real                    a_rightEdge );
      
      void extractSurfaceChargeLo( FArrayBox&              a_sigma,
                             const Box&                    a_bdry_box,
		             const List<JustinsParticle>&  a_list,
                             const int                     a_bdry_dir,
                             const Real                    a_leftEdge,
                             const bool                    a_intermediate_advance );
      
      void extractSurfaceChargeHi( FArrayBox&              a_sigma,
                             const Box&                    a_bdry_box,
		             const List<JustinsParticle>&  a_list,
                             const int                     a_bdry_dir,
                             const Real                    a_rightEdge,
                             const bool                    a_intermediate_advance );
      
      void remove_Lo( List<JustinsParticle>&  a_list,
                const int                     a_bdry_dir,
                const Real                    a_leftEdge );
      
      void remove_Hi( List<JustinsParticle>&  a_list,
                const int                     a_bdry_dir,
                const Real                    a_rightEdge );
      
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      Vector<int> m_bc_binary;
      Vector<int> m_periodic_bc_binary;
   
   private:
      
      string m_species_name;
      Real m_species_mass;
      int m_species_charge;
      int m_interp_order;

      const DomainGrid& m_mesh;
   
      int m_verbosity;

      Vector<List<JustinsParticle>>  m_outflow_list_vector;
      Vector<List<JustinsParticle>>  m_inflow_list_vector;
  
      // boundary probes used for diagnostics
      RealVect m_delta_MassOut_lo;   // delta mass of outflow particles at lo side 
      RealVect m_delta_MassOut_hi;   // delta mass of outflow particles at hi side
      RealVect m_delta_MomXOut_lo;   // delta X-momentum of outflow particles at lo side 
      RealVect m_delta_MomXOut_hi;   // delta X-momentum of outflow particles at hi side
      RealVect m_delta_MomYOut_lo;   // delta Y-momentum of outflow particles at lo side 
      RealVect m_delta_MomYOut_hi;   // delta Y-momentum of outflow particles at hi side
      RealVect m_delta_MomZOut_lo;   // delta Z-momentum of outflow particles at lo side 
      RealVect m_delta_MomZOut_hi;   // delta Z-momentum of outflow particles at hi side
      RealVect m_delta_EnergyOut_lo; // delta energy of outflow particles at lo side 
      RealVect m_delta_EnergyOut_hi; // delta energy of outflow particles at hi side
      // 
      RealVect m_delta_MassIn_lo;   // delta mass of inflow particles at lo side 
      RealVect m_delta_MassIn_hi;   // delta mass of inflow particles at hi side
      RealVect m_delta_MomXIn_lo;   // delta X-momentum of inflow particles at lo side 
      RealVect m_delta_MomXIn_hi;   // delta X-momentum of inflow particles at hi side
      RealVect m_delta_MomYIn_lo;   // delta Y-momentum of inflow particles at lo side 
      RealVect m_delta_MomYIn_hi;   // delta Y-momentum of inflow particles at hi side
      RealVect m_delta_MomZIn_lo;   // delta Z-momentum of inflow particles at lo side 
      RealVect m_delta_MomZIn_hi;   // delta Z-momentum of inflow particles at hi side
      RealVect m_delta_EnergyIn_lo; // delta energy of inflow particles at lo side 
      RealVect m_delta_EnergyIn_hi; // delta energy of outflow particles at hi side

      IntVect m_isInflowBC_lo = IntVect::Zero;
      IntVect m_isInflowBC_hi = IntVect::Zero;     
      IntVect m_isOutflowBC_lo = IntVect::Zero;
      IntVect m_isOutflowBC_hi = IntVect::Zero;     

      InflowBCPtrVect m_inflow_bc_ptr_vect;

      void parseParameters( ParmParse&  a_pp,
                      const CodeUnits&  a_units );
   
      void printParameters() const;

};


#include "NamespaceFooter.H"

#endif
