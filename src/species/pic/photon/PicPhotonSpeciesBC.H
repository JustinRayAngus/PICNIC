#ifndef _PIC_PHOTON_SPECIES_BC_H_
#define _PIC_PHOTON_SPECIES_BC_H_

#include "DomainGrid.H"
#include "ParmParse.H"
#include "List.H"
#include "ParticleData.H"
#include "PhotonParticle.H"
#include "ParticleBC.H"
#include "CodeUnits.H"

#include "MeshInterp.H"

#include "NamespaceHeader.H"

class PicPhotonSpeciesBC
{
   public:
   
   PicPhotonSpeciesBC( const std::string&  a_species_name,
                       const DomainGrid&   a_mesh,
                       const CodeUnits&    a_units,
                       const int           a_verbosity );

   ~PicPhotonSpeciesBC();
   
   void apply( List<PhotonParticle>&  a_outcast_list,
         const bool                    a_intermediate_advance,
         const Real                    a_time );
   
   void removeOutflowParticles();

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
   
   const IntVect& getIsOutflowBC_lo() const { return m_isOutflowBC_lo; }
   const IntVect& getIsOutflowBC_hi() const { return m_isOutflowBC_hi; }

   protected:

      void enforcePeriodic( List<PhotonParticle>&  a_list,
                      const int                     a_dir,
                      const Real                    a_leftEdge,
                      const Real                    a_rightEdge );

      void axis( List<PhotonParticle>&  a_list,
           const int                     a_dir );

      void symmetry_Lo( List<PhotonParticle>&  a_list,
                  const int                     a_dir,
                  const Real                    a_leftEdge );

      void symmetry_Hi( List<PhotonParticle>&  a_list,
                  const int                     a_dir,
                  const Real                    a_rightEdge );

      void outflow_Lo( List<PhotonParticle>&  a_list,
                       List<PhotonParticle>&  a_outflow_list,
                 const int                     a_dir,
                 const Real                    a_leftEdge );

      void outflow_Hi( List<PhotonParticle>&  a_list,
                       List<PhotonParticle>&  a_outflow_list,
                 const int                     a_dir,
                 const Real                    a_rightEdge );

      void remove_Lo( List<PhotonParticle>&  a_list,
                const int                     a_bdry_dir,
                const Real                    a_leftEdge );
      
      void remove_Hi( List<PhotonParticle>&  a_list,
                const int                     a_bdry_dir,
                const Real                    a_rightEdge );
      
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      Vector<int> m_bc_binary;
      Vector<int> m_periodic_bc_binary;

   private:

      string m_species_name;
      int m_interp_order;

      const DomainGrid& m_mesh;

      int m_verbosity;

      Vector<List<PhotonParticle>>  m_outflow_list_vector;
  
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

      IntVect m_isOutflowBC_lo = IntVect::Zero;
      IntVect m_isOutflowBC_hi = IntVect::Zero;

      void parseParameters( ParmParse&  a_pp,
                      const CodeUnits&  a_units );

      void printParameters() const;

};


#include "NamespaceFooter.H"

#endif
