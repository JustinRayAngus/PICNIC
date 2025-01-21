#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MESHINTERP_H_
#define _MESHINTERP_H_

#include "Box.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "Vector.H"
#include "RealVect.H"
#include "MeshInterpF_F.H"
#include "MeshInterpChargeConservingF_F.H"
#include "MeshInterpMassMatrixF_F.H"

#include "NamespaceHeader.H"

typedef enum {number, density, momentum, energy, energyOffDiag, energyFlux} MomentType;
typedef enum {photon_number, photon_density, photon_momentum, photon_energy} PhotonMomentType;
typedef enum {UNKNOWN, NGP, CIC, TSC, W4, CC0, CC1} InterpType;

/**
   This class is used for handling the interaction between the particles and the mesh
   in PIC calculations. It handles depositing particle masses or charges onto the 
   mesh, and interpolating the result of the force solve back to the particle positions.
 **/
class MeshInterp
{
public:

  MeshInterp();

  /// Full constructor 
  MeshInterp(const Box&  a_domain,
             const int a_ghosts,
             const RealVect& a_dx,
             const RealVect& a_domainLeftEdge,
             const RealVect& a_domainRightEdge);

  /// Define function 
  void define(const Box&  a_domain,
              const int a_ghosts,
              const RealVect& a_dx,
              const RealVect& a_domainLeftEdge,
              const RealVect& a_domainRightEdge);

  void setBCcheckLo( const IntVect a_bc_check_lo ) {m_bc_check_lo = a_bc_check_lo;}
  void setBCcheckHi( const IntVect a_bc_check_hi ) {m_bc_check_hi = a_bc_check_hi;}

  /// Deposit the particles in a_particleList onto a_rho.
  /// The type of interpolation to perform is passed in as an argument.
  template <class P>
  void deposit(const List<P>&  a_particleList,
               FArrayBox&      a_rho,
               const InterpType&  a_interpType) const; // JRA const

  template <class P>
  void depositCurrent( FArrayBox&   a_rhoVx,
                       FArrayBox&   a_rhoVy,
                       FArrayBox&   a_rhoVz,
                 const List<P>&     a_particleList,
                 const InterpType&  a_interpType,
                 const Real         a_cnormDt,
                 const int          a_axisymm_car_push,
                 const bool         a_from_explicit_solver ) const;
  
  template <class P>
  void depositMassMatrices( FArrayBox&   a_sigma_xx,
                            FArrayBox&   a_sigma_xy,
                            FArrayBox&   a_sigma_xz,
                            FArrayBox&   a_sigma_yx,
                            FArrayBox&   a_sigma_yy,
                            FArrayBox&   a_sigma_yz,
                            FArrayBox&   a_sigma_zx,
                            FArrayBox&   a_sigma_zy,
                            FArrayBox&   a_sigma_zz,
                            FArrayBox&   a_Jx,
                            FArrayBox&   a_Jy,
                            FArrayBox&   a_Jz,
                      const FArrayBox&   a_Bx,
                      const FArrayBox&   a_By,
                      const FArrayBox&   a_Bz,
                      const Real         a_qovs,
                      const Real         a_alphas,
                      const Real         a_cnormDt,
                      const int          a_inertial_type,
                      const bool         a_anticylic,
                      const List<P>&     a_particleList,
                      const InterpType&  a_interpType ) const;
  
  template <class P>
  void depositVelocity( FArrayBox&   a_rhoV,
                  const int          a_dir,
                  const List<P>&     a_particleList,
                  const InterpType&  a_interpType,
                  const int          a_comp=0 ) const;
  
  /// Interpolate the field given by a_field onto the particles in a_particleList.
  /// a_field must have SpaceDim components. The type of interpolation to perform
  /// is passed in as an argument.
  template <class P>
  void interpolate(List<P>& a_particleList,
                   const FArrayBox& a_field,
                   InterpType& a_interpType);
  
  // Interpolate a_electric_field onto the particles in a_particleList.
  template <class P>
  void interpolateEMfieldsToPart( List<P>&     a_particleList,
                            const FArrayBox&   a_Ex,
                            const FArrayBox&   a_Ey,
                            const FArrayBox&   a_Ez,
                            const FArrayBox&   a_Bx,
                            const FArrayBox&   a_By,
                            const FArrayBox&   a_Bz,
                            const InterpType&  a_interpType,
                            const int          a_blank_B=0 ) const;
  
  template <class P>
  void interpolateEMfieldsToPart( P&           a_particle,
                            const FArrayBox&   a_Ex,
                            const FArrayBox&   a_Ey,
                            const FArrayBox&   a_Ez,
                            const FArrayBox&   a_Bx,
                            const FArrayBox&   a_By,
                            const FArrayBox&   a_Bz,
                            const InterpType&  a_interpType ) const;
  
  template <class P>
  void interpolateEfieldToPart( List<P>&     a_particleList,
                          const int          a_dir,
                          const FArrayBox&   a_field,
                          const InterpType&  a_interpType,
                          const int          a_comp=0) const;
  
  // Interpolate a_magnetic_field onto the particles in a_particleList.
  template <class P>
  void interpolateBfieldToPart( List<P>&     a_particleList,
                          const int          a_dir,
                          const FArrayBox&   a_field,
                          const InterpType&  a_interpType,
                          const int          a_comp=0 ) const;
  
  //
  //
  //

  /// Deposit the particles in a_particleList onto a_rho.
  /// The type of interpolation to perform is passed in as an argument.
  template <class P>
  void moment( FArrayBox&  a_moment,
         const List<P>&    a_particleList,
         const Real&       a_species_mass,
         const MomentType  a_momentType ) const; // JRA const

  /// Deposit photon particles in a_particleList onto a_rho.
  /// The type of interpolation to perform is passed in as an argument.
  template <class P>
  void photonMoment( FArrayBox&  a_moment,
               const List<P>&    a_particleList,
               const PhotonMomentType  a_momentType ) const; // JRA const

   // set the particle weights using a grid-profile for density
   //
   template <class P>
   void setWeightFromGridProfile( List<P>&    a_particleList,
                            const FArrayBox&  a_density,
                            const Real        a_partsPerCell ) const;
  
  //
  //
  //
  
  template <class P>
  void depositCurrent_testing( FArrayBox&   a_rhoVx,
                               FArrayBox&   a_rhoVy,
                               FArrayBox&   a_rhoVz,
                         const List<P>&     a_particleList,
                         const InterpType&  a_interpType ) const;

  template <class P>
  void interpolateEMfieldsToPart_testing( List<P>&     a_particleList,
                              const FArrayBox&   a_Ex,
                              const FArrayBox&   a_Ey,
                              const FArrayBox&   a_Ez,
                              const FArrayBox&   a_Bx,
                              const FArrayBox&   a_By,
                              const FArrayBox&   a_Bz,
                              const InterpType&  a_interpType ) const;

private:
  
  template <class P>
  void interpolateEMfieldsToPart_CIC( List<P>&     a_particleList,
                                const FArrayBox&   a_Ex,
                                const FArrayBox&   a_Ey,
                                const FArrayBox&   a_Ez,
                                const FArrayBox&   a_Bx,
                                const FArrayBox&   a_By,
                                const FArrayBox&   a_Bz ) const;
  
  template <class P>
  void interpolateBfieldsToPart_CIC( List<P>&     a_particleList,
                               const FArrayBox&   a_Bx,
                               const FArrayBox&   a_By,
                               const FArrayBox&   a_Bz ) const;

  template <class P>
  void interpolateEfieldDirToPart_CIC( List<P>&     a_particleList,
                                 const FArrayBox&   a_E,
                                 const int          a_dir,
                                 const int          a_comp=0 ) const;
  
  template <class P>
  void interpolateEMfieldsToPart_TSC( List<P>&     a_particleList,
                                const FArrayBox&   a_Ex,
                                const FArrayBox&   a_Ey,
                                const FArrayBox&   a_Ez,
                                const FArrayBox&   a_Bx,
                                const FArrayBox&   a_By,
                                const FArrayBox&   a_Bz ) const;
  template <class P>
  void interpolateEToPart_CC0( List<P>&    a_particleList,
                         const FArrayBox&  a_Ex ) const;
  
  template <class P>
  void interpolateEToPart_CC0( List<P>&    a_particleList,
                         const FArrayBox&  a_Ex,
                         const FArrayBox&  a_Ey,
                         const FArrayBox&  a_Ez ) const;
  
  template <class P>
  void interpolateEToPart_CC1( List<P>&    a_particleList,
                         const FArrayBox&  a_Ex ) const;
  
  template <class P>
  void interpolateEToPart_CC1( List<P>&    a_particleList,
                         const FArrayBox&  a_Ex,
                         const FArrayBox&  a_Ey ) const;
  
  /// a wrapper function that deposits the contribution of a single 
  //  particle to the specified moment using the 
  /// interpolating functions defined in MeshInterpF.ChF
  void momentParticle( FArrayBox&   a_moment,
                 const RealVect&    a_domainLeftEdge,
                 const RealVect&    a_dx,
                 const RealVect&    a_position,
                 const std::array<Real,3>&   a_velocity,
                 const Real         a_weight,
                 const Real         a_species_mass,
                 const MomentType   a_momentType ) const; // JRA const;
  
  void momentPhotonParticle( FArrayBox&   a_moment,
                       const RealVect&    a_domainLeftEdge,
                       const RealVect&    a_dx,
                       const RealVect&    a_position,
                       const std::array<Real,3>&   a_velocity,
                       const Real         a_weight,
                       const PhotonMomentType  a_momentType ) const; // JRA const;

  void setParticleWeight( Real&       a_weight,
                    const Real&       a_partsPerCell,
                    const FArrayBox&  a_density,
                    const RealVect&   a_domainLeftEdge,
                    const RealVect&   a_dx,
                    const RealVect&   a_position ) const;


  /// a wrapper function that deposits a single particle using the 
  /// interpolating functions defined in MeshInterpF.ChF
  void depositParticle( FArrayBox&   a_rho,
                  const RealVect&    a_domainLeftEdge,
                  const RealVect&    a_dx,
                  const RealVect&    a_position,
                  const Real&        a_weight,
                  const Real&        a_kernal,
                  const IntVect&     a_stag,
                  const InterpType&  a_interpType,
                  const int          a_comp=0) const; // JRA const;
  
  /// a wrapper function that interpolates the field to a single
  /// particle position using the interpolating functions defined in
  /// MeshInterpF.ChF
  void interpolateParticle(RealVect& a_particleField,
                           const FArrayBox& a_field,
                           const RealVect& a_domainLeftEdge,
                           const RealVect& a_dx,
                           const RealVect& a_position,
                           const InterpType& a_interpType);
  
   //
   // JRA
   //
   /*
   std::array<Real,3> interpolateToParticle(
                                  const FArrayBox&   a_Ex,
                                  const FArrayBox&   a_Ey,
                                  const FArrayBox&   a_Ez,
                                  const RealVect&    a_domainLeftEdge,
                                  const RealVect&    a_dx,
                                  const RealVect&    a_position,
                                  const InterpType&  a_interpType ) const;
   */
   Real interpolateToParticle(
                         const FArrayBox&   a_field,
                         const RealVect&    a_domainLeftEdge,
                         const RealVect&    a_dx,
                         const RealVect&    a_position,
                         const IntVect&     a_stag,
                         const InterpType&  a_interpType,
                         const int          a_comp ) const;
  
  Box m_domain;
  int m_ghosts;
  RealVect m_dx;
  RealVect m_domainLeftEdge;
  RealVect m_domainRightEdge;
  IntVect m_bc_check_lo;
  IntVect m_bc_check_hi;
};

#include "NamespaceFooter.H"

// implementation
#include "MeshInterpI.H"

#endif // include guard
