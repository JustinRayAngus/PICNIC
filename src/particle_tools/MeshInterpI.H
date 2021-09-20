#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MESHINTERPI_H_
#define _MESHINTERPI_H_

#include "MeshInterp.H"

#include "NamespaceHeader.H"

template <class P>
void MeshInterp::deposit( const List<P>&  a_particleList,
                          FArrayBox&      a_rho,
                          const InterpType&     a_interpType ) const // JRA const
{
  CH_TIME("MeshInterp::deposit()");

  // user-specified way for interpolating particle density to 
  // grid points in a FArrayBox

  const Real this_kernal = 1.0; 
  const IntVect stag = a_rho.box().type();
  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_w = this_particle.weight();
    depositParticle(a_rho, 
                    m_domainLeftEdge, 
                    m_dx, 
                    this_x, 
                    this_w,
                    this_kernal,
                    stag, 
                    a_interpType);
  }

}

template <class P>
void MeshInterp::depositCurrent( FArrayBox&   a_rhoVx,
                                 FArrayBox&   a_rhoVy,
                                 FArrayBox&   a_rhoVz,
                           const List<P>&     a_particleList,
                           const InterpType&  a_interpType ) const
{
   CH_TIME("MeshInterp::depositCurrent()");

   // user-specified way for interpolating particle velocity density to 
   // grid points in a FArrayBox (used for current density calc)

   for (ListIterator<P> lit(a_particleList); lit; ++lit) {
      const P& this_particle = lit();
      const RealVect& this_x = this_particle.position();
      const std::array<Real,3>& this_vel = this_particle.velocity(); 
      const Real& this_w = this_particle.weight();
      
      switch (a_interpType) {
         case CIC:
         FORT_CIC_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, 0),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(this_x),
                                   CHF_CONST_REAL(this_vel[0]),
                                   CHF_CONST_REAL(this_vel[1]),
                                   CHF_CONST_REAL(this_vel[2]),
                                   CHF_CONST_REAL(this_w) );
         break;
         default:
         MayDay::Error("Invalid interpolation type in MeshInterp::interpolateToParticle.");
      }

   }

}

template <class P>
void MeshInterp::depositCurrent1D( FArrayBox&   a_rhoVx,
                                   FArrayBox&   a_rhoVvirt,
                             const List<P>&     a_particleList,
                             const InterpType&  a_interpType ) const
{
   //CH_TIME("MeshInterp::depositCurrent1D()");

   // user-specified way for interpolating particle velocity density to 
   // grid points in a FArrayBox (used for current density calc)

   for (ListIterator<P> lit(a_particleList); lit; ++lit) {
      const P& this_particle = lit();
      const RealVect& this_x = this_particle.position();
      const std::array<Real,3>& this_vel = this_particle.velocity(); 
      const Real& this_w = this_particle.weight();
      
      switch (a_interpType) {
         case CIC:
         FORT_CIC_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVvirt, 0),
                                   CHF_FRA1(a_rhoVvirt, 1),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(this_x),
                                   CHF_CONST_REAL(this_vel[0]),
                                   CHF_CONST_REAL(this_vel[1]),
                                   CHF_CONST_REAL(this_vel[2]),
                                   CHF_CONST_REAL(this_w) );
         break;
         default:
         MayDay::Error("Invalid interpolation type in MeshInterp::interpolateToParticle.");
      }

   }

}

template <class P>
void MeshInterp::depositVelocity( FArrayBox&   a_rhoV,
                            const int          a_dir,
                            const List<P>&     a_particleList,
                            const InterpType&  a_interpType,
                            const int          a_comp ) const
{
  CH_TIME("MeshInterp::depositVelocity()");

  CH_assert(a_dir<3);
  
  // user-specified way for interpolating particle velocity density to 
  // grid points in a FArrayBox (used for current density calc)

  const IntVect stag = a_rhoV.box().type();
  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    const Real& this_kernal = this_particle.velocity(a_dir); 
    const Real& this_w = this_particle.weight();
    depositParticle(a_rhoV, 
                    m_domainLeftEdge, 
                    m_dx, 
                    this_x, 
                    this_w,
                    this_kernal,
                    stag, 
                    a_interpType,
                    a_comp);
  }

}

template <class P>
void MeshInterp::interpolate(List<P>& a_particleList,
                             const FArrayBox& a_field,
                             InterpType& a_interpType)
{
  for (ListIterator<P> lit(a_particleList); lit; ++lit)
    {
      RealVect particleField(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
      P& this_particle = lit();
      RealVect& this_x = this_particle.position();
      interpolateParticle(particleField,
                          a_field,
                          m_domainLeftEdge,
                          m_dx,
                          this_x,
                          a_interpType);
      this_particle.setAcceleration(particleField);
    }
}

//
//  new functions by JRA
//

template <class P>
void MeshInterp::interpolateEMfieldsToPart( List<P>&     a_particleList,
                                      const FArrayBox&   a_Ex,
                                      const FArrayBox&   a_Ey,
                                      const FArrayBox&   a_Ez,
                                      const FArrayBox&   a_Bx,
                                      const FArrayBox&   a_By,
                                      const FArrayBox&   a_Bz,
                                      const InterpType&  a_interpType ) const
{
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     std::array<Real,3> this_Ep = {0.0,0.0,0.0};
     std::array<Real,3> this_Bp = {0.0,0.0,0.0};
     P& this_particle = lit();
     const RealVect& this_xp = this_particle.position();
     
     switch (a_interpType) {
        case CIC:
        FORT_CIC_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_Ex, 0),
                                     CHF_CONST_FRA1(a_Ey, 0),
                                     CHF_CONST_FRA1(a_Ez, 0),
                                     CHF_CONST_FRA1(a_Bx, 0),
                                     CHF_CONST_FRA1(a_By, 0),
                                     CHF_CONST_FRA1(a_Bz, 0),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xp) );
        break;
        default:
        MayDay::Error("Invalid interpolation type in MeshInterp::interpolateToParticle.");
     }
     this_particle.setElectricField(this_Ep);
     this_particle.setMagneticField(this_Bp);
  }

}

template <class P>
void MeshInterp::interpolateEMfieldsToPart1D( List<P>&     a_particleList,
                                        const FArrayBox&   a_Ex,
                                        const FArrayBox&   a_Evirt,
                                        const FArrayBox&   a_Bx,
                                        const FArrayBox&   a_Bvirt,
                                        const InterpType&  a_interpType ) const
{
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     std::array<Real,3> this_Ep = {0.0,0.0,0.0};
     std::array<Real,3> this_Bp = {0.0,0.0,0.0};
     P& this_particle = lit();
     const RealVect& this_xp = this_particle.position();
   
     switch (a_interpType) {
        case CIC:
        FORT_CIC_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_Ex, 0),
                                     CHF_CONST_FRA1(a_Evirt, 0),
                                     CHF_CONST_FRA1(a_Evirt, 1),
                                     CHF_CONST_FRA1(a_Bx, 0),
                                     CHF_CONST_FRA1(a_Bvirt, 0),
                                     CHF_CONST_FRA1(a_Bvirt, 1),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xp) );
        break;
        default:
        MayDay::Error("Invalid interpolation type in MeshInterp::interpolateToParticle.");
     }
     this_particle.setElectricField(this_Ep);
     this_particle.setMagneticField(this_Bp);
  }

}

template <class P>
void MeshInterp::interpolateEfieldToPart( List<P>&     a_particleList,
                                    const int          a_dir,
                                    const FArrayBox&   a_electric_field,
                                    const InterpType&  a_interpType,
                                    const int          a_comp ) const
{
  
  const IntVect stag = a_electric_field.box().type();
  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     Real this_Ep_dir = 0.0;
     P& this_particle = lit();
     RealVect& this_x = this_particle.position();
     this_Ep_dir = interpolateToParticle(
                            a_electric_field,
                            m_domainLeftEdge,
                            m_dx,
                            this_x,
                            stag, 
                            a_interpType,
                            a_comp);
     this_particle.setElectricField(this_Ep_dir,a_dir);
  }

}

template <class P>
void MeshInterp::interpolateBfieldToPart( List<P>&     a_particleList,
                                    const int          a_dir,
                                    const FArrayBox&   a_magnetic_field,
                                    const InterpType&  a_interpType,
                                    const int          a_comp ) const
{
  
  const IntVect stag = a_magnetic_field.box().type();
  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     Real this_Bp_dir = 0.0;
     P& this_particle = lit();
     RealVect& this_x = this_particle.position();
     this_Bp_dir = interpolateToParticle(
                            a_magnetic_field,
                            m_domainLeftEdge,
                            m_dx,
                            this_x,
                            stag, 
                            a_interpType,
                            a_comp);
     this_particle.setMagneticField(this_Bp_dir,a_dir);
  }

}

template <class P>
void MeshInterp::moment( FArrayBox&  a_moment,
                   const List<P>&    a_particleList,
                   const Real&       a_species_mass,
                   const MomentType  a_momentType ) const // JRA const
{
  //CH_TIME("MeshInterp::moment()");

  if(a_momentType==momentum) CH_assert(a_moment.nComp()==3);
  if(a_momentType==energy) CH_assert(a_moment.nComp()==3);

  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
     const P& this_particle = lit();
     const RealVect& this_x = this_particle.position();
     const std::array<Real,3>& this_v = this_particle.velocity();
     const Real& this_w = this_particle.weight();
     momentParticle( a_moment, 
                     m_domainLeftEdge, 
                     m_dx, 
                     this_x, 
                     this_v, 
                     this_w,
                     a_species_mass,
                     a_momentType );
   }
}

template <class P>
void MeshInterp::setWeightFromGridProfile( List<P>&    a_particleList,
                                     const FArrayBox&  a_density,
                                     const Real        a_partsPerCell ) const
{

  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
     P& this_particle = lit();
     Real& this_w = this_particle.weight();
     const RealVect& this_x = this_particle.position();
     setParticleWeight( this_w,
                        a_partsPerCell, 
                        a_density,
                        m_domainLeftEdge, 
                        m_dx, 
                        this_x );
   }
}


#include "NamespaceFooter.H"

#endif // include guard
