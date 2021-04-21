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
     //interpolateToParticle( this_Ep_dir,
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
     //interpolateToParticle( this_Bp_dir,
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
