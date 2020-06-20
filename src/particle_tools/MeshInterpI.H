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
void MeshInterp::deposit(const List<P>& a_particleList,
                         FArrayBox& a_rho,
                         InterpType& a_interpType) const // JRA const
{
  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
    const P& this_particle = lit();
    const RealVect& this_x = this_particle.position();
    //const Real& this_m = this_particle.mass();
    const Real& this_m = this_particle.weight(); // JRA
    //const Real& this_pos_virt = this_particle.pos_virt();
    //if(!procID()) cout << "this_pos_virt = " << this_pos_virt << endl;
    //if(!procID()) cout << "this_pos_virt = " << this_particle.pos_virt() << endl;
    depositParticle(a_rho, 
                    m_domainLeftEdge, 
                    m_dx, 
                    this_x, 
                    this_m, 
                    a_interpType);
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
void MeshInterp::moment( FArrayBox&  a_moment,
                   const List<P>&    a_particleList,
                   const Real&       a_species_mass,
                   const MomentType  a_momentType ) const // JRA const
{

  if(a_momentType==momentum) CH_assert(a_moment.nComp()==SpaceDim);

  for (ListIterator<P> lit(a_particleList); lit; ++lit)
  {
     const P& this_particle = lit();
     const RealVect& this_x = this_particle.position();
     const RealVect& this_v = this_particle.velocity();
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
