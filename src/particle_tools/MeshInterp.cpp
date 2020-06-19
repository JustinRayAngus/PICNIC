#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MeshInterp.H"
#include "NamespaceHeader.H"

MeshInterp::MeshInterp()
{}

MeshInterp::MeshInterp(const Box&  a_domain,
                       const RealVect& a_dx,
                       const RealVect& a_domainLeftEdge)
{
  m_domain = a_domain;
  m_domainLeftEdge = a_domainLeftEdge;
  m_dx = a_dx;
}

void MeshInterp::define(const Box&  a_domain,
                        const RealVect& a_dx,
                        const RealVect& a_domainLeftEdge)
{
  m_domain = a_domain;
  m_domainLeftEdge = a_domainLeftEdge;
  m_dx = a_dx;
}

void MeshInterp::depositParticle(FArrayBox& a_rho,
                                 const RealVect& a_domainLeftEdge,
                                 const RealVect& a_dx,
                                 const RealVect& a_position,
                                 const Real&     a_mass,
                                 const InterpType a_interpType) const // JRA const

{
  switch (a_interpType)
    {
    case NGP:
      FORT_NGP_DEPOSIT(CHF_FRA1(a_rho, 0),
                       CHF_CONST_REALVECT(a_domainLeftEdge),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_position),
                       CHF_CONST_REAL(a_mass));
      break;
    case CIC:
      FORT_CIC_DEPOSIT(CHF_FRA1(a_rho, 0),
                       CHF_CONST_REALVECT(a_domainLeftEdge),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_position),
                       CHF_CONST_REAL(a_mass));
      break;
    case TSC:
      FORT_TSC_DEPOSIT(CHF_FRA1(a_rho, 0),
                       CHF_CONST_REALVECT(a_domainLeftEdge),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_position),
                       CHF_CONST_REAL(a_mass));
      break;
    case W4:
      FORT_W4_DEPOSIT(CHF_FRA1(a_rho, 0),
                      CHF_CONST_REALVECT(a_domainLeftEdge),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_position),
                      CHF_CONST_REAL(a_mass));
      break;
    default:
      MayDay::Error("Invalid interpolation type in MeshInterp::depositParticle");
    }
}

void MeshInterp::interpolateParticle(RealVect& a_particleField,
                                     const FArrayBox& a_field,
                                     const RealVect& a_domainLeftEdge,
                                     const RealVect& a_dx,
                                     const RealVect& a_position,
                                     const InterpType& a_interpType)
{
  switch (a_interpType)
    {
    case NGP:
      FORT_NGP_INTERPOLATE(CHF_REALVECT(a_particleField),
                           CHF_CONST_FRA(a_field),
                           CHF_CONST_REALVECT(a_domainLeftEdge),
                           CHF_CONST_REALVECT(a_dx),
                           CHF_CONST_REALVECT(a_position));

      break;
    case CIC:
      FORT_CIC_INTERPOLATE(CHF_REALVECT(a_particleField),
                           CHF_CONST_FRA(a_field),
                           CHF_CONST_REALVECT(a_domainLeftEdge),
                           CHF_CONST_REALVECT(a_dx),
                           CHF_CONST_REALVECT(a_position));
      break;
    case TSC:
      FORT_TSC_INTERPOLATE(CHF_REALVECT(a_particleField),
                           CHF_CONST_FRA(a_field),
                           CHF_CONST_REALVECT(a_domainLeftEdge),
                           CHF_CONST_REALVECT(a_dx),
                           CHF_CONST_REALVECT(a_position));
      break;
    case W4:
      FORT_W4_INTERPOLATE(CHF_REALVECT(a_particleField),
                          CHF_CONST_FRA(a_field),
                          CHF_CONST_REALVECT(a_domainLeftEdge),
                          CHF_CONST_REALVECT(a_dx),
                          CHF_CONST_REALVECT(a_position));
      break;
    default:
      MayDay::Error("Invalid interpolation type in MeshInterp::interpolateParticle.");
    }
}

//
//
//

void MeshInterp::momentParticle( FArrayBox&  a_moment,
                           const RealVect&   a_domainLeftEdge,
                           const RealVect&   a_dx,
                           const RealVect&   a_position,
                           const RealVect&   a_velocity,
                           const Real&       a_weight, 
                           const Real&       a_species_mass, 
                           const MomentType  a_momentType ) const // JRA const

{ 
  Real kernal;
  switch (a_momentType)
    {
    case density:
      kernal = 1.0;
      FORT_MOMENT_DEPOSIT( CHF_FRA1(a_moment, 0),
                   CHF_CONST_REALVECT(a_domainLeftEdge),
                   CHF_CONST_REALVECT(a_dx),
                   CHF_CONST_REALVECT(a_position),
                   CHF_CONST_REAL(kernal),
                   CHF_CONST_REAL(a_weight) );
      break;
    case momentum:
      // CH_assert(a_moment.nComp()==SpaceDim); // done at higher level
      for(int dir=0; dir<SpaceDim; dir++) {
         kernal = a_species_mass*a_velocity[dir];
         FORT_MOMENT_DEPOSIT( CHF_FRA1(a_moment, dir),
                      CHF_CONST_REALVECT(a_domainLeftEdge),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_position),
                      CHF_CONST_REAL(kernal),
                      CHF_CONST_REAL(a_weight) );
      }
      break;
    case energy:
      kernal = 0.0;
      for(int dir=0; dir<SpaceDim; dir++) {
         kernal = kernal + a_velocity[dir]*a_velocity[dir];
      }
      kernal = a_species_mass*kernal/2.0;
      FORT_MOMENT_DEPOSIT( CHF_FRA1(a_moment, 0),
                   CHF_CONST_REALVECT(a_domainLeftEdge),
                   CHF_CONST_REALVECT(a_dx),
                   CHF_CONST_REALVECT(a_position),
                   CHF_CONST_REAL(kernal),
                   CHF_CONST_REAL(a_weight) );
      break;
    case heatFlux:
      MayDay::Error("heatFlux type in MeshInterp::momentParticle not defined yet");
      break;
    default:
      MayDay::Error("Invalid moment type in MeshInterp::momentParticle");
    }
}

#include "NamespaceFooter.H"
