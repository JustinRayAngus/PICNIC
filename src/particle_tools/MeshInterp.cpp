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
                       const int  a_ghosts,
                       const RealVect& a_dx,
                       const RealVect& a_domainLeftEdge,
                       const RealVect& a_domainRightEdge)
{
  m_domain = a_domain;
  m_domainLeftEdge = a_domainLeftEdge;
  m_domainRightEdge = a_domainRightEdge;
  m_ghosts = a_ghosts;
  m_dx = a_dx;
  m_bc_check_lo = IntVect::Zero;
  m_bc_check_hi = IntVect::Zero;
}

void MeshInterp::define(const Box&  a_domain,
                        const int  a_ghosts,
                        const RealVect& a_dx,
                        const RealVect& a_domainLeftEdge,
                        const RealVect& a_domainRightEdge)
{
  m_domain = a_domain;
  m_domainLeftEdge = a_domainLeftEdge;
  m_domainRightEdge = a_domainRightEdge;
  m_ghosts = a_ghosts;
  m_dx = a_dx;
  m_bc_check_lo = IntVect::Zero;
  m_bc_check_hi = IntVect::Zero;
}

void MeshInterp::depositParticle( FArrayBox&   a_rho,
                            const RealVect&    a_domainLeftEdge,
                            const RealVect&    a_dx,
                            const RealVect&    a_position,
                            const Real&        a_weight,
                            const Real&        a_kernal,
                            const IntVect&     a_stag,
                            const InterpType&  a_interpType,
                            const int          a_comp ) const // JRA const

{
  switch (a_interpType)
    {
    case NGP:
      FORT_NGP_DEPOSIT(CHF_FRA1(a_rho, 0),
                       CHF_CONST_REALVECT(a_domainLeftEdge),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_position),
                       CHF_CONST_REAL(a_weight));
      break;
    case CIC:
      FORT_CIC_DEPOSIT(CHF_FRA1(a_rho, a_comp),
                       CHF_CONST_REALVECT(a_domainLeftEdge),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_position),
                       CHF_CONST_INTVECT(a_stag),
                       CHF_CONST_REAL(a_kernal),
                       CHF_CONST_REAL(a_weight));
      break;
    case TSC:
      FORT_TSC_DEPOSIT(CHF_FRA1(a_rho, a_comp),
                       CHF_CONST_REALVECT(a_domainLeftEdge),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_CONST_REALVECT(a_position),
                       CHF_CONST_INTVECT(a_stag),
                       CHF_CONST_REAL(a_kernal),
                       CHF_CONST_REAL(a_weight));
      break;
    case W4:
      FORT_W4_DEPOSIT(CHF_FRA1(a_rho, 0),
                      CHF_CONST_REALVECT(a_domainLeftEdge),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_position),
                      CHF_CONST_REAL(a_weight));
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

Real MeshInterp::interpolateToParticle(
                                  const FArrayBox&   a_field,
                                  const RealVect&    a_domainLeftEdge,
                                  const RealVect&    a_dx,
                                  const RealVect&    a_position,
                                  const IntVect&     a_stag,
                                  const InterpType&  a_interpType,
                                  const int          a_comp ) const
{
   Real a_Fp_dir = 0.0;
   switch (a_interpType) {
      case CIC:
      FORT_CIC_INTERPOLATE_DIR( CHF_REAL(a_Fp_dir),
                                CHF_CONST_FRA1(a_field, a_comp),
                                CHF_CONST_REALVECT(a_domainLeftEdge),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_REALVECT(a_position), 
                                CHF_CONST_INTVECT(a_stag) );
      break;
      default:
      MayDay::Error("Invalid interpolation type in MeshInterp::interpolateToParticle.");
   }
   return a_Fp_dir;

}

//
//
//

void MeshInterp::momentParticle( FArrayBox&   a_moment,
                           const RealVect&    a_domainLeftEdge,
                           const RealVect&    a_dx,
                           const RealVect&    a_position,
                           const std::array<Real,3>&   a_velocity,
                           const Real&        a_weight, 
                           const Real&        a_species_mass, 
                           const MomentType&  a_momentType ) const // JRA const

{ 
  //CH_TIME("MeshInterp::momentParticle()"); // timer here seems to affect time spent here...
  
  Real kernal, gammap=1.0;
  Real kernal0, kernal1, kernal2;
  switch (a_momentType)
    {
    case number:
      FORT_COUNT_DEPOSIT( CHF_FRA1(a_moment, 0),
                   CHF_CONST_REALVECT(a_domainLeftEdge),
                   CHF_CONST_REALVECT(a_dx),
                   CHF_CONST_REALVECT(a_position) );
      break;
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
      kernal0 = a_species_mass*a_velocity[0];
      kernal1 = a_species_mass*a_velocity[1];
      kernal2 = a_species_mass*a_velocity[2];
      FORT_MOMENT3V_DEPOSIT( CHF_FRA(a_moment),
                      CHF_CONST_REALVECT(a_domainLeftEdge),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_position),
                      CHF_CONST_REAL(kernal0),
                      CHF_CONST_REAL(kernal1),
                      CHF_CONST_REAL(kernal2),
                      CHF_CONST_REAL(a_weight) );
      break;
    case energy:
#ifdef RELATIVISTIC_PARTICLES
      for (int n=0; n<3; n++) gammap += a_velocity[n]*a_velocity[n];
      gammap = sqrt(gammap);
#endif
      kernal = a_species_mass/(gammap + 1.0);
      kernal0 = kernal*a_velocity[0]*a_velocity[0];
      kernal1 = kernal*a_velocity[1]*a_velocity[1];
      kernal2 = kernal*a_velocity[2]*a_velocity[2];
      FORT_MOMENT3V_DEPOSIT( CHF_FRA(a_moment),
                   CHF_CONST_REALVECT(a_domainLeftEdge),
                   CHF_CONST_REALVECT(a_dx),
                   CHF_CONST_REALVECT(a_position),
                   CHF_CONST_REAL(kernal0),
                   CHF_CONST_REAL(kernal1),
                   CHF_CONST_REAL(kernal2),
                   CHF_CONST_REAL(a_weight) );
      break;
    case heatFlux:
      MayDay::Error("heatFlux type in MeshInterp::momentParticle not defined yet");
      break;
    default:
      MayDay::Error("Invalid moment type in MeshInterp::momentParticle");
    }
}

void MeshInterp::setParticleWeight( Real&       a_weight,
                              const Real&       a_partsPerCell,
                              const FArrayBox&  a_density,
                              const RealVect&   a_domainLeftEdge,
                              const RealVect&   a_dx,
                              const RealVect&   a_position ) const
{ 
   CH_TIME("MeshInterp::setParticleWeight()");

   FORT_SET_PARTICLE_WEIGHT( CHF_REAL(a_weight),
                             CHF_CONST_REAL(a_partsPerCell), 
                             CHF_CONST_FRA1(a_density,0),
                             CHF_CONST_REALVECT(a_domainLeftEdge), 
                             CHF_CONST_REALVECT(a_dx), 
                             CHF_CONST_REALVECT(a_position) );
}


#include "NamespaceFooter.H"
