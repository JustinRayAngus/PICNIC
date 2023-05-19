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
#include <array>

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
                           const InterpType&  a_interpType,
                           const int          a_cyl_car,
                           const bool         a_from_explicit_solver ) const
{
   CH_TIME("MeshInterp::depositCurrent()");

   // user-specified way for interpolating particle velocity density to 
   // grid points in a FArrayBox (used for current density calc)
    
   const int z_comp = a_rhoVz.nComp()-1;

   for (ListIterator<P> lit(a_particleList); lit; ++lit) {
      
      const RealVect& xpold = lit().position_old();
      const RealVect& xp = lit().position();
      const std::array<Real,3>& up = lit().velocity(); 
      const Real& wp = lit().weight();
   
#ifdef RELATIVISTIC_PARTICLES
      Real gammap = 1.0;
      if(a_from_explicit_solver) {
         gammap += up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
         gammap = sqrt(gammap);
      }
      else {
         std::array<Real,3> upnew;
         Real gbsq_old, gbsq_new;
         const std::array<Real,3>& upold = lit().velocity_old();
         for(int n=0; n<3; n++) upnew[n] = 2.0*up[n] - upold[n];
         gbsq_old = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
         gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
         gammap = 0.5*(sqrt(1.0 + gbsq_old) + sqrt(1.0 + gbsq_new));
      }
      const Real wpog = wp/gammap;
#else
      const Real wpog = wp;
#endif

      std::array<Real,3> up0 = up;
      if(a_cyl_car) { // transform up from CAR to CYL
         int dirpr = 0;
         int dirpth = 1;
#if CH_SPACEDIM==2
         dirpth = 2;
#endif
         const Real costhp = lit().position_virt(0);
         const Real sinthp = lit().position_virt(1);
         const Real upx = up[dirpr];
         const Real upy = up[dirpth];
         up0[dirpr]  =  costhp*upx + sinthp*upy; 
         up0[dirpth] = -sinthp*upx + costhp*upy;
      }
      
      switch (a_interpType) {
         case CIC:
         FORT_CIC_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, z_comp),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(xp),
                                   CHF_CONST_REAL(up0[0]),
                                   CHF_CONST_REAL(up0[1]),
                                   CHF_CONST_REAL(up0[2]),
                                   CHF_CONST_REAL(wpog) );
         break;
         case TSC:
         FORT_TSC_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, z_comp),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(xp),
                                   CHF_CONST_REAL(up0[0]),
                                   CHF_CONST_REAL(up0[1]),
                                   CHF_CONST_REAL(up0[2]),
                                   CHF_CONST_REAL(wpog) );
         break;
         case CC0: // zero-order charge conserving scheme
#if CH_SPACEDIM == 1
         FORT_CC0_1D_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, z_comp),
                                   CHF_CONST_INTVECT(m_bc_check_lo),
                                   CHF_CONST_INTVECT(m_bc_check_hi),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_domainRightEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(xpold),
                                   CHF_CONST_REALVECT(xp),
                                   CHF_CONST_REAL(up0[0]),
                                   CHF_CONST_REAL(up0[1]),
                                   CHF_CONST_REAL(up0[2]),
                                   CHF_CONST_REAL(wpog) );
#endif
#if CH_SPACEDIM == 2
         FORT_CC0_2D_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, z_comp),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(xpold),
                                   CHF_CONST_REALVECT(xp),
                                   CHF_CONST_REAL(up0[0]),
                                   CHF_CONST_REAL(up0[1]),
                                   CHF_CONST_REAL(up0[2]),
                                   CHF_CONST_REAL(wpog) );
#endif
         break;
         case CC1: // first-order charge conserving scheme
#if CH_SPACEDIM == 1
         FORT_CC1_1D_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, z_comp),
                                   CHF_CONST_INT(m_ghosts),
                                   CHF_CONST_INTVECT(m_bc_check_lo),
                                   CHF_CONST_INTVECT(m_bc_check_hi),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_domainRightEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(xpold),
                                   CHF_CONST_REALVECT(xp),
                                   CHF_CONST_REAL(up0[0]),
                                   CHF_CONST_REAL(up0[1]),
                                   CHF_CONST_REAL(up0[2]),
                                   CHF_CONST_REAL(wpog) );
#endif
#if CH_SPACEDIM == 2
         FORT_CC1_2D_DEPOSIT_CURRENT( CHF_FRA1(a_rhoVx, 0),
                                   CHF_FRA1(a_rhoVy, 0),
                                   CHF_FRA1(a_rhoVz, z_comp),
                                   CHF_CONST_INT(m_ghosts),
                                   CHF_CONST_REALVECT(m_domainLeftEdge),
                                   CHF_CONST_REALVECT(m_dx),
                                   CHF_CONST_REALVECT(xpold),
                                   CHF_CONST_REALVECT(xp),
                                   CHF_CONST_REAL(up0[0]),
                                   CHF_CONST_REAL(up0[1]),
                                   CHF_CONST_REAL(up0[2]),
                                   CHF_CONST_REAL(wpog) );
#endif
         break;
         default:
         MayDay::Error("Invalid interpolation type in MeshInterp::depositCurrent()");
      }

   }

}

template <class P>
void MeshInterp::depositMassMatrices( FArrayBox&   a_sigma_xx,
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
                                const int          a_inert_type,
                                const List<P>&     a_particleList,
                                const InterpType&  a_interpType ) const
{
   const int z_comp = a_Jz.nComp()-1; 
   for (ListIterator<P> lit(a_particleList); lit; ++lit) {

      const RealVect& xpbar = lit().position();
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upbar = lit().velocity(); 
      const std::array<Real,3>& upold = lit().velocity_old(); 
      const Real wp = lit().weight()*a_qovs;

      std::array<Real,3> upbar0 = upbar;
      Real inert_fact = 0.0;
      if(a_inert_type==1) {
         inert_fact = lit().position_virt(0); // dtheta
#if CH_SPACEDIM==2
         cout << "JRA: the anti-cyclic nature for cyl RZ has not yet been";
         cout << " accounted for in the mass matrix deposits" << endl;
#endif
      }
      else if(a_inert_type==2) {
         const Real costhp = lit().position_virt(0);
         const Real sinthp = lit().position_virt(1);
         inert_fact = costhp*costhp + sinthp*sinthp; // not equal to unity in general
         int dirpr = 0;
         int dirpth = 1;
#if CH_SPACEDIM==2
         dirpth = 2;
         cout << "JRA: the anti-cyclic nature for cyl RZ has not yet been";
         cout << " accounted for in the mass matrix deposits" << endl;
#endif
         const Real upx = upbar[dirpr];
         const Real upy = upbar[dirpth];
         upbar0[dirpr]  =  costhp*upx + sinthp*upy;
         upbar0[dirpth] = -sinthp*upx + costhp*upy;
      }
      
      switch (a_interpType) {
         case CIC:
         FORT_CIC_DEPOSIT_MASS_MATRIX( CHF_FRA(a_sigma_xx),
                                       CHF_FRA(a_sigma_xy),
                                       CHF_FRA(a_sigma_xz),
                                       CHF_FRA(a_sigma_yx),
                                       CHF_FRA(a_sigma_yy),
                                       CHF_FRA(a_sigma_yz),
                                       CHF_FRA(a_sigma_zx),
                                       CHF_FRA(a_sigma_zy),
                                       CHF_FRA(a_sigma_zz),
                                       CHF_FRA1(a_Jx, 0),
                                       CHF_FRA1(a_Jy, 0),
                                       CHF_FRA1(a_Jz, z_comp),
                                 CHF_CONST_FRA1(a_Bx, 0),
                                 CHF_CONST_FRA1(a_By, 0),
                                 CHF_CONST_FRA1(a_Bz, z_comp),
                                 CHF_CONST_REALVECT(m_domainLeftEdge),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REALVECT(xpold),
                                 CHF_CONST_REALVECT(xpbar),
                                 CHF_CONST_REAL(upold[0]),
                                 CHF_CONST_REAL(upold[1]),
                                 CHF_CONST_REAL(upold[2]),
                                 CHF_CONST_REAL(upbar0[0]),
                                 CHF_CONST_REAL(upbar0[1]),
                                 CHF_CONST_REAL(upbar0[2]),
                                 CHF_CONST_REAL(a_alphas),
                                 CHF_CONST_INT(a_inert_type),
                                 CHF_CONST_REAL(inert_fact),
                                 CHF_CONST_REAL(wp) );
         break;
         case TSC:
         FORT_TSC_DEPOSIT_MASS_MATRIX( CHF_FRA(a_sigma_xx),
                                       CHF_FRA(a_sigma_xy),
                                       CHF_FRA(a_sigma_xz),
                                       CHF_FRA(a_sigma_yx),
                                       CHF_FRA(a_sigma_yy),
                                       CHF_FRA(a_sigma_yz),
                                       CHF_FRA(a_sigma_zx),
                                       CHF_FRA(a_sigma_zy),
                                       CHF_FRA(a_sigma_zz),
                                       CHF_FRA1(a_Jx, 0),
                                       CHF_FRA1(a_Jy, 0),
                                       CHF_FRA1(a_Jz, z_comp),
                                 CHF_CONST_FRA1(a_Bx, 0),
                                 CHF_CONST_FRA1(a_By, 0),
                                 CHF_CONST_FRA1(a_Bz, z_comp),
                                 CHF_CONST_REALVECT(m_domainLeftEdge),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REALVECT(xpold),
                                 CHF_CONST_REALVECT(xpbar),
                                 CHF_CONST_REAL(upold[0]),
                                 CHF_CONST_REAL(upold[1]),
                                 CHF_CONST_REAL(upold[2]),
                                 CHF_CONST_REAL(upbar0[0]),
                                 CHF_CONST_REAL(upbar0[1]),
                                 CHF_CONST_REAL(upbar0[2]),
                                 CHF_CONST_REAL(a_alphas),
                                 CHF_CONST_INT(a_inert_type),
                                 CHF_CONST_REAL(inert_fact),
                                 CHF_CONST_REAL(wp) );
         break;
         case CC0:
#if CH_SPACEDIM == 1
         FORT_CC0_1D_DEPOSIT_MASS_MATRIX( CHF_FRA(a_sigma_xx),
                                       CHF_FRA(a_sigma_xy),
                                       CHF_FRA(a_sigma_xz),
                                       CHF_FRA(a_sigma_yx),
                                       CHF_FRA(a_sigma_yy),
                                       CHF_FRA(a_sigma_yz),
                                       CHF_FRA(a_sigma_zx),
                                       CHF_FRA(a_sigma_zy),
                                       CHF_FRA(a_sigma_zz),
                                       CHF_FRA1(a_Jx, 0),
                                       CHF_FRA1(a_Jy, 0),
                                       CHF_FRA1(a_Jz, z_comp),
                                 CHF_CONST_FRA1(a_Bx, 0),
                                 CHF_CONST_FRA1(a_By, 0),
                                 CHF_CONST_FRA1(a_Bz, z_comp),
                                 CHF_CONST_INTVECT(m_bc_check_lo),
                                 CHF_CONST_INTVECT(m_bc_check_hi),
                                 CHF_CONST_REALVECT(m_domainLeftEdge),
                                 CHF_CONST_REALVECT(m_domainRightEdge),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REALVECT(xpold),
                                 CHF_CONST_REALVECT(xpbar),
                                 CHF_CONST_REAL(upold[0]),
                                 CHF_CONST_REAL(upold[1]),
                                 CHF_CONST_REAL(upold[2]),
                                 CHF_CONST_REAL(upbar0[0]),
                                 CHF_CONST_REAL(upbar0[1]),
                                 CHF_CONST_REAL(upbar0[2]),
                                 CHF_CONST_REAL(a_alphas),
                                 CHF_CONST_INT(a_inert_type),
                                 CHF_CONST_REAL(inert_fact),
                                 CHF_CONST_REAL(wp) );
#endif
         break;
         case CC1:
#if CH_SPACEDIM == 1
         FORT_CC1_1D_DEPOSIT_MASS_MATRIX( CHF_FRA(a_sigma_xx),
                                       CHF_FRA(a_sigma_xy),
                                       CHF_FRA(a_sigma_xz),
                                       CHF_FRA(a_sigma_yx),
                                       CHF_FRA(a_sigma_yy),
                                       CHF_FRA(a_sigma_yz),
                                       CHF_FRA(a_sigma_zx),
                                       CHF_FRA(a_sigma_zy),
                                       CHF_FRA(a_sigma_zz),
                                       CHF_FRA1(a_Jx, 0),
                                       CHF_FRA1(a_Jy, 0),
                                       CHF_FRA1(a_Jz, z_comp),
                                 CHF_CONST_FRA1(a_Bx, 0),
                                 CHF_CONST_FRA1(a_By, 0),
                                 CHF_CONST_FRA1(a_Bz, z_comp),
                                 CHF_CONST_INTVECT(m_bc_check_lo),
                                 CHF_CONST_INTVECT(m_bc_check_hi),
                                 CHF_CONST_REALVECT(m_domainLeftEdge),
                                 CHF_CONST_REALVECT(m_domainRightEdge),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REALVECT(xpold),
                                 CHF_CONST_REALVECT(xpbar),
                                 CHF_CONST_REAL(upold[0]),
                                 CHF_CONST_REAL(upold[1]),
                                 CHF_CONST_REAL(upold[2]),
                                 CHF_CONST_REAL(upbar0[0]),
                                 CHF_CONST_REAL(upbar0[1]),
                                 CHF_CONST_REAL(upbar0[2]),
                                 CHF_CONST_REAL(a_alphas),
                                 CHF_CONST_INT(a_inert_type),
                                 CHF_CONST_REAL(inert_fact),
                                 CHF_CONST_REAL(wp) );
#endif
#if CH_SPACEDIM == 2
         FORT_CC1_2D_DEPOSIT_MASS_MATRIX( CHF_FRA(a_sigma_xx),
                                       CHF_FRA(a_sigma_xy),
                                       CHF_FRA(a_sigma_xz),
                                       CHF_FRA(a_sigma_yx),
                                       CHF_FRA(a_sigma_yy),
                                       CHF_FRA(a_sigma_yz),
                                       CHF_FRA(a_sigma_zx),
                                       CHF_FRA(a_sigma_zy),
                                       CHF_FRA(a_sigma_zz),
                                       CHF_FRA1(a_Jx, 0),
                                       CHF_FRA1(a_Jy, 0),
                                       CHF_FRA1(a_Jz, z_comp),
                                 CHF_CONST_FRA1(a_Bx, 0),
                                 CHF_CONST_FRA1(a_By, 0),
                                 CHF_CONST_FRA1(a_Bz, z_comp),
                                 CHF_CONST_REALVECT(m_domainLeftEdge),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REALVECT(xpold),
                                 CHF_CONST_REALVECT(xpbar),
                                 CHF_CONST_REAL(upold[0]),
                                 CHF_CONST_REAL(upold[1]),
                                 CHF_CONST_REAL(upold[2]),
                                 CHF_CONST_REAL(upbar0[0]),
                                 CHF_CONST_REAL(upbar0[1]),
                                 CHF_CONST_REAL(upbar0[2]),
                                 CHF_CONST_REAL(a_alphas),
                                 CHF_CONST_INT(a_inert_type),
                                 CHF_CONST_REAL(inert_fact),
                                 CHF_CONST_REAL(wp) );
#endif
         break;
         default:
         MayDay::Error("Invalid interpolation type in MeshInterp::depositMassMatrices()");
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
                                      const FArrayBox&   a_E0,
                                      const FArrayBox&   a_E1,
                                      const FArrayBox&   a_E2,
                                      const FArrayBox&   a_B0,
                                      const FArrayBox&   a_B1,
                                      const FArrayBox&   a_B2,
                                      const InterpType&  a_interpType,
                                      const int          a_blank_B ) const
{
 
  const int z_comp = a_E2.nComp()-1;
 
  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     std::array<Real,3> this_Ep = {0.0,0.0,0.0};
     std::array<Real,3> this_Bp = {0.0,0.0,0.0};
     P& this_particle = lit();
     const RealVect& this_xpold = this_particle.position_old();
     const RealVect& this_xp = this_particle.position();
     
     switch (a_interpType) {
        case CIC:
        FORT_CIC_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xp) );
        break;
        case TSC:
        FORT_TSC_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xp) );
        break;
        case CC0:
#if CH_SPACEDIM == 1
        FORT_CC0_1D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_INTVECT(m_bc_check_lo),
                                     CHF_CONST_INTVECT(m_bc_check_hi),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_domainRightEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
#if CH_SPACEDIM == 2
        FORT_CC0_2D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
        break;
        case CC1:
#if CH_SPACEDIM == 1
        FORT_CC1_1D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_INT(m_ghosts),
                                     CHF_CONST_INTVECT(m_bc_check_lo),
                                     CHF_CONST_INTVECT(m_bc_check_hi),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_domainRightEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
#if CH_SPACEDIM == 2
        FORT_CC1_2D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_INT(m_ghosts),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
        break;
        default:
        MayDay::Error("Invalid interpolation type in MeshInterp::interpolateEMfieldsToPart()");
     }
     this_particle.setElectricField(this_Ep);
     if(!a_blank_B) this_particle.setMagneticField(this_Bp);
  }

}

template <class P>
void MeshInterp::interpolateEMfieldsToPart( P&           a_particle,
                                      const FArrayBox&   a_E0,
                                      const FArrayBox&   a_E1,
                                      const FArrayBox&   a_E2,
                                      const FArrayBox&   a_B0,
                                      const FArrayBox&   a_B1,
                                      const FArrayBox&   a_B2,
                                      const InterpType&  a_interpType ) const
{
 
  const int z_comp = a_E2.nComp()-1;
 
  std::array<Real,3> this_Ep = {0.0,0.0,0.0};
  std::array<Real,3> this_Bp = {0.0,0.0,0.0};
  const RealVect& this_xpold = a_particle.position_old();
  const RealVect& this_xp = a_particle.position();
     
  switch (a_interpType) {
     case CIC:
     FORT_CIC_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                  CHF_REAL(this_Ep[1]),
                                  CHF_REAL(this_Ep[2]),
                                  CHF_REAL(this_Bp[0]),
                                  CHF_REAL(this_Bp[1]),
                                  CHF_REAL(this_Bp[2]),
                                  CHF_CONST_FRA1(a_E0, 0),
                                  CHF_CONST_FRA1(a_E1, 0),
                                  CHF_CONST_FRA1(a_E2, z_comp),
                                  CHF_CONST_FRA1(a_B0, 0),
                                  CHF_CONST_FRA1(a_B1, 0),
                                  CHF_CONST_FRA1(a_B2, z_comp),
                                  CHF_CONST_REALVECT(m_domainLeftEdge),
                                  CHF_CONST_REALVECT(m_dx),
                                  CHF_CONST_REALVECT(this_xp) );
     break;
     case TSC:
     FORT_TSC_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                  CHF_REAL(this_Ep[1]),
                                  CHF_REAL(this_Ep[2]),
                                  CHF_REAL(this_Bp[0]),
                                  CHF_REAL(this_Bp[1]),
                                  CHF_REAL(this_Bp[2]),
                                  CHF_CONST_FRA1(a_E0, 0),
                                  CHF_CONST_FRA1(a_E1, 0),
                                  CHF_CONST_FRA1(a_E2, z_comp),
                                  CHF_CONST_FRA1(a_B0, 0),
                                  CHF_CONST_FRA1(a_B1, 0),
                                  CHF_CONST_FRA1(a_B2, z_comp),
                                  CHF_CONST_REALVECT(m_domainLeftEdge),
                                  CHF_CONST_REALVECT(m_dx),
                                  CHF_CONST_REALVECT(this_xp) );
     break;
     case CC0:
#if CH_SPACEDIM == 1
     FORT_CC0_1D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_INTVECT(m_bc_check_lo),
                                     CHF_CONST_INTVECT(m_bc_check_hi),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_domainRightEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
#if CH_SPACEDIM == 2
     FORT_CC0_2D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
     break;
     case CC1:
#if CH_SPACEDIM == 1
     FORT_CC1_1D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_INT(m_ghosts),
                                     CHF_CONST_INTVECT(m_bc_check_lo),
                                     CHF_CONST_INTVECT(m_bc_check_hi),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_domainRightEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
#if CH_SPACEDIM == 2
     FORT_CC1_2D_INTERPOLATE_FIELDS( CHF_REAL(this_Ep[0]),
                                     CHF_REAL(this_Ep[1]),
                                     CHF_REAL(this_Ep[2]),
                                     CHF_REAL(this_Bp[0]),
                                     CHF_REAL(this_Bp[1]),
                                     CHF_REAL(this_Bp[2]),
                                     CHF_CONST_FRA1(a_E0, 0),
                                     CHF_CONST_FRA1(a_E1, 0),
                                     CHF_CONST_FRA1(a_E2, z_comp),
                                     CHF_CONST_FRA1(a_B0, 0),
                                     CHF_CONST_FRA1(a_B1, 0),
                                     CHF_CONST_FRA1(a_B2, z_comp),
                                     CHF_CONST_INT(m_ghosts),
                                     CHF_CONST_REALVECT(m_domainLeftEdge),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REALVECT(this_xpold),
                                     CHF_CONST_REALVECT(this_xp) );
#endif
      break;
      default:
      MayDay::Error("Invalid interpolation type in MeshInterp::interpolateEMfieldsToPart()");
   }
   a_particle.setElectricField(this_Ep);
   a_particle.setMagneticField(this_Bp);

}

template <class P>
void MeshInterp::getFieldDerivativesAtPart( RealVect&    a_dExdy,
                                            RealVect&    a_dEydy,
                                            RealVect&    a_dEzdy,
                                            RealVect&    a_dBxdy,
                                            RealVect&    a_dBydy,
                                            RealVect&    a_dBzdy,
                                      const P&           a_particle,
                                      const FArrayBox&   a_E0,
                                      const FArrayBox&   a_E1,
                                      const FArrayBox&   a_E2,
                                      const FArrayBox&   a_B0,
                                      const FArrayBox&   a_B1,
                                      const FArrayBox&   a_B2,
                                      const InterpType&  a_interpType ) const
{
 
  const int z_comp = a_E2.nComp()-1;
 
  const RealVect& this_xpold = a_particle.position_old();
  const RealVect& this_xp = a_particle.position();

  for (int dir=0; dir<SpaceDim; dir++) {
    a_dExdy[dir] = 0.0;
    a_dEydy[dir] = 0.0;
    a_dEzdy[dir] = 0.0;
    a_dBxdy[dir] = 0.0;
    a_dBydy[dir] = 0.0;
    a_dBzdy[dir] = 0.0;
  }
     
  switch (a_interpType) {
     case CIC:
     FORT_CIC_FIELD_DERIVATIVES( CHF_REALVECT(a_dExdy),
                                 CHF_REALVECT(a_dEydy),
                                 CHF_REALVECT(a_dEzdy),
                                 CHF_REALVECT(a_dBxdy),
                                 CHF_REALVECT(a_dBydy),
                                 CHF_REALVECT(a_dBzdy),
                                 CHF_CONST_FRA1(a_E0, 0),
                                 CHF_CONST_FRA1(a_E1, 0),
                                 CHF_CONST_FRA1(a_E2, z_comp),
                                 CHF_CONST_FRA1(a_B0, 0),
                                 CHF_CONST_FRA1(a_B1, 0),
                                 CHF_CONST_FRA1(a_B2, z_comp),
                                 CHF_CONST_REALVECT(m_domainLeftEdge),
                                 CHF_CONST_REALVECT(m_dx),
                                 CHF_CONST_REALVECT(this_xp) );
     break;
     //case TSC:
     //break;
     case CC0:
     FORT_CC0_1D_FIELD_DERIVATIVES( CHF_REALVECT(a_dExdy),
                                    CHF_REALVECT(a_dEydy),
                                    CHF_REALVECT(a_dEzdy),
                                    CHF_REALVECT(a_dBxdy),
                                    CHF_REALVECT(a_dBydy),
                                    CHF_REALVECT(a_dBzdy),
                                    CHF_CONST_FRA1(a_E0, 0),
                                    CHF_CONST_FRA1(a_E1, 0),
                                    CHF_CONST_FRA1(a_E2, z_comp),
                                    CHF_CONST_FRA1(a_B0, 0),
                                    CHF_CONST_FRA1(a_B1, 0),
                                    CHF_CONST_FRA1(a_B2, z_comp),
                                    CHF_CONST_INTVECT(m_bc_check_lo),
                                    CHF_CONST_INTVECT(m_bc_check_hi),
                                    CHF_CONST_REALVECT(m_domainLeftEdge),
                                    CHF_CONST_REALVECT(m_domainRightEdge),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REALVECT(this_xpold),
                                    CHF_CONST_REALVECT(this_xp) );
     break;
     case CC1:
#if CH_SPACEDIM == 1
     FORT_CC1_1D_FIELD_DERIVATIVES( CHF_REALVECT(a_dExdy),
                                    CHF_REALVECT(a_dEydy),
                                    CHF_REALVECT(a_dEzdy),
                                    CHF_REALVECT(a_dBxdy),
                                    CHF_REALVECT(a_dBydy),
                                    CHF_REALVECT(a_dBzdy),
                                    CHF_CONST_FRA1(a_E0, 0),
                                    CHF_CONST_FRA1(a_E1, 0),
                                    CHF_CONST_FRA1(a_E2, z_comp),
                                    CHF_CONST_FRA1(a_B0, 0),
                                    CHF_CONST_FRA1(a_B1, 0),
                                    CHF_CONST_FRA1(a_B2, z_comp),
                                    CHF_CONST_INTVECT(m_bc_check_lo),
                                    CHF_CONST_INTVECT(m_bc_check_hi),
                                    CHF_CONST_REALVECT(m_domainLeftEdge),
                                    CHF_CONST_REALVECT(m_domainRightEdge),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REALVECT(this_xpold),
                                    CHF_CONST_REALVECT(this_xp) );
#endif
#if CH_SPACEDIM == 2
#endif
      break;
      default:
      MayDay::Error("Invalid interpolation type in MeshInterp::getFieldDerivativesAtPart()");
   }

}

template <class P>
void MeshInterp::interpolateBfieldsToPart1D( List<P>&     a_particleList,
                                       const FArrayBox&   a_Bx,
                                       const FArrayBox&   a_Bvirt,
                                       const InterpType&  a_interpType ) const
{
  
  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     std::array<Real,3> this_Bp = {0.0,0.0,0.0};
     P& this_particle = lit();
     const RealVect& this_xp = this_particle.position();
   
     switch (a_interpType) {
        case CIC:
        FORT_CIC_INTERPOLATE_MAGNETIC_FIELD( CHF_REAL(this_Bp[0]),
                                             CHF_REAL(this_Bp[1]),
                                             CHF_REAL(this_Bp[2]),
                                             CHF_CONST_FRA1(a_Bx, 0),
                                             CHF_CONST_FRA1(a_Bvirt, 0),
                                             CHF_CONST_FRA1(a_Bvirt, 1),
                                             CHF_CONST_REALVECT(m_domainLeftEdge),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_REALVECT(this_xp) );
        break;
        default:
        MayDay::Error("Invalid interpolation type in MeshInterp::interpolateBfieldsToPart1D()");
     }
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

/////////////
//
//   non-fortran functions
//
/////////////

template <class P>
void MeshInterp::depositInflowCurrent( FArrayBox&   a_rhoVx,
                                       FArrayBox&   a_rhoVy,
                                       FArrayBox&   a_rhoVz,
                                 const List<P>&     a_particleList,
                                 const int          a_bdry_dir, 
                                 const int          a_bdry_side,
                                 const Real         a_cnormDt ) const
{
   CH_TIME("MeshInterp::depositInflowCurrent()");
   
   // inflow current is deposited via time-orbit-averged NGP for in-plane J
   // and time-orbit-averaged CIC for out-of-plane J

   // Note that this function is only called from the implicit solver.

#if CH_SPACEDIM == 1
   
   const int z_comp = a_rhoVz.nComp()-1;
   
   Real X0;
   if(a_bdry_side==0) X0 = m_domainLeftEdge[a_bdry_dir];
   else X0 = m_domainRightEdge[a_bdry_dir];

   IntVect index1, ii;
   RealVect xpnew, xpbar1;
   Real cnormDt0, cnormDt1; 
   Real volume = m_dx.product();
   Real w0, l0, orbit_fact;

   for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     const Real& wp = lit().weight();
     const RealVect& xpold = lit().position_old();
     const RealVect& xpbar = lit().position();
     const std::array<Real,3>& upold = lit().velocity_old(); 
     const std::array<Real,3>& upbar = lit().velocity(); 

     if(upbar[a_bdry_dir]==0.0) continue;
      
     Real gammap_old = 1.0;
     Real gammap = 1.0;
#ifdef RELATIVISTIC_PARTICLES
     std::array<Real,3> upnew;
     Real gbsq_old, gbsq_new;
     for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
     gbsq_old = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
     gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
     gammap_old = sqrt(1.0 + gbsq_old);
     gammap = 0.5*(gammap_old + sqrt(1.0+gbsq_new));
#endif

     cnormDt0 = (X0-xpold[a_bdry_dir])/(upold[a_bdry_dir]/gammap_old);
     cnormDt1 = a_cnormDt - cnormDt0;
     
     for (int dir=0; dir<SpaceDim; dir++) {
       xpnew[dir] = 2.0*xpbar[dir] - xpold[dir];
       if(dir==a_bdry_dir) xpbar1[dir] = (xpnew[dir] + X0)/2.0;
       else xpbar1[dir] = xpbar[dir];
       index1[dir] = floor((xpbar1[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
     }
     
     // this method assumes that the new position is no more than half a grid
     // spacing away from the boundary
     if(a_bdry_side==0) {CH_assert(xpnew[a_bdry_dir]<X0+m_dx[a_bdry_dir]/2.0);} 
     else {CH_assert(xpnew[a_bdry_dir]>X0-m_dx[a_bdry_dir]/2.0);}
      
     // orbit-averaged NGP for in-plane J
     orbit_fact = cnormDt1/a_cnormDt;
     a_rhoVx(index1,0) += orbit_fact*wp*upbar[0]/gammap/volume;

     // orbit-averaged CIC for out-of-plane J
     ii = index1;
     l0 = xpbar1[0] - (m_domainLeftEdge[0] + ii[0]*m_dx[0]);
     w0 = 1.0 - abs(l0/m_dx[0]);
     a_rhoVy(ii,0) += orbit_fact*w0*wp*upbar[1]/gammap/volume;
     a_rhoVz(ii,z_comp) += orbit_fact*w0*wp*upbar[2]/gammap/volume;
     
     w0 = 1.0 - w0;
     ii[0] += 1;
     a_rhoVy(ii,0) += orbit_fact*w0*wp*upbar[1]/gammap/volume;
     a_rhoVz(ii,z_comp) += orbit_fact*w0*wp*upbar[2]/gammap/volume;

   }
#endif
}

template <class P>
void MeshInterp::interpolateEMfieldsToInflowParts( List<P>&     a_particleList,
                                             const int          a_bdry_dir, 
                                             const int          a_bdry_side,
                                             const FArrayBox&   a_Ex,
                                             const FArrayBox&   a_Ey,
                                             const FArrayBox&   a_Ez,
                                             const FArrayBox&   a_Bx,
                                             const FArrayBox&   a_By,
                                             const FArrayBox&   a_Bz ) const
{
   CH_TIME("MeshInterp::interpolateEMfieldsToInflowParts()");

#if CH_SPACEDIM == 1
  
   const int z_comp = a_Ez.nComp()-1;
   
   Real X0;
   if(a_bdry_side==0) X0 = m_domainLeftEdge[a_bdry_dir];
   else X0 = m_domainRightEdge[a_bdry_dir];
  
   IntVect index1, ii;
   RealVect xpnew, xpbar1;
   std::array<Real,3> Ep, Bp;
   Real l0, w0;

   for (ListIterator<P> lit(a_particleList); lit; ++lit) {
     
     Ep = {0.0,0.0,0.0};
     Bp = {0.0,0.0,0.0};

     P& this_particle = lit();
     const RealVect& xpold = this_particle.position_old();
     const RealVect& xpbar = this_particle.position();
     
     for (int dir=0; dir<SpaceDim; dir++) {
       xpnew[dir] = 2.0*xpbar[dir] - xpold[dir];
       if(dir==a_bdry_dir) xpbar1[dir] = (xpnew[dir] + X0)/2.0;
       else xpbar1[dir] = xpbar[dir];
       index1[dir] = floor((xpbar1[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
     }

     // this method assumes that the new position is no more than half a grid
     // spacing away from the boundary
     if(a_bdry_side==0) {CH_assert(xpnew[a_bdry_dir]<X0+m_dx[a_bdry_dir]/2.0);} 
     else {CH_assert(xpnew[a_bdry_dir]>X0-m_dx[a_bdry_dir]/2.0);}
     
     // NGP for in-plane E
     Ep[0] += a_Ex(index1,0);
     
     // CIC for out-of-plane E and in-plane B
     ii = index1;
     l0 = xpbar1[0] - (m_domainLeftEdge[0] + ii[0]*m_dx[0]);
     w0 = 1.0 - abs(l0/m_dx[0]);
     Bp[0] += w0*a_Bx(ii,0);
     Ep[1] += w0*a_Ey(ii,0);
     Ep[2] += w0*a_Ez(ii,z_comp);
     
     w0 = 1.0 - w0;
     ii[0] += 1;
     Bp[0] += w0*a_Bx(ii,0);
     Ep[1] += w0*a_Ey(ii,0);
     Ep[2] += w0*a_Ez(ii,z_comp);
     
     // CIC for out-of-plane B
     ii = index1;
     l0 = xpbar1[0] - (m_domainLeftEdge[0] + ii[0]*m_dx[0]/2.0);
     w0 = 1.0 - abs(l0/m_dx[0]);
     Bp[1] += w0*a_By(ii,0);
     Bp[2] += w0*a_Bz(ii,z_comp);
     
     w0 = 1.0 - w0;
     ii = index1;
     ii[a_bdry_dir] += 2*a_bdry_side - 1;
     Bp[1] += w0*a_By(ii,0);
     Bp[2] += w0*a_Bz(ii,z_comp);
      
     this_particle.setElectricField(Ep);
     this_particle.setMagneticField(Bp);

   }
#endif

}

template <class P>
void MeshInterp::depositCurrent_testing( FArrayBox&   a_rhoVx,
                                         FArrayBox&   a_rhoVy,
                                         FArrayBox&   a_rhoVz,
                                   const List<P>&     a_particleList,
                                   const InterpType&  a_interpType ) const
{
   CH_TIME("MeshInterp::depositCurrent()");

   switch (a_interpType) {
      case CIC:
         //depositCurrent_CIC( a_particleList,
         //                    a_rhoVx, a_rhoVy, a_rhoVz );
      break;
      case TSC:
         //depositCurrent_TSC( a_particleList,
         //                    a_rhoVx, a_rhoVy, a_rhoVz );
      break;
      case CC0: // zero-order charge conserving scheme
         CH_assert(SpaceDim<3);
         ////  in-plane J computed using charge-conserving method
         //if(SpaceDim==1) depositCurrent_CC0( a_particleList, a_rhoVx );
         //if(SpaceDim==2) depositCurrent_CC0( a_particleList, a_rhoVx, a_rhoVy );
      break;
      case CC1: // first-order charge conserving scheme
         CH_assert(SpaceDim<3);
         ////  in-plane J computed using charge-conserving method
         //if(SpaceDim==1) depositCurrent_CC1( a_particleList, a_rhoVx );
         //if(SpaceDim==2) depositCurrent_CC1( a_particleList, a_rhoVx, a_rhoVy );
      break;
      default:
      MayDay::Error("Invalid interpolation type in MeshInterp::depositCurrent_testing()");
   }

}

template <class P>
void MeshInterp::interpolateEMfieldsToPart_testing( List<P>&     a_particleList,
                                        const FArrayBox&   a_Ex,
                                        const FArrayBox&   a_Ey,
                                        const FArrayBox&   a_Ez,
                                        const FArrayBox&   a_Bx,
                                        const FArrayBox&   a_By,
                                        const FArrayBox&   a_Bz,
                                        const InterpType&  a_interpType ) const
{
  CH_TIME("MeshInterp::interpolateEMfieldsToPart_testing()");

   switch (a_interpType) {
      case CIC:
         interpolateEMfieldsToPart_CIC( a_particleList,
                                        a_Ex, a_Ey, a_Ez,
                                        a_Bx, a_By, a_Bz );
      break;
      case TSC:
         interpolateEMfieldsToPart_TSC( a_particleList,
                                        a_Ex, a_Ey, a_Ez,
                                        a_Bx, a_By, a_Bz );
      break;
      case CC0:
         CH_assert(SpaceDim<3);
         if(SpaceDim==1) {
            interpolateEMfieldsToPart_CIC( a_particleList,
                                           a_Ex, a_Ey, a_Ez,
                                           a_Bx, a_By, a_Bz );
            interpolateEToPart_CC0( a_particleList, a_Ex );
         }
         if(SpaceDim==2) {
            interpolateBfieldsToPart_CIC( a_particleList,
                                          a_Bx, a_By, a_Bz );
            //interpolateEfieldDirToPart_CIC( a_particleList, a_Ez, 2, 0 );
            interpolateEToPart_CC0( a_particleList, a_Ex, a_Ey, a_Ez );
         }
      break;
      case CC1:
         CH_assert(SpaceDim<3);
         interpolateEMfieldsToPart_CIC( a_particleList,
                                        a_Ex, a_Ey, a_Ez,
                                        a_Bx, a_By, a_Bz );
         if(SpaceDim==1) interpolateEToPart_CC1( a_particleList, a_Ex );
         //if(SpaceDim==2) interpolateEToPart_CC1( a_particleList, a_Ex, a_Ey );
      break;
      default:
      MayDay::Error("Invalid interpolation type in MeshInterp::interpolateEMfieldsToPart_testing");
   }

}

template <class P>
void MeshInterp::interpolateEMfieldsToPart_CIC( List<P>&     a_particleList,
                                          const FArrayBox&   a_Ex,
                                          const FArrayBox&   a_Ey,
                                          const FArrayBox&   a_Ez,
                                          const FArrayBox&   a_Bx,
                                          const FArrayBox&   a_By,
                                          const FArrayBox&   a_Bz ) const
{
  CH_TIME("MeshInterp::interpolateEMfieldsToPart_CIC()");
 
  std::array<Real,3> Ep, Bp;
  IntVect index, index_stag, ii, ii_stag; 
  IntVect ii_Ex, ii_Ey, ii_Ez, ii_Bx, ii_By, ii_Bz; 
  RealVect w0_Ex, w0_Ey, w0_Ez, w0_Bx, w0_By, w0_Bz; 
  RealVect delta, delta_stag, w0, w0_stag;

  const int z_comp = a_Ez.nComp()-1;
  
  const int order = 1;
  const RealVect Xoffset = m_domainLeftEdge + 0.5*order*m_dx;
  const RealVect Xoffset_stag = Xoffset - 0.5*m_dx;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xp = this_particle.position();
     
     Ep = {0.0,0.0,0.0};
     Bp = {0.0,0.0,0.0};
 
     // get index for left most position that contains particle cloud  
     for (int dir=0; dir<SpaceDim; dir++) {
        index[dir] = floor((xp[dir] - Xoffset[dir])/m_dx[dir]);
        index_stag[dir] = floor((xp[dir] - Xoffset_stag[dir])/m_dx[dir]);
     }     
     
     for (int ll=0; ll<2; ll++) {
        ii[0] = index[0] + ll;
        ii_stag[0] = index_stag[0] + ll;

        delta[0] = (m_domainLeftEdge[0] + (ii[0] + 0.5)*m_dx[0] - xp[0])/m_dx[0]; 
        delta_stag[0] = (m_domainLeftEdge[0] + ii_stag[0]*m_dx[0] - xp[0])/m_dx[0];
        w0[0] = 1.0 - abs(delta[0]);
        w0_stag[0] = 1.0 - abs(delta_stag[0]);
           
        ii_Ex[0] = ii[0];
        ii_Ey[0] = ii_stag[0];
        ii_Ez[0] = ii_stag[0];

        w0_Ex[0] = w0[0];
        w0_Ey[0] = w0_stag[0];
        w0_Ez[0] = w0_stag[0];
           
        ii_Bx[0] = ii_stag[0];
        ii_By[0] = ii[0];
        ii_Bz[0] = ii[0];
        
        w0_Bx[0] = w0_stag[0];
        w0_By[0] = w0[0];
        w0_Bz[0] = w0[0];
#if CH_SPACEDIM>=2
        for (int mm=0; mm<2; mm++) {
           ii[1] = index[1] + mm;
           ii_stag[1] = index_stag[1] + mm;

           delta[1] = (m_domainLeftEdge[1] + (ii[1] + 0.5)*m_dx[1] - xp[1])/m_dx[1]; 
           delta_stag[1] = (m_domainLeftEdge[1] + ii_stag[1]*m_dx[1] - xp[1])/m_dx[1];
           w0[1] = 1.0 - abs(delta[1]);
           w0_stag[1] = 1.0 - abs(delta_stag[1]);

           ii_Ex[1] = ii_stag[1];
           ii_Ey[1] = ii[1];
           ii_Ez[1] = ii_stag[1];

           w0_Ex[1] = w0_stag[1];
           w0_Ey[1] = w0[1];
           w0_Ez[1] = w0_stag[1];
           
           ii_Bx[1] = ii[1];
           ii_By[1] = ii_stag[1];
           ii_Bz[1] = ii[1];
        
           w0_Bx[1] = w0[1];
           w0_By[1] = w0_stag[1];
           w0_Bz[1] = w0[1];
#endif
#if CH_SPACEDIM==3
           for (int nn=0; nn<2; nn++) {
              ii[2] = index[2] + nn;
              ii_stag[2] = index_stag[2] + nn;

              delta[2] = (m_domainLeftEdge[2] + (ii[2] + 0.5)*m_dx[2] - xp[2])/m_dx[2]; 
              delta_stag[2] = (m_domainLeftEdge[2] + ii_stag[2]*m_dx[2] - xp[2])/m_dx[2];
              w0[2] = 1.0 - abs(delta[2]);
              w0_stag[2] = 1.0 - abs(delta_stag[2]);

              ii_Ex[2] = ii_stag[2];
              ii_Ey[2] = ii_stag[2];
              ii_Ez[2] = ii[2];

              w0_Ex[2] = w0_stag[2];
              w0_Ey[2] = w0_stag[2];
              w0_Ez[2] = w0[2];
              
              ii_Bx[2] = ii[2];
              ii_By[2] = ii[2];
              ii_Bz[2] = ii_stag[2];
        
              w0_Bx[2] = w0[2];
              w0_By[2] = w0[2];
              w0_Bz[2] = w0_stag[2];
#endif

              Ep[0] = Ep[0] + w0_Ex.product()*a_Ex(ii_Ex,0);
              Ep[1] = Ep[1] + w0_Ey.product()*a_Ey(ii_Ey,0);
              Ep[2] = Ep[2] + w0_Ez.product()*a_Ez(ii_Ez,z_comp);

              Bp[0] = Bp[0] + w0_Bx.product()*a_Bx(ii_Bx,0);
              Bp[1] = Bp[1] + w0_By.product()*a_By(ii_By,0);
              Bp[2] = Bp[2] + w0_Bz.product()*a_Bz(ii_Bz,z_comp);

#if CH_SPACEDIM==3
           }    
#endif
#if CH_SPACEDIM>=2
        }
#endif
     }

     this_particle.setElectricField(Ep);
     this_particle.setMagneticField(Bp);
  }

}

template <class P>
void MeshInterp::interpolateBfieldsToPart_CIC( List<P>&     a_particleList,
                                         const FArrayBox&   a_Bx,
                                         const FArrayBox&   a_By,
                                         const FArrayBox&   a_Bz ) const
{
  CH_TIME("MeshInterp::interpolateBfieldsToPart_CIC()");
 
  std::array<Real,3> Bp;
  IntVect index, index_stag, ii, ii_stag; 
  IntVect ii_Bx, ii_By, ii_Bz; 
  RealVect w0_Bx, w0_By, w0_Bz; 
  RealVect delta, delta_stag, w0, w0_stag;

  const int z_comp = a_Bz.nComp()-1;
  
  const int order = 1;
  const RealVect Xoffset = m_domainLeftEdge + 0.5*order*m_dx;
  const RealVect Xoffset_stag = Xoffset - 0.5*m_dx;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xp = this_particle.position();
     
     Bp = {0.0,0.0,0.0};
 
     // get index for left most position that contains particle cloud  
     for (int dir=0; dir<SpaceDim; dir++) {
        index[dir] = floor((xp[dir] - Xoffset[dir])/m_dx[dir]);
        index_stag[dir] = floor((xp[dir] - Xoffset_stag[dir])/m_dx[dir]);
     }     
     
     for (int ll=0; ll<2; ll++) {
        ii[0] = index[0] + ll;
        ii_stag[0] = index_stag[0] + ll;

        delta[0] = (m_domainLeftEdge[0] + (ii[0] + 0.5)*m_dx[0] - xp[0])/m_dx[0]; 
        delta_stag[0] = (m_domainLeftEdge[0] + ii_stag[0]*m_dx[0] - xp[0])/m_dx[0];
        w0[0] = 1.0 - abs(delta[0]);
        w0_stag[0] = 1.0 - abs(delta_stag[0]);
           
        ii_Bx[0] = ii_stag[0];
        ii_By[0] = ii[0];
        ii_Bz[0] = ii[0];
        
        w0_Bx[0] = w0_stag[0];
        w0_By[0] = w0[0];
        w0_Bz[0] = w0[0];
#if CH_SPACEDIM>=2
        for (int mm=0; mm<2; mm++) {
           ii[1] = index[1] + mm;
           ii_stag[1] = index_stag[1] + mm;

           delta[1] = (m_domainLeftEdge[1] + (ii[1] + 0.5)*m_dx[1] - xp[1])/m_dx[1]; 
           delta_stag[1] = (m_domainLeftEdge[1] + ii_stag[1]*m_dx[1] - xp[1])/m_dx[1];
           w0[1] = 1.0 - abs(delta[1]);
           w0_stag[1] = 1.0 - abs(delta_stag[1]);

           ii_Bx[1] = ii[1];
           ii_By[1] = ii_stag[1];
           ii_Bz[1] = ii[1];
        
           w0_Bx[1] = w0[1];
           w0_By[1] = w0_stag[1];
           w0_Bz[1] = w0[1];
#endif
#if CH_SPACEDIM==3
           for (int nn=0; nn<2; nn++) {
              ii[2] = index[2] + nn;
              ii_stag[2] = index_stag[2] + nn;

              delta[2] = (m_domainLeftEdge[2] + (ii[2] + 0.5)*m_dx[2] - xp[2])/m_dx[2]; 
              delta_stag[2] = (m_domainLeftEdge[2] + ii_stag[2]*m_dx[2] - xp[2])/m_dx[2];
              w0[2] = 1.0 - abs(delta[2]);
              w0_stag[2] = 1.0 - abs(delta_stag[2]);

              ii_Bx[2] = ii[2];
              ii_By[2] = ii[2];
              ii_Bz[2] = ii_stag[2];
        
              w0_Bx[2] = w0[2];
              w0_By[2] = w0[2];
              w0_Bz[2] = w0_stag[2];
#endif

              Bp[0] = Bp[0] + w0_Bx.product()*a_Bx(ii_Bx,0);
              Bp[1] = Bp[1] + w0_By.product()*a_By(ii_By,0);
              Bp[2] = Bp[2] + w0_Bz.product()*a_Bz(ii_Bz,z_comp);

#if CH_SPACEDIM==3
           }    
#endif
#if CH_SPACEDIM>=2
        }
#endif
     }

     this_particle.setMagneticField(Bp);
  }

}

template <class P>
void MeshInterp::interpolateEfieldDirToPart_CIC( List<P>&     a_particleList,
                                           const FArrayBox&   a_E,
                                           const int          a_dir,
                                           const int          a_comp ) const
{
  CH_TIME("MeshInterp::interpolateEfieldsDirToPart_CIC()");
 
  IntVect index, ii; 
  RealVect delta, w0;
  Real E0;
   
  const IntVect stag = a_E.box().type();

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xp = this_particle.position();
     
     E0 = 0.0;
 
     // get index for left most position that contains particle cloud  
     for (int dir=0; dir<SpaceDim; dir++) {
        index[dir] = floor((xp[dir] - m_domainLeftEdge[dir] - 0.5*m_dx[dir]*(1-stag[dir]))/m_dx[dir]);
     }     
     
     for (int ll=0; ll<2; ll++) {
        ii[0] = index[0] + ll;
        delta[0] = (m_domainLeftEdge[0] + (ii[0] + 0.5*(1-stag[0]))*m_dx[0] - xp[0])/m_dx[0];
        w0[0] = 1.0 - abs(delta[0]);
#if CH_SPACEDIM>=2
        for (int mm=0; mm<2; mm++) {
           ii[1] = index[1] + mm;
           delta[1] = (m_domainLeftEdge[1] + (ii[1] + 0.5*(1-stag[1]))*m_dx[1] - xp[1])/m_dx[1]; 
           w0[1] = 1.0 - abs(delta[1]);
#endif
#if CH_SPACEDIM==3
           for (int nn=0; nn<2; nn++) {
              ii[2] = index[2] + nn;
              delta[2] = (m_domainLeftEdge[2] + (ii[2] + 0.5*(1-stag[2]))*m_dx[2] - xp[2])/m_dx[2]; 
              w0[2] = 1.0 - abs(delta[2]);
#endif

              E0 += w0.product()*a_E(ii,a_comp);

#if CH_SPACEDIM==3
           }    
#endif
#if CH_SPACEDIM>=2
        }
#endif
     }

     this_particle.setElectricField(E0,a_dir);
  }

}

template <class P>
void MeshInterp::interpolateEMfieldsToPart_TSC( List<P>&     a_particleList,
                                          const FArrayBox&   a_Ex,
                                          const FArrayBox&   a_Ey,
                                          const FArrayBox&   a_Ez,
                                          const FArrayBox&   a_Bx,
                                          const FArrayBox&   a_By,
                                          const FArrayBox&   a_Bz ) const
{
  CH_TIME("MeshInterp::interpolateEMfieldsToPart_TSC()");
 
  std::array<Real,3> Ep, Bp;
  IntVect index, index_stag, ii, ii_stag; 
  IntVect ii_Ex, ii_Ey, ii_Ez, ii_Bx, ii_By, ii_Bz; 
  RealVect w0_Ex, w0_Ey, w0_Ez, w0_Bx, w0_By, w0_Bz; 
  RealVect delta, delta_stag, w0, w0_stag;

  const int z_comp = a_Ez.nComp()-1;

  const int order = 2;
  const RealVect Xoffset = m_domainLeftEdge + 0.5*order*m_dx;
  const RealVect Xoffset_stag = Xoffset - 0.5*m_dx;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xp = this_particle.position();
     
     Ep = {0.0,0.0,0.0};
     Bp = {0.0,0.0,0.0};
 
     // get index for left most position that contains particle cloud  
     for (int dir=0; dir<SpaceDim; dir++) {
        index[dir] = floor((xp[dir] - Xoffset[dir])/m_dx[dir]);
        index_stag[dir] = floor((xp[dir] - Xoffset_stag[dir])/m_dx[dir]);
     }     
     
     for (int l=0; l<order+1; l++) {
        ii[0] = index[0] + l;
        ii_stag[0] = index_stag[0] + l;

        delta[0] = (m_domainLeftEdge[0] + (ii[0] + 0.5)*m_dx[0] - xp[0])/m_dx[0]; 
        delta_stag[0] = (m_domainLeftEdge[0] + ii_stag[0]*m_dx[0] - xp[0])/m_dx[0];

        if(abs(delta[0])<0.5) {
           w0[0] = 0.75 - pow(delta[0],2);
           w0_stag[0] = 0.75 - pow(delta_stag[0],2);
        }
        else {
           w0[0] = 0.5*pow(1.5 - abs(delta[0]),2);
           w0_stag[0] = 0.5*pow(1.5 - abs(delta_stag[0]),2);
        }
           
        ii_Ex[0] = ii[0];
        ii_Ey[0] = ii_stag[0];
        ii_Ez[0] = ii_stag[0];

        w0_Ex[0] = w0[0];
        w0_Ey[0] = w0_stag[0];
        w0_Ez[0] = w0_stag[0];
           
        ii_Bx[0] = ii_stag[0];
        ii_By[0] = ii[0];
        ii_Bz[0] = ii[0];
        
        w0_Bx[0] = w0_stag[0];
        w0_By[0] = w0[0];
        w0_Bz[0] = w0[0];
#if CH_SPACEDIM>=2
        for (int m=0; m<order+1; m++) {
           ii[1] = index[1] + m;
           ii_stag[1] = index_stag[1] + m;

           delta[1] = (m_domainLeftEdge[1] + (ii[1] + 0.5)*m_dx[1] - xp[1])/m_dx[1]; 
           delta_stag[1] = (m_domainLeftEdge[1] + ii_stag[1]*m_dx[1] - xp[1])/m_dx[1];

           if(abs(delta[1])<0.5) {
              w0[1] = 0.75 - pow(delta[1],2);
              w0_stag[1] = 0.75 - pow(delta_stag[1],2);
           }
           else {
              w0[1] = 0.5*pow(1.5 - abs(delta[1]),2);
              w0_stag[1] = 0.5*pow(1.5 - abs(delta_stag[1]),2);
           }

           ii_Ex[1] = ii_stag[1];
           ii_Ey[1] = ii[1];
           ii_Ez[1] = ii_stag[1];

           w0_Ex[1] = w0_stag[1];
           w0_Ey[1] = w0[1];
           w0_Ez[1] = w0_stag[1];
           
           ii_Bx[1] = ii[1];
           ii_By[1] = ii_stag[1];
           ii_Bz[1] = ii[1];
        
           w0_Bx[1] = w0[1];
           w0_By[1] = w0_stag[1];
           w0_Bz[1] = w0[1];
#endif
#if CH_SPACEDIM==3
           for (int n=0; n<order+1; n++) {
              ii[2] = index[2] + n;
              ii_stag[2] = index_stag[2] + n;

              delta[2] = (m_domainLeftEdge[2] + (ii[2] + 0.5)*m_dx[2] - xp[2])/m_dx[2]; 
              delta_stag[2] = (m_domainLeftEdge[2] + ii_stag[2]*m_dx[2] - xp[2])/m_dx[2];
           
              if(abs(delta[2])<0.5) {
                 w0[2] = 0.75 - pow(delta[2],2);
                 w0_stag[2] = 0.75 - pow(delta_stag[2],2);
              }
              else {
                 w0[2] = 0.5*pow(1.5 - abs(delta[2]),2);
                 w0_stag[2] = 0.5*pow(1.5 - abs(delta_stag[2]),2);
              }

              ii_Ex[2] = ii_stag[2];
              ii_Ey[2] = ii_stag[2];
              ii_Ez[2] = ii[2];

              w0_Ex[2] = w0_stag[2];
              w0_Ey[2] = w0_stag[2];
              w0_Ez2[] = w0[2];
              
              ii_Bx[2] = ii[2];
              ii_By[2] = ii[2];
              ii_Bz[2] = ii_stag[2];
        
              w0_Bx[2] = w0[2];
              w0_By[2] = w0[2];
              w0_Bz[2] = w0_stag[2];
#endif

              Ep[0] += w0_Ex.product()*a_Ex(ii_Ex,0);
              Ep[1] += w0_Ey.product()*a_Ey(ii_Ey,0);
              Ep[2] += w0_Ez.product()*a_Ez(ii_Ez,z_comp);

              Bp[0] += w0_Bx.product()*a_Bx(ii_Bx,0);
              Bp[1] += w0_By.product()*a_By(ii_By,0);
              Bp[2] += w0_Bz.product()*a_Bz(ii_Bz,z_comp);

#if CH_SPACEDIM==3
           }    
#endif
#if CH_SPACEDIM>=2
        }
#endif
     }

     this_particle.setElectricField(Ep);
     this_particle.setMagneticField(Bp);
  }

}

template <class P>
void MeshInterp::interpolateEToPart_CC0( List<P>&     a_particleList,
                                          const FArrayBox&   a_Ex ) const
{
  CH_TIME("MeshInterp::interpolateEToPart_CC0()");
 
  IntVect index_stag_old, index_stag_new; 
  IntVect ii, ii_next, sign, cell_crossings;
  RealVect xpnew, dXp, dXp_sub, xpold0, xpnew0, Xcell;
  Real w0_Ex, Epx; 
  int num_segments;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xpold = this_particle.position_old();
     const RealVect& xpbar = this_particle.position();

     // compute future particle position     
     xpnew = 2.0*xpbar - xpold;
     dXp = xpnew - xpold;
     sign = IntVect::Unit;

     Epx = 0.0;
 
     // get index for left most position for old and new 
     num_segments = 1;
     for (int dir=0; dir<SpaceDim; dir++) {
        index_stag_old[dir] = floor((xpold[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
        index_stag_new[dir] = floor((xpnew[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
        
        if(index_stag_new[dir]<index_stag_old[dir]) sign[dir] = -1;
        cell_crossings[dir] = abs(index_stag_new[dir]-index_stag_old[dir]);
        num_segments = num_segments + cell_crossings[dir];
        Xcell[dir] = m_domainLeftEdge[dir] + (index_stag_old[dir] + 0.5*(1.0-sign[dir]))*m_dx[dir]; 
     }

     xpold0 = xpold;
     ii_next = index_stag_old;

     for (int nn=0; nn<num_segments; nn++) {

        ii = ii_next;
 
        if(nn==num_segments-1) {
           xpnew0 = xpnew;
           dXp_sub = xpnew0 - xpold0;
        } 
        else {
           ii_next[0] = ii[0] + sign[0];
           Xcell[0] = Xcell[0] + sign[0]*m_dx[0];
           xpnew0[0] = Xcell[0];
           dXp_sub[0] = xpnew0[0] - xpold0[0];
        }

        if(num_segments>1) w0_Ex = dXp_sub[0]/dXp[0];
        else w0_Ex = 1.0;

        Epx += w0_Ex*a_Ex(ii,0);
       
        xpold0 = xpnew0;

     }

     this_particle.setElectricField(Epx,0);
  }

}

template <class P>
void MeshInterp::interpolateEToPart_CC0( List<P>&    a_particleList,
                                   const FArrayBox&  a_Ex,
                                   const FArrayBox&  a_Ey,
                                   const FArrayBox&  a_Ez ) const
{
  CH_TIME("MeshInterp::interpolateEToPart_CC0()");
 
  IntVect index, index_stag_old, index_stag_new; 
  IntVect ii, ii_next, sign, cell_crossings;
  IntVect ii_up0, ii_up1;
  RealVect xpnew, dXp, dXp_sub, xpold0, xpnew0, Xcell, delta, w0;
  Real weight, Epx, Epy, dXp_sub02, slope, slope_inv; 
  int num_segments, i0, i1;

  const IntVect stag = a_Ez.box().type();
  const int z_comp = a_Ez.nComp()-1;
  Real Epz;

  ii_up0 = IntVect::Zero;
  ii_up1 = IntVect::Zero;
  ii_up0[0] = 1;
  ii_up1[1] = 1;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xpold = this_particle.position_old();
     const RealVect& xpbar = this_particle.position();

     // compute future particle position     
     xpnew = 2.0*xpbar - xpold;
     dXp = xpnew - xpold;
     sign = IntVect::Unit;

     Epx = 0.0;
     Epy = 0.0;

     // compute slopes assuming 2D
     i0 = 0;
     i1 = SpaceDim - 1;
     slope = dXp[i1]/dXp[i0];
     slope_inv = 1.0/slope;
 
     // get index for left most position for old and new 
     num_segments = 1;
     for (int dir=0; dir<SpaceDim; dir++) {
        index_stag_old[dir] = floor((xpold[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
        index_stag_new[dir] = floor((xpnew[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
        
        if(index_stag_new[dir]<index_stag_old[dir]) sign[dir] = -1;
        cell_crossings[dir] = abs(index_stag_new[dir]-index_stag_old[dir]);
        num_segments = num_segments + cell_crossings[dir];
        Xcell[dir] = m_domainLeftEdge[dir] + (index_stag_old[dir] + 0.5*(1.0-sign[dir]))*m_dx[dir]; 
     }

     xpold0 = xpold;
     ii_next = index_stag_old;

     for (int nn=0; nn<num_segments; nn++) {

        ii = ii_next;
 
        if(nn==num_segments-1) {
           xpnew0 = xpnew;
           dXp_sub = xpnew0 - xpold0;
        } 
        else if(cell_crossings[i0]==0) {
           ii_next[i1] = ii[i1] + sign[i1];
           Xcell[i1] = Xcell[i1] + sign[i1]*m_dx[i1];
           xpnew0[i1] = Xcell[i1];
           dXp_sub[i1] = xpnew0[i1] - xpold0[i1];
           dXp_sub[i0] = slope_inv*dXp_sub[i1];
           xpnew0[i0] = xpold0[i0] + dXp_sub[i0];
        }
        else if(cell_crossings[i1]==0) {
           ii_next[i0] = ii[i0] + sign[i0];
           Xcell[i0] = Xcell[i0] + sign[i0]*m_dx[i0];
           xpnew0[i0] = Xcell[i0];
           dXp_sub[i0] = xpnew0[i0] - xpold0[i0];
           dXp_sub[i1] = slope*dXp_sub[i0];
           xpnew0[i1] = xpold0[i1] + dXp_sub[i1];
        }
        else {
           xpnew0[i0] = Xcell[i0] + sign[i0]*m_dx[i0];
           xpnew0[i1] = Xcell[i1] + sign[i1]*m_dx[i1];
           dXp_sub[i0] = xpnew0[i0] - xpold0[i0];
           dXp_sub[i1] = xpnew0[i1] - xpold0[i1];
           dXp_sub02 = slope_inv*dXp_sub[i1];
          
           if(abs(dXp_sub[i0])<abs(dXp_sub02)) {
              dXp_sub[i1] = slope*dXp_sub[i0];
              xpnew0[i1] = xpold0[i1] + dXp_sub[i1];
              Xcell[i0] = xpnew0[i0];
              ii_next[i0] = ii[i0] + sign[i0];
              cell_crossings[i0] = cell_crossings[i0] - 1;
           }
           else {
              dXp_sub[i0] = slope_inv*dXp_sub[i1];
              xpnew0[i0] = xpold0[i0] + dXp_sub[i0];
              Xcell[i1] = xpnew0[i1];
              ii_next[i1] = ii[i1] + sign[i1];
              cell_crossings[i1] = cell_crossings[i1] - 1;
           }
        }

        for (int dir=0; dir<SpaceDim; dir++) {
           delta[dir] = ( (xpold0[dir] + xpnew0[dir])/2.0 
                      -   (m_domainLeftEdge[dir] + ii[dir]*m_dx[dir]) )/m_dx[dir];
        }

        // interpolate Ex to particle
        if(num_segments>1) w0[i0] = dXp_sub[i0]/dXp[i0];
        else w0[i0] = 1.0;
        w0[i1] = 1.0 - delta[i1];

        weight = w0.product();
        Epx += weight*a_Ex(ii,0);

        w0[i1] = 1.0 - w0[i1];
        weight = w0.product();
        Epx += weight*a_Ex(ii+ii_up1,0);
       
        // interpolate Ey to particle
        w0[i0] = 1.0 - delta[i0];
        if(num_segments>1) w0[i1] = dXp_sub[i1]/dXp[i1];
        else w0[i1] = 1.0;

        weight = w0.product();
        Epy += weight*a_Ey(ii,0);

        w0[i0] = 1.0 - w0[i0];
        weight = w0.product();
        Epy += weight*a_Ey(ii+ii_up0,0);

        xpold0 = xpnew0;
     
        //
        //
        // use CIC for other components
        //
        //

        for (int dir=0; dir<SpaceDim; dir++) {
           index[dir] = floor( (xpbar[dir] - m_domainLeftEdge[dir] 
                                        - 0.5*m_dx[dir]*(1-stag[dir]))/m_dx[dir] );
        }     

        Epz = 0.0;
        for (int ll=0; ll<2; ll++) {
           ii[0] = index[0] + ll;
           delta[0] = (m_domainLeftEdge[0] + (ii[0] + 0.5*(1-stag[0]))*m_dx[0] - xpbar[0])/m_dx[0];
           w0[0] = 1.0 - abs(delta[0]);
#if CH_SPACEDIM>=2
           for (int mm=0; mm<2; mm++) {
              ii[1] = index[1] + mm;
              delta[1] = (m_domainLeftEdge[1] + (ii[1] + 0.5*(1-stag[1]))*m_dx[1] - xpbar[1])/m_dx[1]; 
              w0[1] = 1.0 - abs(delta[1]);
#endif
#if CH_SPACEDIM==3
              for (int nn=0; nn<2; nn++) {
                 ii[2] = index[2] + nn;
                 delta[2] = (m_domainLeftEdge[2] + (ii[2] + 0.5*(1-stag[2]))*m_dx[2] - xpbar[2])/m_dx[2]; 
                 w0[2] = 1.0 - abs(delta[2]);
#endif
                    Epz += w0.product()*a_Ez(ii,z_comp);
#if CH_SPACEDIM==3
              }    
#endif
#if CH_SPACEDIM>=2
           }
#endif
        }

     }

     this_particle.setElectricField(Epx,0);
     this_particle.setElectricField(Epy,1);
     this_particle.setElectricField(Epz,2);

  }

}

template <class P>
void MeshInterp::interpolateEToPart_CC1( List<P>&     a_particleList,
                                          const FArrayBox&   a_Ex ) const
{
  CH_TIME("MeshInterp::interpolateEToPart_CC1()");
 
  IntVect index_old, index_new; 
  IntVect ii, ii_next, sign, cell_crossings;
  RealVect xpnew, dXp, dXp_sub, xpold0, xpnew0, Xcell, xpbar0;
  Real w0_Ex, Epx, l0, w0, seg_factor; 
  int num_segments;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xpold = this_particle.position_old();
     const RealVect& xpbar = this_particle.position();

     // compute future particle position     
     xpnew = 2.0*xpbar - xpold;
     dXp = xpnew - xpold;
     sign = IntVect::Unit;

     Epx = 0.0;
 
     // get index for left most position for old and new 
     num_segments = 1;
     for (int dir=0; dir<SpaceDim; dir++) {
        index_old[dir] = floor((xpold[dir] - m_domainLeftEdge[dir] - 0.5*m_dx[dir])/m_dx[dir]);
        index_new[dir] = floor((xpnew[dir] - m_domainLeftEdge[dir] - 0.5*m_dx[dir])/m_dx[dir]);
        
        if(index_new[dir]<index_old[dir]) sign[dir] = -1;
        cell_crossings[dir] = abs(index_new[dir]-index_old[dir]);
        num_segments = num_segments + cell_crossings[dir];
        Xcell[dir] = m_domainLeftEdge[dir] + (index_old[dir] + 0.5*(1.0-sign[dir]) + 0.5)*m_dx[dir]; 
     }

     xpold0 = xpold;
     ii_next = index_old;

     for (int nn=0; nn<num_segments; nn++) {

        ii = ii_next;
 
        if(nn==num_segments-1) {
           xpnew0 = xpnew;
           dXp_sub = xpnew0 - xpold0;
        } 
        else {
           ii_next[0] = ii[0] + sign[0];
           Xcell[0] = Xcell[0] + sign[0]*m_dx[0];
           xpnew0[0] = Xcell[0];
           dXp_sub[0] = xpnew0[0] - xpold0[0];
        }

        xpbar0[0] = 0.5*(xpnew0[0] + xpold0[0]);
        l0 = m_domainLeftEdge[0] + ii[0]*m_dx[0] + 0.5*m_dx[0] - xpbar0[0];
        w0 = 1.0 - abs(l0/m_dx[0]);

        if(num_segments>1) seg_factor = dXp_sub[0]/dXp[0];
        else seg_factor = 1.0;

        w0_Ex = w0*seg_factor;
        Epx += w0_Ex*a_Ex(ii,0);
        
        w0_Ex = (1.0-w0)*seg_factor;
        ii[0] = ii[0] + 1;
        Epx += w0_Ex*a_Ex(ii,0);
       
        xpold0 = xpnew0;

     }

     this_particle.setElectricField(Epx,0);
  }

}

template <class P>
void MeshInterp::interpolateEToPart_CC1( List<P>&    a_particleList,
                                   const FArrayBox&  a_Ex,
                                   const FArrayBox&  a_Ey ) const
{
  CH_TIME("MeshInterp::interpolateEToPart_CC1()");

  //
  //  JUST A COPY OF _CC0 in 2D. 
  //
 
  IntVect index_stag_old, index_stag_new; 
  IntVect ii, ii_next, sign, cell_crossings;
  IntVect ii_up0, ii_up1;
  RealVect xpnew, dXp, dXp_sub, xpold0, xpnew0, Xcell, delta;
  Real weight, w0, w1, Epx, Epy, dXp_sub02, slope, slope_inv; 
  int num_segments, i0, i1;

  ii_up0 = IntVect::Zero;
  ii_up1 = IntVect::Zero;
  ii_up0[0] = 1;
  ii_up1[1] = 1;

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {

     P& this_particle = lit();
     const RealVect& xpold = this_particle.position_old();
     const RealVect& xpbar = this_particle.position();

     // compute future particle position     
     xpnew = 2.0*xpbar - xpold;
     dXp = xpnew - xpold;
     sign = IntVect::Unit;

     Epx = 0.0;
     Epy = 0.0;

     // compute slopes assuming 2D
     i0 = 0;
     i1 = SpaceDim - 1;
     slope = dXp[i1]/dXp[i0];
     slope_inv = 1.0/slope;
 
     // get index for left most position for old and new 
     num_segments = 1;
     for (int dir=0; dir<SpaceDim; dir++) {
        index_stag_old[dir] = floor((xpold[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
        index_stag_new[dir] = floor((xpnew[dir] - m_domainLeftEdge[dir])/m_dx[dir]);
        
        if(index_stag_new[dir]<index_stag_old[dir]) sign[dir] = -1;
        cell_crossings[dir] = abs(index_stag_new[dir]-index_stag_old[dir]);
        num_segments = num_segments + cell_crossings[dir];
        Xcell[dir] = m_domainLeftEdge[dir] + (index_stag_old[dir] + 0.5*(1.0-sign[dir]))*m_dx[dir]; 
     }

     xpold0 = xpold;
     ii_next = index_stag_old;

     for (int nn=0; nn<num_segments; nn++) {

        ii = ii_next;
 
        if(nn==num_segments-1) {
           xpnew0 = xpnew;
           dXp_sub = xpnew0 - xpold0;
        } 
        else if(cell_crossings[i0]==0) {
           ii_next[i1] = ii[i1] + sign[i1];
           Xcell[i1] = Xcell[i1] + sign[i1]*m_dx[i1];
           xpnew0[i1] = Xcell[i1];
           dXp_sub[i1] = xpnew0[i1] - xpold0[i1];
           dXp_sub[i0] = slope_inv*dXp_sub[i1];
           xpnew0[i0] = xpold0[i0] + dXp_sub[i0];
        }
        else if(cell_crossings[i1]==0) {
           ii_next[i0] = ii[i0] + sign[i0];
           Xcell[i0] = Xcell[i0] + sign[i0]*m_dx[i0];
           xpnew0[i0] = Xcell[i0];
           dXp_sub[i0] = xpnew0[i0] - xpold0[i0];
           dXp_sub[i1] = slope*dXp_sub[i0];
           xpnew0[i1] = xpold0[i1] + dXp_sub[i1];
        }
        else {
           xpnew0[i0] = Xcell[i0] + sign[i0]*m_dx[i0];
           xpnew0[i1] = Xcell[i1] + sign[i1]*m_dx[i1];
           dXp_sub[i0] = xpnew0[i0] - xpold0[i0];
           dXp_sub[i1] = xpnew0[i1] - xpold0[i1];
           dXp_sub02 = slope_inv*dXp_sub[i1];
          
           if(abs(dXp_sub[i0])<abs(dXp_sub02)) {
              dXp_sub[i1] = slope*dXp_sub[i0];
              xpnew0[i1] = xpold0[i1] + dXp_sub[i1];
              Xcell[i0] = xpnew0[i0];
              ii_next[i0] = ii[i0] + sign[i0];
              cell_crossings[i0] = cell_crossings[i0] - 1;
           }
           else {
              dXp_sub[i0] = slope_inv*dXp_sub[i1];
              xpnew0[i0] = xpold0[i0] + dXp_sub[i0];
              Xcell[i1] = xpnew0[i1];
              ii_next[i1] = ii[i1] + sign[i1];
              cell_crossings[i1] = cell_crossings[i1] - 1;
           }
        }

        for (int dir=0; dir<SpaceDim; dir++) {
           delta[dir] = 0.5*(xpold0[dir] + xpnew0[dir] - (m_domainLeftEdge[dir] + ii[dir]*m_dx[dir]))/m_dx[dir];
        }

        if(num_segments>1) w0 = dXp_sub[0]/dXp[0];
        else w0 = 1.0;
        w1 = 1.0 - delta[1];

        weight = w0*w1;
        Epx += weight*a_Ex(ii,0);

        weight = w0*(1.0-w1);
        Epx += weight*a_Ex(ii+ii_up0,0);
       
        w0 = 1.0 - delta[0];
        if(num_segments>1) w1 = dXp_sub[1]/dXp[1];
        else w1 = 1.0;

        weight = w0*w1;
        Epy += weight*a_Ey(ii,0);

        weight = (1.0-w0)*w1;
        Epy += weight*a_Ey(ii+ii_up1,0);

        xpold0 = xpnew0;

     }

     this_particle.setElectricField(Epx,i0);
     this_particle.setElectricField(Epy,i1);

  }

}


#include "NamespaceFooter.H"

#endif // include guard
