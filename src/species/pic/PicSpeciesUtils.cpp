
#include "PicSpeciesUtils.H"
#include "PicnicConstants.H"
#include <iostream>

#include "NamespaceHeader.H"
   
void 
PicSpeciesUtils::applyForces( List<JustinsParticle>&  a_pList,
                        const Real                    a_fnorm,
                        const Real                    a_cnormDt, 
                        const bool                    a_byHalfDt,
                        const bool                    a_anticyclic )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForces()");
   
   // compute upbar=(up^n + up^{n+1})/2 using Boris method
   // upbar = upn + dt/2*q/m*(Ep(xpbar) + upbar x Bp(xpbar)/gammap)
   // if byHalfDt is false, convert upbar to upnew: up -> 2*up - upold

   Real bp0, bp1, bp2, denom;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;
   
   Real alpha = a_fnorm*a_cnormDt/2.0;

   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   std::array<int,3> dirp = {0,1,2};
   if(a_anticyclic) {
      dirp[1] = 2;
      dirp[2] = 1;
   }

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      std::array<Real,3>& up = lit().velocity(); // gamma*beta
      const std::array<Real,3>& upold = lit().velocity_old();
      const std::array<Real,3>& Ep = lit().electric_field();
      const std::array<Real,3>& Bp = lit().magnetic_field();
   
      // add half acceleration to old velocity
      vm0 = upold[dirp[0]] + alpha*Ep[dirp[0]];
      vm1 = upold[dirp[1]] + alpha*Ep[dirp[1]];
      vm2 = upold[dirp[2]] + alpha*Ep[dirp[2]];
      
      // scale Bp by alpha
      bp0 = alpha*Bp[dirp[0]];
      bp1 = alpha*Bp[dirp[1]];
      bp2 = alpha*Bp[dirp[2]];
      
#ifdef RELATIVISTIC_PARTICLES
      Real gammap = sqrt(1.0 + vm0*vm0 + vm1*vm1 + vm2*vm2);
      bp0 /= gammap;
      bp1 /= gammap;
      bp2 /= gammap;
#endif
   
      denom = 1.0 + bp0*bp0 + bp1*bp1 + bp2*bp2;

      // define vpr = vm + vm x bp
      vpr0 = vm0 + vm1*bp2 - vm2*bp1;
      vpr1 = vm1 + vm2*bp0 - vm0*bp2;
      vpr2 = vm2 + vm0*bp1 - vm1*bp0;

      // compute upbar = vm + (vpr x bp)/denom (energy conserving to machine precision)
      up[dirp[0]] = vm0 + (vpr1*bp2 - vpr2*bp1)/denom;
      up[dirp[1]] = vm1 + (vpr2*bp0 - vpr0*bp2)/denom;
      up[dirp[2]] = vm2 + (vpr0*bp1 - vpr1*bp0)/denom;

      // compute upbar with 2nd formulation (kinda works)
      //up[dirp[0]] = (vm0*denom + vpr1*bp2 - vpr2*bp1)/denom;
      //up[dirp[1]] = (vm1*denom + vpr2*bp0 - vpr0*bp2)/denom;
      //up[dirp[2]] = (vm2*denom + vpr0*bp1 - vpr1*bp0)/denom;

      // compute upbar with 3rd formulation (not energy conserving to machine precision)
      //Real vmbp = vm0*bp0 + vm1*bp1 + vm2*bp2;
      //up[dirp[0]] = (vm0 + vm1*bp2 - vm2*bp1 + vmbp*bp0)/denom;
      //up[dirp[1]] = (vm1 + vm2*bp0 - vm0*bp2 + vmbp*bp1)/denom;
      //up[dirp[2]] = (vm2 + vm0*bp1 - vm1*bp0 + vmbp*bp2)/denom;
      
      // compute upbar with 4th formulation (not energy conserving to machine precision)
      //Real vmbp = vm0*bp0 + vm1*bp1 + vm2*bp2;
      //up[dirp[0]] = vm0 + (vm0*(1.0-denom) + vm1*bp2 - vm2*bp1 + vmbp*bp0)/denom;
      //up[dirp[1]] = vm1 + (vm1*(1.0-denom) + vm2*bp0 - vm0*bp2 + vmbp*bp1)/denom;
      //up[dirp[2]] = vm2 + (vm2*(1.0-denom) + vm0*bp1 - vm1*bp0 + vmbp*bp2)/denom;
      
      // compute upbar with 5th formulation (works)
      //Real vmbp = vm0*bp0 + vm1*bp1 + vm2*bp2;
      //up[dirp[0]] = vm0 + (-vm0*(bp0*bp0 + bp1*bp1 + bp2*bp2) + vm1*bp2 - vm2*bp1 + vmbp*bp0)/denom;
      //up[dirp[1]] = vm1 + (-vm1*(bp0*bp0 + bp1*bp1 + bp2*bp2) + vm2*bp0 - vm0*bp2 + vmbp*bp1)/denom;
      //up[dirp[2]] = vm2 + (-vm2*(bp0*bp0 + bp1*bp1 + bp2*bp2) + vm0*bp1 - vm1*bp0 + vmbp*bp2)/denom;
      
      // compute upbar with 6th formulation (works)
      //up[dirp[0]] = vm0 + (-vm0*(bp1*bp1 + bp2*bp2) + vm1*bp2 - vm2*bp1 + (vm1*bp1 + vm2*bp2)*bp0)/denom;
      //up[dirp[1]] = vm1 + (-vm1*(bp0*bp0 + bp2*bp2) + vm2*bp0 - vm0*bp2 + (vm0*bp0 + vm2*bp2)*bp1)/denom;
      //up[dirp[2]] = vm2 + (-vm2*(bp0*bp0 + bp1*bp1) + vm0*bp1 - vm1*bp0 + (vm0*bp0 + vm1*bp1)*bp2)/denom;
      
      if(!a_byHalfDt) {
         up[dirp[0]] = 2.0*up[dirp[0]] - upold[dirp[0]];
         up[dirp[1]] = 2.0*up[dirp[1]] - upold[dirp[1]];
         up[dirp[2]] = 2.0*up[dirp[2]] - upold[dirp[2]];
      }

   } // end loop over particle list

}

void 
PicSpeciesUtils::applyForces_CYL_CYL( List<JustinsParticle>&  a_pList,
                                const Real                    a_fnorm,
                                const Real                    a_cnormDt,
                                const bool                    a_byHalfDt,
                                const bool                    a_anticyclic )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForces_CYL_CYL()");
   
   // advance up from step n to step n+1 using modified Boris algorithm
   // that uses an iterative method to include the inertia term directly.
   // Only for cylindrical so far. In cylindrical, the equations of motion
   // can be expressed same as they are for Cartesian with the transformation
   // Bpz -> Bpz + vpthbar*mp/qp/rp. So one can still use the Boris algorithm,
   // but iterations are needed in order for the nonlinear term to converge.
   //
   // upbarr  = uprn  + dt/2*qp/mp*(Epr  + (upthbar*Bpz_mod - upzbar*Bpth)/gp)
   // upbarth = upthn + dt/2*qp/mp*(Epth + (upzbar*Bpr      - uprbar*Bpz_mod)/gp)
   // upbarz  = upzn  + dt/2*qp/mp*(Epz  + (uprbar*Bpth     - upthbar*Bpr)/gp)
   // with Bpz_mod = Bpz + upthbar*mp/qp/rp
   //
   // if byHalfDt = false, convert upbar to up: up -> (up + upold)/2
   //
   // Note that 1) convergence is not needed for exact energy conservervation
   //       and 2) 2nd order accurate after 2 iterations

   // advance velocities using Boris algorithm
   Real bp0, bp1, bp2, denom, gammap;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;

   Real alpha = a_fnorm*a_cnormDt/2.0;

   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   std::array<int,3> dirp = {0,1,2};
   if(a_anticyclic) {
      dirp[1] = 2;
      dirp[2] = 1;
   }

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
      const std::array<Real,3>& Ep = lit().electric_field();
      const std::array<Real,3>& Bp = lit().magnetic_field();
    
      const RealVect& xpold = lit().position_old();
  
      // add half acceleration to old velocity
      vm0 = upold[dirp[0]] + alpha*Ep[dirp[0]];
      vm1 = upold[dirp[1]] + alpha*Ep[dirp[1]];
      vm2 = upold[dirp[2]] + alpha*Ep[dirp[2]];
      
      // scale Bp by alpha
      bp0 = alpha*Bp[dirp[0]];
      bp1 = alpha*Bp[dirp[1]];

#ifdef RELATIVISTIC_PARTICLES
      gammap = sqrt(1.0 + vm0*vm0 + vm1*vm1 + vm2*vm2);
      bp0 /= gammap;
      bp1 /= gammap;
#else
      gammap = 1.0;
#endif
      
      Real dtheta = lit().position_virt(0);
      bool set_dtheta = false; 
      if(dtheta==0.0) { // first-order predictor for dtheta
         dtheta = a_cnormDt/2.0*upold[dirp[1]]/xpold[dirp[0]]/gammap;
         set_dtheta = true;
      }
   
      // add inertia term to bpz
      bp2 = alpha*Bp[dirp[2]]/gammap + sin(dtheta);
      denom = 1.0 + bp0*bp0 + bp1*bp1 + bp2*bp2;

      // define vpr = vm + vm x bp
      vpr0 = vm0 + vm1*bp2 - vm2*bp1;
      vpr1 = vm1 + vm2*bp0 - vm0*bp2;
      vpr2 = vm2 + vm0*bp1 - vm1*bp0;

      // add another half acceleration
      up[dirp[0]] = vm0 + (vpr1*bp2 - vpr2*bp1)/denom;
      up[dirp[1]] = vm1 + (vpr2*bp0 - vpr0*bp2)/denom;
      up[dirp[2]] = vm2 + (vpr0*bp1 - vpr1*bp0)/denom;
      
      if(set_dtheta) { // 2nd-order corrector for dtheta
         Real rpbar = xpold[dirp[0]] + a_cnormDt/2.0*up[dirp[0]];
         dtheta = a_cnormDt/2.0*up[dirp[1]]/rpbar/gammap;
	 std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
         position_virt[0] = dtheta;
      }

      if(!a_byHalfDt) {
         up[dirp[0]] = 2.0*up[dirp[0]] - upold[dirp[0]];
         up[dirp[1]] = 2.0*up[dirp[1]] - upold[dirp[1]];
         up[dirp[2]] = 2.0*up[dirp[2]] - upold[dirp[2]];
      }

   } // end loop over particle list

}

void 
PicSpeciesUtils::applyForces_SPH_SPH( List<JustinsParticle>&  a_pList,
                                const Real                    a_fnorm,
                                const Real                    a_cnormDt,
                                const bool                    a_byHalfDt )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForces_SPH_SPH()");
   
   // advance up from step n to step n+1 using modified Boris algorithm
   // for spherical geometry with predictor-corrector for dtheta and dphi
   //
   // upbarr  = uprn  + dt/2*qp/mp*(Epr  + (upthbar*Bpph_mod - upphbar*Bpth_mod)/gp)
   // upbarth = upthn + dt/2*qp/mp*(Epth + (upphbar*Bpr_mod  - uprbar*Bpph_mod)/gp)
   // upbarph = upphn + dt/2*qp/mp*(Epph + (uprbar*Bpth_mod - upthbar*Bpr_mod)/gp)
   // with Bpr_mod  = Bpr  + mp/qp*sin(phi)*dtheta/dt
   //      Bpth_mod = Bpth - mp/qp*dphi/dt
   //      Bpph_mod = Bpph + mp/qp*cos(phi)*dtheta/dt
   //
   // if byHalfDt = false, convert upbar to up: up -> (up + upold)/2

   // advance velocities using Boris algorithm
   Real bp0, bp1, bp2, denom, gammap;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;

   Real alpha = a_fnorm*a_cnormDt/2.0;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
      const std::array<Real,3>& Ep = lit().electric_field();
      const std::array<Real,3>& Bp = lit().magnetic_field();
    
      const RealVect& xpold = lit().position_old();
  
      // add half acceleration to old velocity
      vm0 = upold[0] + alpha*Ep[0];
      vm1 = upold[1] + alpha*Ep[1];
      vm2 = upold[2] + alpha*Ep[2];
      
#ifdef RELATIVISTIC_PARTICLES
      gammap = sqrt(1.0 + vm0*vm0 + vm1*vm1 + vm2*vm2);
#else
      gammap = 1.0;
#endif
      
      // scale Bp by alpha/gammap
      bp0 = alpha*Bp[0]/gammap;
      bp1 = alpha*Bp[1]/gammap;
      bp2 = alpha*Bp[2]/gammap;

      
      Real dtheta = lit().position_virt(0);
      Real dphi = lit().position_virt(1);
      bool set_dtheta = false; 
      if(dtheta==0.0) { // first-order predictor for dtheta and dphi
         dtheta = a_cnormDt/2.0*upold[1]/xpold[0]/gammap;
         dphi = a_cnormDt/2.0*upold[2]/xpold[0]/gammap;
	 dtheta = std::sin(dtheta); // avoid large angles
	 dphi = std::sin(dphi);     // avoid large angles
         set_dtheta = true;
      }
   
      // add inertia terms to Bp
      bp0 += std::sin(dphi)*dtheta;
      bp1 -= dphi;
      bp2 += std::cos(dphi)*dtheta;
      denom = 1.0 + bp0*bp0 + bp1*bp1 + bp2*bp2;

      // define vpr = vm + vm x bp
      vpr0 = vm0 + vm1*bp2 - vm2*bp1;
      vpr1 = vm1 + vm2*bp0 - vm0*bp2;
      vpr2 = vm2 + vm0*bp1 - vm1*bp0;

      // add another half acceleration
      up[0] = vm0 + (vpr1*bp2 - vpr2*bp1)/denom;
      up[1] = vm1 + (vpr2*bp0 - vpr0*bp2)/denom;
      up[2] = vm2 + (vpr0*bp1 - vpr1*bp0)/denom;
      
      if(set_dtheta) { // 2nd-order corrector for dtheta
         Real rpbar = xpold[0] + a_cnormDt/2.0*up[0];
         dphi = a_cnormDt/2.0*up[2]/rpbar/gammap;
	 dphi = std::sin(dphi);     // avoid large angles
         dtheta = a_cnormDt/2.0*up[1]/rpbar/gammap/std::cos(dphi);
	 dtheta = std::sin(dtheta); // avoid large angles
         std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
         position_virt[0] = dtheta;
         position_virt[1] = dphi;
      }

      if(!a_byHalfDt) {
         up[0] = 2.0*up[0] - upold[0];
         up[1] = 2.0*up[1] - upold[1];
         up[2] = 2.0*up[2] - upold[2];
      }

   } // end loop over particle list

}
      

void 
PicSpeciesUtils::applyForces_CYL_HYB( List<JustinsParticle>&  a_pList,
                               const Real                    a_fnorm,
                               const Real                    a_cnormDt,
                               const bool                    a_anticyclic )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForces_CYL_HYB()");
   
   // advance velocities using Boris algorithm
   Real bp0, bp1, bp2, denom, gammap;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;

   Real alpha = a_fnorm*a_cnormDt/2.0;

   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   std::array<int,3> dirp = {0,1,2};
   if(a_anticyclic) {
      dirp[1] = 2;
      dirp[2] = 1;
   }

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
      const std::array<Real,3>& Ep = lit().electric_field();
      const std::array<Real,3>& Bp = lit().magnetic_field();
    
      // add half acceleration to old velocity
      vm0 = upold[dirp[0]] + alpha*Ep[dirp[0]];
      vm1 = upold[dirp[1]] + alpha*Ep[dirp[1]];
      vm2 = upold[dirp[2]] + alpha*Ep[dirp[2]];
      
      // scale Bp by alpha
      bp0 = alpha*Bp[dirp[0]];
      bp1 = alpha*Bp[dirp[1]];
      bp2 = alpha*Bp[dirp[2]];

#ifdef RELATIVISTIC_PARTICLES
      gammap = sqrt(1.0 + vm0*vm0 + vm1*vm1 + vm2*vm2);
      bp0 /= gammap;
      bp1 /= gammap;
      bp2 /= gammap;
#else
      gammap = 1.0;
#endif
      denom = 1.0 + bp0*bp0 + bp1*bp1 + bp2*bp2;
      
      // redefine vm0 and vm1 using modified definitions for upoldr and upoldth
      const Real dthp = lit().position_virt(0);
      const Real costhp = std::cos(dthp);
      const Real sinthp = std::sin(dthp);
      
      const Real upoldr_2  = costhp*upold[dirp[0]] + sinthp*upold[dirp[1]]; 
      const Real upoldth_2 = -sinthp*upold[dirp[0]] + costhp*upold[dirp[1]]; 
      vm0 = upoldr_2 + alpha*Ep[dirp[0]];
      vm1 = upoldth_2 + alpha*Ep[dirp[1]];
 
      // define vpr = vm + vm x bp
      vpr0 = vm0 + vm1*bp2 - vm2*bp1;
      vpr1 = vm1 + vm2*bp0 - vm0*bp2;
      vpr2 = vm2 + vm0*bp1 - vm1*bp0;

      // add another half acceleration
      up[dirp[0]] = vm0 + (vpr1*bp2 - vpr2*bp1)/denom;
      up[dirp[1]] = vm1 + (vpr2*bp0 - vpr0*bp2)/denom;
      up[dirp[2]] = vm2 + (vpr0*bp1 - vpr1*bp0)/denom;
      
   } // end loop over particle list

}

void 
PicSpeciesUtils::applyForces_SPH_HYB( List<JustinsParticle>&  a_pList,
                                const Real                    a_fnorm,
                                const Real                    a_cnormDt )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForces_SPH_HYB()");
      
   // In the HYBRID method, the velocity advance is done in CAR (to avoid needing to resolve
   // inertial forces), but the position advance is done in SPH (for charge conservation).
   // The straightforward approach is to 1) convert the velocity and fields from SPH to CAR,
   // 2) update the velocity in CAR, 3) convert that velocity back to SPH to do the position
   // update. However, we do it in the following more compact way. We keep the velocity and fields
   // in SPH (no conversion to CAR), do the advance as if the velocity and fields were in CAR, 
   // but use a modified upold vector that is formed by 1) transforming upold from SPH to CAR
   // using old values for theta (Pi/2) and phi (0), and then transforming this CAR upold back to
   // SPH using the time-centered values for theta and phi. 
   
   // advance velocities using Boris algorithm
   Real bp0, bp1, bp2, denom, gammap;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;

   Real thp, php, costhp, sinthp, cosphp, sinphp;
   Real upxold, upyold, upzold;
   Real upoldr_2, upoldth_2, upoldph_2;

   Real alpha = a_fnorm*a_cnormDt/2.0;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
      const std::array<Real,3>& Ep = lit().electric_field();
      const std::array<Real,3>& Bp = lit().magnetic_field();
    
      // add half acceleration to old velocity
      vm0 = upold[0] + alpha*Ep[0];
      vm1 = upold[1] + alpha*Ep[1];
      vm2 = upold[2] + alpha*Ep[2];
      
      // scale Bp by alpha
      bp0 = alpha*Bp[0];
      bp1 = alpha*Bp[1];
      bp2 = alpha*Bp[2];

#ifdef RELATIVISTIC_PARTICLES
      gammap = sqrt(1.0 + vm0*vm0 + vm1*vm1 + vm2*vm2);
      bp0 /= gammap;
      bp1 /= gammap;
      bp2 /= gammap;
#else
      gammap = 1.0;
#endif
      denom = 1.0 + bp0*bp0 + bp1*bp1 + bp2*bp2;
      
      thp = lit().position_virt(0);
      php = lit().position_virt(1);
      costhp = std::cos(thp);
      sinthp = std::sin(thp);
      cosphp = std::cos(php);
      sinphp = std::sin(php);
      
      // compute CAR upold from SPH upold using theta = 0.0 and phi = 0.0
      upxold = upold[0];
      upyold = upold[1];
      upzold = upold[2];

      // transform CAR upold to SPH upold using time-centered values for angles
      upoldr_2  =  cosphp*(costhp*upxold + sinthp*upyold) + sinphp*upzold; 
      upoldth_2 = -sinthp*upxold + costhp*upyold;
      upoldph_2 = -sinphp*(costhp*upxold + sinthp*upyold) + cosphp*upzold; 
      
      // redefine vm using modified definitions for upold
      vm0 = upoldr_2  + alpha*Ep[0];
      vm1 = upoldth_2 + alpha*Ep[1];
      vm2 = upoldph_2 + alpha*Ep[2];
 
      // define vpr = vm + vm x bp
      vpr0 = vm0 + vm1*bp2 - vm2*bp1;
      vpr1 = vm1 + vm2*bp0 - vm0*bp2;
      vpr2 = vm2 + vm0*bp1 - vm1*bp0;

      // add another half acceleration
      up[0] = vm0 + (vpr1*bp2 - vpr2*bp1)/denom;
      up[1] = vm1 + (vpr2*bp0 - vpr0*bp2)/denom;
      up[2] = vm2 + (vpr0*bp1 - vpr1*bp0)/denom;
      
   } // end loop over particle list

}

void PicSpeciesUtils::stepNormTransfer_CYL_CAR( List<JustinsParticle>&  a_pList,
                                        List<JustinsParticle>&  a_temp_pList,
                                  const RealVect&               a_dX,
                                  const Real                    a_cnormHalfDt, 
                                  const Real                    a_rtol, 
                                  const bool                    a_reverse )
{
   CH_TIME("PicSpeciesUtils::stepNormTransfer_CYL_CAR()");

   // Check for convergence of xpbar-xpold and vpbar(xpbar)*cnormDt/2.
   // Note that vpbar has to be updated after xpbar before coming here.
   // a_reverse = false: if not converged, move particle from pList to temp_pList
   // a_reverse = true: if converged, move particle from temp_pList to pList

   Real dxp0, rel_diff_max;
   RealVect dxp;
   Real xpnew, ypnew, rpnew;
   
   std::array<int,3> dirp = {0,1,2};
#if CH_SPACEDIM==2
   dirp[1] = 2;
   dirp[2] = 1;
#endif

   if(a_reverse) { // if converged, temp_pList ==> pList, else update position 
      
      ListIterator<JustinsParticle> lit(a_temp_pList);
      for(lit.begin(); lit.ok();) {
         
         RealVect& xpbar = lit().position();
	 const RealVect& xpold = lit().position_old();
         const std::array<Real,3>& upbar = lit().velocity();
	 std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();

         // compute rpnew using updated value of CAR positions	 
         xpnew = xpold[0] + 2.0*a_cnormHalfDt*upbar[0];
         ypnew = 2.0*a_cnormHalfDt*upbar[dirp[1]];
         rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew);

         dxp0 = xpbar[0]-xpold[0];
         dxp[0] = (rpnew - xpold[0])/2.0;
         rel_diff_max = std::abs(dxp0-dxp[0])/a_dX[0];
#if CH_SPACEDIM==2
         dxp0 = xpbar[1]-xpold[1];
         dxp[1]  = upbar[dirp[2]]*a_cnormHalfDt;
         Real rel_diff_dir = std::abs(dxp0-dxp[1])/a_dX[1];
         rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
#endif             

         if (rel_diff_max < a_rtol) a_pList.transfer(lit);
         else {
            for(int dir=0; dir<SpaceDim; dir++) xpbar[dir] = xpold[dir] + dxp[dir];
	    position_virt[0] = (xpnew + xpold[0])/2.0;
	    position_virt[1] = ypnew/2.0;
            ++lit;
         }

      }

   }
   else { // if not converged, update position and send to temp_pList
 
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok();) {
         
         RealVect& xpbar = lit().position();
         const RealVect& xpold = lit().position_old();
         const std::array<Real,3>& upbar = lit().velocity();
	 std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();

         // compute rpnew using updated value of CAR positions	 
         xpnew = xpold[0] + 2.0*a_cnormHalfDt*upbar[0];
         ypnew = 2.0*a_cnormHalfDt*upbar[dirp[1]];
         rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew);
 
         dxp0 = xpbar[0]-xpold[0];
         dxp[0] = (rpnew - xpold[0])/2.0;
         rel_diff_max = std::abs(dxp0-dxp[0])/a_dX[0];
#if CH_SPACEDIM==2
         dxp0 = xpbar[1]-xpold[1];
         dxp[1]  = upbar[dirp[2]]*a_cnormHalfDt;
         Real rel_diff_dir = std::abs(dxp0-dxp[1])/a_dX[1];
         rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
#endif   
          
         if (rel_diff_max >= a_rtol) {
            for(int dir=0; dir<SpaceDim; dir++) xpbar[dir] = xpold[dir] + dxp[dir];
	    position_virt[0] = (xpnew + xpold[0])/2.0;
	    position_virt[1] = ypnew/2.0;
            a_temp_pList.transfer(lit);
         }
         else ++lit;

      }

   }

}

#if CH_SPACEDIM==1
void PicSpeciesUtils::stepNormTransfer_SPH_CAR( List<JustinsParticle>&  a_pList,
                                        List<JustinsParticle>&  a_temp_pList,
                                  const RealVect&               a_dX,
                                  const Real                    a_cnormHalfDt, 
                                  const Real                    a_rtol, 
                                  const bool                    a_reverse )
{
   CH_TIME("PicSpeciesUtils::stepNormTransfer_SPH_CAR()");

   // Check for convergence of rpbar-xpold and vrpbar(rpbar)*cnormDt/2.

   Real drp0, drp1, rel_diff_max;
   Real xpnew, ypnew, zpnew, rpnew;   

   if(a_reverse) { // if converged, temp_pList ==> pList, else update position 
      
      ListIterator<JustinsParticle> lit(a_temp_pList);
      for(lit.begin(); lit.ok();) {
 
         RealVect& rpbar = lit().position();
         const RealVect& xpold = lit().position_old();
         const std::array<Real,3>& upbar = lit().velocity();

	 std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();

         // compute rpnew using updated value of CAR positions	 
         xpnew = xpold[0] + 2.0*a_cnormHalfDt*upbar[0];
         ypnew = 2.0*a_cnormHalfDt*upbar[1];
         zpnew = 2.0*a_cnormHalfDt*upbar[2];
         rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew + zpnew*zpnew);
         
         drp0 = rpbar[0] - xpold[0];
         drp1 = (rpnew - xpold[0])/2.0;
         rel_diff_max = std::abs(drp0-drp1)/a_dX[0];
         if (rel_diff_max < a_rtol) a_pList.transfer(lit);
         else {
            rpbar[0] = (rpnew + xpold[0])/2.0;
	    position_virt[0] = (xpnew + xpold[0])/2.0;
	    position_virt[1] = ypnew/2.0;
	    position_virt[2] = zpnew/2.0;
            ++lit;
         }

      }

   }
   else { // if not converged, update position and send to temp_pList
 
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok();) {
         
         RealVect& rpbar = lit().position();
         const RealVect& xpold = lit().position_old();
         const std::array<Real,3>& upbar = lit().velocity();
         
	 std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
         
         xpnew = xpold[0] + 2.0*a_cnormHalfDt*upbar[0];
         ypnew = 2.0*a_cnormHalfDt*upbar[1];
         zpnew = 2.0*a_cnormHalfDt*upbar[2];
         rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew + zpnew*zpnew);
         
         drp0 = rpbar[0] - xpold[0];
         drp1 = (rpnew - xpold[0])/2.0;
         rel_diff_max = std::abs(drp0-drp1)/a_dX[0];
         if (rel_diff_max >= a_rtol) {
            rpbar[0] = (rpnew + xpold[0])/2.0;
	    position_virt[0] = (xpnew + xpold[0])/2.0;
	    position_virt[1] = ypnew/2.0;
	    position_virt[2] = zpnew/2.0;
            a_temp_pList.transfer(lit);
         }
         else ++lit;

      }

   }

}
#endif

void PicSpeciesUtils::stepNormTransfer( List<JustinsParticle>&  a_pList,
                                        List<JustinsParticle>&  a_temp_pList,
                                  const RealVect&               a_dX,
                                  const Real                    a_cnormHalfDt, 
                                  const Real                    a_rtol, 
                                  const bool                    a_reverse )
{
   CH_TIME("PicSpeciesUtils::stepNormTransfer()");

   // Check for convergence of xpbar-xpold and vpbar(xpbar)*cnormDt/2.
   // Note that vpbar has to be updated after xpbar before coming here.
   // a_reverse = false: if not converged, move particle from pList to temp_pList
   // a_reverse = true: if converged, move particle from temp_pList to pList

   Real dxp0, rel_diff_dir, rel_diff_max;
   RealVect dxp;
#ifdef RELATIVISTIC_PARTICLES
   Real gammap, gbsq_new, gbsq_old;
   std::array<Real,3> upnew;
#endif
      
   if(a_reverse) { // if converged, temp_pList ==> pList, else update position 
      
      ListIterator<JustinsParticle> li(a_temp_pList);
      for(li.begin(); li.ok();) {
 
         RealVect& xpbar = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& upbar = li().velocity();

#ifdef RELATIVISTIC_PARTICLES
         const std::array<Real,3>& upold = li().velocity_old();
         for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
         gbsq_old = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
         gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
         gammap = 0.5*( sqrt(1.0 + gbsq_old) + sqrt(1.0 + gbsq_new) );
#endif         

         rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xpbar[dir]-xpold[dir];
#ifdef RELATIVISTIC_PARTICLES
            dxp[dir]  = upbar[dir]/gammap*a_cnormHalfDt;
#else
            dxp[dir]  = upbar[dir]*a_cnormHalfDt;
#endif
            rel_diff_dir = std::abs(dxp0-dxp[dir])/a_dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max < a_rtol) a_pList.transfer(li);
         else {
            for(int dir=0; dir<SpaceDim; dir++) xpbar[dir] = xpold[dir] + dxp[dir];
            ++li;
         }

      }

   }
   else { // if not converged, update position and send to temp_pList
 
      ListIterator<JustinsParticle> li(a_pList);
      for(li.begin(); li.ok();) {
         
         RealVect& xpbar = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& upbar = li().velocity();
         
#ifdef RELATIVISTIC_PARTICLES
         const std::array<Real,3>& upold = li().velocity_old();
         for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
         gbsq_old = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
         gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
         gammap = 0.5*( sqrt(1.0 + gbsq_old) + sqrt(1.0 + gbsq_new) );
#endif         

         rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xpbar[dir]-xpold[dir];
#ifdef RELATIVISTIC_PARTICLES
            dxp[dir]  = upbar[dir]/gammap*a_cnormHalfDt;
#else
            dxp[dir]  = upbar[dir]*a_cnormHalfDt;
#endif
            rel_diff_dir = std::abs(dxp0-dxp[dir])/a_dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max >= a_rtol) {
            for(int dir=0; dir<SpaceDim; dir++) xpbar[dir] = xpold[dir] + dxp[dir];
            a_temp_pList.transfer(li);
         }
         else ++li;

      }

   }

}

void PicSpeciesUtils::stepNormTransferInflow( List<JustinsParticle>&  a_pList,
                                              List<JustinsParticle>&  a_temp_pList,
                                        const RealVect&               a_dX,
                                        const Real                    a_Xbdry,
                                        const int                     a_bdry_dir,
                                        const Real                    a_cnormDt,
                                        const Real                    a_rtol,
                                        const bool                    a_reverse )
{
   
   Real cnormDt0, cnormDt1;
   Real dxp10, dxp1, rel_diff_dir;
 
   Real gammap_old, gammap;
#ifdef RELATIVISTIC_PARTICLES
   Real gbsq_new;
   std::array<Real,3> upnew;
#endif
     
   RealVect X0;
   X0[a_bdry_dir] = a_Xbdry;

   if(a_reverse) { // if converged, temp_pList ==> pList 
      
      ListIterator<JustinsParticle> li(a_temp_pList);
      for(li.begin(); li.ok();) {
         
         const RealVect& xpbar = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& upbar = li().velocity();
         const std::array<Real,3>& upold = li().velocity_old();

         // put reflected particle back
         if(upbar[a_bdry_dir]==0.0) {
            a_pList.transfer(li);
            continue;
         }
        
         gammap_old = 1.0;
#ifdef RELATIVISTIC_PARTICLES
         gammap_old += upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
         gammap_old = sqrt(gammap_old);
#endif

         cnormDt0 = (X0[a_bdry_dir]-xpold[a_bdry_dir])*gammap_old/upold[a_bdry_dir];
         for(int dir=0; dir<SpaceDim; dir++) {
            if(dir==a_bdry_dir) continue;
            X0[dir] = xpold[dir] + upold[dir]/gammap_old*cnormDt0;
         }
         cnormDt1 = a_cnormDt - cnormDt0;
         
#ifdef RELATIVISTIC_PARTICLES
         for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
         gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
         gammap = 0.5*(gammap_old + sqrt(1.0 + gbsq_new));
#else
         gammap = 1.0;
#endif

         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp10 = 2.0*xpbar[dir]-xpold[dir]-X0[dir];
            dxp1  = upbar[dir]/gammap*cnormDt1;
            rel_diff_dir = std::abs(dxp10-dxp1)/a_dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max < a_rtol) a_pList.transfer(li);
         else ++li;

      }

   }
   else { // if not converged, pList ==> temp_pList
 
      ListIterator<JustinsParticle> li(a_pList);
      for(li.begin(); li.ok();) {

         const RealVect& xpbar = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& upbar = li().velocity();
         const std::array<Real,3>& upold = li().velocity_old();
         
         gammap_old = 1.0;
#ifdef RELATIVISTIC_PARTICLES
         gammap_old += upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
         gammap_old = sqrt(gammap_old);
#endif

         cnormDt0 = (X0[a_bdry_dir]-xpold[a_bdry_dir])*gammap_old/upold[a_bdry_dir];
         for(int dir=0; dir<SpaceDim; dir++) {
            if(dir==a_bdry_dir) continue;
            X0[dir] = xpold[dir] + upold[dir]/gammap_old*cnormDt0;
         }
         cnormDt1 = a_cnormDt - cnormDt0;

#ifdef RELATIVISTIC_PARTICLES
         for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
         gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
         gammap = 0.5*(gammap_old + sqrt(1.0 + gbsq_new));
#else
         gammap = 1.0;
#endif

         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp10 = 2.0*xpbar[dir]-xpold[dir]-X0[dir];
            dxp1  = upbar[dir]/gammap*cnormDt1;
            rel_diff_dir = std::abs(dxp10-dxp1)/a_dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max >= a_rtol) a_temp_pList.transfer(li);
         else ++li;

      }

   }

}

#include "NamespaceFooter.H"

