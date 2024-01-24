
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

#include "NamespaceFooter.H"

