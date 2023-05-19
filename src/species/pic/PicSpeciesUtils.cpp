
#include "PicSpeciesUtils.H"
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
PicSpeciesUtils::applyForcesAxisymm( List<JustinsParticle>&  a_pList,
                               const Real                    a_fnorm,
                               const Real                    a_cnormDt,
                               const bool                    a_byHalfDt,
                               const bool                    a_anticyclic )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForcesAxisymm()");
   
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
      if(dtheta==0.0) {
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
      
      if(set_dtheta) {
         Real rpbar = xpold[dirp[0]] + a_cnormDt/2.0*up[dirp[0]];
         dtheta = a_cnormDt/2.0*up[dirp[1]]/rpbar/gammap;
         std::array<Real,2>& position_virt = lit().position_virt();
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
PicSpeciesUtils::applyForcesHYBRID( List<JustinsParticle>&  a_pList,
                              const Real                    a_fnorm,
                              const Real                    a_cnormDt,
                              const bool                    a_anticyclic )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpeciesUtils::applyForcesHYBRID()");
   
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
      
      const Real upoldr_2 = costhp*upold[dirp[0]] + sinthp*upold[dirp[1]]; 
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

void PicSpeciesUtils::stepNormTransferCYLCAR( List<JustinsParticle>&  a_pList,
                                        List<JustinsParticle>&  a_temp_pList,
                                  const RealVect&               a_dX,
                                  const Real                    a_cnormHalfDt, 
                                  const Real                    a_rtol, 
                                  const bool                    a_reverse )
{
   CH_TIME("PicSpeciesUtils::stepNormTransferCYLCAR()");

   // Check for convergence of xpbar-xpold and vpbar(xpbar)*cnormDt/2.
   // Note that vpbar has to be updated after xpbar before coming here.
   // a_reverse = false: if not converged, move particle from pList to temp_pList
   // a_reverse = true: if converged, move particle from temp_pList to pList

   Real dxp0, rel_diff_max;
   RealVect dxp;
      
   if(a_reverse) { // if converged, temp_pList ==> pList, else update position 
      
      ListIterator<JustinsParticle> lit(a_temp_pList);
      for(lit.begin(); lit.ok();) {
 
         RealVect& xpbar = lit().position();
         const RealVect& xpold = lit().position_old();
         const std::array<Real,2>& position_virt = lit().position_virt();
         
         const Real costhp = position_virt[0];
         const Real sinthp = position_virt[1];

         const Real xpbar0 = costhp*xpbar[0];
         const Real ypbar0 = sinthp*xpbar[0];
         
         const Real xpnew = 2.0*xpbar0 - xpold[0];
         const Real ypnew = 2.0*ypbar0;
         const Real rpnew = sqrt(xpnew*xpnew + ypnew*ypnew);
 
         dxp0 = xpbar[0]-xpold[0];
         dxp[0] = (rpnew + xpold[0])/2.0 - xpold[0];
         rel_diff_max = std::abs(dxp0-dxp[0])/a_dX[0];
#if CH_SPACEDIM==2
         dxp0 = xpbar[1]-xpold[1];
         const std::array<Real,3>& upbar = lit().velocity();
         dxp[1]  = upbar[1]*a_cnormHalfDt;
         Real rel_diff_dir = std::abs(dxp0-dxp[1])/a_dX[1];
         rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
#endif             

         if (rel_diff_max < a_rtol) a_pList.transfer(lit);
         else {
            for(int dir=0; dir<SpaceDim; dir++) xpbar[dir] = xpold[dir] + dxp[dir];
            ++lit;
         }

      }

   }
   else { // if not converged, update position and send to temp_pList
 
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok();) {
         
         RealVect& xpbar = lit().position();
         const RealVect& xpold = lit().position_old();
         const std::array<Real,2>& position_virt = lit().position_virt();
         
         const Real costhp = position_virt[0];
         const Real sinthp = position_virt[1];

         const Real xpbar0 = costhp*xpbar[0];
         const Real ypbar0 = sinthp*xpbar[0];
         
         const Real xpnew = 2.0*xpbar0 - xpold[0];
         const Real ypnew = 2.0*ypbar0;
         const Real rpnew = sqrt(xpnew*xpnew + ypnew*ypnew);
 
         dxp0 = xpbar[0]-xpold[0];
         dxp[0] = (rpnew + xpold[0])/2.0 - xpold[0];
         rel_diff_max = std::abs(dxp0-dxp[0])/a_dX[0];
#if CH_SPACEDIM==2
         dxp0 = xpbar[1]-xpold[1];
         const std::array<Real,3>& upbar = lit().velocity();
         dxp[1]  = upbar[1]*a_cnormHalfDt;
         Real rel_diff_dir = std::abs(dxp0-dxp[1])/a_dX[1];
         rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
#endif   
          
         if (rel_diff_max >= a_rtol) {
            for(int dir=0; dir<SpaceDim; dir++) xpbar[dir] = xpold[dir] + dxp[dir];
            a_temp_pList.transfer(lit);
         }
         else ++lit;

      }

   }

}

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

