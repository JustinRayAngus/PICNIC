
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
   
   // compute vpbar=(vp^n + vp^{n+1})/2 using Boris method
   // vpbar = vpn + dt/2*q/m*(Ep(xpbar) + vpbar x Bp(xpbar))
   // if byHalfDt is false, convert vpbar to vpnew: vp -> 2*vp - vpold

   // advance velocities using Boris algorithm
   Real bp0, bp1, bp2, denom;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;
   
   Real alpha = a_fnorm*a_cnormDt/2.0;

   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   std::array<int,3> dirp = {0,1,2};
   if(a_anticyclic) {
      dirp[1] = 2;
      dirp[2] = 1;
   }

   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {

      std::array<Real,3>& vp = li().velocity(); // actually beta
      const std::array<Real,3>& vpold = li().velocity_old(); // actually beta
      const std::array<Real,3>& Ep = li().electric_field();
      const std::array<Real,3>& Bp = li().magnetic_field();
    
      // add half acceleration to old velocity
      vm0 = vpold[dirp[0]] + alpha*Ep[dirp[0]];
      vm1 = vpold[dirp[1]] + alpha*Ep[dirp[1]];
      vm2 = vpold[dirp[2]] + alpha*Ep[dirp[2]];
   
      // scale Bp by alpha
      bp0 = alpha*Bp[dirp[0]];
      bp1 = alpha*Bp[dirp[1]];
      bp2 = alpha*Bp[dirp[2]];
      denom = 1.0 + bp0*bp0 + bp1*bp1 + bp2*bp2;

      // define vpr = vm + vm x bp
      vpr0 = vm0 + vm1*bp2 - vm2*bp1;
      vpr1 = vm1 + vm2*bp0 - vm0*bp2;
      vpr2 = vm2 + vm0*bp1 - vm1*bp0;

      // compute vpbar = vm + (vpr x bp)/denom (energy conserving to machine precision)
      vp[dirp[0]] = vm0 + (vpr1*bp2 - vpr2*bp1)/denom;
      vp[dirp[1]] = vm1 + (vpr2*bp0 - vpr0*bp2)/denom;
      vp[dirp[2]] = vm2 + (vpr0*bp1 - vpr1*bp0)/denom;

      // compute vpbar with 2nd formulation (kinda works)
      //vp[dirp[0]] = (vm0*denom + vpr1*bp2 - vpr2*bp1)/denom;
      //vp[dirp[1]] = (vm1*denom + vpr2*bp0 - vpr0*bp2)/denom;
      //vp[dirp[2]] = (vm2*denom + vpr0*bp1 - vpr1*bp0)/denom;

      // compute vpbar with 3rd formulation (not energy conserving to machine precision)
      //Real vmbp = vm0*bp0 + vm1*bp1 + vm2*bp2;
      //vp[dirp[0]] = (vm0 + vm1*bp2 - vm2*bp1 + vmbp*bp0)/denom;
      //vp[dirp[1]] = (vm1 + vm2*bp0 - vm0*bp2 + vmbp*bp1)/denom;
      //vp[dirp[2]] = (vm2 + vm0*bp1 - vm1*bp0 + vmbp*bp2)/denom;
      
      // compute vpbar with 4th formulation (not energy conserving to machine precision)
      //Real vmbp = vm0*bp0 + vm1*bp1 + vm2*bp2;
      //vp[dirp[0]] = vm0 + (vm0*(1.0-denom) + vm1*bp2 - vm2*bp1 + vmbp*bp0)/denom;
      //vp[dirp[1]] = vm1 + (vm1*(1.0-denom) + vm2*bp0 - vm0*bp2 + vmbp*bp1)/denom;
      //vp[dirp[2]] = vm2 + (vm2*(1.0-denom) + vm0*bp1 - vm1*bp0 + vmbp*bp2)/denom;
      
      // compute vpbar with 5th formulation (works)
      //Real vmbp = vm0*bp0 + vm1*bp1 + vm2*bp2;
      //vp[dirp[0]] = vm0 + (-vm0*(bp0*bp0 + bp1*bp1 + bp2*bp2) + vm1*bp2 - vm2*bp1 + vmbp*bp0)/denom;
      //vp[dirp[1]] = vm1 + (-vm1*(bp0*bp0 + bp1*bp1 + bp2*bp2) + vm2*bp0 - vm0*bp2 + vmbp*bp1)/denom;
      //vp[dirp[2]] = vm2 + (-vm2*(bp0*bp0 + bp1*bp1 + bp2*bp2) + vm0*bp1 - vm1*bp0 + vmbp*bp2)/denom;
      
      // compute vpbar with 6th formulation (works)
      //vp[dirp[0]] = vm0 + (-vm0*(bp1*bp1 + bp2*bp2) + vm1*bp2 - vm2*bp1 + (vm1*bp1 + vm2*bp2)*bp0)/denom;
      //vp[dirp[1]] = vm1 + (-vm1*(bp0*bp0 + bp2*bp2) + vm2*bp0 - vm0*bp2 + (vm0*bp0 + vm2*bp2)*bp1)/denom;
      //vp[dirp[2]] = vm2 + (-vm2*(bp0*bp0 + bp1*bp1) + vm0*bp1 - vm1*bp0 + (vm0*bp0 + vm1*bp1)*bp2)/denom;
      
      if(!a_byHalfDt) {
         vp[dirp[0]] = 2.0*vp[dirp[0]] - vpold[dirp[0]];
         vp[dirp[1]] = 2.0*vp[dirp[1]] - vpold[dirp[1]];
         vp[dirp[2]] = 2.0*vp[dirp[2]] - vpold[dirp[2]];
      }

   } // end loop over particle list

}

void 
PicSpeciesUtils::applyForcesAxisymm( List<JustinsParticle>&  a_pList,
                               const Real                    a_fnorm,
                               const Real                    a_cnormDt,
                               const bool                    a_byHalfDt,
                               const bool                    a_anticyclic,
                               const int                     a_iter_max )
{
   if(a_pList.length()==0) return;   
   
   // advance vp from step n to step n+1 using modified Boris algorithm
   // that uses an iterative method to include the inertia term directly.
   // Only for cylindrical so far. In cylindrical, the equations of motion
   // can be expressed same as they are for Cartesian with the transformation
   // Bpz -> Bpz + vpthbar*mp/qp/rp. So one can still use the Boris algorithm,
   // but iterations are needed in order for the nonlinear term to converge.
   //
   // vpr  = vprn  + dt*qp/mp*(Epr  + vpthbar*Bpz - vpzbar*Bpth) + dt*vpthbar^2*mp/rp
   //      = vprn  + dt*qp/mp*(Epr  + vpthbar*Bpz_mod - vpzbar*Bpth)
   // vpth = vpthn + dt*qp/mp*(Epth + vpzbar*Bpr - vprbar*Bpz) - dt*vpthbar*vprbar*mp/rp
   //      = vpthn + dt*qp/mp*(Epth + vpzbar*Bpr - vprbar*Bpz_mod)
   // vpz = vpzn + dt*qp/mp*(Epz + vprbar*Bpth - vpthbar*Bpr)
   // with Bpz_mod = Bpz + vpthbar*mp/qp/rp
   // if byHalfDt is true, convert to vpbar: vp -> (vp + vpold)/2
   //
   // Note that 1) convergence is not needed for exact energy conservervation
   //       and 2) 2nd order accurate after 2 iterations

   // advance velocities using Boris algorithm
   Real t0, t1, t2, denom, s0, s1, s2;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;
   Real vpl0, vpl1, vpl2;

   Real alpha = a_fnorm*a_cnormDt/2.0;

   // iterative method for axisymmetric 
   Real moqr = 0.0; // mp/qp/rp
   Real vth_bar = 0.0;
   Real vth_bar_old, dVth;
   const Real step_tol = 1.0e-12;

   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   std::array<int,3> dirp = {0,1,2};
   if(a_anticyclic) {
      dirp[1] = 2;
      dirp[2] = 1;
   }

   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {

      std::array<Real,3>& vp = li().velocity(); // actually beta
      const std::array<Real,3>& vpold = li().velocity_old(); // actually beta
      const std::array<Real,3>& Ep = li().electric_field();
      const std::array<Real,3>& Bp = li().magnetic_field();
    
      const RealVect& xp = li().position();
      vth_bar = vp[dirp[1]];
      if(std::abs(xp[0])>0.0) moqr = 1.0/a_fnorm/xp[0];

      // add half acceleration to old velocity
      vm0 = vpold[dirp[0]] + alpha*Ep[dirp[0]];
      vm1 = vpold[dirp[1]] + alpha*Ep[dirp[1]];
      vm2 = vpold[dirp[2]] + alpha*Ep[dirp[2]];
   
      int iter = 0;
      while (iter<a_iter_max) {

         // define rotation vectors t and s
         t0 = alpha*Bp[dirp[0]];
         t1 = alpha*Bp[dirp[1]];
         t2 = alpha*(Bp[dirp[2]] + moqr*vth_bar);
         denom = 1.0 + t0*t0 + t1*t1 + t2*t2;
         s0 = 2.0*t0/denom;
         s1 = 2.0*t1/denom;
         s2 = 2.0*t2/denom;

         // define vpr = vm + vm x t
         vpr0 = vm0 + vm1*t2 - vm2*t1;
         vpr1 = vm1 + vm2*t0 - vm0*t2;
         vpr2 = vm2 + vm0*t1 - vm1*t0;

         // rotate (define vplus = vminus + vprime x s)
         vpl0 = vm0 + vpr1*s2 - vpr2*s1;
         vpl1 = vm1 + vpr2*s0 - vpr0*s2;
         vpl2 = vm2 + vpr0*s1 - vpr1*s0;

         // add another half acceleration
         vp[dirp[0]] = vpl0 + alpha*Ep[dirp[0]];
         vp[dirp[1]] = vpl1 + alpha*Ep[dirp[1]];
         vp[dirp[2]] = vpl2 + alpha*Ep[dirp[2]];

         vth_bar_old = vth_bar;
         vth_bar = (vp[dirp[1]] + vpold[dirp[1]])/2.0;
         dVth = std::abs(vth_bar-vth_bar_old);
         if(dVth < step_tol) break;
         //cout << "JRA iter = " << iter << endl;
         //cout << "JRA dVth = " << dVth << endl;
         ++iter;

      }

      if(a_byHalfDt) {
         vp[dirp[0]] = (vp[dirp[0]] + vpold[dirp[0]])/2.0;
         vp[dirp[1]] = (vp[dirp[1]] + vpold[dirp[1]])/2.0;
         vp[dirp[2]] = (vp[dirp[2]] + vpold[dirp[2]])/2.0;
      }

   } // end loop over particle list

}

void PicSpeciesUtils::stepNormTransfer( List<JustinsParticle>&  a_pList,
                                        List<JustinsParticle>&  a_temp_pList,
                                  const RealVect&               a_dX,
                                  const Real                    a_cnormHalfDt, 
                                  const Real                    a_rtol, 
                                  const bool                    a_reverse )
{
   // Check for convergence of xpbar-xpold and vpbar(xpbar)*cnormDt/2.
   // Note that vpbar has to be updated after xpbar before coming here.
   // a_reverse = false: if not converged, move particle from pList to temp_pList
   // a_reverse = true: if converged, move particle from temp_pList to pList

   Real dxp0, dxp, rel_diff_dir;
         
   if(a_reverse) { // if converged, temp_pList ==> pList 
      
      ListIterator<JustinsParticle> li(a_temp_pList);
      for(li.begin(); li.ok();) {
         const RealVect& xpbar = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& vpbar = li().velocity();
         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xpbar[dir]-xpold[dir];
            dxp  = vpbar[dir]*a_cnormHalfDt;
            rel_diff_dir = std::abs(dxp0-dxp)/a_dX[dir];
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
         const std::array<Real,3>& vpbar = li().velocity();
         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xpbar[dir]-xpold[dir];
            dxp  = vpbar[dir]*a_cnormHalfDt;
            rel_diff_dir = std::abs(dxp0-dxp)/a_dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max >= a_rtol) a_temp_pList.transfer(li);
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
         
   RealVect X0;
   X0[a_bdry_dir] = a_Xbdry;

   if(a_reverse) { // if converged, temp_pList ==> pList 
      
      ListIterator<JustinsParticle> li(a_temp_pList);
      for(li.begin(); li.ok();) {
         const RealVect& xpbar = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& vpbar = li().velocity();

         // put reflected particle back
         if(vpbar[a_bdry_dir]==0.0) {
            a_pList.transfer(li);
            continue;
         }
        
         const std::array<Real,3>& vpold = li().velocity_old();
         cnormDt0 = (X0[a_bdry_dir]-xpold[a_bdry_dir])/vpold[a_bdry_dir];
         for(int dir=0; dir<SpaceDim; dir++) {
            if(dir==a_bdry_dir) continue;
            X0[dir] = xpold[dir] + vpold[dir]*cnormDt0;
         }
         cnormDt1 = a_cnormDt - cnormDt0;

         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp10 = 2.0*xpbar[dir]-xpold[dir]-X0[dir];
            dxp1  = vpbar[dir]*cnormDt1;
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
         const std::array<Real,3>& vpbar = li().velocity();
         
         const std::array<Real,3>& vpold = li().velocity_old();
         cnormDt0 = (X0[a_bdry_dir]-xpold[a_bdry_dir])/vpold[a_bdry_dir];
         for(int dir=0; dir<SpaceDim; dir++) {
            if(dir==a_bdry_dir) continue;
            X0[dir] = xpold[dir] + vpold[dir]*cnormDt0;
         }
         cnormDt1 = a_cnormDt - cnormDt0;

         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp10 = 2.0*xpbar[dir]-xpold[dir]-X0[dir];
            dxp1  = vpbar[dir]*cnormDt1;
            rel_diff_dir = std::abs(dxp10-dxp1)/a_dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max >= a_rtol) a_temp_pList.transfer(li);
         else ++li;
      }

   }

}

#include "NamespaceFooter.H"

