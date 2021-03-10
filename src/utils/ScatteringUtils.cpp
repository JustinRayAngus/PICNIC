#include "ScatteringUtils.H"
//#include "ScatteringUtilsF_F.H"

#include "Constants.H"
#include "MathUtils.H"
#include <cmath>
#include <cstdlib>
#include "CH_Timer.H"

#include "NamespaceHeader.H"
      
void ScatteringUtils::computeDeltaU( std::array<Real,3>&  a_deltaU, 
                               const std::array<Real,3>&  a_vp1,
                               const std::array<Real,3>&  a_vp2,
                               const Real&                a_theta )
{
   
   CH_TIME("ScatteringUtils::computeDeltaU()");

   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];
   Real u = sqrt(ux*ux + uy*uy + uz*uz);
   Real uperp = sqrt(ux*ux + uy*uy);

   // set cos and sin of passed polar angle
   Real costh = cos(a_theta);
   Real sinth = sin(a_theta);
   
   // set random azimuthal angle
   Real phi = Constants::TWOPI*MathUtils::rand();
   Real cosphi = cos(phi);
   Real sinphi = sin(phi);

   // define deltaU
   a_deltaU[0] = ux*uz/uperp*sinth*cosphi - uy*u/uperp*sinth*sinphi - ux*(1.-costh);   
   a_deltaU[1] = uy*uz/uperp*sinth*cosphi + ux*u/uperp*sinth*sinphi - uy*(1.-costh);   
   a_deltaU[2] = -uperp*sinth*cosphi - uz*(1.-costh);   
   
   /*   
   Real OneMiCosth = 1.-costh;
   Real sinthOvUperp = sinth/uperp;

   // define deltaU
   a_deltaU[0] = ux*uz*sinthOvUperp*cosphi - uy*u*sinthOvUperp*sinphi - ux*OneMiCosth;   
   a_deltaU[1] = uy*uz*sinthOvUperp*cosphi + ux*u*sinthOvUperp*sinphi - uy*OneMiCosth;   
   a_deltaU[2] = -uperp*sinth*cosphi - uz*OneMiCosth;   
   */

}

void ScatteringUtils::computeDeltaU( std::array<Real,3>&  a_deltaU,
                               const Real&                ux, 
                               const Real&                uy, 
                               const Real&                uz, 
                               const Real&                a_theta )
{
   
   CH_TIME("ScatteringUtils::computeDeltaU()");

   Real u = sqrt(ux*ux + uy*uy + uz*uz);
   Real uperp = sqrt(ux*ux + uy*uy);

   // set cos and sin of passed polar angle
   Real costh = cos(a_theta);
   Real sinth = sin(a_theta);
   
   // set random azimuthal angle
   Real phi = Constants::TWOPI*MathUtils::rand();
   Real cosphi = cos(phi);
   Real sinphi = sin(phi);

   // define deltaU
   a_deltaU[0] = ux*uz/uperp*sinth*cosphi - uy*u/uperp*sinth*sinphi - ux*(1.-costh);   
   a_deltaU[1] = uy*uz/uperp*sinth*cosphi + ux*u/uperp*sinth*sinphi - uy*(1.-costh);   
   a_deltaU[2] = -uperp*sinth*cosphi - uz*(1.-costh);   

}


