#include "ScatteringUtils.H"
//#include "ScatteringUtilsF_F.H"

#include "Constants.H"
#include "MathUtils.H"
#include <cmath>
#include <cstdlib>
#include "CH_Timer.H"

#include "NamespaceHeader.H"
      
void ScatteringUtils::computeDeltaU( std::array<Real,3>&  a_deltaU, 
                               const Real  a_vx1, 
                               const Real  a_vy1,
                               const Real  a_vz1,
                               const Real  a_vx2,
                               const Real  a_vy2,
                               const Real  a_vz2, 
                               const Real  a_th )
{
   
   CH_TIME("ScatteringUtils::computeDeltaU()");

   Real ux = a_vx1-a_vx2;
   Real uy = a_vy1-a_vy2;
   Real uz = a_vz1-a_vz2;
   Real u = sqrt(ux*ux + uy*uy + uz*uz);
   Real uperp = sqrt(ux*ux + uy*uy);

   // set cos and sin of passed polar angle
   Real costh = cos(a_th);
   Real sinth = sin(a_th);
   
   // set random azimuthal angle
   Real phi = Constants::TWOPI*MathUtils::rand();
   Real cosphi = cos(phi);
   Real sinphi = sin(phi);

   // define deltaU
   a_deltaU[0] = ux/uperp*uz*sinth*cosphi - uy*uperp*u*sinth*sinphi - ux*(1.-costh);   
   a_deltaU[1] = uy/uperp*uz*sinth*cosphi + ux*uperp*u*sinth*sinphi - uy*(1.-costh);   
   a_deltaU[2] = -uperp*sinth*cosphi - uz*(1.-costh);   

}
