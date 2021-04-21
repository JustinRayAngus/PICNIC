#include "ParticleUtils.H"
#include "ParticleUtilsF_F.H"

#include "NamespaceHeader.H"
      
void ParticleUtils::borisPusher( std::array<Real,3>&  a_vp,
                           const std::array<Real,3>&  a_vpold,
                           const std::array<Real,3>&  a_Ep,
                           const std::array<Real,3>&  a_Bp,
                           const Real&                a_fnorm_const,
                           const Real&                a_cnormDt ) 
{
   CH_TIME("ParticleUtils::borisPusher()");

   const Real zeroValue = 0.0;
   FORT_BORIS_PUSHER( CHF_REAL(a_vp[0]),
                      CHF_REAL(a_vp[1]),
                      CHF_REAL(a_vp[2]),
                      CHF_CONST_REAL(a_vpold[0]),
                      CHF_CONST_REAL(a_vpold[1]),
                      CHF_CONST_REAL(a_vpold[2]),
                      CHF_CONST_REAL(a_Ep[0]),
                      CHF_CONST_REAL(a_Ep[1]),
                      CHF_CONST_REAL(a_Ep[2]),
                      CHF_CONST_REAL(a_Bp[0]),
                      CHF_CONST_REAL(a_Bp[1]),
                      CHF_CONST_REAL(a_Bp[2]),
                      CHF_CONST_REAL(a_fnorm_const),
                      CHF_CONST_REAL(a_cnormDt) );

}

void ParticleUtils::computeDeltaU( std::array<Real,3>&  a_deltaU, 
                             const std::array<Real,3>&  a_vp1,
                             const std::array<Real,3>&  a_vp2,
                             const Real&                a_theta,
                             const Real&                a_phi )
{  
   CH_TIME("ParticleUtils::computeDeltaU() v1");

   // compute relative velocity vector
   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];

   // set cos and sin of polar angle
   Real costh = cos(a_theta);
   Real sinth = sin(a_theta);
   
   // set cos and sin of azimuthal angle
   Real cosphi = cos(a_phi);
   Real sinphi = sin(a_phi);
   
   computeDeltaU(a_deltaU,ux,uy,uz,costh,sinth,cosphi,sinphi);
   
}

void ParticleUtils::computeDeltaU( std::array<Real,3>&  a_deltaU,
                             const Real&                ux, 
                             const Real&                uy, 
                             const Real&                uz, 
                             const Real&                a_theta,
                             const Real&                a_phi )
{
   //CH_TIME("ParticleUtils::computeDeltaU() v2");

   // set cos and sin of polar angle
   Real costh = cos(a_theta);
   Real sinth = sin(a_theta);
   
   // set cos and sin of azimuthal angle
   Real cosphi = cos(a_phi);
   Real sinphi = sin(a_phi);
   
   computeDeltaU(a_deltaU,ux,uy,uz,costh,sinth,cosphi,sinphi);

}

void ParticleUtils::computeDeltaU( std::array<Real,3>&  a_deltaU,
                             const Real&                ux, 
                             const Real&                uy, 
                             const Real&                uz, 
                             const Real&                costh,
                             const Real&                sinth,
                             const Real&                cosphi,
                             const Real&                sinphi )
{
   CH_TIME("ParticleUtils::computeDeltaU() v3");

   // define u and uperp
   Real u = sqrt(ux*ux + uy*uy + uz*uz);
   Real uperp = sqrt(ux*ux + uy*uy);

   // define deltaU
   a_deltaU[0] = ux*uz/uperp*sinth*cosphi - uy*u/uperp*sinth*sinphi - ux*(1.-costh);   
   a_deltaU[1] = uy*uz/uperp*sinth*cosphi + ux*u/uperp*sinth*sinphi - uy*(1.-costh);   
   a_deltaU[2] = -uperp*sinth*cosphi - uz*(1.-costh);   

   /*
   FORT_COMPUTE_DELTAU( CHF_REAL(a_deltaU[0]),
                        CHF_REAL(a_deltaU[1]),
                        CHF_REAL(a_deltaU[2]),
                        CHF_CONST_REAL(ux),
                        CHF_CONST_REAL(uy),
                        CHF_CONST_REAL(uz),
                        CHF_CONST_REAL(costh),
                        CHF_CONST_REAL(sinth),
                        CHF_CONST_REAL(cosphi),
                        CHF_CONST_REAL(sinphi) );
   */

}

#include "NamespaceFooter.H"


