#include "ParticleUtils.H"

#include "NamespaceHeader.H"

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

#include "NamespaceFooter.H"


