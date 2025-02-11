#include "ScatteringUtils.H"

#include "NamespaceHeader.H"

void ScatteringUtils::computeDeltaU( std::array<Real,3>&  a_deltaU, 
                               const std::array<Real,3>&  a_vp1,
                               const std::array<Real,3>&  a_vp2,
                               const Real                 a_costh,
                               const Real                 a_sinth,
                               const Real                 a_phi )
{  

   // compute relative velocity vector
   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];

   // set cos and sin of azimuthal angle
   Real cosphi = cos(a_phi);
   Real sinphi = sin(a_phi);
   
   computeDeltaU(a_deltaU,ux,uy,uz,a_costh,a_sinth,cosphi,sinphi);
   
}

void ScatteringUtils::computeDeltaU( std::array<Real,3>&  a_deltaU,
                               const Real                 ux, 
                               const Real                 uy, 
                               const Real                 uz, 
                               const Real                 a_costh,
                               const Real                 a_sinth,
                               const Real                 a_phi )
{
   // set cos and sin of azimuthal angle
   Real cosphi = cos(a_phi);
   Real sinphi = sin(a_phi);   
   computeDeltaU(a_deltaU,ux,uy,uz,a_costh,a_sinth,cosphi,sinphi);
}

void ScatteringUtils::computeDeltaU( std::array<Real,3>&  a_deltaU, 
                               const std::array<Real,3>&  a_vp1,
                               const std::array<Real,3>&  a_vp2,
                               const Real                 a_potential,
                               const Real                 a_mu,
                               const Real                 a_costh,
                               const Real                 a_sinth,
                               const Real                 a_phi )
{  

   // compute pre-collision relative velocity vector
   Real ux = a_vp1[0]-a_vp2[0];
   Real uy = a_vp1[1]-a_vp2[1];
   Real uz = a_vp1[2]-a_vp2[2];

   // set cos and sin of azimuthal angle
   Real cosphi = cos(a_phi);
   Real sinphi = sin(a_phi);
      
   // define u, uprime, and uperp
   Real u = sqrt(ux*ux + uy*uy + uz*uz);
   Real uprime = u*u - 2.0*a_potential/a_mu;
   CH_assert(uprime>=0.0);
   uprime = sqrt(uprime);
   
   // define deltaU
   Real uperp = sqrt(ux*ux + uy*uy);
   if(uperp==0.0) {
      a_deltaU[0] = a_sinth*cosphi*uprime;
      a_deltaU[1] = a_sinth*sinphi*uprime;   
      a_deltaU[2] = a_costh*uprime - u;
   }
   else {
      a_deltaU[0] = (ux*uz/uperp*a_sinth*cosphi - uy*u/uperp*a_sinth*sinphi + ux*a_costh)*uprime/u - ux;
      a_deltaU[1] = (uy*uz/uperp*a_sinth*cosphi + ux*u/uperp*a_sinth*sinphi + uy*a_costh)*uprime/u - uy;   
      a_deltaU[2] = (-uperp*a_sinth*cosphi + uz*a_costh)*uprime/u - uz;
   }
   
}

void ScatteringUtils::LorentzTransform( Real&                gammapst,
                                        std::array<Real,3>&  upst,
                                  const Real                 gammap,
                                  const std::array<Real,3>&  up,
                                  const Real                 gammacm,
                                  const std::array<Real,3>&  ucm,
                                  const bool                 reverse )
{
    // Relativistic LorentzTransform to st frame (u = gamma*beta)

    long double ucmdotup = ucm[0]*up[0] + ucm[1]*up[1] + ucm[2]*up[2];
    if (reverse) { // same as using -ucm
        gammapst = gammacm*gammap + ucmdotup;
        long double upst_fact = (1.0/(1.0+gammacm)*ucmdotup + gammap);
        for(int n=0; n<3; n++) { upst[n] = up[n] + upst_fact*ucm[n]; }
    }
    else {
        gammapst = gammacm*gammap - ucmdotup;
        long double upst_fact = (1.0/(1.0+gammacm)*ucmdotup - gammap);
        for(int n=0; n<3; n++) { upst[n] = up[n] + upst_fact*ucm[n]; }
    }

}

Real ScatteringUtils::semilogInterp( const std::vector<Real>&  a_X,
                                     const std::vector<Real>&  a_Y,
                                     const Real                a_X0,
                                     const int                 a_index )
{
     
   Real Y0;
   CH_assert(a_index<a_X.size()-1);

   Real log10X0  = log10(a_X0);       
   Real log10Xup = log10(a_X[a_index]);       
   Real log10Xdn = log10(a_X[a_index+1]);       

   Y0 = ( a_Y[a_index+1]*(log10X0 - log10Xdn) 
      +   a_Y[a_index]*(log10Xup - log10X0) )
      / ( log10Xup - log10Xdn );

   return Y0;

}

Real ScatteringUtils::loglogInterp( const std::vector<Real>&  a_X,
                                    const std::vector<Real>&  a_Y,
                                    const Real                a_X0,
                                    const int                 a_index )
{
     
   Real Y0;
   CH_assert(a_index<a_X.size()-1);

   Real log10X0  = log10(a_X0);       
   Real log10Xup = log10(a_X[a_index]);       
   Real log10Xdn = log10(a_X[a_index+1]);       

   Y0 = ( log10(a_Y[a_index+1])*(log10X0 - log10Xdn) 
      +   log10(a_Y[a_index])*(log10Xup - log10X0) )
      / ( log10Xup - log10Xdn );

   return pow(10.0,Y0);

}
   

#include "NamespaceFooter.H"


