#include "PicSpecies.H"
#include <cmath>
#include "PicnicConstants.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"
#include "BinFabFactory.H"

#include "CH_HDF5.H"
#include "ParticleIO.H"

#include "MathUtils.H"
#include "SpaceUtils.H"
#include "PicSpeciesUtils.H"

#include "NamespaceHeader.H"

PicSpecies::PicSpecies( ParmParse&         a_ppspc,
                        const int          a_species,
                        const string&      a_name,
                        const DomainGrid&  a_mesh )
   : m_species(a_species),
     m_mass(),
     m_charge(),
     m_Uint(0.0),
     m_push_type(PLANAR),
     m_species_bc(NULL),
     m_interpRhoToGrid(CIC),
     m_interpJToGrid(CIC),
     m_interpEToParts(CIC),
     m_motion(true),
     m_forces(true),
     m_scatter(false),
     m_write_all_part_comps(false),
     m_charge_conserving_deposit(false),
     m_mesh(a_mesh),
     m_meshInterp(nullptr)
{
   m_name = a_name;

   std::string push_type_str = "PLANAR";
   if(m_mesh.axisymmetric()) {
   
      const string& geom_type = m_mesh.geomType();

      if(geom_type=="sph_R") { // 1D spherical
         push_type_str = "SPH_SPH";
         a_ppspc.query("push_type",push_type_str);
         if(push_type_str=="SPH_SPH") {m_push_type = SPH_SPH;}
         else if(push_type_str=="SPH_HYB") {m_push_type = SPH_HYB;}
         else if(push_type_str=="SPH_CAR") {m_push_type = SPH_CAR;}
         else {
            if(!procID()) { 
               cout << "EXIT FAILURE!!!" << endl;
               cout << "push_type = " << push_type_str << " for species = " << m_species
                    << " is not a valid push_type" << endl;
               cout << "Valid options for " << geom_type << " geometry are "
                    << "SPH_SPH, SPH_HYB, SPH_CAR" << endl;
            }
            exit(EXIT_FAILURE);
         }
      }
      else { // 1D or 2D cylindrical
         push_type_str = "CYL_CYL";
         a_ppspc.query("push_type",push_type_str);
         if(push_type_str=="CYL_BORIS") {m_push_type = CYL_BORIS;}
         else if(push_type_str=="CYL_CYL") {m_push_type = CYL_CYL;}
         else if(push_type_str=="CYL_HYBRID") {m_push_type = CYL_HYB;} // legacy
         else if(push_type_str=="CYL_HYB") {m_push_type = CYL_HYB;}
         else if(push_type_str=="CYL_CAR") {m_push_type = CYL_CAR;}
         else {
            if(!procID()) { 
               cout << "EXIT FAILURE!!!" << endl;
               cout << "push_type = " << push_type_str << " for species = " << m_species 
                    << " is not a valid push_type" << endl;
               cout << "Valid options for " << geom_type << " geometry are "
                    << "CYL_BORIS, CYL_CYL, CYL_HYB, CYL_CAR" << endl;
            }
            exit(EXIT_FAILURE);
         }
      }

   }

   createMeshInterp();
  
   std::string interp_type_N="CIC";
   a_ppspc.query( "interp_type_N", interp_type_N );
   //if(interp_type_N=="NGP") m_interpRhoToGrid = NGP;
   if(interp_type_N=="CIC") m_interpRhoToGrid = CIC;
   if(interp_type_N=="TSC") m_interpRhoToGrid = TSC;

   std::string interp_type_J="CIC";
   a_ppspc.query( "interp_type_J", interp_type_J );
   if(interp_type_J=="CIC") m_interpJToGrid = CIC;
   if(interp_type_J=="TSC") m_interpJToGrid = TSC;
   if(interp_type_J=="CC0") m_interpJToGrid = CC0;
   if(interp_type_J=="CC1") m_interpJToGrid = CC1;
   
   std::string interp_type_E="CIC";
   a_ppspc.query( "interp_type_E", interp_type_E );
   if(interp_type_E=="CIC") m_interpEToParts = CIC;
   if(interp_type_E=="TSC") m_interpEToParts = TSC;
   if(interp_type_E=="CC0") m_interpEToParts = CC0;
   if(interp_type_E=="CC1") m_interpEToParts = CC1;
 
   if(m_interpJToGrid==CC0 || m_interpJToGrid==CC1) m_charge_conserving_deposit = true;

   // flags for energy/charge conservation at physical boundaries
   m_interp_bc_check = false;
   m_suborbit_inflowJ = false;
   m_use_suborbit_model = false;
   m_suborbit_fast_particles = false;
   if(m_charge_conserving_deposit) m_interp_bc_check = true;
   a_ppspc.query( "interp_bc_check", m_interp_bc_check );
   a_ppspc.query( "suborbit_inflow_J", m_suborbit_inflowJ );
   a_ppspc.query( "use_suborbit_model", m_use_suborbit_model );
   if(m_use_suborbit_model) a_ppspc.query("suborbit_testing",m_suborbit_testing);
   a_ppspc.query( "suborbit_fast_particles", m_suborbit_fast_particles );
   if(m_suborbit_fast_particles) {  m_use_suborbit_model = true; }

   a_ppspc.get( "mass", m_mass );
   a_ppspc.get( "charge", m_charge ); 
   a_ppspc.query( "potential", m_Uint );
   a_ppspc.query( "motion", m_motion );
   a_ppspc.query( "forces", m_forces );
   a_ppspc.query( "scatter", m_scatter );
   a_ppspc.query( "write_all_comps", m_write_all_part_comps );

   if (!m_mesh.axisymmetric()) {
      if (m_charge==0) { m_forces = false; }
   }

   if (!procID()) {
      cout << "  name = " << m_name << endl;
      cout << "  number = " << m_species << endl;
      cout << "  mass/me   = " << m_mass << endl;
      cout << "  charge/qe = " << m_charge << endl;
      cout << "  Uint [eV] = " << m_Uint << endl;
      cout << "  interpolate N scheme = " << interp_type_N << endl;
      cout << "  interpolate J scheme = " << interp_type_J << endl;
      cout << "  interpolate E/B scheme = " << interp_type_E << endl;
      cout << "  interpolate bc check = " << (m_interp_bc_check?"true":"false") << endl;
      cout << "  suborbit inflow J    = " << (m_suborbit_inflowJ?"true":"false") << endl;
      cout << "  use suborbit model   = " << (m_use_suborbit_model?"true":"false") << endl;
      cout << "  suborbit fast parts  = " << (m_suborbit_fast_particles?"true":"false") << endl;
      if (m_suborbit_testing) {
         cout << "  suborbit testing = true (fixed suborbit iter_max = 10)" << endl;
      }
      cout << "  push type = " << push_type_str << endl;
      cout << "  motion = " << m_motion << endl;
      cout << "  forces = " << m_forces << endl;
      cout << "  scatter = " << m_scatter << endl << endl;
   }

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   const int ghosts(m_mesh.ghosts());

   m_stable_dt = meshSpacing[0];

   // initialize the member LevelDatas
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   m_Nppc.define(grids,1,ghostVect);
   m_density.define(grids,1,ghostVect);
   m_momentum.define(grids,3,ghostVect);
   m_energy.define(grids,3,ghostVect); // direction-dependent energy
   m_energyOffDiag.define(grids,3,ghostVect);
   m_energyFlux.define(grids,3,ghostVect);

   m_currentDensity.define(grids,1,ghostVect);
   m_inflowJ.define(grids,1,ghostVect);
   m_suborbitJ.define(grids,1,ghostVect);
#if CH_SPACEDIM<3
   m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect);
   m_inflowJ_virtual.define(grids,3-SpaceDim,ghostVect);
   m_suborbitJ_virtual.define(grids,3-SpaceDim,ghostVect);
#endif
   m_chargeDensity.define(grids,1,ghostVect);
   m_chargeDensity_faces.define(grids,1,ghostVect);
   m_chargeDensity_nodes.define(grids,1,ghostVect);
   m_surfaceCharge_nodes.define(grids,1,ghostVect);

   m_temperature.define(grids,3,ghostVect);
   m_velocity.define(grids,3,ghostVect);
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_Nppc[dit].setVal(0.0);
      m_density[dit].setVal(0.0);
      m_momentum[dit].setVal(0.0);
      m_energy[dit].setVal(0.0);
      m_energyOffDiag[dit].setVal(0.0);
      m_energyFlux[dit].setVal(0.0);
      m_temperature[dit].setVal(0.0);
      m_velocity[dit].setVal(0.0);
      m_surfaceCharge_nodes[dit].getFab().setVal(0.0);
   } 

   // each box has to be square with fixedBoxSize length to use ParticleData()
   int fixedBoxSize;
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      Box thisbox( grids[dit] ); 
      fixedBoxSize = thisbox.size(0);
      for (int dir=1; dir<SpaceDim; ++dir) CH_assert(thisbox.size(dir)==fixedBoxSize);
      break;
   } 
   
   // ParticleData<T> behaves similar to LevelData<FArrayBox>
   // ghosts?
   m_data.define(grids, domain, fixedBoxSize,
                 meshSpacing, meshOrigin);
   
   m_data_suborbit.define(grids, domain, fixedBoxSize,
                          meshSpacing, meshOrigin);
  
   // In order to initialize a LevelData<BinFab<T>>, first need to define a
   // BinFabFactory. The factor has a "create" function to define "new" pointers
   // to the BinFab that lives at the box level... Not sure if I need to use this
   // or If I can just do what I'm doing below....
   //
   BinFabFactory<JustinsParticlePtr> bfptrFactory(meshSpacing, meshOrigin);
   m_data_binfab_ptr.define(grids, 1, 0*IntVect::Unit, bfptrFactory);

}

PicSpecies::~PicSpecies()
{
   if(m_species_bc!=NULL) {
      delete m_species_bc;
      m_species_bc = NULL;
   }
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
      m_meshInterp = NULL;
   }
}

void PicSpecies::applyForces( List<JustinsParticle>&  a_pList,
                        const Real   a_cnormDt, 
                        const bool   a_byHalfDt )
{
   CH_TIME("PicSpecies::applyForces()");
    
   if(a_pList.length()==0) return;   

   if(m_push_type==CYL_CAR) { 
      
      // transform Ep and Bp from CYL to CAR
      
      std::array<int,3> dirp = {0,1,2};
      if(m_mesh.anticyclic()) {
         dirp[1] = 2;
         dirp[2] = 1;
      }
   
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         
         std::array<Real,3>& Ep = lit().electric_field();
         std::array<Real,3>& Bp = lit().magnetic_field();
         
         const Real xpbar = lit().position_virt(0);
         const Real ypbar = lit().position_virt(1);
         const Real rpbar = lit().position(0); // (rpnew + rpold)/2.0

         const Real costhp = xpbar/rpbar;
         const Real sinthp = ypbar/rpbar;

         Real Epr  = Ep[dirp[0]];
         Real Epth = Ep[dirp[1]];
         Real Bpr  = Bp[dirp[0]];
         Real Bpth = Bp[dirp[1]];
         Ep[dirp[0]] = costhp*Epr - sinthp*Epth; 
         Ep[dirp[1]] = sinthp*Epr + costhp*Epth; 
         Bp[dirp[0]] = costhp*Bpr - sinthp*Bpth; 
         Bp[dirp[1]] = sinthp*Bpr + costhp*Bpth; 
      
      }

   }
#if CH_SPACEDIM==1
   else if(m_push_type==SPH_CAR) { 
      
      // transform Ep and Bp from SPH to CAR
      
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         
         std::array<Real,3>& Ep = lit().electric_field();
         std::array<Real,3>& Bp = lit().magnetic_field();
         
         const Real xpbar = lit().position_virt(0);
         const Real ypbar = lit().position_virt(1);
         const Real zpbar = lit().position_virt(2);
         
         // compute sines and cosines
         const Real rpbar  = lit().position(0); // (rpnew + rpold)/2.0
         // using below for rpbar is not energy/charge conserving. Also it
         // is 2nd order accurate for single particle pusher test with Er=const
         // Using above gives 4th order. I don't know why...
         //const Real rpbar = std::sqrt(xpbar*xpbar + ypbar*ypbar + zpbar*zpbar);
         const Real rppol = std::sqrt(xpbar*xpbar + ypbar*ypbar);
         const Real costhp = xpbar/rppol;
         const Real sinthp = ypbar/rppol;
         const Real cosphp = rppol/rpbar;
         const Real sinphp = zpbar/rpbar;

         // tranform Ep and Bp from SPH to CAR
         const Real Epr  = Ep[0];
         const Real Epth = Ep[1];
         const Real Epph = Ep[2];
         const Real Bpr  = Bp[0];
         const Real Bpth = Bp[1];
         const Real Bpph = Bp[2];
         Ep[0] = costhp*(cosphp*Epr - sinphp*Epph) - sinthp*Epth; 
         Ep[1] = sinthp*(cosphp*Epr - sinphp*Epph) + costhp*Epth; 
         Ep[2] = sinphp*Epr + cosphp*Epph; 
         Bp[0] = costhp*(cosphp*Bpr - sinphp*Bpph) - sinthp*Bpth; 
         Bp[1] = sinthp*(cosphp*Bpr - sinphp*Bpph) + costhp*Bpth; 
         Bp[2] = sinphp*Bpr + cosphp*Bpph; 
      
      }

   }
#endif
   
   if(m_push_type==CYL_CYL) {
      PicSpeciesUtils::applyForces_CYL_CYL( a_pList, m_fnorm_const, a_cnormDt,
                                           a_byHalfDt, m_mesh.anticyclic() );
   }
   else if(m_push_type==SPH_SPH) {
      PicSpeciesUtils::applyForces_SPH_SPH( a_pList, m_fnorm_const, a_cnormDt,
                                            a_byHalfDt );
   }
   else if(m_push_type==CYL_HYB) {
      PicSpeciesUtils::applyForces_CYL_HYB( a_pList, m_fnorm_const, a_cnormDt,
                                          m_mesh.anticyclic() );
   }
   else if(m_push_type==SPH_HYB) {
      PicSpeciesUtils::applyForces_SPH_HYB( a_pList, m_fnorm_const, a_cnormDt );
   }
   else {
      PicSpeciesUtils::applyForces( a_pList, m_fnorm_const, a_cnormDt, 
                                    a_byHalfDt, m_mesh.anticyclic() ); 
   }
   
   if(m_push_type==CYL_HYB) { // compute thpbar and convert upbar from CAR to CYL
      
      std::array<int,3> dirp = {0,1,2};
      if(m_mesh.anticyclic()) {
         dirp[1] = 2;
         dirp[2] = 1;
      }

      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         
         std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
         Real thpbar = position_virt[0];
         //if(thpbar!=0.0) break;
         if(thpbar==0.0) {
            const RealVect& xpold = lit().position_old();
            std::array<Real,3>& upbar = lit().velocity();
            
            Real costhp = std::cos(thpbar);
            Real sinthp = std::sin(thpbar);

            // convert upbar from CYL to CAR
            const Real upr  = upbar[0];
            const Real upth = upbar[dirp[1]];
            const Real upx  = costhp*upr - sinthp*upth;
            const Real upy  = sinthp*upr + costhp*upth;

            // update thpbar
            const Real xpbar = xpold[0] + a_cnormDt/2.0*upx;
            const Real ypbar = a_cnormDt/2.0*upy;
            thpbar = std::atan2(ypbar,xpbar);
            position_virt[0] = thpbar;

            // update cosines and sines
            costhp = std::cos(thpbar);
            sinthp = std::sin(thpbar);

            // convert upbar from CAR to CYL
            upbar[dirp[0]] =  costhp*upx + sinthp*upy; 
            upbar[dirp[1]] = -sinthp*upx + costhp*upy;
         }

      }

   }
   else if(m_push_type==SPH_HYB) { // compute thpbar, phpbar, and convert upbar from CAR to SPH
      
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         
         std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
         Real thpbar = position_virt[0];
         Real phpbar = position_virt[1];
         if(thpbar==0.0 && phpbar==0.0) {
            const RealVect& xpold = lit().position_old();
            std::array<Real,3>& upbar = lit().velocity();
            
            Real costhp = std::cos(thpbar);
            Real sinthp = std::sin(thpbar);
            Real cosphp = std::cos(phpbar);
            Real sinphp = std::sin(phpbar);
           
            // convert upbar from SPH to CAR
            const Real upr  = upbar[0];
            const Real upth = upbar[1];
            const Real upph = upbar[2];
            const Real upx  =  cosphp*costhp*upr - sinthp*upth - sinphp*costhp*upph; 
            const Real upy  =  cosphp*sinthp*upr + costhp*upth - sinphp*sinthp*upph; 
            const Real upz  =  sinphp*upr + cosphp*upph; 

            // update phpbar and thpbar
            const Real xpbar = xpold[0] + a_cnormDt/2.0*upx;
            const Real ypbar = a_cnormDt/2.0*upy;
            const Real zpbar = a_cnormDt/2.0*upz;
            const Real rpbar = std::sqrt(xpbar*xpbar + ypbar*ypbar + zpbar*zpbar);
            if(rpbar>0.0) phpbar = std::asin(zpbar/rpbar);
            else phpbar = 0.0;
            thpbar = std::atan2(ypbar,xpbar);
            position_virt[0] = thpbar;
            position_virt[1] = phpbar;

            // update cosines and sines
            costhp = std::cos(thpbar);
            sinthp = std::sin(thpbar);
            cosphp = std::cos(phpbar);
            sinphp = std::sin(phpbar);

            // convert upbar from CAR to SPH
            upbar[0] =  cosphp*(costhp*upx + sinthp*upy) + sinphp*upz; 
            upbar[1] = -sinthp*upx + costhp*upy;
            upbar[2] = -sinphp*(costhp*upx + sinthp*upy) + cosphp*upz;
         }

      }

   }

}

void PicSpecies::advancePositionsExplicit( const Real  a_full_dt,
                                           const bool  a_half_step )
{
   if(!m_motion) return;
   CH_TIME("PicSpecies::advancePositionsExplicit()");
   
   const Real cnormDt = m_cvac_norm*a_full_dt;
   Real dt_factor = 1.0;
   if(a_half_step) dt_factor = 0.5; 

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      advancePositionsExplicit( pList, cnormDt*dt_factor );
   }

}

void PicSpecies::advancePositionsExplicit( List<JustinsParticle>&  a_pList,
                                     const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositionsExplicit() from particles");
   
   // xp^{n+1} = xp^{n} + vp*cnormDt;
   // vpbar = upbar/gamma, gamma = sqrt(1 + upbar^2)

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      RealVect& xp = lit().position();
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& up = lit().velocity();
#ifdef RELATIVISTIC_PARTICLES
      Real gammap = sqrt(1.0 + up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
      for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + up[dir]/gammap*a_cnormDt;
#else
      for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + up[dir]*a_cnormDt;
#endif
   }

}

void PicSpecies::advancePositionsImplicit( const Real  a_full_dt )
{
   if(!m_motion) return;
   CH_TIME("PicSpecies::advancePositionsImplicit()");
  
   const Real cnormDt = m_cvac_norm*a_full_dt;

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      advancePositionsImplicit( pList, cnormDt );
   }

}

void PicSpecies::advancePositionsImplicit( List<JustinsParticle>&  a_pList,
                                     const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositionsImplicit() from particles");
   
   // advance particle positions by half dt

   // xp^{n+1} = xp^{n} + vpbar*cnormDt;
   // vpbar = (up^{n+1}+up^{n}/(gamma^{n+1} + gamma^{n})
   // gamma^{n+1} = sqrt(1 + up^{n+1})
   // gamma^{n} = sqrt(1 + up^{n})
   // xpbar = (xp^{n+1} + xp^{n})/2.0

#ifdef RELATIVISTIC_PARTICLES
   Real gammap, gbsq_old, gbsq_new;
   std::array<Real,3> upnew;
#endif
      
   if(m_push_type==CYL_CAR) { advancePositions_CYL_CAR( a_pList, a_cnormDt ); }
#if CH_SPACEDIM==1
   else if(m_push_type==SPH_CAR) { advancePositions_SPH_CAR( a_pList, a_cnormDt ); }
#endif
   else {
   
      const Real cnormHalfDt = a_cnormDt*0.5;
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         RealVect& xp = lit().position();
         const RealVect& xpold = lit().position_old();
         const std::array<Real,3>& upbar = lit().velocity();
#ifdef RELATIVISTIC_PARTICLES
         const std::array<Real,3>& upold = lit().velocity_old();
         for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
         gbsq_old = upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
         gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
         gammap = 0.5*(sqrt(1.0 + gbsq_old) + sqrt(1.0 + gbsq_new));
         for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + upbar[dir]/gammap*cnormHalfDt;
#else
         for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + upbar[dir]*cnormHalfDt;
#endif
      }

   }

}

void PicSpecies::advancePositions_CYL_CAR( List<JustinsParticle>&  a_pList,
                                     const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositions_CYL_CAR() from particles");
   
   // advance particle positions by half dt using CAR method in CYL geometry

   std::array<int,3> dirp = {0,1,2};
   if(m_mesh.anticyclic()) {
      dirp[1] = 2;
      dirp[2] = 1;
   }
   
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      
      RealVect& xpbar = lit().position();
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upbar = lit().velocity();
      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
         
      const Real xpnew = xpold[0] + upbar[0]*a_cnormDt;
      const Real ypnew = upbar[dirp[1]]*a_cnormDt;
      const Real rpnew = sqrt(xpnew*xpnew + ypnew*ypnew);

      // update CAR xpbar and ypbar      
      position_virt[0] = (xpnew + xpold[0])/2.0;
      position_virt[1] = ypnew/2.0;

      // update CYL rpbar and zpbar
      xpbar[0] = (rpnew + xpold[0])/2.0;
#if CH_SPACEDIM==2
      xpbar[1] = xpold[1] + upbar[dirp[2]]*a_cnormDt/2.0;
#endif

   }

}

#if CH_SPACEDIM==1
void PicSpecies::advancePositions_SPH_CAR( List<JustinsParticle>&  a_pList,
                                     const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositions_SPH_CAR() from particles");
   
   // advance particle positions by half dt using CAR method in SPH geometry
   
   Real xpnew, ypnew, zpnew;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upbar = lit().velocity();
      xpnew = xpold[0] + a_cnormDt*upbar[0];
      ypnew = a_cnormDt*upbar[1];
      zpnew = a_cnormDt*upbar[2];

      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      position_virt[0] = (xpnew + xpold[0])/2.0;
      position_virt[1] = ypnew/2.0;
      position_virt[2] = zpnew/2.0;

      RealVect& rpbar = lit().position();
      const Real rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew + zpnew*zpnew);
      rpbar[0] = (rpnew + xpold[0])/2.0;

   }

}
#endif

void PicSpecies::stepNormTransfer( List<JustinsParticle>&  a_in_pList,
                                   List<JustinsParticle>&  a_out_pList,
                             const Real                    a_cnormDt, 
                             const bool                    a_reverse )
{
   CH_TIME("PicSpecies::stepNormTransfer()");

   // Check for convergence of xpbar-xpold and vpbar(xpbar)*cnormDt/2.
   // Note that vpbar has to be updated after xpbar before coming here.
   // a_reverse = true: if converged, move particle to out_pList
   // a_reverse = false: if not converged, move particle to out_pList
      
   if(m_push_type==CYL_CAR) { 
      stepNormTransfer_CYL_CAR( a_in_pList, a_out_pList, a_cnormDt, a_reverse );
   }
#if CH_SPACEDIM==1
   else if(m_push_type==SPH_CAR) {
      stepNormTransfer_SPH_CAR( a_in_pList, a_out_pList, a_cnormDt, a_reverse );
   }
#endif
   else {

      Real dxp0, rel_diff_dir, rel_diff_max;
      RealVect dxp;
#ifdef RELATIVISTIC_PARTICLES
      Real gammap, gbsq_new, gbsq_old;
      std::array<Real,3> upnew;
#endif

      const Real cnormHalfDt = 0.5*a_cnormDt;
      const RealVect& dX = m_mesh.getdX();
      
      ListIterator<JustinsParticle> li(a_in_pList);
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

         // compute step error in the particle position
         rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xpbar[dir]-xpold[dir];
#ifdef RELATIVISTIC_PARTICLES
            dxp[dir]  = upbar[dir]/gammap*cnormHalfDt;
#else
            dxp[dir]  = upbar[dir]*cnormHalfDt;
#endif
            rel_diff_dir = std::abs(dxp0-dxp[dir])/dX[dir];
            rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
         }

         // check for convergence (or lack thereof)
         if(a_reverse) { // if converged, transfer to out_pList
            if(rel_diff_max < m_rtol) { a_out_pList.transfer(li); }
            else { // update particle position and iterator
               for(int dir=0; dir<SpaceDim; dir++) { xpbar[dir] = xpold[dir] + dxp[dir]; }
               ++li;
            }
         }
         else { // update position and, if not converged, transfer to out_pList
            for(int dir=0; dir<SpaceDim; dir++) { xpbar[dir] = xpold[dir] + dxp[dir]; }
            if(rel_diff_max >= m_rtol) { 
               a_out_pList.transfer(li);
            }
            else { ++li; }
         }

      }

   }

}

void PicSpecies::stepNormTransfer_CYL_CAR( List<JustinsParticle>&  a_in_pList,
                                           List<JustinsParticle>&  a_out_pList,
                                     const Real                    a_cnormDt,
                                     const bool                    a_reverse )
{
   CH_TIME("PicSpecies::stepNormTransfer_CYL_CAR()");

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
   const Real cnormHalfDt = 0.5*a_cnormDt;
#endif

   const RealVect& dX = m_mesh.getdX();
      
   ListIterator<JustinsParticle> lit(a_in_pList);
   for(lit.begin(); lit.ok();) {
         
      RealVect& xpbar = lit().position();
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upbar = lit().velocity();
      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();

      // compute rpnew using updated value of CAR positions
      xpnew = xpold[0] + a_cnormDt*upbar[0];
      ypnew = a_cnormDt*upbar[dirp[1]];
      rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew);

      dxp0 = xpbar[0]-xpold[0];
      dxp[0] = (rpnew - xpold[0])/2.0;
      rel_diff_max = std::abs(dxp0-dxp[0])/dX[0];
#if CH_SPACEDIM==2
      dxp0 = xpbar[1]-xpold[1];
      dxp[1]  = upbar[dirp[2]]*cnormHalfDt;
      Real rel_diff_dir = std::abs(dxp0-dxp[1])/dX[1];
      rel_diff_max = std::max(rel_diff_max,rel_diff_dir);
#endif 

      // check for convergence (or lack thereof)
      if(a_reverse) { // if converged, transfer to out_pList
         if(rel_diff_max < m_rtol) { a_out_pList.transfer(lit); }
         else {
            // update particle position and iterator 
            for(int dir=0; dir<SpaceDim; dir++) { xpbar[dir] = xpold[dir] + dxp[dir]; }
            position_virt[0] = (xpnew + xpold[0])/2.0;
            position_virt[1] = ypnew/2.0;
            ++lit;
         }
      }
      else { // update position and, if not converged, transfer to out_pList
         for(int dir=0; dir<SpaceDim; dir++) { xpbar[dir] = xpold[dir] + dxp[dir]; }
         position_virt[0] = (xpnew + xpold[0])/2.0;
         position_virt[1] = ypnew/2.0;
         if(rel_diff_max >= m_rtol) { 
            a_out_pList.transfer(lit);
         }
         else { ++lit; }
      }

   }

}

#if CH_SPACEDIM==1
void PicSpecies::stepNormTransfer_SPH_CAR( List<JustinsParticle>&  a_in_pList,
                                           List<JustinsParticle>&  a_out_pList,
                                     const Real                    a_cnormDt,
                                     const bool                    a_reverse )
{
   CH_TIME("PicSpecies::stepNormTransfer_SPH_CAR()");

   // Check for convergence of rpbar-xpold and vrpbar(rpbar)*cnormDt/2.

   Real drp0, drp1, rel_diff_max;
   Real xpnew, ypnew, zpnew, rpnew;   
   
   const RealVect& dX = m_mesh.getdX();

   ListIterator<JustinsParticle> lit(a_in_pList);
   for(lit.begin(); lit.ok();) {
 
      RealVect& rpbar = lit().position();
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upbar = lit().velocity();

      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();

      // compute rpnew using updated value of CAR positions
      xpnew = xpold[0] + a_cnormDt*upbar[0];
      ypnew = a_cnormDt*upbar[1];
      zpnew = a_cnormDt*upbar[2];
      rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew + zpnew*zpnew);
         
      // compute step error in the particle position
      drp0 = rpbar[0] - xpold[0];
      drp1 = (rpnew - xpold[0])/2.0;
      rel_diff_max = std::abs(drp0-drp1)/dX[0];

      // check for convergence (or lack thereof)
      if(a_reverse) { 
         // if converged, transfer to out_pList
         if(rel_diff_max < m_rtol) { a_out_pList.transfer(lit); }
         else { // update particle position and iterator
            rpbar[0] = (rpnew + xpold[0])/2.0;
            position_virt[0] = (xpnew + xpold[0])/2.0;
            position_virt[1] = ypnew/2.0;
            position_virt[2] = zpnew/2.0;
            ++lit;
         }
      }
      else { // update position and, if not converged, transfer to out_pList
         rpbar[0] = (rpnew + xpold[0])/2.0;
         position_virt[0] = (xpnew + xpold[0])/2.0;
         position_virt[1] = ypnew/2.0;
         position_virt[2] = zpnew/2.0;
         if(rel_diff_max >= m_rtol) { 
            a_out_pList.transfer(lit); 
         }
         else { ++lit; }
      }

   }

}
#endif

void PicSpecies::transferFastParticles()
{
   if(!m_suborbit_fast_particles) { return; }
   if(m_interpJToGrid != CC1) { return; }
   CH_TIME("PicSpecies::transferFastParticles()");

   // This function is called from preRHSOp::setMassMatrices() prior to 
   // call to species->accumulateMassMatrices(). It is used to handle
   // high-energy particles with more cell crossings than permitted 
   // by the CC1 scheme for the mass matrices. These particles are 
   // sent to the suborbit container and are not included in the MM.
      
   const RealVect& dX = m_mesh.getdX();
   const int max_crossings = m_mesh.ghosts()-SpaceDim;
   int index_old, index_new, num_crossings;
   RealVect xpnew;

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      
      // get list of particles in main container
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      // get list for suborbit container
      ListBox<JustinsParticle>& box_list_suborbit = m_data_suborbit[dit];
      List<JustinsParticle>& pList_suborbit = box_list_suborbit.listItems();
   
      // loop over main container and search for fast particles
      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok();) {
 
         bool fastParticle = false;
         const RealVect& xpbar = lit().position();
         const RealVect& xpold = lit().position_old();
         for (int dir=0; dir<SpaceDim; dir++) { 
            xpnew[dir] = 2.0*xpbar[dir] - xpold[dir];
            index_old = std::floor((xpold[dir] - 0.5*dX[dir])/dX[dir]);
            index_new = std::floor((xpnew[dir] - 0.5*dX[dir])/dX[dir]);
            num_crossings = std::abs(index_new - index_old);
            if ( num_crossings > max_crossings ) {
               fastParticle = true;
               cout << "Notice: fast particle found with num_crossings = "
                    << num_crossings << " in dir = " << dir << endl;
               const Box grown_box = grow(BL.get(dit),m_mesh.ghosts());
               const int box_imin = grown_box.smallEnd(dir);
               const int box_imax = grown_box.bigEnd(dir);
               CH_assert(index_new<box_imax-1);
               CH_assert(index_new>=box_imin);
            }
         }
         if (fastParticle) {
            lit().setNumSubOrbits( 2 );
            pList_suborbit.transfer(lit);
         }
         else { ++lit; }

      }

   }

}

void PicSpecies::advanceInflowPartToBdry( RealVect&            a_xpold,
                                          Real&                a_cnormDt_sub,
                                    const Real                 a_cnormDt,
                                    const std::array<Real,3>&  a_upold,
                                    const int                  a_bdry_dir,
                                    const int                  a_bdry_side )
{
   CH_TIME("PicSpecies::advanceInflowPartToBdry()");
         
   // free-streaming advance of inflow particle to the physical boundary

   // set the boundary position
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   Real X0;
   if (a_bdry_side==0) { X0 = Xmin[a_bdry_dir]; }
   else { X0 = Xmax[a_bdry_dir]; } 
   
   Real gammap_old = 1.0; 
#ifdef RELATIVISTIC_PARTICLES
   gammap_old += a_upold[0]*a_upold[0] + a_upold[1]*a_upold[1] + a_upold[2]*a_upold[2];
   gammap_old = std::sqrt(gammap_old);
#endif

   // compute Dt for bring particle to the boundary
   const Real cnormDt0 = (X0-a_xpold[a_bdry_dir])/(a_upold[a_bdry_dir]/gammap_old);

   // update the old position to be at the boundary
   for (int dir=0; dir<SpaceDim; dir++) {
      if (dir==a_bdry_dir) { a_xpold[dir] = X0; }
      else { a_xpold[dir] = a_xpold[dir] + a_upold[dir]/gammap_old*cnormDt0; }
   }

   // set the remaining time
   a_cnormDt_sub = a_cnormDt - cnormDt0;
   CH_assert(a_cnormDt_sub>0.0);
      
}

void PicSpecies::advancePositions_2ndHalf()
{
   CH_TIME("PicSpecies::advancePositions_2ndHalf()");
   
   if(!m_motion) return;
   if(m_push_type==CYL_CAR) return;
   if(m_push_type==SPH_CAR) return;
    
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      advancePositions_2ndHalf(pList);
   }

}

void PicSpecies::advancePositions_2ndHalf( List<JustinsParticle>&  a_pList )
{
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      const RealVect& xpold = lit().position_old();
      RealVect& xp = lit().position();
      for(int dir=0; dir<SpaceDim; dir++) {
         xp[dir] = 2.0*xp[dir] - xpold[dir];
      }
   }
}

void PicSpecies::checkForNAN( List<JustinsParticle>&  a_pList,
                        const std::string&  a_string )
{
   CH_TIME("PicSpecies::checkForNAN()");
 
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         bool is_nan = false;
         const std::array<Real,3>& vp = lit().velocity();
         const std::array<Real,3>& vpold = lit().velocity();
         const RealVect& xpold = lit().position_old();
         const RealVect& xp = lit().position();
         //const std::array<Real,2>& position_virt = lit().position_virt();
         const std::array<Real,3>& Ep = lit().electric_field();
         const std::array<Real,3>& Bp = lit().magnetic_field();
    
         for(int dir = 0; dir<SpaceDim; dir++) {
            if(xp[dir]!=xp[dir]) { is_nan = true; }
         }
         //if( position_virt[0]!=position_virt[0] ) { is_nan = true; }
         if( (vp[0]!=vp[0]) || (vp[1]!=vp[1]) || (vp[2]!=vp[2]) ) { is_nan = true; }
         if( (Ep[0]!=Ep[0]) || (Ep[1]!=Ep[1]) || (Ep[2]!=Ep[2]) ) { is_nan = true; }
         if( (Bp[0]!=Bp[0]) || (Bp[1]!=Bp[1]) || (Bp[2]!=Bp[2]) ) { is_nan = true; }
         if(is_nan) {
            cout << "pic species = " << m_species << endl;
            //cout << "position virt= " << position_virt[0] << " " << position_virt[1] << endl;
            for(int dir=0; dir<SpaceDim; dir++) {cout << "old position = " << xpold[dir] << endl;}
            for(int dir=0; dir<SpaceDim; dir++) {cout << "position = " << xp[dir] << endl;}
            cout << "velocity = " << vp[0] << " " << vp[1] << " " << vp[2] << endl;
            cout << "old velocity = " << vpold[0] << " " << vpold[1] << " " << vpold[2] << endl;
            cout << "Ep = " << Ep[0] << " " << Ep[1] << " " << Ep[2] << endl;
            cout << "Bp = " << Bp[0] << " " << Bp[1] << " " << Bp[2] << endl;
            cout << a_string << endl;
            exit(EXIT_FAILURE);
         }
      }
   
}

void PicSpecies::checkForNAN( const std::string&  a_string )
{
   CH_TIME("PicSpecies::checkForNAN()");
 
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      checkForNAN(pList, a_string);
   }
   
}

void PicSpecies::applyBCs( const bool  a_intermediate_advance,
                           const Real  a_time )
{
   CH_TIME("PicSpecies::applyBCs()");
 
   if(!m_motion) return;
 
   // gather outcast particles
   m_data.gatherOutcast();
   //if(!a_intermediate_advance) m_data.gatherOutcast(); // JRA, testing

   // apply BCs to outcasts that are also out of bounds
   List<JustinsParticle>& outcast_list = m_data.outcast();
   m_species_bc->apply( outcast_list, m_surfaceCharge_nodes,
                        a_intermediate_advance, a_time );

   // remap the outcasts
   m_data.remapOutcast();
   if(!m_data.isClosed()) {
      List<JustinsParticle>& outcast_list = m_data.outcast();
      cout << "m_data is not closed in PicSpecies::applyBCs() " << endl;
      cout << "procID() = " << procID() << endl;
      cout << "species = " << m_species << endl;
      cout << "outcast_list.length() = " << outcast_list.length() << endl;
      ListIterator<JustinsParticle> lit(outcast_list);
      for(lit.begin(); lit.ok(); ++lit) {
         cout << "position = " << lit().position() << endl;
         cout << "velocity = " << lit().velocity()[0] << endl;
         cout << "position old = " << lit().position_old() << endl;
         cout << "velocity old = " << lit().velocity_old()[0] << endl;
      }
      exit(EXIT_FAILURE);
   } 
   
}

void PicSpecies::advanceVelocities( const Real  a_full_dt, 
                                    const bool  a_half_step )
{
   if(!m_forces) return;
   CH_TIME("PicSpecies::advanceVelocities()");
   
   const Real cnormDt = a_full_dt*m_cvac_norm;
 
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      applyForces(pList, cnormDt, a_half_step);

   }

}

void PicSpecies::advanceVelocities_2ndHalf( const Real  a_dt )
{
   if(!m_forces) return;
   CH_TIME("PicSpecies::advanceVelocities_2ndHalf()");
   
   if(m_push_type==CYL_HYB) {

      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         advanceVelocities_CYL_HYB_2ndHalf( a_dt, pList );
      }

   }
   else if(m_push_type==CYL_CAR) {
   
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         advanceParticles_CYL_CAR_2ndHalf( pList );
      }

   }
#if CH_SPACEDIM==1
   else if(m_push_type==SPH_HYB) {

      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         advanceVelocities_SPH_HYB_2ndHalf( a_dt, pList );
      }
   }
   else if(m_push_type==SPH_CAR) {
   
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         advanceParticles_SPH_CAR_2ndHalf( pList );
      }

   }
#endif
   else {

      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
      
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         advanceVelocities_2ndHalf( pList );

      }
      
   }

}

void PicSpecies::advanceVelocities_2ndHalf( List<JustinsParticle>&  a_pList )
{
   CH_TIME("PicSpecies::advanceVelocities_2ndHalf() from pList");

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
         
      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
       
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
   
   }

}

void PicSpecies::advanceVelocities_CYL_HYB_2ndHalf( const Real  a_dt,
                                        List<JustinsParticle>&  a_pList )
{
   CH_TIME("PicSpecies::advanceVelocities_CYL_HYB_2ndHalf()");
   
   // convert particle velocities from CYL at t_{n+1/2} to CYL at t_{n+1}
   // 1) use thetapbar to convert time-centered up from CYL to CAR
   // 2) use upbar in CAR to compute xpnew and ypnew ==>  thetapnew
   // 3) convert CAR upbar to new
   // 4) use thetapnew to convert new up from CAR to CYL

   const Real cnormDt = a_dt*m_cvac_norm;
         
   std::array<int,3> dirp = {0,1,2};
   if(m_mesh.anticyclic()) {
      dirp[1] = 2;
      dirp[2] = 1;
   }
 
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
         
      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
         
      const RealVect& xpold = lit().position_old();
            
      const std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      const Real thetabar = position_virt[0];
      Real cos_thp = std::cos(thetabar);
      Real sin_thp = std::sin(thetabar);
         
      // convert upbar from CYL to CAR
      Real upr  = up[dirp[0]];
      Real upth = up[dirp[1]];
      up[dirp[0]] = cos_thp*upr - sin_thp*upth; 
      up[dirp[1]] = sin_thp*upr + cos_thp*upth;
            
      // compute theta at t^{n+1}
      const Real xpnew = xpold[dirp[0]] + cnormDt*up[dirp[0]];
      const Real ypnew = cnormDt*up[dirp[1]];
      const Real thetanew = std::atan2(ypnew,xpnew);
      cos_thp = std::cos(thetanew);
      sin_thp = std::sin(thetanew);
            
      // convert CAR upbar from time-centered to new
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
            
      // convert upnew from CAR to CYL
      Real upx = up[dirp[0]];
      Real upy = up[dirp[1]];
      up[dirp[0]] =  cos_thp*upx + sin_thp*upy; 
      up[dirp[1]] = -sin_thp*upx + cos_thp*upy;
            
   }

}

void PicSpecies::advanceVelocities_SPH_HYB_2ndHalf( const Real  a_dt,
                                        List<JustinsParticle>&  a_pList )
{
   CH_TIME("PicSpecies::advanceVelocities_SPH_HYB_2ndHalf()");
   
   // convert particle velocities from SPH at t_{n+1/2} to CYL at t_{n+1}
   // 1) use thpbar and phpbar to convert time-centered up from SPH to CAR
   // 2) use upbar in CAR to compute xpnew, ypnew and zpnew ==>  thpnew and phpnew
   // 3) convert CAR upbar to new
   // 4) use thpnew and phpnew to convert new up from CAR to SPH

   const Real cnormDt = a_dt*m_cvac_norm;
         
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
         
      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
         
      const RealVect& xpold = lit().position_old();
            
      const std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      const Real thpbar = position_virt[0];
      const Real phpbar = position_virt[1];
      Real costhp = std::cos(thpbar);
      Real sinthp = std::sin(thpbar);
      Real cosphp = std::cos(phpbar);
      Real sinphp = std::sin(phpbar);
      
      // convert upbar from SPH to CAR
      Real upr  = up[0];
      Real upth = up[1];
      Real upph = up[2];
      up[0] = costhp*(cosphp*upr - sinphp*upph) - sinthp*upth; 
      up[1] = sinthp*(cosphp*upr - sinphp*upph) + costhp*upth;
      up[2] = sinphp*upr + cosphp*upph;
            
      // compute theta at t^{n+1}
      const Real xpnew = xpold[0] + cnormDt*up[0];
      const Real ypnew = cnormDt*up[1];
      const Real zpnew = cnormDt*up[2];
      const Real rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew + zpnew*zpnew);
      const Real thpnew = std::atan2(ypnew,xpnew);
      const Real phpnew = std::asin(zpnew/rpnew);
      costhp = std::cos(thpnew);
      sinthp = std::sin(thpnew);
      cosphp = std::cos(phpnew);
      sinphp = std::sin(phpnew);
            
      // convert CAR upbar from time-centered to new
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
            
      // convert upnew from CAR to SPH
      Real upx = up[0];
      Real upy = up[1];
      Real upz = up[2];
      up[0] =  cosphp*(costhp*upx + sinthp*upy) + sinphp*upz; 
      up[1] = -sinthp*upx + costhp*upy;
      up[2] = -sinphp*(costhp*upx + sinthp*upy) + cosphp*upz;
            
   }

}

void PicSpecies::advanceParticles_CYL_CAR_2ndHalf( List<JustinsParticle>&  a_pList )
{
   CH_TIME("PicSpecies::advanceParticles_CYL_CAR_2ndHalf()");
   
   // convert particle positions and velocities from CAR at t_{n+1/2} to CAR at t_{n+1}
   // convert new time velocities from CAR to CYL
   // rotate position vector back to theta = 0
         
   Real xpnew, ypnew, rpnew;
   Real costhp, sinthp, upx, upy;

   std::array<int,3> dirp = {0,1,2};
   if(m_mesh.anticyclic()) {
      dirp[1] = 2;
      dirp[2] = 1;
   }
   
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      RealVect& xp = lit().position();           // time-centered {rp,zp}
      std::array<Real,3>& up = lit().velocity(); // time-centered CAR
            
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upold = lit().velocity_old();

      // compute new-time CAR position from time-centered CAR position
      const std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      xpnew = 2.0*position_virt[0] - xpold[0];
      ypnew = 2.0*position_virt[1];
      rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew);
            
      // convert CAR velocity from time-centered to new time
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
            
      // rotate up new from CAR to CYL
      costhp = xpnew/rpnew;
      sinthp = ypnew/rpnew;
      upx = up[dirp[0]];
      upy = up[dirp[1]];
      up[dirp[0]] =  costhp*upx + sinthp*upy; 
      up[dirp[1]] = -sinthp*upx + costhp*upy;
  
      xp[0] = rpnew;
#if CH_SPACEDIM==2 // update z-position for 2D
      xp[1] = 2.0*xp[1] - xpold[1];
#endif

   }
 
}

#if CH_SPACEDIM==1
void PicSpecies::advanceParticles_SPH_CAR_2ndHalf( List<JustinsParticle>&  a_pList )
{
   CH_TIME("PicSpecies::advanceParticles_SPH_CAR_2ndHalf()");

   // convert particle positions and velocities from CAR at t_{n+1/2} to CAR at t_{n+1}
   // convert new time velocities from CAR to SPH
   // rotate position vector back to theta = 0 and phi = 0

   Real xpnew, ypnew, zpnew;      
   Real rpnew, rppol, costhp, sinthp, cosphp, sinphp;
   Real upx, upy, upz;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      RealVect& xp = lit().position();           // time-centered SPH position
      std::array<Real,3>& up = lit().velocity(); // time-centered CAR velocity
            
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upold = lit().velocity_old();

      // compute new-time CAR position from time-centered CAR position
      const std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      xpnew = 2.0*position_virt[0] - xpold[0];
      ypnew = 2.0*position_virt[1];
      zpnew = 2.0*position_virt[2];
            
      // convert CAR velocity from time-centered to new time
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
            
      // rotate up new from CAR to SPH
      rpnew = std::sqrt(xpnew*xpnew + ypnew*ypnew + zpnew*zpnew);
      rppol = std::sqrt(xpnew*xpnew + ypnew*ypnew);
      costhp = xpnew/rppol;
      sinthp = ypnew/rppol;
      cosphp = rppol/rpnew;
      sinphp = zpnew/rpnew;
      upx = up[0];
      upy = up[1];
      upz = up[2];
      up[0] =  cosphp*(costhp*upx + sinthp*upy) + sinphp*upz; 
      up[1] = -sinthp*upx + costhp*upy;
      up[2] = -sinphp*(costhp*upx + sinthp*upy) + cosphp*upz; 
      
      xp[0] = rpnew;
   
   }
 
}
#endif

#if CH_SPACEDIM<3
void PicSpecies::rebaseVirtualPositions()
{
   if(!m_mesh.axisymmetric()) return;
    
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      rebaseVirtualPositions( pList );
   }

}

void PicSpecies::rebaseVirtualPositions( List<JustinsParticle>&  a_pList )
{
   CH_TIME("PicSpecies::rebaseVirtualPositions() from pList");

   const int nComps = 4-CH_SPACEDIM; 

   if(m_push_type==CYL_CAR || m_push_type==SPH_CAR) {
      
      // virtual position is xpbar, ypbar, and (if SPH_CAR) zpbar
      // rotate CAR position vector back to theta = 0 and phi = 0
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         const RealVect& xp = lit().position();
         std::array<Real,nComps>& position_virt = lit().position_virt();
         position_virt[0] = xp[0];
         for (int n=1; n<nComps; n++) position_virt[n] = 0.0;
      }

   }
   else {
      
      // virtual position is theta, and (if SPH_CAR) phi
      // rotate theta and phi back to zero
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         std::array<Real,nComps>& position_virt = lit().position_virt();
         for (int n=0; n<nComps; n++) position_virt[n] = 0.0;
      }

   }
 
}
#endif

void PicSpecies::applyInertialForces( const Real  a_dt, 
                                      const bool  a_use_avg_velocity,
                                      const bool  a_update_positions,
                                      const bool  a_half_positions )
{
   if(m_push_type!=CYL_BORIS) return;
   CH_TIME("PicSpecies::applyInertialForces()");

   // See Delzanno JCP 2013

   const Real cnormDt = a_dt*m_cvac_norm;
   
   Real rp, beta_x, beta_y, x_car, y_car, beta_r, beta_th;
   
   int th_dir = 1;
   if(m_mesh.anticyclic()) th_dir = 2;
 
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {

         std::array<Real,3>& betap = lit().velocity();
         const RealVect& xpold = lit().position_old();

         // correct vr and vth velocities to account for inertial forces
         if(a_use_avg_velocity) {
            const std::array<Real,3>& betap_old = lit().velocity_old();
            beta_x = (betap[0] + betap_old[0])/2.0;
            beta_y = (betap[th_dir] + betap_old[th_dir])/2.0;
         }
         else {
            beta_x = betap[0];
            beta_y = betap[th_dir];
         }

         x_car = xpold[0] + beta_x*cnormDt;
         y_car = beta_y*cnormDt;
         rp = sqrt(x_car*x_car + y_car*y_car);

         if(rp > 0.0) { // rotate velocity vector from x,y to r,th
            beta_r = (x_car*betap[0] + y_car*betap[th_dir])/rp;
            beta_th = (x_car*betap[th_dir] - y_car*betap[0])/rp;
            betap[0] = beta_r;
            betap[th_dir] = beta_th;
            //betap[0] = (rp - xpold[0])/cnormDt;  // JRA, NO! Not energy conservative.
         }

         if(a_update_positions) {
            RealVect& xp = lit().position();
            if(a_half_positions) xp[0] = (rp + xpold[0])/2.0;
            else xp[0] = rp;
         }      

      }

   }

}

void PicSpecies::advanceParticles( const EMFields&  a_emfields,
                                   const Real       a_dt )
{
   CH_TIME("PicSpecies::advanceParticles()");
   
   // xpbar = xpn + dt/2*upbar/(0.5*(gammap_new + gammap_old))
   // upbar = upn + dt/2*q/m*(Ep(xpbar) + upbar/gammap_bar x Bp(xpbar))
   
   if (m_iter_order_swap) {
      advancePositionsImplicit( a_dt );
   }
   interpolateFieldsToParticles( a_emfields );
   addExternalFieldsToParticles( a_emfields ); 
   advanceVelocities( a_dt, true );
   if (!m_iter_order_swap) {
      advancePositionsImplicit( a_dt );
   }

}

void PicSpecies::advanceParticlesIteratively( const EMFields&  a_emfields,
                                              const Real       a_dt )
{
   CH_TIME("PicSpecies::advanceParticlesIteratively()");
   
   // Picard method for coupled half dt advance of particle positions and velocities
   // xpbar = xpn + dt/2*upbar/(0.5*(gammap_new + gammap_old))
   // upbar = upn + dt/2*q/m*(Ep(xpbar) + upbar/gammap_bar x Bp(xpbar))

   if(m_iter_max==0 || !m_motion || !m_forces || m_charge==0) {
      advanceParticles( a_emfields, a_dt );
      return;
   }
   
   const Real cnormDt = a_dt*m_cvac_norm;

   const LevelData<EdgeDataBox>& Efield = a_emfields.getFilteredElectricField();
   const LevelData<FluxBox>& Bfield = a_emfields.getMagneticField();
   const LevelData<NodeFArrayBox>& Efield_virt = a_emfields.getVirtualElectricField();
   const LevelData<FArrayBox>& Bfield_virt = a_emfields.getVirtualMagneticField();
              
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      const EdgeDataBox& Efield_inPlane = Efield[dit];
      const FluxBox& Bfield_inPlane = Bfield[dit];
      const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
      const FArrayBox& Bfield_virtual = Bfield_virt[dit];

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      m_num_parts_its += pList.length();
    
      // initial push to update particle velocities
      interpolateFieldsToParticles( pList, 
                                    Efield_inPlane, Bfield_inPlane, 
                                    Efield_virtual, Bfield_virtual ); 
      addExternalFieldsToParticles( pList, a_emfields ); 
      applyForces(pList, cnormDt, true);
      m_num_apply_its += pList.length();
     
      // update particle positions and transfer un-converged particles to a temp list
      List<JustinsParticle> temp_pList;
      const int pListLengthStart = pList.length();
      stepNormTransfer( pList, temp_pList, cnormDt, false );

      // loop over temp list and move back to main list when converged     
      int iter(1);
      while(temp_pList.length()>0) {
 
         // update particle velocities
         interpolateFieldsToParticles( temp_pList, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual ); 
         addExternalFieldsToParticles( temp_pList, a_emfields ); 
         applyForces(temp_pList, cnormDt, true);
         m_num_apply_its += temp_pList.length();
         
         // send back to pList if converged, else update particle positions
         stepNormTransfer( temp_pList, pList, cnormDt, true );

         if(temp_pList.length()==0) break;
         
         if(iter>=m_iter_max) {
         
            int num_not_converged = temp_pList.length();
            if (m_verbose_particles) {
               cout << "JRA: Picard for particles: iter = " << iter << endl;
               cout << "     num not converged = " << num_not_converged << endl;
               cout << "     m_name = " << m_name << endl;
               cout << "     cnormDt = " << cnormDt << endl;
               Box cell_box = BL.get(dit);
               inspectParticles( temp_pList, false, cell_box,
                                 Efield_inPlane, Bfield_inPlane,
                                 Efield_virtual, Bfield_virtual );
            }
            break;

         }
         iter += 1;

      }

      if( m_use_suborbit_model ) { // move non-converged particles to sub-orbit container
         ListBox<JustinsParticle>& box_list_suborbit = m_data_suborbit[dit];
         List<JustinsParticle>& pList_suborbit = box_list_suborbit.listItems();
         ListIterator<JustinsParticle> lit(temp_pList);
         for(lit.begin(); lit.ok();) { 
            lit().setNumSubOrbits( 2 );
            pList_suborbit.transfer(lit);
         }
      }
      else { // put back particles that could not converge...     
         ListIterator<JustinsParticle> lit(temp_pList);
         for(lit.begin(); lit.ok();) { pList.transfer(lit); }
         const int pListLengthEnd = pList.length();      
         if(pListLengthEnd!=pListLengthStart) exit(EXIT_FAILURE);
      }

   } // end loop over boxes on this proc
 
}

void PicSpecies::mergeSubOrbitParticles()
{
   CH_TIME("PicSpecies::mergeSubOrbitParticles()");
   if(!m_use_suborbit_model) return;

   // this function is called at the end of the pic advance prior to
   // converting xp and vp from time-centered to new time values

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
    
      ListBox<JustinsParticle>& box_list_suborbit = m_data_suborbit[dit];
      List<JustinsParticle>& pList_suborbit = box_list_suborbit.listItems();
      
      if(pList_suborbit.length()==0) continue;

      // return the suborbit particle back into the main list
      ListIterator<JustinsParticle> lit(pList_suborbit);
      for(lit.begin(); lit.ok();) { 
         lit().setNumSubOrbits( 0 );
         pList.transfer(lit); 
      }

   }

}

void PicSpecies::removeOutflowParticles()
{
   CH_TIME("PicSpecies::removeOutflowParticles()");
   
   if(!m_motion) return;
    
   // remove the outcasts that are out of bounds
   m_species_bc->removeOutflowParticles( m_surfaceCharge_nodes );
   
   if(!m_data.isClosed()) {
      List<JustinsParticle>& outcast_list = m_data.outcast();
      cout << "procID() = " << procID() << endl;
      cout << "outcast_list.length() = " << outcast_list.length() << endl;
      ListIterator<JustinsParticle> lit(outcast_list);
      for(lit.begin(); lit.ok(); ++lit) {
         cout << "position = " << lit().position() << endl;
         cout << "velocity = " << lit().velocity()[0] << endl;
         cout << "position old = " << lit().position_old() << endl;
         cout << "velocity old = " << lit().velocity_old()[0] << endl;
      }
   } 
   CH_assert(m_data.isClosed());
   
}

void PicSpecies::createInflowParticles( const Real  a_time, 
                                        const Real  a_dt  )
{
   CH_TIME("PicSpecies::createInflowParticles()");
   if(!m_motion) return;

   const Real cnormDt = a_dt*m_cvac_norm;
   m_species_bc->createInflowParticles( a_time, cnormDt, m_data );
}

void PicSpecies::injectInflowParticles()
{
   CH_TIME("PicSpecies::injectInflowParticles()");
   if( m_suborbit_inflowJ ) { return; }
   else { m_species_bc->injectInflowParticles( m_data, m_surfaceCharge_nodes ); }
}

void PicSpecies::updateOldParticlePositions()
{
   CH_TIME("PicSpecies::updateOldParticlePositions()");
   
   if(!m_motion) return;
    
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {

         RealVect& xpold = lit().position_old();
         const RealVect& xp = lit().position();
         for(int dir=0; dir<SpaceDim; dir++) xpold[dir] = xp[dir];

      }
      
   }

}

void PicSpecies::updateOldParticleVelocities()
{
   if(!m_forces) return;
   CH_TIME("PicSpecies::updateOldParticleVelocities()");
    
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {
         const std::array<Real,3>& up = lit().velocity();
         lit().setOldVelocity(up);
      }
      
   }

}

void PicSpecies::setStableDt()
{
   CH_TIME("PicSpecies::setStableDt()");
   
   // set the stable time step based on particles crossing a grid cell
   const RealVect& dX(m_mesh.getdX());
   Real maxDtinv = 0.0;
   Real gammap, thisDtinv;

   // Each proc loops over its own boxes
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      //ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = m_data[dit].listItems();
      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {

         const std::array<Real,3>& up = lit().velocity();
         gammap = 1.0;
#ifdef RELATIVISTIC_PARTICLES
         gammap += up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
         gammap = sqrt(gammap);
#endif
         for(int dir=0; dir<SpaceDim; dir++) {
            thisDtinv = abs(up[dir]/gammap)/dX[dir];
            maxDtinv = Max(thisDtinv,maxDtinv);
         }

      }
      
   }

   // update stable time step based on particle courant
   Real local_stable_dt = 1.0/maxDtinv/m_cvac_norm;
   Real stable_dt = local_stable_dt;
#ifdef CH_MPI
   MPI_Allreduce( &local_stable_dt, &stable_dt, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD ); 
#endif
   m_stable_dt = stable_dt; 

}

void PicSpecies::binTheParticles()
{
   CH_TIME("PicSpecies::binTheParticles()");
   
   ///////////////////////////////////////////////////////////
   //
   //   fill BinFab container with pointers to particles
   //
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for (dit.begin(); dit.ok(); ++dit) { // loop over boxes
      
      //m_data_binfab_ptr[dit].reBin(); 
      List<JustinsParticle>& pList = m_data[dit].listItems();
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = m_data_binfab_ptr[dit];
      
      // clear the binfab_ptr container
      const Box gridBox = BL.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices
         const IntVect ig = gbit();
         List<JustinsParticlePtr>& cell_pList_ptr = thisBinFab_ptr(ig,0);
         cell_pList_ptr.clear();
      }

      // refill the binfab_ptr container
      ListIterator<JustinsParticle> li(pList);
      CH_XD::List<JustinsParticlePtr> pListPtr;
      for(li.begin(); li.ok(); ++li) {
         JustinsParticlePtr particlePtr(li());
         pListPtr.append(particlePtr);
      }
      thisBinFab_ptr.addItems(pListPtr);
   }
}

void PicSpecies::initialize( const CodeUnits&    a_units,
                             const Real          a_time,
                             const std::string&  a_restart_file_name )
{
   // set the normalization constant for the equation of motion
   const Real Escale = a_units.getScale(a_units.ELECTRIC_FIELD);
   const Real Xscale = a_units.getScale(a_units.LENGTH);
   const Real qom    = m_charge/m_mass*Constants::QE/Constants::ME;
   const Real cvacSq = Constants::CVAC*Constants::CVAC;
   m_fnorm_const = qom/cvacSq*Escale*Xscale;

   m_cvac_norm = a_units.CvacNorm();

   if(a_restart_file_name.empty()) { initializeFromInputFile(a_units); }
   else {
      initializeFromRestartFile( a_time, a_restart_file_name );
      
      // need to set the particle electric and magnetic fields vectors to zero,
      // as they are not saved in the checkpoint files, and are thus uninitialized
      // when loading particles from a restart file. I have ran into issues where
      // a NAN is present for neutral particles that causes issues. In the future,
      // the particles will not own Ep and Bp, and these lines can be removed.
      const std::array<Real,3> zero_vect = {0.0,0.0,0.0};
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         ListIterator<JustinsParticle> lit(pList);
         for(lit.begin(); lit.ok(); ++lit) {
            lit().setElectricField( zero_vect );
            lit().setMagneticField( zero_vect );
         }
      }
   }

   if(m_push_type==CYL_CAR) {
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         ListIterator<JustinsParticle> lit(pList);
         for(lit.begin(); lit.ok(); ++lit) {
            const Real xp0 = lit().position(0);
#if CH_SPACEDIM==1
            std::array<Real,3> position_virt = {xp0,0.0,0.0};
#elif CH_SPACEDIM==2
            std::array<Real,2> position_virt = {xp0,0.0};
#endif
            lit().setPositionVirt( position_virt );
         }
      }
   }
   
#if CH_SPACEDIM==1
   if(m_push_type==SPH_CAR) {
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         ListIterator<JustinsParticle> lit(pList);
         for(lit.begin(); lit.ok(); ++lit) {
            const Real xp0 = lit().position(0);
            std::array<Real,3> position_virt = {xp0,0.0,0.0};
            lit().setPositionVirt( position_virt );
         }
      }
   }
#endif
   
   int totalParticleCount = m_data.numParticles();
   if(!procID()) {
      cout << "Finished initializing pic species " << m_name  << endl;
      cout << "total particles  = " << totalParticleCount << endl << endl;
   }
   
   // define BCs object for this species
   int verbosity = 1; 
   m_species_bc = new PicSpeciesBC( m_name, m_mass, m_charge, m_interpRhoToGrid, 
                                    m_mesh, a_units, verbosity );
      
   if(m_interp_bc_check) { // only used for charge-conserving schemes so far
      const IntVect& is_outflow_bc_lo = m_species_bc->getIsOutflowBC_lo(); 
      const IntVect& is_outflow_bc_hi = m_species_bc->getIsOutflowBC_hi(); 
      m_meshInterp->setBCcheckLo(is_outflow_bc_lo);
      m_meshInterp->setBCcheckHi(is_outflow_bc_hi);
   }
   
   // check for perturbing initial positions
   if(a_restart_file_name.empty()) {
      ParmParse pp_spc( m_name.c_str() );
      bool perturb_positions = false;
      pp_spc.query("perturb_positions",perturb_positions);
      if(perturb_positions) perturbPositions();
  }

}

void PicSpecies::initializeFromInputFile( const CodeUnits&  a_units )
{
   if(!procID()) cout << "Initializing pic species " << m_name  << " from input file..." << endl;
   int verbosity=0;
   
   // initialize boundary probes to zero   
   for (int dir=0; dir<SpaceDim; dir++) {
      m_MassOut_lo[dir] = 0.0;
      m_MassOut_hi[dir] = 0.0;
      m_MomXOut_lo[dir] = 0.0;
      m_MomXOut_hi[dir] = 0.0;
      m_MomYOut_lo[dir] = 0.0;
      m_MomYOut_hi[dir] = 0.0;
      m_MomZOut_lo[dir] = 0.0;
      m_MomZOut_hi[dir] = 0.0;
      m_EnergyOut_lo[dir] = 0.0;
      m_EnergyOut_hi[dir] = 0.0;
      //
      m_MassIn_lo[dir] = 0.0;
      m_MassIn_hi[dir] = 0.0;
      m_MomXIn_lo[dir] = 0.0;
      m_MomXIn_hi[dir] = 0.0;
      m_MomYIn_lo[dir] = 0.0;
      m_MomYIn_hi[dir] = 0.0;
      m_MomZIn_lo[dir] = 0.0;
      m_MomZIn_hi[dir] = 0.0;
      m_EnergyIn_lo[dir] = 0.0;
      m_EnergyIn_hi[dir] = 0.0;
   }
   
   // get some mesh info
   const RealVect& dX(m_mesh.getdX());
   const LevelData<FArrayBox>& Xcc(m_mesh.getXcc());

   // loop over all ICs for this species
   std::string spcIC("IC." + m_name);
   ParmParse ppspcIC( spcIC.c_str() );

   int IC_count = 0;
   bool more_ICs = true;
   while(more_ICs) {
      
      if(!ppspcIC.contains("parts_per_cell")) break;
      std::vector<int> partsPerCellstd(SpaceDim);
      ppspcIC.getarr("parts_per_cell",partsPerCellstd,0,SpaceDim);
      IntVect partsPerCell(IntVect::Zero);
      for (int dir=0; dir<SpaceDim; ++dir) {
         partsPerCell[dir] = partsPerCellstd[dir];
         CH_assert( partsPerCell[dir]>0 );
      }
   
      // parse the spatial range for the initial profile
      RealVect sXmin = m_mesh.getXmin();
      RealVect sXmax = m_mesh.getXmax();
   
      ppspcIC.query("X_min", sXmin[0]);
      ppspcIC.query("X_max", sXmax[0]);
      if(SpaceDim==2) {
         ppspcIC.query("Z_min", sXmin[1]);
         ppspcIC.query("Z_max", sXmax[1]);
      }
      if(SpaceDim==3) {
         ppspcIC.query("Y_min", sXmin[1]);
         ppspcIC.query("Y_max", sXmax[1]);
         ppspcIC.query("Z_min", sXmin[2]);
         ppspcIC.query("Z_max", sXmax[2]);
      }
      
      bool cylindrical_profile = false;
      RealVect cyl_base = IntVect::Zero;
      if(SpaceDim==2) {
          ppspcIC.query("cylindrical_profile",cylindrical_profile);
          cyl_base[0] = sXmin[0];
          cyl_base[1] = (sXmax[1]+sXmin[1])/2.0;
      }
      
      // parse the initial profiles of moments to construct
      // initial particle positions and velocities
      GridFunctionFactory  gridFactory;
      const Real this_time = 0.0;
   
      // set density profile from ICs
      const std::string spcdenIC(spcIC + ".density");
      ParmParse ppdenIC( spcdenIC.c_str() );
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppdenIC,m_mesh,verbosity);
      gridFunction->assign( m_density, m_mesh, this_time );

      // set temperature profiles from ICs
      const DisjointBoxLayout& grids(m_mesh.getDBL());
      LevelData<FArrayBox> tempProfile;
      tempProfile.define(grids,1,m_temperature.ghostVect());

      const std::string spctemp0IC(spcIC + ".temperature_0");
      ParmParse pptemp0IC( spctemp0IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionTemp0 = gridFactory.create(pptemp0IC,m_mesh,verbosity);
      gridFunctionTemp0->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_temperature[dit].copy(tempProfile[dit],0,0,1);
      }
   
      const std::string spctemp1IC(spcIC + ".temperature_1");
      ParmParse pptemp1IC( spctemp1IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionTemp1 = gridFactory.create(pptemp1IC,m_mesh,verbosity);
      gridFunctionTemp1->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_temperature[dit].copy(tempProfile[dit],0,1,1);
      }
   
      const std::string spctemp2IC(spcIC + ".temperature_2");
      ParmParse pptemp2IC( spctemp2IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionTemp2 = gridFactory.create(pptemp2IC,m_mesh,verbosity);
      gridFunctionTemp2->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_temperature[dit].copy(tempProfile[dit],0,2,1);
      }
   
      // set mean velocity profiles from ICs
      const std::string spcvel0IC(spcIC + ".velocity_0");
      ParmParse ppvel0IC( spcvel0IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionVel0 = gridFactory.create(ppvel0IC,m_mesh,verbosity);
      gridFunctionVel0->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_velocity[dit].copy(tempProfile[dit],0,0,1);
      }
   
      const std::string spcvel1IC(spcIC + ".velocity_1");
      ParmParse ppvel1IC( spcvel1IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionVel1 = gridFactory.create(ppvel1IC,m_mesh,verbosity);
      gridFunctionVel1->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_velocity[dit].copy(tempProfile[dit],0,1,1);
      }
   
      const std::string spcvel2IC(spcIC + ".velocity_2");
      ParmParse ppvel2IC( spcvel2IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionVel2 = gridFactory.create(ppvel2IC,m_mesh,verbosity);
      gridFunctionVel2->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_velocity[dit].copy(tempProfile[dit],0,2,1);
      }
      
      int totalPartsPerCell = partsPerCell.product();
      const Real volume_scale = m_mesh.getVolumeScale();
      Real cellVolume = volume_scale*dX.product();
      const Real numDen_scale = a_units.getScale(a_units.NUMBER_DENSITY);
      Real pWeight = 0.0; 
      
      // query for casting weight to float 
      bool float_weight = false;
      ppspcIC.query("use_float_for_weights",float_weight);
      
      // query for uniform particle weight option
      bool uniform_particle_weights = false;
      Real pWeight_fixed = 0.0; 
      ppspcIC.query("uniform_particle_weights",uniform_particle_weights);
      if(uniform_particle_weights) {
         pWeight_fixed = minimumParticleWeight( m_density, sXmin, sXmax,
                                                numDen_scale, cellVolume, totalPartsPerCell );
         ppspcIC.query( "fixed_weight", pWeight_fixed );
         if(float_weight) pWeight_fixed = (float)pWeight_fixed;
         if(!procID()) {
            cout << std::setprecision(16) << std::scientific << endl;
            cout << "using fixed particle weight for species = " << m_species << endl;
            cout << "pWeight_fixed  = " << pWeight_fixed << endl;
            cout << std::setprecision(8) << std::defaultfloat << endl;
         }
      }

      // query for scaling of Nppc with jacobian option
      bool nppc_jacobian_scaling = false;
      if(!nppc_jacobian_scaling) { ppspcIC.query("nppc_jacobian_scaling",nppc_jacobian_scaling); }
      const int Nppc0 = totalPartsPerCell;
      Real Ja0 = 0.0;
      Real nppc_jacobian_scaling_power = 1.0;
      if(nppc_jacobian_scaling) { 
         Ja0 = minimumJacobian();
         ppspcIC.query("nppc_jacobian_scaling.power",nppc_jacobian_scaling_power);
         if(!procID()) {
            cout << "nppc_jacobian_scaling = true" << endl;
            cout << "Ja0   = " << Ja0 << endl;
            cout << "power = " << nppc_jacobian_scaling_power << endl;
         }
      }

   ////////////////////////////////////////////////////////////////////////////////////////

      // create sub-box and dX for particles
      Box partSubBox(IntVect::Zero, partsPerCell-IntVect::Unit);
      RealVect dXpart = dX;
      if(totalPartsPerCell>1) {
         for(int dir=0; dir<SpaceDim; dir++) dXpart[dir] /= partsPerCell[dir];
      }

      // loop over boxes and set the initial particle values (pos., vel., weight)
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      
      uint64_t ID;
      bool set_ID_from_IC_count = false;
      ppspcIC.query("set_ID_from_IC_count",set_ID_from_IC_count);
      if(set_ID_from_IC_count) ID = IC_count;
      else ID = procID()*512 + 1; // hack for testing purposes
   
      const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  

      for (dit.begin(); dit.ok(); ++dit) {

         CH_XD::List<JustinsParticle> thisList;
            
         // loop over grid indices
         const Box gridBox = BL.get(dit);
         BoxIterator gbit(gridBox);
         for(gbit.begin(); gbit.ok(); ++gbit) {

            const IntVect ig = gbit(); // grid index
               
            Real local_density = m_density[dit].get(ig,0);
            if(local_density<=0.0) continue;
            
            Real local_Jacobian = Jacobian[dit].get(ig,0);
            std::array<Real,3> local_temperature;
            std::array<Real,3> local_velocity;
            for(int dir=0; dir<3; dir++) {
               local_temperature[dir] = m_temperature[dit].get(ig,dir);
               local_velocity[dir] = m_velocity[dit].get(ig,dir);
            }

            if(uniform_particle_weights) {
               pWeight = pWeight_fixed;
               totalPartsPerCell = round(numDen_scale*local_density*local_Jacobian*cellVolume/pWeight_fixed); 
               if(totalPartsPerCell>0) {
#if CH_SPACEDIM==1
                  partsPerCell[0] = totalPartsPerCell;
#elif CH_SPACEDIM==2
                  if(dX[0]<=dX[1]) {
                     partsPerCell[0] = ceil(sqrt((Real)totalPartsPerCell*dX[0]/dX[1]));
                     partsPerCell[1] = ceil((Real)totalPartsPerCell/(Real)partsPerCell[0]);
                  }
                  else {
                     partsPerCell[1] = ceil(sqrt((Real)totalPartsPerCell*dX[1]/dX[0]));
                     partsPerCell[0] = ceil((Real)totalPartsPerCell/(Real)partsPerCell[1]);
                  }
                  CH_assert(partsPerCell.product()>=totalPartsPerCell);
#else
                  cout << "JRA: uniform_particle_weights not implemented for 3D" << endl;
                  exit(EXIT_FAILURE);
#endif
                  for(int dir=0; dir<SpaceDim; dir++) dXpart[dir] = dX[dir]/partsPerCell[dir];
                  Box partSubBox_new(IntVect::Zero, partsPerCell-IntVect::Unit);
                  partSubBox = partSubBox_new;
               }
            }
            else if(nppc_jacobian_scaling) {
               totalPartsPerCell = round(std::pow(local_Jacobian/Ja0,nppc_jacobian_scaling_power))*Nppc0;
               pWeight = numDen_scale*local_density*local_Jacobian*cellVolume/(Real)totalPartsPerCell;
               if(float_weight) pWeight = (float)pWeight;
               if(totalPartsPerCell>0) {
#if CH_SPACEDIM==1
                  partsPerCell[0] = totalPartsPerCell;
#elif CH_SPACEDIM==2
                  if(dX[0]<=dX[1]) {
                     partsPerCell[0] = ceil(sqrt((Real)totalPartsPerCell*dX[0]/dX[1]));
                     partsPerCell[1] = ceil((Real)totalPartsPerCell/(Real)partsPerCell[0]);
                  }
                  else {
                     partsPerCell[1] = ceil(sqrt((Real)totalPartsPerCell*dX[1]/dX[0]));
                     partsPerCell[0] = ceil((Real)totalPartsPerCell/(Real)partsPerCell[1]);
                  }
                  CH_assert(partsPerCell.product()>=totalPartsPerCell);
#else
                  cout << "VG: nppc_jacobian_scaling not implemented for 3D" << endl;
                  exit(EXIT_FAILURE);
#endif                  
                  for(int dir=0; dir<SpaceDim; dir++) dXpart[dir] = dX[dir]/partsPerCell[dir];
                  Box partSubBox_new(IntVect::Zero, partsPerCell-IntVect::Unit);
                  partSubBox = partSubBox_new;
               }
            }
            else {
               pWeight = numDen_scale*local_density*local_Jacobian*cellVolume/(Real)totalPartsPerCell;
               if(float_weight) pWeight = (float)pWeight;
            }
            if(totalPartsPerCell<=0) continue;        

            RealVect local_Xcc; 
            for(int dir=0; dir<SpaceDim; dir++) local_Xcc[dir] = Xcc[dit].get(ig,dir);

            // loop over subgrid corresponding to where particles are placed in each grid cell    
            Real V0 = sqrt(Constants::QE/Constants::ME); // ele thermal speed at 1eV [m/s]
            BoxIterator pbit(partSubBox);
            int pCount = 0;
            for(pbit.begin(); pbit.ok(); ++pbit) {

               pCount += 1; // needed for uniform_particle_weights = true
               if(pCount>totalPartsPerCell) break;

               // set particle position uniformly on grid
               RealVect Xpart = local_Xcc - 0.5*dX;
               bool part_outside = false;
               IntVect ipg = pbit();
              
               if(cylindrical_profile) {
                  for(int dir=0; dir<SpaceDim; dir++) {
                     Xpart[dir] += (ipg[dir] + 0.5)*dXpart[dir];
                  }
                  Real Rpart = std::pow(Xpart[0]-cyl_base[0],2) + std::pow(Xpart[1]-cyl_base[1],2);
                  Rpart = std::sqrt(Rpart);
                  if(Rpart>1.0) { part_outside = true; }
               }
               else{
                  for(int dir=0; dir<SpaceDim; dir++) {
                     Xpart[dir] += (ipg[dir] + 0.5)*dXpart[dir];
                     if(Xpart[dir]<sXmin[dir] || Xpart[dir]>sXmax[dir]) { 
                        part_outside=true;
                        break;
                     }
                  }
               }


               if(part_outside) continue;

               // initialize particle velocities by randomly sampling a maxwellian
               std::array<Real,3> BetaPart = {0,0,0};
               Real thisVT;
               for(int dir=0; dir<3; dir++) { 
                  thisVT = V0*sqrt(local_temperature[dir]/m_mass); // [m/s]
                  BetaPart[dir] = (thisVT*MathUtils::randn() + local_velocity[dir])/Constants::CVAC;
               }
            
               // create this particle
               JustinsParticle particle(pWeight, Xpart, BetaPart);
               particle.setID(ID);
               if(!set_ID_from_IC_count) ID = ID + 1;

               // append particle to the list
               thisList.append(particle);
            
            }

         }

         // add particles destructively to this ListBox. 
         // those left behind are outcasts.
         m_data[dit].addItemsDestructive(thisList);
      
         // add the outcast list to be re-distributed later. 
         m_data.outcast().catenate(thisList); // no iterator? outside loop?

      }
      m_data.remapOutcast();
      CH_assert(m_data.isClosed());

      // check for more ICs
      IC_count = IC_count + 1;
      stringstream spcIC_ss;
      spcIC_ss << "IC_" << IC_count << "." << m_name; 
      spcIC = spcIC_ss.str();

      ppspcIC = spcIC.c_str();
      if(!ppspcIC.contains("parts_per_cell")) {
         IC_count = IC_count;
         more_ICs = false;
      }
  
   } // end while loop

   //
   //   check for individual particle creations
   //
 
   std::string spc_part(m_name + ".particle.0");
   ParmParse ppspc_part( spc_part.c_str() );
   if(ppspc_part.contains("weight") && procID()==0) {

      int part_count = 0;
      List<JustinsParticle> species_pList;
      while(ppspc_part.contains("weight")) {

         Real pWeight;
         RealVect Xpart;
         std::array<Real,3> BetaPart;
  
         // parse particle information
         ppspc_part.get("weight",pWeight);
         //
         std::vector<Real> Xpart_vect(SpaceDim);
         ppspc_part.getarr("position",Xpart_vect,0,SpaceDim);
         for (int n=0; n<SpaceDim; n++) Xpart[n] = Xpart_vect.at(n);
         //
         std::vector<Real> Vpart_vect(3);
         ppspc_part.getarr("velocity",Vpart_vect,0,3);
         for (int n=0; n<3; n++) BetaPart[n] = Vpart_vect.at(n)/Constants::CVAC;
               
         // create this particle and append it to the list
         JustinsParticle particle(pWeight, Xpart, BetaPart);
         uint64_t ID = part_count;
         particle.setID(ID);
         species_pList.append(particle);
       
         // update particle count and redefine ppspc_part
         part_count = part_count + 1;
         stringstream spc_part_ss;
         spc_part_ss << m_name << ".particle." << part_count; 
         spc_part = spc_part_ss.str();
         ppspc_part = spc_part.c_str();

      } // end while loop
      m_data.outcast().catenate(species_pList);
  
   }

   m_data.remapOutcast();    // put the particles in the correct box
   m_data.outcast().clear(); // remove particles created that may be out of bounds
   CH_assert(m_data.isClosed());

}

void PicSpecies::perturbPositions()
{

   if(!procID()) cout << "perturb positions for species: " << m_name << endl;

   std::string spc_pert(m_name + ".perturb_positions");
   ParmParse pp_pert( spc_pert.c_str() );

   std::string type;
   pp_pert.query("type",type);
   CH_assert(type == "sin" || type == "cos");
   
   // get the mode amplitude and number 
   Real amplitude;
   pp_pert.get("amplitude",amplitude);
   std::vector<int> temp(SpaceDim,0);
   pp_pert.getarr( "mode", temp, 0, SpaceDim );
   IntVect mode = IntVect( temp );
   if(!procID()) cout << "amplitude = " << amplitude << endl;
   if(!procID()) cout << "mode = " << mode << endl;
   
   const RealVect Xmin = m_mesh.getXmin();
   const RealVect Xmax = m_mesh.getXmax();
   const RealVect L = Xmax-Xmin;
   const Real twoPi = Constants::TWOPI;
 
   // loop over boxes  
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
   
      // loop over particle list and perturb positions
      ListIterator<JustinsParticle> li(pList);
      if(type=="sin") {
         for(li.begin(); li.ok(); ++li) {
            RealVect& xp = li().position();
            RealVect& xpold = li().position_old();
            for (int dir=0; dir<SpaceDim; dir++) {
               Real arg = twoPi*mode[dir]*xp[dir]/L[dir];
               xp[dir] = xp[dir] + amplitude*L[dir]*sin(arg);
               xpold[dir] = xp[dir];
            }
         }
      }
      else {
         for(li.begin(); li.ok(); ++li) {
            RealVect& xp = li().position();
            RealVect& xpold = li().position_old();
            for (int dir=0; dir<SpaceDim; dir++) {
               Real arg = twoPi*mode[dir]*xp[dir]/L[dir];
               xp[dir] = xp[dir] + amplitude*L[dir]*cos(arg);
               xpold[dir] = xp[dir];
            }
         }
      }
      
   }
   const Real dummy_time = 0.0;
   applyBCs(false,dummy_time);

}

void PicSpecies::initializeFromRestartFile( const Real          a_time,
                                            const std::string&  a_restart_file_name )
{
   if(!procID()) cout << "Initializing pic species " << m_name  << " from restart file..." << endl;

#ifdef CH_USE_HDF5
   HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );
   
   //  read cummulative boundary diagnostic data
   HDF5HeaderData header;
   
   // set the species group
   std::stringstream s;
   s << "species" << m_species; 
   if(!procID()) cout << s.str() << endl;
   const std::string group_name = std::string( s.str() + "_data");
   handle.setGroup(group_name);
   
   header.readFromFile( handle );
   
   // initialize boundary probes to zero
   const ProblemDomain& domain(m_mesh.getDomain()); 
   for(int dir=0; dir<SpaceDim; dir++) {         
      if(domain.isPeriodic(dir)) continue;
      std::string dirstr = std::to_string(dir);
      m_MassOut_lo[dir] = header.m_real["massOut_lo"+dirstr];
      m_MassOut_hi[dir] = header.m_real["massOut_hi"+dirstr];
      m_MomXOut_lo[dir] = header.m_real["momXOut_lo"+dirstr];
      m_MomXOut_hi[dir] = header.m_real["momXOut_hi"+dirstr];
      m_MomYOut_lo[dir] = header.m_real["momYOut_lo"+dirstr];
      m_MomYOut_hi[dir] = header.m_real["momYOut_hi"+dirstr];
      m_MomZOut_lo[dir] = header.m_real["momZOut_lo"+dirstr];
      m_MomZOut_hi[dir] = header.m_real["momZOut_hi"+dirstr];
      m_EnergyOut_lo[dir] = header.m_real["energyOut_lo"+dirstr];
      m_EnergyOut_hi[dir] = header.m_real["energyOut_hi"+dirstr];
   
      m_MassIn_lo[dir] = header.m_real["massIn_lo"+dirstr];
      m_MassIn_hi[dir] = header.m_real["massIn_hi"+dirstr];
      m_MomXIn_lo[dir] = header.m_real["momXIn_lo"+dirstr];
      m_MomXIn_hi[dir] = header.m_real["momXIn_hi"+dirstr];
      m_MomYIn_lo[dir] = header.m_real["momYIn_lo"+dirstr];
      m_MomYIn_hi[dir] = header.m_real["momYIn_hi"+dirstr];
      m_MomZIn_lo[dir] = header.m_real["momZIn_lo"+dirstr];
      m_MomZIn_hi[dir] = header.m_real["momZIn_hi"+dirstr];
      m_EnergyIn_lo[dir] = header.m_real["energyIn_lo"+dirstr];
      m_EnergyIn_hi[dir] = header.m_real["energyIn_hi"+dirstr];
   }
      
   if( m_charge != 0 ) { // read in the cumulative surface charge diagnostic
      int read_error;
      LevelData<NodeFArrayBox> sigma_temp;
      const DisjointBoxLayout& grids(m_mesh.getDBL());
      read_error = read(handle, sigma_temp, "surface_charge", grids);
      if(read_error==1) {
         if(!procID()) {
            cout << "Notice: surface_charge not found in checkpoint file "
                 << a_restart_file_name << " for group = " << group_name << endl;
            cout << "...most likely an old restart file..." << endl;
         }
      }
      else{
        sigma_temp.copyTo(m_surfaceCharge_nodes);  
#if CH_SPACEDIM>1
        postWriteSurfaceCharge();
#endif
      }
   }

   // read in the particle data
   readParticlesFromHDF( handle, m_data, "particles" );
   handle.close();

   m_data.remapOutcast();
   CH_assert(m_data.isClosed());

   zeroFields();

#else
      MayDay::Error("restart only defined with hdf5");
#endif

}

void PicSpecies::setNppc() const
{
    
   const DisjointBoxLayout& grids = m_Nppc.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_Nppc = m_Nppc[dit];
      box_Nppc.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = number;
      m_meshInterp->moment( box_Nppc,
                            box_list.listItems(),
                            m_mass,
                            thisMoment );
   }
   m_Nppc.exchange();
     
}

void PicSpecies::setNumberDensity() const
{
   CH_TIME("PicSpecies::setNumberDensity()");
    
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_density.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density[dit];
      box_rho.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = density;
      m_meshInterp->moment( box_rho,
                            box_list.listItems(),
                            m_mass,
                            thisMoment );
 
      const FArrayBox& box_Ja = Jacobian[dit];
      box_rho.divide(box_Ja,0,0,1);
      box_rho.divide(volume_scale);

   }
   m_density.exchange();
     
}

void PicSpecies::setMomentumDensity() const
{
   CH_TIME("PicSpecies::setMomentumDensity()");
    
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_momentum.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_mom = m_momentum[dit];
      box_mom.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = momentum;
      m_meshInterp->moment( box_mom,
                            box_list.listItems(),
                            m_mass/volume_scale,
                            thisMoment );
      
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_mom.nComp(); ++n) { 
         box_mom.divide(box_Ja,0,n,1);
      }
 
   }
   m_momentum.exchange(); 
     
}

void PicSpecies::setEnergyDensity() const
{
   CH_TIME("PicSpecies::setEnergyDensity()");
    
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_energy.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& box_ene = m_energy[dit];
      box_ene.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = energy;
      m_meshInterp->moment( box_ene,
                            box_list.listItems(),
                            m_mass/volume_scale,
                            thisMoment ); 
      
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_ene.nComp(); ++n) { 
         box_ene.divide(box_Ja,0,n,1);
      }

   }
   m_energy.exchange(); 
     
}

void PicSpecies::setEnergyOffDiag() const
{
   CH_TIME("PicSpecies::setEnergyOffDiag()");
    
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_energyOffDiag.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& box_eneOF = m_energyOffDiag[dit];
      box_eneOF.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = energyOffDiag;
      m_meshInterp->moment( box_eneOF,
                            box_list.listItems(),
                            m_mass/volume_scale,
                            thisMoment ); 
      
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_eneOF.nComp(); ++n) { 
         box_eneOF.divide(box_Ja,0,n,1);
      }

   }
   m_energyOffDiag.exchange(); 
     
}

void PicSpecies::setEnergyDensityFlux() const
{
   CH_TIME("PicSpecies::setEnergyDensityFlux()");
    
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids = m_energyFlux.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& box_eneF = m_energyFlux[dit];
      box_eneF.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = energyFlux;
      m_meshInterp->moment( box_eneF,
                            box_list.listItems(),
                            m_mass/volume_scale,
                            thisMoment ); 
      
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_eneF.nComp(); ++n) { 
         box_eneF.divide(box_Ja,0,n,1);
      }

   }
   m_energyFlux.exchange(); 
     
}

void PicSpecies::setNumberDensityFromBinFab() const
{
   CH_TIME("PicSpecies::setNumberDensityFromBinFab()");
   
   JustinsParticle* this_part_ptr = NULL;  
    
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   const Real dV_mapped = m_mesh.getMappedCellVolume();
   const Real volume_scale = m_mesh.getVolumeScale();
   const Real dV_phys = dV_mapped*volume_scale;
   const Real kernal = 1.0/dV_phys;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density[dit];
      box_rho.setVal(0.0);
      
      const BinFab<JustinsParticlePtr>& thisBinFab = m_data_binfab_ptr[dit];
      
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices

         const IntVect ig = gbit(); // grid index
         const List<JustinsParticlePtr>& cell_pList = thisBinFab(ig,0);
         
         Real rho0 = 0.0;
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { // loop over particles in this grid cell
            JustinsParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            const Real& wp = this_part_ptr->weight();
            rho0 += wp;
         }
         rho0 *= kernal;
         box_rho.set(ig,0,rho0);

      }
 
      const FArrayBox& box_Ja = Jacobian[dit];
      box_rho.divide(box_Ja,0,0,1);

   }

   this_part_ptr = NULL;
   delete this_part_ptr;
     
}

void PicSpecies::setMomentumDensityFromBinFab() const
{
   CH_TIME("PicSpecies::setMomentumDensityFromBinFab()");
   
   JustinsParticle* this_part_ptr = NULL;  
    
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   const Real dV_mapped = m_mesh.getMappedCellVolume();
   const Real volume_scale = m_mesh.getVolumeScale();
   const Real dV_phys = dV_mapped*volume_scale;
   const Real kernal = m_mass/dV_phys;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_mom = m_momentum[dit];
      box_mom.setVal(0.0);

      const BinFab<JustinsParticlePtr>& thisBinFab = m_data_binfab_ptr[dit];
      
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices

         const IntVect ig = gbit(); // grid index
         const List<JustinsParticlePtr>& cell_pList = thisBinFab(ig,0);
         
         Real mom0 = 0.0;
         Real mom1 = 0.0;
         Real mom2 = 0.0;
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { // loop over particles in this grid cell
            JustinsParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            const Real& wp = this_part_ptr->weight();
            const std::array<Real,3>& up = this_part_ptr->velocity();
            mom0 += wp*up[0];
            mom1 += wp*up[1];
            mom2 += wp*up[2];
         }
         mom0 *= kernal;
         mom1 *= kernal;
         mom2 *= kernal;
         box_mom.set(ig,0,mom0);
         box_mom.set(ig,1,mom1);
         box_mom.set(ig,2,mom2);

      }
 
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<3; ++n) box_mom.divide(box_Ja,0,n,1);

   }

   this_part_ptr = NULL;
   delete this_part_ptr;
     
}

void PicSpecies::setEnergyDensityFromBinFab() const
{
   CH_TIME("PicSpecies::setEnergyDensityFromBinFab()");
   
   JustinsParticle* this_part_ptr = NULL;  
    
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   const Real dV_mapped = m_mesh.getMappedCellVolume();
   const Real volume_scale = m_mesh.getVolumeScale();
   const Real dV_phys = dV_mapped*volume_scale;
   const Real kernal = 0.5*m_mass/dV_phys;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_ene = m_energy[dit];
      box_ene.setVal(0.0);

      const BinFab<JustinsParticlePtr>& thisBinFab = m_data_binfab_ptr[dit];
      
      const Box gridBox = grids.get(dit);
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices

         const IntVect ig = gbit(); // grid index
         const List<JustinsParticlePtr>& cell_pList = thisBinFab(ig,0);
         
         Real ene0 = 0.0;
         Real ene1 = 0.0;
         Real ene2 = 0.0;
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) { // loop over particles in this grid cell
            JustinsParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            const Real& wp = this_part_ptr->weight();
            const std::array<Real,3>& up = this_part_ptr->velocity();
            ene0 += wp*up[0]*up[0];
            ene1 += wp*up[1]*up[1];
            ene2 += wp*up[2]*up[2];
         }
         ene0 *= kernal;
         ene1 *= kernal;
         ene2 *= kernal;
         box_ene.set(ig,0,ene0);
         box_ene.set(ig,1,ene1);
         box_ene.set(ig,2,ene2);

      }
 
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<3; ++n) box_ene.divide(box_Ja,0,n,1);

   }

   this_part_ptr = NULL;
   delete this_part_ptr;
     
}

void PicSpecies::setChargeDensity()
{
   CH_TIME("PicSpecies::setChargeDensity()");
    
   CH_assert(m_charge != 0);
      
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   DataIterator dit = grids.dataIterator();

   for(dit.begin(); dit.ok(); ++dit) {
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      // deposit on cells 
      FArrayBox& this_rho = m_chargeDensity[dit];
      this_rho.setVal(0.0);
      
      m_meshInterp->deposit( box_list.listItems(), 
                             this_rho,
                             m_interpRhoToGrid );
 
      this_rho.mult(m_charge/volume_scale); 
      const FArrayBox& box_Ja = Jacobian[dit];
      this_rho.divide(box_Ja,0,0,1);
    
   }
   
   // apply BC to Rho (only used for symmetry BC right now)
   m_species_bc->applyToRho(m_chargeDensity);
   
   // add ghost cells to valid cells
   LDaddOp<FArrayBox> addOp;
   m_chargeDensity.exchange(m_chargeDensity.interval(), m_mesh.reverseCopier(), addOp);
   //m_chargeDensity.exchange(); // needed if more than 1 box per proccesor. bug?

}

void PicSpecies::setChargeDensityOnFaces()
{
   CH_TIME("PicSpecies::setChargeDensityOnFaces()");
    
   CH_assert(m_charge != 0);
      
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();

   for(dit.begin(); dit.ok(); ++dit) {
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      //  deposit on faces 
      for( int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_rho = m_chargeDensity_faces[dit][dir];
         this_rho.setVal(0.0);
      
         m_meshInterp->deposit( box_list.listItems(), 
                                this_rho,
                                m_interpRhoToGrid ); 
         this_rho.mult(m_charge/volume_scale); 
      }
      
   }
   
   // apply BC to Rho (only used for symmetry BC right now)
   m_species_bc->applyToRho(m_chargeDensity_faces);
   
   // add ghost cells to valid cells
   LDaddFaceOp<FluxBox> addFaceOp;
   m_chargeDensity_faces.exchange(m_chargeDensity_faces.interval(), m_mesh.reverseCopier(), addFaceOp);
   //SpaceUtils::exchangeFluxBox(m_chargeDensity_faces); // needed if more than 1 box per proccesor. bug?
   
   // divide by Jacobian after exchange (corrected Jacobian does not have ghosts)
   const LevelData<FluxBox>& Jacobian = m_mesh.getCorrectedJfc();  
   for(dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_rho = m_chargeDensity_faces[dit][dir];
         this_rho.divide(Jacobian[dit][dir],0,0,1);
      }
   }

}

void PicSpecies::setChargeDensityOnNodes( const bool  a_use_filtering )
{
   CH_TIME("PicSpecies::setChargeDensityOnNodes()");
    
   CH_assert(m_charge != 0);
      
   const Real volume_scale = m_mesh.getVolumeScale();
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();

   for(dit.begin(); dit.ok(); ++dit) {
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      //  deposit on nodes 
      FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
      this_rho.setVal(0.0);
      
      m_meshInterp->deposit( box_list.listItems(), 
                             this_rho,
                             m_interpRhoToGrid ); 
      this_rho.mult(m_charge/volume_scale); 
      
   }
   
   // apply BC to Rho (only used for symmetry BC right now)
   m_species_bc->applyToRho(m_chargeDensity_nodes);
      
   // add charge density deposited on ghost cells to valid cells
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_chargeDensity_nodes.exchange(m_chargeDensity_nodes.interval(), 
                                  m_mesh.reverseCopier(), addNodeOp);
   SpaceUtils::exchangeNodeFArrayBox(m_chargeDensity_nodes); // needed if more than 1 box per proccesor. bug?
                                                             // also needed if only 1 proc in Z for 2D RZ..?
   
   // divide by Jacobian after exchange (corrected Jacobian does not have ghosts)
   const LevelData<NodeFArrayBox>& Jacobian = m_mesh.getCorrectedJnc();  
   for(dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
      this_rho.divide(Jacobian[dit].getFab(),0,0,1);
   }
   
   if(a_use_filtering) {
      SpaceUtils::exchangeNodeFArrayBox(m_chargeDensity_nodes);
      m_species_bc->applyToRhoInGhosts(m_chargeDensity_nodes);
      SpaceUtils::applyBinomialFilter(m_chargeDensity_nodes);
   }
   
}

void PicSpecies::setCurrentDensity( const Real  a_dt,
                                    const bool  a_from_explicit_solver )
{
   CH_TIME("PicSpecies::setCurrentDensity()");
    
   CH_assert(m_charge != 0);
   
   int axisymm_car_push = 0;
   if(m_push_type==CYL_CAR) axisymm_car_push = 1;
   else if(m_push_type==SPH_CAR) axisymm_car_push = 2;
 
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();
   
   const Real cnormDt = a_dt*m_cvac_norm;
   
   // deposit particle current
   for(dit.begin(); dit.ok(); ++dit) {
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      const List<JustinsParticle>& pList = box_list.listItems();
      
      EdgeDataBox& J_inPlane = m_currentDensity[dit];
      for (int dir=0; dir<SpaceDim; ++dir) J_inPlane[dir].setVal(0.0);
#if CH_SPACEDIM<3
      FArrayBox& J_virtual = m_currentDensity_virtual[dit].getFab();
      J_virtual.setVal(0.0);
#endif
      m_meshInterp->depositCurrent( J_inPlane[0],
#if CH_SPACEDIM>=2
                                    J_inPlane[1],
#else
                                    J_virtual,
#endif
#if CH_SPACEDIM==3
                                    J_inPlane[2],
#else
                                    J_virtual,
#endif
                                    pList,
                                    m_interpJToGrid,
                                    cnormDt,
                                    axisymm_car_push,
                                    a_from_explicit_solver );

   }

   // contribute inflow/outflow current to J if being called from explicit solver.
   // This is needed for charge conserving deposits. Note that J from inflow/outflow 
   // particles for iterative implicit solvers is handled differently.
   if(a_from_explicit_solver) {
      m_species_bc->depositInflowOutflowJ( m_currentDensity, m_currentDensity_virtual,
                                          *m_meshInterp, m_interpJToGrid, cnormDt, true );
   }
   
   // multiply by charge/volume
   const Real volume_scale = m_mesh.getVolumeScale();
   for(dit.begin(); dit.ok(); ++dit) {
#if CH_SPACEDIM<3
      m_currentDensity_virtual[dit].getFab().mult(m_charge/volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_currentDensity[dit][dir].mult(m_charge/volume_scale);
      }
   }

   // apply BCs (which are done in the EMfields class), doing addOp exchange, 
   // and dividing by the corrected Jacobian are now done by caller of this function

}

void PicSpecies::advanceInflowParticlesAndSetJ( const EMFields&  a_emfields,
                                                const Real       a_dt,
                                                const bool       a_from_emjacobian )
{
   CH_TIME("PicSpecies::advanceInflowParticlesAndSetJ()");
   if(!m_suborbit_inflowJ) { return; }
   
   SpaceUtils::zero( m_inflowJ );
   SpaceUtils::zero( m_inflowJ_virtual );
   
   Vector<List<JustinsParticle>>& inflow_list_vect = m_species_bc->getInflowListVect();
             
   const LevelData<EdgeDataBox>&   Efield = a_emfields.getFilteredElectricField();
   const LevelData<FluxBox>&       Bfield = a_emfields.getMagneticField();
   const LevelData<NodeFArrayBox>& Efield_virt = a_emfields.getVirtualElectricField();
   const LevelData<FArrayBox>&     Bfield_virt = a_emfields.getVirtualMagneticField();
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
   
      const IntVect& is_inflow_bc_lo = m_species_bc->getIsInflowBC_lo();
      const IntVect& is_inflow_bc_hi = m_species_bc->getIsInflowBC_hi();
   
      bool is_inflow_bdry = false;
      if(bdry_side==0 && is_inflow_bc_lo[bdry_dir]==1) { is_inflow_bdry = true; }
      else if(bdry_side==1 && is_inflow_bc_hi[bdry_dir]==1) { is_inflow_bdry = true; }
      if( !is_inflow_bdry ) { continue; }

      List<JustinsParticle>& pList = inflow_list_vect[b];
      if(pList.length()==0) { continue; }
   
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const EdgeDataBox& Efield_inPlane = Efield[interior_dit];
         const FluxBox&     Bfield_inPlane = Bfield[interior_dit];
         const FArrayBox&   Efield_virtual = Efield_virt[interior_dit].getFab();
         const FArrayBox&   Bfield_virtual = Bfield_virt[interior_dit];
      
         EdgeDataBox& local_J = m_inflowJ[interior_dit]; 
         FArrayBox& local_Jv  = m_inflowJ_virtual[interior_dit].getFab(); 

         advanceSubOrbitParticlesAndSetJ( pList, local_J, local_Jv, a_emfields,
                                          Efield_inPlane, Bfield_inPlane,
                                          Efield_virtual, Bfield_virtual,
                                          bdry_dir, bdry_side,
                                          a_dt, a_from_emjacobian );

         // multiply by charge/volume
         const Real volume_scale = m_mesh.getVolumeScale();
#if CH_SPACEDIM<3
         local_Jv.mult(m_charge/volume_scale);
#endif
         for (int dir=0; dir<SpaceDim; ++dir) {
           local_J[dir].mult(m_charge/volume_scale);
         }

      }

   }


}

void PicSpecies::advanceSubOrbitParticlesAndSetJ( const EMFields&  a_emfields,
                                                  const Real       a_dt,
                                                  const bool       a_from_emjacobian )
{
   CH_TIME("PicSpecies::advanceSubOrbitParticlesAndSetJ()");

   SpaceUtils::zero( m_suborbitJ );
   SpaceUtils::zero( m_suborbitJ_virtual );
   
   const int dummy_bdry_dir = -1;
   const int dummy_bdry_side = -1;

   const LevelData<EdgeDataBox>& Efield = a_emfields.getFilteredElectricField();
   const LevelData<FluxBox>&     Bfield = a_emfields.getMagneticField();
   const LevelData<NodeFArrayBox>& Efield_virt = a_emfields.getVirtualElectricField();
   const LevelData<FArrayBox>&     Bfield_virt = a_emfields.getVirtualMagneticField();
              
   const BoxLayout& BL = m_data_suborbit.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      const EdgeDataBox& Efield_inPlane = Efield[dit];
      const FluxBox&     Bfield_inPlane = Bfield[dit];
      const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
      const FArrayBox& Bfield_virtual = Bfield_virt[dit];
      
      EdgeDataBox& local_J = m_suborbitJ[dit]; 
      FArrayBox& local_Jv  = m_suborbitJ_virtual[dit].getFab(); 

      ListBox<JustinsParticle>& box_list = m_data_suborbit[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      if(pList.length()==0) { continue; }

      advanceSubOrbitParticlesAndSetJ( pList, local_J, local_Jv, a_emfields,
                                       Efield_inPlane, Bfield_inPlane,
                                       Efield_virtual, Bfield_virtual,
                                       dummy_bdry_dir, dummy_bdry_side,
                                       a_dt, a_from_emjacobian );

      // multiply by charge/volume
      const Real volume_scale = m_mesh.getVolumeScale();
#if CH_SPACEDIM<3
      m_suborbitJ_virtual[dit].getFab().mult(m_charge/volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_suborbitJ[dit][dir].mult(m_charge/volume_scale);
      }

   }

}

void PicSpecies::advanceSubOrbitParticlesAndSetJ( List<JustinsParticle>&  a_pList,
                                                  EdgeDataBox&            a_local_J,
                                                  FArrayBox&              a_local_Jv,
                                            const EMFields&     a_emfields,
                                            const EdgeDataBox&  a_Efield_inPlane,
                                            const FluxBox&      a_Bfield_inPlane,
                                            const FArrayBox&    a_Efield_virtual,
                                            const FArrayBox&    a_Bfield_virtual,
                                            const int           a_bdry_dir,
                                            const int           a_bdry_side,
                                            const Real          a_dt,
                                            const bool          a_from_emjacobian )
{
   CH_TIME("PicSpecies::advanceSubOrbitParticlesAndSetJ (from box)");

   int axisymm_car_push = 0;
   if (m_push_type==CYL_CAR) { axisymm_car_push = 1; }
   else if (m_push_type==SPH_CAR) { axisymm_car_push = 2; }
   
   int iter_min = 0;
   int iter_max = m_iter_max;
   if (m_suborbit_testing) { iter_max = std::max(10,iter_max); } // used for CYL pusher test
   if (a_from_emjacobian) { iter_max += iter_max; } // increase iters for buffer

   const Real cnormDt = a_dt*m_cvac_norm;
   bool is_inflow_list = false;
   if (a_bdry_dir>=0 && a_bdry_side>=0) { is_inflow_list = true; }

   IntVect is_outflow_bc_lo = IntVect::Zero;
   IntVect is_outflow_bc_hi = IntVect::Zero;
   if (m_interp_bc_check) { // JRA testing bug I think occurs when a the initial position
                            // of a suborbit is past physical boundary
      is_outflow_bc_lo = m_species_bc->getIsOutflowBC_lo();
      is_outflow_bc_hi = m_species_bc->getIsOutflowBC_hi();
   }

   // define local containers to store J from a single particle
   EdgeDataBox this_Jp(a_local_J.box(),a_local_J.nComp());
   FArrayBox this_Jpv(a_local_Jv.box(),a_local_Jv.nComp());
      
   List<JustinsParticle> temp_pList, finished_pList;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok();) { 
         
      const int num_part_suborbits = lit().numSubOrbits();
      int num_suborbits = num_part_suborbits;
      Real cnormDt_sub = cnormDt/num_suborbits;
   
      // zero containers to store J from a single particle
      this_Jp.setVal(0.0);
      this_Jpv.setVal(0.0);

      // save the old position and velocity
      const RealVect& xpold = lit().position_old();
      const RealVect xpold0_save = xpold;
      const std::array<Real,3>& vpold = lit().velocity_old();
      const std::array<Real,3> vpold0 = vpold;

      RealVect xpold0 = xpold0_save;

      if(is_inflow_list) { // advance particle to bdry and update cnormDt_sub
         advanceInflowPartToBdry( xpold0, cnormDt_sub, cnormDt, 
                                  vpold0, a_bdry_dir, a_bdry_side );
         cnormDt_sub /= num_suborbits;
         lit().setOldPosition(xpold0);
      }

      // set xpbar and vpbar to old values for initial guess of first suborbit
      lit().setPosition(xpold0);
      lit().setVelocity(vpold0);
#if CH_SPACEDIM<3
      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      if(m_push_type==CYL_CYL || m_push_type==CYL_HYB) {
         position_virt[0] = 0.0; // dtheta
      }
      else if(m_push_type==SPH_SPH || m_push_type==SPH_HYB) {
         position_virt[0] = 0.0; // dtheta
         position_virt[1] = 0.0; // dphi
      }
      else if(m_push_type==CYL_CAR) {
         position_virt[0] = xpold0[0]; // xpbar
         position_virt[1] = 0.0;      // ypbar
      }
#if CH_SPACEDIM==1
      else if(m_push_type==SPH_CAR) {
         position_virt[0] = xpold0[0]; // xpbar
         position_virt[1] = 0.0;      // ypbar
         position_virt[2] = 0.0;      // zpbar
      }
#endif
#endif

      // put this particle in a list all by itself
      List<JustinsParticle> single_pList;
      single_pList.transfer(lit);

      // loop over suborbits for this particle
      for (int nv = 0; nv<num_suborbits; nv++) {

         int iter(0);
         while(single_pList.length()>0) {
            interpolateFieldsToParticles( single_pList, a_Efield_inPlane, a_Bfield_inPlane, 
                                          a_Efield_virtual, a_Bfield_virtual ); 
            addExternalFieldsToParticles( single_pList, a_emfields ); 
            applyForces(single_pList, cnormDt_sub, true);
         
            if(iter>=iter_min) { stepNormTransfer( single_pList, temp_pList, cnormDt_sub, true ); }
            else { advancePositionsImplicit( single_pList, cnormDt_sub ); }

            // check for reflected inflow particle and move it to the finalized list.
            // These guys cause a lot of issues when not doing this...
            if(is_inflow_list && single_pList.length()>0) {

               ListIterator<JustinsParticle> lit_single(single_pList);
               lit_single.begin();
               RealVect& xpold = lit_single().position_old();
               std::array<Real,3>& vpbar = lit_single().velocity();
               Real xpnew0 = xpold[a_bdry_dir] + vpbar[a_bdry_dir]*cnormDt_sub;

               // check for reflected particle
               const RealVect& Xmin(m_mesh.getXmin());
               const RealVect& Xmax(m_mesh.getXmax());
               if( (a_bdry_side==0 && xpnew0<Xmin[a_bdry_dir]) ||
                   (a_bdry_side==1 && xpnew0>Xmax[a_bdry_dir]) ) {
                  for(int n=0; n<3; n++) { vpbar[n] = vpold0[n]; }
                  vpbar[a_bdry_dir] = 0.0;
                  lit_single().setPosition(xpold0_save);
                  lit_single().setOldPosition(xpold0_save);
                  lit_single().setOldVelocity(vpold0);
                  lit_single().setNumSubOrbits(1);
                  this_Jp.setVal(0.0);
                  this_Jpv.setVal(0.0);
                  finished_pList.transfer(lit_single);
               }
            }
         
            if(single_pList.length()==0) { break; }

            iter += 1;
            if( iter >= iter_max ) { // this suborbit did not converge

               ListIterator<JustinsParticle> lit_single(single_pList);
               lit_single.begin();

               if( !a_from_emjacobian ) { // increase # suborbits at start again
                  cout << "JRA: iter_max reached for suborbit nv = " << nv + 1 << " of "
                       << num_suborbits << endl;
                  num_suborbits++;
                  cout << "     increasing suborbits to " << num_suborbits << endl;
                  cout << "     old position = " << xpold0 << endl;
                  lit_single().setPosition(xpold0);
                  lit_single().setOldPosition(xpold0);
                  lit_single().setVelocity(vpold0);
                  lit_single().setOldVelocity(vpold0);
                  lit_single().setNumSubOrbits( num_suborbits );
                  if(is_inflow_list) { cnormDt_sub = cnormDt_sub*(num_suborbits-1)/num_suborbits; }
                  else { cnormDt_sub = cnormDt/num_suborbits; }
                  cnormDt_sub = cnormDt/num_suborbits;
                  this_Jp.setVal(0.0);
                  this_Jpv.setVal(0.0);
#if CH_SPACEDIM<3
                  if(m_mesh.axisymmetric()) { rebaseVirtualPositions( single_pList ); }
#endif
                  nv = -1;
                  break;
               }
               else {
                  cout << "JRA: iter_max reached for suborbit nv = " << nv + 1 << " of "
                       << num_suborbits << " during linear stage... not ideal..." << endl;
                  temp_pList.transfer(lit_single);
                  break;
               }

            }
         
         }
         if(nv==-1) { continue; }
        
         if(temp_pList.length()==0) { break; } // reflected inflow particle

         // 3) deposit particles J for this suborbit
         m_meshInterp->depositCurrent( this_Jp[0],
#if CH_SPACEDIM>=2
                                       this_Jp[1],
#else
                                       this_Jpv,
#endif
#if CH_SPACEDIM==3
                                       this_Jp[2],
#else
                                       this_Jpv,
#endif
                                       temp_pList,
                                       m_interpJToGrid,
                                       cnormDt_sub,
                                       axisymm_car_push,
                                       false );

         // 4) convert particle quantities from time-centered to new
         if(m_push_type==CYL_HYB) {
            advanceVelocities_CYL_HYB_2ndHalf( a_dt/num_suborbits, temp_pList );
            advancePositions_2ndHalf( temp_pList );
         }
         else if(m_push_type==CYL_CAR) {
            advanceParticles_CYL_CAR_2ndHalf( temp_pList );
         }
#if CH_SPACEDIM==1
         else if(m_push_type==SPH_HYB) {
            advanceVelocities_SPH_HYB_2ndHalf( a_dt/num_suborbits, temp_pList );
            advancePositions_2ndHalf( temp_pList );
         }
         else if(m_push_type==SPH_CAR) {
            advanceParticles_SPH_CAR_2ndHalf( temp_pList );
         }
#endif
         else {
            advanceVelocities_2ndHalf( temp_pList );
            advancePositions_2ndHalf( temp_pList );
         }
#if CH_SPACEDIM<3
         if(m_mesh.axisymmetric()) { rebaseVirtualPositions( temp_pList ); }
#endif

         // 5) reset xpold and vpold to suborbit xpnew and vpnew
         ListIterator<JustinsParticle> lit_temp(temp_pList);
         for(lit_temp.begin(); lit_temp.ok();) {
            RealVect& xpold = lit_temp().position_old();
            std::array<Real,3>& vpold = lit_temp().velocity_old();
            if(nv==num_suborbits-1) {
               xpold = xpold0_save;
               vpold = vpold0;
               if(is_inflow_list) { // need to convert xp and vp from new to bar
                                    // for inflow particles
                  RealVect& xp = lit_temp().position();
                  std::array<Real,3>& vp = lit_temp().velocity();
                  for (int dir=0; dir<SpaceDim; dir++) {
                     xp[dir] = (xp[dir]+xpold[dir])/2.0;
                  }
                  for (int n=0; n<3; n++) {
                     vp[n] = (vp[n] + vpold[n])/2.0;
                  }
               }
               finished_pList.transfer(lit_temp);
            }
            else {
               const RealVect& xp = lit_temp().position();
               bool past_bdry = false;
               if (m_interp_bc_check) { // need to check if xp out of bounds
                  for (int dir=0; dir<SpaceDim; dir++) {
                     if (is_outflow_bc_hi[dir]) {
                        const RealVect& Xmax(m_mesh.getXmax());
                        if (xp[dir] >= Xmax[dir]) { past_bdry = true; }
                     }
                     else if (is_outflow_bc_lo[dir]) {
                        const RealVect& Xmin(m_mesh.getXmax());
                        if (xp[dir] <= Xmin[dir]) { past_bdry = true; }
                     }
                  }
               }

               if (past_bdry) {
                  xpold = xpold0_save;
                  vpold = vpold0;
                  finished_pList.transfer(lit_temp);
               }
               else { // update old values and put back in single pList
                  const std::array<Real,3>& vp = lit_temp().velocity();
                  xpold = xp;
                  vpold = vp;
                  single_pList.transfer(lit_temp);
               }
            }
         }

      } // end loop over suborbits for this particle

      // divide particle J by number of suborbits and add it to the total
      for(int dir=0; dir<SpaceDim; dir++) {
         if(is_inflow_list) { this_Jp[dir].mult(cnormDt_sub/cnormDt); }
         else { this_Jp[dir].divide(num_suborbits); }
         a_local_J[dir].plus(this_Jp[dir]);
      }
      if(is_inflow_list) { this_Jpv.mult(cnormDt_sub/cnormDt); }
      else { this_Jpv.divide(num_suborbits); }
      a_local_Jv.plus(this_Jpv);

   } // end loop over pList particles

   // transfer the particles back to the main pList owned by m_data_suborbit
   ListIterator<JustinsParticle> lit_final(finished_pList);
   for(lit_final.begin(); lit_final.ok();) { a_pList.transfer(lit_final); }

}

void PicSpecies::accumulateMassMatrices( LevelData<EdgeDataBox>&    a_sigma_xx, 
                                         LevelData<EdgeDataBox>&    a_sigma_xy,
                                         LevelData<EdgeDataBox>&    a_sigma_xz,
#if CH_SPACEDIM==1
                                         LevelData<NodeFArrayBox>&  a_sigma_yx,
                                         LevelData<NodeFArrayBox>&  a_sigma_yy,
                                         LevelData<NodeFArrayBox>&  a_sigma_yz,
#endif
                                         LevelData<NodeFArrayBox>&  a_sigma_zx,
                                         LevelData<NodeFArrayBox>&  a_sigma_zy,
                                         LevelData<NodeFArrayBox>&  a_sigma_zz,  
                                         LevelData<EdgeDataBox>&    a_J0,
                                         LevelData<NodeFArrayBox>&  a_J0v,
                                   const EMFields&                  a_emfields,
                                   const Real                       a_dt ) const
{
   CH_TIME("PicSpecies::accumulateMassMatrices()");
   if( m_charge==0 ) return;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();
   
   //
   //
   //

   const LevelData<FluxBox>& Bfield = a_emfields.getMagneticField();
   const LevelData<FArrayBox>& Bfield_virt = a_emfields.getVirtualMagneticField();

   const Real cnormDt = a_dt*m_cvac_norm;
   const Real alphas = m_fnorm_const*cnormDt/2.0;
   const Real volume_scale = m_mesh.getVolumeScale();
   const Real qovs = m_charge/volume_scale;

   int inert_type = 0;
   if(m_push_type==CYL_CYL) inert_type = 1;
   else if(m_push_type==CYL_HYB) inert_type = 0; // zero is correct here
   else if(m_push_type==CYL_CAR) inert_type = 2;
   else if(m_push_type==SPH_SPH) inert_type = 3;
   else if(m_push_type==SPH_HYB) inert_type = 0;
   else if(m_push_type==SPH_CAR) inert_type = 4;

   for(dit.begin(); dit.ok(); ++dit) {
      
      const FluxBox& Bfield_inPlane = Bfield[dit];
      const FArrayBox& Bfield_virtual = Bfield_virt[dit];   
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      const List<JustinsParticle>& pList = box_list.listItems();
      
      m_meshInterp->depositMassMatrices( a_sigma_xx[dit][0],
                                         a_sigma_xy[dit][0],
                                         a_sigma_xz[dit][0],
#if CH_SPACEDIM==1
                                         a_sigma_yx[dit].getFab(),
                                         a_sigma_yy[dit].getFab(),
                                         a_sigma_yz[dit].getFab(),
                                         a_sigma_zx[dit].getFab(),
                                         a_sigma_zy[dit].getFab(),
                                         a_sigma_zz[dit].getFab(),
                                         a_J0[dit][0],
                                         a_J0v[dit].getFab(),
                                         a_J0v[dit].getFab(),
                                         Bfield_inPlane[0],
                                         Bfield_virtual,
                                         Bfield_virtual,
#elif CH_SPACEDIM==2
                                         a_sigma_xy[dit][1], // yx
                                         a_sigma_xx[dit][1], // yy
                                         a_sigma_xz[dit][1], // yz
                                         a_sigma_zx[dit].getFab(),
                                         a_sigma_zy[dit].getFab(),
                                         a_sigma_zz[dit].getFab(),
                                         a_J0[dit][0],
                                         a_J0[dit][1],
                                         a_J0v[dit].getFab(),
                                         Bfield_inPlane[0],
                                         Bfield_inPlane[1],
                                         Bfield_virtual,
#endif
                                         qovs,
                                         alphas,
                                         cnormDt,
                                         inert_type,
                                         m_mesh.anticyclic(),
                                         pList,
                                         m_interpJToGrid );

   }
  
}

void PicSpecies::interpolateEfieldToParticles( const EMFields&  a_emfields )
{
  if(!m_forces) return;
  if(m_charge==0.0) return;
  CH_TIME("PicSpecies::interpolateEfieldToParticles()");
  
  const int blank_B = 1;
 
  const DisjointBoxLayout& grids(m_mesh.getDBL());
   
  const LevelData<EdgeDataBox>& Efield = a_emfields.getFilteredElectricField();
  const LevelData<FluxBox>& Bfield = a_emfields.getMagneticField();
#if CH_SPACEDIM < 3
  const LevelData<NodeFArrayBox>& Efield_virt = a_emfields.getVirtualElectricField();
  const LevelData<FArrayBox>& Bfield_virt = a_emfields.getVirtualMagneticField();
#endif    

  for(DataIterator dit(grids); dit.ok(); ++dit) {

    const EdgeDataBox& Efield_inPlane = Efield[dit];
    const FluxBox& Bfield_inPlane = Bfield[dit];
#if CH_SPACEDIM < 3
    const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
    const FArrayBox& Bfield_virtual = Bfield_virt[dit];
#endif
         
    ListBox<JustinsParticle>& box_list = m_data[dit];
    List<JustinsParticle>& pList = box_list.listItems();
   
    m_meshInterp->interpolateEMfieldsToPart( pList,
                                             Efield_inPlane[0],
#if CH_SPACEDIM >= 2
                                             Efield_inPlane[1],
#else
                                             Efield_virtual,
#endif
                                             Efield_virtual,
                                             Bfield_inPlane[0],
#if CH_SPACEDIM >= 2
                                             Bfield_inPlane[1],
#else
                                             Bfield_virtual,
#endif
                                             Bfield_virtual,
                                             m_interpEToParts,
                                             blank_B );
   
  }

}

void PicSpecies::interpolateFieldsToParticles( const EMFields&  a_emfields )
{
  if(!m_forces) return;
  if(m_charge==0.0) return;
  CH_TIME("PicSpecies::interpolateFieldsToParticles()");
   
  const DisjointBoxLayout& grids(m_mesh.getDBL());
   
  const LevelData<EdgeDataBox>& Efield = a_emfields.getFilteredElectricField();
  const LevelData<FluxBox>& Bfield = a_emfields.getMagneticField();

#if CH_SPACEDIM == 3
     
  for(DataIterator dit(grids); dit.ok(); ++dit) {

     const EdgeDataBox& Efield_inPlane = Efield[dit];
     const FluxBox& Bfield_inPlane = Bfield[dit];
         
     ListBox<JustinsParticle>& box_list = m_data[dit];
     List<JustinsParticle>& pList = box_list.listItems();
   
     interpolateFieldsToParticles( pList, 
                                   Efield_inPlane, Bfield_inPlane ); 

  }

#else

  const LevelData<NodeFArrayBox>& Efield_virt = a_emfields.getVirtualElectricField();
  const LevelData<FArrayBox>& Bfield_virt = a_emfields.getVirtualMagneticField();
    
  // should really loop over BL boxes from m_data. For now it is same as grids
  for(DataIterator dit(grids); dit.ok(); ++dit) {

     const EdgeDataBox& Efield_inPlane = Efield[dit];
     const FluxBox& Bfield_inPlane = Bfield[dit];
     const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
     const FArrayBox& Bfield_virtual = Bfield_virt[dit];
         
     ListBox<JustinsParticle>& box_list = m_data[dit];
     List<JustinsParticle>& pList = box_list.listItems();
   
     interpolateFieldsToParticles( pList, 
                                   Efield_inPlane, Bfield_inPlane, 
                                   Efield_virtual, Bfield_virtual ); 

  }

#endif

}

void PicSpecies::interpolateFieldsToParticles( List<JustinsParticle>&  a_pList,
                                         const EdgeDataBox&  a_Efield_inPlane,
                                         const FluxBox&      a_Bfield_inPlane )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   CH_TIME("PicSpecies::interpolateFieldsToParticles() from particle list");
   
   m_meshInterp->interpolateEMfieldsToPart( a_pList,
                                            a_Efield_inPlane[0],
                                            a_Efield_inPlane[1],
                                            a_Efield_inPlane[2],
                                            a_Bfield_inPlane[0],
                                            a_Bfield_inPlane[1],
                                            a_Bfield_inPlane[2],
                                            m_interpEToParts );
   
}

void PicSpecies::interpolateFieldsToParticles( List<JustinsParticle>&  a_pList,
                                         const EdgeDataBox&  a_Efield_inPlane,
                                         const FluxBox&      a_Bfield_inPlane,
                                         const FArrayBox&    a_Efield_virtual,
                                         const FArrayBox&    a_Bfield_virtual )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   CH_TIME("PicSpecies::interpolateFieldsToParticles() from particle list");

//#define NEW_EM_INTERP_METHOD
#ifdef NEW_EM_INTERP_METHOD   
   m_meshInterp->interpolateEMfieldsToPart_testing( a_pList,
#else
   m_meshInterp->interpolateEMfieldsToPart( a_pList,
#endif
                                            a_Efield_inPlane[0],
#if CH_SPACEDIM > 1
                                            a_Efield_inPlane[1],
#else
                                            a_Efield_virtual,
#endif
                                            a_Efield_virtual,
                                            a_Bfield_inPlane[0],
#if CH_SPACEDIM > 1
                                            a_Bfield_inPlane[1],
#else
                                            a_Bfield_virtual,
#endif
                                            a_Bfield_virtual,
                                            m_interpEToParts );

}

void PicSpecies::interpolateFieldsToParticle( JustinsParticle&  a_particle,
                                        const EdgeDataBox&      a_Efield_inPlane,
                                        const FluxBox&          a_Bfield_inPlane,
                                        const FArrayBox&        a_Efield_virtual,
                                        const FArrayBox&        a_Bfield_virtual )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   CH_TIME("PicSpecies::interpolateFieldsToParticle() single particle");

   m_meshInterp->interpolateEMfieldsToPart( a_particle,
                                            a_Efield_inPlane[0],
#if CH_SPACEDIM > 1
                                            a_Efield_inPlane[1],
#else
                                            a_Efield_virtual,
#endif
                                            a_Efield_virtual,
                                            a_Bfield_inPlane[0],
#if CH_SPACEDIM > 1
                                            a_Bfield_inPlane[1],
#else
                                            a_Bfield_virtual,
#endif
                                            a_Bfield_virtual,
                                            m_interpEToParts );

}

void PicSpecies::addExternalFieldsToParticles( const EMFields&  a_emfields )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   if(!a_emfields.externalFields()) return;
   CH_TIME("PicSpecies::addExternalFieldsToParticles()");
   
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      addExternalFieldsToParticles( pList, a_emfields );

   }
   
}

void PicSpecies::addExternalFieldsToParticles( List<JustinsParticle>&  a_pList, 
                                         const EMFields&               a_emfields )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   if(!a_emfields.externalFields()) return;
   CH_TIME("PicSpecies::addExternalFieldsToParticles() from particle list");
   
   std::array<Real,3> extE;
   std::array<Real,3> extB;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      // get external field values at particle position 
      const RealVect& Xp = lit().position();
      extE = a_emfields.getExternalE(Xp);
      extB = a_emfields.getExternalB(Xp);

      // add the external fields to the particle field arrays
      std::array<Real,3>& Ep = lit().electric_field();
      std::array<Real,3>& Bp = lit().magnetic_field();
      for (int n=0; n<3; n++) {
         Ep[n] += extE[n];
         Bp[n] += extB[n];
      }

   }

}
  
void PicSpecies::inspectBinFab( const LevelData<BinFab<JustinsParticlePtr>>&  a_binfab_ptr)
{
   JustinsParticle* this_part_ptr = NULL;  
   const DisjointBoxLayout& grids = a_binfab_ptr.disjointBoxLayout();
   const int thisProcID = 0;   

   DataIterator dit(grids);
   for (dit.begin(); dit.ok(); ++dit) { // loop over boxes

      //CH_XD::List<JustinsParticle> thisList;
      const BinFab<JustinsParticlePtr>& thisBinFab = a_binfab_ptr[dit];
     
      const Box gridBox = grids.get(dit);
      int thisNumBox = thisBinFab.numItems( gridBox );
      if(procID()==thisProcID) {
         cout << "JRA: gridBox = " << gridBox << endl;
         cout << "JRA: thisNumBox = " << thisNumBox << endl;
      }

      int comp=0;
      int thisNumCell;
      int num=0;
      BoxIterator gbit(gridBox);
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices

         const IntVect ig = gbit(); // grid index
         const List<JustinsParticlePtr>& cell_pList = thisBinFab(ig,comp);
         thisNumCell = cell_pList.length();
         
         ListIterator<JustinsParticlePtr> lit(cell_pList);
         for (lit.begin(); lit.ok(); ++lit) // loop over particles in this grid cell
         {
            JustinsParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            //const uint64_t& this_ID = this_part_ptr->ID();
            const RealVect& this_x = this_part_ptr->position();
            if(procID()==thisProcID && num==0) {
               cout << "JRA: thisNumCell = " << thisNumCell << endl;
               //cout << "JRA: ID = " << this_ID << endl;
               cout << "JRA: position = " << this_x << endl;
            }
         }
         num=1;
      
      }

   }
   
   this_part_ptr = NULL;
   delete this_part_ptr;

}

void PicSpecies::createMeshInterp()
{

   const ProblemDomain& domain(m_mesh.getDomain()); 
   const int ghosts(m_mesh.ghosts());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshLower(m_mesh.getXmin());
   const RealVect& meshUpper(m_mesh.getXmax());
  
   if(m_meshInterp!=NULL) delete m_meshInterp;
   m_meshInterp = static_cast<MeshInterp*> (new MeshInterp( domain.domainBox(),
                                                            ghosts, meshSpacing,
                                                            meshLower, meshUpper ));

}

void PicSpecies::globalMoments(std::vector<Real>&  a_global_moments) const
{
   CH_TIME("PicSpecies::globalMoments()");

   Real mass_local = 0.0;
   Real momX_local = 0.0;
   Real momY_local = 0.0;
   Real momZ_local = 0.0;
   Real energy_local = 0.0;
   
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      const ListBox<JustinsParticle>& box_list = m_data[dit];
      const List<JustinsParticle>& pList = box_list.listItems();

      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok(); ++lit) {

         const Real& wp = lit().weight();
         mass_local += wp;         
         
         const std::array<Real,3>& up = lit().velocity();
         momX_local += wp*up[0];
         momY_local += wp*up[1];
         momZ_local += wp*up[2];

#ifdef RELATIVISTIC_PARTICLES
         Real gbsq = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
         Real gammap = sqrt(1.0 + gbsq);
         energy_local += wp*gbsq*2.0/(gammap + 1.0);
#else
         energy_local += wp*(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
#endif
      }
 
   }
   mass_local *= m_mass*Constants::ME; // [kg]
   momX_local *= m_mass*Constants::ME*Constants::CVAC; // [kg-m/s]
   momY_local *= m_mass*Constants::ME*Constants::CVAC; // [kg-m/s]
   momZ_local *= m_mass*Constants::ME*Constants::CVAC; // [kg-m/s]
   energy_local *= 0.5*m_mass*Constants::ME*Constants::CVAC*Constants::CVAC; // [Joules]

   std::vector<Real> local_moments;
   local_moments.push_back(mass_local);
   local_moments.push_back(momX_local);
   local_moments.push_back(momY_local);
   local_moments.push_back(momZ_local);
   local_moments.push_back(energy_local);

   a_global_moments.resize(local_moments.size());
#ifdef CH_MPI
   MPI_Allreduce( local_moments.data(),
                  a_global_moments.data(),
                  local_moments.size(),
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
   a_global_moments = local_moments;
#endif

}

Real PicSpecies::max_wpdt(const CodeUnits&  a_units,
                          const Real&       a_dt)
{
   Real max_density_local = 0.0;

   setNumberDensity();
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      Box grid_box( grids[dit] ); 
      Real box_max = m_density[dit].max(grid_box,0);
      max_density_local = Max(max_density_local,box_max);
   } 

   Real max_density_global = max_density_local;
#ifdef CH_MPI
   MPI_Allreduce( &max_density_local,
                  &max_density_global,
                  1,
                  MPI_CH_REAL,
                  MPI_MAX,
                  MPI_COMM_WORLD );
#endif
  
  Real max_wpdt = a_units.wpNorm()*sqrt(max_density_global/m_mass)*a_dt;
  return max_wpdt;

}

void PicSpecies::bdryMoments(std::vector<Real>&  a_bdry_moments)
{
   CH_TIME("PicSpecies::bdryMoments()");
   
   const RealVect& delta_MassOut_lo = m_species_bc->getDeltaMassOut_lo();
   const RealVect& delta_MassOut_hi = m_species_bc->getDeltaMassOut_hi();
   const RealVect& delta_MomXOut_lo = m_species_bc->getDeltaMomXOut_lo();
   const RealVect& delta_MomXOut_hi = m_species_bc->getDeltaMomXOut_hi();
   const RealVect& delta_MomYOut_lo = m_species_bc->getDeltaMomYOut_lo();
   const RealVect& delta_MomYOut_hi = m_species_bc->getDeltaMomYOut_hi();
   const RealVect& delta_MomZOut_lo = m_species_bc->getDeltaMomZOut_lo();
   const RealVect& delta_MomZOut_hi = m_species_bc->getDeltaMomZOut_hi();
   const RealVect& delta_EnergyOut_lo = m_species_bc->getDeltaEnergyOut_lo();
   const RealVect& delta_EnergyOut_hi = m_species_bc->getDeltaEnergyOut_hi();
   
   const RealVect& delta_MassIn_lo = m_species_bc->getDeltaMassIn_lo();
   const RealVect& delta_MassIn_hi = m_species_bc->getDeltaMassIn_hi();
   const RealVect& delta_MomXIn_lo = m_species_bc->getDeltaMomXIn_lo();
   const RealVect& delta_MomXIn_hi = m_species_bc->getDeltaMomXIn_hi();
   const RealVect& delta_MomYIn_lo = m_species_bc->getDeltaMomYIn_lo();
   const RealVect& delta_MomYIn_hi = m_species_bc->getDeltaMomYIn_hi();
   const RealVect& delta_MomZIn_lo = m_species_bc->getDeltaMomZIn_lo();
   const RealVect& delta_MomZIn_hi = m_species_bc->getDeltaMomZIn_hi();
   const RealVect& delta_EnergyIn_lo = m_species_bc->getDeltaEnergyIn_lo();
   const RealVect& delta_EnergyIn_hi = m_species_bc->getDeltaEnergyIn_hi();
   
   std::vector<Real> local_delta;
   const Real mass_fact = m_mass*Constants::ME;
   const Real momentum_fact = mass_fact*Constants::CVAC;
   const Real energy_fact = momentum_fact*Constants::CVAC;
   for (int dir=0; dir<SpaceDim; ++dir) {
      local_delta.push_back(mass_fact*delta_MassOut_lo[dir]); // [kg]
      local_delta.push_back(mass_fact*delta_MassOut_hi[dir]); // [kg]
      local_delta.push_back(momentum_fact*delta_MomXOut_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomXOut_hi[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomYOut_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomYOut_hi[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomZOut_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomZOut_hi[dir]); // [kg-m/s]
      local_delta.push_back(energy_fact*delta_EnergyOut_lo[dir]); // [Joules]
      local_delta.push_back(energy_fact*delta_EnergyOut_hi[dir]); // [Joules]
      //
      local_delta.push_back(mass_fact*delta_MassIn_lo[dir]); // [kg]
      local_delta.push_back(mass_fact*delta_MassIn_hi[dir]); // [kg]
      local_delta.push_back(momentum_fact*delta_MomXIn_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomXIn_hi[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomYIn_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomYIn_hi[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomZIn_lo[dir]); // [kg-m/s]
      local_delta.push_back(momentum_fact*delta_MomZIn_hi[dir]); // [kg-m/s]
      local_delta.push_back(energy_fact*delta_EnergyIn_lo[dir]); // [Joules]
      local_delta.push_back(energy_fact*delta_EnergyIn_hi[dir]); // [Joules]
   }
   
   a_bdry_moments.resize(local_delta.size());
#ifdef CH_MPI
   MPI_Allreduce( local_delta.data(),
                  a_bdry_moments.data(),
                  local_delta.size(),
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
   a_bdry_moments = local_delta;
#endif
   
   // add previous running sums to bdry moments
   for (int dir=0; dir<SpaceDim; ++dir) {
      a_bdry_moments.at(0+20*dir) += m_MassOut_lo[dir];
      a_bdry_moments.at(1+20*dir) += m_MassOut_hi[dir];
      a_bdry_moments.at(2+20*dir) += m_MomXOut_lo[dir];
      a_bdry_moments.at(3+20*dir) += m_MomXOut_hi[dir];
      a_bdry_moments.at(4+20*dir) += m_MomYOut_lo[dir];
      a_bdry_moments.at(5+20*dir) += m_MomYOut_hi[dir];
      a_bdry_moments.at(6+20*dir) += m_MomZOut_lo[dir];
      a_bdry_moments.at(7+20*dir) += m_MomZOut_hi[dir];
      a_bdry_moments.at(8+20*dir) += m_EnergyOut_lo[dir];
      a_bdry_moments.at(9+20*dir) += m_EnergyOut_hi[dir];
      //
      a_bdry_moments.at(10+20*dir) += m_MassIn_lo[dir];
      a_bdry_moments.at(11+20*dir) += m_MassIn_hi[dir];
      a_bdry_moments.at(12+20*dir) += m_MomXIn_lo[dir];
      a_bdry_moments.at(13+20*dir) += m_MomXIn_hi[dir];
      a_bdry_moments.at(14+20*dir) += m_MomYIn_lo[dir];
      a_bdry_moments.at(15+20*dir) += m_MomYIn_hi[dir];
      a_bdry_moments.at(16+20*dir) += m_MomZIn_lo[dir];
      a_bdry_moments.at(17+20*dir) += m_MomZIn_hi[dir];
      a_bdry_moments.at(18+20*dir) += m_EnergyIn_lo[dir];
      a_bdry_moments.at(19+20*dir) += m_EnergyIn_hi[dir];
   }

   // update totals and reset deltas to zero   
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_MassOut_lo[dir] = a_bdry_moments.at(0+20*dir);
      m_MassOut_hi[dir] = a_bdry_moments.at(1+20*dir);
      m_MomXOut_lo[dir] = a_bdry_moments.at(2+20*dir);
      m_MomXOut_hi[dir] = a_bdry_moments.at(3+20*dir);
      m_MomYOut_lo[dir] = a_bdry_moments.at(4+20*dir);
      m_MomYOut_hi[dir] = a_bdry_moments.at(5+20*dir);
      m_MomZOut_lo[dir] = a_bdry_moments.at(6+20*dir);
      m_MomZOut_hi[dir] = a_bdry_moments.at(7+20*dir);
      m_EnergyOut_lo[dir] = a_bdry_moments.at(8+20*dir);
      m_EnergyOut_hi[dir] = a_bdry_moments.at(9+20*dir);
      //
      m_MassIn_lo[dir] = a_bdry_moments.at(10+20*dir);
      m_MassIn_hi[dir] = a_bdry_moments.at(11+20*dir);
      m_MomXIn_lo[dir] = a_bdry_moments.at(12+20*dir);
      m_MomXIn_hi[dir] = a_bdry_moments.at(13+20*dir);
      m_MomYIn_lo[dir] = a_bdry_moments.at(14+20*dir);
      m_MomYIn_hi[dir] = a_bdry_moments.at(15+20*dir);
      m_MomZIn_lo[dir] = a_bdry_moments.at(16+20*dir);
      m_MomZIn_hi[dir] = a_bdry_moments.at(17+20*dir);
      m_EnergyIn_lo[dir] = a_bdry_moments.at(18+20*dir);
      m_EnergyIn_hi[dir] = a_bdry_moments.at(19+20*dir);
   }
   m_species_bc->zeroDeltas();

}

void PicSpecies::picardParams(std::vector<Real>&  a_solver_params)
{

   Real avg_pi = avg_picard_its();   
   Real max_pi = max_avg_picard_its();

   a_solver_params.clear();
   a_solver_params.push_back(avg_pi);
   a_solver_params.push_back(max_pi);
   
   // reset the running counts
   m_num_parts_its = 0;
   m_num_apply_its = 0;

}

Real PicSpecies::avg_picard_its()
{
   // Return the average number of picard iterations over all particles for
   // this species. This is meant to represent an average number of picard 
   // its per particle per newton iteration per time step, but its actually an average
   // taken over the number of steps since the last call to this function.

   uint64_t num_parts_local = m_num_parts_its;
   uint64_t num_apply_local = m_num_apply_its;
   uint64_t num_parts_global = m_num_parts_its;
   uint64_t num_apply_global = m_num_apply_its;
#ifdef CH_MPI
   MPI_Allreduce( &num_parts_local, &num_parts_global,
                  1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( &num_apply_local, &num_apply_global,
                  1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD );
#endif 

   Real avg_its = 0.0;
   if(num_parts_global>0) avg_its = (Real)num_apply_global/(Real)num_parts_global;
   
   return avg_its;
   
}

Real PicSpecies::max_avg_picard_its()
{

   // returns the maximum average number of picard iterations amongst the processors

   Real avg_its_local = 0.0;
   if(m_num_parts_its>0) avg_its_local = (Real)m_num_apply_its/(Real)m_num_parts_its;

   Real max_avg_its = avg_its_local;
#ifdef CH_MPI
   MPI_Allreduce( &avg_its_local, &max_avg_its,
                  1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif 

   return max_avg_its;

}

Real PicSpecies::minimumParticleWeight( const LevelData<FArrayBox>&  a_density,
                                        const RealVect&              a_sXmin,
                                        const RealVect&              a_sXmax,
                                        const Real                   a_numDen_scale,
                                        const Real                   a_cell_volume,
                                        const int                    a_partsPerCell ) const
{
   Real min_pWeight_local = DBL_MAX;

   const DisjointBoxLayout& grids(a_density.getBoxes());
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   const LevelData<FluxBox>& Xfc = m_mesh.getXfc();  

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      Box grid_box( grids[dit] );
      BoxIterator gbit(grid_box);
      for(gbit.begin(); gbit.ok(); ++gbit) {

         const IntVect ig = gbit(); // grid index
          
         // only consider cells that completely contain the species profile  
         RealVect local_Xfc_lo, local_Xfc_hi; 
         for(int dir=0; dir<SpaceDim; dir++) {
            local_Xfc_lo[dir] = Xfc[dit][dir].get(ig,dir);
            IntVect ig_hi = ig;
            ig_hi[dir] += 1;
            local_Xfc_hi[dir] = Xfc[dit][dir].get(ig_hi,dir);
         }
         bool cell_contains_profile = true;
         for (int dir=0; dir<SpaceDim; dir++) {
            if(local_Xfc_lo[dir]<a_sXmin[dir] || local_Xfc_hi[dir]>a_sXmax[dir]) {
               cell_contains_profile = false;
               break;
            }
         }
         if(!cell_contains_profile) continue;

         // compute the local mass for this cell         
         Real local_density = a_density[dit].get(ig,0);
         if(local_density<=0.0) continue;
         Real local_Jacobian = Jacobian[dit].get(ig,0);
         Real local_pWeight  = a_numDen_scale*local_density*local_Jacobian*a_cell_volume/(Real)a_partsPerCell; 
         min_pWeight_local = Min(min_pWeight_local,local_pWeight);
      }

   } 

   Real min_pWeight_global = min_pWeight_local;
#ifdef CH_MPI
   MPI_Allreduce( &min_pWeight_local,
                  &min_pWeight_global,
                  1,
                  MPI_CH_REAL,
                  MPI_MIN,
                  MPI_COMM_WORLD );
#endif
  
  return min_pWeight_global;

}

Real PicSpecies::minimumJacobian( ) const
{
   Real min_Jacobian_local = DBL_MAX;

   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      Box grid_box( grids[dit] );
      BoxIterator gbit(grid_box);
      for(gbit.begin(); gbit.ok(); ++gbit) {
         const IntVect ig = gbit(); // grid index
         Real local_Jacobian = Jacobian[dit].get(ig,0);
         min_Jacobian_local = std::min(min_Jacobian_local,local_Jacobian);
      }

   } 

   Real min_Jacobian_global = min_Jacobian_local;
#ifdef CH_MPI
   MPI_Allreduce( &min_Jacobian_local,
                  &min_Jacobian_global,
                  1,
                  MPI_CH_REAL,
                  MPI_MIN,
                  MPI_COMM_WORLD );
#endif
  
  return min_Jacobian_global;

}

void PicSpecies::inspectParticles( List<JustinsParticle>&  a_pList ) const
{
            
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      cout << std::setprecision(16) << std::scientific << endl;
      cout << "position bar = " << lit().position() << endl;
      cout << "position old = " << lit().position_old() << endl;
      cout << "position new = " << 2.0*lit().position()-lit().position_old() << endl;
      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      cout << "pos_virt0 =  " << position_virt[0] << endl;
      if (position_virt.size()>1) { cout << "pos_virt1 =  " << position_virt[1] << endl; }
      if (position_virt.size()>2) { cout << "pos_virt2 =  " << position_virt[2] << endl; }
      cout << "vpbar = " << lit().velocity()[0] << ", " << lit().velocity()[1] << ", " 
                         << lit().velocity()[2] << endl;
      cout << "vpold = " << lit().velocity_old()[0] << ", " << lit().velocity_old()[1] << ", " 
                         << lit().velocity_old()[2] << endl;
         std::array<Real,3>& Ep = lit().electric_field();
         std::array<Real,3>& Bp = lit().magnetic_field();
      cout << "Ep = " << Ep[0] << ", " << Ep[1] << ", " << Ep[2] << endl;
      cout << "Bp = " << Bp[0] << ", " << Bp[1] << ", " << Bp[2] << endl;
      cout << std::setprecision(8) << std::defaultfloat << endl;
      cout << endl;
         
   }
      
}

void PicSpecies::inspectParticles( List<JustinsParticle>&  a_pList,
                             const bool          a_print_fields,
                             const Box&          a_cell_box,
                             const EdgeDataBox&  a_Efield_inPlane,
                             const FluxBox&      a_Bfield_inPlane,
                             const FArrayBox&    a_Efield_virtual,
                             const FArrayBox&    a_Bfield_virtual ) const
{
            
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      cout << std::setprecision(16) << std::scientific << endl;
      cout << "position bar = " << lit().position() << endl;
      cout << "position old = " << lit().position_old() << endl;
      cout << "position new = " << 2.0*lit().position()-lit().position_old() << endl;
      std::array<Real,4-CH_SPACEDIM>& position_virt = lit().position_virt();
      cout << "pos_virt0 =  " << position_virt[0] << endl;
      if (position_virt.size()>1) { cout << "pos_virt1 =  " << position_virt[1] << endl; }
      if (position_virt.size()>2) { cout << "pos_virt2 =  " << position_virt[2] << endl; }
      cout << "vpbar = " << lit().velocity()[0] << ", " << lit().velocity()[1] << ", " 
                         << lit().velocity()[2] << endl;
      cout << "vpold = " << lit().velocity_old()[0] << ", " << lit().velocity_old()[1] << ", " 
                         << lit().velocity_old()[2] << endl << endl;
      std::array<Real,3>& Ep = lit().electric_field();
      std::array<Real,3>& Bp = lit().magnetic_field();
      cout << "Ep = " << Ep[0] << ", " << Ep[1] << ", " << Ep[2] << endl;
      cout << "Bp = " << Bp[0] << ", " << Bp[1] << ", " << Bp[2] << endl;
      cout << std::setprecision(8) << std::defaultfloat << endl;
   }
      
   if(a_print_fields) {

      //Box cell_box = a_Efield_inPlane[0].box();
      Box cell_box = a_cell_box;
      cell_box.grow(1);
      Box node_box = surroundingNodes(cell_box);
      cout << " fnorm = " << m_fnorm_const << endl;
      cout << " cell_box = " << cell_box << endl;
      cout << " node_box = " << node_box << endl;
      BoxIterator gbit(cell_box);
      BoxIterator gnbit(node_box);

      for(gbit.begin(); gbit.ok(); ++gbit) {
         const IntVect ig = gbit(); // grid index
         cout << a_Efield_inPlane[0].get(ig,0) << ", ";
         for (int dir=0; dir<SpaceDim; dir++) {
            cout << "E[dir=" << dir << "] = " << a_Efield_inPlane[dir].get(ig,0) << ", ";
         }
      }
      cout << endl;

      for(gnbit.begin(); gnbit.ok(); ++gnbit) {
         const IntVect ig = gnbit(); // grid index
         for (int n=0; n<a_Efield_virtual.nComp(); n++) {
            cout << "Ev(n=" << n << ") = " <<  a_Efield_virtual.get(ig,n) << ", ";
         }
      }
      cout << endl;
      
      for(gnbit.begin(); gnbit.ok(); ++gnbit) {
         const IntVect ig = gnbit(); // grid index
         for (int dir=0; dir<SpaceDim; dir++) {
            cout << "B[dir=" << dir << "] = " << a_Bfield_inPlane[dir].get(ig,0) << ", ";
         }
      }
      cout << endl;

      for(gbit.begin(); gbit.ok(); ++gbit) {
         const IntVect ig = gbit(); // grid index
         for (int n=0; n<a_Bfield_virtual.nComp(); n++) {
            cout << "Bv(n=" << n << ") = " <<  a_Bfield_virtual.get(ig,n) << ", ";
         }
      }
      cout << endl;

   }

}

bool PicSpecies::isSpecies( const string&  a_name ) const
{
   if(name() == a_name) return true;
   return false;
}


#include "NamespaceFooter.H"

