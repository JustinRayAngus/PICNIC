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

   if(m_mesh.axisymmetric()) {

      // get the push type for axisymmetric geometries
      m_push_type = CYL_CYL;
      std::string push_type = "CYL_CYL";
      a_ppspc.query("push_type",push_type);
      if(push_type=="CYL_BORIS") m_push_type = CYL_BORIS;
      if(push_type=="CYL_HYBRID") m_push_type = CYL_HYBRID;
      if(push_type=="CYL_CAR") m_push_type = CYL_CAR;

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
   if(m_charge_conserving_deposit) m_interp_bc_check = true;
   a_ppspc.query( "interp_bc_check", m_interp_bc_check );
   a_ppspc.query( "suborbit_inflow_J", m_suborbit_inflowJ );
   a_ppspc.query( "use_suborbit_model", m_use_suborbit_model );
   if(m_use_suborbit_model) a_ppspc.query("suborbit_testing",m_suborbit_testing);

   a_ppspc.get( "mass", m_mass );
   a_ppspc.get( "charge", m_charge ); 
   a_ppspc.query( "potential", m_Uint );
   a_ppspc.query( "motion", m_motion );
   a_ppspc.query( "forces", m_forces );
   a_ppspc.query( "scatter", m_scatter );
   a_ppspc.query( "write_all_comps", m_write_all_part_comps );

   if(!m_mesh.axisymmetric()) {
      if(m_charge==0) m_forces = false;
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
      cout << "  interpolate bc check = " << m_interp_bc_check << endl;
      cout << "  suborbit inflow J = " << m_suborbit_inflowJ << endl;
      cout << "  use suborbit model = " << m_use_suborbit_model << endl;
      if(m_suborbit_testing) {
	 cout << "  suborbit testing = true (fixed suborbit iter_max = 10)" << endl;
      }
      if(m_push_type==CYL_BORIS)  cout << "  push type = CYL_BORIS" << endl;
      if(m_push_type==CYL_CYL)    cout << "  push type = CYL_CYL" << endl;
      if(m_push_type==CYL_HYBRID) cout << "  push type = CYL_HYBRID" << endl;
      if(m_push_type==CYL_CAR) cout << "  push type = CYL_CAR" << endl;
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
         
         const Real costhp = lit().position_virt(0);
         const Real sinthp = lit().position_virt(1);

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
   
   if(m_push_type==CYL_CYL) {
      PicSpeciesUtils::applyForcesAxisymm( a_pList, m_fnorm_const, a_cnormDt,
                                           a_byHalfDt, m_mesh.anticyclic() );
   }
   else if(m_push_type==CYL_HYBRID) {
      PicSpeciesUtils::applyForcesHYBRID( a_pList, m_fnorm_const, a_cnormDt,
                                          m_mesh.anticyclic() );
   }
   else {
      PicSpeciesUtils::applyForces( a_pList, m_fnorm_const, a_cnormDt, 
                                    a_byHalfDt, m_mesh.anticyclic() ); 
   }
   
   if(m_push_type==CYL_HYBRID) { // compute thetabar and convert upbar from CAR to CYL
      
      std::array<int,3> dirp = {0,1,2};
      if(m_mesh.anticyclic()) {
         dirp[1] = 2;
         dirp[2] = 1;
      }

      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         
         std::array<Real,2>& position_virt = lit().position_virt();
         Real thetabar = position_virt[0];
         //if(thetabar!=0.0) break;
         if(thetabar==0.0) {
            const RealVect& xpold = lit().position_old();
            std::array<Real,3>& upbar = lit().velocity();
            
            const Real xpbar = xpold[dirp[0]] + a_cnormDt/2.0*upbar[dirp[0]];
            const Real ypbar = a_cnormDt/2.0*upbar[dirp[1]];
            thetabar = std::atan2(ypbar,xpbar);
            position_virt[0] = thetabar;
            Real costhp = std::cos(thetabar);
            Real sinthp = std::sin(thetabar);

            const Real upx = upbar[dirp[0]];
            const Real upy = upbar[dirp[1]];
            upbar[dirp[0]] =  costhp*upx + sinthp*upy; 
            upbar[dirp[1]] = -sinthp*upx + costhp*upy;
         }

      }

   }
   
   if(m_push_type==CYL_CAR) { // transform upbar from CAR to CYL
      
      std::array<int,3> dirp = {0,1,2};
      if(m_mesh.anticyclic()) {
         dirp[1] = 2;
         dirp[2] = 1;
      }
      
      ListIterator<JustinsParticle> lit(a_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         
         const RealVect& xpold = lit().position_old();
         const std::array<Real,3>& upbar = lit().velocity();
         std::array<Real,2>& position_virt = lit().position_virt();
         
         const RealVect& xpbar = lit().position();
        
         const Real xpbar0 = xpold[0] + a_cnormDt/2.0*upbar[0];
         const Real ypbar0 = a_cnormDt/2.0*upbar[dirp[1]];

         const Real costhp = xpbar0/xpbar[0];
         const Real sinthp = ypbar0/xpbar[0];
         
         // update costhp and sinthp
         position_virt[0] = costhp;
         position_virt[1] = sinthp; 
    
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

void PicSpecies::advancePositionsImplicit( const Real  a_full_dt,
                                           const bool  a_half_step )
{
   if(!m_motion) return;
   CH_TIME("PicSpecies::advancePositionsImplicit()");
   
   const Real cnormDt = m_cvac_norm*a_full_dt;
   Real dt_factor = 1.0;
   if(a_half_step) dt_factor = 0.5; 

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      if(m_push_type==CYL_CAR) advancePositionsCYLCAR( pList, cnormDt );
      else advancePositionsImplicit( pList, cnormDt*dt_factor );
   }

}

void PicSpecies::advancePositionsImplicit( List<JustinsParticle>&  a_pList,
                                     const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositionsImplicit() from particles");
   
   // xp^{n+1} = xp^{n} + vp*cnormDt;
   // vp = (up^{n+1}+up^{n}/(gamma^{n+1} + gamma^{n})
   // gamma^{n+1} = sqrt(1 + up^{n+1})
   // gamma^{n} = sqrt(1 + up^{n})

#ifdef RELATIVISTIC_PARTICLES
   Real gammap, gbsq_old, gbsq_new;
   std::array<Real,3> upnew;
#endif
   
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
      for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + upbar[dir]/gammap*a_cnormDt;
#else
      for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + upbar[dir]*a_cnormDt;
#endif
   }

}

void PicSpecies::advancePositionsCYLCAR( List<JustinsParticle>&  a_pList,
                                   const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositionsCYLCAR() from particles");
   
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      
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
 
      xpbar[0] = (rpnew + xpold[0])/2.0;
#if CH_SPACEDIM==2
      const std::array<Real,3>& upbar = lit().velocity();
      xpbar[1] = xpold[1] + upbar[1]*a_cnormDt/2.0;
#endif

   }

}

void PicSpecies::advanceInflowPositions( List<JustinsParticle>&  a_pList,
                                   const int                     a_bdry_dir,
                                   const int                     a_bdry_side,
                                   const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advanceInflowPositions() from particle list");
         
   // Note that this function is only called from implicit solver 

   // set the boundary position
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   RealVect X0;
   if(a_bdry_side==0) X0[a_bdry_dir] = Xmin[a_bdry_dir];
   else X0[a_bdry_dir] = Xmax[a_bdry_dir]; 
   
   Real cnormDt0, cnormDt1;
   Real gammap_old, gammap;
#ifdef RELATIVISTIC_PARTICLES
   Real gbsq_new;
   std::array<Real,3> upnew;
#endif

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      RealVect& xp = lit().position();
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upold = lit().velocity_old();

      gammap_old = 1.0; 
#ifdef RELATIVISTIC_PARTICLES
      gammap_old += upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
      gammap_old = sqrt(gammap_old);
#endif

      cnormDt0 = (X0[a_bdry_dir]-xpold[a_bdry_dir])/(upold[a_bdry_dir]/gammap_old);
      for(int dir=0; dir<SpaceDim; dir++) {
         if(dir==a_bdry_dir) continue;
         X0[dir] = xpold[dir] + upold[dir]/gammap_old*cnormDt0;
      }
      cnormDt1 = a_cnormDt - cnormDt0;
      CH_assert(cnormDt1>0.0);
      const std::array<Real,3>& upbar = lit().velocity();
      
      gammap = 1.0; 
#ifdef RELATIVISTIC_PARTICLES
      for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
      gbsq_new = upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
      gammap = 0.5*(gammap_old + sqrt(1.0 + gbsq_new));
#endif

      for(int dir=0; dir<SpaceDim; dir++) {
         xp[dir] = X0[dir] + upbar[dir]/gammap*cnormDt1;
         xp[dir] = (xpold[dir] + xp[dir])/2.0; // convert to xpbar
      }
      if(upbar[a_bdry_dir]==0.0) xp[a_bdry_dir] = xpold[a_bdry_dir];
   }
      
}

void PicSpecies::advancePositions_2ndHalf()
{
   CH_TIME("PicSpecies::advancePositions_2ndHalf()");
   
   if(!m_motion) return;
   if(m_push_type==CYL_CAR) return;
    
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

void PicSpecies::applyBCs( const bool  a_intermediate_advance )
{
   CH_TIME("PicSpecies::applyBCs()");
 
   if(!m_motion) return;
 
   // gather outcast particles
   m_data.gatherOutcast();
   //if(!a_intermediate_advance) m_data.gatherOutcast(); // JRA, testing

   // apply BCs to outcasts that are also out of bounds
   List<JustinsParticle>& outcast_list = m_data.outcast();
   Real dummy_time = 0.0;
   m_species_bc->apply(outcast_list,a_intermediate_advance,dummy_time);

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
   
   if(m_push_type==CYL_HYBRID) {

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

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
         
      std::array<Real,3>& up = lit().velocity();
      const std::array<Real,3>& upold = lit().velocity_old();
       
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
   
      if(m_push_type==CYL_CYL) {
         std::array<Real,2>& position_virt = lit().position_virt();
         position_virt[0] = 0.0;
      }

   }

}

void PicSpecies::advanceVelocities_CYL_HYB_2ndHalf( const Real  a_dt,
		                        List<JustinsParticle>&  a_pList )
{
   
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
            
      std::array<Real,2>& position_virt = lit().position_virt();
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
            
      // reset theta back to zero
      position_virt[0] = 0.0;

   }

}

void PicSpecies::advanceParticles_CYL_CAR_2ndHalf( List<JustinsParticle>&  a_pList )
{
   
   // convert particle positions from CYL at t_{n+1/2} to CYL at t_{n+1}
   // convert particle velocities from CAR at t_{n+1/2} to CYL at t_{n+1}
         
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

      std::array<Real,2>& position_virt = lit().position_virt();
      long double costhp = position_virt[0];
      long double sinthp = position_virt[1];
 
      // rotate time-centered up from CYL to CAR
      //const long double upr = up[dirp[0]];
      //const long double upth = up[dirp[1]];
      //const long double ratio = 1.0/(costhp*costhp + sinthp*sinthp);
      //up[dirp[0]] = ratio*(costhp*upr - sinthp*upth); 
      //up[dirp[1]] = ratio*(sinthp*upr + costhp*upth);
            
      // compute rp new
      const Real xpbar0 = costhp*xp[0]; // time-centered xp
      const Real ypbar0 = sinthp*xp[0]; // time-centered yp
      const Real xpnew = 2.0*xpbar0 - xpold[0];
      const Real ypnew = 2.0*ypbar0;
      const Real rpnew = sqrt(xpnew*xpnew + ypnew*ypnew);
            
      // convert CAR up from time-centered to new
      up[0] = 2.0*up[0] - upold[0];
      up[1] = 2.0*up[1] - upold[1];
      up[2] = 2.0*up[2] - upold[2];
            
      // rotate up new from CAR to CYL
      costhp = xpnew/rpnew;
      sinthp = ypnew/rpnew;
      const Real upx = up[dirp[0]];
      const Real upy = up[dirp[1]];
      up[dirp[0]] =  costhp*upx + sinthp*upy; 
      up[dirp[1]] = -sinthp*upx + costhp*upy;
  
      // convert time-centered {rp,zp} to new
      xp[0] = rpnew;
      //xp[0] = 2.0*xp[0] - xpold[0];
#if CH_SPACEDIM==2
      xp[1] = 2.0*xp[1] - xpold[1];
#endif
      // rebase virtual costhp and sinthp back to theta = 0
      position_virt[0] = 1.0;
      position_virt[1] = 0.0;
         
   }
 
}

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

void PicSpecies::advanceParticles( const ElectroMagneticFields&  a_em_fields,
                                   const Real                    a_dt )
{
   CH_TIME("PicSpecies::advanceParticles()");
   
   // xpbar = xpn + dt/2*upbar/(0.5*(gammap_new + gammap_old))
   // upbar = upn + dt/2*q/m*(Ep(xpbar) + upbar/gammap_bar x Bp(xpbar))
   
   if (m_iter_order_swap) {
      advancePositionsImplicit( a_dt, true );
   }
   interpolateFieldsToParticles( a_em_fields );
   addExternalFieldsToParticles( a_em_fields ); 
   advanceVelocities( a_dt, true );
   if (!m_iter_order_swap) {
      advancePositionsImplicit( a_dt, true );
   }

}

void PicSpecies::advanceParticlesIteratively( const ElectroMagneticFields&  a_em_fields,
                                              const Real                    a_dt )
{
   CH_TIME("PicSpecies::advanceParticlesIteratively()");
   
   // Picard method for coupled half dt advance of particle positions and velocities
   // xpbar = xpn + dt/2*upbar/(0.5*(gammap_new + gammap_old))
   // upbar = upn + dt/2*q/m*(Ep(xpbar) + upbar/gammap_bar x Bp(xpbar))

   if(m_iter_max==0 || !m_motion || !m_forces || m_charge==0) {
      advanceParticles( a_em_fields, a_dt );
      return;
   }
   
   const Real cnormDt = a_dt*m_cvac_norm;
   const Real cnormHalfDt = 0.5*cnormDt;

   const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
   const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();
   const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
   const LevelData<FArrayBox>& Bfield_virt = a_em_fields.getVirtualMagneticField();
              
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
    
      // update xpbar and upbar 
      if(m_iter_order_swap) {
         if(m_push_type==CYL_CAR) advancePositionsCYLCAR( pList, cnormDt );
         else advancePositionsImplicit( pList, cnormHalfDt );
      }
      interpolateFieldsToParticles( pList, 
                                    Efield_inPlane, Bfield_inPlane, 
                                    Efield_virtual, Bfield_virtual ); 
      addExternalFieldsToParticles( pList, a_em_fields ); 
      applyForces(pList, cnormDt, true);
      m_num_apply_its += pList.length();
     
      // for high precision, sometime need to do at least two force applys
      // (trying to do this with rtol gets too close to machine precision)
      if(m_iter_min_two) { 
         if(m_push_type==CYL_CAR) advancePositionsCYLCAR( pList, cnormDt );
         else advancePositionsImplicit( pList, cnormHalfDt ); 
         interpolateFieldsToParticles( pList, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual ); 
         addExternalFieldsToParticles( pList, a_em_fields ); 
         applyForces(pList, cnormDt, true);
         m_num_apply_its += pList.length();
      }
       
      // transfer not converged particles to a temp list
      List<JustinsParticle> temp_pList;
      const int pListLengthStart = pList.length();
      if(m_push_type==CYL_CAR) {
         PicSpeciesUtils::stepNormTransferCYLCAR( pList, temp_pList, m_mesh.getdX(),
                                            cnormHalfDt, m_rtol, false );
      }
      else{
         PicSpeciesUtils::stepNormTransfer( pList, temp_pList, m_mesh.getdX(),
                                         cnormHalfDt, m_rtol, false );
      }

      // loop over temp list and move back to main list when converged     
      int iter(1);
      while(temp_pList.length()>0) {
 
         interpolateFieldsToParticles( temp_pList, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual ); 
         addExternalFieldsToParticles( temp_pList, a_em_fields ); 
         applyForces(temp_pList, cnormDt, true);
         m_num_apply_its += temp_pList.length();
         
         if(m_push_type==CYL_CAR) {   
            PicSpeciesUtils::stepNormTransferCYLCAR( pList, temp_pList, m_mesh.getdX(),
                                               cnormHalfDt, m_rtol, true );
         }
         else {
            PicSpeciesUtils::stepNormTransfer( pList, temp_pList, m_mesh.getdX(),
                                               cnormHalfDt, m_rtol, true );
         }

         if(temp_pList.length()==0) break;
         
         if(iter>90) {
            cout << endl;
            cout << "JRA: Picard for particles: iter = " << iter << endl;
            cout << "JRA: m_name = " << m_name << endl;
            cout << "JRA: cnormDt = " << cnormDt << endl;
            Box cell_box = BL.get(dit);
            inspectParticles( temp_pList, cell_box,
                              Efield_inPlane, Bfield_inPlane,
                              Efield_virtual, Bfield_virtual );
         }

         if(iter>=m_iter_max) {
         
            if(m_use_suborbit_model) break;

            int num_not_converged = temp_pList.length();
#if CH_SPACEDIM==1
            // use Newton method to update particles that won't converge with Picard
            if(m_push_type!=CYL_CAR && m_push_type!=CYL_HYBRID) {
              num_not_converged -= newtonParticleUpdate( temp_pList, cnormDt, a_em_fields,
                                                         Efield_inPlane, Bfield_inPlane, 
                                                         Efield_virtual, Bfield_virtual );
            }
#endif
            if ( num_not_converged>0 ) {
               cout << "JRA: Picard for particles: iter = " << iter << endl;
               cout << "JRA: num not converged = " << num_not_converged << endl;
               cout << "JRA: m_name = " << m_name << endl;
               cout << "JRA: cnormDt = " << cnormDt << endl;
               ListIterator<JustinsParticle> lit(temp_pList);
               for(lit.begin(); lit.ok(); ++lit) {
                  cout << "position bar = " << lit().position() << endl;
                  cout << "position old = " << lit().position_old() << endl;
                  cout << "position new = " << 2.0*lit().position()-lit().position_old() << endl;
                  cout << "theta =        " << lit().position_virt(0) << endl;
                  cout << "velocity0 = " << lit().velocity()[0] << endl;
                  cout << "velocity1 = " << lit().velocity()[1] << endl;
                  cout << "velocity2 = " << lit().velocity()[2] << endl;
                  cout << "velocity0 old = " << lit().velocity_old()[0] << endl;
                  cout << "velocity1 old = " << lit().velocity_old()[1] << endl;
                  cout << "velocity2 old = " << lit().velocity_old()[2] << endl;
               }
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

void PicSpecies::advanceInflowParticlesIteratively( const ElectroMagneticFields&  a_em_fields,
                                                    const Real                    a_dt )
{
   CH_TIME("PicSpecies::advanceInflowParticlesIteratively()");
   if(!m_suborbit_inflowJ) return;
   
   // Picard method for coupled half dt advance of particle positions and velocities
   // (inflow particles only. sub-orbits used)
   // xpbar = xpn + dt/2*upbar/(0.5*(gammap_new + gammap_old))
   // upbar = upn + dt/2*q/m*(Ep(xpbar) + upbar/gammap_bar x Bp(xpbar))
   
   const Real cnormDt = a_dt*m_cvac_norm;

   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect& dX(m_mesh.getdX());

   Vector<List<JustinsParticle>>& inflow_list_vect = m_species_bc->getInflowListVect();
             
   const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
   const LevelData<FluxBox>&     Bfield = a_em_fields.getMagneticField();
   const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
   const LevelData<FArrayBox>&     Bfield_virt = a_em_fields.getVirtualMagneticField();
   
   const BoundaryBoxLayoutPtrVect& bdry_layout = m_mesh.getBoundaryLayout();
   for (int b(0); b<bdry_layout.size(); b++) {
      
      const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
      const DisjointBoxLayout& bdry_grids( this_bdry_layout.disjointBoxLayout() );
      const int bdry_dir = this_bdry_layout.dir();
      const int bdry_side(this_bdry_layout.side());
   
      // compute the boundary location
      Real Xbdry;
      if(bdry_side==0) Xbdry = Xmin[bdry_dir];
      else Xbdry = Xmax[bdry_dir]; 

      List<JustinsParticle>& pList = inflow_list_vect[b];
      if(pList.length()==0) continue;      
   
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

         const DataIndex& interior_dit( this_bdry_layout.dataIndex(dit) );
         const EdgeDataBox& Efield_inPlane = Efield[interior_dit];
         const FluxBox&     Bfield_inPlane = Bfield[interior_dit];
         const FArrayBox&   Efield_virtual = Efield_virt[interior_dit].getFab();
         const FArrayBox&   Bfield_virtual = Bfield_virt[interior_dit];

         applyForcesToInflowParticles( pList, a_em_fields, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual, 
                                       bdry_dir, bdry_side, cnormDt ); 
         advanceInflowPositions( pList, bdry_dir, bdry_side, cnormDt ); 
         applyForcesToInflowParticles( pList, a_em_fields, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual, 
                                       bdry_dir, bdry_side, cnormDt ); 

         // transfer not converged particles to a temp list
         List<JustinsParticle> temp_pList;
         const int pListLengthStart = pList.length();      
         PicSpeciesUtils::stepNormTransferInflow( pList, temp_pList, dX, Xbdry,
                                                  bdry_dir, cnormDt, m_rtol, false );
   
         // loop over temp list and move back to main list when converged     
         int iter(1);
         while(temp_pList.length()>0) {
            advanceInflowPositions( temp_pList, bdry_dir, bdry_side, cnormDt ); 
            applyForcesToInflowParticles( temp_pList, a_em_fields, 
                                          Efield_inPlane, Bfield_inPlane, 
                                          Efield_virtual, Bfield_virtual, 
                                          bdry_dir, bdry_side, cnormDt ); 
            
            PicSpeciesUtils::stepNormTransferInflow( pList, temp_pList, dX, Xbdry,
                                                     bdry_dir, cnormDt, m_rtol, true );
            if(temp_pList.length()==0) break;
            
            if(iter>=m_iter_max) {

               int num_not_converged = temp_pList.length();
               if ( num_not_converged>0 ) {
                  cout << "JRA: Picard for inflow particles: iter = " << iter << endl;
                  cout << "JRA: num not converged = " << num_not_converged << endl;
                  ListIterator<JustinsParticle> lit(temp_pList);
                  for(lit.begin(); lit.ok(); ++lit) {
                     cout << "position = " << lit().position() << endl;
                     cout << "velocity = " << lit().velocity()[0] << endl;
                     cout << "position old = " << lit().position_old() << endl;
                     cout << "velocity old = " << lit().velocity_old()[0] << endl;
                  }
               }
               break;

            }
            iter += 1;

         }

         // put back particles that could not converge...     
         ListIterator<JustinsParticle> lit(temp_pList);
         for(lit.begin(); lit.ok();) pList.transfer(lit);
        
         // make sure we got all the particles
         const int pListLengthEnd = pList.length();      
         if(pListLengthEnd!=pListLengthStart) exit(EXIT_FAILURE);

      }
   }

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
   m_species_bc->removeOutflowParticles();
   
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
   if(m_suborbit_inflowJ) return;
   else m_species_bc->injectInflowParticles( m_data );
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

   if(a_restart_file_name.empty()) initializeFromInputFile(a_units);
   else initializeFromRestartFile( a_time, a_restart_file_name );
   
   if(m_push_type==CYL_CAR) {
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      for(dit.begin(); dit.ok(); ++dit) {
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         ListIterator<JustinsParticle> lit(pList);
         for(lit.begin(); lit.ok(); ++lit) {
            std::array<Real,2> position_virt = {1.0,0.0};
            lit().setPositionVirt( position_virt );
         }
      }
   }
   
   int totalParticleCount = m_data.numParticles();
   if(!procID()) {
      cout << "Finished initializing pic species " << m_name  << endl;
      cout << "total particles  = " << totalParticleCount << endl << endl;
   }
   
   // define BCs object for this species
   int verbosity = 1; 
   m_species_bc = new PicSpeciesBC( m_name, m_mass, m_mesh, a_units, verbosity );
      
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
         ppspcIC.get("Z_min", sXmin[1]);
         ppspcIC.get("Z_max", sXmax[1]);
      }
      if(SpaceDim==3) {
         ppspcIC.get("Y_min", sXmin[1]);
         ppspcIC.get("Y_max", sXmax[1]);
         ppspcIC.get("Z_min", sXmin[2]);
         ppspcIC.get("Z_max", sXmax[2]);
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
                  if(SpaceDim==1) partsPerCell[0] = totalPartsPerCell;
                  if(SpaceDim==2) {
                     if(dX[0]<=dX[1]) {
                        partsPerCell[0] = ceil(sqrt((Real)totalPartsPerCell*dX[0]/dX[1]));
                        partsPerCell[1] = ceil((Real)totalPartsPerCell/(Real)partsPerCell[0]);
                     }
                     else {
                        partsPerCell[1] = ceil(sqrt((Real)totalPartsPerCell*dX[1]/dX[0]));
                        partsPerCell[0] = ceil((Real)totalPartsPerCell/(Real)partsPerCell[1]);
                     }
                     CH_assert(partsPerCell.product()>=totalPartsPerCell);
                  }
                  if(SpaceDim==3) {
                     cout << "JRA: uniform_particle_weights not implemented for 3D" << endl;
                     exit(EXIT_FAILURE);
                  }
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
               for(int dir=0; dir<SpaceDim; dir++) {
                  Xpart[dir] += (ipg[dir] + 0.5)*dXpart[dir];
                  if(Xpart[dir]<sXmin[dir] || Xpart[dir]>sXmax[dir]) { 
                     part_outside=true;
                     break;
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
   applyBCs(false);

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

void PicSpecies::setChargeDensityOnNodes()
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
   //SpaceUtils::exchangeNodeFArrayBox(m_chargeDensity_nodes); // needed if more than 1 box per proccesor. bug?
   
   // divide by Jacobian after exchange (corrected Jacobian does not have ghosts)
   const LevelData<NodeFArrayBox>& Jacobian = m_mesh.getCorrectedJnc();  
   for(dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
      this_rho.divide(Jacobian[dit].getFab(),0,0,1);
   }
   
}

void PicSpecies::setCurrentDensity( const bool  a_from_explicit_solver )
{
   CH_TIME("PicSpecies::setCurrentDensity()");
    
   CH_assert(m_charge != 0);
   
   bool cyl_car = false;
   if(m_push_type==CYL_CAR) cyl_car = true;
 
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();
   
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
                                    cyl_car,
                                    a_from_explicit_solver );

   }

   // contribute inflow/outflow current to J if being called from explicit solver.
   // This is needed for charge conserving deposits. Note that J from inflow/outflow 
   // particles for iterative implicit solvers is handled differently.
   if(a_from_explicit_solver) {
      m_species_bc->depositInflowOutflowJ( m_currentDensity, m_currentDensity_virtual,
                                          *m_meshInterp, m_interpJToGrid, true );
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

   // apply BC to J (only used for symmetry BC right now)
   m_species_bc->applyToJ(m_currentDensity,m_currentDensity_virtual);
   
   // add ghost cells to valid cells
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   m_currentDensity.exchange(m_currentDensity.interval(), m_mesh.reverseCopier(), addEdgeOp);
   //SpaceUtils::exchangeEdgeDataBox(m_currentDensity); // needed if more than 1 box per proccesor. bug?
  
   // divide by Jacobian after doing exchange (Corrected Jacobian does not have ghosts) 
   const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();  
   for(dit.begin(); dit.ok(); ++dit) {
      EdgeDataBox& J_inPlane = m_currentDensity[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         J_inPlane[dir].divide(Jec[dit][dir],0,0,1);
      }
   }

#if CH_SPACEDIM<3
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_currentDensity_virtual.exchange(m_currentDensity_virtual.interval(), 
                                     m_mesh.reverseCopier(), addNodeOp);
   //SpaceUtils::exchangeNodeFArrayBox(m_currentDensity_virtual,m_mesh); // needed if more than 1 box per proccesor. bug?
   
   const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();  
   for(dit.begin(); dit.ok(); ++dit) {
      FArrayBox& J_virtual = m_currentDensity_virtual[dit].getFab();
      for (int comp=0; comp<J_virtual.nComp(); ++comp) {
         J_virtual.divide(Jnc[dit].getFab(),0,comp,1);
      }
   }
#endif
   
}

void PicSpecies::setInflowJ( const Real  a_dt )
{
   CH_TIME("PicSpecies::setInflowJ()");
 
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();

   // set containers to zero
   for(dit.begin(); dit.ok(); ++dit) {
      EdgeDataBox& J_inPlane = m_inflowJ[dit];
      for (int dir=0; dir<SpaceDim; ++dir) J_inPlane[dir].setVal(0.0);
#if CH_SPACEDIM<3
      FArrayBox& J_virtual = m_inflowJ_virtual[dit].getFab();
      J_virtual.setVal(0.0);
#endif
   }

   // deposit inflow J
   const Real cnormDt = a_dt*m_cvac_norm;
   m_species_bc->depositInflowJ( m_inflowJ, m_inflowJ_virtual, 
                                *m_meshInterp, cnormDt );

   // multiply by charge/volume
   const Real volume_scale = m_mesh.getVolumeScale();
   for(dit.begin(); dit.ok(); ++dit) {
#if CH_SPACEDIM<3
      m_inflowJ_virtual[dit].getFab().mult(m_charge/volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_inflowJ[dit][dir].mult(m_charge/volume_scale);
      }
   }

   // apply BC to J (only used for symmetry BC right now)
   m_species_bc->applyToJ(m_inflowJ,m_inflowJ_virtual);
   
   // add ghost cells to valid cells
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   m_inflowJ.exchange(m_inflowJ.interval(), m_mesh.reverseCopier(), addEdgeOp);
  
   // divide by Jacobian after doing exchange (Corrected Jacobian does not have ghosts) 
   const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();  
   for(dit.begin(); dit.ok(); ++dit) {
      EdgeDataBox& J_inPlane = m_inflowJ[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         J_inPlane[dir].divide(Jec[dit][dir],0,0,1);
      }
   }

#if CH_SPACEDIM<3
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_inflowJ_virtual.exchange(m_inflowJ_virtual.interval(), 
                              m_mesh.reverseCopier(), addNodeOp);
   
   const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();  
   for(dit.begin(); dit.ok(); ++dit) {
      FArrayBox& J_virtual = m_inflowJ_virtual[dit].getFab();
      for (int comp=0; comp<J_virtual.nComp(); ++comp) {
         J_virtual.divide(Jnc[dit].getFab(),0,comp,1);
      }
   }
#endif
  
}

void PicSpecies::advanceSubOrbitParticlesAndSetJ( const ElectroMagneticFields&  a_em_fields,
		                                  const Real                    a_dt,
                                                  const bool                    a_from_emjacobian )
{
   CH_TIME("PicSpecies::advanceSubOrbitParticlesAndSetJ()");

   bool cyl_car = false;
   if(m_push_type==CYL_CAR) cyl_car = true;
 
   int iter_max = m_iter_max;
   if(m_suborbit_testing) iter_max = std::max(10,iter_max); // used for CYL pusher test
   if(a_from_emjacobian) iter_max += iter_max; // increase iters for buffer

   SpaceUtils::zero( m_suborbitJ );
   SpaceUtils::zero( m_suborbitJ_virtual );

   const Real cnormDt = a_dt*m_cvac_norm;

   const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
   const LevelData<FluxBox>&     Bfield = a_em_fields.getMagneticField();
   const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
   const LevelData<FArrayBox>&     Bfield_virt = a_em_fields.getVirtualMagneticField();
              
   const BoxLayout& BL = m_data_suborbit.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      const EdgeDataBox& Efield_inPlane = Efield[dit];
      const FluxBox&     Bfield_inPlane = Bfield[dit];
      const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
      const FArrayBox& Bfield_virtual = Bfield_virt[dit];

      // define references to suborbitJ containers on this box      
      EdgeDataBox& local_J = m_suborbitJ[dit]; 
      FArrayBox& local_Jv  = m_suborbitJ_virtual[dit].getFab(); 
     
      // define local containers to store J from a single particle
      EdgeDataBox this_Jp(local_J.box(),local_J.nComp());
      FArrayBox this_Jpv(local_Jv.box(),local_Jv.nComp()); // is this correct?
      
      // Loop over each particle in the list and advance/deposit for each suborbit
      ListBox<JustinsParticle>& box_list = m_data_suborbit[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      List<JustinsParticle> temp_pList, finished_pList;

      ListIterator<JustinsParticle> lit(pList);
      for(lit.begin(); lit.ok();) { 
         
     	 const int num_part_suborbits = lit().numSubOrbits();
         int num_suborbits = num_part_suborbits;
         Real cnormDt_sub = cnormDt/num_suborbits;
   
	 // zero containers to store J from a single particle
	 this_Jp.setVal(0.0);
	 this_Jpv.setVal(0.0);

	 // save the old position and velocity
         const RealVect& xpold = lit().position_old();
	 const RealVect xpold0 = xpold;
	 const std::array<Real,3>& vpold = lit().velocity_old();
         const std::array<Real,3> vpold0 = vpold;

	 // set xpbar and vpbar to old values for initial guess of first suborbit
	 lit().setPosition(xpold);
	 lit().setVelocity(vpold);
         if(m_push_type==CYL_CYL || m_push_type==CYL_HYBRID) { // reset dthetap to zero
            std::array<Real,2>& position_virt = lit().position_virt();
            position_virt[0] = 0.0;
	 }
         if(m_push_type==CYL_CAR) { // reset costhp = 1, sinthp = 0
            std::array<Real,2>& position_virt = lit().position_virt();
            position_virt[0] = 1.0;
            position_virt[1] = 0.0;
	 }

	 // put this particle in a list all by itself
         List<JustinsParticle> single_pList;
	 single_pList.transfer(lit);
	 
	 // loop over suborbits for this particle
	 for (int nv = 0; nv<num_suborbits; nv++) {

            int iter(0);
            while(single_pList.length()>0) {
               interpolateFieldsToParticles( single_pList, Efield_inPlane, Bfield_inPlane, 
                                             Efield_virtual, Bfield_virtual ); 
               addExternalFieldsToParticles( single_pList, a_em_fields ); 
               applyForces(single_pList, cnormDt_sub, true);
         
	       if(m_push_type==CYL_CAR) {   
                  PicSpeciesUtils::stepNormTransferCYLCAR( temp_pList, single_pList, m_mesh.getdX(),
                                         cnormDt_sub/2.0, m_rtol, true );
               }
               else {
                  PicSpeciesUtils::stepNormTransfer( temp_pList, single_pList, m_mesh.getdX(),
                                         cnormDt_sub/2.0, m_rtol, true );
               }
         
	       if(single_pList.length()==0) break;

               iter += 1;
	       if( iter >= iter_max ) { // this suborbit did not converge

                  ListIterator<JustinsParticle> lit_single(single_pList);
	          lit_single.begin();
		  if( !a_from_emjacobian ) { // increase # suborbits at start again
		     cout << "JRA: iter_max reached for suborbit nv = " << nv + 1 << " of " 
		          << num_suborbits << endl;
                     num_suborbits++;
		     cout << "     increasing suborbits to " << num_suborbits << endl;
		     lit_single().setPosition(xpold0);
		     lit_single().setOldPosition(xpold0);
		     lit_single().setVelocity(vpold0);
		     lit_single().setOldVelocity(vpold0);
                     lit_single().setNumSubOrbits( num_suborbits );
                     cnormDt_sub = cnormDt/num_suborbits;
	             this_Jp.setVal(0.0);
	             this_Jpv.setVal(0.0);
                     if(m_push_type==CYL_CYL || m_push_type==CYL_HYBRID) { // reset dthetap to zero
                        std::array<Real,2>& position_virt = lit_single().position_virt();
                        position_virt[0] = 0.0;
                     }
                     if(m_push_type==CYL_CAR) { // reset costhp = 1, sinthp = 0
                        std::array<Real,2>& position_virt = lit_single().position_virt();
                        position_virt[0] = 1.0;
                        position_virt[1] = 0.0;
	             }
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
	    if(nv==-1) continue;

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
                                          cyl_car,
                                          false );

	    // 4) convert particle quantities from time-centered to new
            if(m_push_type==CYL_HYBRID) {
	       advanceVelocities_CYL_HYB_2ndHalf( a_dt/num_suborbits, temp_pList );
	       advancePositions_2ndHalf( temp_pList );
	    } 
	    else if(m_push_type==CYL_CAR) {
               advanceParticles_CYL_CAR_2ndHalf( temp_pList );
	    }
	    else {
	       advanceVelocities_2ndHalf( temp_pList );
	       advancePositions_2ndHalf( temp_pList );
	    }

	    // 5) reset xpold and vpold to suborbit xpnew and vpnew
            ListIterator<JustinsParticle> lit_temp(temp_pList);
            for(lit_temp.begin(); lit_temp.ok();) {
               RealVect& xpold = lit_temp().position_old();
	       std::array<Real,3>& vpold = lit_temp().velocity_old();
	       if(nv==num_suborbits-1) {
		  xpold = xpold0;
	          vpold = vpold0;
	          finished_pList.transfer(lit_temp);
	       }
               else { // update old values and put back in single pList
                  const RealVect& xp = lit_temp().position();
	          const std::array<Real,3>& vp = lit_temp().velocity();
	          xpold = xp;
    	          vpold = vp;
	          single_pList.transfer(lit_temp);
	       }
	    }

	 } // end loop over suborbits for this particle

	 // divide particle J by number of suborbits and add it to the total
	 for(int dir=0; dir<SpaceDim; dir++) {
	    this_Jp[dir].divide(num_suborbits);
	    local_J[dir].plus(this_Jp[dir]);
	 }
	 this_Jpv.divide(num_suborbits);
	 local_Jv.plus(this_Jpv);

      } // end loop over pList particles
	 
      // transfer the particles back to the main pList owned by m_data_suborbit
      ListIterator<JustinsParticle> lit_final(finished_pList);
      for(lit_final.begin(); lit_final.ok();) pList.transfer(lit_final);
      
   } // end loop over boxes

   //////////////////////////////////////////////////////////////////////
   //
   //          multiply by constants, set BCs, and exchange
   //
   //////////////////////////////////////////////////////////////////////


   // multiply by charge/volume
   const Real volume_scale = m_mesh.getVolumeScale();
   for(dit.begin(); dit.ok(); ++dit) {
#if CH_SPACEDIM<3
      m_suborbitJ_virtual[dit].getFab().mult(m_charge/volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_suborbitJ[dit][dir].mult(m_charge/volume_scale);
      }
   }

   // apply BC to J (only used for symmetry BC right now)
   m_species_bc->applyToJ(m_suborbitJ,m_suborbitJ_virtual);
   
   // add ghost cells to valid cells
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   m_suborbitJ.exchange(m_suborbitJ.interval(), m_mesh.reverseCopier(), addEdgeOp);
  
   // divide by Jacobian after doing exchange (Corrected Jacobian does not have ghosts) 
   const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();  
   for(dit.begin(); dit.ok(); ++dit) {
      EdgeDataBox& J_inPlane = m_suborbitJ[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         J_inPlane[dir].divide(Jec[dit][dir],0,0,1);
      }
   }

#if CH_SPACEDIM<3
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_suborbitJ_virtual.exchange(m_suborbitJ_virtual.interval(), 
                              m_mesh.reverseCopier(), addNodeOp);
   
   const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();  
   for(dit.begin(); dit.ok(); ++dit) {
      FArrayBox& J_virtual = m_suborbitJ_virtual[dit].getFab();
      for (int comp=0; comp<J_virtual.nComp(); ++comp) {
         J_virtual.divide(Jnc[dit].getFab(),0,comp,1);
      }
   }
#endif

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
                                   const ElectroMagneticFields&     a_em_fields,
                                   const Real                       a_dt ) const
{
   CH_TIME("PicSpecies::accumulateMassMatrices()");
   if( m_charge==0 ) return;

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();
   
   //
   //
   //

   const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();
   const LevelData<FArrayBox>& Bfield_virt = a_em_fields.getVirtualMagneticField();

   const Real cnormDt = a_dt*m_cvac_norm;
   const Real alphas = m_fnorm_const*cnormDt/2.0;
   const Real volume_scale = m_mesh.getVolumeScale();
   const Real qovs = m_charge/volume_scale;

   int inert_type = 0;
   if(m_push_type==CYL_CYL) inert_type = 1;
   if(m_push_type==CYL_HYBRID) inert_type = 0; // zero is correct here
   if(m_push_type==CYL_CAR) inert_type = 2;

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
                                         inert_type,
                                         pList,
                                         m_interpJToGrid );

   }
  
}

void PicSpecies::interpolateEfieldToParticles( const ElectroMagneticFields&  a_em_fields )
{
  if(!m_forces) return;
  if(m_charge==0.0) return;
  CH_TIME("PicSpecies::interpolateEfieldToParticles()");
  
  const int blank_B = 1;
 
  const DisjointBoxLayout& grids(m_mesh.getDBL());
   
  const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
  const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();
#if CH_SPACEDIM < 3
  const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
  const LevelData<FArrayBox>& Bfield_virt = a_em_fields.getVirtualMagneticField();
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

void PicSpecies::interpolateFieldsToParticles( const ElectroMagneticFields&  a_em_fields )
{
  if(!m_forces) return;
  if(m_charge==0.0) return;
  CH_TIME("PicSpecies::interpolateFieldsToParticles()");
   
  const DisjointBoxLayout& grids(m_mesh.getDBL());
   
  const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
  const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();

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

  const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
  const LevelData<FArrayBox>& Bfield_virt = a_em_fields.getVirtualMagneticField();
    
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

void PicSpecies::applyForcesToInflowParticles( List<JustinsParticle>&  a_pList,
                                         const ElectroMagneticFields&  a_em_fields,
                                         const EdgeDataBox&            a_Efield_inPlane,
                                         const FluxBox&                a_Bfield_inPlane,
                                         const FArrayBox&              a_Efield_virtual,
                                         const FArrayBox&              a_Bfield_virtual,
                                         const int                     a_bdry_dir, 
                                         const int                     a_bdry_side,
                                         const Real                    a_cnormDt )
{
   if(a_pList.length()==0) return;   
   CH_TIME("PicSpecies::applyForcesToInflowParticles()");
   
   // set the boundary position
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   Real X0;
   if(a_bdry_side==0) X0 = Xmin[a_bdry_dir];
   else X0 = Xmax[a_bdry_dir]; 
      
   // gather the electric and magnetic fields at the particle
   m_meshInterp->interpolateEMfieldsToInflowParts( a_pList,
                                            a_bdry_dir,
                                            a_bdry_side,
                                            a_Efield_inPlane[0],
#if CH_SPACEDIM >= 2
                                            a_Efield_inPlane[1],
#else
                                            a_Efield_virtual,
#endif
                                            a_Efield_virtual,
                                            a_Bfield_inPlane[0],
#if CH_SPACEDIM >= 2
                                            a_Bfield_inPlane[1],
#else
                                            a_Bfield_virtual,
#endif
                                            a_Bfield_virtual );
   
   // add external fields
   addExternalFieldsToParticles( a_pList, a_em_fields ); 
   
   // advance the velocities to upbar
   Real cnormDt0, cnormDt1;
   Real gammap_old;

   List<JustinsParticle> temp_pList;
   const int pList_length0 = a_pList.length();
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok();) {

      // compute time step for sub-orbit
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upold = lit().velocity_old();

      gammap_old = 1.0; 
#ifdef RELATIVISTIC_PARTICLES
      gammap_old += upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
      gammap_old = sqrt(gammap_old);
#endif

      cnormDt0 = (X0-xpold[a_bdry_dir])*gammap_old/upold[a_bdry_dir];
      cnormDt1 = a_cnormDt - cnormDt0;
      CH_assert(cnormDt1>0.0);
         
      // put this particle in a list all by itself
      List<JustinsParticle> single_pList;
      single_pList.transfer(lit);
   
      // Boris push this particle
      applyForces(single_pList, cnormDt1, true);
      
      // check for reflection and move this particle to the temp list
      ListIterator<JustinsParticle> lit_single(single_pList);
      for(lit_single.begin(); lit_single.ok();) {
         std::array<Real,3>& upbar = lit_single().velocity();
         int sign_fact = 1 - 2*a_bdry_side;
         if(sign_fact*upbar[a_bdry_dir]<=0.0) {
            for(int n=0; n<3; n++) upbar[n] = upold[n];
            upbar[a_bdry_dir] = 0.0;
         }
         temp_pList.transfer(lit_single);
      }

   }
   
   // transfer particles from the temp list back to the passed list 
   ListIterator<JustinsParticle> lit_temp(temp_pList);
   for(lit_temp.begin(); lit_temp.ok();) a_pList.transfer(lit_temp);
   CH_assert(a_pList.length()==pList_length0);

}

void PicSpecies::getFieldDerivativesAtParticle( RealVect&         a_dExdy,
                                                RealVect&         a_dEydy,
                                                RealVect&         a_dEzdy,
                                                RealVect&         a_dBxdy,
                                                RealVect&         a_dBydy,
                                                RealVect&         a_dBzdy,
                                          const JustinsParticle&  a_particle,
                                          const EdgeDataBox&      a_Efield_inPlane,
                                          const FluxBox&          a_Bfield_inPlane,
                                          const FArrayBox&        a_Efield_virtual,
                                          const FArrayBox&        a_Bfield_virtual )
{
   CH_TIME("PicSpecies::getFieldDerivativesAtParticle() single particle");

   m_meshInterp->getFieldDerivativesAtPart( a_dExdy, a_dEydy, a_dEzdy,
                                            a_dBxdy, a_dBydy, a_dBzdy,
                                            a_particle,
                                            a_Efield_inPlane[0],
#if CH_SPACEDIM >= 2
                                            a_Efield_inPlane[1],
#else
                                            a_Efield_virtual,
#endif
                                            a_Efield_virtual,
                                            a_Bfield_inPlane[0],
#if CH_SPACEDIM >= 2
                                            a_Bfield_inPlane[1],
#else
                                            a_Bfield_virtual,
#endif
                                            a_Bfield_virtual,
                                            m_interpEToParts );


}

int PicSpecies::newtonParticleUpdate( List<JustinsParticle>&  a_pList,
                                const Real          a_cnormDt, 
                                const ElectroMagneticFields&  a_em_fields,
                                const EdgeDataBox&  a_Efield_inPlane,
                                const FluxBox&      a_Bfield_inPlane,
                                const FArrayBox&    a_Efield_virtual,
                                const FArrayBox&    a_Bfield_virtual )
{
   CH_TIME("PicSpecies::newtonParticleUpdate()");

   // Newton method to get self-consistent particle position and velocity. 
   // Sometimes Picard iteration method is slow or it fails.
   // I've only seen failures with CC0 scheme at cell crossings where 
   // oscillating solutions can occur. Note that the CC1 scheme reduces 
   // to CC0 near physical boundaries.
   //
   // Only implemented for 1D planar far. 
   // What about axisymmetric ?
   // What about 2D ?
   
   CH_assert(SpaceDim==1);
   CH_assert(m_interpEToParts!=TSC); // derivative calc for TSC not implemented yet
   
   int not_converged_count = 0;   
   
   Real *yp = new Real[m_newton_maxits+1]; 
   Real fp, dfpdy;
   RealVect dExdy, dEydy, dEzdy;
   RealVect dBxdy, dBydy, dBzdy;
   RealVect dbxdy, dbydy, dbzdy;
   
   const Real alphas = m_fnorm_const*a_cnormDt/2.0;
   
   // get boundary info neeed to properly treat particles at 
   // inflow/outflow boundaries
   const IntVect& is_outflow_bc_lo = m_species_bc->getIsOutflowBC_lo(); 
   const IntVect& is_outflow_bc_hi = m_species_bc->getIsOutflowBC_hi(); 
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect& dX(m_mesh.getdX());
 
   Real denom, dbpsq;
   Real gammap_hat, gammap_bar;
#ifdef RELATIVISTIC_PARTICLES
   Real gammap_old, gammap_new, dgpbdy, dgphdy;
   std::array<Real,3> upnew;
#endif
   std::array<Real,3> dupdy, wp, bp, dwpdy, dbpdy;
   Real wpdotbp, dwpdotbp;
            
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      
      const RealVect& xpold = lit().position_old();
      const std::array<Real,3>& upold = lit().velocity_old();
      RealVect& xpbar = lit().position();
      std::array<Real,3>& upbar = lit().velocity();
      std::array<Real,3>& Ep = lit().electric_field();
      std::array<Real,3>& Bp = lit().magnetic_field();
      
#ifdef RELATIVISTIC_PARTICLES
      gammap_old = 1.0 + upold[0]*upold[0] + upold[1]*upold[1] + upold[2]*upold[2];
      gammap_old = sqrt(gammap_old);
      for(int n=0; n<3; n++) upnew[n] = 2.0*upbar[n] - upold[n];
      gammap_new = 1.0 + upnew[0]*upnew[0] + upnew[1]*upnew[1] + upnew[2]*upnew[2];
      gammap_new = sqrt(gammap_new);
      gammap_hat = 0.5*(gammap_old + gammap_new);
#else
      gammap_hat = 1.0;
#endif

      RealVect xpnew;
      
      //  use current value of upbar to set initial value for y = xpnew-xpold 
      for(int dir=0; dir<SpaceDim; dir++) {
         yp[0] = a_cnormDt*upbar[dir]/gammap_hat;
      }
       
      // use newton method to get self-consistent upbar and xpbar
      int n_last = 0, ng_last = 0;
      Real dyp, rel_error;
      Real yp0 = yp[0];

      // guess loop is experimental! 
      Real guess_factor = 1.0;
      //Real damping_factor = 1.0;
      bool converged = false;
      bool is_bdry_part = false;
      for (int ng=0; ng<m_newton_num_guess; ng++) {
      
         guess_factor = 1.0 + 0.1*ng;
         if(ng>19) guess_factor = 1.0 + 0.01*(ng-19);
         if(ng>29) guess_factor = 1.0 + 0.001*(ng-29);
               
         if(is_bdry_part && xpold[0]<Xmin[0]) {
            //if(ng<10) damping_factor = 1.0 - 0.1*ng;
            if(ng>1) { // not ideal to change xpold, but its ok for inflow
                   // particles because it won't affect charge-conservation
                   // and convergence issues are typically only observed
                   // for low-velocity particles that barely cross the bdry.
                   // Still a WIP....
                   // tried using damp_fact>1 and it did help some cases,
                   // but it can also lead to issues via an unphysical xpnew
               //damping_factor = 1.0;
               guess_factor = 1.0;
               RealVect& xpold_2 = lit().position_old();
               Real dxpold = (xpold[0]-Xmin[0])/pow((double)ng,2);
               xpold_2[0] = Xmin[0] + dxpold; // updates xpold[0]
               yp0 = max(Xmin[0] - xpold[0],yp[m_newton_maxits]);
            }
         }
         
         yp[0] = yp0*guess_factor;

         for (int n=0; n<m_newton_maxits; n++) { // Newton iterations
              
            // this is a hack, but it seems to work much better 
            // than using a damping factor 
            if(is_bdry_part && xpold[0]<Xmin[0]) {
               if(ng==1) {
               RealVect& xpold_2 = lit().position_old();
               xpold_2[0] = xpnew[0] - upbar[0]/gammap_hat*a_cnormDt;
               }
            }       
 
            // update xpnew and xpbar 
            xpnew[0] = xpold[0] + yp[n];
            xpbar[0] = (xpold[0] + xpnew[0])/2.0;
         
            // redefine initial guess for boundary crosssing 
            // particle to ensure that it crosses the boundary
            if(is_outflow_bc_lo[0] && !is_bdry_part) {
               if(xpnew[0]<Xmin[0] || xpold[0]<Xmin[0]) {
                  is_bdry_part = true;
                  yp0 = Xmin[0] - xpold[0];
               }
            }
            if(is_outflow_bc_hi[0] && !is_bdry_part) {
               if(xpnew[0]>Xmax[0] || xpold[0]>Xmax[0]) {
                  is_bdry_part = true;
                  yp0 = Xmax[0] - xpold[0];
               }
            }

            // update Ep and Bp using current value of xpbar     
            interpolateFieldsToParticle( lit(), 
                                         a_Efield_inPlane, a_Bfield_inPlane, 
                                         a_Efield_virtual, a_Bfield_virtual );
            if(a_em_fields.externalFields()) {
               std::array<Real,3> extE, extB;
               extE = a_em_fields.getExternalE(xpbar);
               extB = a_em_fields.getExternalB(xpbar);
               for (int n=0; n<3; n++) {
                  Ep[n] += extE[n];
                  Bp[n] += extB[n];
               }
            }

            // compute upbar[0]
      
            // add half acceleration to old velocity
            wp[0] = upold[0] + alphas*Ep[0];
            wp[1] = upold[1] + alphas*Ep[1];
            wp[2] = upold[2] + alphas*Ep[2];
      
            // compute relativistic factor for Lorentz force
            gammap_bar = 1.0;
#ifdef RELATIVISTIC_PARTICLES
            gammap_bar += wp[0]*wp[0] + wp[1]*wp[1] + wp[2]*wp[2];
            gammap_bar = sqrt(gammap_bar);
#endif
            if(m_push_type==CYL_HYBRID) {
           
               const Real dthp = lit().position_virt(0);
               const Real costhp = std::cos(dthp);
               const Real sinthp = std::sin(dthp);
      
               const Real upoldr_2 = costhp*upold[0] + sinthp*upold[1]; 
               const Real upoldth_2 = -sinthp*upold[0] + costhp*upold[1]; 
            
               wp[0] = upoldr_2 + alphas*Ep[0];
               wp[1] = upoldth_2 + alphas*Ep[1];

            } 
   
            // scale Bp by alpha/gammap_bar
            bp[0] = alphas*Bp[0]/gammap_bar;
            bp[1] = alphas*Bp[1]/gammap_bar;
            bp[2] = alphas*Bp[2]/gammap_bar;
            if(m_push_type==CYL_CYL) { 
               const Real dtheta = lit().position_virt(0);
               bp[2] = bp[2] + sin(dtheta);
            }
            denom = 1.0 + bp[0]*bp[0] + bp[1]*bp[1] + bp[2]*bp[2];
            wpdotbp = wp[0]*bp[0] + wp[1]*bp[1] + wp[2]*bp[2];

            //
            //
            //

            // update upbar
            upbar[0] = ( wp[0] + wp[1]*bp[2] - wp[2]*bp[1] + wpdotbp*bp[0] )/denom;
            upbar[1] = ( wp[1] + wp[2]*bp[0] - wp[0]*bp[2] + wpdotbp*bp[1] )/denom;
            upbar[2] = ( wp[2] + wp[0]*bp[1] - wp[1]*bp[0] + wpdotbp*bp[2] )/denom;
  
            // update dEp/dy and dBp/dy using current value of xpbar     
            getFieldDerivativesAtParticle( dExdy, dEydy, dEzdy, dBxdy, dBydy, dBzdy,
                                           lit(), a_Efield_inPlane, a_Bfield_inPlane, 
                                           a_Efield_virtual, a_Bfield_virtual ); 
            
            int dir0 = 0;
            if(a_em_fields.externalFields()) {
               std::array<Real,3> extdEdX, extdBdX;
               extdEdX = a_em_fields.getExternaldEdX(xpbar,dir0);
               dExdy[0] += 0.5*extdEdX[0];
               dEydy[0] += 0.5*extdEdX[1];
               dEzdy[0] += 0.5*extdEdX[2];
               extdBdX = a_em_fields.getExternaldBdX(xpbar,dir0);
               dBxdy[0] += 0.5*extdBdX[0];
               dBydy[0] += 0.5*extdBdX[1];
               dBzdy[0] += 0.5*extdBdX[2];
            }
            
            // compute dwpdy
            dwpdy[0] = alphas*dExdy[dir0];
            dwpdy[1] = alphas*dEydy[dir0];
            dwpdy[2] = alphas*dEzdy[dir0];

            // compute dbpdy
#ifdef RELATIVISTIC_PARTICLES
            dgpbdy = (wp[0]*dwpdy[0] + wp[1]*dwpdy[1] + wp[2]*dwpdy[2] )/gammap_bar;
            dbpdy[0] = (alphas*dBxdy[dir0] - bp[0]*dgpbdy)/gammap_bar;
            dbpdy[1] = (alphas*dBydy[dir0] - bp[1]*dgpbdy)/gammap_bar;
            dbpdy[2] = (alphas*dBzdy[dir0] - bp[2]*dgpbdy)/gammap_bar;
#else
            dbpdy[0] = alphas*dBxdy[dir0]/gammap_bar;
            dbpdy[1] = alphas*dBydy[dir0]/gammap_bar;
            dbpdy[2] = alphas*dBzdy[dir0]/gammap_bar;
#endif
            
            // compute dupdy
            dwpdotbp = dwpdy[0]*bp[0] + dwpdy[1]*bp[1] + dwpdy[2]*bp[2]
                     + wp[0]*dbpdy[0] + wp[1]*dbpdy[1] + wp[2]*dbpdy[2];
            dbpsq = 2.0*(dbpdy[0]*bp[0] + dbpdy[0]*bp[1] + dbpdy[0]*bp[2]);

            dupdy[0] = ( dwpdy[0] + dwpdy[1]*bp[2] - dwpdy[2]*bp[1] 
                                  + wp[1]*dbpdy[2] - wp[2]*dbpdy[1]
                     +   dwpdotbp*bp[0] + wpdotbp*dbpdy[0]
                     -   upbar[0]*dbpsq )/denom;
            
            dupdy[1] = ( dwpdy[1] + dwpdy[2]*bp[0] - dwpdy[0]*bp[2] 
                                  + wp[2]*dbpdy[0] - wp[0]*dbpdy[2]
                     +   dwpdotbp*bp[1] + wpdotbp*dbpdy[1]
                     -   upbar[1]*dbpsq )/denom;
            
            dupdy[2] = ( dwpdy[2] + dwpdy[0]*bp[1] - dwpdy[1]*bp[0] 
                                  + wp[0]*dbpdy[1] - wp[1]*dbpdy[0]
                     +   dwpdotbp*bp[2] + wpdotbp*dbpdy[2]
                     -   upbar[2]*dbpsq )/denom;

            // compute fp and dfp/dy
            fp = yp[n] - a_cnormDt*upbar[0]/gammap_hat;
#ifdef RELATIVISTIC_PARTICLES
            dgphdy = (upnew[0]*dupdy[0] + upnew[1]*dupdy[1] + upnew[2]*dupdy[2])/gammap_new;
            dfpdy = 1.0 - a_cnormDt/gammap_hat*(dupdy[0] - upbar[0]*dgphdy/gammap_hat);
#else
            dfpdy = 1.0 - a_cnormDt*dupdy[0]/gammap_hat;
#endif

            // update yp with Newton correction
            dyp = -fp/dfpdy;
            yp[n+1] = yp[n] + dyp;
            //yp[n+1] = yp[n] + dyp*damping_factor;
         
            // compute relative tolerance
            n_last = n;
            rel_error = abs(dyp/dX[0]);
            if(rel_error<m_rtol) {
               converged = true;
               break;
            }

         } // end Newton iterations
      
         ng_last = ng;
         if(converged) break;

      } // end initial guess loop
 

      if(n_last>std::round(m_newton_maxits/2)) cout << "particle newton iter = " << n_last << endl;
      if(ng_last>1) cout << "particle guess iter = " << ng_last << endl;
      //if(ng_last>1) cout << "damping factor = " << damping_factor << endl;
      if(!converged) {
      
         not_converged_count += 1;
         cout << "JRA: newton not converged " << endl;
         cout << "rel_error = " << rel_error << endl;
         
         // compute data needed to run local newton solve in matlab
         // to explore particles that won't converge
         IntVect index_old_stag, index_new_stag;
         IntVect index_lo_stag, index_up_stag, index_lo, index_up;
         Real Ex0, Ex1, By0, By1, Ez0, Ez1, Ez2;
         for(int dir=0; dir<SpaceDim; dir++) {
            xpnew[dir] = xpold[dir] + yp0;
            xpbar[dir] = (xpold[dir] + xpnew[dir])/2.0;
            index_old_stag[dir] = floor((xpold[dir] - Xmin[dir])/dX[dir]);
            index_new_stag[dir] = floor((xpnew[dir] - Xmin[dir])/dX[dir]);
            index_lo[dir] = floor((xpbar[dir] - Xmin[dir] - 0.5*dX[dir])/dX[dir]);
            index_up[dir] = index_lo[dir]+1;
            index_lo_stag[dir] = floor((xpbar[dir] - Xmin[dir])/dX[dir]);
            index_up_stag[dir] = index_lo_stag[dir]+1;
            Ex0 = a_Efield_inPlane[dir](index_lo,0);
            Ex1 = a_Efield_inPlane[dir](index_up,0);
            By0 = a_Bfield_virtual(index_lo,0);
            By1 = a_Bfield_virtual(index_up,0);
            Ez0 = a_Efield_virtual(-IntVect::Unit,1);
            Ez1 = a_Efield_virtual(IntVect::Zero,1);
            Ez2 = a_Efield_virtual(IntVect::Unit,1);
         }
         cout << std::setprecision(16) << std::scientific << endl;
         cout << "index_old = " << index_old_stag[0] << endl;
         cout << "index_new = " << index_new_stag[0] << endl;
         cout << "index_lo = " << index_lo[0] << endl;
         cout << "index_lo_stag = " << index_lo_stag[0] << endl;
         cout << "xpold = " << xpold[0] << endl;
         cout << "xpbar = " << xpbar[0] << endl;
         //cout << "cnormDt = " << a_cnormDt << endl;
         //cout << "alphas = " << alphas << endl;
         //cout << "yp0 = " << yp0 << ";" << endl;
         //cout << "xpold = " << xpold[0] << ";" << endl;
         //cout << "upold_x = " << upold[0] << ";" << endl;
         //cout << "upold_y = " << upold[1] << ";" << endl;
         //cout << "upold_z = " << upold[2] << ";" << endl;
         //cout << "Ex0 = " << Ex0 << ";" << endl;
         //cout << "Ex1 = " << Ex1 << ";" << endl;
         //cout << "Ez0 = " << Ez0 << ";" << endl;
         //cout << "Ez1 = " << Ez1 << ";" << endl;
         //cout << "Ez2 = " << Ez2 << ";" << endl;
         //cout << "By0 = " << By0 << ";" << endl;
         //cout << "By1 = " << By1 << ";" << endl;
         cout << std::setprecision(3) << std::scientific << endl;

      }
     
   } // end loop over particle list
   
   // only in-plane vpbar is set above. apply Forces
   // to update all components of vpbar
   applyForces(a_pList, a_cnormDt, true);
  
   delete[] yp;

   int num_converged = a_pList.length()-not_converged_count;
   return num_converged;
 
}

void PicSpecies::addExternalFieldsToParticles( const ElectroMagneticFields&  a_em_fields )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   if(!a_em_fields.externalFields()) return;
   CH_TIME("PicSpecies::addExternalFieldsToParticles()");
   
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      addExternalFieldsToParticles( pList, a_em_fields );

   }
   
}

void PicSpecies::addExternalFieldsToParticles( List<JustinsParticle>&  a_pList, 
                                         const ElectroMagneticFields&  a_em_fields )
{
   if(!m_forces) return;
   if(m_charge==0.0) return;
   if(!a_em_fields.externalFields()) return;
   CH_TIME("PicSpecies::addExternalFieldsToParticles() from particle list");
   
   std::array<Real,3> extE;
   std::array<Real,3> extB;

   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {

      // get external field values at particle position 
      const RealVect& Xp = lit().position();
      extE = a_em_fields.getExternalE(Xp);
      extB = a_em_fields.getExternalB(Xp);

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

void PicSpecies::inspectParticles( List<JustinsParticle>&  a_pList,
                             const Box&          a_cell_box,
                             const EdgeDataBox&  a_Efield_inPlane,
                             const FluxBox&      a_Bfield_inPlane,
                             const FArrayBox&    a_Efield_virtual,
                             const FArrayBox&    a_Bfield_virtual ) const
{
            
   ListIterator<JustinsParticle> lit(a_pList);
   for(lit.begin(); lit.ok(); ++lit) {
      cout << "position bar = " << lit().position() << endl;
      cout << "position old = " << lit().position_old() << endl;
      cout << "position new = " << 2.0*lit().position()-lit().position_old() << endl;
      cout << "theta =        " << lit().position_virt(0) << endl;
      cout << "vpbar = " << lit().velocity()[0] << ", " << lit().velocity()[1] << ", " << lit().velocity()[2] << endl;
      cout << "vpold = " << lit().velocity_old()[0] << ", " << lit().velocity_old()[1] << ", " << lit().velocity_old()[2] << endl;
   }
      
   //Box cell_box = a_Efield_inPlane[0].box();
   Box cell_box = a_cell_box;
   cell_box.grow(1);
   Box node_box = surroundingNodes(cell_box);
   cout << " fnorm = " << m_fnorm_const << endl;
   cout << " cell_box = " << cell_box << endl;
   cout << " node_box = " << node_box << endl;
   BoxIterator gbit(cell_box);
   BoxIterator gnbit(node_box);

   cout << "Er = ";
   for(gbit.begin(); gbit.ok(); ++gbit) {
      const IntVect ig = gbit(); // grid index
      cout << a_Efield_inPlane[0].get(ig,0) << ", ";
   }
   cout << endl;

   cout << "Eth = ";
   for(gnbit.begin(); gnbit.ok(); ++gnbit) {
      const IntVect ig = gnbit(); // grid index
      cout << a_Efield_virtual.get(ig,0) << ", ";
   }
   cout << endl;
            
   cout << "Ez = ";
   for(gnbit.begin(); gnbit.ok(); ++gnbit) {
      const IntVect ig = gnbit(); // grid index
      cout << a_Efield_virtual.get(ig,1) << ", ";
   }
   cout << endl;
            
   cout << "Br = ";
   for(gnbit.begin(); gnbit.ok(); ++gnbit) {
      const IntVect ig = gnbit(); // grid index
      cout << a_Bfield_inPlane[0].get(ig,0) << ", ";
   }
   cout << endl;

   cout << "Bth = ";
   for(gbit.begin(); gbit.ok(); ++gbit) {
      const IntVect ig = gbit(); // grid index
      cout << a_Bfield_virtual.get(ig,0) << ", ";
   }
   cout << endl;
            
   cout << "Bz = ";
   for(gbit.begin(); gbit.ok(); ++gbit) {
      const IntVect ig = gbit(); // grid index
      cout << a_Bfield_virtual.get(ig,1) << ", ";
   }
   cout << endl;

}

bool PicSpecies::isSpecies( const string&  a_name ) const
{
   if(name() == a_name) return true;
   return false;
}


#include "NamespaceFooter.H"

