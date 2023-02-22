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
     m_use_axisymmetric_boris(false), 
     m_axisymmetric_iter_max(1), 
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
      m_use_axisymmetric_boris = true;
      a_ppspc.query("use_axisymmetric_boris",m_use_axisymmetric_boris);
      if(!m_use_axisymmetric_boris) {
         m_axisymmetric_iter_max = 2;
         a_ppspc.query("axisymmetric_iter_max", m_axisymmetric_iter_max);
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
   m_deposit_bdry_J = false;
   m_suborbit_inflowJ = false;
   if(m_charge_conserving_deposit) {
      m_interp_bc_check = true;
      m_deposit_bdry_J = true;
   }
   a_ppspc.query( "interp_bc_check", m_interp_bc_check );
   a_ppspc.query( "deposit_outflow_J", m_deposit_bdry_J );
   a_ppspc.query( "suborbit_inflow_J", m_suborbit_inflowJ );

   a_ppspc.get( "mass", m_mass );
   a_ppspc.get( "charge", m_charge ); 
   a_ppspc.query( "potential", m_Uint );
   a_ppspc.query( "motion", m_motion );
   a_ppspc.query( "forces", m_forces );
   a_ppspc.query( "scatter", m_scatter );
   a_ppspc.query( "write_all_comps", m_write_all_part_comps );

   if(m_charge==0) m_forces = false;

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
      cout << "  deposit outflow J = " << m_deposit_bdry_J << endl;
      cout << "  suborbit inflow J = " << m_suborbit_inflowJ << endl;
      if(m_mesh.axisymmetric()) cout << "  axisym Boris = " << m_use_axisymmetric_boris << endl;
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
   
   m_density.define(grids,1,ghostVect);
   m_density_binfab.define(grids,1,IntVect::Zero);
   m_momentum.define(grids,3,ghostVect);
   m_energy.define(grids,3,ghostVect); // direction-dependent energy

   m_currentDensity.define(grids,1,ghostVect);
   m_inflowJ.define(grids,1,ghostVect);
#if CH_SPACEDIM<3
   m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect);
   m_inflowJ_virtual.define(grids,3-SpaceDim,ghostVect);
#endif
   m_chargeDensity.define(grids,1,ghostVect);
   m_chargeDensity_faces.define(grids,1,ghostVect);
   m_chargeDensity_nodes.define(grids,1,ghostVect);

   m_temperature.define(grids,3,ghostVect);
   m_velocity.define(grids,3,ghostVect);
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_density[dit].setVal(0.0);
      m_density_binfab[dit].setVal(0.0);
      m_momentum[dit].setVal(0.0);
      m_energy[dit].setVal(0.0);
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

void PicSpecies::repositionOutcastsAndApplyForces( const ElectroMagneticFields&  a_em_fields,
                                                   const Real                    a_dt, 
                                                   const bool                    a_byHalfDt )
{
   CH_TIME("PicSpecies::repositionOutcastsAndApplyForces()");
    
   if(!m_forces) return;

   // this function is only called by the iterative implicit solvers for particles that cross
   // a physical boundary with the outflow bc during the iterations. Some may end up being 
   // outcasts in the end, but to be sure the particles are placed back on the boundary and
   // their velocity is updated prior to the next position advance. For each iteration,
   // if they are out of the boundary after the position advance then they do not contribute
   // to the current
   //
   // Note that calling remapOutcast() is expensive - even if there are no outcast particles.
   // Only call if needed.
   
   const Real cnormDt = a_dt*m_cvac_norm;

   // reposition the outcasts that are out of bounds to be just inside the domain boundaries
   List<JustinsParticle>& outcast_pList = m_data.outcast();
   m_species_bc->repositionOutflowParticles(outcast_pList);

   // interpolate fields to particles and push them 
   interpolateFieldsToOutcasts( a_em_fields );
   addExternalFieldsToParticles( outcast_pList, a_em_fields ); 
   applyForces(outcast_pList, cnormDt, a_byHalfDt);

   const int local_closed = m_data.isClosed();
   int global_closed = local_closed;
#ifdef CH_MPI
   MPI_Allreduce( &local_closed, &global_closed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD ); 
#endif

   // remap the outcasts
   if(global_closed==0) {
      m_data.remapOutcast();
      CH_assert(m_data.isClosed());
   }
 
}

void PicSpecies::applyForces( List<JustinsParticle>&  a_pList,
                        const Real   a_cnormDt, 
                        const bool   a_byHalfDt )
{
   CH_TIME("PicSpecies::applyForces()");
    
   if(a_pList.length()==0) return;   
   
   if(m_mesh.axisymmetric() && !m_use_axisymmetric_boris) {
      PicSpeciesUtils::applyForcesAxisymm( a_pList, m_fnorm_const, a_cnormDt,
                                           a_byHalfDt, m_mesh.anticyclic(),
                                           m_axisymmetric_iter_max );
   }
   else {
      PicSpeciesUtils::applyForces( a_pList, m_fnorm_const, a_cnormDt, 
                                    a_byHalfDt, m_mesh.anticyclic() ); 
   }

}

void PicSpecies::advancePositions( const Real  a_full_dt,
                                   const bool  a_half_step )
{
   CH_TIME("PicSpecies::advancePositions()");
  
   if(!m_motion) return;

   const Real cnormDt = m_cvac_norm*a_full_dt;
   Real dt_factor = 1.0;
   if(a_half_step) dt_factor = 0.5; 

   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      // loop over particles in this box and advance
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      advancePositions( pList, cnormDt*dt_factor ); 
      
   }

}

void PicSpecies::advancePositions( List<JustinsParticle>&  a_pList,
                             const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advancePositions() from particle list");
         
   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {
      RealVect& xp = li().position();
      const RealVect& xpold = li().position_old();
      const std::array<Real,3>& betap = li().velocity();
      for(int dir=0; dir<SpaceDim; dir++) xp[dir] = xpold[dir] + betap[dir]*a_cnormDt;
   }
      
}

void PicSpecies::advanceInflowPositions( List<JustinsParticle>&  a_pList,
                                   const int                     a_bdry_dir,
                                   const int                     a_bdry_side,
                                   const Real                    a_cnormDt )
{
   CH_TIME("PicSpecies::advanceInflowPositions() from particle list");
         
   // set the boundary position
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   RealVect X0;
   if(a_bdry_side==0) X0[a_bdry_dir] = Xmin[a_bdry_dir];
   else X0[a_bdry_dir] = Xmax[a_bdry_dir]; 
   
   Real cnormDt0, cnormDt1;

   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {
      RealVect& xp = li().position();
      const RealVect& xpold = li().position_old();
      const std::array<Real,3>& vpold = li().velocity_old();
      cnormDt0 = (X0[a_bdry_dir]-xpold[a_bdry_dir])/vpold[a_bdry_dir];
      for(int dir=0; dir<SpaceDim; dir++) {
         if(dir==a_bdry_dir) continue;
         X0[dir] = xpold[dir] + vpold[dir]*cnormDt0;
      }
      cnormDt1 = a_cnormDt - cnormDt0;
      CH_assert(cnormDt1>0.0);
      const std::array<Real,3>& vpbar = li().velocity();
      for(int dir=0; dir<SpaceDim; dir++) {
         xp[dir] = X0[dir] + vpbar[dir]*cnormDt1;
         xp[dir] = (xpold[dir] + xp[dir])/2.0; // convert to xpbar
      }
      if(vpbar[a_bdry_dir]==0.0) xp[a_bdry_dir] = xpold[a_bdry_dir];
   }
      
}

void PicSpecies::advancePositions_2ndHalf()
{
   CH_TIME("PicSpecies::advancePositions_2ndHalf()");
   
   if(!m_motion) return;
    
   // Each proc loops over its own boxes
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         const RealVect& xpold = li().position_old();
         RealVect& xp = li().position();

         // update particle position
         for(int dir=0; dir<SpaceDim; dir++) {
            xp[dir] = 2.0*xp[dir] - xpold[dir];
         }

      } // end loop over particle list
      
   } // end loop over boxes

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
   CH_assert(m_data.isClosed()); // make sure all particles are accounted for
   
}

void PicSpecies::advanceVelocities( const Real  a_full_dt, 
                                    const bool  a_half_step )
{
   CH_TIME("PicSpecies::advanceVelocities()");

   if(!m_forces) return;
   
   const Real cnormDt = a_full_dt*m_cvac_norm;
 
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      applyForces(pList, cnormDt, a_half_step);

   }

}

void PicSpecies::advanceVelocities_2ndHalf()
{
   CH_TIME("PicSpecies::advanceVelocities_2ndHalf()");
         
   if(!m_forces) return;
   
   // convert velocities from t_{n+1/2} to t_{n+1}
 
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         std::array<Real,3>& vp = li().velocity(); // actually beta
         const std::array<Real,3>& vpold = li().velocity_old(); // actually beta

         vp[0] = 2.0*vp[0] - vpold[0];
         vp[1] = 2.0*vp[1] - vpold[1];
         vp[2] = 2.0*vp[2] - vpold[2];
      
      }
      
   }

}

void PicSpecies::applyInertialForces( const Real  a_dt, 
                                      const bool  a_use_avg_velocity,
                                      const bool  a_update_positions,
                                      const bool  a_half_positions )
{
   CH_TIME("PicSpecies::applyInertialForces()");

   // See Delzanno JCP 2013
   
   if(!m_mesh.axisymmetric() || !m_use_axisymmetric_boris) return;

   const Real cnormDt = a_dt*m_cvac_norm;
   
   Real rp, beta_x, beta_y, x_car, y_car, beta_r, beta_th;
   
   int th_dir = 1;
   if(m_mesh.anticyclic()) th_dir = 2;
 
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         std::array<Real,3>& betap = li().velocity();
         const RealVect& xpold = li().position_old();

         // correct vr and vth velocities to account for inertial forces
         if(a_use_avg_velocity) {
            const std::array<Real,3>& betap_old = li().velocity_old();
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
            RealVect& xp = li().position();
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
   
   // xpbar = xpn + dt/2*vpbar
   // vpbar = vpn + dt/2*q/m*(Ep(xpbar) + vpbar x Bp(xpbar))
   
   if(!m_deposit_bdry_J) { 
      repositionOutcastsAndApplyForces( a_em_fields, a_dt, true );
   }
           
   if (m_iter_order_swap) {
      advancePositions( a_dt, true );
      if(!m_deposit_bdry_J) applyBCs(true);
   }
   interpolateFieldsToParticles( a_em_fields );
   addExternalFieldsToParticles( a_em_fields ); 
   advanceVelocities( a_dt, true );
   if (!m_iter_order_swap) {
      advancePositions( a_dt, true );
      if(!m_deposit_bdry_J) applyBCs(true);
   }

}

void PicSpecies::advanceParticlesIteratively( const ElectroMagneticFields&  a_em_fields,
                                              const Real                    a_dt )
{
   CH_TIME("PicSpecies::advanceParticlesIteratively()");
   
   // Picard method for coupled half dt advance of particle positions and velocities
   //
   // xpbar = xpn + dt/2*vpbar
   // vpbar = vpn + dt/2*q/m*(Ep(xpbar) + vpbar x Bp(xpbar))

   if(m_iter_max==0 || !m_motion || !m_forces) {
      advanceParticles( a_em_fields, a_dt );
      return;
   }
   
   const Real cnormDt = a_dt*m_cvac_norm;
   const Real cnormHalfDt = 0.5*cnormDt;

   if(!m_deposit_bdry_J) { 
      repositionOutcastsAndApplyForces( a_em_fields, a_dt, true );
   }

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
    
      // update xpbar and vpbar 
      if(m_iter_order_swap) advancePositions( pList, cnormHalfDt ); 
      interpolateFieldsToParticles( pList, 
                                    Efield_inPlane, Bfield_inPlane, 
                                    Efield_virtual, Bfield_virtual ); 
      addExternalFieldsToParticles( pList, a_em_fields ); 
      applyForces(pList, cnormDt, true);
      m_num_apply_its += pList.length();
     
      // for high precision, sometime need to do at least two force applys
      // (trying to do this with rtol gets too close to machine precision)
      if(m_iter_min_two) { 
         advancePositions( pList, cnormHalfDt ); 
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
      PicSpeciesUtils::stepNormTransfer( pList, temp_pList, m_mesh.getdX(),
                                         cnormHalfDt, m_rtol, false );
   
      // loop over temp list and move back to main list when converged     
      int iter(1);
      while(temp_pList.length()>0) {
 
         advancePositions( temp_pList, cnormHalfDt ); 
         interpolateFieldsToParticles( temp_pList, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual ); 
         addExternalFieldsToParticles( temp_pList, a_em_fields ); 
         applyForces(temp_pList, cnormDt, true);
         m_num_apply_its += temp_pList.length();
            
         PicSpeciesUtils::stepNormTransfer( pList, temp_pList, m_mesh.getdX(),
                                            cnormHalfDt, m_rtol, true );
   
         if(temp_pList.length()==0) break;
            
         if(iter>=m_iter_max) {

            int num_not_converged = temp_pList.length();

#if CH_SPACEDIM==1
            // use Newton method to update particles that won't converge with Picard
            num_not_converged -= newtonParticleUpdate( temp_pList, cnormDt,
                                                       Efield_inPlane, Bfield_inPlane, 
                                                       Efield_virtual, Bfield_virtual );
#endif
            if ( num_not_converged>0 ) {
               cout << "JRA: Picard for particles: iter = " << iter << endl;
               cout << "JRA: num not converged = " << num_not_converged << endl;
            }
            break;

         }
         iter += 1;

      }

      // put back particles that could not converge...     
      ListIterator<JustinsParticle> li(temp_pList);
      for(li.begin(); li.ok();) pList.transfer(li);
        
      // call advance positions last to enforce charge conservation ?
      //advancePositions( pList, cnormHalfDt );
         
      // make sure we got all the particles
      const int pListLengthEnd = pList.length();      
      if(pListLengthEnd!=pListLengthStart) exit(EXIT_FAILURE);
               
   } // end loop over boxes on this proc
 
   // JRA, need to figure out better way to handle particles crossing
   // outlet boundaries so that I can remove this BC call here.
   // 
   // If call to applyBCs() is made, then outflow particles will be moved to an outflow list. 
   // They will contribute to J in setCurrentDensity() with m_deposit_bdry_J = true. But 
   // they will not contribute to J if using mass matrices. I initially tried to add function
   // to deposit outflow particles to mass matrices, but for iterative solver I also need
   // functions to iteratively advance the particles in the outflow list. The easier (and 
   // more computationally efficient) thing to do is to never put them in the outflow list
   // to begin with during the intermediate iterative stage. Then the outflow particles
   // are retained in the standard particle list and they contribute to J and MM same as
   // all other particles. Thus, if m_deposit_bdry_J is true, then we do not call applyBCs(). 
   //
   // This flag is confusing, but it is meant to be a temporary flag. End goal is to never 
   // call applyBCs during iterative half advance in any situation.
   //
   if(!m_deposit_bdry_J) applyBCs(true);
 
}

void PicSpecies::advanceInflowParticlesIteratively( const ElectroMagneticFields&  a_em_fields,
                                                    const Real                    a_dt )
{
   CH_TIME("PicSpecies::advanceInflowParticlesIteratively()");
   if(!m_suborbit_inflowJ) return;
   
   // Picard method for coupled half dt advance of particle positions and velocities
   // (inflow particles only)
   //
   // xpbar = xpn + dt/2*vpbar
   // vpbar = vpn + dt/2*q/m*(Ep(xpbar) + vpbar x Bp(xpbar))
   
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
         ListIterator<JustinsParticle> li(temp_pList);
         for(li.begin(); li.ok();) pList.transfer(li);
        
         // make sure we got all the particles
         const int pListLengthEnd = pList.length();      
         if(pListLengthEnd!=pListLengthStart) exit(EXIT_FAILURE);

         // advanceInflowPositions( pList, bdry_dir, bdry_side, cnormDt ); 

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
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         RealVect& xpold = li().position_old();
         const RealVect& xp = li().position();

         for(int dir=0; dir<SpaceDim; dir++) {
            xpold[dir] = xp[dir];
         }

      } // end loop over particle list
      
   } // end loop over boxes

}

void PicSpecies::updateOldParticleVelocities()
{
   CH_TIME("PicSpecies::updateOldParticleVelocities()");
   
   if(!m_forces) return;
    
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         const std::array<Real,3>& vp = li().velocity();  // actually beta
         li().setOldVelocity(vp);

      } // end loop over particle list
      
   } // end loop over boxes

}

void PicSpecies::setStableDt()
{
   CH_TIME("PicSpecies::setStableDt()");
   
   // set the stable time step based on particles crossing a grid cell
   const RealVect& dX(m_mesh.getdX());
   Real maxDtinv = 0.0;
   Real thisDtinv;

   // Each proc loops over its own boxes
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      //ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = m_data[dit].listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         const std::array<Real,3>&  v = li().velocity(); // actually beta
         for(int dir=0; dir<SpaceDim; dir++) {
            thisDtinv = abs(v[dir])/dX[dir];
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
   m_volume_scale = a_units.getScale(a_units.VOLUME);

   if(a_restart_file_name.empty()) initializeFromInputFile(a_units);
   else initializeFromRestartFile( a_time, a_restart_file_name );
   
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
      Real cellVolume = m_volume_scale*dX.product();
      const Real numDen_scale = a_units.getScale(a_units.NUMBER_DENSITY);
      Real pWeight = 0.0; 
      
      // query for uniform particle weight option
      bool uniform_particle_weights = false;
      Real pWeight_fixed = 0.0; 
      ppspcIC.query("uniform_particle_weights",uniform_particle_weights);
      if(uniform_particle_weights) {
         pWeight_fixed = minimumParticleWeight( m_density, sXmin, sXmax,
                                                numDen_scale, cellVolume, totalPartsPerCell );
      }
      
      // query for casting weight to float 
      bool float_weight = false;
      ppspcIC.query("use_float_for_weights",float_weight);

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
      uint64_t ID = procID()*512 + 1; // hack for testing purposes
      //Real ID = 0.;
   
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
               ID = ID + 1;
 
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

#else
      MayDay::Error("restart only defined with hdf5");
#endif

}

void PicSpecies::setNumberDensity()
{
   CH_TIME("PicSpecies::setNumberDensity()");
    
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
      box_rho.divide(m_volume_scale);

   }
   m_density.exchange();
     
}

void PicSpecies::setNumberDensityFromBinFab()
{
   CH_TIME("PicSpecies::setNumberFromBinFab()");
   
   JustinsParticle* this_part_ptr = NULL;  
    
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   const Real dV_mapped = m_mesh.getMappedCellVolume();
   const Real dV_phys = dV_mapped*m_volume_scale;

   const DisjointBoxLayout& grids = m_density_binfab.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density_binfab[dit];
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
            const Real& this_w = this_part_ptr->weight();
            rho0 += this_w;
         }
         box_rho.set(ig,0,rho0);

      }
 
      const FArrayBox& box_Ja = Jacobian[dit];
      box_rho.divide(box_Ja,0,0,1);
      box_rho.divide(dV_phys);

   }

   this_part_ptr = NULL;
   delete this_part_ptr;
     
}

void PicSpecies::setMomentumDensity()
{
   CH_TIME("PicSpecies::setMomentumDensity()");
    
   const DisjointBoxLayout& grids = m_momentum.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_mom = m_momentum[dit];
      box_mom.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = momentum;
      m_meshInterp->moment( box_mom,
                            box_list.listItems(),
                            m_mass/m_volume_scale,
                            thisMoment );
      
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_mom.nComp(); ++n) { 
         box_mom.divide(box_Ja,0,n,1);
      }
 
   }
   m_momentum.exchange(); 
     
}

void PicSpecies::setEnergyDensity()
{
   CH_TIME("PicSpecies::setEnergyDensity()");
    
   const DisjointBoxLayout& grids = m_energy.disjointBoxLayout();
   const LevelData<FArrayBox>& Jacobian = m_mesh.getJcc();  
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& box_ene = m_energy[dit];
      box_ene.setVal(0.0);
      for (auto n=0; n<m_energy.nComp(); n++) {
         SpaceUtils::setVal(box_ene,0.0,n);
      }
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = energy;
      m_meshInterp->moment( box_ene,
                            box_list.listItems(),
                            m_mass/m_volume_scale,
                            thisMoment ); 
      
      const FArrayBox& box_Ja = Jacobian[dit];
      for (auto n=0; n<box_ene.nComp(); ++n) { 
         box_ene.divide(box_Ja,0,n,1);
      }

   }
   m_energy.exchange(); 
     
}

void PicSpecies::setChargeDensity()
{
   CH_TIME("PicSpecies::setChargeDensity()");
    
   CH_assert(m_charge != 0);
      
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
 
      this_rho.mult(m_charge/m_volume_scale); 
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
         this_rho.mult(m_charge/m_volume_scale); 
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
      this_rho.mult(m_charge/m_volume_scale); 
      
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
                                    m_interpJToGrid );

   }

   // contribute inflow/outflow current to J if being called from explicit solver.
   // This is needed for charge conserving deposits. Note that J from inflow/outflow 
   // particles for iterative implicit solvers is handled differently.
   if(a_from_explicit_solver && m_deposit_bdry_J) {
      m_species_bc->depositInflowOutflowJ( m_currentDensity, m_currentDensity_virtual,
                                          *m_meshInterp, m_interpJToGrid );
   }
   
   // multiply by charge/volume
   for(dit.begin(); dit.ok(); ++dit) {
#if CH_SPACEDIM<3
      m_currentDensity_virtual[dit].getFab().mult(m_charge/m_volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_currentDensity[dit][dir].mult(m_charge/m_volume_scale);
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
   for(dit.begin(); dit.ok(); ++dit) {
#if CH_SPACEDIM<3
      m_inflowJ_virtual[dit].getFab().mult(m_charge/m_volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_inflowJ[dit][dir].mult(m_charge/m_volume_scale);
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
                                         LevelData<EdgeDataBox>&    a_Jtilde,
                                         LevelData<NodeFArrayBox>&  a_Jtildev,
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
                                         a_Jtilde[dit][0],
                                         a_Jtildev[dit].getFab(),
                                         a_Jtildev[dit].getFab(),
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
                                         a_Jtilde[dit][0],
                                         a_Jtilde[dit][1],
                                         a_Jtildev[dit].getFab(),
                                         Bfield_inPlane[0],
                                         Bfield_inPlane[1],
                                         Bfield_virtual,
#endif
                                         m_charge,
                                         alphas,
                                         pList,
                                         m_interpJToGrid );

   }
  
}

void PicSpecies::interpolateEfieldToParticles( const ElectroMagneticFields&  a_em_fields )
{
  CH_TIME("PicSpecies::interpolateEfieldToParticles()");
   
  if(!m_forces) return;
  
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
  CH_TIME("PicSpecies::interpolateFieldsToParticles()");
   
  if(!m_forces) return;
   
  const DisjointBoxLayout& grids(m_mesh.getDBL());
   
  const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
  const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();

  if(SpaceDim==3) {
     
    for(DataIterator dit(grids); dit.ok(); ++dit) {

      const EdgeDataBox& Efield_inPlane = Efield[dit];
      const FluxBox& Bfield_inPlane = Bfield[dit];
         
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
   
      interpolateFieldsToParticles( pList, 
                                    Efield_inPlane, Bfield_inPlane ); 

    }

  }
  else {

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

  }

}

void PicSpecies::interpolateFieldsToOutcasts( const ElectroMagneticFields&  a_em_fields )
{
   CH_TIME("PicSpecies::interpolateFieldsToOutcasts()");
   
   if(!m_forces) return;
   if(m_data.isClosed()) return;  
 
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   // CAUTION!!! This only works when each processor has only 1 box
   // I should create a box level outcast list that lives in ParticleData     
   
   List<JustinsParticle>& pList = m_data.outcast();
   
   const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
   const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();

#if CH_SPACEDIM == 3
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      const EdgeDataBox& Efield_inPlane = Efield[dit];
      const FluxBox& Bfield_inPlane = Bfield[dit];
      interpolateFieldsToParticles( pList, 
                                    Efield_inPlane, Bfield_inPlane );
   }
#else
   const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
   const LevelData<FArrayBox>& Bfield_virt = a_em_fields.getVirtualMagneticField();

   for(DataIterator dit(grids); dit.ok(); ++dit) {

      const EdgeDataBox& Efield_inPlane = Efield[dit];
      const FluxBox& Bfield_inPlane = Bfield[dit];
      const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
      const FArrayBox& Bfield_virtual = Bfield_virt[dit];
      
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
   CH_TIME("PicSpecies::interpolateFieldsToParticles() from particle list");

//#define NEW_EM_INTERP_METHOD
#ifdef NEW_EM_INTERP_METHOD   
   m_meshInterp->interpolateEMfieldsToPart_testing( a_pList,
#else
   m_meshInterp->interpolateEMfieldsToPart( a_pList,
#endif
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

void PicSpecies::interpolateFieldsToParticle( JustinsParticle&  a_particle,
                                        const EdgeDataBox&      a_Efield_inPlane,
                                        const FluxBox&          a_Bfield_inPlane,
                                        const FArrayBox&        a_Efield_virtual,
                                        const FArrayBox&        a_Bfield_virtual )
{
   CH_TIME("PicSpecies::interpolateFieldsToParticle() single particle");

   m_meshInterp->interpolateEMfieldsToPart( a_particle,
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
   CH_TIME("PicSpecies::applyForcesToInflowParticles()");
   if(a_pList.length()==0) return;   
   
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
   
   // advance the velocities to vpbar
   Real cnormDt0, cnormDt1;
   List<JustinsParticle> temp_pList;

   const int pList_length0 = a_pList.length();
   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok();) {

      // compute time step for sub-orbit
      const RealVect& xpold = li().position_old();
      const std::array<Real,3>& vpold = li().velocity_old(); // beta
      cnormDt0 = (X0-xpold[a_bdry_dir])/vpold[a_bdry_dir];
      cnormDt1 = a_cnormDt - cnormDt0;
      CH_assert(cnormDt1>0.0);
         
      // put this particle in a list all by itself
      List<JustinsParticle> single_pList;
      single_pList.transfer(li);
   
      // Boris push this particle
      applyForces(single_pList, cnormDt1, true);
      
      // check for reflection and move this particle to the temp list
      ListIterator<JustinsParticle> li_single(single_pList);
      for(li_single.begin(); li_single.ok();) {
         std::array<Real,3>& vpbar = li_single().velocity();
         int sign_fact = 1 - 2*a_bdry_side;
         if(sign_fact*vpbar[a_bdry_dir]<=0.0) {
            for(int n=0; n<3; n++) vpbar[n] = vpold[n];
            vpbar[a_bdry_dir] = 0.0;
         }
         temp_pList.transfer(li_single);
      }

      //std::array<Real,3>& vp = li().velocity();
      //const std::array<Real,3>& Ep = li().electric_field();
      //Real alpha = m_fnorm_const*cnormDt1/2.0;
      //vp[0] = vpold[0] + alpha*Ep[0];
      //vp[1] = vpold[1];
      //vp[2] = vpold[2];

   }
   
   // transfer particles from the temp list back to the passed list 
   ListIterator<JustinsParticle> li_temp(temp_pList);
   for(li_temp.begin(); li_temp.ok();) a_pList.transfer(li_temp);
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
   CH_TIME("PicSpecies::interpolateFieldsToParticle() single particle");

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
   // Only implemented for 1D so far. 
   // What about external fields ?
   // What about relativisitic ?
   // What about 2D ?
   
   CH_assert(SpaceDim==1);
   CH_assert(m_interpEToParts!=TSC); // derivative calc for TSC not implemented yet
      
   const int n_max = 20; 
   int not_converged_count = 0;   

   std::array<Real,n_max+1> yp; 
   Real fp, dvpdy, dfpdy;
   RealVect dExdy, dEydy, dEzdy;
   RealVect dBxdy, dBydy, dBzdy;
   
   const Real alphas = m_fnorm_const*a_cnormDt/2.0;
   const Real alphasSq = alphas*alphas;
   const Real alphasCu = alphas*alphasSq;
   
   // get boundary info neeed to properly treat particles at 
   // inflow/outflow boundaries
   const IntVect& is_outflow_bc_lo = m_species_bc->getIsOutflowBC_lo(); 
   const IntVect& is_outflow_bc_hi = m_species_bc->getIsOutflowBC_hi(); 
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect& dX(m_mesh.getdX());
 
   Real Bpsq, gammaB;
   Real BpXvpn_x, BpXEp_x, Bpdotvpn, BpdotEp;
   Real dBpXvpn_x, dBpdotvpn, dBpXEp_x, dBpdotEp, dBsq; 

   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {
      
      const RealVect& xpold = li().position_old();
      const std::array<Real,3>& vpold = li().velocity_old();
      
      // references to particle data that will be updated in this loop
      RealVect& xpbar = li().position();
      std::array<Real,3>& vpbar = li().velocity();
      std::array<Real,3>& Ep = li().electric_field();
      std::array<Real,3>& Bp = li().magnetic_field();

      RealVect xpnew;
      
      //  use current value of vpbar to set initial value for y = xpnew-xpold 
      for(int dir=0; dir<SpaceDim; dir++) {
         yp[0] = a_cnormDt*vpbar[dir];
      }
       
      // use newton method to get self-consistent vpbar and xpbar
      int n_last = 0, ng_last = 0;
      Real dyp, rel_error;
      Real yp0 = yp[0];

      const int n_guess = 30;
      Real guess_factor = 1.0;
      //Real damping_factor = 1.0;
      bool converged = false;
      bool is_bdry_part = false;
      for (int ng=0; ng<n_guess; ng++) {
      
         guess_factor = 1.0 + 0.1*ng;
         if(ng>19) guess_factor = 1.0 + 0.01*(ng-19);
               
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
               RealVect& xpold_2 = li().position_old();
               Real dxpold = (xpold[0]-Xmin[0])/pow((double)ng,2);
               xpold_2[0] = Xmin[0] + dxpold; // updates xpold[0]
               yp0 = max(Xmin[0] - xpold[0],yp[n_max-1]);
            }
         }
         
         yp[0] = yp0*guess_factor;

         for (int n=0; n<n_max; n++) { // Newton iterations
              
            // this is a hack, but it seems to work much better 
            // than using a damping factor 
            if(is_bdry_part && xpold[0]<Xmin[0]) {
               if(ng==1) {
               RealVect& xpold_2 = li().position_old();
               xpold_2[0] = xpnew[0] - vpbar[0]*a_cnormDt;
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
            interpolateFieldsToParticle( li(), 
                                         a_Efield_inPlane, a_Bfield_inPlane, 
                                         a_Efield_virtual, a_Bfield_virtual ); 
      
            // compute vpbar[0]
            Bpsq = Bp[0]*Bp[0] + Bp[1]*Bp[1] + Bp[2]*Bp[2];
            gammaB = 1.0 + alphasSq*Bpsq;      
            Bpdotvpn = Bp[0]*vpold[0] + Bp[1]*vpold[1] + Bp[2]*vpold[2];
            BpdotEp = Bp[0]*Ep[0] + Bp[1]*Ep[1] + Bp[2]*Ep[2];
            BpXvpn_x = Bp[1]*vpold[2] - Bp[2]*vpold[1];
            BpXEp_x = Bp[1]*Ep[2] - Bp[2]*Ep[1];
  
            vpbar[0] = ( vpold[0] + alphas*Ep[0] 
                       - alphas*BpXvpn_x  + alphasSq*Bpdotvpn*Bp[0]
                       - alphasSq*BpXEp_x + alphasCu*BpdotEp*Bp[0] )/gammaB;
  
            // update dEp/dy and dBp/dy using current value of xpbar     
            getFieldDerivativesAtParticle( dExdy, dEydy, dEzdy, dBxdy, dBydy, dBzdy,
                                           li(), a_Efield_inPlane, a_Bfield_inPlane, 
                                           a_Efield_virtual, a_Bfield_virtual ); 

            // compute dvpbar[0]/dy
            dBpdotvpn = dBxdy[0]*vpold[0] + dBydy[0]*vpold[1] + dBzdy[0]*vpold[2];
            dBsq = 2.0*(dBxdy[0]*Bp[0] + dBydy[0]*Bp[1] + dBzdy[0]*Bp[2]);
            dBpdotEp  = dBxdy[0]*Ep[0] + dBydy[0]*Ep[1] + dBzdy[0]*Ep[2] 
                      + Bp[0]*dExdy[0] + Bp[1]*dEydy[0] + Bp[2]*dEzdy[0];
            dBpXvpn_x = dBydy[0]*vpold[2] - dBzdy[0]*vpold[1];
            dBpXEp_x  = dBydy[0]*Ep[2] + Bp[1]*dEzdy[0] - dBzdy[0]*Ep[1] - Bp[2]*dEydy[0];

            dvpdy = ( alphas*dExdy[0]
                    - alphas*dBpXvpn_x  + alphasSq*(dBpdotvpn*Bp[0] + Bpdotvpn*dBxdy[0])
                    - alphasSq*dBpXEp_x + alphasCu*(dBpdotEp*Bp[0] + BpdotEp*dBxdy[0])
                    - vpbar[0]*alphasSq*dBsq )/gammaB;
            //dvpdy = alphas*dExdy[0]/gammaB;

            // compute fp and dfp/dy
            fp = yp[n] - a_cnormDt*vpbar[0];
            dfpdy = 1.0 - a_cnormDt*dvpdy;

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
 
      if(n_last>10) cout << "particle newton iter = " << n_last << endl;
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
         cout << "cnormDt = " << a_cnormDt << endl;
         cout << "alphas = " << alphas << endl;
         cout << "yp0 = " << yp0 << ";" << endl;
         cout << "xpold = " << xpold[0] << ";" << endl;
         cout << "vpold_x = " << vpold[0] << ";" << endl;
         cout << "vpold_y = " << vpold[1] << ";" << endl;
         cout << "vpold_z = " << vpold[2] << ";" << endl;
         cout << "Ex0 = " << Ex0 << ";" << endl;
         cout << "Ex1 = " << Ex1 << ";" << endl;
         cout << "Ez0 = " << Ez0 << ";" << endl;
         cout << "Ez1 = " << Ez1 << ";" << endl;
         cout << "Ez2 = " << Ez2 << ";" << endl;
         cout << "By0 = " << By0 << ";" << endl;
         cout << "By1 = " << By1 << ";" << endl;
         cout << std::setprecision(3) << std::scientific << endl;

      }
      
   } // end loop over particle list
   
   // only in-plane vpbar is set above. apply Forces
   // to update all components of vpbar
   applyForces(a_pList, a_cnormDt, true);
   
   // calling advancePositions last will enforce charge conservation
   // for non-converged particles ? .. mass matrix...         
   //advancePositions( a_pList, a_cnormDt/2.0 ); 

   int num_converged = a_pList.length()-not_converged_count;
   return num_converged;
 
}

void PicSpecies::addExternalFieldsToParticles( const ElectroMagneticFields&  a_em_fields )
{
   CH_TIME("PicSpecies::addExternalFieldsToParticles()");
   
   if(!m_forces) return;
   if(!a_em_fields.externalFields()) return;
   
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
   CH_TIME("PicSpecies::addExternalFieldsToParticles() from particle list");
   
   if(!m_forces) return;
   if(!a_em_fields.externalFields()) return;
   
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

void PicSpecies::numberDensity( LevelData<FArrayBox>&  a_rho )
{
   CH_TIME("PicSpecies::numberDensity()");
 
   setNumberDensity();
   const DisjointBoxLayout& grids = a_rho.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      a_rho[dit].copy(m_density[dit]);   
   }
  
}

void PicSpecies::momentumDensity( LevelData<FArrayBox>&  a_mom )
{
   CH_TIME("PicSpecies::momentumDensity()");
 
   setMomentumDensity();
   const DisjointBoxLayout& grids = a_mom.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      a_mom[dit].copy(m_momentum[dit]);   
   }
  
}

void PicSpecies::energyDensity( LevelData<FArrayBox>&  a_ene )
{
   CH_TIME("PicSpecies::energyDensity()");
 
   setEnergyDensity();
   const DisjointBoxLayout& grids = a_ene.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      a_ene[dit].copy(m_energy[dit]);   
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

      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         const Real& wp = li().weight();
         mass_local += wp;         
         
         const std::array<Real,3>& vp = li().velocity(); // actually beta
         momX_local += wp*vp[0];         
         momY_local += wp*vp[1];         
         momZ_local += wp*vp[2];         

         energy_local += wp*(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);

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

bool PicSpecies::isSpecies( const string&  a_name ) const
{
   if(name() == a_name) return true;
   return false;
}


#include "NamespaceFooter.H"

