#include "PicSpecies.H"
#include <cmath>
#include "Constants.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"
#include "BinFabFactory.H"

#include "CH_HDF5.H"
#include "ParticleIO.H"

#include "MathUtils.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

PicSpecies::PicSpecies( ParmParse&         a_ppspc,
                        const int          a_species,
                        const string&      a_name,
                        const MeshInterp&  a_meshInterp,
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
     m_meshInterp(a_meshInterp)
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
   if(SpaceDim<3) {
      m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect);
   }
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
}

void PicSpecies::repositionInflowParticles()
{
   CH_TIME("PicSpecies::repositionInflowParticles()");
    
   // this function is only called by the iterative implicit solvers for particles that
   // inflow during the iterative half particle advance. Need to remove them from m_data
   // and place them back in the inflow_list in order to treate them appropriately
   // for the next particle advance 

   m_species_bc->repositionInflowParticles( m_data );
   
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

   // advance velocities using Boris algorithm
   Real t0, t1, t2, denom, s0, s1, s2;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;
   Real vpl0, vpl1, vpl2;
   
   Real alpha = m_fnorm_const*a_cnormDt/2.0;
   //if(!procID()) {
      //cout << std::setprecision(10) << std::scientific << endl;
      //cout << "m_fnorm_const = " << m_fnorm_const << endl;
      //cout << "cnormDt = " << a_cnormDt << endl;
   //}

   // iterative method for axisymmetric 
   Real moqr = 0.0; // mp/qp/rp
   Real vth_bar = 0.0;
   Real vth_bar_old, dVth;
   const Real  step_tol = 1.0e-12;
   
   // particle arrays are stored as {X,Z,Y} for anticyclic (i.e. cyl_RZ)
   std::array<int,3> dirp = {0,1,2};
   if(m_mesh.anticyclic()) {
      dirp[1] = 2;
      dirp[2] = 1;
   }

   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {

      std::array<Real,3>& vp = li().velocity(); // actually beta
      const std::array<Real,3>& vpold = li().velocity_old(); // actually beta
      const std::array<Real,3>& Ep = li().electric_field();
      const std::array<Real,3>& Bp = li().magnetic_field();
    
      if(m_mesh.axisymmetric() && !m_use_axisymmetric_boris) {
         RealVect& xp = li().position();
         vth_bar = vp[dirp[1]];
         if(abs(xp[0])>0.0) moqr = 1.0/m_fnorm_const/xp[0];
      }

      // add half acceleration to old velocity
      vm0 = vpold[dirp[0]] + alpha*Ep[dirp[0]];
      vm1 = vpold[dirp[1]] + alpha*Ep[dirp[1]];
      vm2 = vpold[dirp[2]] + alpha*Ep[dirp[2]];
   
      int iter = 0;
      while (iter<m_axisymmetric_iter_max) {

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

         if(!m_mesh.axisymmetric() || m_use_axisymmetric_boris) break;
         vth_bar_old = vth_bar;
         vth_bar = (vp[dirp[1]] + vpold[dirp[1]])/2.0;
         dVth = abs(vth_bar-vth_bar_old);
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
   
   //applyBCs(a_half_step);

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
                                   const Real                    a_dt,
                                   const bool                    a_iter_order_swap )
{
   CH_TIME("PicSpecies::advanceParticles()");
   
   // xbar = xn + dt/2*vbar
   // vbar = vn + dt/2*q/m*(E(xbar) + vbar x B(xbar))
   
   const Real cnormDt = a_dt*m_cvac_norm;
   const Real cnormHalfDt = 0.5*cnormDt;
           
   if (!a_iter_order_swap) {
      advancePositions( a_dt, true );
      applyBCs(true);
   }
   interpolateFieldsToParticles( a_em_fields );
   addExternalFieldsToParticles( a_em_fields ); 
   advanceVelocities( a_dt, true );
   if (a_iter_order_swap) {
      advancePositions( a_dt, true );
      applyBCs(true);
   }

}

void PicSpecies::advanceParticlesIteratively( const ElectroMagneticFields&  a_em_fields,
                                              const Real                    a_dt,
                                              const bool                    a_iter_order_swap,
                                              const int                     a_iter_max )
{
   CH_TIME("PicSpecies::advanceParticlesIteratively()");
   
   // Picard method for coupled half dt advance of particle positions and velocities
   //
   // xbar = xn + dt/2*vbar
   // vbar = vn + dt/2*q/m*(E(xbar) + vbar x B(xbar))
   
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
      
      // update xpbar and vpbar 
      if(m_forces && a_iter_order_swap) {
         interpolateFieldsToParticles( pList, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual ); 
         addExternalFieldsToParticles( pList, a_em_fields ); 
         applyForces(pList, cnormDt, true);
      }
      if(m_motion) advancePositions( pList, cnormHalfDt ); 
      if(m_forces) { // interpolate fields to particles and push by dt/2
         interpolateFieldsToParticles( pList, 
                                       Efield_inPlane, Bfield_inPlane, 
                                       Efield_virtual, Bfield_virtual ); 
         addExternalFieldsToParticles( pList, a_em_fields ); 
         applyForces(pList, cnormDt, true);
      }
       
      // no need to iterative if either motion or forces are false
      if(m_motion && m_forces) {

         // transfer not converged particles to a temp list
         List<JustinsParticle> temp_pList;
         const int pListLengthStart = pList.length();      
         stepNormParticleTransfer(pList,temp_pList,cnormHalfDt,false);
   
         // loop over temp list and move back to main list when converged     
         int iter(1);
         while(temp_pList.length()>0) {
            advancePositions( temp_pList, cnormHalfDt ); 
            interpolateFieldsToParticles( temp_pList, 
                                          Efield_inPlane, Bfield_inPlane, 
                                          Efield_virtual, Bfield_virtual ); 
            addExternalFieldsToParticles( temp_pList, a_em_fields ); 
            applyForces(temp_pList, cnormDt, true);
            
            stepNormParticleTransfer(pList,temp_pList,cnormHalfDt,true);
            if (temp_pList.length()==0) break;
            
            // new method for updating particles that won't converge with Picard
            if ( iter >= a_iter_max && m_interpJToGrid == CC0 ) {
               //newtonParticleUpdate( temp_pList, cnormDt,
               //                      Efield_inPlane, Bfield_inPlane, 
               //                      Efield_virtual, Bfield_virtual );

               // set particle position to old position
               ListIterator<JustinsParticle> li(temp_pList);
               for(li.begin(); li.ok(); ++li) {
                  RealVect& xpbar = li().position();
                  const RealVect& xpold = li().position_old();
                  xpbar = xpold;
               }
            
               // compute vbar using old position
               interpolateFieldsToParticles( temp_pList, 
                                             Efield_inPlane, Bfield_inPlane, 
                                             Efield_virtual, Bfield_virtual ); 
               addExternalFieldsToParticles( temp_pList, a_em_fields ); 
               applyForces(temp_pList, cnormDt, true);
           
               // update position
               cellCrossingPositionUpdate( temp_pList, cnormDt);
               
               // compute vbar using new particle position
               interpolateFieldsToParticles( temp_pList, 
                                             Efield_inPlane, Bfield_inPlane, 
                                             Efield_virtual, Bfield_virtual ); 
               addExternalFieldsToParticles( temp_pList, a_em_fields ); 
               applyForces(temp_pList, cnormDt, true);

            }

            if ( iter >= a_iter_max ) {
               cout << "JRA: Picard for particles: iter = " << iter << endl;
               cout << "temp_pList.length() = " << temp_pList.length() << endl;
               break;
            }
            
            iter += 1;

         }

         // call this before putting back non-converged parts...
         // This will break exact charge-conservation, but solver will fail
         // if I don't do this... hopefully the issue is only with NGP method...
         //
         // Maybe don't do it at all... it doesn't seem to help with getting
         // charge-conservation, so there must be something up in JFNK...?
         // For now, we don't do it so smoke tests with mass_matricies pass. 
         //if(m_motion && a_iter_order_swap) advancePositions( pList, cnormHalfDt ); 
          
         // put back particles that could not converge...     
         ListIterator<JustinsParticle> li(temp_pList);
         for(li.begin(); li.ok();) pList.transfer(li);
         
         // make sure we got all the particles
         const int pListLengthEnd = pList.length();      
         if(pListLengthEnd!=pListLengthStart) exit(EXIT_FAILURE);

      } // end if doing iterations
               
   }
 
   // JRA, need to figure out better way to handle particles crossing
   // outlet boundaries so that I can remove this BC call here.
   applyBCs(true);
 
}

void PicSpecies::stepNormParticleTransfer( List<JustinsParticle>&  a_pList,
                                           List<JustinsParticle>&  a_temp_pList,
                                     const Real                    a_cnormHalfDt, 
                                     const bool                    a_reverse )
{
   CH_TIME("PicSpecies::stepNormParticleTransfer()");

   // need to do an extra check for periodic BCs to get dxp0 correct   
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect Lbox = Xmax-Xmin;

   const RealVect& dX(m_mesh.getdX());
   Real dxp0, dxp, rel_diff_dir;
   const Real rel_tol = 1.0e-12;
         
   if(a_reverse) { // temp_pList ==> pList 
      
      ListIterator<JustinsParticle> li(a_temp_pList);
      for(li.begin(); li.ok();) {
         const RealVect& xp = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& vp = li().velocity();
         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xp[dir]-xpold[dir];
            dxp  = vp[dir]*a_cnormHalfDt;
            rel_diff_dir = abs(dxp0-dxp)/dX[dir];
            rel_diff_max = max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max < rel_tol) a_pList.transfer(li);
         else ++li;
      }

   }
   else { //pList ==> temp_pList
 
      ListIterator<JustinsParticle> li(a_pList);
      for(li.begin(); li.ok();) {
         const RealVect& xp = li().position();
         const RealVect& xpold = li().position_old();
         const std::array<Real,3>& vp = li().velocity();
         Real rel_diff_max = 0.0;
         for(int dir=0; dir<SpaceDim; dir++) {
            dxp0 = xp[dir]-xpold[dir];
            dxp  = vp[dir]*a_cnormHalfDt;
            rel_diff_dir = abs(dxp0-dxp)/dX[dir];
            rel_diff_max = max(rel_diff_max,rel_diff_dir);
         }
         if (rel_diff_max >= rel_tol) a_temp_pList.transfer(li);
         else ++li;
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

void PicSpecies::createInflowParticles( const Real&  a_time, 
                                        const Real&  a_dt  )
{
   CH_TIME("PicSpecies::createInflowParticles()");
   
   if(!m_motion) return;

   const Real cnormDt = a_dt*m_cvac_norm;
   m_species_bc->createInflowParticles(a_time,cnormDt);
    
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

void PicSpecies::setStableDt( const CodeUnits&  a_units )
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
   Real local_stable_dt = 1.0/maxDtinv/a_units.CvacNorm();
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

   if(a_restart_file_name.empty()) {
      initializeFromInputFile(a_units);
   }
   else {
      initializeFromRestartFile(a_time,a_restart_file_name);
   }
   
   int totalParticleCount = m_data.numParticles();
   if(!procID()) {
      cout << "Finished initializing pic species " << m_name  << endl;
      cout << "total particles  = " << totalParticleCount << endl << endl;
   }
   
   // define BCs object for this species
   int verbosity = 1; 
   m_species_bc = new PicSpeciesBC( m_name, m_mass, m_mesh, a_units, verbosity );

}

void PicSpecies::initializeFromInputFile( const CodeUnits&  a_units )
{
   if(!procID()) cout << "Initializing pic species " << m_name  << " from input file..." << endl;
   int verbosity=0;
   
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

void PicSpecies::initializeFromRestartFile( const Real          a_time,
                                            const std::string&  a_restart_file_name )
{
   if(!procID()) cout << "Initializing pic species " << m_name  << " from restart file..." << endl;

#ifdef CH_USE_HDF5
   HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );

   // read in the particle data
   std::stringstream s;
   s << "species" << m_species; 
   if(!procID()) cout << s.str() << endl;
   const std::string group_name = std::string( s.str() + "_data");
   handle.setGroup(group_name);
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
      m_meshInterp.moment( box_rho,
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
      m_meshInterp.moment( box_mom,
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
      m_meshInterp.moment( box_ene,
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
      
      m_meshInterp.deposit( box_list.listItems(), 
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
      
         m_meshInterp.deposit( box_list.listItems(), 
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
      
      m_meshInterp.deposit( box_list.listItems(), 
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
   
   // divide by Jacobian after exchange (corrected Jacobian does not have ghosts)
   const LevelData<NodeFArrayBox>& Jacobian = m_mesh.getCorrectedJnc();  
   for(dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_rho = m_chargeDensity_nodes[dit].getFab();
      this_rho.divide(Jacobian[dit].getFab(),0,0,1);
   }
   
}

void PicSpecies::setCurrentDensity()
{
   CH_TIME("PicSpecies::setCurrentDensity()");
    
   CH_assert(m_charge != 0);
      
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();

   for(dit.begin(); dit.ok(); ++dit) {
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      const List<JustinsParticle>& pList = box_list.listItems();
      
      EdgeDataBox& J_inPlane = m_currentDensity[dit];
      for (int dir=0; dir<SpaceDim; ++dir) J_inPlane[dir].setVal(0.0);
#if CH_SPACEDIM < 3
      FArrayBox& J_virtual = m_currentDensity_virtual[dit].getFab();
      J_virtual.setVal(0.0);
#endif
      m_meshInterp.depositCurrent( J_inPlane[0],
#if CH_SPACEDIM >= 2
                                   J_inPlane[1],
#else
                                   J_virtual,
#endif
#if CH_SPACEDIM == 3
                                   J_inPlane[2],
#else
                                   J_virtual,
#endif
                                   pList,
                                   m_interpJToGrid );
#if CH_SPACEDIM < 3
      J_virtual.mult(m_charge/m_volume_scale);
#endif
      for (int dir=0; dir<SpaceDim; ++dir) {
         J_inPlane[dir].mult(m_charge/m_volume_scale);
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

   if(SpaceDim<3) {
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
   }
   
}

void PicSpecies::accumulateMassMatrices( LevelData<EdgeDataBox>&    a_sigma_xx, 
                                         LevelData<EdgeDataBox>&    a_sigma_xy,
                                         LevelData<EdgeDataBox>&    a_sigma_xz,
                                         LevelData<NodeFArrayBox>&  a_sigma_yx,
                                         LevelData<NodeFArrayBox>&  a_sigma_yy,
                                         LevelData<NodeFArrayBox>&  a_sigma_yz,
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
   CH_assert(SpaceDim==1); // only works for 1D so far...

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
      
      m_meshInterp.depositMassMatrices( a_sigma_xx[dit][0],
                                        a_sigma_xy[dit][0],
                                        a_sigma_xz[dit][0],
                                        a_sigma_yx[dit].getFab(),
                                        a_sigma_yy[dit].getFab(),
                                        a_sigma_yz[dit].getFab(),
                                        a_sigma_zx[dit].getFab(),
                                        a_sigma_zy[dit].getFab(),
                                        a_sigma_zz[dit].getFab(),
                                        a_Jtilde[dit][0],
                                        a_Jtildev[dit].getFab(),
                                        Bfield_inPlane[0],
                                        Bfield_virtual,
                                        m_charge,
                                        alphas,
                                        pList,
                                        m_interpJToGrid );

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

   if(SpaceDim==3) {
     
     for(DataIterator dit(grids); dit.ok(); ++dit) {

       const EdgeDataBox& Efield_inPlane = Efield[dit];
       const FluxBox& Bfield_inPlane = Bfield[dit];
      
       interpolateFieldsToParticles( pList, 
                                     Efield_inPlane, Bfield_inPlane );
     }

   }
   else {
 
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

  }
   
}

void PicSpecies::interpolateFieldsToParticles( List<JustinsParticle>&  a_pList,
                                         const EdgeDataBox&  a_Efield_inPlane,
                                         const FluxBox&      a_Bfield_inPlane )
{
   CH_TIME("PicSpecies::interpolateFieldsToParticles() from particle list");
   
   m_meshInterp.interpolateEMfieldsToPart( a_pList,
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
   m_meshInterp.interpolateEMfieldsToPart_testing( a_pList,
#else
   m_meshInterp.interpolateEMfieldsToPart( a_pList,
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

void PicSpecies::newtonParticleUpdate( List<JustinsParticle>&  a_pList,
                                 const Real          a_cnormDt, 
                                 const EdgeDataBox&  a_Efield_inPlane,
                                 const FluxBox&      a_Bfield_inPlane,
                                 const FArrayBox&    a_Efield_virtual,
                                 const FArrayBox&    a_Bfield_virtual )
{
   CH_TIME("PicSpecies::newtonParticleUpdate()");

   // Sometimes Picard iteration method fails for particles crossing cells
   // when using the charge-conserving interpolation method. Need to do 
   // Newton method to get self-consistent particle position and velocity. 
   //
   // Only implemented for 1D electrostatic so far. 
   // What about external fields ?
   // What about vpxBp ?
   // What about relativisitic ?
   // What about 2D ?
   
   CH_assert(SpaceDim==1);
   
   RealVect xpnew, Xcell;
   IntVect index_old, index_new;
   Real E0, E1;
       
   const RealVect& dX(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   
   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {
      RealVect& xp = li().position();
      const RealVect& xpold = li().position_old();
      std::array<Real,3>& vpbar = li().velocity();
      const std::array<Real,3>& vpold = li().velocity_old();
      std::array<Real,3>& Ep = li().electric_field();

      //  use old position to compute Ep -> vpbar -> xpnew 
      for(int dir=0; dir<SpaceDim; dir++) {
         index_old[dir] = floor((xpold[dir] - meshOrigin[dir])/dX[dir]);
      }
      E0 = a_Efield_inPlane[0](index_old,0);
      
      Ep[0] = E0;
      vpbar[0] = vpold[0] + a_cnormDt/2.0*m_fnorm_const*Ep[0];
      xpnew[0] = xpold[0] + a_cnormDt*vpbar[0];
      
      for(int dir=0; dir<SpaceDim; dir++) {
         index_new[dir] = floor((xpnew[dir] - meshOrigin[dir])/dX[dir]);
      }
      E1 = a_Efield_inPlane[0](index_new,0);
     
      // set the value for X between xpold and xpnew
      for(int dir=0; dir<SpaceDim; dir++) {
         int index_max = max(index_old[dir],index_new[dir]);
         Xcell[dir] = meshOrigin[dir] + dX[dir]*index_max;
      }
  
      bool print_stuff = true;
      if(print_stuff) {
         cout << std::setprecision(10) << std::scientific << endl;
         cout << "JRA: xpold = " << xpold[0] << endl;
         cout << "     Xcell = " << Xcell[0] << endl;
         cout << "     xpnew = " << xpnew[0] << endl;
         cout << "     vpold = " << vpold[0] << endl;
         cout << "     vpbar = " << vpbar[0] << endl;
         cout << "     E0 = " << E0 << endl;
         cout << "     E1 = " << E1 << endl;
         cout << "     index_old = " << index_old[0] << endl;
         cout << "     index_new = " << index_new[0] << endl;
         cout << std::setprecision(2) << std::scientific << endl;
      }

      if(index_new[0]==index_old[0]) { // shouldn't happen, but it does...
         // I think its because I update position first using old vbar.
         // If I updated vbar first, then this should never occur
         xp[0] = (xpnew[0] + xpold[0])/2.0;
         continue;
      }

      // use Newton method to get xpnew, vpbar, and Ep
      xpnew[0] = Xcell[0];
      Real dXp = xpnew[0] - xpold[0];
      Real fun, Jac, alpha, dX_seg0, dX_seg1;
      
      dX_seg0 = Xcell[0] - xpold[0];
      alpha = a_cnormDt*a_cnormDt/2.0*m_fnorm_const*dX_seg0*(E0 - E1);
      
      int n_max = 30; 
      int n_last = 0;
      Real dXp0 = dXp;
      Real rel_tol = 1.0e-12;
      Real rel_error = 1.0;
      for (int n=0; n<n_max; n++) {

         dX_seg1 = xpnew[0] - Xcell[0];
         Ep[0] = dX_seg0/dXp*E0 + dX_seg1/dXp*E1;

         vpbar[0] = vpold[0] + a_cnormDt/2.0*m_fnorm_const*Ep[0];
         
         fun = dXp - a_cnormDt*vpbar[0];
         Jac = 1.0 + alpha/dXp/dXp;

         dXp = dXp - fun/Jac;         
         xpnew[0] = xpold[0] + dXp;
    
         rel_error = abs(dXp-dXp0)/dX[0];
         cout << " rel_error = " << rel_error << endl;
         n_last = n;
         if(rel_error<rel_tol) break;
         dXp0 = dXp;
      }
      cout << " n = " << n_last << endl;
      if(n_last==n_max-1) {
         cout << " newton not converged " << endl;
      }
      dX_seg1 = xpnew[0] - Xcell[0];
      Ep[0] = dX_seg0/dXp*E0 + dX_seg1/dXp*E1;
      vpbar[0] = vpold[0] + a_cnormDt/2.0*m_fnorm_const*Ep[0];
      xp[0] = (xpnew[0] + xpold[0])/2.0;
      //xp[0] = xpold[0] + a_cnormDt/2.0*vpbar[0];
      
      if(n_last>=2) { // JRA, hack for testing purposes ...
                      // It seems that all particles that give issues are those
                      // that live just on one side of a cell and just barely
                      // cross over into the other cell....
         Ep[0] = E0;
         vpbar[0] = vpold[0] + a_cnormDt/2.0*m_fnorm_const*Ep[0];
         xpnew[0] = Xcell[0];
         xp[0] = (xpnew[0] + xpold[0])/2.0;
      }

   }

}

void PicSpecies::cellCrossingPositionUpdate( List<JustinsParticle>&  a_pList,
                                       const Real          a_cnormDt ) 
{
   CH_TIME("PicSpecies::cellCrossingPositionUpdate()");

   // Sometimes Picard iteration method fails for particles crossing cells
   // when using the charge-conserving interpolation method. This always seems
   // to occur when a particle starts just on one side of a cell face and cross
   // just to the other. See Chen and Chacon 2015 paper below Eq. 43.  
   //
   // In this situation, we force the new position to be right on the cell face.
   //
   
   RealVect xpnew, Xcell;
   IntVect index_old, index_new;
   IntVect cell_crossings;
      
   Real epsilon = 1.0e-12;
   Real rel_tol = 1.0e-4;
       
   const RealVect& dX(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   
   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {

      // set xpnew using vpbar computed with old position
      RealVect& xpbar = li().position();
      const RealVect& xpold = li().position_old();
      const std::array<Real,3>& vpbar = li().velocity();      
      for(int dir=0; dir<SpaceDim; dir++) {
         xpnew[dir] = xpold[dir] + a_cnormDt*vpbar[dir];
      }

      // compute index for old and new particle positions
      for(int dir=0; dir<SpaceDim; dir++) {
         index_old[dir] = floor((xpold[dir] - meshOrigin[dir])/dX[dir]);
         index_new[dir] = floor((xpnew[dir] - meshOrigin[dir])/dX[dir]);
      }
      cell_crossings = index_new - index_old;
     
      // set the value for Xcell between xpold and xpnew
      for(int dir=0; dir<SpaceDim; dir++) {
         int shift = 0;
         if(cell_crossings[dir]<0) shift = 1;
         Xcell[dir] = meshOrigin[dir] + dX[dir]*(index_new[dir]+shift);  
      }

      // need to force position to be on cell face in dir
      // (don't know how to do this reliably for SpaceDim>1)
      // How to determine which dir cell crossing is an issue?
      // Using dXp = xpnew-Xcell with tolerance may not work well 
      // if not using Newton method to determing true particle orbit....
      
      // For now, just use position based on forces from xpold position
      // to update position (xpnew = xpold before call here)
      Real dXp;
      RealVect sign = IntVect::Unit;
      for(int dir=0; dir<SpaceDim; dir++) {
         dXp = xpnew[dir]-Xcell[dir];
         if(cell_crossings[dir]!=0 && abs(dXp)<rel_tol*dX[dir]) {
            if(cell_crossings[dir]<0) sign[dir] = -1;
            xpnew[dir] = Xcell[dir] + dX[dir]*epsilon*sign[dir];
         }
         xpbar[dir] = (xpnew[dir] + xpold[dir])/2.0;
      }

   }

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

