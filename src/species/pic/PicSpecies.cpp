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
#include "ParticleUtils.H"

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
     m_interpToGrid(CIC),
     m_species_bc(NULL),
     m_motion(true),
     m_forces(true),
     m_scatter(false),
     m_write_all_part_comps(false),
     m_mesh(a_mesh),
     m_meshInterp(a_meshInterp)
{
   m_name = a_name;
   
   string interp_type_parts;
   a_ppspc.query( "interp_type_parts", interp_type_parts );
   //a_ppspc.query( "interp_type_fields", interp_type_fields );
   if(interp_type_parts=="CIC") {
      m_interpToGrid = CIC;
   }
   CH_assert(m_interpToGrid==CIC);
 
   a_ppspc.get( "mass", m_mass );
   a_ppspc.get( "charge", m_charge );
   a_ppspc.query( "Uint", m_Uint );
   a_ppspc.query( "motion", m_motion );
   a_ppspc.query( "forces", m_forces );
   a_ppspc.query( "scatter", m_scatter );
   a_ppspc.query( "write_all_comps", m_write_all_part_comps );

   if ( procID() == 0 ) {
      cout << "  name = " << m_name << endl;
      cout << "  number = " << m_species << endl;
      cout << "  mass/me   = " << m_mass << endl;
      cout << "  charge/qe = " << m_charge << endl;
      cout << "  Uint [eV] = " << m_Uint << endl;
      cout << "  interp to grid type = " << interp_type_parts << endl;
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
                                                   const Real&                   a_cnormDt, 
                                                   const bool&                   a_byHalfDt )
{
   CH_TIME("PicSpecies::repositionOutcastsAndApplyForces()");
    
   // this function is only called by the iterative implicit solvers for particles that cross
   // a physical boundary with the outflow bc during the iterations. Some may end up being 
   // outcasts in the end, but to be sure the particles are placed back on the boundary and
   // their velocity is updated prior to the next position advance. For each iteration,
   // if they are out of the boundary after the position advance then they do not contribute
   // to the current

   // reposition the outcasts that are out of bounds to be just inside the domain boundaries
   List<JustinsParticle>& outcast_pList = m_data.outcast();
   m_species_bc->repositionOutflowParticles(outcast_pList);

   // interpolate the fields to the particles 
   interpolateFieldsToOutcasts( a_em_fields );

   // advance velocities using Boris algorithm
   applyForces(outcast_pList, a_cnormDt, a_byHalfDt);     

   // finally, remap the outcasts
   m_data.remapOutcast();
   CH_assert(m_data.isClosed()); // make sure all particles are accounted for
   
}

void PicSpecies::applyForces( List<JustinsParticle>&  a_pList,
                        const Real&  a_cnormDt, 
                        const bool&  a_byHalfDt )
{
   CH_TIME("PicSpecies::applyForces()");
   
   // advance velocities using Boris algorithm
   Real t0, t1, t2, denom, s0, s1, s2;
   Real vm0, vm1, vm2, vpr0, vpr1, vpr2;
   Real vpl0, vpl1, vpl2;
   
   Real alpha = m_fnorm_const*a_cnormDt/2.0;
 
   ListIterator<JustinsParticle> li(a_pList);
   for(li.begin(); li.ok(); ++li) {

      std::array<Real,3>& vp = li().velocity(); // actually beta
      const std::array<Real,3>& vpold = li().velocity_old(); // actually beta
      const std::array<Real,3>& Ep = li().electric_field();
      const std::array<Real,3>& Bp = li().magnetic_field();
    
      //ParticleUtils::borisPusher(vp,vpold,Ep,Bp,m_fnorm_const,a_cnormDt);
   
      t0 = alpha*Bp[0];
      t1 = alpha*Bp[1];
      t2 = alpha*Bp[2];
      denom = 1.0 + t0*t0 + t1*t1 + t2*t2;
      s0 = 2.0*t0/denom;
      s1 = 2.0*t1/denom;
      s2 = 2.0*t2/denom;

      // add half acceleration to old velocity
      vm0 = vpold[0] + alpha*Ep[0];
      vm1 = vpold[1] + alpha*Ep[1];
      vm2 = vpold[2] + alpha*Ep[2];

      // define vpr = vm + vm x t
      vpr0 = vm0 + vm1*t2 - vm2*t1;
      vpr1 = vm1 + vm2*t0 - vm0*t2;
      vpr2 = vm2 + vm0*t1 - vm1*t0;

      // rotate (define vplus = vminus + vprime x s)
      vpl0 = vm0 + vpr1*s2 - vpr2*s1;
      vpl1 = vm1 + vpr2*s0 - vpr0*s2;
      vpl2 = vm2 + vpr0*s1 - vpr1*s0;

      // add another half acceleration
      vp[0] = vpl0 + alpha*Ep[0];
      vp[1] = vpl1 + alpha*Ep[1];
      vp[2] = vpl2 + alpha*Ep[2];

      if(a_byHalfDt) {
         vp[0] = (vp[0] + vpold[0])/2.0;
         vp[1] = (vp[1] + vpold[1])/2.0;
         vp[2] = (vp[2] + vpold[2])/2.0;
      }

   } // end loop over particle list

}

void PicSpecies::advancePositions( const Real&  a_cnormDt,
                                   const bool&  a_intermediate_advance )
{
   CH_TIME("PicSpecies::advancePositions()");
    
   // Each proc loops over its own boxes
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      // loop over particles in this box and advance
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         RealVect& xp = li().position();
         const RealVect& xpold = li().position_old();
         std::array<Real,3>& vp = li().velocity(); // actually beta

         // update particle position
         for(int dir=0; dir<SpaceDim; dir++) {
            xp[dir] = xpold[dir] + vp[dir]*a_cnormDt;
         }

      } // end loop over particle list
      
   } // end loop over boxes

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

void PicSpecies::advancePositions_2ndHalf()
{
   CH_TIME("PicSpecies::advancePositions_2ndHalf()");
    
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

   // gather outcast particles
   m_data.gatherOutcast();

   // apply BCs to outcasts that are also out of bounds
   List<JustinsParticle>& outcast_list = m_data.outcast();
   Real dummy_time = 0.0;
   bool intermediate_advance = false;
   m_species_bc->apply(outcast_list,intermediate_advance,dummy_time);
   
   // remap the outcasts
   m_data.remapOutcast();
   CH_assert(m_data.isClosed()); // make sure all particles are accounted for
   
}

void PicSpecies::advanceVelocities( const Real&  a_cnormDt, 
                                    const bool&  a_byHalfDt )
{
   CH_TIME("PicSpecies::advanceVelocities()");
              
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();

      applyForces(pList, a_cnormDt, a_byHalfDt);

   }

}

void PicSpecies::advanceVelocities_2ndHalf()
{
   CH_TIME("PicSpecies::advanceVelocities_2ndHalf()");
   
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

void PicSpecies::removeOutflowParticles()
{
   CH_TIME("PicSpecies::removeOutflowParticles()");
    
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
   
   m_species_bc->createInflowParticles(a_time,a_dt);
    
}

void PicSpecies::updateOldParticlePositions()
{
   CH_TIME("PicSpecies::updateOldParticlePositions()");
    
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

   // update stable time step
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

   int IC_count = 1;
   bool more_ICs = true;
   while(more_ICs) {
      
      std::vector<int> partsPerCellstd(SpaceDim);
      ppspcIC.getarr("parts_per_cell",partsPerCellstd,0,SpaceDim);
      IntVect partsPerCell; // convert std::vector<int> to IntVect
      for (int dir=0; dir<SpaceDim; ++dir) {
         partsPerCell[dir] = 0;
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
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppdenIC,verbosity);
      gridFunction->assign( m_density, m_mesh, this_time );

      // set temperature profiles from ICs
      const DisjointBoxLayout& grids(m_mesh.getDBL());
      LevelData<FArrayBox> tempProfile;
      tempProfile.define(grids,1,m_temperature.ghostVect());

      const std::string spctemp0IC(spcIC + ".temperature_0");
      ParmParse pptemp0IC( spctemp0IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionTemp0 = gridFactory.create(pptemp0IC,verbosity);
      gridFunctionTemp0->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_temperature[dit].copy(tempProfile[dit],0,0,1);
      }
   
      const std::string spctemp1IC(spcIC + ".temperature_1");
      ParmParse pptemp1IC( spctemp1IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionTemp1 = gridFactory.create(pptemp1IC,verbosity);
      gridFunctionTemp1->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_temperature[dit].copy(tempProfile[dit],0,1,1);
      }
   
      const std::string spctemp2IC(spcIC + ".temperature_2");
      ParmParse pptemp2IC( spctemp2IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionTemp2 = gridFactory.create(pptemp2IC,verbosity);
      gridFunctionTemp2->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_temperature[dit].copy(tempProfile[dit],0,2,1);
      }
   
      // set mean velocity profiles from ICs
      const std::string spcvel0IC(spcIC + ".velocity_0");
      ParmParse ppvel0IC( spcvel0IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionVel0 = gridFactory.create(ppvel0IC,verbosity);
      gridFunctionVel0->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_velocity[dit].copy(tempProfile[dit],0,0,1);
      }
   
      const std::string spcvel1IC(spcIC + ".velocity_1");
      ParmParse ppvel1IC( spcvel1IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionVel1 = gridFactory.create(ppvel1IC,verbosity);
      gridFunctionVel1->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_velocity[dit].copy(tempProfile[dit],0,1,1);
      }
   
      const std::string spcvel2IC(spcIC + ".velocity_2");
      ParmParse ppvel2IC( spcvel2IC.c_str() );
      RefCountedPtr<GridFunction> gridFunctionVel2 = gridFactory.create(ppvel2IC,verbosity);
      gridFunctionVel2->assign( tempProfile, m_mesh, this_time );
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         m_velocity[dit].copy(tempProfile[dit],0,2,1);
      }

   ////////////////////////////////////////////////////////////////////////////////////////

      int totalPartsPerCell = 1;
      Real cellVolume = m_volume_scale;
      for (int dir=0; dir<SpaceDim; dir++)
      {
         totalPartsPerCell *= partsPerCell[dir];
         cellVolume = cellVolume*dX[dir];
      }
      Real pWeight = 0.0; 
   
      // create sub-box and dX for particles
      const Box partSubBox(IntVect::Zero, partsPerCell-IntVect::Unit);
      BoxIterator pbit(partSubBox);
   
      RealVect dXpart = dX;
      if(totalPartsPerCell > 1 ) {
         for(int dir=0; dir<SpaceDim; dir++) dXpart[dir] /= partsPerCell[dir];
      }

      // loop over boxes and set the initial particle values (pos., vel., weight)
      const BoxLayout& BL = m_data.getBoxes();
      DataIterator dit(BL);
      uint64_t ID = procID()*512 + 1; // hack for testing purposes
      //Real ID = 0.;
      Real numDen_scale = a_units.getScale(a_units.NUMBER_DENSITY); // multiply input by this to get density in SI
      for (dit.begin(); dit.ok(); ++dit) {

         CH_XD::List<JustinsParticle> thisList;
     
         // loop over grid indices
         const Box gridBox = BL.get(dit);
         BoxIterator gbit(gridBox);
         for(gbit.begin(); gbit.ok(); ++gbit) {

            const IntVect ig = gbit(); // grid index
            Real local_density = m_density[dit].get(ig,0);
            pWeight = numDen_scale*local_density*cellVolume/(Real)totalPartsPerCell; 
        
            RealVect local_Xcc; 
            std::array<Real,3> local_temperature;
            std::array<Real,3> local_velocity;
            for(int dir=0; dir<SpaceDim; dir++) {
               local_Xcc[dir] = Xcc[dit].get(ig,dir);
            }
            for(int dir=0; dir<3; dir++) {
               local_temperature[dir] = m_temperature[dit].get(ig,dir);
               local_velocity[dir] = m_velocity[dit].get(ig,dir);
            }
     
            // loop over subgrid corresponding to where particles are placed in each grid cell    
            Real V0 = sqrt(Constants::QE/Constants::ME); // ele thermal speed at 1eV [m/s]
            for(pbit.begin(); pbit.ok(); ++pbit) {
            
               // set particle position uniformly on grid
               RealVect Xpart = local_Xcc - 0.5*dX;
               IntVect ipg = pbit();
               for(int dir=0; dir<SpaceDim; dir++) {
                  Xpart[dir] += (ipg[dir] + 0.5)*dXpart[dir];
               }

               // check to see if this particle position is inside specified region
               if(Xpart[0]<sXmin[0] || Xpart[0]>sXmax[0]) continue;
               if(SpaceDim==2) {
                  if(Xpart[1]<sXmin[1] || Xpart[1]>sXmax[1]) continue;
               }
               if(SpaceDim==3) {
                  if(Xpart[2]<sXmin[2] || Xpart[2]>sXmax[2]) continue;
               }

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

         // finally, add particles destructively to this ListBox. Those that are
         // left behind are outcasts.
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
         IC_count = IC_count - 1;
         more_ICs = false;
      }
  
   } // end while loop
   
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
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density[dit];
      box_rho.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = density;
      m_meshInterp.moment( box_rho,
                           box_list.listItems(),
                           m_mass,
                           thisMoment );
 
      box_rho.divide(m_volume_scale);

   }
   m_density.exchange(); // causes ERROR: corrupted double-linked list at code exit!!!!   
     
}

void PicSpecies::setMomentumDensity()
{
   CH_TIME("PicSpecies::setMomentumDensity()");
    
   //CH_assert(m_data.isClosed());
    
   const DisjointBoxLayout& grids = m_momentum.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_mom = m_momentum[dit];
      box_mom.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = momentum;
      m_meshInterp.moment( box_mom,
                           box_list.listItems(),
                           m_mass/m_volume_scale,
                           thisMoment ); 

   }
   m_momentum.exchange(); 
     
}

void PicSpecies::setEnergyDensity()
{
   CH_TIME("PicSpecies::setEnergyDensity()");
    
   //CH_assert(m_data.isClosed());
    
   const DisjointBoxLayout& grids = m_energy.disjointBoxLayout();
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

   }
   m_energy.exchange(); 
     
}

void PicSpecies::setChargeDensity()
{
   CH_TIME("PicSpecies::setChargeDensity()");
    
   //CH_assert(m_data.isClosed());
   CH_assert(m_charge != 0);
      
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   DataIterator dit = grids.dataIterator();

   for(dit.begin(); dit.ok(); ++dit) {
      
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      // deposit on cells 
      FArrayBox& this_rho = m_chargeDensity[dit];
      this_rho.setVal(0.0);
      
      m_meshInterp.deposit( box_list.listItems(), 
                            this_rho,
                            m_interpToGrid ); 
      this_rho.mult(m_charge/m_volume_scale); 
    
   }
   
   // add ghost cells to valid cells
   LDaddOp<FArrayBox> addOp;
   m_chargeDensity.exchange(m_chargeDensity.interval(), m_mesh.reverseCopier(), addOp);
   //m_chargeDensity.exchange(); // needed if more than 1 box per proccesor. bug?

}

void PicSpecies::setChargeDensityOnFaces()
{
   CH_TIME("PicSpecies::setChargeDensityOnFaces()");
    
   //CH_assert(m_data.isClosed());
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
                               m_interpToGrid ); 
         this_rho.mult(m_charge/m_volume_scale); 
      }
      
   }
   
   // add ghost cells to valid cells
   LDaddFaceOp<FluxBox> addFaceOp;
   m_chargeDensity_faces.exchange(m_chargeDensity_faces.interval(), m_mesh.reverseCopier(), addFaceOp);
   //SpaceUtils::exchangeFluxBox(m_chargeDensity_faces); // needed if more than 1 box per proccesor. bug?

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
                            m_interpToGrid ); 
      this_rho.mult(m_charge/m_volume_scale); 
      
   }
      
   // add charge density deposited on ghost cells to valid cells
   LDaddNodeOp<NodeFArrayBox> addNodeOp;
   m_chargeDensity_nodes.exchange(m_chargeDensity_nodes.interval(), 
                                  m_mesh.reverseCopier(), addNodeOp);
   
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
      if(SpaceDim==3) {
         m_meshInterp.depositCurrent( J_inPlane[0],
                                      J_inPlane[1],
                                      J_inPlane[2],
                                      pList,
                                      m_interpToGrid );
      }
      else {
         FArrayBox& J_virtual = m_currentDensity_virtual[dit].getFab();
         J_virtual.setVal(0.0);
         if(SpaceDim==2) {
            m_meshInterp.depositCurrent( J_inPlane[0],
                                         J_inPlane[1],
                                         J_virtual,
                                         pList,
                                         m_interpToGrid );
         }
         else {
            m_meshInterp.depositCurrent1D( J_inPlane[0],
                                           J_virtual,
                                           pList,
                                           m_interpToGrid );
         }
         J_virtual.mult(m_charge/m_volume_scale);
      }
      for (int dir=0; dir<SpaceDim; ++dir) J_inPlane[dir].mult(m_charge/m_volume_scale);
   
   }
   
   // add ghost cells to valid cells
   LDaddEdgeOp<EdgeDataBox> addEdgeOp;
   m_currentDensity.exchange(m_currentDensity.interval(), m_mesh.reverseCopier(), addEdgeOp);
   //SpaceUtils::exchangeEdgeDataBox(m_currentDensity); // needed if more than 1 box per proccesor. bug?
   
   if(SpaceDim<3) {
      LDaddNodeOp<NodeFArrayBox> addNodeOp;
      m_currentDensity_virtual.exchange(m_currentDensity_virtual.interval(), 
                                        m_mesh.reverseCopier(), addNodeOp);
      //SpaceUtils::exchangeNodeFArrayBox(m_currentDensity_virtual,m_mesh); // needed if more than 1 box per proccesor. bug?
   }
   
}

void PicSpecies::interpolateFieldsToParticles( const ElectroMagneticFields&  a_em_fields )
{
   CH_TIME("PicSpecies::interpolateFieldsToParticles()");
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   
   const LevelData<EdgeDataBox>& Efield = a_em_fields.getElectricField();
   const LevelData<FluxBox>& Bfield = a_em_fields.getMagneticField();

   if(SpaceDim==3) {
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         const EdgeDataBox& Efield_inPlane = Efield[dit];
         const FluxBox& Bfield_inPlane = Bfield[dit];
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         
         m_meshInterp.interpolateEMfieldsToPart( pList,
                                                 Efield_inPlane[0],
                                                 Efield_inPlane[1],
                                                 Efield_inPlane[2],
                                                 Bfield_inPlane[0],
                                                 Bfield_inPlane[1],
                                                 Bfield_inPlane[2],
                                                 m_interpToGrid );
      }
   
   }
   else {

      const LevelData<NodeFArrayBox>& Efield_virt = a_em_fields.getVirtualElectricField();
      const LevelData<FArrayBox>& Bfield_virt = a_em_fields.getVirtualMagneticField();
      
      for(DataIterator dit(grids); dit.ok(); ++dit) {

         const EdgeDataBox& Efield_inPlane = Efield[dit];
         const FluxBox& Bfield_inPlane = Bfield[dit];
         ListBox<JustinsParticle>& box_list = m_data[dit];
         List<JustinsParticle>& pList = box_list.listItems();
         
         const FArrayBox& Efield_virtual = Efield_virt[dit].getFab();
         const FArrayBox& Bfield_virtual = Bfield_virt[dit];
         if(SpaceDim==2) {
            m_meshInterp.interpolateEMfieldsToPart( pList,
                                                    Efield_inPlane[0],
                                                    Efield_inPlane[1],
                                                    Efield_virtual,
                                                    Bfield_inPlane[0],
                                                    Bfield_inPlane[1],
                                                    Bfield_virtual,
                                                    m_interpToGrid );
         }  
         else {
            m_meshInterp.interpolateEMfieldsToPart1D( pList,
                                                      Efield_inPlane[0],
                                                      Efield_virtual,
                                                      Bfield_inPlane[0],
                                                      Bfield_virtual,
                                                      m_interpToGrid );
         }
      }

   }

}

void PicSpecies::interpolateFieldsToOutcasts( const ElectroMagneticFields&  a_em_fields )
{
   CH_TIME("PicSpecies::interpolateFieldsToOutcasts()");
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
         
         m_meshInterp.interpolateEMfieldsToPart( pList,
                                                 Efield_inPlane[0],
                                                 Efield_inPlane[1],
                                                 Efield_inPlane[2],
                                                 Bfield_inPlane[0],
                                                 Bfield_inPlane[1],
                                                 Bfield_inPlane[2],
                                                 m_interpToGrid );
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
         if(SpaceDim==2) {
            m_meshInterp.interpolateEMfieldsToPart( pList,
                                                    Efield_inPlane[0],
                                                    Efield_inPlane[1],
                                                    Efield_virtual,
                                                    Bfield_inPlane[0],
                                                    Bfield_inPlane[1],
                                                    Bfield_virtual,
                                                    m_interpToGrid );
         }  
         else {
            m_meshInterp.interpolateEMfieldsToPart1D( pList,
                                                      Efield_inPlane[0],
                                                      Efield_virtual,
                                                      Bfield_inPlane[0],
                                                      Bfield_virtual,
                                                      m_interpToGrid );
         }
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

/*
void PicSpecies::createMeshInterp()
{
   CH_TIME("PicSpecies::createMeshInterp()");

   // get some mesh information
   //
   //DomainGrid* mesh = DomainGrid::mesh;
   const ProblemDomain& domain(m_mesh.getDomain()); 
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
  
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
      m_meshInterp = NULL;
   }
   // create/set the meshInterp object pointer
   //m_meshInterp = new MeshInterp( domain.domainBox(),
   //                               meshSpacing,
   //                               meshOrigin );
   m_meshInterp = static_cast<MeshInterp*> (new MeshInterp( domain.domainBox(),
                                                            meshSpacing,
                                                            meshOrigin  ));
   
   //m_meshInterp = static_cast<RefCountedPtr<MeshInterp>> (new MeshInterp( domain.domainBox(),
   //                                                                       meshSpacing,
   //                                                                       meshOrigin  ));
   //m_isMeshInterpSet = true;

}
*/ 

bool PicSpecies::isSpecies( const string&  a_name ) const
{
   if(name() == a_name) return true;
   return false;
}


#include "NamespaceFooter.H"

