
#include "PicSpecies.H"
#include <array>
#include <cmath>
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"
#include "BinFabFactory.H"

#include "MathUtils.H"

#include "NamespaceHeader.H"

PicSpecies::PicSpecies( ParmParse&         a_ppspc,
                        const string&      a_name,
                        const MeshInterp&  a_meshInterp,
                        //MeshInterp*  a_meshInterp,
                        const DomainGrid&  a_mesh )
   : m_mass(1.0),
     m_Uint(0.0),
     m_charge(1.0),
     m_motion(true),
     m_forces(true),
     m_mesh(a_mesh),
     //m_meshInterp(NULL)
     m_meshInterp(a_meshInterp)
{
   // maybe parse some stuff
   m_name = a_name;
 
   //createMeshInterp();
 
   a_ppspc.query( "mass", m_mass );
   a_ppspc.query( "charge", m_charge );

   if ( procID() == 0 ) {
      cout << "name = " << m_name << endl;
      cout << "mass = " << m_mass << endl;
      cout << "charge = " << m_charge << endl;
   }

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   const int ghosts(m_mesh.ghosts());
   //RealVect particleOrigin = RealVect(D_DECL(0.0,0.0,0.0));
   //if(!procID()) cout << "particleOrigin = " << particleOrigin << endl; 

   // set initial piston position
   //
   const RealVect& Xmax(m_mesh.getXmax());
   rpiston = Xmax[0];
   m_stable_dt = meshSpacing[0]/abs(vpiston); // initialize stable time step

   // initialize the member LevelDatas
   //
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   m_density.define(grids,1,ghostVect);
   m_momentum.define(grids,SpaceDim,ghostVect);
   m_energy.define(grids,1,ghostVect);
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
   // I need a better way to ensure this and get this parameter here!
   //
   int fixedBoxSize;
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      Box thisbox( grids[dit] ); 
      fixedBoxSize = thisbox.size(0);
      break;
   } 
   
   // ParticleData<T> behaves similar to LevelData<FArrayBox>
   // ghosts?
   //
   m_data.define(grids, domain, fixedBoxSize,
                 meshSpacing, meshOrigin);
   
   
   // In order to initialize a LevelData<BinFab<T>>, first need to define a
   // BinFabFactory. The factor has a "create" function to define "new" pointers
   // to the BinFab that lives at the box level... Not sure if I need to use this
   // or If I can just do what I'm doing below....
   //
   BinFabFactory<JustinsParticle> bfFactory(meshSpacing, meshOrigin);
   m_data_binfab.define(grids, 1, 1*IntVect::Unit, bfFactory);
   
   // Do I need to loop over the grids and define each BinFab<T> in the
   // LevelData? What is the "create" pointer in BinFabFactor<> used for?
   // Do I need to use that?
   //
   // How do I define a DataIndex?
   //
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      //const DataIndex& dataIndex( const DataIterator& a_dit ) const;
      //BinFab<T>*  thisbinfab_pointer = 0;
      //thisbinfab_pointer = bfFactor.create( grids[dit], 1, thisDataIndex );
      //BinFab<T>*  thisbinfab = bfFactory.create( grids[dit], 1, thisDataIndex ); 
      //BinFab<JustinsParticle>& thisbinfab = m_data_binfab[dit];
      //thisbinfab.define(grids[dit], meshSpacing, meshOrigin);
   }
   

}


PicSpecies::~PicSpecies()
{
   /*
   if(m_meshInterp!=NULL) {
      delete m_meshInterp;
      m_meshInterp = NULL;
   } 
   */
}

void PicSpecies::advancePositions( const Real& a_dt )
{
   CH_TIME("PicParticle::advancePositions()");
    
   // get some Grid info
   //
   //const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& dX(m_mesh.getdX());
   //const RealVect& meshOrigin(m_mesh.getXmin());
   //const int ghosts(m_mesh.ghosts());
   
   // compute domain extent and get total num parts and sim volume
   //
   const IntVect domainDimensions = domain.size();
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect Lbox = Xmax - Xmin; 
   
   Real maxDtinv = abs(vpiston)/dX[0];

   // Each proc loops over its own boxes, setting the initial
   // particle positions in each cell.
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      // //ListBox<Particle>& box_list = m_data[dit];
      // List<Particle>& pList = m_data[dit].listItems();
      // ListIterator<Particle> li(pList);
      List<JustinsParticle>& pList = m_data[dit].listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         RealVect&  x = li().position();
         RealVect&  v = li().velocity();
         //RealVect&  a = li().acceleration();

         // update particle position
         x += v*a_dt;

         // check to see if particle crosses piston boundary
         // or axis boundary
         //
         if(x[0]>=rpiston) { // specular refletion
            Real dr = x[0]-rpiston;
            x[0] = rpiston-dr;
            //v[0] = -2.0*abs(v[0]+vpiston);
            v[0] = -(v[0]-vpiston) + vpiston;
         }
         if(x[0]<=Xmin[0]) { // symmetry
            x[0] = -x[0];
            v[0] = abs(v[0]);
         }
         
         // set periodic boundary conditons
         //
         for(int dir=0; dir<SpaceDim; dir++) {
            if( x[dir]<Xmin[dir] && domain.isPeriodic(dir) ) {
               x[dir] = x[dir] + Lbox[dir];
            }
            if( x[dir]>Xmax[dir] && domain.isPeriodic(dir) ) {
               x[dir] = x[dir] - Lbox[dir];
            }
         }

         // set max dt for stability
         //
         Real thisDtinv = abs(v[0])/dX[0];
         maxDtinv = Max(thisDtinv,maxDtinv);
         for(int dir=1; dir<SpaceDim; dir++) {
            thisDtinv = abs(v[dir])/dX[dir];
            maxDtinv = Max(thisDtinv,maxDtinv);
         }

      }
   
   }
   m_data.gatherOutcast();
   m_data.remapOutcast();
   CH_assert(m_data.isClosed());

   // update piston position
   rpiston = rpiston + vpiston*a_dt;
   if(rpiston<dX[0]) {
      rpiston = Xmin[0] + dX[0];
      vpiston = -vpiston;
      if(!procID()) cout << "piston too close to axis !!! " << endl;
      exit(EXIT_FAILURE);
   }

   // update stable time step
   //
   Real local_stable_dt = 1.0/maxDtinv;
   Real stable_dt = local_stable_dt;
#ifdef CH_MPI
   MPI_Allreduce( &local_stable_dt, &stable_dt, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD ); 
#endif
   m_stable_dt = stable_dt; 

}

void PicSpecies::initialize()
{
   // initilize the particle position and velocities
   //

   // set ICs for this species 
   // 
   const std::string spcIC("IC." + m_name);
   ParmParse ppspcIC( spcIC.c_str() );
   std::vector<int> partsPerCellstd;
   ppspcIC.getarr("parts_per_cell",partsPerCellstd,0,SpaceDim);
   IntVect partsPerCell; // convert std::vector<int> to IntVect
   for (int dir=0; dir<SpaceDim; ++dir) {
      partsPerCell[dir] = 0;
      partsPerCell[dir] = partsPerCellstd[dir];
      CH_assert( partsPerCell[dir]>0 );
   }

   // parse the initial profiles of moments to construct
   // initial particle positions and velocities
   //
   GridFunctionFactory  gridFactory;
   const Real this_time = 0.0;
   
   // set density profile from ICs
   //
   const std::string spcdenIC("IC." + m_name + ".density");
   ParmParse ppdenIC( spcdenIC.c_str() );
   RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppdenIC,1);
   gridFunction->assign( m_density, m_mesh, this_time );

   // set temperature profiles from ICs
   //
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   LevelData<FArrayBox> tempProfile;
   tempProfile.define(grids,1,m_temperature.ghostVect());

   const std::string spctemp0IC("IC." + m_name + ".temperature_0");
   ParmParse pptemp0IC( spctemp0IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionTemp0 = gridFactory.create(pptemp0IC,1);
   gridFunctionTemp0->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(tempProfile[dit],0,0,1);
   }
   
   const std::string spctemp1IC("IC." + m_name + ".temperature_1");
   ParmParse pptemp1IC( spctemp1IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionTemp1 = gridFactory.create(pptemp1IC,1);
   gridFunctionTemp1->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(tempProfile[dit],0,1,1);
   }
   
   const std::string spctemp2IC("IC." + m_name + ".temperature_2");
   ParmParse pptemp2IC( spctemp2IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionTemp2 = gridFactory.create(pptemp2IC,1);
   gridFunctionTemp2->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(tempProfile[dit],0,2,1);
   }
   
   // set mean velocity profiles from ICs
   //
   const std::string spcvel0IC("IC." + m_name + ".velocity_0");
   ParmParse ppvel0IC( spcvel0IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionVel0 = gridFactory.create(ppvel0IC,1);
   gridFunctionVel0->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_velocity[dit].copy(tempProfile[dit],0,0,1);
   }
   
   const std::string spcvel1IC("IC." + m_name + ".velocity_1");
   ParmParse ppvel1IC( spcvel1IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionVel1 = gridFactory.create(ppvel1IC,1);
   gridFunctionVel1->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_velocity[dit].copy(tempProfile[dit],0,1,1);
   }
   
   const std::string spcvel2IC("IC." + m_name + ".velocity_2");
   ParmParse ppvel2IC( spcvel2IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionVel2 = gridFactory.create(ppvel2IC,1);
   gridFunctionVel2->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_velocity[dit].copy(tempProfile[dit],0,2,1);
   }

   if(!procID()) {
      cout << "Initializing pic species " << m_name  << endl;
   }

   ////////////////////////////////////////////////////////////////////////////////////////

   // get some mesh info
   //
   const RealVect& dX(m_mesh.getdX());
   const LevelData<FArrayBox>& Xcc(m_mesh.getXcc());
   
   int totalPartsPerCell = 1;
   Real cellVolume = 1.0;
   for (int dir=0; dir<SpaceDim; dir++)
   {
      totalPartsPerCell *= partsPerCell[dir];
      cellVolume = cellVolume*dX[dir];
   }
   Real pWeight = 0.0; 
   
   // create sub-box and dX for particles
   //
   const Box partSubBox(IntVect::Zero, partsPerCell-IntVect::Unit);
   BoxIterator pbit(partSubBox);
   
   RealVect dXpart = dX;
   if(totalPartsPerCell > 1 ) {
      for(int dir=0; dir<SpaceDim; dir++) dXpart[dir] /= partsPerCell[dir];
   }

   // loop over boxes and set the initial
   // particle values (pos., vel., weight)
   //
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for (dit.begin(); dit.ok(); ++dit) {

      CH_XD::List<JustinsParticle> thisList;
     
      // loop over grid indices
      //
      const Box gridBox = BL.get(dit);
      BoxIterator gbit(gridBox);
      for(gbit.begin(); gbit.ok(); ++gbit) {

         const IntVect ig = gbit(); // grid index
         Real local_density = m_density[dit].get(ig,0);
         pWeight = local_density*cellVolume/(Real)totalPartsPerCell; 
        
         RealVect local_Xcc; 
         RealVect local_temperature;
         RealVect local_velocity;
         for(int dir=0; dir<SpaceDim; dir++) {
            local_Xcc[dir] = Xcc[dit].get(ig,dir);
            local_temperature[dir] = m_temperature[dit].get(ig,dir);
            local_velocity[dir] = m_velocity[dit].get(ig,dir);
         }
     
         // loop over subgrid corresponding to where particles are 
         // placed in each grid cell    
         //
         for(pbit.begin(); pbit.ok(); ++pbit) {
            
            // set particle position uniformly on grid
            //
            RealVect Xpart = local_Xcc - 0.5*dX;
            IntVect ipg = pbit();
            for(int dir=0; dir<SpaceDim; dir++) {
               Xpart[dir] += (ipg[dir] + 0.5)*dXpart[dir];
            }

            // initialize particle velocities by randomly sampling 
            // a maxwellian
            //
            RealVect Vpart = RealVect::Zero;
            Real thisRand, thisVT;
            for(int dir=0; dir<SpaceDim; dir++) { 
               thisRand = MathUtils::rand();
               Vpart[dir] = MathUtils::errorinv(2.0*thisRand-1.0);
               thisVT = 4.19e5*sqrt(local_temperature[dir]/m_mass); // [m/s]
               Vpart[dir] = Vpart[dir]*sqrt(2.0)*thisVT + local_velocity[dir];
            } 
            
            // create this particle and append it to the list
            //
            JustinsParticle particle(pWeight, Xpart, Vpart);
 
            // set velocities in virtual directions
            //
            if(SpaceDim<3){
               Real thisVpart, local_temp_virt, local_vel_virt;
               for(int dir=0; dir<3-SpaceDim; dir++) { 
                  thisRand = MathUtils::rand();
                  thisVpart = MathUtils::errorinv(2.0*thisRand-1.0);
                  local_temp_virt = m_temperature[dit].get(ig,SpaceDim+dir);
                  local_vel_virt = m_velocity[dit].get(ig,SpaceDim+dir);
                  thisVT = 4.19e5*sqrt(local_temp_virt/m_mass); // [m/s]
                  thisVpart = thisVpart*sqrt(2.0)*thisVT + local_vel_virt;
                  particle.setVelocityVirt(thisVpart, dir);
               }  
            }
            
            // append particle to the list
            //
            thisList.append(particle);
 
         }

      }

      m_data_binfab[dit].addItems(thisList); // JRA testing binfab

      // finally, add particles destructively to this ListBox. Those that are
      // left behind are outcasts.
      m_data[dit].addItemsDestructive(thisList);
      
      // If all the particles were contained on the Box that 
      // created them, thisList would now be empty. But, 
      // because we perturb the initial particle distribution,
      // particles can move off the box that created them. We 
      // add these to the outcast list to be re-distributed later. 
      m_data.outcast().catenate(thisList); // no iterator? outside loop?

   }
   //const int numOutcast = m_PNew.numOutcast();
   m_data.remapOutcast();
   CH_assert(m_data.isClosed());
   
   if(!procID()) {
      cout << "Finished initializing pic species " << m_name  << endl;
   }

   /*
   //  this code is depricated...remove here and in MeshInterp soon!!!
   //
   const DisjointBoxLayout& grids = m_density.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      if(!procID()) cout << "JRA: grids[dit] = " << grids[dit] << endl;
      const FArrayBox& this_density = m_density[dit];
      ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = box_list.listItems();
      m_meshInterp.setWeightFromGridProfile( pList,
                                             this_density,
                                             totalPartsPerCell ); 
   }
   */

}

void PicSpecies::setNumberDensity()
{
   CH_TIME("PicParticle::setNumberDensity()");
    
   CH_assert(m_data.isClosed());
    
   const DisjointBoxLayout& grids = m_density.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_rho = m_density[dit];
      box_rho.setVal(0.0);
      //const ListBox<Particle>& box_list = m_data[dit];
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      //const List<Particle>& pList = m_data[dit].listItems();
      /*
      int interpFlag = 0;
      InterpType interpMethod = (InterpType)interpFlag;
      //m_meshInterp->deposit( box_list.listItems(),
      m_meshInterp.deposit( box_list.listItems(),
                            box_rho,
                            interpMethod ); 
      // NOTE that for m_meshInterp being const ref, I have to remove const
      // qualifiers in 4 places in MeshInterp files. !!!!!!!!!!!!!!!!!!!
      */
      MomentType thisMoment = density;
      m_meshInterp.moment( box_rho,
                           box_list.listItems(),
                           m_mass,
                           thisMoment ); 

   }
   m_density.exchange(); // causes ERROR: corrupted double-linked list at code exit!!!!   
     
}

void PicSpecies::setMomentumDensity()
{
   CH_TIME("PicParticle::setMomentumDensity()");
    
   CH_assert(m_data.isClosed());
    
   const DisjointBoxLayout& grids = m_momentum.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_mom = m_momentum[dit];
      box_mom.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = momentum;
      m_meshInterp.moment( box_mom,
                           box_list.listItems(),
                           m_mass,
                           thisMoment ); 

   }
   m_momentum.exchange(); 
     
}

void PicSpecies::setEnergyDensity()
{
   CH_TIME("PicParticle::setEnergyDensity()");
    
   CH_assert(m_data.isClosed());
    
   const DisjointBoxLayout& grids = m_energy.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& box_ene = m_energy[dit];
      box_ene.setVal(0.0);
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      MomentType thisMoment = energy;
      m_meshInterp.moment( box_ene,
                           box_list.listItems(),
                           m_mass,
                           thisMoment ); 

   }
   m_energy.exchange(); 
     
}

void PicSpecies::numberDensity( LevelData<FArrayBox>&  a_rho )
{
   CH_TIME("PicParticle::numberDensity()");
 
   setNumberDensity();
   const DisjointBoxLayout& grids = a_rho.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      a_rho[dit].copy(m_density[dit]);   
   }
  
}

void PicSpecies::momentumDensity( LevelData<FArrayBox>&  a_mom )
{
   CH_TIME("PicParticle::momentumDensity()");
 
   setMomentumDensity();
   const DisjointBoxLayout& grids = a_mom.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      a_mom[dit].copy(m_momentum[dit]);   
   }
  
}

void PicSpecies::energyDensity( LevelData<FArrayBox>&  a_ene )
{
   CH_TIME("PicParticle::energyDensity()");
 
   setEnergyDensity();
   const DisjointBoxLayout& grids = a_ene.disjointBoxLayout();
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      a_ene[dit].copy(m_energy[dit]);   
   }
  
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

