
#include "PicSpecies.H"
#include <array>
#include <cmath>
#include "Constants.H"
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
     m_scatter(false),
     m_mesh(a_mesh),
     //m_meshInterp(NULL)
     m_meshInterp(a_meshInterp)
{
   // maybe parse some stuff
   m_name = a_name;
 
   //createMeshInterp();
 
   a_ppspc.query( "mass", m_mass );
   a_ppspc.query( "charge", m_charge );
   a_ppspc.query( "motion", m_motion );
   a_ppspc.query( "forces", m_forces );
   a_ppspc.query( "scatter", m_scatter );

   if ( procID() == 0 ) {
      cout << "  name = " << m_name << endl;
      cout << "  mass = " << m_mass << endl;
      cout << "  charge = " << m_charge << endl;
      cout << "  motion = " << m_motion << endl;
      cout << "  forces = " << m_forces << endl;
      cout << "  scatter = " << m_scatter << endl;
   }

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   const int ghosts(m_mesh.ghosts());
   //RealVect particleOrigin = RealVect(D_DECL(0.0,0.0,0.0));
   //if(!procID()) cout << "particleOrigin = " << particleOrigin << endl; 

   m_stable_dt = meshSpacing[0]; // initialize stable time step

   // initialize the member LevelDatas
   //
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   m_density.define(grids,1,ghostVect);
   m_momentum.define(grids,3,ghostVect);
   m_energy.define(grids,3,ghostVect); // direction-depenent energy
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
   BinFabFactory<JustinsParticlePtr> bfptrFactory(meshSpacing, meshOrigin);
   m_data_binfab_ptr.define(grids, 1, 0*IntVect::Unit, bfptrFactory);

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
   
   // Each proc loops over its own boxes
   //
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      //ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = m_data[dit].listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         RealVect&  x = li().position();
         std::array<Real,3>&  v = li().velocity();

         // update particle position
         //x += v*a_dt;
         for(int dir=0; dir<SpaceDim; dir++) {
            x[dir] += v[dir]*a_dt;
         }

         // set particle boundary conditions here for now
         //
         // What I should do is create a list of particles that cross the
         // domain boundaries, and then apply the BCs to that list elsewhere
         //
         for(int dir=0; dir<SpaceDim; dir++) {
            if( domain.isPeriodic(dir) ) {
               if( x[dir]<Xmin[dir] ) x[dir] = x[dir] + Lbox[dir];
               if( x[dir]>=Xmax[dir] ) x[dir] = x[dir] - Lbox[dir];
            }
            else { // symmetry BCs
               if(x[dir]<=Xmin[dir]) {
                  x[dir] = 2.*Xmin[dir] - x[dir];
                  v[dir] = -v[dir];
               }
               if(x[dir]>=Xmax[dir]) {
                  x[dir] = 2.*Xmax[dir] - x[dir];
                  v[dir] = -v[dir];
               }
            }
         }
         //
         //
         //////////////////////////////////////////////////

      } // end loop over particle list
      
   } // end loop over boxes

   m_data.gatherOutcast();
   m_data.remapOutcast();
   CH_assert(m_data.isClosed());
   
}


void PicSpecies::setStableDt()
{
   CH_TIME("PicParticle::setStableDt()");
   
   // set the stable time step based on particles crossing a grid cell
   //
   const RealVect& dX(m_mesh.getdX());
   Real maxDtinv = 0.0;
   Real thisDtinv;

   // Each proc loops over its own boxes
   //
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for(dit.begin(); dit.ok(); ++dit) {

      //ListBox<JustinsParticle>& box_list = m_data[dit];
      List<JustinsParticle>& pList = m_data[dit].listItems();
      ListIterator<JustinsParticle> li(pList);
      for(li.begin(); li.ok(); ++li) {

         std::array<Real,3>&  v = li().velocity();
         for(int dir=0; dir<SpaceDim; dir++) {
            thisDtinv = abs(v[dir])/dX[dir];
            maxDtinv = Max(thisDtinv,maxDtinv);
         }

      }
      
   }

   // update stable time step
   //
   Real local_stable_dt = 1.0/maxDtinv;
   Real stable_dt = local_stable_dt;
#ifdef CH_MPI
   MPI_Allreduce( &local_stable_dt, &stable_dt, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD ); 
#endif
   m_stable_dt = stable_dt; 
   //if(!procID()) cout << "m_particle_dt = " << m_stable_dt << endl;

}

void PicSpecies::binTheParticles()
{
   CH_TIME("PicParticle::binTheParticles()");
   
   ///////////////////////////////////////////////////////////
   //
   //   fill BinFab container with pointers to particles
   //
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);
   for (dit.begin(); dit.ok(); ++dit) { // loop over boxes
      //m_data_binfab_ptr[dit].reBin(); 
      const Box gridBox = BL.get(dit);
      BoxIterator gbit(gridBox);
      List<JustinsParticle>& pList = m_data[dit].listItems();
      
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = m_data_binfab_ptr[dit];
      
      // clear the binfab_ptr container
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
         pListPtr.append(particlePtr); // JRA PTR
      }
      thisBinFab_ptr.addItems(pListPtr); // JRA testing binfab_ptr
   }
}

void PicSpecies::testParticleShuffling( const Real& a_dt )
{
   CH_TIME("PicParticle::testParticleShuffling()");
  
    JustinsParticle* this_part_ptr = NULL;  
   
   // fill BinFab container with pointers to particles
   //
   binTheParticles();
   
   // loop over lists in each cell and test shuffle
   //
   const DisjointBoxLayout& grids = m_data_binfab_ptr.disjointBoxLayout();
   DataIterator ditg(grids);
   int thisNumCell;
   int verbosity=0; // using this as a verbosity flag
   for (ditg.begin(); ditg.ok(); ++ditg) { // loop over boxes
   
      const Box gridBox = grids.get(ditg);
      BoxIterator gbit(gridBox);
      
      BinFab<JustinsParticlePtr>& thisBinFab_ptr = m_data_binfab_ptr[ditg];
      std::vector<JustinsParticlePtr> shuffled_parts;
      for (gbit.begin(); gbit.ok(); ++gbit) { // loop over grid indices
         
         const IntVect ig = gbit(); // grid index
         List<JustinsParticlePtr>& cell_pList = thisBinFab_ptr(ig,0);
         thisNumCell = cell_pList.length();
         if(thisNumCell < 2) break;        
         if(!procID() && verbosity) cout << "JRA: thisNumCell = " << thisNumCell << endl;
 
         //cell_pList.shuffle();
         //thisNumCell = cell_pList.length();
         
         ListIterator<JustinsParticlePtr> lit(cell_pList);
           
         // copy the iterators to a vector
         //std::vector<JustinsParticlePtr> shuffled_parts;
         shuffled_parts.clear();
         shuffled_parts.reserve(thisNumCell);
         for (lit.begin(); lit.ok(); ++lit) shuffled_parts.push_back(lit());
        
         // randomly choose two unique elements from the vector
         int random_index1, random_index2; 
         for (auto n=0; n<thisNumCell; n++) {
            random_index1 = std::rand() % thisNumCell;  
            random_index2 = std::rand() % thisNumCell;
            while(random_index2==random_index1) random_index2 = std::rand() % thisNumCell;  
            if(procID()==0 && verbosity) {
               cout << "JRA random_index1 = " << random_index1 << endl;
               cout << "JRA random_index2 = " << random_index2 << endl;
            }
         }

         // randomly shuffle the entire vector of iterators 
         std::random_shuffle(shuffled_parts.begin(), shuffled_parts.end());
         
         // one way to loop over vector
         for (auto it = shuffled_parts.begin(); it != shuffled_parts.end(); ++it) {  
            JustinsParticlePtr& this_particle = (*it);
            this_part_ptr = this_particle.getPointer();
            const uint64_t& this_ID = this_part_ptr->ID();
            if(procID()==0 && verbosity) {
               cout << "JRA (shuffled): ID = " << this_ID << endl;
            }
         }
         
         // another way to loop over vector
         for (auto n=0; n<shuffled_parts.size(); n++) {  
            JustinsParticlePtr& this_particle = shuffled_parts[n];
            this_part_ptr = this_particle.getPointer();
            const uint64_t& this_ID = this_part_ptr->ID();
            if(procID()==0 && verbosity) {
               cout << "JRA (shuffled 2): ID = " << this_ID << endl;
            }
         }
        
         // loop over the particles in this cell using the iterator 
         for (lit.begin(); lit.ok(); ++lit) {
            JustinsParticlePtr& this_particle = lit();
            this_part_ptr = this_particle.getPointer();
            const uint64_t& this_ID = this_part_ptr->ID();
            const RealVect& this_x = this_part_ptr->position();
            if(procID()==0 && verbosity) {
               cout << "JRA: ID = " << this_ID << endl;
               cout << "JRA: position = " << this_x << endl;
            }
         }
         verbosity=0;
      }

   }

   this_part_ptr = NULL;
   delete this_part_ptr;
   //
   //
   ////////////////////////////////// 
   
}

void PicSpecies::initialize()
{
   // initilize the particle position and velocities
   //
   if(!procID()) {
      cout << "Initializing pic species " << m_name  << "..." << endl;
   }
   int verbosity=0; // using this as a verbosity flag

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
   RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppdenIC,verbosity);
   gridFunction->assign( m_density, m_mesh, this_time );

   // set temperature profiles from ICs
   //
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   LevelData<FArrayBox> tempProfile;
   tempProfile.define(grids,1,m_temperature.ghostVect());

   const std::string spctemp0IC("IC." + m_name + ".temperature_0");
   ParmParse pptemp0IC( spctemp0IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionTemp0 = gridFactory.create(pptemp0IC,verbosity);
   gridFunctionTemp0->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(tempProfile[dit],0,0,1);
   }
   
   const std::string spctemp1IC("IC." + m_name + ".temperature_1");
   ParmParse pptemp1IC( spctemp1IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionTemp1 = gridFactory.create(pptemp1IC,verbosity);
   gridFunctionTemp1->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(tempProfile[dit],0,1,1);
   }
   
   const std::string spctemp2IC("IC." + m_name + ".temperature_2");
   ParmParse pptemp2IC( spctemp2IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionTemp2 = gridFactory.create(pptemp2IC,verbosity);
   gridFunctionTemp2->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(tempProfile[dit],0,2,1);
   }
   
   // set mean velocity profiles from ICs
   //
   const std::string spcvel0IC("IC." + m_name + ".velocity_0");
   ParmParse ppvel0IC( spcvel0IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionVel0 = gridFactory.create(ppvel0IC,verbosity);
   gridFunctionVel0->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_velocity[dit].copy(tempProfile[dit],0,0,1);
   }
   
   const std::string spcvel1IC("IC." + m_name + ".velocity_1");
   ParmParse ppvel1IC( spcvel1IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionVel1 = gridFactory.create(ppvel1IC,verbosity);
   gridFunctionVel1->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_velocity[dit].copy(tempProfile[dit],0,1,1);
   }
   
   const std::string spcvel2IC("IC." + m_name + ".velocity_2");
   ParmParse ppvel2IC( spcvel2IC.c_str() );
   RefCountedPtr<GridFunction> gridFunctionVel2 = gridFactory.create(ppvel2IC,verbosity);
   gridFunctionVel2->assign( tempProfile, m_mesh, this_time );
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_velocity[dit].copy(tempProfile[dit],0,2,1);
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
   uint64_t ID = procID()*512 + 1; // hack for testing purposes
   //Real ID = 0.;
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
         std::array<Real,3> local_temperature;
         std::array<Real,3> local_velocity;
         for(int dir=0; dir<SpaceDim; dir++) {
            local_Xcc[dir] = Xcc[dit].get(ig,dir);
         }
         for(int dir=0; dir<3; dir++) {
            local_temperature[dir] = m_temperature[dit].get(ig,dir);
            local_velocity[dir] = m_velocity[dit].get(ig,dir);
         }
     
         // loop over subgrid corresponding to where particles are 
         // placed in each grid cell    
         //
         Real V0 = sqrt(Constants::QE/Constants::ME); // ele thermal speed at 1eV [m/s]
         //Real beta0 = V0/Constants::CVAC;
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
            std::array<Real,3> Vpart = {0,0,0};
            Real thisRand, thisVT;
            for(int dir=0; dir<3; dir++) { 
               thisRand = MathUtils::rand();
               Vpart[dir] = MathUtils::errorinv(2.0*thisRand-1.0);
               thisVT = V0*sqrt(local_temperature[dir]/m_mass); // [m/s]
               Vpart[dir] = Vpart[dir]*sqrt(2.0)*thisVT + local_velocity[dir];
            } 
            
            // create this particle and append it to the list
            //
            JustinsParticle particle(pWeight, Xpart, Vpart);
            particle.setID(ID);
            ID = ID + 1;
 
            // append particle to the list
            //
            thisList.append(particle);
            
         }

      }

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
      cout << "Finished initializing pic species " << m_name  << endl << endl;
   }

   // JRA, testing what's in binfab...
   //
   //testParticleShuffling( 0.0 );
   //inspectBinFab( m_data_binfab_ptr );

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
      const ListBox<JustinsParticle>& box_list = m_data[dit];
      
      //const List<JustinsParticle>& pList = box_list.listItems();
      //if(procID()==0 || procID()==8) cout << "JRA: pList.length() = " << pList.length() << endl;
     
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

