
#include "PicSpecies.H"
#include <array>
#include <cmath>
#include "BoxIterator.H"
#include "ProblemDomain.H"

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
   m_velocity.define(grids,1,ghostVect);
   m_temperature.define(grids,1,ghostVect);
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      m_density[dit].setVal(0.0);
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
   Real density, meanVelocity, temperature;
   ppspcIC.query( "density", density );
   ppspcIC.query( "meanVelocity", meanVelocity );
   ppspcIC.query( "temperature", temperature );

   if(!procID()) {
      cout << "PicSpecies::inititialize() Setting initial conditions for species " << endl;
      cout << "density = " << density << endl;
      cout << "meanVelocity = " << meanVelocity << endl;
      cout << "temperature = " << temperature << endl;
      cout << "particles per cell = " << partsPerCell[0] << endl;
   }

   // get some mesh info
   //
   //const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   //const RealVect& meshOrigin(m_mesh.getXmin());
   //const int ghosts(m_mesh.ghosts());
   
   // compute domain extent and get total num parts and sim volume
   //
   const IntVect domainDimensions = domain.size();
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect Lbox = Xmax - Xmin; 
   
   int numParticles = 1;
   Real volume = 1.0;
   for (int dir=0; dir<SpaceDim; dir++)
   {
      numParticles *= partsPerCell[dir]*domainDimensions[dir];
      volume *= Lbox[dir];
   }
   Real pWeight = density*volume/(Real)numParticles; 
   
   if(!procID()) {
      cout << "JRA: numParticles = " << numParticles << endl;
      cout << "JRA: volume = " << volume << endl << endl;
      cout << "JRA: pWeight = " << pWeight << endl << endl;
   }
   

   // Each proc loops over its own boxes, setting the initial
   // particle positions in each cell.
   const BoxLayout& BL = m_data.getBoxes();
   DataIterator dit(BL);

   for (dit.begin(); dit.ok(); ++dit) {

      //CH_XD::List<Particle> thisList;
      CH_XD::List<JustinsParticle> thisList;

      // refine the box so we can iterate over cell centers
      // and get the right number of particles per cell
      // 
      RealVect dXpart = meshSpacing;
      const Box thisBox = BL.get(dit);
      Box partBox = thisBox;
      
      const int totalPartsPerCell = partsPerCell[0]*partsPerCell[1];
      if(totalPartsPerCell > 1 ) {
         partBox.refine(partsPerCell);
         for(int dir=0; dir<SpaceDim; dir++) dXpart[dir] /= partsPerCell[dir];
      }

      BoxIterator bit(partBox);
      for (bit.begin(); bit.ok(); ++bit) {

          // index of this cell
          IntVect iv = bit();
          
          RealVect Xpart = Xmin;
          RealVect Vpart = 0.0*Xmin;  // init all velocity to zero
          for (int dir=0; dir<SpaceDim; dir++) {
             Xpart[dir] += (iv[dir] + 0.5)*dXpart[dir];
          }

          // enforce periodic boundary conditions
          // check for perodic direction ? 
          // domain.isPeriodic[dir]
          //
          for(int dir=0; dir<SpaceDim; dir++) {

             if(Xpart[dir]<Xmin[dir]) {
                Xpart[dir] += Lbox[dir];
             }
             if(Xpart[dir]>Xmax[dir]) {
                Xpart[dir] -= Lbox[dir];
             }

          }

          /////////////////////////////////////////////////////

          // initialize particle velocities by randomly sampling 
          // a maxwellian
          //
          double rand = MathUtils::rand();
          Vpart[1] = MathUtils::errorinv(2.0*rand-1.0);
          
          /////////////////////////////////////////////////////

          // create this particle and append to list
          //Particle particle(pWeight, Xpart, Vpart);
          JustinsParticle particle(pWeight, Xpart, Vpart);
          thisList.append(particle);
          //if(!procID()) cout << "pos_virt = " << particle.pos_virt() << endl; 
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
   //const RealVect meshSpace = m_PNew.meshSpacing();
   //const RealVect thisOrigin = m_PNew.origin();
   //const ProblemDomain& thisPhysDomain= m_PNew.physDomain();
   m_data.remapOutcast();
   CH_assert(m_data.isClosed());

   setNumberDensity();

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
      
      int interpFlag = 0;
      InterpType interpMethod = (InterpType)interpFlag;
      //m_meshInterp->deposit( box_list.listItems(),
      m_meshInterp.deposit( box_list.listItems(),
                            box_rho,
                            interpMethod ); 
      // NOTE that for m_meshInterp being const ref, I have to remove const
      // qualifiers in 4 places in MeshInterp files. !!!!!!!!!!!!!!!!!!!

   }
   m_density.exchange(); // causes ERROR: corrupted double-linked list at code exit!!!!   
     
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

