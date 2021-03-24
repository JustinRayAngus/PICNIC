#include "dataFileIO.H"
#include "CH_Timer.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "BoxIterator.H"
#include "DomainGrid.H"
#include "ParticleIO.H"

#include "NamespaceHeader.H"


dataFileIO::dataFileIO( ParmParse&   a_pp,
                  const DomainGrid&  a_mesh )
   :
     m_plot_particle_moments(false),
     m_mesh(a_mesh),
     m_verbosity(0)
{

   // parse stuff from input about output directory names
   // and different plotting options for different types of data
   //
   
   writeMeshDataFile();

}

void dataFileIO::writeMeshDataFile()
{
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   if(!procID()) cout << "writing mesh data file ..." << endl << endl;
   const ProblemDomain& domain(m_mesh.getDomain());

#ifdef CH_USE_HDF5

   const std::string plotFileName = "mesh.h5";
   HDF5Handle handle( plotFileName.c_str(), HDF5Handle::CREATE );
  
   handle.setGroup("/");

   //
   // write the default header
   //
  
   // set the component names
   HDF5HeaderData header;

   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[SpaceDim];
   coords[0] = 'x';
   if(SpaceDim==2) coords[1] = 'z';
   if(SpaceDim==3) {
      coords[1] = 'y';
      coords[2] = 'z';
   }
   int numMeshComps = SpaceDim;
   for (int dir=0; dir<SpaceDim; dir++)
   {
      sprintf(field_name, "grid_%c", coords[dir]);
      vectNames.push_back(field_name);
   }

   for (int i = 0; i < numMeshComps; ++i)
   {
     sprintf(comp_name, "component_%d", i);
     header.m_string[comp_name] = vectNames[i];
   }
   header.m_int["num_components"] = numMeshComps;

   header.writeToFile(handle);

   //
   // write the cell centered grid data
   //
   
   const LevelData<FArrayBox>& gridOnCells(m_mesh.getXcc());
   
   const std::string group1Name= std::string("cell_centered_grid");
   handle.setGroup(group1Name);

   header.clear();
   header.m_int["num_components"] = gridOnCells.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);

   
   write(handle, gridOnCells.boxLayout());
   write(handle, gridOnCells, "data", gridOnCells.ghostVect());
   
   //
   // write the face centered grid data
   //
   
   const LevelData<FluxBox>& gridOnFaces(m_mesh.getXfc());

   const std::string group2Name = std::string("face_centered_grid");
   handle.setGroup(group2Name);
   
   header.clear();
   header.m_int["is_fluxbox"] = 1;
   header.m_int["num_components"] = gridOnFaces.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, gridOnFaces.boxLayout());
   write(handle, gridOnFaces, "data", gridOnFaces.ghostVect());
   
   //
   // write the edge centered grid data
   //
   // Some issue with with writing edge data box...
   //
 
   /*
   const LevelData<EdgeDataBox>& gridOnEdges(m_mesh.getXec());
   
   //const DisjointBoxLayout& grids = gridOnEdges.getBoxes();
   //if(!procID()) cout << "gridOnEdges.nComp() = " << gridOnEdges.nComp() << endl;
   //for(DataIterator dit(grids); dit.ok(); ++dit) {
   //   cout << "procID() = " << procID() << endl;
   //   cout << "box() " << grids.get(dit) << endl;
   //   for(int dir=0; dir<SpaceDim; dir++) {
   //      cout << "dir =  " << dir << " box = " << gridOnEdges[dit][dir].box() << endl;
   //   }
   //}
   
   const std::string group3Name = std::string("edge_centered_grid");
   handle.setGroup(group3Name);
   
   header.clear();
   header.m_int["is_edgebox"] = 1;
   header.m_int["num_components"] = gridOnEdges.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, gridOnEdges.boxLayout());
   write(handle, gridOnEdges, "data", gridOnEdges.ghostVect());
   */
 
   // close the handle
   handle.close();
   
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif

}



inline
std::string dirPrefix( const std::string& a_prefix,
                       const std::string& a_diag_name )
{
   std::string dir_prefix( a_prefix );
#if 1  // warning, OS dependencies, will not work on all platforms
   //std::string iter_str( a_prefix + "_" + a_diag_name + "_plots" );
   std::string iter_str( a_prefix + "_data" );
#ifdef CH_MPI
   if (procID() == 0) {
#endif
      // only works the first time, subsequent failure is normal and expected
      mkdir( iter_str.c_str(), 0777 );
#ifdef CH_MPI
   }
#endif
   dir_prefix = std::string( iter_str + "/" );
#endif
   return dir_prefix;
}


inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const int a_cur_step )
{  
   std::string dir_prefix( dirPrefix( a_prefix, a_diag_name ) );
   char buffer[100];
   sprintf( buffer, "%04d.", a_cur_step );
   //std::string filename( dir_prefix + a_prefix + "." + a_diag_name + buffer );
   std::string filename( dir_prefix + a_diag_name + buffer + "h5" );
   return filename;
}


inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const std::string& a_species_name,
                          const int a_cur_step,
                          const int a_species_index )
{
   std::string dir_prefix( dirPrefix( a_prefix, a_diag_name ) );
   char buffer0[10];
   char buffer[20];
   sprintf( buffer0, ".%d.", a_species_index);
   sprintf( buffer, "%04d", a_cur_step );
   std::string filename( dir_prefix + a_prefix + buffer0 + a_species_name + "." + a_diag_name + buffer + "." );

   return filename;
}

void dataFileIO::writeParticleDataFile( const PicSpecies&  a_picSpecies,
                                        const int          a_species,
                                        const int          a_cur_step,
                                        const double       a_cur_time )
{
   CH_TIME("dataFileIO::writeParticleDataFile()");
   
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   // get references to particle data   
   const ParticleData<JustinsParticle>& a_Pdata = a_picSpecies.partData();
   
   // It is the job of the caller to maker sure the moments are set
   const LevelData<FArrayBox>& density  = a_picSpecies.getNumberDensity();
   const LevelData<FArrayBox>& momentum = a_picSpecies.getMomentumDensity();
   const LevelData<FArrayBox>& energy   = a_picSpecies.getEnergyDensity();
   const LevelData<FArrayBox>& current  = a_picSpecies.getTestDeposit();
   
   if(!procID()) {
      cout << "writing particle data file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   //const RealVect& meshSpacing(m_mesh.getdX());
   
   ///////////// 
   //
   //   write the particles
   //
   //std::string plotFileNameParts = "parts.h5";
   stringstream s;
   s << "species" << a_species; 
   std::string prefix = s.str();
   //std::string prefix = "particle";
   std::string plotFileNameParts( plotFileName( prefix,
                                                "parts",
                                                a_cur_step ) );

   HDF5Handle handleParts( plotFileNameParts.c_str(), HDF5Handle::CREATE );
   
   handleParts.setGroup("/");

   // write the header stuff
   //
   HDF5HeaderData headerParts;
   
   // set the component names
   //
   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';

   // first write header for mesh data
   //
   int numMeshComps = density.nComp() + momentum.nComp() + energy.nComp();
   numMeshComps = numMeshComps + current.nComp();
   vectNames.push_back("density");
   for (int dir=0; dir<momentum.nComp(); dir++) {
      sprintf(field_name, "momentumDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<energy.nComp(); dir++) {
      sprintf(field_name, "energyDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<current.nComp(); dir++) {
      sprintf(field_name, "currentDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for(int i=0; i<numMeshComps; i++) {
      sprintf(comp_name,"component_%d", i);
      headerParts.m_string[comp_name] = vectNames[i];
   }
   headerParts.m_int["num_components"] = numMeshComps;


   // now write header for particles
   //
   vectNames.clear();
   vectNames.push_back("particle_weight");
   for (int dir=0; dir<SpaceDim; dir++) {
      sprintf(field_name, "particle_position_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   if(SpaceDim<3) {
     for(int i=0; i<3-SpaceDim; i++) {
        sprintf(field_name, "virtual_particle_position_%c", coords[i]);
        vectNames.push_back(field_name);
     }
   }
   for (int dir=0; dir<3; dir++) {
      sprintf(field_name, "particle_velocity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   if(a_picSpecies.writeAll()) {
      for (int dir=0; dir<3; dir++) {
         sprintf(field_name, "particle_electric_field_%c", coords[dir]);
         vectNames.push_back(field_name);
      }
   }
   vectNames.push_back("particle_ID");

   int numComps = vectNames.size();
   for (int i=0; i<numComps; ++i) {
      sprintf(comp_name, "particle_component_%d", i);
      headerParts.m_string[comp_name] = vectNames[i];
   }
   headerParts.m_int["numPartComps"] = numComps;
 
   headerParts.writeToFile(handleParts);


   //
   ///////////////////////////////////// 


   char levelStr[20];
   const int this_level = 0;
   sprintf(levelStr,"%d",this_level);
   const std::string label = std::string("level_") + levelStr;
   handleParts.setGroup(label);  


   int totalParticleCount = a_Pdata.numParticles();
   headerParts.m_int["num_particles"] = totalParticleCount;
   headerParts.m_real["mass"] = a_picSpecies.mass();
   headerParts.m_int["charge"] = a_picSpecies.charge();
   headerParts.m_real["Uint"] = a_picSpecies.Uint();
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_box["prob_domain"] = domain.domainBox();
   
   // write the header 
   headerParts.writeToFile(handleParts);
 
   // write the particles
   write(handleParts, grids);
   writeParticlesToHDF(handleParts, a_Pdata, "particles", a_picSpecies.writeAll() );


   // write the moment data
   LevelData<FArrayBox> momentData;
   momentData.define(grids,numMeshComps,density.ghostVect());
   int compData = 0;
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      momentData[dit].copy(density[dit],0,compData,1);
      compData = compData + 1;
      for (int dir=0; dir<momentum.nComp(); dir++) { 
         momentData[dit].copy(momentum[dit],dir,compData,1);
         compData = compData + 1;
      } 
      for (int dir=0; dir<energy.nComp(); dir++) { 
         momentData[dit].copy(energy[dit],dir,compData,1);
         compData = compData + 1;
      } 
      for (int dir=0; dir<current.nComp(); dir++) { 
         momentData[dit].copy(current[dit],dir,compData,1);
         compData = compData + 1;
      } 
   }
   write(handleParts, momentData.boxLayout());
   write(handleParts, momentData, "data", density.ghostVect());

   handleParts.close();
   
   if(!procID()) cout << "finished writing particle data file" << endl << endl;

}

void dataFileIO::writeBinFabDataFile( const LevelData<BinFab<JustinsParticle>>&  a_Pdata,
                                      const int                      a_cur_step,
                                      const double                   a_cur_time )
{
   CH_TIME("dataFileIO::writeBinFabDataFile()");
   if(!procID()) {
      cout << "writing BinFab particle data file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   //const RealVect& meshSpacing(m_mesh.getdX());
   
   ///////////// 
   //
   //   write the particles
   //
   //std::string plotFileNameParts = "parts.h5";
   std::string prefix = "binfab";
   std::string plotFileNameParts( plotFileName( prefix,
                                                "parts",
                                                a_cur_step ) );

   HDF5Handle handleParts( plotFileNameParts.c_str(), HDF5Handle::CREATE );
   
   handleParts.setGroup("/");
   
   // There is no function to write BinFab???????????
   //write(handleParts, grids);
   //writeParticlesToHDF(handleParts, a_Pdata, "particles");

   handleParts.close();
   
   if(!procID()) cout << "finished writing binfab particle data file" << endl << endl;

}


#include "NamespaceFooter.H"
