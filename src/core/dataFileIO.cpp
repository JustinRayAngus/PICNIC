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

   if(!procID()) cout << "writing mesh data file ..." << endl << endl;

#ifdef CH_USE_HDF5

   //char iter_str[100];
   const std::string plotFileName = "mesh.h5";
   HDF5Handle handle( plotFileName.c_str(), HDF5Handle::CREATE );
  
   handle.setGroup("/");

   // write the header stuff
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
   //handle.close();

   // now write the data
   //

   // setup the level string
   char levelStr[20];
   const int this_level = 0;
   sprintf(levelStr,"%d",this_level);
   const std::string label = std::string("level_") + levelStr;
  
   handle.setGroup(label);
   //header.m_real["dx"] = m_dX[1]; 

   const ProblemDomain& domain(m_mesh.getDomain());
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);

   // write the grid data
   //
   //DomainGrid* mesh = DomainGrid::mesh;
   //DisjointBoxLayout& grids = mesh->m_grids;
   //const DisjointBoxLayout& grids(m_mesh.getDBL());
   const LevelData<FArrayBox>& outputData(m_mesh.getXcc());
   
   //LevelData<FArrayBox> outputData;
   //outputData.define(grids, SpaceDim, IntVect::Zero);   
   //m_DomainGrid.getXphys(outputData);
   //m_fieldNew.copyTo(Interval(0, SpaceDim - 1), outputData, Interval(0, SpaceDim - 1));

   write(handle, outputData.boxLayout());
   write(handle, outputData, "data", outputData.ghostVect());
   
   // now write some more data
   //
   /*
   const LevelData<FluxBox>& outputData2(m_mesh.getXfc());
   const std::string label1 = std::string("level_1");
   handle.setGroup(label1);
   header.m_int["num_components"] = outputData2.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   */
   //LevelData<FArrayBox> outputData2comp1;
   //outputData2comp1.define(grids,1,IntVect::Unit);
   //for(DataIterator dit(grids); dit.ok(); ++dit) {
   //   outputData2comp1[dit].copy(outputData2[dit][0],0,0,1); 
   //}
   //write(handle, outputData2comp1.boxLayout());
   //write(handle, outputData2comp1, "data",IntVect::Unit);
   //write(handle, outputData2.boxLayout());
   //write(handle, outputData2,"data",outputData2.ghostVect());
   
   // close the handle
   //
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

void dataFileIO::writeParticleDataFile( const ParticleData<JustinsParticle>&  a_Pdata,
                                        const LevelData<FArrayBox>&    a_density, 
                                        const LevelData<FArrayBox>&    a_momentum, 
                                        const LevelData<FArrayBox>&    a_energy, 
                                        const Real                     a_mass, 
                                        const int                      a_cur_step,
                                        const double                   a_cur_time )
{
   CH_TIME("dataFileIO::writeParticleDataFile()");
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
   std::string prefix = "particle";
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
   coords[0] = 'x';
   coords[1] = 'y';
   coords[2] = 'z';

   // first write header for mesh data
   //
   int numMeshComps = a_density.nComp() + a_momentum.nComp() + a_energy.nComp();
   vectNames.push_back("density");
   for (int dir=0; dir<a_momentum.nComp(); dir++) {
      sprintf(field_name, "momentumDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<a_energy.nComp(); dir++) {
      sprintf(field_name, "energyDensity_%c", coords[dir]);
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
   for (int dir=0; dir<SpaceDim; dir++) {
      sprintf(field_name, "particle_position_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   if(SpaceDim<3) {
     for(int i=0; i<3-SpaceDim; i++) {
        sprintf(field_name, "virtual_particle_position_%c", i);
        vectNames.push_back(field_name);
     }
   }
   for (int dir=0; dir<3; dir++) {
      sprintf(field_name, "particle_velocity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<3; dir++) {
      sprintf(field_name, "particle_acceleration_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   vectNames.push_back("particle_weight");
   vectNames.push_back("particle_ID");

   int numComps = vectNames.size();
   for (int i=0; i<numComps; ++i) {
      sprintf(comp_name, "particle_component_%d", i);
      headerParts.m_string[comp_name] = vectNames[i];
   }
 
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
   headerParts.m_real["mass"] = a_mass;
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_box["prob_domain"] = domain.domainBox();
   
   // write the header 
   headerParts.writeToFile(handleParts);
 
   // write the particles
   write(handleParts, grids);
   writeParticlesToHDF(handleParts, a_Pdata, "particles");


   // write the moment data
   LevelData<FArrayBox> momentData;
   momentData.define(grids,numMeshComps,a_density.ghostVect());
   int compData = 0;
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      momentData[dit].copy(a_density[dit],0,compData,1);
      compData = compData + 1;
      for (int dir=0; dir<a_momentum.nComp(); dir++) { 
         momentData[dit].copy(a_momentum[dit],dir,compData,1);
         compData = compData + 1;
      } 
      for (int dir=0; dir<a_energy.nComp(); dir++) { 
         momentData[dit].copy(a_energy[dit],dir,compData,1);
         compData = compData + 1;
      } 
   }
   write(handleParts, momentData.boxLayout());
   write(handleParts, momentData, "data", a_density.ghostVect());

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
