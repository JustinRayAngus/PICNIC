#include "dataFileIO.H"
#include "CH_Timer.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "BoxIterator.H"
#include "DomainGrid.H"
#include "ParticleIO.H"

#include "SpaceUtils.H"

#include "NamespaceHeader.H"


dataFileIO::dataFileIO( ParmParse&   a_ppsys,
                  const DomainGrid&  a_mesh,
                  const CodeUnits&   a_units )
   :
     m_plot_particle_moments(false),
     m_mesh(a_mesh),
     m_units(a_units),
     m_verbosity(0)
{

   // parse stuff from input about output directory names
   // and different plotting options for different types of data
   //
  
   m_nodeDataTest = false;
   a_ppsys.query("node_data_test",m_nodeDataTest);
   if(!procID()) cout << "node_data_test = " << m_nodeDataTest << endl; 

   writeMeshDataFile();

}

inline
std::string dirPrefix( const std::string& a_prefix )
{
   std::string dir_prefix( a_prefix );
#if 1  // warning, OS dependencies, will not work on all platforms
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
                          const std::string& a_diag_name )
{  
   std::string dir_prefix( dirPrefix( a_prefix ) );
   std::string filename( dir_prefix + a_diag_name + ".h5" );
   return filename;
}

inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const int a_cur_step )
{  
   std::string dir_prefix( dirPrefix( a_prefix ) );
   char buffer[100];
   sprintf( buffer, "%04d.", a_cur_step );
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
   std::string dir_prefix( dirPrefix( a_prefix ) );
   char buffer0[10];
   char buffer[20];
   sprintf( buffer0, ".%d.", a_species_index);
   sprintf( buffer, "%04d", a_cur_step );
   std::string filename( dir_prefix + a_prefix + buffer0 + a_species_name + "." + a_diag_name + buffer + "." );

   return filename;
}

void dataFileIO::writeMeshDataFile()
{
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   if(!procID()) cout << "writing mesh data file ..." << endl << endl;
   const ProblemDomain& domain(m_mesh.getDomain());

#ifdef CH_USE_HDF5

   //const std::string plotFileName = "mesh_data/mesh.h5";
   //HDF5Handle handle( plotFileName.c_str(), HDF5Handle::CREATE );
  
   std::string prefix = "mesh";
   std::string plotFileNameMesh( plotFileName( prefix,
                                               "mesh" ) );
   HDF5Handle handle( plotFileNameMesh.c_str(), HDF5Handle::CREATE );
   
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
   header.m_real["length_scale_SI"] = m_units.getScale(m_units.LENGTH);

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
    
   const LevelData<EdgeDataBox>& gridOnEdges(m_mesh.getXec());
   
   const std::string group3Name = std::string("edge_centered_grid");
   handle.setGroup(group3Name);
   
   header.clear();
   header.m_int["is_edgebox"] = 1;
   header.m_int["num_components"] = gridOnEdges.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, gridOnEdges.boxLayout());
   write(handle, gridOnEdges, "data", gridOnEdges.ghostVect());
   
   //
   // write the node centered grid data
   //
    
   const LevelData<NodeFArrayBox>& gridOnNodes(m_mesh.getXnc());
   
   const std::string group4Name = std::string("node_centered_grid");
   handle.setGroup(group4Name);
   
   header.clear();
   header.m_int["is_nodebox"] = 1;
   header.m_int["num_components"] = gridOnNodes.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, gridOnNodes.boxLayout());
   write(handle, gridOnNodes, "data", gridOnNodes.ghostVect());

   //
   // write node centered data test
   //
   //
   if(m_nodeDataTest) {

      const DisjointBoxLayout& grids = gridOnNodes.disjointBoxLayout();
      LevelData<NodeFArrayBox> nodeTestData(grids,1,gridOnNodes.ghostVect());
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         FArrayBox& this_nodeData = nodeTestData[dit].getFab();
         this_nodeData.setVal(procID());
      }
      //nodeTestData.exchange();   
      SpaceUtils::exchangeNodeFArrayBox(nodeTestData,m_mesh);   

      const std::string group5Name = std::string("node_centered_test");
      handle.setGroup(group5Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = nodeTestData.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, nodeTestData.boxLayout());
      write(handle, nodeTestData, "data", nodeTestData.ghostVect());

   }

   //
   // close the handle
   //

   handle.close();
   
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif

}

void dataFileIO::writeElectroMagneticFieldsDataFile( const ElectroMagneticFields&  a_emfield,
                                                     const int                     a_cur_step, 
                                                     const double                  a_cur_time )
{
   CH_TIME("dataFileIO::writeElectroMagneticFieldsDataFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   if(!procID()) cout << "writing EandM fields data file ..." << endl << endl;
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   
   std::string base_dir = dirPrefix("mesh");
   stringstream s;
   s << base_dir << "field"; 
   std::string prefix = s.str();
   //std::string prefix = "mesh_data/field";
   std::string plotFileNameFields( plotFileName( prefix,
                                                 "fields",
                                                 a_cur_step ) );
   
   HDF5Handle handle( plotFileNameFields.c_str(), HDF5Handle::CREATE );
   
   //
   // write the default header
   //
  
   handle.setGroup("/");
   HDF5HeaderData header;

   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';
   
   for (int dir=0; dir<SpaceDim; dir++)
   {
      sprintf(field_name, "magneticField_%c", coords[dir]);
      vectNames.push_back(field_name);
   }

   for (int i = 0; i < vectNames.size(); ++i)
   {
     sprintf(comp_name, "component_%d", i);
     header.m_string[comp_name] = vectNames[i];
   }
   header.m_int["num_components"] = 1;
   
   header.writeToFile(handle);
  
   //////////////////
  
   const std::string groupName = std::string("field_data");
   handle.setGroup(groupName);  
   
   header.m_int["step_number"] = a_cur_step;
   header.m_real["time"] = a_cur_time;
   header.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
   header.m_real["electric_field_scale_SI"] = m_units.getScale(m_units.ELECTRIC_FIELD);
   header.m_real["magnetic_field_scale_SI"] = m_units.getScale(m_units.MAGNETIC_FIELD);
   header.m_box["prob_domain"] = domain.domainBox();
   
   header.writeToFile(handle);

   //
   // write the magnetic field data
   //
   
   const LevelData<FluxBox>& Bfield  = a_emfield.getMagneticField();

   const std::string group2Name = std::string("magnetic_field");
   handle.setGroup(group2Name);
   
   header.clear();
   header.m_int["is_fluxbox"] = 1;
   header.m_int["num_components"] = Bfield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, Bfield.boxLayout());
   write(handle, Bfield, "data", Bfield.ghostVect());
   
   //
   // write the virtual magnetic field data
   //

   if(SpaceDim<3) {
   
      const LevelData<FArrayBox>& Bfield_virt  = a_emfield.getVirtualMagneticField();

      const std::string group3Name = std::string("virtual_magnetic_field");
      handle.setGroup(group3Name);
   
      header.clear();
      header.m_int["is_cellbox"] = 1;
      header.m_int["num_components"] = Bfield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, Bfield_virt.boxLayout());
      write(handle, Bfield_virt, "data", Bfield_virt.ghostVect());

   }
   
   //
   // write the electric field data
   //
   
   const LevelData<EdgeDataBox>& Efield  = a_emfield.getElectricField();

   const std::string group4Name = std::string("electric_field");
   handle.setGroup(group4Name);
   
   header.clear();
   header.m_int["is_edgebox"] = 1;
   header.m_int["num_components"] = Efield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, Efield.boxLayout());
   write(handle, Efield, "data", Efield.ghostVect());
   
   //
   // write the virtual electric field data
   //

   if(SpaceDim<3) {
   
      const LevelData<NodeFArrayBox>& Efield_virt  = a_emfield.getVirtualElectricField();

      const std::string group5Name = std::string("virtual_electric_field");
      handle.setGroup(group5Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = Efield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, Efield_virt.boxLayout());
      write(handle, Efield_virt, "data", Efield_virt.ghostVect());

   }
   
   //
   // write the current density data
   //
   
   const LevelData<EdgeDataBox>& Jfield  = a_emfield.getCurrentDensity();

   const std::string group6Name = std::string("current_density");
   handle.setGroup(group6Name);
   
   header.clear();
   header.m_int["is_edgebox"] = 1;
   header.m_int["num_components"] = Jfield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(handle);
   
   write(handle, Jfield.boxLayout());
   write(handle, Jfield, "data", Jfield.ghostVect());
   
   //
   // write the virtual current density data
   //

   if(SpaceDim<3) {
   
      const LevelData<NodeFArrayBox>& Jfield_virt  = a_emfield.getVirtualCurrentDensity();

      const std::string group7Name = std::string("virtual_current_density");
      handle.setGroup(group7Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = Jfield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, Jfield_virt.boxLayout());
      write(handle, Jfield_virt, "data", Jfield_virt.ghostVect());

   }
   
   //
   // close the handle
   //

   handle.close();
   
   if(!procID()) cout << "finished writing field data file" << endl << endl;


}      

void dataFileIO::writeNeutralSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                              const int          a_species,
                                              const int          a_cur_step,
                                              const double       a_cur_time )
{
   CH_TIME("dataFileIO::writeNeutralSpeciesDataFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   if(!procID()) {
      cout << "writing neutral species data file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }

   // write the species particle data file
   writeSpeciesParticleFile( a_picSpecies, a_species, a_cur_step, a_cur_time );
   
   // write the species moment data file
   writeSpeciesMomentsFile( a_picSpecies, a_species, a_cur_step, a_cur_time, false, false );
   
   if(!procID()) cout << "finished writing neutral species data file" << endl << endl;

}

void dataFileIO::writeChargedSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                              const int          a_species,
                                              const int          a_cur_step,
                                              const double       a_cur_time,
                                              const bool         a_write_charge_density,
                                              const bool         a_write_current_density )
{
   CH_TIME("dataFileIO::writeChargedSpeciesDataFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   if(!procID()) {
      cout << "writing charged species data file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }

   // write the species particle data file
   writeSpeciesParticleFile( a_picSpecies, a_species, a_cur_step, a_cur_time );
   
   // write the species moment data file
   writeSpeciesMomentsFile( a_picSpecies, a_species, a_cur_step, a_cur_time, 
                            a_write_charge_density, a_write_current_density );
   
   if(!procID()) cout << "finished writing charged species data file" << endl << endl;

}

void dataFileIO::writeSpeciesParticleFile( const PicSpecies&  a_picSpecies,
                                           const int          a_species,
                                           const int          a_cur_step,
                                           const double       a_cur_time )
{
   CH_TIME("dataFileIO::writeSpeciesParticleFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   // get references to particle data   
   const ParticleData<JustinsParticle>& a_Pdata = a_picSpecies.partData();
   
   int verbosity = 0;
   if(!procID() && verbosity) {
      cout << "writing species " << a_species << " particle file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   
   // create the species particle file
   std::string base_dir = dirPrefix("particle");
   stringstream s;
   s << base_dir << "species" << a_species; 
   std::string prefix = s.str();
   std::string plotFileNameParts( plotFileName( prefix,
                                                "parts",
                                                a_cur_step ) );

   HDF5Handle handleParts( plotFileNameParts.c_str(), HDF5Handle::CREATE );
   
   handleParts.setGroup("/");

   // write the header stuff
   HDF5HeaderData headerParts;
   
   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';

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
   //

   const std::string groupName = std::string("species_data");
   handleParts.setGroup(groupName);  

   int totalParticleCount = a_Pdata.numParticles();
   headerParts.m_int["num_particles"] = totalParticleCount;
   headerParts.m_real["mass"] = a_picSpecies.mass();
   headerParts.m_int["charge"] = a_picSpecies.charge();
   headerParts.m_real["Uint"] = a_picSpecies.Uint();
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
   //headerParts.m_real["number_density_scale_SI"] = m_units.getScale(m_units.NUMBER_DENSITY);
   headerParts.m_box["prob_domain"] = domain.domainBox();
   
   // write the header 
   headerParts.writeToFile(handleParts);
 
   // write the particles
   write(handleParts, grids);
   writeParticlesToHDF(handleParts, a_Pdata, "particles", a_picSpecies.writeAll() );

   //
   // close the handle
   //

   handleParts.close();
   
   if(!procID() && verbosity) {
      cout << "finished writing species " << a_species << " particle file" << endl << endl;
   }

}

void dataFileIO::writeSpeciesMomentsFile( const PicSpecies&  a_picSpecies,
                                          const int          a_species,
                                          const int          a_cur_step,
                                          const double       a_cur_time,
                                          const bool         a_write_charge_density,
                                          const bool         a_write_current_density )
{
   CH_TIME("dataFileIO::writeSpeciesMomentsFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   // It is the job of the caller to maker sure the moments are set
   const LevelData<FArrayBox>& density  = a_picSpecies.getNumberDensity();
   const LevelData<FArrayBox>& momentum = a_picSpecies.getMomentumDensity();
   const LevelData<FArrayBox>& energy   = a_picSpecies.getEnergyDensity();
  
   int verbosity = 0; 
   if(!procID() && verbosity) {
      cout << "writing species moments file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   
   // create the species moment file
   std::string base_dir = dirPrefix("mesh");
   stringstream s;
   s << base_dir << "species" << a_species; 
   std::string prefix = s.str();
   std::string plotFileNameParts( plotFileName( prefix,
                                                "moments",
                                                a_cur_step ) );

   HDF5Handle handleParts( plotFileNameParts.c_str(), HDF5Handle::CREATE );
   
   handleParts.setGroup("/");

   //
   // write the header stuff
   //

   HDF5HeaderData headerParts;
   
   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';

   int numMeshComps = density.nComp() + momentum.nComp() + energy.nComp();
   vectNames.push_back("density");
   for (int dir=0; dir<momentum.nComp(); dir++) {
      sprintf(field_name, "momentumDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<energy.nComp(); dir++) {
      sprintf(field_name, "energyDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for(int i=0; i<numMeshComps; i++) {
      sprintf(comp_name,"component_%d", i);
      headerParts.m_string[comp_name] = vectNames[i];
   }
   headerParts.m_int["num_components"] = numMeshComps;

   headerParts.writeToFile(handleParts);

   //
   ///////////////////////////////////// 
   //

   const std::string groupName = std::string("species_data");
   handleParts.setGroup(groupName);  


   headerParts.m_real["mass"] = a_picSpecies.mass();
   headerParts.m_int["charge"] = a_picSpecies.charge();
   headerParts.m_real["Uint"] = a_picSpecies.Uint();
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
   //headerParts.m_real["number_density_scale_SI"] = m_units.getScale(m_units.NUMBER_DENSITY);
   headerParts.m_box["prob_domain"] = domain.domainBox();
   
   // write the header 
   headerParts.writeToFile(handleParts);
 

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
   }
   write(handleParts, momentData.boxLayout());
   write(handleParts, momentData, "data", density.ghostVect());
   
   if(a_write_charge_density) {
      
      //
      // write the cell centered charge density
      //
   
      const LevelData<FArrayBox>& chargeDensity  = a_picSpecies.getChargeDensity();
      
      const std::string groupRho_Name = std::string("cell_centered_charge_density");
      handleParts.setGroup(groupRho_Name);
   
      headerParts.clear();
      headerParts.m_int["is_cellbox"] = 1;
      headerParts.m_int["num_components"] = chargeDensity.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox(); 
      headerParts.writeToFile(handleParts);
   
      write(handleParts, chargeDensity.boxLayout());
      write(handleParts, chargeDensity, "data", chargeDensity.ghostVect());

      //
      // write the face centered charge density
      //
   
      const LevelData<FluxBox>& chargeDenOnFaces  = a_picSpecies.getChargeDensityOnFaces();

      const std::string groupRhoFaces_Name = std::string("face_centered_charge_density");
      handleParts.setGroup(groupRhoFaces_Name);
   
      headerParts.clear();
      headerParts.m_int["is_fluxbox"] = 1;
      headerParts.m_int["num_components"] = chargeDenOnFaces.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox(); 
      headerParts.writeToFile(handleParts);
   
      write(handleParts, chargeDenOnFaces.boxLayout());
      write(handleParts, chargeDenOnFaces, "data", chargeDenOnFaces.ghostVect());

   }

   if(a_write_current_density) {

      //
      // write the edge centered current density
      //

      const LevelData<EdgeDataBox>& JonEdges  = a_picSpecies.getCurrentDensity();
      const std::string groupJ_Name = std::string("current_density");
      handleParts.setGroup(groupJ_Name);
   
      headerParts.clear();
      headerParts.m_int["is_edgebox"] = 1;
      headerParts.m_int["num_components"] = JonEdges.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox(); 
      headerParts.writeToFile(handleParts);
   
      write(handleParts, JonEdges.boxLayout());
      write(handleParts, JonEdges, "data", JonEdges.ghostVect());

      //
      // write the node centered virtual current density for 1D/2D sims
      //
   
      if(SpaceDim<3) {
         const LevelData<NodeFArrayBox>& JvirtOnNodes  = a_picSpecies.getCurrentDensity_virtual();
         const std::string groupJvirtOnNodes_Name = std::string("virtual_current_density");
         handleParts.setGroup(groupJvirtOnNodes_Name);
   
         headerParts.clear();
         headerParts.m_int["is_nodebox"] = 1;
         headerParts.m_int["num_components"] = JvirtOnNodes.nComp();
         headerParts.m_box["prob_domain"] = domain.domainBox(); 
         headerParts.writeToFile(handleParts);
   
         write(handleParts, JvirtOnNodes.boxLayout());
         write(handleParts, JvirtOnNodes, "data", JvirtOnNodes.ghostVect());
      }

   }
   
   //
   // close the handle
   //

   handleParts.close();
   
   if(!procID() && verbosity) cout << "finished writing species moment data file" << endl << endl;

}
/*
void dataFileIO::writeChargedSpeciesDataFile( const PicSpecies&  a_picSpecies,
                                              const int          a_species,
                                              const int          a_cur_step,
                                              const double       a_cur_time,
                                              const bool         a_writeChargeDensity,
                                              const bool         a_writeCurrentDensity )
{
   CH_TIME("dataFileIO::writeChargedSpeciesDataFile()");
   
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   // get references to particle data   
   const ParticleData<JustinsParticle>& a_Pdata = a_picSpecies.partData();
   
   // It is the job of the caller to maker sure the moments are set
   const LevelData<FArrayBox>& density  = a_picSpecies.getNumberDensity();
   const LevelData<FArrayBox>& momentum = a_picSpecies.getMomentumDensity();
   const LevelData<FArrayBox>& energy   = a_picSpecies.getEnergyDensity();
   // 
   //const LevelData<FArrayBox>& chargeDensity  = a_picSpecies.getChargeDensity();
   
   if(!procID()) {
      cout << "writing charged species data file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   
   ///////////// 
   //
   //   write the particles
   //
   /////////////

   stringstream s;
   s << "species" << a_species; 
   std::string prefix = s.str();
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
   numMeshComps = numMeshComps; // + chargeDensity.nComp();
   vectNames.push_back("density");
   for (int dir=0; dir<momentum.nComp(); dir++) {
      sprintf(field_name, "momentumDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<energy.nComp(); dir++) {
      sprintf(field_name, "energyDensity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   //vectNames.push_back("chargeDensity");
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


   const std::string groupName = std::string("species_data");
   handleParts.setGroup(groupName);  


   int totalParticleCount = a_Pdata.numParticles();
   headerParts.m_int["num_particles"] = totalParticleCount;
   headerParts.m_real["mass"] = a_picSpecies.mass();
   headerParts.m_int["charge"] = a_picSpecies.charge();
   headerParts.m_real["Uint"] = a_picSpecies.Uint();
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
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
      //for (int n=0; n<chargeDensity.nComp(); n++) { 
      //   momentData[dit].copy(chargeDensity[dit],n,compData,1);
      //   compData = compData + 1;
      //} 
   }
   write(handleParts, momentData.boxLayout());
   write(handleParts, momentData, "data", density.ghostVect());
   
   if(a_writeChargeDensity) {
      
      //
      // write the cell centered charge density
      //
   
      const LevelData<FArrayBox>& chargeDensity  = a_picSpecies.getChargeDensity();
      
      const std::string groupRho_Name = std::string("cell_centered_charge_density");
      handleParts.setGroup(groupRho_Name);
   
      headerParts.clear();
      headerParts.m_int["is_cellbox"] = 1;
      headerParts.m_int["num_components"] = chargeDensity.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox(); 
      headerParts.writeToFile(handleParts);
   
      write(handleParts, chargeDensity.boxLayout());
      write(handleParts, chargeDensity, "data", chargeDensity.ghostVect());

      //
      // write the face centered charge density
      //
   
      const LevelData<FluxBox>& chargeDenOnFaces  = a_picSpecies.getChargeDensityOnFaces();

      const std::string groupRhoFaces_Name = std::string("face_centered_charge_density");
      handleParts.setGroup(groupRhoFaces_Name);
   
      headerParts.clear();
      headerParts.m_int["is_fluxbox"] = 1;
      headerParts.m_int["num_components"] = chargeDenOnFaces.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox(); 
      headerParts.writeToFile(handleParts);
   
      write(handleParts, chargeDenOnFaces.boxLayout());
      write(handleParts, chargeDenOnFaces, "data", chargeDenOnFaces.ghostVect());

   }

   if(a_writeCurrentDensity) {

      //
      // write the edge centered current density
      //

      const LevelData<EdgeDataBox>& JonEdges  = a_picSpecies.getCurrentDensity();
      const std::string groupJ_Name = std::string("current_density");
      handleParts.setGroup(groupJ_Name);
   
      headerParts.clear();
      headerParts.m_int["is_edgebox"] = 1;
      headerParts.m_int["num_components"] = JonEdges.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox(); 
      headerParts.writeToFile(handleParts);
   
      write(handleParts, JonEdges.boxLayout());
      write(handleParts, JonEdges, "data", JonEdges.ghostVect());

      //
      // write the node centered virtual current density for 1D/2D sims
      //
   
      if(SpaceDim<3) {
         const LevelData<NodeFArrayBox>& JvirtOnNodes  = a_picSpecies.getCurrentDensity_virtual();
         const std::string groupJvirtOnNodes_Name = std::string("virtual_current_density");
         handleParts.setGroup(groupJvirtOnNodes_Name);
   
         headerParts.clear();
         headerParts.m_int["is_nodebox"] = 1;
         headerParts.m_int["num_components"] = JvirtOnNodes.nComp();
         headerParts.m_box["prob_domain"] = domain.domainBox(); 
         headerParts.writeToFile(handleParts);
   
         write(handleParts, JvirtOnNodes.boxLayout());
         write(handleParts, JvirtOnNodes, "data", JvirtOnNodes.ghostVect());
      }

   }

   //
   // close the handle
   //

   handleParts.close();
   
   if(!procID()) cout << "finished writing charged species data file" << endl << endl;

}
*/

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
