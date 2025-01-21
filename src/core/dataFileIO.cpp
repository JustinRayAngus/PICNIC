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
     m_mesh(a_mesh),
     m_units(a_units),
     m_verbosity(0)
{

   // parse stuff from input about output directory names
   // and different plotting options for different types of data
   //
  
   m_nodeDataTest = false;
   a_ppsys.query("node_data_test",m_nodeDataTest);
   if (!procID()) { cout << "node_data_test = " << m_nodeDataTest << endl; } 

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
   sprintf( buffer, "%06d.", a_cur_step );
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
   sprintf( buffer, "%06d", a_cur_step );
   std::string filename( dir_prefix + a_prefix + buffer0 + a_species_name + "." + a_diag_name + buffer + "." );

   return filename;
}

void dataFileIO::writeMeshDataFile()
{
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   if (!procID()) { cout << "writing mesh data file ..." << endl << endl; }
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
   if (SpaceDim==2) { coords[1] = 'z'; }
   if (SpaceDim==3) {
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
   
   if (m_mesh.writeJacobians()) {
      
      // write the cell centered corrected jacobian
    
      const LevelData<FArrayBox>& JonCells(m_mesh.getJcc());
   
      const std::string group_Jcc_Name = std::string("cell_centered_jacobian");
      handle.setGroup(group_Jcc_Name);
   
      header.clear();
      header.m_int["num_components"] = JonCells.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JonCells.boxLayout());
      write(handle, JonCells, "data", JonCells.ghostVect());
      
      // write the face centered corrected jacobian
    
      const LevelData<FluxBox>& JonFaces(m_mesh.getJfc());
   
      const std::string group_Jfc_Name = std::string("face_centered_jacobian");
      handle.setGroup(group_Jfc_Name);
   
      header.clear();
      header.m_int["is_fluxbox"] = 1;
      header.m_int["num_components"] = JonFaces.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JonFaces.boxLayout());
      write(handle, JonFaces, "data", JonFaces.ghostVect());
      
      // write the edge centered jacobian
    
      const LevelData<EdgeDataBox>& JonEdges(m_mesh.getJec());
   
      const std::string group_Jec_Name = std::string("edge_centered_jacobian");
      handle.setGroup(group_Jec_Name);
   
      header.clear();
      header.m_int["is_edgebox"] = 1;
      header.m_int["num_components"] = JonEdges.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JonEdges.boxLayout());
      write(handle, JonEdges, "data", JonEdges.ghostVect());

      // write the node centered jacobian
    
      const LevelData<NodeFArrayBox>& JonNodes(m_mesh.getJnc());
   
      const std::string group_Jnc_Name = std::string("node_centered_jacobian");
      handle.setGroup(group_Jnc_Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = JonNodes.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JonNodes.boxLayout());
      write(handle, JonNodes, "data", JonNodes.ghostVect());

   }
   
   if (m_mesh.writeCorrectedJacobians()) {
      
      // write the face centered corrected jacobian
    
      const LevelData<FluxBox>& JConFaces(m_mesh.getCorrectedJfc());
   
      const std::string group_JCfc_Name = std::string("face_centered_corrected_jacobian");
      handle.setGroup(group_JCfc_Name);
   
      header.clear();
      header.m_int["is_fluxbox"] = 1;
      header.m_int["num_components"] = JConFaces.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JConFaces.boxLayout());
      write(handle, JConFaces, "data", JConFaces.ghostVect());
      
      // write the edge centered corrected jacobian
    
      const LevelData<EdgeDataBox>& JConEdges(m_mesh.getCorrectedJec());
   
      const std::string group_JCec_Name = std::string("edge_centered_corrected_jacobian");
      handle.setGroup(group_JCec_Name);
   
      header.clear();
      header.m_int["is_edgebox"] = 1;
      header.m_int["num_components"] = JConEdges.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JConEdges.boxLayout());
      write(handle, JConEdges, "data", JConEdges.ghostVect());

      // write the node centered corrected jacobian
    
      const LevelData<NodeFArrayBox>& JConNodes(m_mesh.getCorrectedJnc());
   
      const std::string group_JCnc_Name = std::string("node_centered_corrected_jacobian");
      handle.setGroup(group_JCnc_Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = JConNodes.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JConNodes.boxLayout());
      write(handle, JConNodes, "data", JConNodes.ghostVect());

   }
   
   if (m_mesh.writeMaskedJacobians()) {
      
      // write the face centered masked jacobian
    
      const LevelData<FluxBox>& JMonFaces(m_mesh.getMaskedJfc());
   
      const std::string group_JMfc_Name = std::string("face_centered_masked_jacobian");
      handle.setGroup(group_JMfc_Name);
   
      header.clear();
      header.m_int["is_fluxbox"] = 1;
      header.m_int["num_components"] = JMonFaces.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JMonFaces.boxLayout());
      write(handle, JMonFaces, "data", JMonFaces.ghostVect());
      
      // write the edge centered masked jacobian
    
      const LevelData<EdgeDataBox>& JMonEdges(m_mesh.getMaskedJec());
   
      const std::string group_JMec_Name = std::string("edge_centered_masked_jacobian");
      handle.setGroup(group_JMec_Name);
   
      header.clear();
      header.m_int["is_edgebox"] = 1;
      header.m_int["num_components"] = JMonEdges.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JMonEdges.boxLayout());
      write(handle, JMonEdges, "data", JMonEdges.ghostVect());

      // write the node centered corrected jacobian
    
      const LevelData<NodeFArrayBox>& JMonNodes(m_mesh.getMaskedJnc());
   
      const std::string group_JMnc_Name = std::string("node_centered_masked_jacobian");
      handle.setGroup(group_JMnc_Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = JMonNodes.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(handle);
   
      write(handle, JMonNodes.boxLayout());
      write(handle, JMonNodes, "data", JMonNodes.ghostVect());

   }

   //
   // write node centered data test
   //

   if (m_nodeDataTest) {

      const DisjointBoxLayout& grids = gridOnNodes.disjointBoxLayout();
      LevelData<NodeFArrayBox> nodeTestData(grids,1,gridOnNodes.ghostVect());
   
      for(DataIterator dit(grids); dit.ok(); ++dit) {
         FArrayBox& this_nodeData = nodeTestData[dit].getFab();
         this_nodeData.setVal(procID());
      }
      //nodeTestData.exchange();   
      SpaceUtils::exchangeNodeFArrayBox(nodeTestData);

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

   // close the handle
   handle.close();
   
#else
   MayDay::Error( "plotting only defined with hdf5" );
#endif

}

void dataFileIO::writeEMFields( HDF5Handle&  a_handle,
                          const EMFields&    a_emfield )
{
   CH_TIME("dataFileIO::writeEMFields()");
   
   const ProblemDomain& domain(m_mesh.getDomain());
   HDF5HeaderData header;

   //
   //  write cummulative boundary diagnostic data
   //

   a_handle.setGroup("field_bdry_data");
   //cout << "JRA: a_handle.getGroup() = " << a_handle.getGroup() << endl;

   const RealVect& intSdAdt_lo = a_emfield.getIntSdAdt_lo();
   const RealVect& intSdAdt_hi = a_emfield.getIntSdAdt_hi();

   header.m_real["intSdAdt_lo0"] = intSdAdt_lo[0];
   header.m_real["intSdAdt_hi0"] = intSdAdt_hi[0];
#if CH_SPACEDIM>1
   header.m_real["intSdAdt_lo1"] = intSdAdt_lo[1];
   header.m_real["intSdAdt_hi1"] = intSdAdt_hi[1];
#endif
#if CH_SPACEDIM==3
   header.m_real["intSdAdt_lo2"] = intSdAdt_lo[2];
   header.m_real["intSdAdt_hi2"] = intSdAdt_hi[2];
#endif
   header.writeToFile( a_handle );

   //
   // write the magnetic field data
   //
   
   const LevelData<FluxBox>& Bfield  = a_emfield.getMagneticField();

   const std::string group2Name = std::string("magnetic_field");
   a_handle.setGroup(group2Name);
   
   header.clear();
   header.m_int["is_fluxbox"] = 1;
   header.m_int["num_components"] = Bfield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(a_handle);
   
   write(a_handle, Bfield.boxLayout());
   write(a_handle, Bfield, "data", Bfield.ghostVect());
   
   //
   // write the virtual magnetic field data
   //

   if (SpaceDim<3) {
   
      const LevelData<FArrayBox>& Bfield_virt  = a_emfield.getVirtualMagneticField();

      const std::string group3Name = std::string("virtual_magnetic_field");
      a_handle.setGroup(group3Name);
   
      header.clear();
      header.m_int["is_cellbox"] = 1;
      header.m_int["num_components"] = Bfield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(a_handle);
   
      write(a_handle, Bfield_virt.boxLayout());
      write(a_handle, Bfield_virt, "data", Bfield_virt.ghostVect());

   }
   
   //
   // write the electric field data
   //
   
   const LevelData<EdgeDataBox>& Efield  = a_emfield.getElectricField();

   const std::string group4Name = std::string("electric_field");
   a_handle.setGroup(group4Name);
   
   header.clear();
   header.m_int["is_edgebox"] = 1;
   header.m_int["num_components"] = Efield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(a_handle);
   
   write(a_handle, Efield.boxLayout());
   write(a_handle, Efield, "data", Efield.ghostVect());
   
   //
   // write the virtual electric field data
   //

   if (SpaceDim<3) {
   
      const LevelData<NodeFArrayBox>& Efield_virt  = a_emfield.getVirtualElectricField();

      const std::string group5Name = std::string("virtual_electric_field");
      a_handle.setGroup(group5Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = Efield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(a_handle);
   
      write(a_handle, Efield_virt.boxLayout());
      write(a_handle, Efield_virt, "data", Efield_virt.ghostVect());

   }

   // set group back to default group
   a_handle.setGroup("/");

}

void dataFileIO::writeEMFields_old( HDF5Handle&  a_handle,
                              const EMFields&    a_emfield )
{
   CH_TIME("dataFileIO::writeEMFields_old()");
   
   const ProblemDomain& domain(m_mesh.getDomain());
   HDF5HeaderData header;
   
   //
   // write the magnetic field data
   //
   
   const LevelData<FluxBox>& Bfield  = a_emfield.getMagneticField_old();

   const std::string group2Name = std::string("magnetic_field_old");
   a_handle.setGroup(group2Name);
   
   header.clear();
   header.m_int["is_fluxbox"] = 1;
   header.m_int["num_components"] = Bfield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(a_handle);
   
   write(a_handle, Bfield.boxLayout());
   write(a_handle, Bfield, "data", Bfield.ghostVect());
   
   //
   // write the virtual magnetic field data
   //

   if (SpaceDim<3) {
   
      const LevelData<FArrayBox>& Bfield_virt  = a_emfield.getVirtualMagneticField_old();

      const std::string group3Name = std::string("virtual_magnetic_field_old");
      a_handle.setGroup(group3Name);
   
      header.clear();
      header.m_int["is_cellbox"] = 1;
      header.m_int["num_components"] = Bfield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(a_handle);
   
      write(a_handle, Bfield_virt.boxLayout());
      write(a_handle, Bfield_virt, "data", Bfield_virt.ghostVect());

   }
   
   //
   // write the electric field data
   //
   
   const LevelData<EdgeDataBox>& Efield  = a_emfield.getElectricField_old();

   const std::string group4Name = std::string("electric_field_old");
   a_handle.setGroup(group4Name);
   
   header.clear();
   header.m_int["is_edgebox"] = 1;
   header.m_int["num_components"] = Efield.nComp();
   header.m_box["prob_domain"] = domain.domainBox(); 
   header.writeToFile(a_handle);
   
   write(a_handle, Efield.boxLayout());
   write(a_handle, Efield, "data", Efield.ghostVect());
   
   //
   // write the virtual electric field data
   //

   if (SpaceDim<3) {
   
      const LevelData<NodeFArrayBox>& Efield_virt  = a_emfield.getVirtualElectricField_old();

      const std::string group5Name = std::string("virtual_electric_field_old");
      a_handle.setGroup(group5Name);
   
      header.clear();
      header.m_int["is_nodebox"] = 1;
      header.m_int["num_components"] = Efield_virt.nComp();
      header.m_box["prob_domain"] = domain.domainBox(); 
      header.writeToFile(a_handle);
   
      write(a_handle, Efield_virt.boxLayout());
      write(a_handle, Efield_virt, "data", Efield_virt.ghostVect());

   }

}

void dataFileIO::writeEMFieldsDataFile( const EMFields&             a_emfield,
#ifdef MASS_MATRIX_TEST
                                        const PicSpeciesInterface&  a_pic_species,
#endif
                                        const int                   a_cur_step, 
                                        const double                a_cur_time )
{
   CH_TIME("dataFileIO::writeEMFieldsDataFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   if (!procID()) { cout << "writing emfields data file" << endl; }
   
   const ProblemDomain& domain(m_mesh.getDomain());
   
   std::string base_dir = dirPrefix("mesh");
   stringstream s;
   s << base_dir << "field"; 
   std::string prefix = s.str();
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
   if (m_mesh.anticyclic()) {
      coords[1] = '2';
      coords[2] = '1';
   }
   
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
   //  write the E and B field variables
   //

   writeEMFields(handle,a_emfield);
   
   //
   // write the current density data
   //
  
   if (a_emfield.advance()) {

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
   
      if (SpaceDim<3) {
      
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
   }

#ifdef MASS_MATRIX_TEST
   writeMassMatricesTest( handle, header, a_pic_species );
#endif
   
   //
   // write the div/curl of the fields
   //

   if (a_emfield.writeRho())   { writeEMFieldRho( a_emfield, handle, header ); }
   if (a_emfield.writeSigma()) { writeEMFieldSigma( a_emfield, handle, header ); }
   if (a_emfield.writeExB())   { writeEMFieldExB( a_emfield, handle, header ); }
   if (a_emfield.writeDivs())  { writeEMFieldDivs( a_emfield, handle, header ); }
   if (a_emfield.writeCurls()) { writeEMFieldCurls( a_emfield, handle, header ); }
   if (a_emfield.usePoisson()) { writeEMFieldPotential( a_emfield, handle, header ); }

   // close the handle
   handle.close();

}

#ifdef MASS_MATRIX_TEST
void dataFileIO::writeMassMatricesTest( HDF5Handle&           a_handle,
                                        HDF5HeaderData&       a_header,
                                  const PicSpeciesInterface&  a_pic_species )
{
   CH_TIME("dataFileIO::writeMassMatricesTest()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write J from mass matriex
   //
   const LevelData<EdgeDataBox>& Jtest = a_pic_species.getJtest();
   
   const std::string group6Name = std::string("current_density_test");
   a_handle.setGroup(group6Name);
      
   a_header.clear();
   a_header.m_int["is_edgebox"] = 1;
   a_header.m_int["num_components"] = Jtest.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
      
   write(a_handle, Jtest.boxLayout());
   write(a_handle, Jtest, "data", Jtest.ghostVect());
   
   //
   //
   //

   const LevelData<NodeFArrayBox>& Jtestv  = a_pic_species.getVirtualJtest();

   const std::string groupName = std::string("current_density_virtual_test");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_nodebox"] = 1;
   a_header.m_int["num_components"] = Jtestv.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, Jtestv.boxLayout());
   write(a_handle, Jtestv, "data", Jtestv.ghostVect());

}
#endif

void dataFileIO::writeEMFieldRho( const EMFields&  a_emfield,
                                  HDF5Handle&      a_handle,
                                  HDF5HeaderData&  a_header )
{
   CH_TIME("dataFileIO::writeEMFieldRho()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write the charge density
   //
   
   const LevelData<NodeFArrayBox>& rho  = a_emfield.getChargeDensity();

   const std::string groupName = std::string("charge_density");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_nodebox"] = 1;
   a_header.m_int["num_components"] = rho.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, rho.boxLayout());
   write(a_handle, rho, "data", rho.ghostVect());

}

void dataFileIO::writeEMFieldSigma( const EMFields&  a_emfield,
                                    HDF5Handle&      a_handle,
                                    HDF5HeaderData&  a_header )
{
   CH_TIME("dataFileIO::writeEMFieldSigma()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write the surface charge
   //
   
   const LevelData<NodeFArrayBox>& sigma  = a_emfield.getSurfaceCharge();

   const std::string groupName = std::string("surface_charge");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_nodebox"] = 1;
   a_header.m_int["num_components"] = sigma.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, sigma.boxLayout());
   write(a_handle, sigma, "data", sigma.ghostVect());

}

void dataFileIO::writeEMFieldExB( const EMFields&  a_emfield,
                                  HDF5Handle&      a_handle,
                                  HDF5HeaderData&  a_header )
{
   CH_TIME("dataFileIO::writeEMFieldExB()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write components of ExB for in-plane E
   //
   
   const LevelData<EdgeDataBox>& exb  = a_emfield.getExB();

   const std::string groupName = std::string("exb");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_edgebox"] = 1;
   a_header.m_int["num_components"] = exb.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, exb.boxLayout());
   write(a_handle, exb, "data", exb.ghostVect());
   
   //
   // write components of ExB for virtual E
   //

   if (SpaceDim<3) {
   
      const LevelData<NodeFArrayBox>& evxb = a_emfield.getEvxB();

      const std::string group2Name = std::string("evxb");
      a_handle.setGroup(group2Name);
   
      a_header.clear();
      a_header.m_int["is_nodebox"] = 1;
      a_header.m_int["num_components"] = evxb.nComp();
      a_header.m_box["prob_domain"] = domain.domainBox(); 
      a_header.writeToFile(a_handle);
   
      write(a_handle, evxb.boxLayout());
      write(a_handle, evxb, "data", evxb.ghostVect());

   }

}

void dataFileIO::writeEMFieldDivs( const EMFields&  a_emfield,
                                   HDF5Handle&      a_handle,
                                   HDF5HeaderData&  a_header )
{
   CH_TIME("dataFileIO::writeEMFieldDivs()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write the div of the electric field
   //
   
   const LevelData<NodeFArrayBox>& divE  = a_emfield.getDivE();

   const std::string groupName = std::string("div_electric_field");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_nodebox"] = 1;
   a_header.m_int["num_components"] = divE.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, divE.boxLayout());
   write(a_handle, divE, "data", divE.ghostVect());

   //
   // write the div of the magnetic field
   //
   
   const LevelData<FArrayBox>& divB  = a_emfield.getDivB();

   const std::string group2Name = std::string("div_magnetic_field");
   a_handle.setGroup(group2Name);
   
   a_header.clear();
   a_header.m_int["is_cellbox"] = 1;
   a_header.m_int["num_components"] = divB.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, divB.boxLayout());
   write(a_handle, divB, "data", divB.ghostVect());
 
}

void dataFileIO::writeEMFieldCurls( const EMFields&  a_emfield,
                                    HDF5Handle&      a_handle,
                                    HDF5HeaderData&  a_header )
{
   CH_TIME("dataFileIO::writeEMFieldCurls()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write the curl of the electric field
   //
   
   const LevelData<FluxBox>& curlE  = a_emfield.getCurlE();

   const std::string groupName = std::string("curl_electric_field");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_fluxbox"] = 1;
   a_header.m_int["num_components"] = curlE.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, curlE.boxLayout());
   write(a_handle, curlE, "data", curlE.ghostVect());

   //
   // write the curl of the virtual electric field
   //

   if (SpaceDim<3) {
   
      const LevelData<FArrayBox>& curlE_virt  = a_emfield.getCurlE_virtual();

      const std::string group2Name = std::string("curl_virtual_electric_field");
      a_handle.setGroup(group2Name);
   
      a_header.clear();
      a_header.m_int["is_cellbox"] = 1;
      a_header.m_int["num_components"] = curlE_virt.nComp();
      a_header.m_box["prob_domain"] = domain.domainBox(); 
      a_header.writeToFile(a_handle);
   
      write(a_handle, curlE_virt.boxLayout());
      write(a_handle, curlE_virt, "data", curlE_virt.ghostVect());

   }
   
   //
   // write the curl of the magnetic field
   //
   
   const LevelData<EdgeDataBox>& curlB  = a_emfield.getCurlB();

   const std::string group3Name = std::string("curl_magnetic_field");
   a_handle.setGroup(group3Name);
   
   a_header.clear();
   a_header.m_int["is_edgebox"] = 1;
   a_header.m_int["num_components"] = curlB.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, curlB.boxLayout());
   write(a_handle, curlB, "data", curlB.ghostVect());

   //
   // write the curl of the virtual magnetic field
   //

   if (SpaceDim<3) {
   
      const LevelData<NodeFArrayBox>& curlB_virt  = a_emfield.getCurlB_virtual();

      const std::string group4Name = std::string("curl_virtual_magnetic_field");
      a_handle.setGroup(group4Name);
   
      a_header.clear();
      a_header.m_int["is_nodebox"] = 1;
      a_header.m_int["num_components"] = curlB_virt.nComp();
      a_header.m_box["prob_domain"] = domain.domainBox(); 
      a_header.writeToFile(a_handle);
   
      write(a_handle, curlB_virt.boxLayout());
      write(a_handle, curlB_virt, "data", curlB_virt.ghostVect());

   }
 
}

void dataFileIO::writeEMFieldPotential( const EMFields&  a_emfield,
                                        HDF5Handle&      a_handle,
                                        HDF5HeaderData&  a_header )
{
   CH_TIME("dataFileIO::writeEMFieldPotential()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   const ProblemDomain& domain(m_mesh.getDomain());
   
   //
   // write the potential
   //
   
   const LevelData<NodeFArrayBox>& potential  = a_emfield.getPotential();

   const std::string groupName = std::string("potential");
   a_handle.setGroup(groupName);
   
   a_header.clear();
   a_header.m_int["is_nodebox"] = 1;
   a_header.m_int["num_components"] = potential.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, potential.boxLayout());
   write(a_handle, potential, "data", potential.ghostVect());
   
   //
   // write the rhs
   //
   
   const LevelData<NodeFArrayBox>& rhs  = a_emfield.getRHS();

   const std::string groupName2 = std::string("rhs");
   a_handle.setGroup(groupName2);
   
   a_header.clear();
   a_header.m_int["is_nodebox"] = 1;
   a_header.m_int["num_components"] = rhs.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, rhs.boxLayout());
   write(a_handle, rhs, "data", rhs.ghostVect());
   
   //
   // write the electric field correction
   //
   
   const LevelData<EdgeDataBox>& E_corr  = a_emfield.getElectricFieldCorrection();

   const std::string groupName3 = std::string("electric_field_correction");
   a_handle.setGroup(groupName3);
   
   a_header.clear();
   a_header.m_int["is_edgebox"] = 1;
   a_header.m_int["num_components"] = E_corr.nComp();
   a_header.m_box["prob_domain"] = domain.domainBox(); 
   a_header.writeToFile(a_handle);
   
   write(a_handle, E_corr.boxLayout());
   write(a_handle, E_corr, "data", E_corr.ghostVect());

}

void dataFileIO::writeSpeciesParticleFile( const PicChargedSpecies&  a_picSpecies,
                                           const int          a_species,
                                           const int          a_cur_step,
                                           const double       a_cur_time )
{
   CH_TIME("dataFileIO::writeSpeciesParticleFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   if (!procID()) {
      cout << "writing particle file for species " << a_species << endl;
   }
   
   // get references to particle data   
   const ParticleData<JustinsParticle>& a_Pdata = a_picSpecies.partData();
   
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

   // write the header stuff
   HDF5HeaderData headerParts;
   
   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';
   if (m_mesh.anticyclic()) {
      coords[1] = '2';
      coords[2] = '1';
   }

   vectNames.clear();
   vectNames.push_back("particle_weight");
   for (int dir=0; dir<SpaceDim; dir++) {
      sprintf(field_name, "particle_position_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   if (SpaceDim<3) {
     for(int i=0; i<3-SpaceDim; i++) {
        sprintf(field_name, "virtual_particle_position_%c", coords[i]);
        vectNames.push_back(field_name);
     }
   }
   for (int dir=0; dir<3; dir++) {
      sprintf(field_name, "particle_velocity_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   if (a_picSpecies.writeAll()) {
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
   headerParts.m_real["mass"]  = a_picSpecies.mass();
   headerParts.m_int["charge"] = a_picSpecies.charge();
   headerParts.m_real["Uint"]  = a_picSpecies.Uint();
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
   //headerParts.m_real["number_density_scale_SI"] = m_units.getScale(m_units.NUMBER_DENSITY);
   headerParts.m_box["prob_domain"] = domain.domainBox();
   
   // write the header 
   headerParts.writeToFile(handleParts);
   
   // write the particles
   write(handleParts, grids);
   writeParticlesToHDF(handleParts, a_Pdata, "particles", a_picSpecies.writeAll() ); // JRA, issue for 3D?
   
   handleParts.close();
   
}

void dataFileIO::writePhotonSpeciesParticleFile( const PicPhotonSpecies&  a_photonSpecies,
                                                 const int          a_cur_step,
                                                 const double       a_cur_time )
{
   CH_TIME("dataFileIO::writePhotonSpeciesParticleFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
   const int species_num = a_photonSpecies.species_num();   

   if (!procID()) {
      cout << "writing particle file for photon species " << species_num << endl;
   }
   
   // get references to particle data   
   const ParticleData<PhotonParticle>& a_Pdata = a_photonSpecies.partData();
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   
   // create the species particle file
   std::string base_dir = dirPrefix("particle");
   stringstream s;
   s << base_dir << "species" << species_num; 
   std::string prefix = s.str();
   std::string plotFileNameParts( plotFileName( prefix,
                                                "parts",
                                                a_cur_step ) );

   HDF5Handle handleParts( plotFileNameParts.c_str(), HDF5Handle::CREATE );

   // write the header stuff
   HDF5HeaderData headerParts;
   
   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';
   if (m_mesh.anticyclic()) {
      coords[1] = '2';
      coords[2] = '1';
   }

   vectNames.clear();
   vectNames.push_back("particle_weight");
   for (int dir=0; dir<SpaceDim; dir++) {
      sprintf(field_name, "particle_position_%c", coords[dir]);
      vectNames.push_back(field_name);
   }
   for (int dir=0; dir<3; dir++) {
      sprintf(field_name, "particle_velocity_%c", coords[dir]);
      vectNames.push_back(field_name);
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
   headerParts.m_real["mass"]  = a_photonSpecies.mass();
   headerParts.m_int["charge"] = a_photonSpecies.charge();
   headerParts.m_int["step_number"] = a_cur_step;
   headerParts.m_real["time"] = a_cur_time;
   headerParts.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
   headerParts.m_box["prob_domain"] = domain.domainBox();
   
   // write the header 
   headerParts.writeToFile(handleParts);
   
   // write the particles
   write(handleParts, grids);
   writeParticlesToHDF(handleParts, a_Pdata, "particles", false );
   
   handleParts.close();
   
}

void dataFileIO::writeSpeciesMomentsFile( const PicChargedSpecies&  a_picSpecies,
                                          const int          a_species,
                                          const int          a_cur_step,
                                          const double       a_cur_time,
                                          const bool         a_write_charge_density,
                                          const bool         a_write_surface_charge,
                                          const bool         a_write_current_density,
                                          const bool         a_write_nppc,
                                          const bool         a_write_energy_off_diag,
                                          const bool         a_write_energy_flux )
{
   CH_TIME("dataFileIO::writeSpeciesMomentsFile()");
   // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

   if (!procID()) {
      cout << "writing moments file for species " << a_species << endl;
   }

   // It is the job of the caller to maker sure the moments are set
   const LevelData<FArrayBox>& density  = a_picSpecies.getNumberDensity();
   const LevelData<FArrayBox>& momentum = a_picSpecies.getMomentumDensity();
   const LevelData<FArrayBox>& energy   = a_picSpecies.getEnergyDensity();

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
   // write header information
   //

   HDF5HeaderData headerParts;

   Vector<string> vectNames;
   char field_name[50];
   char comp_name[50];
   char coords[3];
   coords[0] = '0';
   coords[1] = '1';
   coords[2] = '2';
   if (m_mesh.anticyclic()) {
      coords[1] = '2';
      coords[2] = '1';
   }

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
   for(DataIterator dit(grids); dit.ok(); ++dit) {
      int compData = 0;
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

   if (a_write_charge_density) {

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

      //
      // write the node centered charge density
      //

      const LevelData<NodeFArrayBox>& chargeDenOnNodes  = a_picSpecies.getChargeDensityOnNodes();

      const std::string groupRhoNodes_Name = std::string("node_centered_charge_density");
      handleParts.setGroup(groupRhoNodes_Name);

      headerParts.clear();
      headerParts.m_int["is_nodebox"] = 1;
      headerParts.m_int["num_components"] = chargeDenOnNodes.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox();
      headerParts.writeToFile(handleParts);

      write(handleParts, chargeDenOnNodes.boxLayout());
      write(handleParts, chargeDenOnNodes, "data", chargeDenOnNodes.ghostVect());

   }

   if (a_write_surface_charge) {

      //
      // write the node centered surface charge
      //

      const int ghosts = m_mesh.ghosts();
      LevelData<NodeFArrayBox> surfaceCharge(grids,1,IntVect::Unit*ghosts);
      a_picSpecies.getSurfaceChargeOnNodes( surfaceCharge );

      const std::string group_name = std::string("node_centered_surface_charge");
      handleParts.setGroup(group_name);

      headerParts.clear();
      headerParts.m_int["is_nodebox"] = 1;
      headerParts.m_int["num_components"] = surfaceCharge.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox();
      headerParts.writeToFile(handleParts);

      write(handleParts, surfaceCharge.boxLayout());
      write(handleParts, surfaceCharge, "data", surfaceCharge.ghostVect());

   }

   if (a_write_current_density) {

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

      if (SpaceDim<3) {
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

   if (a_write_nppc) {

      //
      // write the macro particle count in each cell
      //

      const LevelData<FArrayBox>& Nppc  = a_picSpecies.getNppc();

      const std::string groupNppc_Name = std::string("cell_centered_nppc");
      handleParts.setGroup(groupNppc_Name);

      headerParts.clear();
      headerParts.m_int["is_cellbox"] = 1;
      headerParts.m_int["num_components"] = Nppc.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox();
      headerParts.writeToFile(handleParts);

      write(handleParts, Nppc.boxLayout());
      write(handleParts, Nppc, "data", Nppc.ghostVect());

   }

   if (a_write_energy_off_diag) {

      //
      // write the off-diagonal components of the energy density tensor
      //

      const LevelData<FArrayBox>& energyOD  = a_picSpecies.getEnergyOffDiag();

      const std::string groupName = std::string("cell_centered_energy_off_diag");
      handleParts.setGroup(groupName);

      headerParts.clear();
      headerParts.m_int["is_cellbox"] = 1;
      headerParts.m_int["num_components"] = energyOD.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox();
      headerParts.writeToFile(handleParts);

      write(handleParts, energyOD.boxLayout());
      write(handleParts, energyOD, "data", energyOD.ghostVect());

   }

   if (a_write_energy_flux) {

      //
      // write the energy density flux
      //

      const LevelData<FArrayBox>& energyFlux  = a_picSpecies.getEnergyDensityFlux();
      
      const std::string groupName = std::string("cell_centered_energy_flux");
      handleParts.setGroup(groupName);

      headerParts.clear();
      headerParts.m_int["is_cellbox"] = 1;
      headerParts.m_int["num_components"] = energyFlux.nComp();
      headerParts.m_box["prob_domain"] = domain.domainBox();
      headerParts.writeToFile(handleParts);

      write(handleParts, energyFlux.boxLayout());
      write(handleParts, energyFlux, "data", energyFlux.ghostVect());

   }

   // close the handle
   handleParts.close();

}
       
void dataFileIO::writePhotonSpeciesMomentsFile( const PicPhotonSpecies&  a_photonSpecies,
                                                const int          a_cur_step, 
                                                const double       a_cur_time )
{
    CH_TIME("dataFileIO::writePhotonSpeciesMomentFile()");
    // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"
   
    const int species_num = a_photonSpecies.species_num();

    if (!procID()) {cout << "writing moments file for photon species " << species_num << endl;}

    // It is the job of the caller to maker sure the moments are set
    const LevelData<FArrayBox>& density  = a_photonSpecies.getNumberDensity();
    const LevelData<FArrayBox>& momentum = a_photonSpecies.getMomentumDensity();
    const LevelData<FArrayBox>& energy   = a_photonSpecies.getEnergyDensity();

    const DisjointBoxLayout& grids(m_mesh.getDBL());
    const ProblemDomain& domain(m_mesh.getDomain());

    // create the species moment file
    std::string base_dir = dirPrefix("mesh");
    stringstream s;
    s << base_dir << "species" << species_num;
    std::string prefix = s.str();
    std::string plotFileNameParts( plotFileName( prefix,
                                                 "moments",
                                                 a_cur_step ) );

    HDF5Handle handleParts( plotFileNameParts.c_str(), HDF5Handle::CREATE );

    handleParts.setGroup("/");

    //
    // write header information
    //

    HDF5HeaderData headerParts;

    Vector<string> vectNames;
    char field_name[50];
    char comp_name[50];
    char coords[3];
    coords[0] = '0';
    coords[1] = '1';
    coords[2] = '2';
    if (m_mesh.anticyclic()) {
       coords[1] = '2';
       coords[2] = '1';
    }

    int numMeshComps = density.nComp() + momentum.nComp() + energy.nComp();
    vectNames.push_back("density");
    for (int dir=0; dir<momentum.nComp(); dir++) {
        sprintf(field_name, "momentumDensity_%c", coords[dir]);
        vectNames.push_back(field_name);
    }
    vectNames.push_back("energyDensity");
    for (int i=0; i<numMeshComps; i++) {
        sprintf(comp_name,"component_%d", i);
        headerParts.m_string[comp_name] = vectNames[i];
    }
    headerParts.m_int["num_components"] = numMeshComps;
    headerParts.writeToFile(handleParts);

    const std::string groupName = std::string("species_data");
    handleParts.setGroup(groupName);

    headerParts.m_real["mass"] = a_photonSpecies.mass();
    headerParts.m_int["charge"] = a_photonSpecies.charge();
    headerParts.m_int["step_number"] = a_cur_step;
    headerParts.m_real["time"] = a_cur_time;
    headerParts.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
    headerParts.m_box["prob_domain"] = domain.domainBox();

    // write the header 
    headerParts.writeToFile(handleParts);

    // write the moment data
    LevelData<FArrayBox> momentData;
    momentData.define(grids,numMeshComps,density.ghostVect());
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        int compData = 0;
        momentData[dit].copy(density[dit],0,compData,1);
        compData = compData + 1;
        for (int dir=0; dir<momentum.nComp(); dir++) {
            momentData[dit].copy(momentum[dit],dir,compData,1);
            compData = compData + 1;
        }
        momentData[dit].copy(energy[dit],0,compData,1);
        compData = compData + 1;
    }
    write(handleParts, momentData.boxLayout());
    write(handleParts, momentData, "data", density.ghostVect());

    //
    // write the macro particle count in each cell
    //

    const LevelData<FArrayBox>& Nppc  = a_photonSpecies.getNppc();

    const std::string groupNppc_Name = std::string("cell_centered_nppc");
    handleParts.setGroup(groupNppc_Name);

    headerParts.clear();
    headerParts.m_int["is_cellbox"] = 1;
    headerParts.m_int["num_components"] = Nppc.nComp();
    headerParts.m_box["prob_domain"] = domain.domainBox();
    headerParts.writeToFile(handleParts);

    write(handleParts, Nppc.boxLayout());
    write(handleParts, Nppc, "data", Nppc.ghostVect());

    // close the handle
    handleParts.close();

}

void dataFileIO::writeFusionProductsDataFile( const ScatteringPtrVect&  a_scattering_ptr_vect,
                                              const int                 a_num_fusion,
                                              const int                 a_cur_step,
                                              const double              a_cur_time )
{
    CH_TIME("dataFileIO::writeFusionProductsDataFile()");
    // See Chombo_3.2/lib/src/BoxTools/CH_HDF5.H for "write"

    if (!procID()) {
        cout << "writing fusion diagnostic file for fusion products" << std::endl;
    }

    // create the species moment file
    std::string base_dir = dirPrefix("mesh");
    stringstream s;
    s << base_dir << "fusion";
    std::string prefix = s.str();
    std::string plotFileNameParts( plotFileName( prefix,
                                                "products",
                                                a_cur_step ) );

    HDF5Handle handle( plotFileNameParts.c_str(), HDF5Handle::CREATE );
    handle.setGroup("/");

    HDF5HeaderData header;

    header.m_int["num_fusion"] = a_num_fusion;
    header.m_int["step_number"] = a_cur_step;
    header.m_real["time"] = a_cur_time;
    header.m_real["time_scale_SI"] = m_units.getScale(m_units.TIME);
    header.writeToFile(handle);

    writeFusionProducts( handle, a_scattering_ptr_vect );

    // close the handle
    handle.close();

}

void dataFileIO::writeFusionProducts( HDF5Handle&         a_handle,
                                const ScatteringPtrVect&  a_scattering_ptr_vect )
{
    CH_TIME("dataFileIO::writeFusionProducts()");

    const ProblemDomain& domain(m_mesh.getDomain());
    HDF5HeaderData header;

    int fusion_num = -1;
    for (int sct=0; sct<a_scattering_ptr_vect.size(); sct++) {
        ScatteringPtr this_scattering(a_scattering_ptr_vect[sct]);
        if (this_scattering->getScatteringType()==FUSION) {

            // Create the group name for this fusion process. DO NOT CHANGE!
            // The same naming convention is used to read in from the checkpoint
            // file on a restart. See ScatteringInterface::readCheckpoint();
            fusion_num += 1;
            stringstream ssfn;
            ssfn << "fusion_data_" << fusion_num;
            const std::string groupName = ssfn.str();
            a_handle.setGroup(groupName);

            const LevelData<FArrayBox>& fusionProducts = this_scattering->getFusionProducts();
            const int num_comps = fusionProducts.nComp();
            const std::string fusion_type = this_scattering->getScatteringSubTypeName();

            header.clear();
            header.m_int["sct_num"]  = sct;
            header.m_string["fusion_type"] = fusion_type;
            header.m_int["species1"] = this_scattering->species1();
            header.m_int["species2"] = this_scattering->species2();
            header.m_real["fusionQ"] = this_scattering->getFusionQ();
            if (num_comps>1) { header.m_real["fusionQb"] = this_scattering->getFusionQb(); }

            header.m_int["is_cellbox"] = 1;
            header.m_int["num_components"] = num_comps;
            header.m_box["prob_domain"] = domain.domainBox();
            header.writeToFile(a_handle);

            write(a_handle, fusionProducts.boxLayout());
            write(a_handle, fusionProducts, "data", fusionProducts.ghostVect());

        }
    }

    // set group back to default group
    a_handle.setGroup("/");

}

void dataFileIO::writeBinFabDataFile( const LevelData<BinFab<JustinsParticle>>&  a_Pdata,
                                      const int                      a_cur_step,
                                      const double                   a_cur_time )
{
   CH_TIME("dataFileIO::writeBinFabDataFile()");
   if (!procID()) {
      cout << "writing BinFab particle data file" << endl;
      cout << "at step = " << a_cur_step << " and time = " << a_cur_time << endl;
      cout << "... " << endl;
   }

   ///////////// 
   //
   //   write the particles
   //
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

   if (!procID()) { cout << "finished writing binfab particle data file" << endl << endl; }

}

void dataFileIO::writeCheckpointEMFields( HDF5Handle&  a_handle,
                                    const int          a_write_old_data,
                                    const EMFields&    a_emfield )
{
   CH_TIME("dataFileIO::writeCheckpointEMFields()");

   writeEMFields( a_handle, a_emfield );
   if ( a_write_old_data > 0 ) { // needed for time schemes with fields staggered in time
      writeEMFields_old( a_handle, a_emfield );
   }

}

void dataFileIO::writeCheckpointFusion( HDF5Handle&         a_handle,
                                  const ScatteringPtrVect&  a_scattering_ptr_vect )
{
   CH_TIME("dataFileIO::writeCheckpointFusion()");
   writeFusionProducts( a_handle, a_scattering_ptr_vect );
}

void dataFileIO::writeCheckpointSpecies( HDF5Handle&           a_handle,
                                   const PicSpeciesInterface&  a_pic_species )
{
   CH_TIME("dataFileIO::writeCheckpointSpecies()");

   const PicChargedSpeciesPtrVect& pic_species_ptr_vect = a_pic_species.getChargedPtrVect();
   for (int sp=0; sp<pic_species_ptr_vect.size(); sp++) {
      const PicChargedSpeciesPtr species = pic_species_ptr_vect[sp];
      writeCheckpointParticles( a_handle, *species, species->species_num() );
   }

   const PicPhotonSpeciesPtrVect& photon_species_ptr_vect = a_pic_species.getPhotonPtrVect();
   for (int sp=0; sp<photon_species_ptr_vect.size(); sp++) {
      const PicPhotonSpeciesPtr species = photon_species_ptr_vect[sp];
      writeCheckpointPhotonParticles( a_handle, *species );
   }

}

void dataFileIO::writeCheckpointParticles( HDF5Handle&  a_handle,
                                     const PicChargedSpecies&  a_picSpecies,
                                     const int          a_species )
{
   CH_TIME("dataFileIO::writeCheckpointParticles()");

   // get references to particle data   
   const ParticleData<JustinsParticle>& a_Pdata = a_picSpecies.partData();

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());

   HDF5HeaderData headerParts;

   stringstream s;
   s << "species" << a_species;
   const std::string groupName = std::string( s.str() + "_data");
   a_handle.setGroup(groupName);

   int totalParticleCount = a_Pdata.numParticles();
   headerParts.m_int["species_number"] = a_species;
   headerParts.m_int["num_particles"]  = totalParticleCount;
   headerParts.m_real["mass"]  = a_picSpecies.mass();
   headerParts.m_int["charge"] = a_picSpecies.charge();
   headerParts.m_real["Uint"]  = a_picSpecies.Uint();
   headerParts.m_box["prob_domain"] = domain.domainBox();

   const RealVect& massOut_lo = a_picSpecies.getMassOut_lo();
   const RealVect& massOut_hi = a_picSpecies.getMassOut_hi();
   const RealVect& momXOut_lo = a_picSpecies.getMomXOut_lo();
   const RealVect& momXOut_hi = a_picSpecies.getMomXOut_hi();
   const RealVect& momYOut_lo = a_picSpecies.getMomYOut_lo();
   const RealVect& momYOut_hi = a_picSpecies.getMomYOut_hi();
   const RealVect& momZOut_lo = a_picSpecies.getMomZOut_lo();
   const RealVect& momZOut_hi = a_picSpecies.getMomZOut_hi();
   const RealVect& energyOut_lo = a_picSpecies.getEnergyOut_lo();
   const RealVect& energyOut_hi = a_picSpecies.getEnergyOut_hi();

   const RealVect& massIn_lo = a_picSpecies.getMassIn_lo();
   const RealVect& massIn_hi = a_picSpecies.getMassIn_hi();
   const RealVect& momXIn_lo = a_picSpecies.getMomXIn_lo();
   const RealVect& momXIn_hi = a_picSpecies.getMomXIn_hi();
   const RealVect& momYIn_lo = a_picSpecies.getMomYIn_lo();
   const RealVect& momYIn_hi = a_picSpecies.getMomYIn_hi();
   const RealVect& momZIn_lo = a_picSpecies.getMomZIn_lo();
   const RealVect& momZIn_hi = a_picSpecies.getMomZIn_hi();
   const RealVect& energyIn_lo = a_picSpecies.getEnergyIn_lo();
   const RealVect& energyIn_hi = a_picSpecies.getEnergyIn_hi();

   for (int dir=0; dir<SpaceDim; dir++) {
      if (domain.isPeriodic(dir)) { continue; }
      std::string dirstr = std::to_string(dir);
      headerParts.m_real["massOut_lo"+dirstr] = massOut_lo[dir];
      headerParts.m_real["massOut_hi"+dirstr] = massOut_hi[dir];
      headerParts.m_real["momXOut_lo"+dirstr] = momXOut_lo[dir];
      headerParts.m_real["momXOut_hi"+dirstr] = momXOut_hi[dir];
      headerParts.m_real["momYOut_lo"+dirstr] = momYOut_lo[dir];
      headerParts.m_real["momYOut_hi"+dirstr] = momYOut_hi[dir];
      headerParts.m_real["momZOut_lo"+dirstr] = momZOut_lo[dir];
      headerParts.m_real["momZOut_hi"+dirstr] = momZOut_hi[dir];
      headerParts.m_real["energyOut_lo"+dirstr] = energyOut_lo[dir];
      headerParts.m_real["energyOut_hi"+dirstr] = energyOut_hi[dir];

      headerParts.m_real["massIn_lo"+dirstr] = massIn_lo[dir];
      headerParts.m_real["massIn_hi"+dirstr] = massIn_hi[dir];
      headerParts.m_real["momXIn_lo"+dirstr] = momXIn_lo[dir];
      headerParts.m_real["momXIn_hi"+dirstr] = momXIn_hi[dir];
      headerParts.m_real["momYIn_lo"+dirstr] = momYIn_lo[dir];
      headerParts.m_real["momYIn_hi"+dirstr] = momYIn_hi[dir];
      headerParts.m_real["momZIn_lo"+dirstr] = momZIn_lo[dir];
      headerParts.m_real["momZIn_hi"+dirstr] = momZIn_hi[dir];
      headerParts.m_real["energyIn_lo"+dirstr] = energyIn_lo[dir];
      headerParts.m_real["energyIn_hi"+dirstr] = energyIn_hi[dir];
   }

   if (a_picSpecies.charge()!=0) {
     // write the cumulative surface charge diagnostic
#if CH_SPACEDIM>1
     a_picSpecies.preWriteSurfaceCharge();
#endif
     const LevelData<NodeFArrayBox>& surfaceCharge  = a_picSpecies.getSurfaceChargeOnNodes();
     write(a_handle, surfaceCharge, "surface_charge", surfaceCharge.ghostVect());
   }

   // write the header 
   headerParts.writeToFile(a_handle);

   // write the particles
   write(a_handle, grids);
   const bool writeAll = true;
   writeParticlesToHDF(a_handle, a_Pdata, "particles", writeAll );

   // set group back to default group
   a_handle.setGroup("/");

}

void dataFileIO::writeCheckpointPhotonParticles( HDF5Handle&     a_handle,
                                           const PicPhotonSpecies&  a_photonSpecies )
{
   CH_TIME("dataFileIO::writeCheckpointPhotonParticles()");

   // get references to particle data   
   const ParticleData<PhotonParticle>& a_Pdata = a_photonSpecies.partData();
   const int species_num = a_photonSpecies.species_num();

   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());

   HDF5HeaderData headerParts;

   stringstream s;
   s << "species" << species_num;
   const std::string groupName = std::string( s.str() + "_data");
   a_handle.setGroup(groupName);

   int totalParticleCount = a_Pdata.numParticles();
   headerParts.m_int["species_number"] = species_num;
   headerParts.m_int["num_particles"]  = totalParticleCount;
   headerParts.m_real["mass"]  = a_photonSpecies.mass();
   headerParts.m_int["charge"] = a_photonSpecies.charge();
   headerParts.m_box["prob_domain"] = domain.domainBox();

   const RealVect& massOut_lo = a_photonSpecies.getMassOut_lo();
   const RealVect& massOut_hi = a_photonSpecies.getMassOut_hi();
   const RealVect& momXOut_lo = a_photonSpecies.getMomXOut_lo();
   const RealVect& momXOut_hi = a_photonSpecies.getMomXOut_hi();
   const RealVect& momYOut_lo = a_photonSpecies.getMomYOut_lo();
   const RealVect& momYOut_hi = a_photonSpecies.getMomYOut_hi();
   const RealVect& momZOut_lo = a_photonSpecies.getMomZOut_lo();
   const RealVect& momZOut_hi = a_photonSpecies.getMomZOut_hi();
   const RealVect& energyOut_lo = a_photonSpecies.getEnergyOut_lo();
   const RealVect& energyOut_hi = a_photonSpecies.getEnergyOut_hi();

   for (int dir=0; dir<SpaceDim; dir++) {
      if (domain.isPeriodic(dir)) { continue; }
      std::string dirstr = std::to_string(dir);
      headerParts.m_real["massOut_lo"+dirstr] = massOut_lo[dir];
      headerParts.m_real["massOut_hi"+dirstr] = massOut_hi[dir];
      headerParts.m_real["momXOut_lo"+dirstr] = momXOut_lo[dir];
      headerParts.m_real["momXOut_hi"+dirstr] = momXOut_hi[dir];
      headerParts.m_real["momYOut_lo"+dirstr] = momYOut_lo[dir];
      headerParts.m_real["momYOut_hi"+dirstr] = momYOut_hi[dir];
      headerParts.m_real["momZOut_lo"+dirstr] = momZOut_lo[dir];
      headerParts.m_real["momZOut_hi"+dirstr] = momZOut_hi[dir];
      headerParts.m_real["energyOut_lo"+dirstr] = energyOut_lo[dir];
      headerParts.m_real["energyOut_hi"+dirstr] = energyOut_hi[dir];
   }

   // write the header 
   headerParts.writeToFile(a_handle);

   // write the particles
   write(a_handle, grids);
   writeParticlesToHDF(a_handle, a_Pdata, "particles", false );

   // set group back to default group
   a_handle.setGroup("/");

}



#include "NamespaceFooter.H"
