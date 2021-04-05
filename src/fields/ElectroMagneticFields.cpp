
#include <array>
#include <cmath>
#include "Constants.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "GridFunctionFactory.H"

#include "ElectroMagneticFields.H"

#include "MathUtils.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

ElectroMagneticFields::ElectroMagneticFields( ParmParse&   a_ppflds,
                                        const DomainGrid&  a_mesh,
                                        const bool&        a_verbosity )
   : m_verbosity(a_verbosity),
     m_mesh(a_mesh)
{
   // parse input file
   a_ppflds.get( "advance", m_advance );

   if ( procID() == 0 ) {
      cout << " advance fields = " << m_advance << endl << endl;
   }
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& meshSpacing(m_mesh.getdX());
   const RealVect& meshOrigin(m_mesh.getXmin());
   const int ghosts(m_mesh.ghosts());

   m_stable_dt = meshSpacing[0]/Constants::CVAC;
   for (int dir=1; dir<SpaceDim; dir++) {
      m_stable_dt = Min(m_stable_dt,meshSpacing[dir]/Constants::CVAC);
   } 

   //
   // initialize the member LevelDatas
   //
   //
   const IntVect ghostVect = ghosts*IntVect::Unit; 
   
   // initialize the magnetic field containers
   m_magneticField.define(grids,1,ghostVect);     // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_magneticField_virtual.define(grids,3-SpaceDim,ghostVect); // FArrayBox
   }

   // initialize the electric field containers
   m_electricField.define(grids,1,ghostVect);     // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_electricField_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   }
   
   // initialize the current density containers
   m_currentDensity.define(grids,1,ghostVect);    // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_currentDensity_virtual.define(grids,3-SpaceDim,ghostVect); // NodeFArrayBox
   }
   
   // initialize the curlE containers
   m_curlE.define(grids,1,IntVect::Zero);     // FluxBox with 1 comp
   if(SpaceDim<3) {
      m_curlE_virtual.define(grids,3-SpaceDim,IntVect::Zero); // FArrayBox
   }
   
   // initialize the curlB containers
   m_curlB.define(grids,1,IntVect::Zero);     // EdgeDataBox with 1 comp
   if(SpaceDim<3) {
      m_curlB_virtual.define(grids,3-SpaceDim,IntVect::Zero); // NodeFArrayBox
   }

}

ElectroMagneticFields::~ElectroMagneticFields()
{
}

void ElectroMagneticFields::initialize()
{
   if(!procID()) cout << "Initializing electromagnetic fields..." << endl;
   
   //
   // parse the initial profiles for the fields
   //
   
   GridFunctionFactory  gridFactory;
   const Real this_time = 0.0;
   
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   LevelData<FArrayBox> temporaryProfile;
   temporaryProfile.define(grids,1,m_magneticField_virtual.ghostVect());
   LevelData<NodeFArrayBox> temporaryNodeProfile;
   temporaryNodeProfile.define(grids,1,m_electricField_virtual.ghostVect());
   
   // set initial magnetic field profile from ICs
   for (int dir=0; dir<SpaceDim; dir++) {
      stringstream sB; 
      sB << "IC.em_fields.magnetic." << dir; 
      ParmParse ppBfieldIC( sB.str().c_str() );
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppBfieldIC,m_verbosity);
      gridFunction->assign( m_magneticField, dir, m_mesh, this_time );
   }
   
   // set initial virtual magnetic field profile from ICs
   if(SpaceDim<3) {
      int numVirt = 3-SpaceDim;
      for( int comp=0; comp<numVirt; comp++) {
         int field_dir = comp + SpaceDim;
         stringstream sBv; 
         sBv << "IC.em_fields.magnetic." << field_dir; 
         ParmParse ppBvfieldIC( sBv.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppBvfieldIC,m_verbosity);
         gridFunction->assign( temporaryProfile, m_mesh, this_time );
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            m_magneticField_virtual[dit].copy(temporaryProfile[dit],0,comp,1);
         }
      }
   }
   
   // set initial electric field profile from ICs
   for (int dir=0; dir<SpaceDim; dir++) {
      stringstream sE; 
      sE << "IC.em_fields.electric." << dir; 
      ParmParse ppEfieldIC( sE.str().c_str() );
      RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppEfieldIC,m_verbosity);
      gridFunction->assign( m_electricField, dir, m_mesh, this_time );
   }
   
   // set initial virtual electric field profile from ICs
   if(SpaceDim<3) {
      int numVirt = 3-SpaceDim;
      for( int comp=0; comp<numVirt; comp++) {
         int field_dir = comp + SpaceDim;
         stringstream sEv; 
         sEv << "IC.em_fields.electric." << field_dir; 
         ParmParse ppEvfieldIC( sEv.str().c_str() );
         RefCountedPtr<GridFunction> gridFunction = gridFactory.create(ppEvfieldIC,m_verbosity);
         gridFunction->assign( temporaryNodeProfile, m_mesh, this_time );
         for(DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& this_Ev  = m_electricField_virtual[dit].getFab();
            this_Ev.copy(temporaryNodeProfile[dit].getFab(),0,comp,1);
         }
      }
   }

   if(!procID()) cout << "Finished initializing electromagnetic fields " << endl << endl;

}

void ElectroMagneticFields::advanceMagneticField( const Real& a_dt )
{
   CH_TIME("ElectroMagneticFields::advanceMagneticField()");
    
   // get some Grid info
   //
   const DisjointBoxLayout& grids(m_mesh.getDBL());
   const ProblemDomain& domain(m_mesh.getDomain());
   const RealVect& dX(m_mesh.getdX());
   //const int ghosts(m_mesh.ghosts());
   
   // compute domain extent and get total num parts and sim volume
   //
   const IntVect domainDimensions = domain.size();
   const RealVect& Xmin(m_mesh.getXmin());
   const RealVect& Xmax(m_mesh.getXmax());
   const RealVect Lbox = Xmax - Xmin; 
   
}

void ElectroMagneticFields::advanceElectricField( const Real& a_dt )
{
   CH_TIME("ElectroMagneticFields::advanceElectricField()");
    
}


#include "NamespaceFooter.H"

