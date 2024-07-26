#include "InsulatorBC.H"
#include "TimeFunctionFactory.H"
#include "InsulatorBCUtils.H"

#include "NamespaceHeader.H"

InsulatorBC::InsulatorBC( const int&         a_bdry_layout_index,
                          const string&      a_bdry_name,
                          const DomainGrid&  a_mesh )
   : m_bdry_layout_index(a_bdry_layout_index)
{

   parseParameters(a_bdry_name);
   defineInsulatorBinaries( a_mesh );

}

void InsulatorBC::defineInsulatorBinaries( const DomainGrid&  a_mesh )
{
   CH_TIME("InsulatorBC::defineInsulatorBinaries()");
 
   const BoundaryBoxLayoutPtrVect& bdry_layout = a_mesh.getBoundaryLayout();
   const DisjointBoxLayout& grids( a_mesh.getDBL() );
   const int num_ghosts = a_mesh.ghosts();

   m_ICbinary_cc.define(grids, 1, IntVect::Unit*num_ghosts);
   m_ICbinary_fc.define(grids, 1, IntVect::Unit*num_ghosts);
   m_ICbinary_ec.define(grids, 1, IntVect::Unit*num_ghosts);
   m_ICbinary_nc.define(grids, 1, IntVect::Unit*num_ghosts);
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      m_ICbinary_cc[dit].setVal(0.0); // initialize to conductor everywhere
      for (int dir=0; dir<SpaceDim; dir++) {
         m_ICbinary_fc[dit][dir].setVal(0.0); // initialize to conductor everywhere
         m_ICbinary_ec[dit][dir].setVal(0.0); // initialize to conductor everywhere
      }
      m_ICbinary_nc[dit].setVal(0.0); // initialize to conductor everywhere
   }
 
   for(int b(0); b<bdry_layout.size(); b++) {
      if(b==m_bdry_layout_index) {
         
         const BoundaryBoxLayout& this_bdry_layout( *(bdry_layout[b]) );
         InsulatorBCUtils::defineInsulatorBinary( m_ICbinary_cc,
                                                  this_bdry_layout,
                                                  m_Xmin_insulator,
                                                  m_Xmax_insulator,
                                                  a_mesh );
            
         InsulatorBCUtils::defineInsulatorBinary( m_ICbinary_fc,
                                                  this_bdry_layout,
                                                  m_Xmin_insulator,
                                                  m_Xmax_insulator,
                                                  a_mesh );
            
         InsulatorBCUtils::defineInsulatorBinary( m_ICbinary_ec,
                                                  this_bdry_layout,
                                                  m_Xmin_insulator,
                                                  m_Xmax_insulator,
                                                  a_mesh );
         
         InsulatorBCUtils::defineInsulatorBinary( m_ICbinary_nc,
                                                  this_bdry_layout,
                                                  m_Xmin_insulator,
                                                  m_Xmax_insulator,
                                                  a_mesh );
      }
   }

}

void
InsulatorBC::parseParameters( const std::string&  a_bdry_name )
{
   
   string prefix0 = "BC.insulator.";
   prefix0 += a_bdry_name;
   ParmParse fpp0( prefix0.c_str() );
   
   fpp0.query("absorbing",m_absorbing);
#if CH_SPACEDIM==1
   fpp0.query("Bv_comp",m_Bv_comp);
#endif

   fpp0.get( "X0_min", m_Xmin_insulator[0] );
   fpp0.get( "X0_max", m_Xmax_insulator[0] );
   if (SpaceDim>1) {        
      fpp0.get( "X1_min", m_Xmin_insulator[1] );
      fpp0.get( "X1_max", m_Xmax_insulator[1] );
   }
   if (SpaceDim>2) {        
      fpp0.get( "X2_min", m_Xmin_insulator[2] );
      fpp0.get( "X2_max", m_Xmax_insulator[2] );
   }

   // create time-function for insulator BC
   string prefixtf( fpp0.prefix() );
   prefixtf += ".time_function";
   ParmParse tfpp( prefixtf.c_str() );
   TimeFunctionFactory  timeFactory;
   m_timeFunction = timeFactory.create(tfpp,1);
         
   if (!procID()) {
      cout << a_bdry_name << " contains insulator : " << endl;
      cout << "  X0_min_insulator = " << m_Xmin_insulator[0] << endl;
      cout << "  X0_max_insulator = " << m_Xmax_insulator[0] << endl;
#if CH_SPACEDIM>1
      cout << "  X1_min_insulator = " << m_Xmin_insulator[1] << endl;
      cout << "  X1_max_insulator = " << m_Xmax_insulator[1] << endl;
#elif CH_SPACEDIM>2
      cout << "  X2_min_insulator = " << m_Xmin_insulator[2] << endl;
      cout << "  X2_max_insulator = " << m_Xmax_insulator[2] << endl;
#endif
      if (m_absorbing) { cout << "  Boundary is an absorbing boundary " << endl; }
      cout << "  Bv_comp = " << m_Bv_comp << endl;
      m_timeFunction->printParameters();
      cout << std::endl;
   }

}

#include "NamespaceFooter.H"

