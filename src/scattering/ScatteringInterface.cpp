
#include "ScatteringInterface.H"
#include "ParmParse.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"

ScatteringInterface::ScatteringInterface( const PicSpeciesInterface&  a_pic_species_intf )
   : m_verbosity(true),
     m_scatterDt(DBL_MAX),
     m_weight_method(PROBABILISTIC),
     m_model_search_count_max(99),
     m_izn_energy_joules(0.0),
     m_exc_energy_joules(0.0)
{

   ParmParse pp_scatter("scattering");
        
   // look for specified method for scattering particles of different weights
   std::string weight_method;
   pp_scatter.query("weight_method",weight_method);
   if(!weight_method.empty()) {
      if(weight_method=="PROBABILISTIC" || weight_method=="probabilistic") {
         m_weight_method = PROBABILISTIC;
      }
      else if(weight_method=="CONSERVATIVE" || weight_method=="conservative") {
         m_weight_method = CONSERVATIVE;
      }
      else {
         cout << "weight_method = " << weight_method << endl;
         MayDay::Error( "ScatteringInterface: invalid weight_method" );
      }
   }

   // Create all scattering objects
   ScatteringFactory scatteringFactory;
   
   ScatteringPtrVect coulomb_ptr_vect;
   ScatteringPtrVect elastic_ptr_vect;
   ScatteringPtrVect inelastic_ptr_vect;
   
   if(!procID()) cout << "ScatteringInterface: Creating scattering objects..." << endl;

   // check for all Coloumb scattering flag
   ParmParse pp_scatterC("scattering.coulomb");
   bool coulomb_all = false;
   pp_scatterC.query("all",coulomb_all);
   if(coulomb_all) {
      const PicSpeciesPtrVect& pic_species_ptr_vect = a_pic_species_intf.getPtrVect();
      createAllCoulomb( coulomb_ptr_vect, pic_species_ptr_vect, pp_scatterC );
   } 

   //bool more_scattering_models(true);
   pp_scatter.query("model_search_count_max",m_model_search_count_max);
   int model_search_count = 0;
   while(model_search_count < m_model_search_count_max) {
      
      stringstream s;
      s << "scattering." << model_search_count;
      ParmParse pp_scatter( s.str().c_str() );
     
      // Create scattering object and add it to the scattering vector
      if(pp_scatter.contains("model")) {
         ScatteringPtr this_scattering = scatteringFactory.create( pp_scatter, m_weight_method, 1 );
         if(this_scattering->getScatteringType()==COULOMB) {
            coulomb_ptr_vect.push_back(this_scattering);
         }
         else if(this_scattering->getScatteringType()==ELASTIC ||
                 this_scattering->getScatteringType()==CHARGE_EXCHANGE) {
            elastic_ptr_vect.push_back(this_scattering);
         }
         else {
            inelastic_ptr_vect.push_back(this_scattering);
         }
      } 
      model_search_count += 1;
      if(model_search_count==m_model_search_count_max && 
         pp_scatter.contains("model")) m_model_search_count_max += 1;

   }
   m_num_coulomb = coulomb_ptr_vect.size();
   m_num_elastic = elastic_ptr_vect.size();
   m_num_inelastic = inelastic_ptr_vect.size();

   // concatenate the different ptr vects to m_scattering_ptr_vect
   m_scattering_ptr_vect = coulomb_ptr_vect;
   m_scattering_ptr_vect.append( elastic_ptr_vect );
   m_scattering_ptr_vect.append( inelastic_ptr_vect );
        
   if(!procID()) {
      cout << "ScatteringInterface: Finished creating " << m_scattering_ptr_vect.size();
      cout << " scattering objects" << endl;
      cout << "num_coulomb = " << m_num_coulomb << endl;
      cout << "num_elastic = " << m_num_elastic << endl;
      cout << "num_inelastic = " << m_num_inelastic << endl;
      cout << endl;
   }

}
  
void 
ScatteringInterface::initialize( const PicSpeciesInterface&  a_pic_species_intf,
                                 const DomainGrid&           a_mesh,
                                 const std::string&          a_restart_file_name )
{
      
   if(!procID()) cout << "Initializing scattering objects..." << endl << endl;
   
   for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) {
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      this_scattering->initialize( a_pic_species_intf, a_mesh );
   }

   if(!procID()) {
      cout << "Finished initializing " << m_scattering_ptr_vect.size();
      cout << " scattering objects" << endl << endl;
   }

   if(!a_restart_file_name.empty() && m_scattering_ptr_vect.size()>0) {

      HDF5Handle handle( a_restart_file_name, HDF5Handle::OPEN_RDONLY );
      readCheckpoint( handle );
      handle.close();

   }

}
  
void 
ScatteringInterface::applyScattering( PicSpeciesInterface&  a_pic_species_intf,
                                const DomainGrid&           a_mesh,
                                const Real                  a_dt_sec )
{
   CH_TIME("ScatteringInterface::applyScattering()");
   
   // do Coulomb first
   for (int sct=0; sct<m_num_coulomb; sct++) {
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      this_scattering->applyScattering( a_pic_species_intf, a_mesh, a_dt_sec );
   }
   int offset = m_num_coulomb;

   // do all other elastic collisions second
   for (int sct=offset; sct<offset+m_num_elastic; sct++) {
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      this_scattering->applyScattering( a_pic_species_intf, a_mesh, a_dt_sec );
   }
   offset += m_num_elastic;
   
   // do inelastic last
   for (int sct=offset; sct<offset+m_num_inelastic; sct++) {
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      this_scattering->applyScattering( a_pic_species_intf, a_mesh, a_dt_sec );
   }
   CH_assert(offset+m_num_inelastic==m_scattering_ptr_vect.size());

}

void ScatteringInterface::setScatteringProbes()
{  
   CH_TIME("ScatteringInterface::setScatteringProbes()");
   
   Real deltaE_izn_local = 0.0;
   Real deltaE_exc_local = 0.0;
   for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) {
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      if(this_scattering->getScatteringType()==IONIZATION) {
         deltaE_izn_local += this_scattering->getDeltaEizn();
         this_scattering->zeroDeltaEizn();
      }
      if(this_scattering->getScatteringType()==MONTE_CARLO_NULL) {
         deltaE_exc_local += this_scattering->getDeltaEexc();
         this_scattering->zeroDeltaEexc();
         deltaE_izn_local += this_scattering->getDeltaEizn();
         this_scattering->zeroDeltaEizn();
      }
   }

   Real deltaE_izn_global = 0.0;
   Real deltaE_exc_global = 0.0;
#ifdef CH_MPI
   MPI_Allreduce( &deltaE_izn_local,
                  &deltaE_izn_global,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
   MPI_Allreduce( &deltaE_exc_local,
                  &deltaE_exc_global,
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD );
#else
   deltaE_izn_global = deltaE_izn_local;   
   deltaE_exc_global = deltaE_exc_local;   
#endif

  m_izn_energy_joules += deltaE_izn_global;
  m_exc_energy_joules += deltaE_exc_global;

}

Real
ScatteringInterface::scatterDt( const PicSpeciesInterface&  a_pic_species_intf )
{
   CH_TIME("System::scatterDt()");

   Real scatterDt_local = DBL_MAX;
   
   // set mean free time for each scattering object
   for (int sct=0; sct<m_scattering_ptr_vect.size(); sct++) {
      ScatteringPtr this_scattering(m_scattering_ptr_vect[sct]);
      this_scattering->setMeanFreeTime( a_pic_species_intf );
      scatterDt_local = Min(scatterDt_local,this_scattering->scatterDt()); // [s]
   }

   Real m_scatterDt = scatterDt_local;
#ifdef CH_MPI
   MPI_Allreduce( &scatterDt_local, &m_scatterDt, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD ); 
#endif
   
   return m_scatterDt;

}
 

void 
ScatteringInterface::writeCheckpoint( HDF5Handle&  a_handle )
{
   HDF5HeaderData header;
   
   //const std::string& orig_group = a_handle.getGroup();

   //const std::string groupName= std::string("scattering");
   //a_handle.setGroup(groupName);
      
   setScatteringProbes();
   header.m_real["izn_energy_joules"] = m_izn_energy_joules;
   header.m_real["exc_energy_joules"] = m_exc_energy_joules;

   header.writeToFile(a_handle);
   
   //a_handle.setGroup(orig_group);

}

void 
ScatteringInterface::readCheckpoint( HDF5Handle&  a_handle )
{
   HDF5HeaderData header;
   
   //const std::string& orig_group = a_handle.getGroup();
   
   //const std::string groupName= std::string("scattering");
   //a_handle.setGroup(groupName);  

   header.readFromFile( a_handle );
   m_izn_energy_joules = header.m_real["izn_energy_joules"];
   m_exc_energy_joules = header.m_real["exc_energy_joules"];
   
   //a_handle.setGroup(orig_group); // doesnt work for some reason

}

void 
ScatteringInterface::createAllCoulomb( ScatteringPtrVect&  a_coulomb_ptr_vect,
                                 const PicSpeciesPtrVect&  a_pic_species_ptr_vect,
                                 const ParmParse&          a_pp_scatterC ) const
{
         
  if(!procID()) { 
    cout << "ScatteringInterface: Creating Coulomb scattering";
    cout << " objects for all charged species pairs..." << endl;
  }
 
  // query for fixed Coulomb logarithm 
  Real Clog = 0.0;
  if(a_pp_scatterC.contains("coulomb_logarithm")) {
     a_pp_scatterC.get( "coulomb_logarithm", Clog ); 
     CH_assert(Clog>=2.0);
  }
  
  // query the Coulomb collision model type       
  std::string model = "TA";
  a_pp_scatterC.query("model",model);

  // query for the angular scattering type used for Coulomb class
  std::string angular_model = "TAKIZUKA";
  a_pp_scatterC.query("angular_scattering",angular_model);
  
  // query for weight method type used for Coulomb class
  std::string weight_method = "PROBABILISTIC";
  a_pp_scatterC.query("weight_method",weight_method);
  
  if(!procID()) { 
    cout << "ScatteringInterface: Coulomb model type:       " << model << endl;
    cout << "ScatteringInterface: angular scattering model: " << angular_model << endl;
    cout << "ScatteringInterface: weight method:            " << weight_method << endl;
  }

  ScatteringFactory scatteringFactory;

  // create like-like species first
  for (int sp1=0; sp1<a_pic_species_ptr_vect.size(); sp1++) {

    const PicSpeciesPtr this_picSpecies(a_pic_species_ptr_vect[sp1]);
    const bool use_scattering = this_picSpecies->scatter();
    const int charge = this_picSpecies->charge();
    if(!use_scattering || charge==0.0) continue;
            
    // Create scattering object and add it to the scattering vector
    ScatteringPtr this_scattering = scatteringFactory.createCoulomb( sp1, sp1, Clog, model, angular_model, 
		                                                     weight_method, 1 );
    a_coulomb_ptr_vect.push_back(this_scattering);
            
  }
  
  // create unlike species second
  for (int sp1=0; sp1<a_pic_species_ptr_vect.size()-1; sp1++) {

    const PicSpeciesPtr this_picSpecies1(a_pic_species_ptr_vect[sp1]);
    const bool use_scattering1 = this_picSpecies1->scatter();
    const int charge1 = this_picSpecies1->charge();
    if(!use_scattering1 || charge1==0.0) continue;
    
    for (int sp2=sp1+1; sp2<a_pic_species_ptr_vect.size(); sp2++) {
    
      PicSpeciesPtr this_picSpecies2(a_pic_species_ptr_vect[sp2]);
      bool use_scattering2 = this_picSpecies2->scatter();
      int charge2 = this_picSpecies2->charge();
      if(!use_scattering2 || charge2==0.0) continue;
               
      // Create scattering object and add it to the scattering vector
      ScatteringPtr this_scattering = scatteringFactory.createCoulomb( sp1, sp2, Clog, model, angular_model, 
		                                                       weight_method, 1 );
      a_coulomb_ptr_vect.push_back(this_scattering);

    }

  }

}

#include "NamespaceFooter.H"

