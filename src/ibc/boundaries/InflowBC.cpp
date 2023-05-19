#include "InflowBC.H"
#include "TimeFunctionFactory.H"
#include "PicnicConstants.H"
#include "MathUtils.H"

#include "NamespaceHeader.H"

InflowBC::InflowBC( const int&             a_bdry_layout_index,
                    const string&          a_bdry_name,
                    const int&             a_bdry_dir,
                    const Side::LoHiSide&  a_bdry_side,
                    const string&          a_species_name,
                    const Real&            a_species_mass,
                    const DomainGrid&      a_mesh,
                    const CodeUnits&       a_units )
   : m_bdry_layout_index(a_bdry_layout_index),
     m_bdry_name(a_bdry_name),
     m_bdry_dir(a_bdry_dir),
     m_bdry_side(a_bdry_side),
     m_impose_neumann_density(false),
     m_species_name(a_species_name)
{

   parseParameters();

   // set the bdry location value and the normal dX
   const RealVect& dX(a_mesh.getdX());
   const RealVect& Xmin(a_mesh.getXmin());
   const RealVect& Xmax(a_mesh.getXmax());

   m_dX = dX;
   if(m_bdry_side==0) m_X_bdry = Xmin[m_bdry_dir];
   if(m_bdry_side==1) m_X_bdry = Xmax[m_bdry_dir];

   // set the particle weight (NOT YET MODIFED FOR AXISYMMETRIC GEOMS)
   Real volume_scale = a_units.getScale(a_units.VOLUME);
   Real numDen_scale = a_units.getScale(a_units.NUMBER_DENSITY);
   Real cellVolume = volume_scale;
   for (int dir=0; dir<SpaceDim; dir++) cellVolume = cellVolume*dX[dir];
   m_weight = numDen_scale*m_density*cellVolume/(Real)m_discrete_samples; 
   
   // predefine direction dependent mean and thermal velocity / speed of light 
   Real V0 = sqrt(Constants::QE/Constants::ME); // ele thermal speed at 1eV [m/s]
   for(int dir=0; dir<3; dir++) { 
      m_gbeta_VT[dir] = V0*sqrt(m_temperature[dir]/a_species_mass)/Constants::CVAC;
      m_gbeta_U[dir] = m_velocity[dir]/Constants::CVAC;
   }
      

}

void InflowBC::apply( List<JustinsParticle>&  a_inflow_pList,
                const List<JustinsParticle>&  a_valid_pList,
                const FArrayBox&              a_Xec,
                const Box&                    a_bdry_box,
                const Real                    a_time,
                const Real                    a_cnormDt )
{
   CH_TIME("InflowBC::apply()");
      
   int num_samples = m_discrete_samples;

   // get the time function values
   if(m_timeFunction!=NULL) {
      double multiplier;
      m_timeFunction->getValue(multiplier,a_time);
      num_samples = int(m_discrete_samples*multiplier);
   }

   if(m_impose_neumann_density) { // assuming uniform weights for now
         
      const Real X0 = m_X_bdry + m_dX[m_bdry_dir];
      int bdry_count = 0;

      // loop over valid particles and find those next to bdry
      ListIterator<JustinsParticle> lit(a_valid_pList);
      for(lit.begin(); lit.ok(); ++lit) {
         const JustinsParticle& this_particle = lit();
         const RealVect& this_xpold = this_particle.position_old();
         if(this_xpold[m_bdry_dir] <= X0) bdry_count += 1;
      }
    
      // modify number of samples to try and maintain a uniform density
      // between the ghost cell and the first cell inside domain
      int diff = num_samples - bdry_count;
      if(abs(diff)>0.2*num_samples) cout << "JRA: inflow diff = " << diff << endl;
      num_samples = num_samples + diff;

   }


   BoxIterator bit(a_bdry_box); // grid indices for boundary box    
   IntVect ib;         // cell index
   RealVect local_Xlo;
   for(bit.begin(); bit.ok(); ++bit) {

      if(num_samples<=0) continue;

      ib = bit();
      for(int dir=0; dir<SpaceDim; dir++) {
         local_Xlo[dir] = a_Xec.get(ib,dir);
      }
   
      for (int n(0); n<num_samples; n++) {

         // initialize particle velocities by randomly sampling a maxwellian
         RealVect Xpart, Xpart_old;
         Xpart_old[m_bdry_dir] = m_X_bdry + (2.0*m_bdry_side - 1.0)*m_dX[m_bdry_dir]*MathUtils::rand();
         std::array<Real,3> up = {0,0,0};
         
#ifdef RELATIVISTIC_PARTICLES
         for(int dir=0; dir<3; dir++) { 
            up[dir] = m_gbeta_VT[dir]*MathUtils::randn() + m_gbeta_U[dir];
         }

         Real gammap = 1.0;
         gammap += up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
         gammap = sqrt(gammap);

         Xpart[m_bdry_dir] = Xpart_old[m_bdry_dir] + up[m_bdry_dir]/gammap*a_cnormDt;
         if( (m_bdry_side==0 && Xpart[m_bdry_dir]>m_X_bdry) ||
             (m_bdry_side==1 && Xpart[m_bdry_dir]<m_X_bdry) ) {

            for(int dir=0; dir<SpaceDim; dir++) { 
               if(dir!=m_bdry_dir) {
                  Xpart_old[dir] = local_Xlo[dir] + MathUtils::rand()*m_dX[dir];
               }
               Xpart[dir] = Xpart_old[dir] + 0.5*up[dir]/gammap*a_cnormDt;
            }
            JustinsParticle particle(m_weight, Xpart, up);
            particle.setOldPosition(Xpart_old);
            uint64_t ID = 0; // need to declare ID to make smoke tests not fail
            particle.setID(ID);
            a_inflow_pList.append(particle);

         }
#else
         // only go further of the particle makes it into the physical domain
         up[m_bdry_dir] = m_gbeta_VT[m_bdry_dir]*MathUtils::randn() + m_gbeta_U[m_bdry_dir];
         Xpart[m_bdry_dir] = Xpart_old[m_bdry_dir] + up[m_bdry_dir]*a_cnormDt;
         if( (m_bdry_side==0 && Xpart[m_bdry_dir]>m_X_bdry) ||
             (m_bdry_side==1 && Xpart[m_bdry_dir]<m_X_bdry) ) {

            for(int dir=0; dir<3; dir++) { 
               if(dir!=m_bdry_dir) {
                  up[dir] = m_gbeta_VT[dir]*MathUtils::randn() + m_gbeta_U[dir];
               }
            }
            for(int dir=0; dir<SpaceDim; dir++) { 
               if(dir!=m_bdry_dir) {
                  Xpart_old[dir] = local_Xlo[dir] + MathUtils::rand()*m_dX[dir];
               }
               Xpart[dir] = Xpart_old[dir] + 0.5*up[dir]*a_cnormDt;
            }
            JustinsParticle particle(m_weight, Xpart, up);
            particle.setOldPosition(Xpart_old);
            uint64_t ID = 0; // need to declare ID to make smoke tests not fail
            particle.setID(ID);
            a_inflow_pList.append(particle);

         }
#endif

      }

   } // end loop over bdry_grid indices

}

void InflowBC::parseParameters()
{
   
   string prefix0 = "BC." + m_species_name + "." + m_bdry_name + ".inflow";
   ParmParse fpp0( prefix0.c_str() );
 
   fpp0.get( "density", m_density ); // [1/m^3]
   fpp0.get( "discrete_samples", m_discrete_samples );
   fpp0.query( "impose_neumann_density", m_impose_neumann_density );
      
   std::vector<Real> temp_vect(3); // [eV] 
   std::vector<Real> vel_vect(3); // [m/s] 
   fpp0.getarr("temperature",temp_vect,0,3);
   fpp0.getarr("velocity",vel_vect,0,3);
   for (int dir=0; dir<3; ++dir) {
      m_temperature[dir] = temp_vect[dir];
      m_velocity[dir] = vel_vect[dir];
   }
   
   // create time-function for density
   string prefixtf( fpp0.prefix() );
   prefixtf += ".time_function";
   ParmParse tfpp( prefixtf.c_str() );
   if(tfpp.contains("type")) {
      TimeFunctionFactory  timeFactory;
      m_timeFunction = timeFactory.create(tfpp,1);
   }
   else m_timeFunction = RefCountedPtr<TimeFunction>(NULL);

}

void InflowBC::printParameters()
{
   if(!procID()) {
      cout << "Inflow for species " << m_species_name << " on boundary " << m_bdry_name << endl;
      cout << "impose neumann density = " << m_impose_neumann_density << endl;
      cout << "density = " << m_density << endl;
      cout << "particle weight = " << m_weight << endl;
      for (int dir=0; dir<3; ++dir) cout << "velocity(dir=" << dir << ") = " << m_velocity[dir] << endl;
      for (int dir=0; dir<3; ++dir) cout << "temperature(dir=" << dir << ") = " << m_temperature[dir] << endl;
      if(m_timeFunction!=NULL) cout << "time function specified for density" << endl;
   }

}

#include "NamespaceFooter.H"

