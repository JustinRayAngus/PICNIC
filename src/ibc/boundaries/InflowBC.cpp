#include "InflowBC.H"
#include "Constants.H"
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

   // set the particle weight
   Real volume_scale = a_units.getScale(a_units.VOLUME);
   Real numDen_scale = a_units.getScale(a_units.NUMBER_DENSITY);
   Real cellVolume = volume_scale;
   for (int dir=0; dir<SpaceDim; dir++) cellVolume = cellVolume*dX[dir];
   m_weight = numDen_scale*m_density*cellVolume/(Real)m_discrete_samples; 
   
   // predefine direction dependent mean and thermal velocity / speed of light 
   Real V0 = sqrt(Constants::QE/Constants::ME); // ele thermal speed at 1eV [m/s]
   for(int dir=0; dir<3; dir++) { 
      m_beta_VT[dir] = V0*sqrt(m_temperature[dir]/a_species_mass)/Constants::CVAC;
      m_beta_U[dir] = m_velocity[dir]/Constants::CVAC;
   }
      

}

void InflowBC::apply( List<JustinsParticle>&  a_pList,
                const FArrayBox&              a_Xec,
                const Box&                    a_bdry_box,
                const Real&                   a_cnormDt )
{
   CH_TIME("InflowBC::apply()");
      
   BoxIterator bit(a_bdry_box); // grid indices for boundary box    
   IntVect ib;         // cell index
   RealVect local_Xlo;
   for(bit.begin(); bit.ok(); ++bit) {

      ib = bit();
      for(int dir=0; dir<SpaceDim; dir++) {
         local_Xlo[dir] = a_Xec.get(ib,dir);
      }
   
      for (int n(0); n<m_discrete_samples; n++) {

         // initialize particle velocities by randomly sampling a maxwellian
         RealVect Xpart, Xpart_old;
         Xpart_old[m_bdry_dir] = m_X_bdry + (2.0*m_bdry_side - 1.0)*m_dX[m_bdry_dir]*MathUtils::rand();
         std::array<Real,3> BetaPart = {0,0,0};
         BetaPart[m_bdry_dir] = m_beta_VT[m_bdry_dir]*MathUtils::randn() + m_beta_U[m_bdry_dir];

         // only go further of the particle makes it into the physical domain
         Xpart[m_bdry_dir] = Xpart_old[m_bdry_dir] + BetaPart[m_bdry_dir]*a_cnormDt;
         if( (m_bdry_side==0 && Xpart[m_bdry_dir]>m_X_bdry) ||
             (m_bdry_side==1 && Xpart[m_bdry_dir]<m_X_bdry) ) {

            for(int dir=0; dir<3; dir++) { 
               if(dir!=m_bdry_dir) {
                  BetaPart[dir] = m_beta_VT[dir]*MathUtils::randn() + m_beta_U[dir];
               }
            }
            for(int dir=0; dir<SpaceDim; dir++) { 
               if(dir!=m_bdry_dir) {
                  Xpart_old[dir] = local_Xlo[dir] + MathUtils::rand()*m_dX[dir];
                  Xpart[dir] = Xpart_old[dir] + BetaPart[dir]*a_cnormDt;
               }
            }
            JustinsParticle particle(m_weight, Xpart, BetaPart);
            particle.setOldPosition(Xpart_old);
            uint64_t ID = 0; // need to declare ID to make smoke tests not fail
            particle.setID(ID);
            a_pList.append(particle);

         }

      }

   } // end loop over bdry_grid indices

}

void InflowBC::parseParameters()
{
   
   string prefix0 = "BC." + m_species_name + "." + m_bdry_name + ".inflow";
   ParmParse fpp0( prefix0.c_str() );
 
   fpp0.get( "density", m_density ); // [1/m^3]
   fpp0.get( "discrete_samples", m_discrete_samples );
      
   std::vector<Real> temp_vect(3); // [eV] 
   std::vector<Real> vel_vect(3); // [m/s] 
   fpp0.getarr("temperature",temp_vect,0,3);
   fpp0.getarr("velocity",vel_vect,0,3);
   for (int dir=0; dir<3; ++dir) {
      m_temperature[dir] = temp_vect[dir];
      m_velocity[dir] = vel_vect[dir];
   }

}

void InflowBC::printParameters()
{
   if(!procID()) {
      cout << "Inflow for species " << m_species_name << " on boundary " << m_bdry_name << endl;
      cout << "density = " << m_density << endl;
      for (int dir=0; dir<3; ++dir) cout << "velocity(dir=" << dir << ") = " << m_velocity[dir] << endl;
      for (int dir=0; dir<3; ++dir) cout << "temperature(dir=" << dir << ") = " << m_temperature[dir] << endl;
   }

}

#include "NamespaceFooter.H"

