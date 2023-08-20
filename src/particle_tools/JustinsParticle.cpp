#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "JustinsParticle.H"
#include "NamespaceHeader.H"

/// default constructor
JustinsParticle::JustinsParticle() :
  BinItem()
{}

JustinsParticle::~JustinsParticle()
{}

JustinsParticle::JustinsParticle( const Real                 a_weight,
                                  const RealVect&            a_position,
                                  const std::array<Real,3>&  a_velocity )
  :
  BinItem(a_position),
  m_weight(a_weight),
  m_velocity(a_velocity)
{
  setOldPosition(a_position);
  setOldVelocity(a_velocity);
  //std::array<Real,3> zero_array = {0,0,0};
  //setElectricField(zero_array);
  //setMagneticField(zero_array);
#if CH_SPACEDIM<3
  std::array<Real,2> thisPosVirt;
  for(int i=0; i<2; i++) thisPosVirt[i] = 0.0 + i*0.0;
  setPositionVirt(thisPosVirt);
#endif
  setID();
}

void JustinsParticle::define( const Real                 a_weight,
                              const RealVect&            a_position,
                              const std::array<Real,3>&  a_velocity )
{
  setID();
  setWeight(a_weight);
  setPosition(a_position);
  setOldPosition(a_position);
  setVelocity(a_velocity);
  setOldVelocity(a_velocity);
  //std::array<Real,3> zero_array = {0,0,0};
  //setElectricField(zero_array);
  //setMagneticField(zero_array);
#if CH_SPACEDIM<3
  std::array<Real,2> thisPosVirt;
  for(int i=0; i<2; i++) thisPosVirt[i] = 0.0 + i*0.0;
  setPositionVirt(thisPosVirt);
#endif
}

void JustinsParticle::setKillTag()
{
  m_kill_tag = 1;
}
  
const int& JustinsParticle::killTag() const
{
  return m_kill_tag;
}

void JustinsParticle::setNumSubOrbits( const int a_num_suborbits )
{
  m_num_suborbits = a_num_suborbits;
}

const int& JustinsParticle::numSubOrbits() const
{
  return m_num_suborbits;
}

void JustinsParticle::setID(const uint64_t a_ID)
//void JustinsParticle::setID(const Real a_ID)
{
  m_ID = a_ID;
}

void JustinsParticle::setID()
{
  m_ID = reinterpret_cast<uint64_t>(this);
}

const uint64_t& JustinsParticle::ID() const
//const Real& JustinsParticle::ID() const
{
  return m_ID;
}

uint64_t& JustinsParticle::ID()
//Real& JustinsParticle::ID()
{
  return m_ID;
}

void JustinsParticle::setWeight(const Real a_weight)
{
  m_weight = a_weight;
}

const Real& JustinsParticle::weight() const
{
  return m_weight;
}

Real& JustinsParticle::weight()
{
  return m_weight;
}

//
// set/get the old position
//

void JustinsParticle::setOldPosition(const RealVect& a_position_old)
{
  m_position_old = a_position_old;
}

void JustinsParticle::setOldPosition(const Real a_position_old, const int a_dimension)
{
  m_position_old[a_dimension] = a_position_old;
}

const RealVect& JustinsParticle::position_old() const
{
  return m_position_old;
}

RealVect& JustinsParticle::position_old()
{
  return m_position_old;
}

Real JustinsParticle::position_old(const int a_dir) const
{
  return m_position_old[a_dir];
}

//
// get/set the velocity
//

void JustinsParticle::setVelocity(const std::array<Real,3>& a_velocity)
{
  m_velocity = a_velocity;
}

void JustinsParticle::setVelocity(const Real& a_velocity,
                            const int   a_dir)
{
  m_velocity[a_dir] = a_velocity;
}

std::array<Real,3>& JustinsParticle::velocity()
{
  return m_velocity;
}

const std::array<Real,3>& JustinsParticle::velocity() const
{
  return m_velocity;
}

Real JustinsParticle::velocity(const int a_dir) const
{
  return m_velocity[a_dir];
}

//
// set/get the old velocity
//

void JustinsParticle::setOldVelocity(const std::array<Real,3>& a_velocity_old)
{
  m_velocity_old = a_velocity_old;
}

std::array<Real,3>& JustinsParticle::velocity_old()
{
  return m_velocity_old;
}

const std::array<Real,3>& JustinsParticle::velocity_old() const
{
  return m_velocity_old;
}

//
// electric_field functions
//

void JustinsParticle::setElectricField( const std::array<Real,3>&  a_electric_field )
{
  m_electric_field = a_electric_field;
}

void JustinsParticle::setElectricField( const Real&  a_electric_field,
                                        const int    a_dir )
{
  m_electric_field[a_dir] = a_electric_field;
}

std::array<Real,3>& JustinsParticle::electric_field()
{
  return m_electric_field;
}

const std::array<Real,3>& JustinsParticle::electric_field() const
{
  return m_electric_field;
}

Real JustinsParticle::electric_field(const int a_dir) const
{
  return m_electric_field[a_dir];
}

//
// magnetic_field functions
//

void JustinsParticle::setMagneticField( const std::array<Real,3>&  a_magnetic_field )
{
  m_magnetic_field = a_magnetic_field;
}

void JustinsParticle::setMagneticField( const Real&  a_magnetic_field,
                                        const int    a_dir )
{
  m_magnetic_field[a_dir] = a_magnetic_field;
}

std::array<Real,3>& JustinsParticle::magnetic_field()
{
  return m_magnetic_field;
}

const std::array<Real,3>& JustinsParticle::magnetic_field() const
{
  return m_magnetic_field;
}

Real JustinsParticle::magnetic_field(const int a_dir) const
{
  return m_magnetic_field[a_dir];
}

#if CH_SPACEDIM<3

// virtual position functions (for 1D/2D sims)

void JustinsParticle::setPositionVirt(const std::array<Real,2>& a_pos_virt)
{
  m_pos_virt = a_pos_virt;
}

void JustinsParticle::setPositionVirt(const Real a_pos_virt,
                                      const int  a_comp)
{
  m_pos_virt[a_comp] = a_pos_virt;
}

std::array<Real,2>& JustinsParticle::position_virt()
{
  return m_pos_virt;
}

const std::array<Real,2>& JustinsParticle::position_virt() const
{
  return m_pos_virt;
}

Real JustinsParticle::position_virt(const int a_comp) const
{
  return m_pos_virt[a_comp];
}

#endif

//////////

bool JustinsParticle::operator == (const JustinsParticle& a_p) const
{
  return ( m_ID        == a_p.m_ID       &&
           m_kill_tag  == a_p.m_kill_tag &&
           m_num_suborbits  == a_p.m_num_suborbits &&
           m_weight    == a_p.m_weight   &&
           m_position  == a_p.m_position &&
           m_position_old  == a_p.m_position_old &&
           m_velocity  == a_p.m_velocity &&
           m_velocity_old  == a_p.m_velocity &&
           m_electric_field  == a_p.m_electric_field &&
#if CH_SPACEDIM<3
           m_magnetic_field  == a_p.m_magnetic_field &&
           m_pos_virt  == a_p.m_pos_virt );
#else
           m_magnetic_field  == a_p.m_magnetic_field );

#endif
}

bool JustinsParticle::operator == (const JustinsParticle* a_p) const
{
  return (*this == *a_p);
}

bool JustinsParticle::operator != (const JustinsParticle& a_p) const
{
  return !(*this == a_p);
}

int JustinsParticle::size() const
{
  return ( BinItem::size() + sizeof(m_weight) + sizeof(m_ID)
                           + sizeof(m_position_old) 
#if CH_SPACEDIM<3
                           + sizeof(m_pos_virt)
#endif
                           //+ sizeof(m_electric_field) + sizeof(m_magnetic_field)
                           + sizeof(m_velocity) + sizeof(m_velocity_old) );
}

int JustinsParticle::sizeOutput() const
{ // don't include fields or kill_tag in output
  return ( BinItem::size() + sizeof(m_weight) + sizeof(m_ID) 
#if CH_SPACEDIM==1
                           + sizeof(m_pos_virt)
#elif CH_SPACEDIM==2
                           + sizeof(m_pos_virt[0])
#endif
                           + sizeof(m_velocity) );
}

// Write a linear (binary) representation of the internal data.
// Assumes that sufficient memory for the buffer has already been
// allocated by the caller.
void JustinsParticle::linearOut(void* buf) const
{
   Real* buffer = (Real*)buf;

   *buffer++ = m_weight;

   for(int i=0; i<SpaceDim; ++i) {
      *buffer++ = m_position[i];
   }
   
   for(int i=0; i<SpaceDim; ++i) {
      *buffer++ = m_position_old[i];
   }
   
#if CH_SPACEDIM<3
   for(int i=0; i<2; i++) {
      *buffer++ = m_pos_virt[i];
   }
#endif
   
   for(int i=0; i<3; i++) {
      *buffer++ = m_velocity[i];
   }
   
   for(int i=0; i<3; i++) {
      *buffer++ = m_velocity_old[i];
   }
   /*
   for(int i=0; i<3; i++) {
      *buffer++ = m_electric_field[i];
   }
   
   for(int i=0; i<3; i++) {
      *buffer++ = m_magnetic_field[i];
   }
   */
   *buffer = m_ID;
   
}

// allocated by the caller.
void JustinsParticle::linearOutOutput(void* buf) const
{
   Real* buffer = (Real*)buf;
   
   *buffer++ = m_weight;

   for(int i=0; i<SpaceDim; ++i) {
      *buffer++ = m_position[i];
   }
   
#if CH_SPACEDIM<3
   for(int i=0; i<3-SpaceDim; i++) {
      *buffer++ = m_pos_virt[i];
   }
#endif
   
   for(int i=0; i<3; i++) {
      *buffer++ = m_velocity[i];
   }
   
   *buffer = m_ID;
   
}

// Read a linear (binary) representation of the internal data.
// Assumes the buffer contains the correct data.
void JustinsParticle::linearIn(void* buf)
{
   Real* buffer = (Real*)buf;
   
   m_weight = *buffer++;
   
   for(int i=0; i<SpaceDim; ++i) {
      m_position[i] = *buffer++;
   }
   
   for(int i=0; i<SpaceDim; ++i) {
      m_position_old[i] = *buffer++;
   }
   
#if CH_SPACEDIM<3
   for(int i=0; i<2; i++) {
      m_pos_virt[i] = *buffer++;
   }
#endif
   
   for(int i=0; i<3; i++) {
      m_velocity[i] = *buffer++;
   }
   
   for(int i=0; i<3; i++) {
      m_velocity_old[i] = *buffer++;
   }
   /*
   for(int i=0; i<3; i++) {
      m_electric_field[i] = *buffer++;
   }
   
   for(int i=0; i<3; i++) {
      m_magnetic_field[i] = *buffer++;
   }
   */
   m_ID = *buffer;

}

std::ostream & operator<<(std::ostream& ostr, const JustinsParticle& p)
{
  ostr << " JustinsParticle : " << std::endl;
  ostr << " weight " << p.weight() << std::endl;
  ostr << " ID " << p.ID() << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << std::endl << " velocity ( ";
  for ( int i=0; i<3; ++i ){ ostr << " " << p.velocity(i); }
  ostr << " ) ";
  ostr << std::endl << " electric_field ( ";
  for ( int i=0; i<3; ++i ){ ostr << " " << p.electric_field(i); }
  ostr << std::endl << " magnetic_field ( ";
  for ( int i=0; i<3; ++i ){ ostr << " " << p.magnetic_field(i); }
  ostr << " ) " << std::endl;
  return ostr;
}

#include "NamespaceFooter.H"
