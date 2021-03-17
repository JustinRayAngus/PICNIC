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

JustinsParticle::JustinsParticle( const Real       a_weight,
                                  const RealVect&  a_position,
                                  const std::array<Real,3>&  a_velocity )
  :
  BinItem(a_position),
  m_weight(a_weight),
  m_velocity(a_velocity)
{
  std::array<Real,3> thisElectricField = {0,0,0};
  setElectricField(thisElectricField);
  if(SpaceDim<3) {
     std::array<Real,3-SpaceDim> thisPosVirt;
     for(int i=0; i<3-SpaceDim; i++) {
        thisPosVirt[i] = 0.0 + i*0.0;
     }
     setPositionVirt(thisPosVirt);
  }
  setID();
}

void JustinsParticle::define( const Real       a_weight,
                              const RealVect&  a_position,
                              const std::array<Real,3>&  a_velocity )
{
  setID();
  setWeight(a_weight);
  setPosition(a_position);
  setVelocity(a_velocity);
  std::array<Real,3> thisElectricField = {0,0,0};
  setElectricField(thisElectricField);
  if(SpaceDim<3) {
     std::array<Real,3-SpaceDim> thisPosVirt;
     for(int i=0; i<3-SpaceDim; i++) {
        thisPosVirt[i] = 0.0 + i*0.0;
     }
     setPositionVirt(thisPosVirt);
  }
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

// velocity functions
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

// electric_field functions
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
// virtual position functions (for 1D/2D sims)
//
void JustinsParticle::setPositionVirt(const std::array<Real,3-SpaceDim>& a_pos_virt)
{
  m_pos_virt = a_pos_virt;
}

void JustinsParticle::setPositionVirt(const Real& a_pos_virt,
                            const int   a_dir)
{
  m_pos_virt[a_dir] = a_pos_virt;
}

std::array<Real,3-SpaceDim>& JustinsParticle::positionVirt()
{
  return m_pos_virt;
}

const std::array<Real,3-SpaceDim>& JustinsParticle::positionVirt() const
{
  return m_pos_virt;
}

Real JustinsParticle::positionVirt(const int a_dir) const
{
  return m_pos_virt[a_dir];
}

//
//////////
//

bool JustinsParticle::operator == (const JustinsParticle& a_p) const
{
  return ( m_ID        == a_p.m_ID       &&
           m_weight    == a_p.m_weight   &&
           m_position  == a_p.m_position &&
           m_velocity  == a_p.m_velocity &&
           m_electric_field  == a_p.m_electric_field &&
           m_pos_virt  == a_p.m_pos_virt );
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
                           + sizeof(m_pos_virt)
                           + sizeof(m_velocity) + sizeof(m_electric_field) );
}

int JustinsParticle::sizeOutput() const
{ // don't include fields in output
  return ( BinItem::size() + sizeof(m_weight) + sizeof(m_ID) 
                           + sizeof(m_pos_virt)
                           + sizeof(m_velocity) );
}

// Write a linear (binary) representation of the internal data.
// Assumes that sufficient memory for the buffer has already been
// allocated by the caller.
void JustinsParticle::linearOut(void* buf) const
{
   Real* buffer = (Real*)buf;

   *buffer++ = m_weight;

   D_TERM6( *buffer++ = m_position[0];,
            *buffer++ = m_position[1];,
            *buffer++ = m_position[2];,
            *buffer++ = m_position[3];,
            *buffer++ = m_position[4];,
            *buffer++ = m_position[5];);
   
   for(int i=0; i<3-SpaceDim; i++) {
      *buffer++ = m_pos_virt[i];
   }
   
   for(int i=0; i<3; i++) {
      *buffer++ = m_velocity[i];
   }

   for(int i=0; i<3; i++) {
      *buffer++ = m_electric_field[i];
   }

   *buffer = m_ID;

}

// allocated by the caller.
void JustinsParticle::linearOutOutput(void* buf) const
{
   Real* buffer = (Real*)buf;
   
   *buffer++ = m_weight;

   D_TERM6( *buffer++ = m_position[0];,
            *buffer++ = m_position[1];,
            *buffer++ = m_position[2];,
            *buffer++ = m_position[3];,
            *buffer++ = m_position[4];,
            *buffer++ = m_position[5];);
   
   for(int i=0; i<3-SpaceDim; i++) {
      *buffer++ = m_pos_virt[i];
   }
   
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
   
   D_TERM6( m_position[0] = *buffer++;,
            m_position[1] = *buffer++;,
            m_position[2] = *buffer++;,
            m_position[3] = *buffer++;,
            m_position[4] = *buffer++;,
            m_position[5] = *buffer++;);
   
   for(int i=0; i<3-SpaceDim; i++) {
      m_pos_virt[i] = *buffer++;
   }
   
   for(int i=0; i<3; i++) {
      m_velocity[i] = *buffer++;
   }

   for(int i=0; i<3; i++) {
      m_electric_field[i] = *buffer++;
   }
   
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
  ostr << " ) " << std::endl;
  return ostr;
}

#include "NamespaceFooter.H"
