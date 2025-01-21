#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PhotonParticle.H"
#include "NamespaceHeader.H"

/// default constructor
PhotonParticle::PhotonParticle() :
  BinItem()
{}

PhotonParticle::~PhotonParticle()
{}

PhotonParticle::PhotonParticle( const Real                 a_weight,
                                const RealVect&            a_position,
                                const std::array<Real,3>&  a_velocity )
  :
  BinItem(a_position),
  m_weight(a_weight),
  m_velocity(a_velocity)
{
  setID();
}

void PhotonParticle::define( const Real                 a_weight,
                             const RealVect&            a_position,
                             const std::array<Real,3>&  a_velocity )
{
  setID();
  setWeight(a_weight);
  setPosition(a_position);
  setVelocity(a_velocity);
}

void PhotonParticle::setKillTag()
{
  m_kill_tag = 1;
}
  
const int& PhotonParticle::killTag() const
{
  return m_kill_tag;
}

void PhotonParticle::setID(const uint64_t a_ID)
{
  m_ID = a_ID;
}

void PhotonParticle::setID()
{
  m_ID = reinterpret_cast<uint64_t>(this);
}

const uint64_t& PhotonParticle::ID() const
{
  return m_ID;
}

uint64_t& PhotonParticle::ID()
{
  return m_ID;
}

void PhotonParticle::setWeight(const Real a_weight)
{
  m_weight = a_weight;
}

const Real& PhotonParticle::weight() const
{
  return m_weight;
}

Real& PhotonParticle::weight()
{
  return m_weight;
}

//
// get/set the velocity
//

void PhotonParticle::setVelocity(const std::array<Real,3>& a_velocity)
{
  m_velocity = a_velocity;
}

void PhotonParticle::setVelocity(const Real& a_velocity,
                            const int   a_dir)
{
  m_velocity[a_dir] = a_velocity;
}

std::array<Real,3>& PhotonParticle::velocity()
{
  return m_velocity;
}

const std::array<Real,3>& PhotonParticle::velocity() const
{
  return m_velocity;
}

Real PhotonParticle::velocity(const int a_dir) const
{
  return m_velocity[a_dir];
}


//////////

bool PhotonParticle::operator == (const PhotonParticle& a_p) const
{
  return ( m_ID        == a_p.m_ID       &&
           m_kill_tag  == a_p.m_kill_tag &&
           m_weight    == a_p.m_weight   &&
           m_position  == a_p.m_position &&
           m_velocity  == a_p.m_velocity );
}

bool PhotonParticle::operator == (const PhotonParticle* a_p) const
{
  return (*this == *a_p);
}

bool PhotonParticle::operator != (const PhotonParticle& a_p) const
{
  return !(*this == a_p);
}

int PhotonParticle::size() const
{
  return ( BinItem::size() + sizeof(m_weight) + sizeof(m_ID)
                           + sizeof(m_velocity) );
}

int PhotonParticle::sizeOutput() const
{
  return ( BinItem::size() + sizeof(m_weight) + sizeof(m_ID) 
                           + sizeof(m_velocity) );
}

// Write a linear (binary) representation of the internal data.
// Assumes that sufficient memory for the buffer has already been
// allocated by the caller.
void PhotonParticle::linearOut(void* buf) const
{
   Real* buffer = (Real*)buf;

   *buffer++ = m_weight;
   for(int i=0; i<SpaceDim; ++i) { *buffer++ = m_position[i]; }
   for(int i=0; i<3; i++) { *buffer++ = m_velocity[i]; }
   *buffer = m_ID;
   
}

// allocated by the caller.
void PhotonParticle::linearOutOutput(void* buf) const
{
   Real* buffer = (Real*)buf;
   
   *buffer++ = m_weight;
   for(int i=0; i<SpaceDim; ++i) { *buffer++ = m_position[i]; }
   for(int i=0; i<3; i++) { *buffer++ = m_velocity[i]; }
   *buffer = m_ID;
   
}

// Read a linear (binary) representation of the internal data.
// Assumes the buffer contains the correct data.
void PhotonParticle::linearIn(void* buf)
{
   Real* buffer = (Real*)buf;
   
   m_weight = *buffer++;
   for(int i=0; i<SpaceDim; ++i) { m_position[i] = *buffer++; }
   for(int i=0; i<3; i++) { m_velocity[i] = *buffer++; }
   m_ID = *buffer;

}

std::ostream & operator<<(std::ostream& ostr, const PhotonParticle& p)
{
  ostr << " PhotonParticle : " << std::endl;
  ostr << " weight " << p.weight() << std::endl;
  ostr << " ID " << p.ID() << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) ";
  ostr << std::endl << " velocity ( ";
  for ( int i=0; i<3; ++i ){ ostr << " " << p.velocity(i); }
  ostr << " ) ";
  return ostr;
}

#include "NamespaceFooter.H"
