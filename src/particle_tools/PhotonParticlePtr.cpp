#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PhotonParticlePtr.H"
#include "NamespaceHeader.H"

/// default constructor
PhotonParticlePtr::PhotonParticlePtr() :
  BinItem()
{}

PhotonParticlePtr::~PhotonParticlePtr()
{
   m_part_ptr = NULL;
   delete m_part_ptr;
}

PhotonParticlePtr::PhotonParticlePtr(PhotonParticle& a_photon_part)
  :
  BinItem(a_photon_part.position()),
  m_part_ptr( &a_photon_part )
{
  //if(!procID()) std::cout << "JRA: contructing PhotonParticlePtr " << std::endl;
  //setPosition(a_photon_part.position());
  //m_part_ptr = &a_photon_part;
}

void PhotonParticlePtr::define(PhotonParticle& a_photon_part)
{
  setPosition(a_photon_part.position());
  m_part_ptr = &a_photon_part;
}

// return the pointer
PhotonParticle* PhotonParticlePtr::getPointer()
{
  return m_part_ptr;
}


bool PhotonParticlePtr::operator == (const PhotonParticlePtr& a_p) const
{
  return ( m_position  == a_p.m_position &&
           m_part_ptr  == a_p.m_part_ptr );
}

bool PhotonParticlePtr::operator == (const PhotonParticlePtr* a_p) const
{
  return (*this == *a_p);
}

bool PhotonParticlePtr::operator != (const PhotonParticlePtr& a_p) const
{
  return !(*this == a_p);
}

int PhotonParticlePtr::size() const
{
  return ( BinItem::size() + sizeof(m_part_ptr) );
}

// Write a linear (binary) representation of the internal data.
// Assumes that sufficient memory for the buffer has already been
// allocated by the caller.
void PhotonParticlePtr::linearOut(void* buf) const
{
   Real* buffer = (Real*)buf;
   D_TERM6( *buffer++ = m_position[0];,
            *buffer++ = m_position[1];,
            *buffer++ = m_position[2];,
            *buffer++ = m_position[3];,
            *buffer++ = m_position[4];,
            *buffer++ = m_position[5];);
   
}


// Read a linear (binary) representation of the internal data.
// Assumes the buffer contains the correct data.
void PhotonParticlePtr::linearIn(void* buf)
{
   Real* buffer = (Real*)buf;
   D_TERM6( m_position[0] = *buffer++;,
            m_position[1] = *buffer++;,
            m_position[2] = *buffer++;,
            m_position[3] = *buffer++;,
            m_position[4] = *buffer++;,
            m_position[5] = *buffer++;);
   
}

std::ostream & operator<<(std::ostream& ostr, const PhotonParticlePtr& p)
{
  ostr << " PhotonParticlePtr : " << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) " << std::endl;
  return ostr;
}

#include "NamespaceFooter.H"
