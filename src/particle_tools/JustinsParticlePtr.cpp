#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "JustinsParticlePtr.H"
#include "NamespaceHeader.H"

/// default constructor
JustinsParticlePtr::JustinsParticlePtr() :
  BinItem()
{}

JustinsParticlePtr::~JustinsParticlePtr()
{
   m_part_ptr = NULL;
   delete m_part_ptr;
}

JustinsParticlePtr::JustinsParticlePtr(JustinsParticle& a_justins_part)
  :
  BinItem(a_justins_part.position()),
  m_part_ptr( &a_justins_part )
{
  //if(!procID()) std::cout << "JRA: contructing JustinsParticlePtr " << std::endl;
  //setPosition(a_justins_part.position());
  //m_part_ptr = &a_justins_part;
}

void JustinsParticlePtr::define(JustinsParticle& a_justins_part)
{
  setPosition(a_justins_part.position());
  m_part_ptr = &a_justins_part;
}

// return the pointer
//
JustinsParticle* JustinsParticlePtr::getPointer()
{
  return m_part_ptr;
}


bool JustinsParticlePtr::operator == (const JustinsParticlePtr& a_p) const
{
  return ( m_position  == a_p.m_position &&
           m_part_ptr  == a_p.m_part_ptr );
}

bool JustinsParticlePtr::operator == (const JustinsParticlePtr* a_p) const
{
  return (*this == *a_p);
}

bool JustinsParticlePtr::operator != (const JustinsParticlePtr& a_p) const
{
  return !(*this == a_p);
}

int JustinsParticlePtr::size() const
{
  return ( BinItem::size() + sizeof(m_part_ptr) );
}

// Write a linear (binary) representation of the internal data.
// Assumes that sufficient memory for the buffer has already been
// allocated by the caller.
void JustinsParticlePtr::linearOut(void* buf) const
{
   Real* buffer = (Real*)buf;
   D_TERM6( *buffer++ = m_position[0];,
            *buffer++ = m_position[1];,
            *buffer++ = m_position[2];,
            *buffer++ = m_position[3];,
            *buffer++ = m_position[4];,
            *buffer++ = m_position[5];);
   
  // Do I need to add ptr to buffer? Can only add real values..
  // *buffer = sizeof(m_part_ptr); // JRA, is this correct?
}


// Read a linear (binary) representation of the internal data.
// Assumes the buffer contains the correct data.
void JustinsParticlePtr::linearIn(void* buf)
{
   Real* buffer = (Real*)buf;
   D_TERM6( m_position[0] = *buffer++;,
            m_position[1] = *buffer++;,
            m_position[2] = *buffer++;,
            m_position[3] = *buffer++;,
            m_position[4] = *buffer++;,
            m_position[5] = *buffer++;);
   
  //m_part_ptr = buffer; // Cant reader ptr from Real valued buffer. What to do?
  // m_part_ptr = *buffer; // dont know how to "read" ptr from buffer
}

std::ostream & operator<<(std::ostream& ostr, const JustinsParticlePtr& p)
{
  ostr << " JustinsParticlePtr : " << std::endl;
  ostr << " position ( ";
  for ( int i=0; i<SpaceDim; ++i ){ ostr << " " << p.position(i); }
  ostr << " ) " << std::endl;
  return ostr;
}

#include "NamespaceFooter.H"
