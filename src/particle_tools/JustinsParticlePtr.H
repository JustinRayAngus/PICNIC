#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _JUSTINSPARTICLEPTR_H_
#define _JUSTINSPARTICLEPTR_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"
#include "JustinsParticle.H"
#include "NamespaceHeader.H"

//  A pointer to a JustinsParticle. 

class JustinsParticlePtr : public BinItem
{
public:

  /// default constructor
  JustinsParticlePtr();

  /// same as define()
  JustinsParticlePtr(JustinsParticle& a_part);

  /// default
  virtual ~JustinsParticlePtr();

  /// for non-default constructor
  void define(JustinsParticle& a_part);

  ///////////////
  //

  virtual bool operator==(const JustinsParticlePtr& a_p) const;
  virtual bool operator==(const JustinsParticlePtr* a_p) const;
  virtual bool operator!=(const JustinsParticlePtr& a_p) const;

  JustinsParticle* getPointer();

  /// {\bf linearization functions}

  /** Returns the size, in number of bytes, of a flat linear
      representation of the data in this object.
   */
  virtual int size() const;

  // return size, in number of bytes, of a comp of the internal data
  //virtual size_t sizeOfComp(const int a_comp) const;

  ///
  /** Write a linear binary representation of the internal data. Assumes that
      sufficient memory for the buffer has already been allocated by
      the caller.
  */
  virtual void linearOut(void* buf) const;

  ///
  /** Read a linear binary representation of the internal data.
      Assumes that the buffer has the correct data.
  */
  virtual void linearIn(void* buf);

  friend std::ostream & operator<<(std::ostream& ostr,const JustinsParticlePtr& p);

protected:

  JustinsParticle* m_part_ptr;

};

std::ostream & operator<<(std::ostream& ostr, const JustinsParticlePtr& p);

#include "NamespaceFooter.H"

#endif // include guard

