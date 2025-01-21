#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PHOTON_PARTICLE_PTR_H_
#define _PHOTON_PARTICLE_PTR_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"
#include "PhotonParticle.H"
#include "NamespaceHeader.H"

//  A pointer to a PhotonParticle. 

class PhotonParticlePtr : public BinItem
{
public:

  /// default constructor
  PhotonParticlePtr();

  /// same as define()
  PhotonParticlePtr(PhotonParticle& a_part);

  /// default
  virtual ~PhotonParticlePtr();

  /// for non-default constructor
  void define(PhotonParticle& a_part);

  ///////////////
  //

  virtual bool operator==(const PhotonParticlePtr& a_p) const;
  virtual bool operator==(const PhotonParticlePtr* a_p) const;
  virtual bool operator!=(const PhotonParticlePtr& a_p) const;

  PhotonParticle* getPointer();

  virtual int size() const;

  virtual void linearOut(void* buf) const;

  virtual void linearIn(void* buf);

  friend std::ostream & operator<<(std::ostream& ostr,const PhotonParticlePtr& p);

protected:

  PhotonParticle* m_part_ptr;

};

std::ostream & operator<<(std::ostream& ostr, const PhotonParticlePtr& p);

#include "NamespaceFooter.H"

#endif // include guard

