#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PHOTONPARTICLE_H_
#define _PHOTONPARTICLE_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>
#include <array>

#include "BinItem.H"
#include "NamespaceHeader.H"

//  Photon particle class. 
//  Includes weight, position, velocity.
//  The position is separated into those corresonding to the 
//  SpaceDim of the build and virtual components

class PhotonParticle : public BinItem
{
public:

  /// default constructor
  PhotonParticle();

  /// same as define()
  PhotonParticle( const Real       a_weight,
                  const RealVect&  a_position,
                  const std::array<Real,3>&  a_velocity = {0,0,0} );

  /// default
  virtual ~PhotonParticle();

  /// for non-default constructor
  void define( const Real       a_weight,
               const RealVect&  a_position,
               const std::array<Real,3>&  a_velocity = {0,0,0} );
               //const RealVect&  a_velocity = RealVect::Zero );
  
  /// set particle ID
  void setID();
  
  void setID(const uint64_t a_ID);

  // retrieve reference ID
  uint64_t& ID();

  // retrieve const. reference ID
  const uint64_t& ID() const;
  
  void setKillTag();
  
  // retrieve reference to kill tag
  const int& killTag() const;
  
  // set particle weight
  void setWeight(const Real a_weight);

  // retrieve reference weight field
  Real& weight();

  // retrieve const. reference weight field
  const Real& weight() const;
  
  // set the old position
  void setOldPosition(const RealVect& a_position_old);

  // set the velocity
  void setVelocity(const std::array<Real,3>& a_velocity);

  // set a_dir component of velocity
  void setVelocity(const Real& a_velocity, const int a_dir);

  // retrieve reference virtual velocity field
  std::array<Real,3>& velocity();

  // retrieve const. reference virtual velocity field
  const std::array<Real,3>& velocity() const;
  
  // retrieve a_dir component of velocity
  Real velocity(const int a_dir) const;

  ///////////////

  virtual bool operator==(const PhotonParticle& a_p) const;
  virtual bool operator==(const PhotonParticle* a_p) const;
  virtual bool operator!=(const PhotonParticle& a_p) const;

  /// {\bf linearization functions}

  /** Returns the size, in number of bytes, of a flat linear
      representation of the data in this object.
   */
  virtual int size() const;

  // return size, in number of bytes, of a comp of the internal data
  //virtual size_t sizeOfComp(const int a_comp) const;
  
  int sizeOutput() const;
  void linearOutOutput(void* buf) const;

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

  friend std::ostream & operator<<(std::ostream& ostr,const PhotonParticle& p);

protected:

  uint64_t m_ID;
  int      m_kill_tag = 0;
  Real     m_weight;
  std::array<Real,3> m_velocity;

};

std::ostream & operator<<(std::ostream& ostr, const PhotonParticle& p);

#include "NamespaceFooter.H"

#endif // include guard

