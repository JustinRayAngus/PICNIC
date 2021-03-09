#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _JUSTINSPARTICLE_H_
#define _JUSTINSPARTICLE_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"
#include "NamespaceHeader.H"

//  A basic particle class. 
//  Includes weight, position, velocity, and acceleration
//  position is separated into those corresonding to the 
//  SpaceDim of the build and virtual components
//

class JustinsParticle : public BinItem
{
public:

  /// default constructor
  JustinsParticle();

  /// same as define()
  JustinsParticle( const Real       a_weight,
                   const RealVect&  a_position,
                   const RealVect&  a_velocity = RealVect::Zero );

  /// default
  virtual ~JustinsParticle();

  /// for non-default constructor
  void define( const Real       a_weight,
               const RealVect&  a_position,
               const RealVect&  a_velocity = RealVect::Zero );
  
  /// set particle ID
  void setID();
  
  void setID(const uint64_t a_ID);
  //void setID(const Real a_ID);

  // retrieve reference ID
  uint64_t& ID();
  //Real& ID();

  // retrieve const. reference ID
  const uint64_t& ID() const;
  //const Real& ID() const;

  /// set particle weight
  void setWeight(const Real a_weight);

  // retrieve reference weight field
  Real& weight();

  // retrieve const. reference weight field
  const Real& weight() const;
  
  /// set velocity field
  void setVelocity(const RealVect& a_velocity);
  void setVelocity(const Real& a_velocity, const int a_dir);

  // retrieve reference velocity field
  RealVect& velocity();

  // retrieve const. reference velocity field
  const RealVect& velocity() const;

  // retrieve vel. comp. in a given direction
  Real velocity(const int a_dir) const;

  /// set acceleration field
  void setAcceleration(const std::array<Real,3>& a_acceleration);

  // retrieve reference acceleration field
  std::array<Real,3>& acceleration();

  // retrieve const. reference acceleration field
  const std::array<Real,3>& acceleration() const;

  // retrieve acc. comp. in a given direction
  Real acceleration(const int a_dir) const;

  //    
  //  define functions for virtual position (for 1D/2D sim)
  //
  
  /// set virtual position field
  void setPositionVirt(const std::array<Real,3-SpaceDim>& a_pos_virt);
  void setPositionVirt(const Real& a_pos_virt, const int a_dir);

  // retrieve reference virtual position field
  std::array<Real,3-SpaceDim>& positionVirt();

  // retrieve const. reference virtual position field
  const std::array<Real,3-SpaceDim>& positionVirt() const;

  // retrieve virt. pos. comp. in a given direction
  Real positionVirt(const int a_dir) const;

  //
  //  define functions for virtual velocity (for 1D/2D sims)
  //

  /// set virtual velocity field
  void setVelocityVirt(const std::array<Real,3-SpaceDim>& a_vel_virt);
  void setVelocityVirt(const Real& a_vel_virt, const int a_dir);

  // retrieve reference virtual velocity field
  std::array<Real,3-SpaceDim>& velocityVirt();

  // retrieve const. reference virtual velocity field
  const std::array<Real,3-SpaceDim>& velocityVirt() const;

  // retrieve virt. vel. comp. in a given direction
  Real velocityVirt(const int a_dir) const;
  
  //
  ///////////////
  //

  virtual bool operator==(const JustinsParticle& a_p) const;
  virtual bool operator==(const JustinsParticle* a_p) const;
  virtual bool operator!=(const JustinsParticle& a_p) const;

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

  friend std::ostream & operator<<(std::ostream& ostr,const JustinsParticle& p);

protected:

  uint64_t m_ID;
  //Real     m_ID;
  Real     m_weight;
  RealVect m_velocity;
  std::array<Real,3> m_acceleration; // virtual acceleration for 1D/2D sims
  std::array<Real,3-SpaceDim> m_pos_virt; // virtual position for 1D/2D sims
  std::array<Real,3-SpaceDim> m_vel_virt; // virtual velocity for 1D/2D sims
};

std::ostream & operator<<(std::ostream& ostr, const JustinsParticle& p);

#include "NamespaceFooter.H"

#endif // include guard

