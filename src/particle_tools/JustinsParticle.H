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
//  Includes weight, position, velocity, and fields.
//  The position is separated into those corresonding to the 
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
                   const std::array<Real,3>&  a_velocity = {0,0,0} );
                   //const RealVect&  a_velocity = RealVect::Zero );

  /// default
  virtual ~JustinsParticle();

  /// for non-default constructor
  void define( const Real       a_weight,
               const RealVect&  a_position,
               const std::array<Real,3>&  a_velocity = {0,0,0} );
               //const RealVect&  a_velocity = RealVect::Zero );
  
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
  

  // set the old position
  void setOldPosition(const RealVect& a_position_old);

  // Change one component of the old position
  void setOldPosition(const Real a_position_old, const int a_dimension);

  // Return a reference to the old position
  RealVect& position_old();

  // Return a constant reference to the old position
  const RealVect& position_old() const;

  // Gets one component of the old position
  Real position_old(const int a_dir) const;

  
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
  

  // set the old velocity
  void setOldVelocity(const std::array<Real,3>& a_velocity_old);
 
  // Return a reference to the old velocity
  std::array<Real,3>& velocity_old();

  // Return a constant reference to the old velocity
  const std::array<Real,3>& velocity_old() const;
  

  /// set electric field
  void setElectricField(const std::array<Real,3>& a_electric_field);
  void setElectricField(const Real& a_electric_field, const int a_dir);

  // retrieve reference electric field
  std::array<Real,3>& electric_field();

  // retrieve const. reference electric field
  const std::array<Real,3>& electric_field() const;

  // retrieve electric field comp. in a given direction
  Real electric_field(const int a_dir) const;

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
  
  int sizeOutput() const; // JRA
  void linearOutOutput(void* buf) const; // JRA

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
  std::array<Real,3> m_velocity;
  std::array<Real,3> m_electric_field;
  //std::array<Real,3> m_magnetic_field;
  std::array<Real,3-SpaceDim> m_pos_virt; // virtual position for 1D/2D sims
 
  RealVect  m_position_old;
  std::array<Real,3> m_velocity_old;

};

std::ostream & operator<<(std::ostream& ostr, const JustinsParticle& p);

#include "NamespaceFooter.H"

#endif // include guard

