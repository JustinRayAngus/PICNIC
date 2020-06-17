/***
 *
 * DomainGrid (mesh) class header file
 *
***/

#ifndef DomainGrid_h
#define DomainGrid_h

#include "ParmParse.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"

#include "NamespaceHeader.H"

class DomainGrid
{

public:

  DomainGrid( ParmParse&          ppgrid,
        const int                 numGhosts,
        const ProblemDomain&      domain,
        const DisjointBoxLayout&  grids );
  
  virtual ~DomainGrid() {;}

  //static DomainGrid* mesh; // pointer to a DomainGrid instance

 
  // functions to return private members
  //
  virtual const ProblemDomain& getDomain() const {return m_domain;}
  
  virtual const DisjointBoxLayout& getDBL() const {return m_grids;}
 
  virtual const RealVect& getdX() const {return m_dX;}
  
  virtual const RealVect& getXmin() const {return m_Xmin;}
  
  virtual const RealVect& getXmax() const {return m_Xmax;}
 
  virtual const bool axisymmetric() const { return m_axisymmetric; }
 
  virtual const int ghosts() const { return m_ghosts; }
  
  virtual const double getMappedCellVolume() const { return m_mapped_cell_volume; }
  
  virtual const RealVect& getMappedFaceArea() const { return m_mapped_face_area; }
  
  virtual const LevelData<FArrayBox>& getXcc() const {return m_Xcc;}
  
  virtual const LevelData<FluxBox>& getXfc() const {return m_Xfc;}

private:
  
  bool m_axisymmetric;
  int m_ghosts;
  
  RealVect m_Xmin, m_Xmax, m_dX;
  double m_mapped_cell_volume;
  RealVect m_mapped_face_area;  
  
  ProblemDomain m_domain;
  DisjointBoxLayout m_grids;
  
  LevelData<FArrayBox> m_Xcc; // coords at cell center
  LevelData<FluxBox> m_Xfc;
  //LevelData<EdgeDataBox> m_edge_centered_real_coords;

  void setRealCoords();

};

#include "NamespaceFooter.H"

#endif
