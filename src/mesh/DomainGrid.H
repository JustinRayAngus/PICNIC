/***
 *
 * DomainGrid (mesh) class header file
 *
***/

#ifndef DomainGrid_h
#define DomainGrid_h

#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"

#include "BoundaryBoxLayout.H"

#include "NamespaceHeader.H"

class DomainGrid
{

public:

  DomainGrid( const ProblemDomain&      domain,
              const DisjointBoxLayout&  grids,
              const int                 numGhosts,
              const Real                length_scale );
  
  virtual ~DomainGrid() {;}

  //static DomainGrid* mesh; // pointer to a DomainGrid instance

 
  // functions to return private members
  //
  virtual const ProblemDomain& getDomain() const {return m_domain;}
  
  virtual const DisjointBoxLayout& getDBL() const {return m_grids;}
  
  virtual const BoundaryBoxLayoutPtrVect& getBoundaryLayout() const {return m_domain_bdry_layout;}
  
  virtual const BoundaryBoxLayoutPtrVect& getPeriodicBoundaryLayout() const {return m_periodic_bdry_layout;}
 
  virtual const RealVect& getdX() const {return m_dX;}
  
  virtual const RealVect& getXmin() const {return m_Xmin;}
  
  virtual const RealVect& getXmax() const {return m_Xmax;}
 
  virtual const bool axisymmetric() const { return m_axisymmetric; }
  
  virtual const bool anticyclic() const { return m_anticyclic; }
  
  virtual const std::string geomType() const { return m_geom_type; }
  
  virtual const bool writeJacobians() const { return m_write_jacobians; }
  
  virtual const bool writeCorrectedJacobians() const { return m_write_corrected_jacobians; }
  
  virtual const bool writeMaskedJacobians() const { return m_write_masked_jacobians; }
 
  virtual const int ghosts() const { return m_ghosts; }
  
  const Copier forwardCopier() const { return m_forwardCopier; }
  
  const Copier reverseCopier() const { return m_reverseCopier; }
  
  const Real getMappedCellVolume() const { return m_mapped_cell_volume; }
  
  const Real getVolumeScale() const { return m_volume_scale; }
  
  const Real getAreaScale() const { return m_area_scale; }
  
  virtual const RealVect& getMappedFaceArea() const { return m_mapped_face_area; }
  
  virtual const LevelData<FArrayBox>& getXcc() const {return m_Xcc;}
  
  virtual const LevelData<FluxBox>& getXfc() const {return m_Xfc;}
  
  virtual const LevelData<EdgeDataBox>& getXec() const {return m_Xec;}
  
  virtual const LevelData<NodeFArrayBox>& getXnc() const {return m_Xnc;}
  
  virtual const LevelData<FArrayBox>& getJcc() const {return m_Jcc;}
  
  virtual const LevelData<FluxBox>& getJfc() const {return m_Jfc;}
  
  virtual const LevelData<EdgeDataBox>& getJec() const {return m_Jec;}
  
  virtual const LevelData<NodeFArrayBox>& getJnc() const {return m_Jnc;}
  
  virtual const LevelData<FluxBox>& getMaskedJfc() const {return m_masked_Jfc;}
  
  virtual const LevelData<EdgeDataBox>& getMaskedJec() const {return m_masked_Jec;}
  
  virtual const LevelData<NodeFArrayBox>& getMaskedJnc() const {return m_masked_Jnc;}
  
  virtual const LevelData<NodeFArrayBox>& getNodeMask() const {return m_node_mask;}
  
  virtual const LevelData<FluxBox>& getCorrectedJfc() const {return m_corrected_Jfc;}
  
  virtual const LevelData<EdgeDataBox>& getCorrectedJec() const {return m_corrected_Jec;}
  
  virtual const LevelData<NodeFArrayBox>& getCorrectedJnc() const {return m_corrected_Jnc;}


private:
  
  bool m_axisymmetric;
  bool m_anticyclic;
  bool m_write_jacobians;
  bool m_write_corrected_jacobians;
  bool m_write_masked_jacobians;
  int m_ghosts;
 
  Copier m_forwardCopier; 
  Copier m_reverseCopier; 
      
  std::string m_geom_type; 
 
  enum VOLUME_CORRECTION {
     NONE,
     VERBONCOEUR,
     CONSERVATIVE,
  };
      
  VOLUME_CORRECTION m_volume_correction;

  RealVect m_Xmin, m_Xmax, m_dX;
  Real m_mapped_cell_volume;
  RealVect m_mapped_face_area;  
  Real m_volume_scale;
  Real m_area_scale;
  
  ProblemDomain m_domain;
  DisjointBoxLayout m_grids;
  BoundaryBoxLayoutPtrVect m_domain_bdry_layout;
  BoundaryBoxLayoutPtrVect m_periodic_bdry_layout;
  
  LevelData<FArrayBox> m_Xcc;     // coords at cell center
  LevelData<FluxBox> m_Xfc;       // coords at cell faces
  LevelData<EdgeDataBox> m_Xec;   // coords at cell edges
  LevelData<NodeFArrayBox> m_Xnc; // coords at cell nodes
  
  LevelData<FArrayBox> m_Jcc;     // Jacobian at cell center
  LevelData<FluxBox> m_Jfc;       // Jacobian at cell faces
  LevelData<EdgeDataBox> m_Jec;   // Jacobian at cell edges
  LevelData<NodeFArrayBox> m_Jnc; // Jacobian at cell nodes
  
  // masked Jacobians for nodal containers with 1/2 multiple on shared
  // nodes. These are for taking global volume moments.
  LevelData<FluxBox> m_masked_Jfc;
  LevelData<EdgeDataBox> m_masked_Jec;
  LevelData<NodeFArrayBox> m_masked_Jnc;
  LevelData<NodeFArrayBox> m_node_mask;
  
  // corrected Jacobian needed for correct charge and current 
  // density accumulation in curvilinear coords (see Verboncoeur)
  LevelData<FluxBox> m_corrected_Jfc;
  LevelData<EdgeDataBox> m_corrected_Jec;
  LevelData<NodeFArrayBox> m_corrected_Jnc;

  void setRealCoords();
  
  void setJacobian();

  void displayJacobianData();
  
  void defineBoundaryBoxLayout();

};

#include "NamespaceFooter.H"

#endif
