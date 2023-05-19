#ifndef _GLOBALDOF_H_
#define _GLOBALDOF_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "SpaceUtils.H"

/* Fortran ordering */
#define _BoxIndexDOF_(p,i,N,ndims) \
  { \
    p = i[ndims-1]; \
    for (int c=ndims-2; c>=0; c--) p = (p*N[c]) + i[c];\
  }

#include "NamespaceHeader.H"

template <typename T>
class GlobalDOF
{
  public:

    GlobalDOF<T>()
    {
      m_is_defined = false;
      m_mpi_offset = -1;
    }

    GlobalDOF<T>( const int           a_offset,
                  const int           a_mpi_offset,
                  const LevelData<T>& a_data )
    {
      define( a_offset, a_mpi_offset, a_data );
    }

    virtual ~GlobalDOF<T>() {}

    inline bool isDefined() const
    {
      return m_is_defined;
    }

    inline int mpiOffset() const 
    { 
      return m_mpi_offset; 
    }
    
    inline const LevelData<T>& dofDataLHS() const 
    { 
      CH_assert(isDefined());
      return m_gdofs_lhs; 
    }

    inline const LevelData<T>& dofDataRHS() const 
    { 
      CH_assert(isDefined());
      return m_gdofs_rhs; 
    }

    virtual void define(  const int           a_offset,
                          const int           a_mpi_offset,
                          const LevelData<T>& a_data )
    {
      CH_assert(!isDefined());

      auto grids( a_data.disjointBoxLayout() );
      auto phys_domain( grids.physDomain() );
      const int ncomp(a_data.nComp());
      auto gvec(a_data.ghostVect());
      m_gdofs_lhs.define( grids, ncomp, gvec );
      m_gdofs_rhs.define( grids, ncomp, gvec );

      {
        long box_offset(0);
        for (auto dit(a_data.dataIterator()); dit.ok(); ++dit) {
          m_gdofs_lhs[dit].setVal(-1);
          box_offset += setDOF( m_gdofs_lhs[dit], 
                                phys_domain,
                                (a_mpi_offset+a_offset+box_offset), 
                                gvec );
        }
      }

      SpaceUtils::copyLevelData( m_gdofs_lhs, m_gdofs_rhs );
      consistentExchange( m_gdofs_rhs );
    
      m_mpi_offset = a_mpi_offset;
      m_is_defined = true;
      return;
    }

  protected:

    int setDOF( T&, const ProblemDomain&, const int, const IntVect& );

    void consistentExchange( LevelData<T>& );

    inline void setDOF( FArrayBox&  a_dof, 
                        const Box&  a_box,
                        const int   a_offset )
    {
      int ncomp( a_dof.nComp() );

      auto bigEnd(a_box.bigEnd()), smallEnd(a_box.smallEnd());
      auto gridSize = bigEnd; gridSize -= smallEnd; gridSize += 1;

      int comp_offset = 0;
      for (int n(0); n < ncomp; n++) {
        BoxIterator bit(a_box);
        for (bit.begin(); bit.ok(); ++bit) {
          auto index = bit();
          int p; _BoxIndexDOF_( p, (index-smallEnd), gridSize, SpaceDim);
          a_dof(index, n) = p + comp_offset + a_offset;
        }
        comp_offset += a_box.numPts();
      }

      return;
    }

    bool m_is_defined;
    int m_mpi_offset;
    LevelData<T> m_gdofs_lhs, m_gdofs_rhs;

  private:

};

class GlobalDOFEMFields
{
  public:
    
    GlobalDOFEMFields()
    {
      m_is_defined = false;
      m_mpi_offset = -1;
    }

    ~GlobalDOFEMFields() { }

    virtual void define(  const int,
                          const LevelData<FluxBox>&,
                          const LevelData<FArrayBox>&,
                          const LevelData<EdgeDataBox>&,
                          const LevelData<NodeFArrayBox>& );

    virtual void define(  const int,
                          const LevelData<EdgeDataBox>&,
                          const LevelData<NodeFArrayBox>& );

    inline bool isDefined() const { return m_is_defined; }

    inline const LevelData<FluxBox>& dofDataLHSB() const
    {
      CH_assert(isDefined());
      return m_gdof_B.dofDataLHS();
    }

    inline const LevelData<FArrayBox>& dofDataLHSBv() const
    {
      CH_assert( isDefined() && (SpaceDim<3) );
      return m_gdof_Bv.dofDataLHS();
    }

    inline const LevelData<EdgeDataBox>& dofDataLHSE() const
    {
      CH_assert(isDefined());
      return m_gdof_E.dofDataLHS();
    }

    inline const LevelData<NodeFArrayBox>& dofDataLHSEv() const
    {
      CH_assert( isDefined() && (SpaceDim<3) );
      return m_gdof_Ev.dofDataLHS();
    }

    inline const LevelData<FluxBox>& dofDataRHSB() const
    {
      CH_assert(isDefined());
      return m_gdof_B.dofDataRHS();
    }

    inline const LevelData<FArrayBox>& dofDataRHSBv() const
    {
      CH_assert( isDefined() && (SpaceDim<3) );
      return m_gdof_Bv.dofDataRHS();
    }

    inline const LevelData<EdgeDataBox>& dofDataRHSE() const
    {
      CH_assert(isDefined());
      return m_gdof_E.dofDataRHS();
    }

    inline const LevelData<NodeFArrayBox>& dofDataRHSEv() const
    {
      CH_assert( isDefined() && (SpaceDim<3) );
      return m_gdof_Ev.dofDataRHS();
    }

  protected:

    bool m_is_defined;
    int m_mpi_offset;

    GlobalDOF<FluxBox> m_gdof_B;
    GlobalDOF<FArrayBox> m_gdof_Bv;

    GlobalDOF<EdgeDataBox> m_gdof_E;
    GlobalDOF<NodeFArrayBox> m_gdof_Ev;

  private:

};

#include "NamespaceFooter.H"

#endif
