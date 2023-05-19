#include <vector>
#include "GlobalDOF.H"

#include "NamespaceHeader.H"

template <>
int GlobalDOF<FArrayBox>::setDOF( FArrayBox&            a_dof,
                                  const ProblemDomain&,
                                  const int             a_offset,
                                  const IntVect&        a_gvec )
{
  auto ncomp( a_dof.nComp() );
  auto box( grow( a_dof.box(), -1*a_gvec ) );

  setDOF( a_dof, box, a_offset );

  return box.numPts()*ncomp;
}

template <>
int GlobalDOF<FluxBox>::setDOF( FluxBox&              a_dof,
                                const ProblemDomain&  a_phys_domain,
                                const int             a_offset,
                                const IntVect&        a_gvec )
{
  auto ncomp( a_dof.nComp() );
  int offset(0);

  for (int dir = 0; dir < SpaceDim; dir++) {

    auto box( grow( a_dof[dir].box(), -1*a_gvec ) );
    int idir_bdry = a_phys_domain.domainBox().bigEnd(dir);
    if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
  
    setDOF( a_dof[dir], box, a_offset+offset );
    offset += box.numPts()*ncomp;
  }

  return offset;
}

template <>
int GlobalDOF<EdgeDataBox>::setDOF( EdgeDataBox&          a_dof,
                                    const ProblemDomain&  a_phys_domain,
                                    const int             a_offset,
                                    const IntVect&        a_gvec )
{
  auto ncomp( a_dof.nComp() );
  int offset(0);

  for (int dir = 0; dir < SpaceDim; dir++) {

    auto box( grow( a_dof[dir].box(), -1*a_gvec ) );
    for (int adir=0; adir<SpaceDim; ++adir) {
       if (adir != dir) {
          int idir_bdry = a_phys_domain.domainBox().bigEnd(adir);
          if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
       }
    }
  
    setDOF( a_dof[dir], box, a_offset+offset );
    offset += box.numPts()*ncomp;
  }

  return offset;
}

template <>
int GlobalDOF<NodeFArrayBox>::setDOF( NodeFArrayBox&        a_dof,
                                      const ProblemDomain&  a_phys_domain,
                                      const int             a_offset,
                                      const IntVect&        a_gvec )
{
  auto ncomp( a_dof.nComp() );
  auto box( surroundingNodes( grow( a_dof.box(), -1*a_gvec ) ) );
  for (int dir=0; dir<SpaceDim; ++dir) {
     int idir_bdry = a_phys_domain.domainBox().bigEnd(dir);
     if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
  }

  setDOF( a_dof.getFab(), box, a_offset );

  return box.numPts()*ncomp;
}

template <>
void GlobalDOF<FArrayBox>::consistentExchange( LevelData<FArrayBox>& a_data )
{
  a_data.exchange();
}

template <>
void GlobalDOF<FluxBox>::consistentExchange( LevelData<FluxBox>& a_data )
{
  SpaceUtils::exchangeFluxBox( a_data );
}

template <>
void GlobalDOF<EdgeDataBox>::consistentExchange( LevelData<EdgeDataBox>& a_data )
{
  SpaceUtils::exchangeEdgeDataBox( a_data );
}

template <>
void GlobalDOF<NodeFArrayBox>::consistentExchange( LevelData<NodeFArrayBox>& a_data )
{
  SpaceUtils::exchangeNodeFArrayBox( a_data );
}

void GlobalDOFEMFields::define( const int                       a_vec_size,
                                const LevelData<FluxBox>&       a_B,
                                const LevelData<FArrayBox>&     a_Bv,
                                const LevelData<EdgeDataBox>&   a_E,
                                const LevelData<NodeFArrayBox>& a_Ev )
{
  CH_assert(!isDefined());

  int rank, nproc;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  rank = 0;
  nproc = 1;
#endif

  /* Compute the MPI offset */
  {
    std::vector<int> local_npts(nproc, 0);
    local_npts[rank] = a_vec_size;
#ifdef CH_MPI
    MPI_Allreduce(  MPI_IN_PLACE,
                    local_npts.data(),
                    nproc,
                    MPI_INT,
                    MPI_SUM,
                    MPI_COMM_WORLD );
#endif
    m_mpi_offset = 0;
    for (int i=0; i<rank; i++) m_mpi_offset += local_npts[i];
  }

  int local_offset(0);
  m_gdof_B.define( local_offset, m_mpi_offset, a_B );
  local_offset += SpaceUtils::nDOF( a_B );
  if (SpaceDim < 3) {
    m_gdof_Bv.define( local_offset, m_mpi_offset, a_Bv );
    local_offset += SpaceUtils::nDOF( a_Bv );
  }
  m_gdof_E.define( local_offset, m_mpi_offset, a_E );
  local_offset += SpaceUtils::nDOF( a_E );
  if (SpaceDim < 3) {
    m_gdof_Ev.define( local_offset, m_mpi_offset, a_Ev );
    local_offset += SpaceUtils::nDOF( a_Ev );
  }

  if (local_offset != a_vec_size) {
    printf("Error in GlobalDOF::define() on rank %d", rank);
    printf(": mismatch in vector size and offset counting.\n");
  }
  CH_assert(local_offset == a_vec_size);

  m_is_defined = true;
  return;
}

void GlobalDOFEMFields::define( const int                       a_vec_size,
                                const LevelData<EdgeDataBox>&   a_E,
                                const LevelData<NodeFArrayBox>& a_Ev )
{
  CH_assert(!isDefined());

  int rank, nproc;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  rank = 0;
  nproc = 1;
#endif

  /* Compute the MPI offset */
  {
    std::vector<int> local_npts(nproc, 0);
    local_npts[rank] = a_vec_size;
#ifdef CH_MPI
    MPI_Allreduce(  MPI_IN_PLACE,
                    local_npts.data(),
                    nproc,
                    MPI_INT,
                    MPI_SUM,
                    MPI_COMM_WORLD );
#endif
    m_mpi_offset = 0;
    for (int i=0; i<rank; i++) m_mpi_offset += local_npts[i];
  }

  int local_offset(0);
  m_gdof_E.define( local_offset, m_mpi_offset, a_E );
  local_offset += SpaceUtils::nDOF( a_E );
  if (SpaceDim < 3) {
    m_gdof_Ev.define( local_offset, m_mpi_offset, a_Ev );
    local_offset += SpaceUtils::nDOF( a_Ev );
  }

  if (local_offset != a_vec_size) {
    printf("Error in GlobalDOF::define() on rank %d", rank);
    printf(": mismatch in vector size and offset counting.\n");
  }
  CH_assert(local_offset == a_vec_size);

  m_is_defined = true;
  return;
}

#include "NamespaceFooter.H"
