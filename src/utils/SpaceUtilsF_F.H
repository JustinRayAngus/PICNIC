#ifndef _SPACEUTILSF_F_H_
#define _SPACEUTILSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDBWENOFACE 
#define GUARDBWENOFACE 
// Prototype for Fortran procedure BWENOFACE ...
//
void FORTRAN_NAME( BWENOFACE ,bwenoface )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_BWENOFACE FORTRAN_NAME( inlineBWENOFACE, inlineBWENOFACE)
#define FORTNT_BWENOFACE FORTRAN_NAME( BWENOFACE, bwenoface)

inline void FORTRAN_NAME(inlineBWENOFACE, inlineBWENOFACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_BWENOFACE");
 FORTRAN_NAME( BWENOFACE ,bwenoface )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDBWENOFACE 

#ifndef GUARDWENO5FACE 
#define GUARDWENO5FACE 
// Prototype for Fortran procedure WENO5FACE ...
//
void FORTRAN_NAME( WENO5FACE ,weno5face )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_CONST_FRA(smoothFn)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_WENO5FACE FORTRAN_NAME( inlineWENO5FACE, inlineWENO5FACE)
#define FORTNT_WENO5FACE FORTRAN_NAME( WENO5FACE, weno5face)

inline void FORTRAN_NAME(inlineWENO5FACE, inlineWENO5FACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_CONST_FRA(smoothFn)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_WENO5FACE");
 FORTRAN_NAME( WENO5FACE ,weno5face )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_CONST_FRA(smoothFn)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDWENO5FACE 

#ifndef GUARDTVDFACE 
#define GUARDTVDFACE 
// Prototype for Fortran procedure TVDFACE ...
//
void FORTRAN_NAME( TVDFACE ,tvdface )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(limiter)
      ,CHFp_CONST_INT(idir) );

#define FORT_TVDFACE FORTRAN_NAME( inlineTVDFACE, inlineTVDFACE)
#define FORTNT_TVDFACE FORTRAN_NAME( TVDFACE, tvdface)

inline void FORTRAN_NAME(inlineTVDFACE, inlineTVDFACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(limiter)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_TVDFACE");
 FORTRAN_NAME( TVDFACE ,tvdface )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(limiter)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDTVDFACE 

#ifndef GUARDUW1FACE 
#define GUARDUW1FACE 
// Prototype for Fortran procedure UW1FACE ...
//
void FORTRAN_NAME( UW1FACE ,uw1face )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_UW1FACE FORTRAN_NAME( inlineUW1FACE, inlineUW1FACE)
#define FORTNT_UW1FACE FORTRAN_NAME( UW1FACE, uw1face)

inline void FORTRAN_NAME(inlineUW1FACE, inlineUW1FACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_UW1FACE");
 FORTRAN_NAME( UW1FACE ,uw1face )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDUW1FACE 

#ifndef GUARDUW1C2FACE 
#define GUARDUW1C2FACE 
// Prototype for Fortran procedure UW1C2FACE ...
//
void FORTRAN_NAME( UW1C2FACE ,uw1c2face )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(cellFun)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(limType)
      ,CHFp_CONST_INT(idir) );

#define FORT_UW1C2FACE FORTRAN_NAME( inlineUW1C2FACE, inlineUW1C2FACE)
#define FORTNT_UW1C2FACE FORTRAN_NAME( UW1C2FACE, uw1c2face)

inline void FORTRAN_NAME(inlineUW1C2FACE, inlineUW1C2FACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA(cellFun)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(limType)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_UW1C2FACE");
 FORTRAN_NAME( UW1C2FACE ,uw1c2face )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA(cellFun)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(limType)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDUW1C2FACE 

#ifndef GUARDUW3FACE 
#define GUARDUW3FACE 
// Prototype for Fortran procedure UW3FACE ...
//
void FORTRAN_NAME( UW3FACE ,uw3face )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_UW3FACE FORTRAN_NAME( inlineUW3FACE, inlineUW3FACE)
#define FORTNT_UW3FACE FORTRAN_NAME( UW3FACE, uw3face)

inline void FORTRAN_NAME(inlineUW3FACE, inlineUW3FACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_UW3FACE");
 FORTRAN_NAME( UW3FACE ,uw3face )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDUW3FACE 

#ifndef GUARDUW5FACE 
#define GUARDUW5FACE 
// Prototype for Fortran procedure UW5FACE ...
//
void FORTRAN_NAME( UW5FACE ,uw5face )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_UW5FACE FORTRAN_NAME( inlineUW5FACE, inlineUW5FACE)
#define FORTNT_UW5FACE FORTRAN_NAME( UW5FACE, uw5face)

inline void FORTRAN_NAME(inlineUW5FACE, inlineUW5FACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_UW5FACE");
 FORTRAN_NAME( UW5FACE ,uw5face )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDUW5FACE 

#ifndef GUARDQUICKFACE 
#define GUARDQUICKFACE 
// Prototype for Fortran procedure QUICKFACE ...
//
void FORTRAN_NAME( QUICKFACE ,quickface )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_QUICKFACE FORTRAN_NAME( inlineQUICKFACE, inlineQUICKFACE)
#define FORTNT_QUICKFACE FORTRAN_NAME( QUICKFACE, quickface)

inline void FORTRAN_NAME(inlineQUICKFACE, inlineQUICKFACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_QUICKFACE");
 FORTRAN_NAME( QUICKFACE ,quickface )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDQUICKFACE 

#ifndef GUARDC2FACETOCELL 
#define GUARDC2FACETOCELL 
// Prototype for Fortran procedure C2FACETOCELL ...
//
void FORTRAN_NAME( C2FACETOCELL ,c2facetocell )(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(facePhi)
      ,CHFp_FRA1(cellPhi) );

#define FORT_C2FACETOCELL FORTRAN_NAME( inlineC2FACETOCELL, inlineC2FACETOCELL)
#define FORTNT_C2FACETOCELL FORTRAN_NAME( C2FACETOCELL, c2facetocell)

inline void FORTRAN_NAME(inlineC2FACETOCELL, inlineC2FACETOCELL)(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(facePhi)
      ,CHFp_FRA1(cellPhi) )
{
 CH_TIMELEAF("FORT_C2FACETOCELL");
 FORTRAN_NAME( C2FACETOCELL ,c2facetocell )(
      CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(facePhi)
      ,CHFt_FRA1(cellPhi) );
}
#endif  // GUARDC2FACETOCELL 

#ifndef GUARDFACE_SCALAR_TO_CELL 
#define GUARDFACE_SCALAR_TO_CELL 
// Prototype for Fortran procedure FACE_SCALAR_TO_CELL ...
//
void FORTRAN_NAME( FACE_SCALAR_TO_CELL ,face_scalar_to_cell )(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_FRA(cellPhi) );

#define FORT_FACE_SCALAR_TO_CELL FORTRAN_NAME( inlineFACE_SCALAR_TO_CELL, inlineFACE_SCALAR_TO_CELL)
#define FORTNT_FACE_SCALAR_TO_CELL FORTRAN_NAME( FACE_SCALAR_TO_CELL, face_scalar_to_cell)

inline void FORTRAN_NAME(inlineFACE_SCALAR_TO_CELL, inlineFACE_SCALAR_TO_CELL)(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA(facePhi)
      ,CHFp_FRA(cellPhi) )
{
 CH_TIMELEAF("FORT_FACE_SCALAR_TO_CELL");
 FORTRAN_NAME( FACE_SCALAR_TO_CELL ,face_scalar_to_cell )(
      CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA(facePhi)
      ,CHFt_FRA(cellPhi) );
}
#endif  // GUARDFACE_SCALAR_TO_CELL 

#ifndef GUARDC2_EDGES_TO_EDGES 
#define GUARDC2_EDGES_TO_EDGES 
// Prototype for Fortran procedure C2_EDGES_TO_EDGES ...
//
void FORTRAN_NAME( C2_EDGES_TO_EDGES ,c2_edges_to_edges )(
      CHFp_BOX(edgeBox)
      ,CHFp_CONST_FRA1(Ein)
      ,CHFp_FRA1(Eout)
      ,CHFp_CONST_INT(idir) );

#define FORT_C2_EDGES_TO_EDGES FORTRAN_NAME( inlineC2_EDGES_TO_EDGES, inlineC2_EDGES_TO_EDGES)
#define FORTNT_C2_EDGES_TO_EDGES FORTRAN_NAME( C2_EDGES_TO_EDGES, c2_edges_to_edges)

inline void FORTRAN_NAME(inlineC2_EDGES_TO_EDGES, inlineC2_EDGES_TO_EDGES)(
      CHFp_BOX(edgeBox)
      ,CHFp_CONST_FRA1(Ein)
      ,CHFp_FRA1(Eout)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_C2_EDGES_TO_EDGES");
 FORTRAN_NAME( C2_EDGES_TO_EDGES ,c2_edges_to_edges )(
      CHFt_BOX(edgeBox)
      ,CHFt_CONST_FRA1(Ein)
      ,CHFt_FRA1(Eout)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDC2_EDGES_TO_EDGES 

#ifndef GUARDC2_NODES_TO_EDGES 
#define GUARDC2_NODES_TO_EDGES 
// Prototype for Fortran procedure C2_NODES_TO_EDGES ...
//
void FORTRAN_NAME( C2_NODES_TO_EDGES ,c2_nodes_to_edges )(
      CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA(ncFin)
      ,CHFp_FRA(ecFout) );

#define FORT_C2_NODES_TO_EDGES FORTRAN_NAME( inlineC2_NODES_TO_EDGES, inlineC2_NODES_TO_EDGES)
#define FORTNT_C2_NODES_TO_EDGES FORTRAN_NAME( C2_NODES_TO_EDGES, c2_nodes_to_edges)

inline void FORTRAN_NAME(inlineC2_NODES_TO_EDGES, inlineC2_NODES_TO_EDGES)(
      CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA(ncFin)
      ,CHFp_FRA(ecFout) )
{
 CH_TIMELEAF("FORT_C2_NODES_TO_EDGES");
 FORTRAN_NAME( C2_NODES_TO_EDGES ,c2_nodes_to_edges )(
      CHFt_BOX(edgeBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA(ncFin)
      ,CHFt_FRA(ecFout) );
}
#endif  // GUARDC2_NODES_TO_EDGES 

#ifndef GUARDC2_NODES_TO_CELLS 
#define GUARDC2_NODES_TO_CELLS 
// Prototype for Fortran procedure C2_NODES_TO_CELLS ...
//
void FORTRAN_NAME( C2_NODES_TO_CELLS ,c2_nodes_to_cells )(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_FRA(ncFin)
      ,CHFp_FRA(ccFout) );

#define FORT_C2_NODES_TO_CELLS FORTRAN_NAME( inlineC2_NODES_TO_CELLS, inlineC2_NODES_TO_CELLS)
#define FORTNT_C2_NODES_TO_CELLS FORTRAN_NAME( C2_NODES_TO_CELLS, c2_nodes_to_cells)

inline void FORTRAN_NAME(inlineC2_NODES_TO_CELLS, inlineC2_NODES_TO_CELLS)(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_FRA(ncFin)
      ,CHFp_FRA(ccFout) )
{
 CH_TIMELEAF("FORT_C2_NODES_TO_CELLS");
 FORTRAN_NAME( C2_NODES_TO_CELLS ,c2_nodes_to_cells )(
      CHFt_BOX(cellBox)
      ,CHFt_CONST_FRA(ncFin)
      ,CHFt_FRA(ccFout) );
}
#endif  // GUARDC2_NODES_TO_CELLS 

#ifndef GUARDC2_CELLS_TO_NODES 
#define GUARDC2_CELLS_TO_NODES 
// Prototype for Fortran procedure C2_CELLS_TO_NODES ...
//
void FORTRAN_NAME( C2_CELLS_TO_NODES ,c2_cells_to_nodes )(
      CHFp_BOX(nodeBox)
      ,CHFp_CONST_FRA(ccFin)
      ,CHFp_FRA(ncFout) );

#define FORT_C2_CELLS_TO_NODES FORTRAN_NAME( inlineC2_CELLS_TO_NODES, inlineC2_CELLS_TO_NODES)
#define FORTNT_C2_CELLS_TO_NODES FORTRAN_NAME( C2_CELLS_TO_NODES, c2_cells_to_nodes)

inline void FORTRAN_NAME(inlineC2_CELLS_TO_NODES, inlineC2_CELLS_TO_NODES)(
      CHFp_BOX(nodeBox)
      ,CHFp_CONST_FRA(ccFin)
      ,CHFp_FRA(ncFout) )
{
 CH_TIMELEAF("FORT_C2_CELLS_TO_NODES");
 FORTRAN_NAME( C2_CELLS_TO_NODES ,c2_cells_to_nodes )(
      CHFt_BOX(nodeBox)
      ,CHFt_CONST_FRA(ccFin)
      ,CHFt_FRA(ncFout) );
}
#endif  // GUARDC2_CELLS_TO_NODES 

#ifndef GUARDEDGE_GRAD_AT_CELLS 
#define GUARDEDGE_GRAD_AT_CELLS 
// Prototype for Fortran procedure EDGE_GRAD_AT_CELLS ...
//
void FORTRAN_NAME( EDGE_GRAD_AT_CELLS ,edge_grad_at_cells )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dX)
      ,CHFp_CONST_FRA1(PhiOnEdges)
      ,CHFp_FRA1(gradPhi) );

#define FORT_EDGE_GRAD_AT_CELLS FORTRAN_NAME( inlineEDGE_GRAD_AT_CELLS, inlineEDGE_GRAD_AT_CELLS)
#define FORTNT_EDGE_GRAD_AT_CELLS FORTRAN_NAME( EDGE_GRAD_AT_CELLS, edge_grad_at_cells)

inline void FORTRAN_NAME(inlineEDGE_GRAD_AT_CELLS, inlineEDGE_GRAD_AT_CELLS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dX)
      ,CHFp_CONST_FRA1(PhiOnEdges)
      ,CHFp_FRA1(gradPhi) )
{
 CH_TIMELEAF("FORT_EDGE_GRAD_AT_CELLS");
 FORTRAN_NAME( EDGE_GRAD_AT_CELLS ,edge_grad_at_cells )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dX)
      ,CHFt_CONST_FRA1(PhiOnEdges)
      ,CHFt_FRA1(gradPhi) );
}
#endif  // GUARDEDGE_GRAD_AT_CELLS 

#ifndef GUARDEDGE_GRAD_AT_NODES 
#define GUARDEDGE_GRAD_AT_NODES 
// Prototype for Fortran procedure EDGE_GRAD_AT_NODES ...
//
void FORTRAN_NAME( EDGE_GRAD_AT_NODES ,edge_grad_at_nodes )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dX)
      ,CHFp_CONST_FRA1(PhiOnEdges)
      ,CHFp_FRA1(gradPhi) );

#define FORT_EDGE_GRAD_AT_NODES FORTRAN_NAME( inlineEDGE_GRAD_AT_NODES, inlineEDGE_GRAD_AT_NODES)
#define FORTNT_EDGE_GRAD_AT_NODES FORTRAN_NAME( EDGE_GRAD_AT_NODES, edge_grad_at_nodes)

inline void FORTRAN_NAME(inlineEDGE_GRAD_AT_NODES, inlineEDGE_GRAD_AT_NODES)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dX)
      ,CHFp_CONST_FRA1(PhiOnEdges)
      ,CHFp_FRA1(gradPhi) )
{
 CH_TIMELEAF("FORT_EDGE_GRAD_AT_NODES");
 FORTRAN_NAME( EDGE_GRAD_AT_NODES ,edge_grad_at_nodes )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dX)
      ,CHFt_CONST_FRA1(PhiOnEdges)
      ,CHFt_FRA1(gradPhi) );
}
#endif  // GUARDEDGE_GRAD_AT_NODES 

#ifndef GUARDC2CELL 
#define GUARDC2CELL 
// Prototype for Fortran procedure C2CELL ...
//
void FORTRAN_NAME( C2CELL ,c2cell )(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(edgePhi)
      ,CHFp_FRA1(cellPhi) );

#define FORT_C2CELL FORTRAN_NAME( inlineC2CELL, inlineC2CELL)
#define FORTNT_C2CELL FORTRAN_NAME( C2CELL, c2cell)

inline void FORTRAN_NAME(inlineC2CELL, inlineC2CELL)(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(edgePhi)
      ,CHFp_FRA1(cellPhi) )
{
 CH_TIMELEAF("FORT_C2CELL");
 FORTRAN_NAME( C2CELL ,c2cell )(
      CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(edgePhi)
      ,CHFt_FRA1(cellPhi) );
}
#endif  // GUARDC2CELL 

#ifndef GUARDEDGE_SCALAR_TO_CELL 
#define GUARDEDGE_SCALAR_TO_CELL 
// Prototype for Fortran procedure EDGE_SCALAR_TO_CELL ...
//
void FORTRAN_NAME( EDGE_SCALAR_TO_CELL ,edge_scalar_to_cell )(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA(edgePhi)
      ,CHFp_FRA(cellPhi) );

#define FORT_EDGE_SCALAR_TO_CELL FORTRAN_NAME( inlineEDGE_SCALAR_TO_CELL, inlineEDGE_SCALAR_TO_CELL)
#define FORTNT_EDGE_SCALAR_TO_CELL FORTRAN_NAME( EDGE_SCALAR_TO_CELL, edge_scalar_to_cell)

inline void FORTRAN_NAME(inlineEDGE_SCALAR_TO_CELL, inlineEDGE_SCALAR_TO_CELL)(
      CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA(edgePhi)
      ,CHFp_FRA(cellPhi) )
{
 CH_TIMELEAF("FORT_EDGE_SCALAR_TO_CELL");
 FORTRAN_NAME( EDGE_SCALAR_TO_CELL ,edge_scalar_to_cell )(
      CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA(edgePhi)
      ,CHFt_FRA(cellPhi) );
}
#endif  // GUARDEDGE_SCALAR_TO_CELL 

#ifndef GUARDC2EDGE 
#define GUARDC2EDGE 
// Prototype for Fortran procedure C2EDGE ...
//
void FORTRAN_NAME( C2EDGE ,c2edge )(
      CHFp_FRA(edgePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_C2EDGE FORTRAN_NAME( inlineC2EDGE, inlineC2EDGE)
#define FORTNT_C2EDGE FORTRAN_NAME( C2EDGE, c2edge)

inline void FORTRAN_NAME(inlineC2EDGE, inlineC2EDGE)(
      CHFp_FRA(edgePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_C2EDGE");
 FORTRAN_NAME( C2EDGE ,c2edge )(
      CHFt_FRA(edgePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_BOX(edgeBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDC2EDGE 

#ifndef GUARDUW1EDGE 
#define GUARDUW1EDGE 
// Prototype for Fortran procedure UW1EDGE ...
//
void FORTRAN_NAME( UW1EDGE ,uw1edge )(
      CHFp_FRA(edgePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_UW1EDGE FORTRAN_NAME( inlineUW1EDGE, inlineUW1EDGE)
#define FORTNT_UW1EDGE FORTRAN_NAME( UW1EDGE, uw1edge)

inline void FORTRAN_NAME(inlineUW1EDGE, inlineUW1EDGE)(
      CHFp_FRA(edgePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_UW1EDGE");
 FORTRAN_NAME( UW1EDGE ,uw1edge )(
      CHFt_FRA(edgePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(edgeBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDUW1EDGE 

#ifndef GUARDC2FACE 
#define GUARDC2FACE 
// Prototype for Fortran procedure C2FACE ...
//
void FORTRAN_NAME( C2FACE ,c2face )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_C2FACE FORTRAN_NAME( inlineC2FACE, inlineC2FACE)
#define FORTNT_C2FACE FORTRAN_NAME( C2FACE, c2face)

inline void FORTRAN_NAME(inlineC2FACE, inlineC2FACE)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_C2FACE");
 FORTRAN_NAME( C2FACE ,c2face )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDC2FACE 

#ifndef GUARDEXTRAP_FOR_CC_OPS 
#define GUARDEXTRAP_FOR_CC_OPS 
// Prototype for Fortran procedure extrap_for_cc_ops ...
//
void FORTRAN_NAME( EXTRAP_FOR_CC_OPS ,extrap_for_cc_ops )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) );

#define FORT_EXTRAP_FOR_CC_OPS FORTRAN_NAME( inlineEXTRAP_FOR_CC_OPS, inlineEXTRAP_FOR_CC_OPS)
#define FORTNT_EXTRAP_FOR_CC_OPS FORTRAN_NAME( EXTRAP_FOR_CC_OPS, extrap_for_cc_ops)

inline void FORTRAN_NAME(inlineEXTRAP_FOR_CC_OPS, inlineEXTRAP_FOR_CC_OPS)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) )
{
 CH_TIMELEAF("FORT_EXTRAP_FOR_CC_OPS");
 FORTRAN_NAME( EXTRAP_FOR_CC_OPS ,extrap_for_cc_ops )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(order)
      ,CHFt_BOX(facebox)
      ,CHFt_BOX(interiorbox)
      ,CHFt_FRA(array) );
}
#endif  // GUARDEXTRAP_FOR_CC_OPS 

#ifndef GUARDEXTRAP_FOR_FC_OPS 
#define GUARDEXTRAP_FOR_FC_OPS 
// Prototype for Fortran procedure extrap_for_fc_ops ...
//
void FORTRAN_NAME( EXTRAP_FOR_FC_OPS ,extrap_for_fc_ops )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) );

#define FORT_EXTRAP_FOR_FC_OPS FORTRAN_NAME( inlineEXTRAP_FOR_FC_OPS, inlineEXTRAP_FOR_FC_OPS)
#define FORTNT_EXTRAP_FOR_FC_OPS FORTRAN_NAME( EXTRAP_FOR_FC_OPS, extrap_for_fc_ops)

inline void FORTRAN_NAME(inlineEXTRAP_FOR_FC_OPS, inlineEXTRAP_FOR_FC_OPS)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(order)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorbox)
      ,CHFp_FRA(array) )
{
 CH_TIMELEAF("FORT_EXTRAP_FOR_FC_OPS");
 FORTRAN_NAME( EXTRAP_FOR_FC_OPS ,extrap_for_fc_ops )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(order)
      ,CHFt_BOX(facebox)
      ,CHFt_BOX(interiorbox)
      ,CHFt_FRA(array) );
}
#endif  // GUARDEXTRAP_FOR_FC_OPS 

#ifndef GUARDFACE_CENTERED_GRAD_COMPONENT 
#define GUARDFACE_CENTERED_GRAD_COMPONENT 
// Prototype for Fortran procedure face_centered_grad_component ...
//
void FORTRAN_NAME( FACE_CENTERED_GRAD_COMPONENT ,face_centered_grad_component )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(grad) );

#define FORT_FACE_CENTERED_GRAD_COMPONENT FORTRAN_NAME( inlineFACE_CENTERED_GRAD_COMPONENT, inlineFACE_CENTERED_GRAD_COMPONENT)
#define FORTNT_FACE_CENTERED_GRAD_COMPONENT FORTRAN_NAME( FACE_CENTERED_GRAD_COMPONENT, face_centered_grad_component)

inline void FORTRAN_NAME(inlineFACE_CENTERED_GRAD_COMPONENT, inlineFACE_CENTERED_GRAD_COMPONENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(grad) )
{
 CH_TIMELEAF("FORT_FACE_CENTERED_GRAD_COMPONENT");
 FORTRAN_NAME( FACE_CENTERED_GRAD_COMPONENT ,face_centered_grad_component )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(var)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_INT(order)
      ,CHFt_FRA1(grad) );
}
#endif  // GUARDFACE_CENTERED_GRAD_COMPONENT 

#ifndef GUARDCELL_CENTERED_GRAD_COMPONENT 
#define GUARDCELL_CENTERED_GRAD_COMPONENT 
// Prototype for Fortran procedure cell_centered_grad_component ...
//
void FORTRAN_NAME( CELL_CENTERED_GRAD_COMPONENT ,cell_centered_grad_component )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(grad) );

#define FORT_CELL_CENTERED_GRAD_COMPONENT FORTRAN_NAME( inlineCELL_CENTERED_GRAD_COMPONENT, inlineCELL_CENTERED_GRAD_COMPONENT)
#define FORTNT_CELL_CENTERED_GRAD_COMPONENT FORTRAN_NAME( CELL_CENTERED_GRAD_COMPONENT, cell_centered_grad_component)

inline void FORTRAN_NAME(inlineCELL_CENTERED_GRAD_COMPONENT, inlineCELL_CENTERED_GRAD_COMPONENT)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_CONST_INT(order)
      ,CHFp_FRA1(grad) )
{
 CH_TIMELEAF("FORT_CELL_CENTERED_GRAD_COMPONENT");
 FORTRAN_NAME( CELL_CENTERED_GRAD_COMPONENT ,cell_centered_grad_component )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(var)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_CONST_INT(order)
      ,CHFt_FRA1(grad) );
}
#endif  // GUARDCELL_CENTERED_GRAD_COMPONENT 

#ifndef GUARDFACE_INTERPOLATE 
#define GUARDFACE_INTERPOLATE 
// Prototype for Fortran procedure face_interpolate ...
//
void FORTRAN_NAME( FACE_INTERPOLATE ,face_interpolate )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(order)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_FRA1(face_var) );

#define FORT_FACE_INTERPOLATE FORTRAN_NAME( inlineFACE_INTERPOLATE, inlineFACE_INTERPOLATE)
#define FORTNT_FACE_INTERPOLATE FORTRAN_NAME( FACE_INTERPOLATE, face_interpolate)

inline void FORTRAN_NAME(inlineFACE_INTERPOLATE, inlineFACE_INTERPOLATE)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(order)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_FRA1(face_var) )
{
 CH_TIMELEAF("FORT_FACE_INTERPOLATE");
 FORTRAN_NAME( FACE_INTERPOLATE ,face_interpolate )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(order)
      ,CHFt_CONST_FRA1(var)
      ,CHFt_FRA1(face_var) );
}
#endif  // GUARDFACE_INTERPOLATE 

#ifndef GUARDSECOND_ORDER_EXTRAPOLATION 
#define GUARDSECOND_ORDER_EXTRAPOLATION 
// Prototype for Fortran procedure second_order_extrapolation ...
//
void FORTRAN_NAME( SECOND_ORDER_EXTRAPOLATION ,second_order_extrapolation )(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_BOX(srcbox)
      ,CHFp_BOX(dstbox)
      ,CHFp_FRA(array) );

#define FORT_SECOND_ORDER_EXTRAPOLATION FORTRAN_NAME( inlineSECOND_ORDER_EXTRAPOLATION, inlineSECOND_ORDER_EXTRAPOLATION)
#define FORTNT_SECOND_ORDER_EXTRAPOLATION FORTRAN_NAME( SECOND_ORDER_EXTRAPOLATION, second_order_extrapolation)

inline void FORTRAN_NAME(inlineSECOND_ORDER_EXTRAPOLATION, inlineSECOND_ORDER_EXTRAPOLATION)(
      CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_BOX(srcbox)
      ,CHFp_BOX(dstbox)
      ,CHFp_FRA(array) )
{
 CH_TIMELEAF("FORT_SECOND_ORDER_EXTRAPOLATION");
 FORTRAN_NAME( SECOND_ORDER_EXTRAPOLATION ,second_order_extrapolation )(
      CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_BOX(srcbox)
      ,CHFt_BOX(dstbox)
      ,CHFt_FRA(array) );
}
#endif  // GUARDSECOND_ORDER_EXTRAPOLATION 

#ifndef GUARDCOPY 
#define GUARDCOPY 
// Prototype for Fortran procedure copy ...
//
void FORTRAN_NAME( COPY ,copy )(
      CHFp_BOX(gridbox)
      ,CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src) );

#define FORT_COPY FORTRAN_NAME( inlineCOPY, inlineCOPY)
#define FORTNT_COPY FORTRAN_NAME( COPY, copy)

inline void FORTRAN_NAME(inlineCOPY, inlineCOPY)(
      CHFp_BOX(gridbox)
      ,CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src) )
{
 CH_TIMELEAF("FORT_COPY");
 FORTRAN_NAME( COPY ,copy )(
      CHFt_BOX(gridbox)
      ,CHFt_FRA(dst)
      ,CHFt_CONST_FRA(src) );
}
#endif  // GUARDCOPY 

#ifndef GUARDVECTOR_NORM 
#define GUARDVECTOR_NORM 
// Prototype for Fortran procedure vector_norm ...
//
void FORTRAN_NAME( VECTOR_NORM ,vector_norm )(
      CHFp_BOX(gridbox)
      ,CHFp_FRA1(dst)
      ,CHFp_CONST_FRA(src) );

#define FORT_VECTOR_NORM FORTRAN_NAME( inlineVECTOR_NORM, inlineVECTOR_NORM)
#define FORTNT_VECTOR_NORM FORTRAN_NAME( VECTOR_NORM, vector_norm)

inline void FORTRAN_NAME(inlineVECTOR_NORM, inlineVECTOR_NORM)(
      CHFp_BOX(gridbox)
      ,CHFp_FRA1(dst)
      ,CHFp_CONST_FRA(src) )
{
 CH_TIMELEAF("FORT_VECTOR_NORM");
 FORTRAN_NAME( VECTOR_NORM ,vector_norm )(
      CHFt_BOX(gridbox)
      ,CHFt_FRA1(dst)
      ,CHFt_CONST_FRA(src) );
}
#endif  // GUARDVECTOR_NORM 

#ifndef GUARDINSPECT_FARRAYBOX 
#define GUARDINSPECT_FARRAYBOX 
// Prototype for Fortran procedure inspect_farraybox ...
//
void FORTRAN_NAME( INSPECT_FARRAYBOX ,inspect_farraybox )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(FA0) );

#define FORT_INSPECT_FARRAYBOX FORTRAN_NAME( inlineINSPECT_FARRAYBOX, inlineINSPECT_FARRAYBOX)
#define FORTNT_INSPECT_FARRAYBOX FORTRAN_NAME( INSPECT_FARRAYBOX, inspect_farraybox)

inline void FORTRAN_NAME(inlineINSPECT_FARRAYBOX, inlineINSPECT_FARRAYBOX)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA1(FA0) )
{
 CH_TIMELEAF("FORT_INSPECT_FARRAYBOX");
 FORTRAN_NAME( INSPECT_FARRAYBOX ,inspect_farraybox )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA1(FA0) );
}
#endif  // GUARDINSPECT_FARRAYBOX 

#ifndef GUARDINSPECT_FLUXBOX 
#define GUARDINSPECT_FLUXBOX 
// Prototype for Fortran procedure inspect_fluxbox ...
//
void FORTRAN_NAME( INSPECT_FLUXBOX ,inspect_fluxbox )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(Flux)
      ,CHFp_CONST_INT(dir) );

#define FORT_INSPECT_FLUXBOX FORTRAN_NAME( inlineINSPECT_FLUXBOX, inlineINSPECT_FLUXBOX)
#define FORTNT_INSPECT_FLUXBOX FORTRAN_NAME( INSPECT_FLUXBOX, inspect_fluxbox)

inline void FORTRAN_NAME(inlineINSPECT_FLUXBOX, inlineINSPECT_FLUXBOX)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(Flux)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_INSPECT_FLUXBOX");
 FORTRAN_NAME( INSPECT_FLUXBOX ,inspect_fluxbox )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(Flux)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDINSPECT_FLUXBOX 

}

#endif
