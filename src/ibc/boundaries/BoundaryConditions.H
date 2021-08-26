#ifndef  _BOUNDARY_CONDITIONS_H_
#define  _BOUNDARY_CONDITIONS_H_

#include "IntVect.H"
#include "Box.H"
#include "LoHiSide.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

namespace BoundaryConditions {

   void setBC( FArrayBox&       a_this_soln,
         const Box&             a_boundary_box,
         const std::string&     a_bc_type,
         const int&             a_dir,
         const Side::LoHiSide&  a_side );

   void setDirichletBC( FArrayBox&             a_this_soln,
                        const IntVect&         a_ghosts,
                        const Box&             a_boundary_box,
                        const FluxBox&         a_face_values,
                        const int&             a_dir,
                        const Side::LoHiSide&  a_side );

}

#include "NamespaceFooter.H"

#endif
