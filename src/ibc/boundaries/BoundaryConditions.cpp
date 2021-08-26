#include "BoundaryConditions.H"
#include "RealVect.H"

#include "BoundaryConditionsF_F.H"

#include "NamespaceHeader.H"

void BoundaryConditions::setBC( FArrayBox&       a_dst,
                          const Box&             a_boundary_box,
                          const std::string&     a_bc_type,
                          const int&             a_dir,
                          const Side::LoHiSide&  a_side )
{
   CH_TIMERS("BoundaryConditions::setBC()");
   
   // determine if a_dst is stag in dir direction
   const IntVect& stag_vect = a_dst.box().type();
   const int& is_stag = stag_vect[a_dir];
   
   const int ISIDE(a_side);
   IntVect ghosts = a_boundary_box.size();
 
   if (a_bc_type == "even" || a_bc_type == "odd") {
      int evenodd = 1;
      if(a_bc_type=="even") evenodd = 1;
      if(a_bc_type=="odd") evenodd = -1;
      FORT_EVENODD_BC( CHF_FRA(a_dst),
                       CHF_BOX(a_boundary_box),
                       CHF_CONST_INT(evenodd),
                       CHF_CONST_INT(is_stag),
                       CHF_CONST_INT(a_dir),
                       CHF_CONST_INT(ISIDE) );
   }
   
   else if (a_bc_type == "symmetry") {
      CH_assert( !is_stag );
      CH_assert( a_dst.nComp()==SpaceDim );
      FORT_SYMMETRY_BC( CHF_FRA(a_dst),
                        CHF_BOX(a_boundary_box),
  	                CHF_CONST_INTVECT(ghosts),
                        CHF_CONST_INT(a_dir),
                        CHF_CONST_INT(ISIDE) );
   }

   else { // default is fill ghost cells based on extrap
      int order = 2;
      FORT_EXTRAP_BC( CHF_FRA(a_dst),
                      CHF_BOX(a_boundary_box),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(ISIDE),
                      CHF_CONST_INT(order) );
   }  

}

void BoundaryConditions::setDirichletBC( FArrayBox&       a_this_soln,
                                   const IntVect&         a_ghosts,
                                   const Box&             a_boundary_box,
                                   const FluxBox&         a_face_values,
                                   const int&             a_dir,
                                   const Side::LoHiSide&  a_side )

{
   const int ISIDE(a_side);
   FORT_DIRICHLET_BC( CHF_FRA(a_this_soln),
                      CHF_BOX(a_boundary_box),
                      CHF_CONST_INTVECT(a_ghosts),
                      CHF_CONST_FRA(a_face_values[a_dir]),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(ISIDE) );
}

#include "NamespaceFooter.H"
