#ifndef _BC_UTILS_H_
#define _BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
#include "Vector.H"

#include "NamespaceHeader.H"

namespace BCUtils {

   inline
   Box getGhostBox( const Box&             a_domain_box,
                    const Box&             a_box,
                    const int&             a_dir,
                    const Side::LoHiSide&  a_side,
                    const int&             a_num_ghosts )
   {
      const Box g_box( adjCellBox( a_box, a_dir, a_side, a_num_ghosts ) );
      const Box g_domain_box( adjCellBox( a_domain_box, a_dir, a_side, a_num_ghosts ) );
      return (g_box & g_domain_box);
   }

   inline
   bool touchesPhysicalBoundary( const Box&             a_domain_box,
                                 const Box&             a_box,
                                 const int&             a_dir,
                                 const Side::LoHiSide&  a_side )
   {
      const Box g_box( adjCellBox( a_box, a_dir, a_side, 1 ) );
      const Box g_domain_box( adjCellBox( a_domain_box, a_dir, a_side, 1 ) );
      return g_box.intersects( g_domain_box );
   }


}

#include "NamespaceFooter.H"

#endif
