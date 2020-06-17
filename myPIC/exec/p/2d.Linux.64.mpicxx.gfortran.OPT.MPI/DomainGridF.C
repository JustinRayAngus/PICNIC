#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine GET_NC_MAPPED_COORDS(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
      integer i,j
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        xi(i,j,0) = i*dx(0)
                  xi(i,j,1) = j*dx(1)
      
      enddo
      enddo
      return
      end
      subroutine GET_CC_MAPPED_COORDS(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
      integer i,j
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        xi(i,j,0) = (i + half)*dx(0)
                  xi(i,j,1) = (j + half)*dx(1)
      
      enddo
      enddo
      return
      end
      subroutine GET_FC_MAPPED_COORDS(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dir
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
      integer i,j
      double precision offset(0:CH_SPACEDIM-1)
      offset(0) = half
                offset(1) = half
      offset(dir) = zero
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        xi(i,j,0) = (i + offset(0))*dx(0)
                  xi(i,j,1) = (j + offset(1))*dx(1)
      
      enddo
      enddo
      return
      end
      subroutine GET_EC_MAPPED_COORDS(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dir
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
      integer i,j
      double precision offset(0:CH_SPACEDIM-1)
      offset(0) = zero
                offset(1) = zero
      offset(dir) = half
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        xi(i,j,0) = (i + offset(0))*dx(0)
                  xi(i,j,1) = (j + offset(1))*dx(1)
      
      enddo
      enddo
      return
      end
