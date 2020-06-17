      subroutine GET_NC_MAPPED_COORDS(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
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
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        xi(i,j,0) = (i + (0.500d0))*dx(0)
                  xi(i,j,1) = (j + (0.500d0))*dx(1)
      enddo
      enddo
      return
      end
      subroutine GET_FC_MAPPED_COORDS(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      double precision offset(0:2 -1)
      offset(0) = (0.500d0)
                offset(1) = (0.500d0)
      offset(dir) = (0.0d0)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        xi(i,j,0) = (i + offset(0))*dx(0)
                  xi(i,j,1) = (j + offset(1))*dx(1)
      enddo
      enddo
      return
      end
      subroutine GET_EC_MAPPED_COORDS(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      double precision offset(0:2 -1)
      offset(0) = (0.0d0)
                offset(1) = (0.0d0)
      offset(dir) = (0.500d0)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        xi(i,j,0) = (i + offset(0))*dx(0)
                  xi(i,j,1) = (j + offset(1))*dx(1)
      enddo
      enddo
      return
      end
