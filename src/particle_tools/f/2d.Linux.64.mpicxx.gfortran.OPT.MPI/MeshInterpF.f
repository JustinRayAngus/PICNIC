      subroutine NGP_DEPOSIT(
     & rho
     & ,irholo0,irholo1
     & ,irhohi0,irhohi1
     & ,left_edge
     & ,dx
     & ,x
     & ,q
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irholo0,irholo1
      integer irhohi0,irhohi1
      REAL*8 rho(
     & irholo0:irhohi0,
     & irholo1:irhohi1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      REAL*8 q
      integer index(0:2 - 1)
      integer idir
      REAL*8 volume
      volume = dx(0)*dx(1)
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir)) / dx(idir))
      enddo
      rho(index(0), index(1)) =
     & rho(index(0), index(1)) + q / volume
      end
      subroutine NGP_INTERPOLATE(
     & particle_field
     & ,field
     & ,ifieldlo0,ifieldlo1
     & ,ifieldhi0,ifieldhi1
     & ,nfieldcomp
     & ,left_edge
     & ,dx
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 particle_field(0:1)
      integer nfieldcomp
      integer ifieldlo0,ifieldlo1
      integer ifieldhi0,ifieldhi1
      REAL*8 field(
     & ifieldlo0:ifieldhi0,
     & ifieldlo1:ifieldhi1,
     & 0:nfieldcomp-1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      integer index(0:2 - 1)
      integer idir
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir)) / dx(idir))
      enddo
      do idir = 0, 2 - 1
        particle_field(idir) = field(index(0), index(1), idir)
      enddo
      end
      subroutine CIC_DEPOSIT(
     & rho
     & ,irholo0,irholo1
     & ,irhohi0,irhohi1
     & ,left_edge
     & ,dx
     & ,x
     & ,q
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irholo0,irholo1
      integer irhohi0,irhohi1
      REAL*8 rho(
     & irholo0:irhohi0,
     & irholo1:irhohi1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      REAL*8 q
      integer index(0:2 - 1)
      integer idir
      REAL*8 volume, l0,l1
      REAL*8 particle_rho
      REAL*8 weight, w0,w1
      integer ii,jj
      volume = dx(0)*dx(1)
      particle_rho = q / volume
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir) - 0.5d0*dx(idir))
     & / dx(idir))
      enddo
      do ii = index(0), index(0) + 1
        l0 = ii*dx(0) + 0.5d0*dx(0) - x(0) + left_edge(0)
        do jj = index(1), index(1) + 1
          l1 = jj*dx(1) + 0.5d0*dx(1) - x(1) + left_edge(1)
            w0 = 1.d0 - abs(l0 / dx(0))
            w1 = 1.d0 - abs(l1 / dx(1))
            weight = w0 *w1
            rho(ii, jj) =
     & rho(ii, jj) + particle_rho * weight
                enddo
              enddo
      end
      subroutine CIC_INTERPOLATE(
     & particle_field
     & ,field
     & ,ifieldlo0,ifieldlo1
     & ,ifieldhi0,ifieldhi1
     & ,nfieldcomp
     & ,left_edge
     & ,dx
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 particle_field(0:1)
      integer nfieldcomp
      integer ifieldlo0,ifieldlo1
      integer ifieldhi0,ifieldhi1
      REAL*8 field(
     & ifieldlo0:ifieldhi0,
     & ifieldlo1:ifieldhi1,
     & 0:nfieldcomp-1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      integer index(0:2 - 1)
      integer idir
      REAL*8 l0,l1
      REAL*8 weight, w0,w1
      integer ii,jj
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir) - 0.5d0*dx(idir))
     & / dx(idir))
      enddo
      do ii = index(0), index(0) + 1
        l0 = ii*dx(0) + 0.5d0*dx(0) - x(0) + left_edge(0)
        do jj = index(1), index(1) + 1
          l1 = jj*dx(1) + 0.5d0*dx(1) - x(1) + left_edge(1)
            w0 = 1.d0 - abs(l0 / dx(0))
            w1 = 1.d0 - abs(l1 / dx(1))
            weight = w0 *w1
            do idir = 0, 2 - 1
              particle_field(idir) = particle_field(idir) +
     & weight * field(ii, jj, idir)
            enddo
                enddo
              enddo
      end
      subroutine TSC_DEPOSIT(
     & rho
     & ,irholo0,irholo1
     & ,irhohi0,irhohi1
     & ,left_edge
     & ,dx
     & ,x
     & ,q
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irholo0,irholo1
      integer irhohi0,irhohi1
      REAL*8 rho(
     & irholo0:irhohi0,
     & irholo1:irhohi1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      REAL*8 q
      integer index(0:2 - 1)
      integer idir
      REAL*8 volume, l0,l1
      REAL*8 particle_rho
      REAL*8 weight, w0,w1
      integer ii,jj
      volume = dx(0)*dx(1)
      particle_rho = q / volume
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir) - 1.0d0*dx(idir))
     & / dx(idir))
      enddo
      do ii = index(0), index(0) + 2
        l0 = ii*dx(0) + 0.5d0*dx(0) - x(0) + left_edge(0)
        do jj = index(1), index(1) + 2
          l1 = jj*dx(1) + 0.5d0*dx(1) - x(1) + left_edge(1)
            if (abs(l0 / dx(0)) .lt. 0.5d0) then
              w0 = 7.5d-1 - (l0 / dx(0))**2.d0
            else
              w0 = 0.5d0 * (1.5d0 - abs(l0 / dx(0)))**2.d0
            endif
            if (abs(l1 / dx(1)) .lt. 0.5d0) then
              w1 = 7.5d-1 - (l1 / dx(1))**2.d0
            else
              w1 = 0.5d0 * (1.5d0 - abs(l1 / dx(1)))**2.d0
            endif
            weight = w0 *w1
            rho(ii, jj) =
     & rho(ii, jj) + particle_rho * weight
                enddo
              enddo
      end
      subroutine TSC_INTERPOLATE(
     & particle_field
     & ,field
     & ,ifieldlo0,ifieldlo1
     & ,ifieldhi0,ifieldhi1
     & ,nfieldcomp
     & ,left_edge
     & ,dx
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 particle_field(0:1)
      integer nfieldcomp
      integer ifieldlo0,ifieldlo1
      integer ifieldhi0,ifieldhi1
      REAL*8 field(
     & ifieldlo0:ifieldhi0,
     & ifieldlo1:ifieldhi1,
     & 0:nfieldcomp-1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      integer index(0:2 - 1)
      integer idir
      REAL*8 volume, l0,l1
      REAL*8 weight, w0,w1
      integer ii,jj
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir) - 1.0d0*dx(idir))
     & / dx(idir))
      enddo
      do ii = index(0), index(0) + 2
        l0 = ii*dx(0) + 0.5d0*dx(0) - x(0) + left_edge(0)
        do jj = index(1), index(1) + 2
          l1 = jj*dx(1) + 0.5d0*dx(1) - x(1) + left_edge(1)
            if (abs(l0 / dx(0)) .lt. 0.5d0) then
              w0 = 7.5d-1 - (l0 / dx(0))**2.d0
            else
              w0 = 0.5d0 * (1.5d0 - abs(l0 / dx(0)))**2.d0
            endif
            if (abs(l1 / dx(1)) .lt. 0.5d0) then
              w1 = 7.5d-1 - (l1 / dx(1))**2.d0
            else
              w1 = 0.5d0 * (1.5d0 - abs(l1 / dx(1)))**2.d0
            endif
            weight = w0 *w1
            do idir = 0, 2 - 1
              particle_field(idir) = particle_field(idir) +
     & weight * field(ii, jj, idir)
            enddo
                enddo
              enddo
      end
      subroutine W4_DEPOSIT(
     & rho
     & ,irholo0,irholo1
     & ,irhohi0,irhohi1
     & ,left_edge
     & ,dx
     & ,x
     & ,q
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irholo0,irholo1
      integer irhohi0,irhohi1
      REAL*8 rho(
     & irholo0:irhohi0,
     & irholo1:irhohi1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      REAL*8 q
      integer index(0:2 - 1)
      integer idir
      REAL*8 volume, l0,l1
      REAL*8 particle_rho
      REAL*8 weight, w0,w1
      integer ii,jj
      volume = dx(0)*dx(1)
      particle_rho = q / volume
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir) - 1.5d0*dx(idir))
     & / dx(idir))
      enddo
      do ii = index(0), index(0) + 3
        l0 = ii*dx(0) + 0.5d0*dx(0) - x(0) + left_edge(0)
        do jj = index(1), index(1) + 3
          l1 = jj*dx(1) + 0.5d0*dx(1) - x(1) + left_edge(1)
            if (abs(l0 / dx(0)) .lt. 1.d0) then
              w0 = 1.d0 - 5.d0*abs(l0 / dx(0))**2.d0 / 2.d0 + 3.d0*abs(l
     &0 / dx(0))**3.d0 / 2.d0
            else
              w0 = 0.5d0 * (2.0d0 - abs(l0 / dx(0)))**2.d0 * (1.d0 - abs
     &(l0 / dx(0)))
            endif
            if (abs(l1 / dx(1)) .lt. 1.d0) then
              w1 = 1.d0 - 5.d0*abs(l1 / dx(1))**2.d0 / 2.d0 + 3.d0*abs(l
     &1 / dx(1))**3.d0 / 2.d0
            else
              w1 = 0.5d0 * (2.0d0 - abs(l1 / dx(1)))**2.d0 * (1.d0 - abs
     &(l1 / dx(1)))
            endif
            weight = w0 *w1
            rho(ii, jj) =
     & rho(ii, jj) + particle_rho * weight
                   enddo
                 enddo
      end
      subroutine W4_INTERPOLATE(
     & particle_field
     & ,field
     & ,ifieldlo0,ifieldlo1
     & ,ifieldhi0,ifieldhi1
     & ,nfieldcomp
     & ,left_edge
     & ,dx
     & ,x
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 particle_field(0:1)
      integer nfieldcomp
      integer ifieldlo0,ifieldlo1
      integer ifieldhi0,ifieldhi1
      REAL*8 field(
     & ifieldlo0:ifieldhi0,
     & ifieldlo1:ifieldhi1,
     & 0:nfieldcomp-1)
      REAL*8 left_edge(0:1)
      REAL*8 dx(0:1)
      REAL*8 x(0:1)
      integer index(0:2 - 1)
      integer idir
      REAL*8 volume, l0,l1
      REAL*8 weight, w0,w1
      integer ii,jj
      do idir = 0, 2 - 1
        index(idir) = floor((x(idir) - left_edge(idir) - 1.5d0*dx(idir))
     & / dx(idir))
      enddo
      do ii = index(0), index(0) + 3
        l0 = ii*dx(0) + 0.5d0*dx(0) - x(0) + left_edge(0)
        do jj = index(1), index(1) + 3
          l1 = jj*dx(1) + 0.5d0*dx(1) - x(1) + left_edge(1)
            if (abs(l0 / dx(0)) .lt. 1.d0) then
              w0 = 1.d0 - 5.d0*abs(l0 / dx(0))**2.d0 / 2.d0 + 3.d0*abs(l
     &0 / dx(0))**3.d0 / 2.d0
            else
              w0 = 0.5d0 * (2.0d0 - abs(l0 / dx(0)))**2.d0 * (1.d0 - abs
     &(l0 / dx(0)))
            endif
            if (abs(l1 / dx(1)) .lt. 1.d0) then
              w1 = 1.d0 - 5.d0*abs(l1 / dx(1))**2.d0 / 2.d0 + 3.d0*abs(l
     &1 / dx(1))**3.d0 / 2.d0
            else
              w1 = 0.5d0 * (2.0d0 - abs(l1 / dx(1)))**2.d0 * (1.d0 - abs
     &(l1 / dx(1)))
            endif
            weight = w0 *w1
            do idir = 0, 2 - 1
              particle_field(idir) = particle_field(idir) +
     & weight * field(ii, jj, idir)
            enddo
                enddo
              enddo
      end
