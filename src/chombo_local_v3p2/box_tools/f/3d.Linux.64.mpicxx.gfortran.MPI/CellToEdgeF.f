      subroutine CELLTOEDGE(
     & cellData
     & ,icellDatalo0,icellDatalo1,icellDatalo2
     & ,icellDatahi0,icellDatahi1,icellDatahi2
     & ,edgeData
     & ,iedgeDatalo0,iedgeDatalo1,iedgeDatalo2
     & ,iedgeDatahi0,iedgeDatahi1,iedgeDatahi2
     & ,iedgeBoxlo0,iedgeBoxlo1,iedgeBoxlo2
     & ,iedgeBoxhi0,iedgeBoxhi1,iedgeBoxhi2
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer icellDatalo0,icellDatalo1,icellDatalo2
      integer icellDatahi0,icellDatahi1,icellDatahi2
      REAL*8 cellData(
     & icellDatalo0:icellDatahi0,
     & icellDatalo1:icellDatahi1,
     & icellDatalo2:icellDatahi2)
      integer iedgeDatalo0,iedgeDatalo1,iedgeDatalo2
      integer iedgeDatahi0,iedgeDatahi1,iedgeDatahi2
      REAL*8 edgeData(
     & iedgeDatalo0:iedgeDatahi0,
     & iedgeDatalo1:iedgeDatahi1,
     & iedgeDatalo2:iedgeDatahi2)
      integer iedgeBoxlo0,iedgeBoxlo1,iedgeBoxlo2
      integer iedgeBoxhi0,iedgeBoxhi1,iedgeBoxhi2
      integer dir
      integer i,j,k
      integer ii,jj,kk
      do k = iedgeBoxlo2,iedgeBoxhi2
      do j = iedgeBoxlo1,iedgeBoxhi1
      do i = iedgeBoxlo0,iedgeBoxhi0
        ii = i-CHF_ID(0,dir)
        jj = j-CHF_ID(1,dir)
        kk = k-CHF_ID(2,dir)
        edgeData(i,j,k) = (0.500d0)*(
     & cellData(ii,jj,kk)
     & + cellData(i,j,k) )
      enddo
      enddo
      enddo
        return
        end
