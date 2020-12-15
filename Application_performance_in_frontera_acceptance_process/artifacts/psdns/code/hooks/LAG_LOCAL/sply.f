!      subroutine sply (bs)
      subroutine sply (bs,uy,vy)
c
#ifdef LAG
c
	use mpilag
      implicit none

!      real(p8) bs(bxistart:bxiend,nby)
!      real(p8) bs(bxisize,nby)
      real(p8) bs(xiind,nby)
!	real(p8), allocatable :: uy(:),vy(:)
	real(p8) uy(0:ny),vy(ny)
c
c
      integer x,y,nyp1
c
!      allocate (uy(0:ny),vy(ny))
c
	
      uy(0)=0.
      vy(ny)=0.
      nyp1=ny+1
c
c form nbx y-splines
c
!      do 200 x=bxistart,bxiend
      do 200 x=1,bxisize
c
      do 50 y=1,ny
      uy(y)=(bs(x,y)-uy(y-1))*p(y,2)
 50   continue
c
      do 60 y=ny-1,1,-1
      vy(y)=q(y,2)*vy(y+1)+uy(y)
 60   continue
c
      bs(x,nyp1)=denom(2)*(bs(x,ny)-vy(1)-vy(ny-1))
c
      do 70 y=ny-1,1,-1
      bs(x,y+1)=t(y,2)*bs(x,nyp1)+vy(y)
 70   continue
c
 200  continue
c
!      do x=bxistart,bxiend
      do x=1,bxisize
         bs(x,1)=bs(x,nyp1)
         bs(x,ny+2)=bs(x,2)
         bs(x,ny+3)=bs(x,3)
      end do
c
!       deallocate (uy,vy)
#endif
      return
      end
