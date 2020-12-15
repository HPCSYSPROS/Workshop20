      subroutine phypdfkrsp (vx,lu,title,a,b,hist0,nddu,binvi)
c
c parallelized variable-interval version
c special choice of sampling bins (suggested by Sreeni)
c search algorithm suggested by Steve Pope
c
c PDF of (vx-a)/b
c
	use mpicom
	implicit none
	include 'intvars'
c
        real(b8) :: vx(nx,zisz,yjsz)
c
      real(b8) hist0(nddu)
      character*(*) title
c
	integer nddu,lu
	real(b8) udum,wdum,a,b
c
      real(b8), allocatable :: hist(:),work(:)
c
      real(b8) binvi(nddu)
c
      character*12 fmt
      data fmt/'(1p,7e11.4)'/
c
c initialize the histogram, with end points already selected
c

      allocate (hist(nddu))
      allocate (work(nx))
c
      call hist1vip (1,0,1,udum,wdum,0,1,nddu,binvi,hist,fmt)
c
      do 20 yp=1,yjsz
      do 20 zp=1,zisz
      do 25 x=1,nx
      work(x)=(vx(x,zp,yp)-a)/b
 25   continue
      call hist1vip (0,0,1,work,wdum,nx,1,nddu,binvi,hist,fmt)
 20   continue
c
c 
c collect tally over all tasks
c
      call MPI_REDUCE (hist,hist0,nddu,mpireal,MPI_SUM,
     1                 0,MPI_COMM_WORLD,ierr)
c
      if (taskid.eq.0) then
      write (lu,211) 0,0,40,fmt
      write (lu,212) title
      call hist1vip (0,lu,1,udum,wdum,0,1,nddu,binvi,hist0,fmt)
	write (6,"('phypdfkrsp:',i5,1p,2e12.4,4x,a40)") nddu,binvi(1),hist0(1),title
      end if
      call hist1vip (-1,0,1,udum,wdum,0,1,nddu,binvi,hist0,fmt)
c
      deallocate (hist,work)
c
 211  format (3i3,3x,a12)
 212  format (a40)
c
      return
      end
