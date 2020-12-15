      subroutine scnpdf (vx,i)
 
#ifndef NOSCALAR
	use comsp
	include 'intvars'
 
! phymom must be called first in order to get the variance
!
! routine for forming sample histogram of scalar values
! in physical space, called from sub. predic
! the data are scaled to produce the pdf of a
! standardised random variable
! currently assuming only 1 value of i will be used
!
! extreme values used are -2. to 2., as in e&p fda-87-12
!
! revised to use call to phymom by P.K. 6/92
 
!     parameter (nddu=101)
      parameter (nddu=81)
 
	real(b8) :: vx(nx,zisz,yjsz)
 
	real(b8), allocatable :: rv(:)
	integer, allocatable :: ist(:)
	save ist
 
      dimension hist(nddu)
      real udum(1),wdum(1),hist0(nddu)
      character*6 head
      character*1 tail
      character*7 fn
      character*12 fmt
 
      data head/'scnpdf'/
      data fmt/'(1p,7e11.4)'/
 
	if (istep.eq.1.and.kinit(i+3).eq.0) go to 99
 
      if (.not.allocated(ist).and.i.eq.1) then
      allocate (ist(ncd))
      ist(:)=0
      end if
 
	allocate (rv(nx))
 

! values of the ith scalar is contained in ur(x,i+3+nc,y,z)
 
      j=i+3+nc
 
! extreme limits are +/- 5 std dev
 
      bldd=-5.
      brdd=5.
!
      call hist1 (1,0,1,udum,wdum,0,1,nddu,bldd,brdd,hist,fmt)
!
! note that the expectation of the scalar is always zero
!
! compute the moments
!
      stddev=sqrt(varce(j))
      denom=1./stddev
 
      do 20 yp=1,yjsz
      do 20 zp=1,xisz
 
      do 25 x=1,nx
 25   rv(x)=vx(x,zp,yp)*denom
      call hist1 (0,0,1,rv,wdum,nx,1,nddu,bldd,brdd,hist,fmt)
 
 20   continue
 
! collect tally over all tasks
 
      call MPI_REDUCE (hist,hist0,nddu,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
 
	deallocate (rv)
 
      if (taskid.ne.0) return
 
! open the file on the first time
 
      lu=luscn+(i-1)*6
      if (ist(i).eq.0) then
      ist(i)=1
        write (tail,"(i1)") i
      fn=head//tail
      call fopen1 (lu,fn,'formatted  ')
      end if
!
      write (lu,211) 0,0,40,fmt,skcof(j),flat(j)
      write (lu,212) i,istep-1,time
      call hist1 (0,lu,1,udum,wdum,0,1,nddu,bldd,brdd,hist0,fmt)
 
 211  format (3i3,3x,a12,5x,'skew,flat=',1p,2e10.3)
 212  format ('std pdf scalar#',i1,' step',i5,' time=',1p,e9.3)
 
#endif
 99   return
      end
