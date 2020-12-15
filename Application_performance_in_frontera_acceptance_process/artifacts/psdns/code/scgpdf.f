	subroutine scgpdf (ifc,ilc,i,yp,grx,gry,grz,ithr,nddu,hist,sumll,sumpp)
 
! this routine is extensively modified (02/06/2007)
c
c Changes by PK Yeung, 12/22/08, to enable histogram bins to
c be set dynamically (reading "gpdflim" from input.sc, same
c for all scalars)

#ifndef NOSCALAR
 
	use comsp
	implicit none
	include 'intvars'
!
! standardized pdf of gradients of i-th scalar, with 
! distinction between components parallel and perpendicular 
! to the mean scalar gradient
!
! values in physical space, called from sub. proc2a
!
! also skewness, flatness factor and superskewness
!
	real(b8) grx(nx,zisz),gry(nx,zisz),grz(nx,zisz)
c
	integer i,k,ic,ifc,ilc,lu,iv
	real(b8) diffus,rms1,rms2,terml1,terml2,terml3,
     1           termp1,termp2,termp3,varce0,skew0,flat0,sskew0
!
	integer nddu,ithr
	real(b8) hist(nddu,nc,2),sumll(6,nc),sumpp(6,nc)
!
      real(b8) mean
c
	real sg(3)
!
	integer icall
	save icall
	data icall/0/
!
c
      real bldd,bldd2,brdd,brdd2,udum(1),wdum(1)
      integer para
	real mean0
!
	real, allocatable :: scgll(:),scgp1(:),scgp2(:)
!
!
      character*12 fmt
!
      real nn

      data fmt/'(1p,7e11.4)'/
!
!
      diffus=viscos/pr(i)
 
!*****************************************************************
! Testing grads to determine which component is aligned with the
! mean scalar gradient.
! Assumes mean gradient in only one direction
 
c     para=0     
      para=1     
      do 150 k=1,3
      if (grad(k,i).ne.0) then
        para=k
      endif
 150  continue
#ifdef ROTD
      if (rrate.gt.0.) para=3
#endif

      if (para.eq.0) return
 
      allocate (scgll(nx),scgp1(nx),scgp2(nx))
 
!*******************************************************************

	sg(1)=1./beta1
	sg(2)=1./beta2
	sg(3)=1./beta3

 
 
      if (ifc.eq.1) then
      do 5 k=1,6
      sumpp(k,i)=0.
 5    sumll(k,i)=0.
c
c
	bldd=-gpdflim ; brdd=gpdflim
	bldd2=-gpdflim ; brdd2=gpdflim
c
      call hist1 (1,0,1,udum,wdum,0,1,nddu,bldd,brdd,
     1            hist(1,i,1),fmt)
      call hist1 (1,0,1,udum,wdum,0,1,nddu,bldd2,brdd2,
     1            hist(1,i,2),fmt)
 
      endif
 
	iv=kt(i+3)-kt(4)+1
 
! Sets arrays for one parallel and two perpendicular.
 
c      do 100 zp=1,zisz
      do 100 zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
 
      if (para.eq.1) then
 
! set histogram limits at +/- 10 standard deviations
 
      rms1=sqrt(scgvar(iv,1))
      rms2=sqrt(.5*(scgvar(iv,2)+scgvar(iv,3)))
      do 11 x=1,nx
      scgll(x)=grx(x,zp)*sg(1)
      scgp1(x)=gry(x,zp)*sg(2)
      scgp2(x)=grz(x,zp)*sg(3)
 11   continue
	
!
      else if (para.eq.2) then
      rms1=sqrt(scgvar(iv,2))
      rms2=sqrt(.5*(scgvar(iv,1)+scgvar(iv,3)))
      do 12 x=1,nx
      scgll(x)=gry(x,zp)*sg(2)
      scgp1(x)=grx(x,zp)*sg(1)
      scgp2(x)=grz(x,zp)*sg(3)
 12   continue
!
      else 
      rms1=sqrt(scgvar(iv,3))
      rms2=sqrt(.5*(scgvar(iv,1)+scgvar(iv,2)))
      do 13 x=1,nx
      scgll(x)=grz(x,zp)*sg(3)
      scgp1(x)=grx(x,zp)*sg(1)
      scgp2(x)=gry(x,zp)*sg(2)
 13   continue
!
      end if
!
!
      bldd=-gpdflim*rms1
      brdd=gpdflim*rms1
      bldd2=-gpdflim*rms2
      brdd2=gpdflim*rms2
!
! loop over the current x-z plane
!
      do 20 x=1,nx
!
      terml1=scgll(x)
      terml2=terml1*terml1
      terml3=terml2*terml1
      sumll(1,i)=sumll(1,i)+terml1
      sumll(2,i)=sumll(2,i)+terml2
      sumll(3,i)=sumll(3,i)+terml3
      sumll(4,i)=sumll(4,i)+terml3*terml1
      sumll(5,i)=sumll(5,i)+terml3*terml2
      sumll(6,i)=sumll(6,i)+terml3*terml3
!
      termp1=scgp1(x)
      termp2=termp1*termp1
      termp3=termp2*termp1
      sumpp(1,i)=sumpp(1,i)+termp1
      sumpp(2,i)=sumpp(2,i)+termp2
      sumpp(3,i)=sumpp(3,i)+termp3
      sumpp(4,i)=sumpp(4,i)+termp3*termp1
      sumpp(5,i)=sumpp(5,i)+termp3*termp2
      sumpp(6,i)=sumpp(6,i)+termp3*termp3
!
      termp1=scgp2(x)
      termp2=termp1*termp1
      termp3=termp2*termp1
      sumpp(1,i)=sumpp(1,i)+termp1
      sumpp(2,i)=sumpp(2,i)+termp2
      sumpp(3,i)=sumpp(3,i)+termp3
      sumpp(4,i)=sumpp(4,i)+termp3*termp1
      sumpp(5,i)=sumpp(5,i)+termp3*termp2
      sumpp(6,i)=sumpp(6,i)+termp3*termp3
 
 20   continue
c
! parallel gradient
 
      call hist1 (0,0,1,scgll,wdum,nx,1,nddu,bldd,brdd,
     1            hist(1,i,1),fmt)
      call hist1 (0,0,1,scgp1,wdum,nx,1,nddu,bldd2,brdd2,
     1            hist(1,i,2),fmt)
      call hist1 (0,0,1,scgp2,wdum,nx,1,nddu,bldd2,brdd2,
     1            hist(1,i,2),fmt)
!
 100  continue
 
	deallocate (scgll,scgp1,scgp2)
c
!
      return
	end
	
	subroutine scgpdf_lc (nddu,hist,gradsum)
c
	use comsp
	implicit none
	include 'intvars'
	integer nddu,jthr
	real(b8) hist(nddu,nc,2,0:num_thr-1)
c
	integer lu,ip,ic,i
	real(b8) nn
	real(b8), allocatable :: hist0(:),hist02(:)
      real bldd,bldd2,brdd,brdd2,udum(1),wdum(1)
      real(b8) mean0,gpmom(6),glmom(6)
	real(b8) sum
	real(b8) varce0,skew0,flat0,sskew0

	real(b8) gradsum(6,nc,2,0:num_thr-1)

      character*6 head1,head2
      character*1 tail
      character*7 fn
      character*12 fmt
!
      data head1/'scgpll'/
      data head2/'scgppp'/
      data fmt/'(1p,7e11.4)'/
!
!
c
c sum over the threads first

	if (num_thr.gt.1) then
c
	do ip=1,2
	do ic=1,nc
	do i=1,nddu
	sum=0.
	do jthr=1,num_thr-1
	sum=sum+hist(i,ic,ip,jthr)
	end do
	hist(i,ic,ip,0)=hist(i,ic,ip,0)+sum
	end do
	end do
	end do
c
	do ip=1,2
	do ic=1,nc
	do i=1,6
	sum=0.
	do jthr=1,num_thr-1
	sum=sum+gradsum(i,ic,ip,jthr)
	end do
	gradsum(i,ic,ip,0)=gradsum(i,ic,ip,0)+sum
	end do
	end do
	end do
c
	end if
!
	if (taskid.eq.0) then
	if (jstep.eq.1) then
	do ic=1,nc
        write (tail,"(i1)") ic
        fn=head1//tail
	call fopen1 (69,fn,'formatted  ')
	close (69)
        fn=head2//tail
	call fopen1 (72,fn,'formatted  ')
	close (72)
	end do
	end if
	end if

! collect tally over all tasks
! note double precision for the moments
 
	allocate (hist0(nddu),hist02(nddu))
c
	do 50 i=1,nc

      call MPI_REDUCE (hist(1,i,1,0),hist0,nddu,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
      call MPI_REDUCE (hist(1,i,2,0),hist02,nddu,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
      call MPI_REDUCE (gradsum(1,i,1,0),glmom,6,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
      call MPI_REDUCE (gradsum(1,i,2,0),gpmom,6,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
 
c
! output by task 0
 
      if (taskid.ne.0) go to 50


        write (tail,"(i1)") i
        fn=head1//tail
	call clreop (69,fn,'formatted  ')
	if (jstep.eq.1) rewind (69)
        fn=head2//tail
	call clreop (72,fn,'formatted  ')
	if (jstep.eq.1) rewind (72)
c
! re-scale to obtain standardized pdf
 
      bldd=-gpdflim
      brdd=gpdflim
      bldd2=-gpdflim
      brdd2=gpdflim
!
! if in parallel direction, writes to one logical unit, else writes 
! to a second.
!
! by definition, the gradient fluctuations have zero expected value
!
! output for parallel component
!
c
c beware of integer overflow: nn is declared to be real
c
c     nn=float(nx)*ny*nz
	nn=nx
	nn=nn*ny
	nn=nn*nz

      glmom(:)=glmom(:)/nn
      mean0=glmom(1)
      varce0=glmom(2)
      skew0=glmom(3)/varce0**1.5
      flat0=glmom(4)/varce0**2
      sskew0=glmom(6)/varce0**3
!
	lu=69

	if (kinit(3+i).eq.0.and.istep.eq.1) then
      write (lu,211) 0,0,40,fmt,0.,0.,0.,0.
      write (lu,212) istep-1,time,0.
	else
      write (lu,211) 0,0,40,fmt,mean0,varce0,skew0,flat0
      write (lu,212) istep-1,time,sskew0
	end if
c
      call hist1 (0,lu,1,udum,wdum,0,1,nddu,bldd,brdd,
     1            hist0,fmt)
!
! output for perpendicular component
!
      gpmom(:)=gpmom(:)/nn/2
      mean0=gpmom(1)
      varce0=gpmom(2)
      skew0=gpmom(3)/varce0**1.5
      flat0=gpmom(4)/varce0**2
      sskew0=gpmom(6)/varce0**3
!
	lu=72
	if (kinit(3+i).eq.0.and.istep.eq.1) then
      write (lu,211) 0,0,40,fmt,0.,0.,0.,0.
      write (lu,212) istep-1,time,0.
	else
      write (lu,211) 0,0,40,fmt,mean0,varce0,skew0,flat0
      write (lu,212) istep-1,time,sskew0
	end if
      call hist1 (0,lu,1,udum,wdum,0,1,nddu,bldd2,brdd2,
     1            hist02,fmt)
!
c
 50	continue
!
	deallocate (hist0,hist02)
!
!
 211  format (3i3,3x,a12,3x,'moments:',1p,4e10.3)
 212  format ('scalar grad. step',i6,' time=',1p,e9.3,
     1         5x,' sskew0=',1p,e10.3)
!
#endif
      end
