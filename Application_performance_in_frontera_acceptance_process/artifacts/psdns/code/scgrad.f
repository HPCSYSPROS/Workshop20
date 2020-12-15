      subroutine scgrad (utx,ux,yp,scgy,ithr)
c
c Fixes made by P.K Yeung, 10/27/09
 
#ifndef NOSCALAR
	use comsp
	implicit none
	include 'intvars'
 
	integer i,j,k,ij
	real factor, sum

! calculate mean covariance matrix of gradients of different scalars
 
! called from sub. proc2a, once for each y plane
 
	real(b8) :: ux(nx,zisz,yjsz,nu)
	real(b8) :: utx(nx,zisz,nut)
 
	integer ithr,jthr
c
 	real(b8) scgy(ncgd,3,0:num_thr-1)
c
	if (nc.eq.0) return
c
c	if (yp.eq.0.and.ithr.eq.0) go to 2
	if (yp.eq.0) then
	if (ithr.eq.0) go to 2
	return
	end if
 
 1    continue
 
! x-gradients, stored in ur(x,i+3,y,z)
c
      do 10 zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
      ij=0
      do 10 i=1,nc
      do 10 j=i,nc
      ij=ij+1
      do 15 x=1,nx
      scgy(ij,1,ithr)=scgy(ij,1,ithr)+ux(x,zp,yp,i+3)*ux(x,zp,yp,j+3)
 15   continue
 10   continue
 
! y-gradients, stored in ur(x,i+3+nc,y,z)
 
      do 20 zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
      ij=0
      do 20 i=1,nc
      do 20 j=i,nc
      ij=ij+1
      do 25 x=1,nx
      scgy(ij,2,ithr)=scgy(ij,2,ithr)+ux(x,zp,yp,i+3+nc)*ux(x,zp,yp,j+3+nc)
 25   continue
 20   continue
 
! z-gradients, stored in utr(x,i,z)
 
      ij=0
      do 30 i=1,nc
      do 30 j=i,nc
      ij=ij+1
      do 30 zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
      do 35 x=1,nx
      scgy(ij,3,ithr)=scgy(ij,3,ithr)+utx(x,zp,i)*utx(x,zp,j)
 35   continue
 30   continue
 
      return
 
 2    continue

	if (num_thr.gt.1) then
        do i=1,3
	do ij=1,nc*(nc+1)/2
        sum=0.
        do jthr=0,num_thr-1
        sum=sum+scgy(ij,i,jthr)
        end do
        scgy(ij,i,0)=sum
        end do
        end do
	end if

 
! ilc=1: form global moments by a reduction operation across all tasks,
!        placing the result in task 0, which proceeds to write info.
 
      call MPI_ALLREDUCE (scgy,scgcov,ncgd*3,mpireal,MPI_SUM,
     1                 MPI_COMM_WORLD,mpierr)
c
! average over all nodal points
 
      factor=1./nx/ny/nz
      do 40 j=1,3
      do 40 i=1,ncgd
      scgcov(i,j)=scgcov(i,j)*factor
 40   continue

        scgcov(:,1)=scgcov(:,1)/b11(2)
        scgcov(:,2)=scgcov(:,2)/b22(2)
        scgcov(:,3)=scgcov(:,3)/b33(2)

      ij=0
      do i=1,nc
      do j=i,nc
      ij=ij+1
      if (i.eq.j) then
      do k=1,3
      scgmsq(i,k)=scgcov(ij,k)
      end do
      end if
      end do
      end do
 
#ifdef TEST
c tests done by PK Yeung, 1/2/09
c
        if (taskid.eq.0) then
        do i=1,nc
        check=2.*viscos/pr(i)*(scgmsq(i,1)+scgmsq(i,2)+scgmsq(i,3))
        if (i.eq.1) write (6,901) istep,kstep,ioflag
 901    format ('scgrad: istep,kstep,ioflag=',3i3)
        write (6,902) i, scgmsq(i,1),scgmsq(i,2),scgmsq(i,3)
        write (6,902) i, (2.*corr(kt(i+3))/taylor(i+3,i+3,j)**2,j=1,3)
        write (6,903) check, scdiss(i)
 902    format ('scgrad: i,scgmsq=',i3,1p,3e14.6)
 903    format ('scgrad: check scalar_diss', 1p,2e14.6)
        end do
        end if
#endif
c

      if (taskid.ne.0) return
 
      call scgout
 
#endif
      return
      end
