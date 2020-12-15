      subroutine udxphy (vxy,udxy,udxzr,nv)
c
#ifdef LAG
      use com
c
      implicit none
	include 'intvars'
c
	integer idir
c
#include "fft_stuff.f"

c
	integer nv
	complex(b8) :: vxy(ny,zjsz*xisz)
	complex(b8) :: udxy(ny,zjsz*xisz,nv)
	real(b8) :: udxzr(nx,zisz,yjsz,nv)
c
        complex(b8), allocatable :: bk2i(:)
	complex(b8) bk1i,bk3i
c
	real sum,checkmsq
c
	integer a,ithr,omp_get_thread_num
c
c start from vxy in Fourier space, take derivative in
c ith direction, and transform to vxzr in physical space
c (without changing 
c
        allocate (bk2i(ny))
c
	do y=1,ny
        bk2i(y)=imagi*b22(2)*ky(y)
	end do
c
	ithr=0

c
!$OMP PARALLEL private (ithr,x,xp,z,bk1i,bk3i) shared(bk1,bk3)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
	ithr=0
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 5 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1

         x = xp + xist-1

        bk1i=imagi*bk1(x,2)
	bk3i=imagi*bk3(z,2)
c
	if (nv.eq.3) then
	do y=1,ny
 	udxy(y,a,1)=bk1i*vxy(y,a)
 	udxy(y,a,2)=bk2i(y)*vxy(y,a)
	udxy(y,a,3)=bk3i*vxy(y,a)
	end do
	else
	do y=1,ny
 	udxy(y,a,1)=bk1i*vxy(y,a)
 	udxy(y,a,2)=bk2i(y)*vxy(y,a)
	end do
	end if
c
	call next_xz(xp,z)
c
 5	continue

!$OMP END PARALLEL
c
	deallocate (bk2i)
c
c transform to physicsl space, store result as udxzr
c
	call kxtran (udxy(1,1,1),udxzr(1,1,1,1),udxzr(1,1,1,1),nv)
c
c 3/15/2011: check the mean-squares
c
	if (jstep.eq.1) then
	sum=0.
	do yp=1,yjsz
	do zp=1,zisz
	do x=1,nx
	sum=sum+udxzr(x,zp,yp,1)**2
	end do
	end do
	end do
	call MPI_REDUCE (sum,checkmsq,1,mpireal,MPI_SUM,0,
     1                   MPI_COMM_WORLD,ierr)
	if (taskid.eq.0) then
	checkmsq=checkmsq/nx/ny/nz
	write (6,*) 'istep,check_msq from udxphy:',istep,checkmsq
	end if
	end if

#endif
      return
      end
