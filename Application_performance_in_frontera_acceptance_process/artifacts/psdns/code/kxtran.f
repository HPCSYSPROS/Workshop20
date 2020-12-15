!  3D FFT inverse transform with 2D domain decomposition
!
!  This version uses MPI_Alltoallv to exchange data between processors
!  In the second step, y planes are sent separately
!  The order of array elements in memory is the same in all stages: (x,z,y) 

! Input: XZYg - comlpex array, with y dimension contained entirely,
!               while x and z are block-distributed among processors in 2D grid
! XZgY - complex array (auxiliary), z dimension contained entirely
!               while x and y are block-distributed among processors in 2D grid
! Output: XgZY - an array of real, x dimension is contained entirely within processors memory  
!               while z and y are block-distributed among processors in 2D grid

! !!! CAUTION: In this version: all arrays occupy the same memory space
!

      subroutine kxtran (XgZY,XZgY,XZYg,nv)
      use com
      implicit none
	include 'intvars'
!
#include "fft_stuff.f"

      complex(b8) XgZY(nxpad/2,zisz,yjsz,nv)
      complex(b8) XZgY(nz,xisz,yjsz,nv)
      complex(b8) XZYg(ny,xisz,zjsz,nv)

      real(b8) factor
	integer :: i,nv
c
        integer iy1,iy2,ithr,omp_get_thread_num

c
! Transform in y dimension for all x and z
!
	XZYg(nyhp,:,:,1:nv)=0.

c transform in y direction and switch into z-lines
c
#ifdef HOMOGENEOUS
         call kxcomm1_trans_cyl(XZYg,XZgY,nv)
#else
         call kxcomm1_trans_pen(XZYg,XZgY)
#endif

	call divide_thr (yjsz,iyp1,iypi,'kxtran: yjsz')

! Transform in z dimension for all x, one y-plane at a time
c using FFTW or all y-planes in one call (ESSL)
!
c
	XZgY(nzhp,:,:,1:nv)=cmplx(0.,0.)
c
        ithr=0
!$OMP PARALLEL private(ithr,iy1)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
#ifdef ESSL
	call escpft (XZgY,1,nz,nz,xisz*yjsz/num_thr,1,ithr)
#endif
c
	do 10 i=1,nv

        iy1=iyp1(ithr)
c
c
#ifdef FFTW
            call fftw_execute(plan2_kx(ithr),XZgY(1,1,iy1,i),XZgY(1,1,iy1,i))
#elif defined ESSL
	  call escfft (XZgY(1,1,iy1,i),1,nz,nz,xisz*yjsz/num_thr,1,ithr)
#endif
c
 10	continue


!$OMP END PARALLEL
c
c
c switch into x-lines and take complex-real transform in x
c
#ifdef KXCOMM2
         call kxcomm2_trans (XZgY(1,1,1,1),XgZY(1,1,1,1),nv)
#else
 	do y=1,yjsz
         call kxcomm2t_trans (XZgY(1,1,1,1),XgZY(1,1,1,1),yjsz,y,nv)
	end do
#endif
c
c
      return
      end
