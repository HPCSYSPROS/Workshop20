      subroutine xktran (Source,buf1,Dest,nv)
      use com
      implicit none
#include "intvars"

#include "fft_stuff.f"

      real(b8) Source(nxpad,zisz,yjsz,nv)
      complex(b8) buf1(nzpad,xisz,yjsz,nv)
#ifdef HOMOGENEOUS
      complex(b8) Dest(nypad,num_al,nv)
#else
      complex(b8) Dest(nypad,zjsz,xisz)
#endif
      complex(b8), allocatable :: buf(:,:),buf2(:,:,:)
      real xnorm,ynorm,znorm,tmpvmaxz,vel
      
      real(b8) factor,factor_yz
      integer i,m,iz,ix,x2,iy,a
	integer :: nv
c
	integer ii
c
        integer iy1,iy2,ithr,omp_get_thread_num

c	call init_work

      factor_yz=1./real(nypad*nzpad,b8)
      factor = factor_yz / nxpad

c take real-to-complex transform in x, then switch to x-lines
c
         call xkcomm1_trans(Source,buf1,nv)
c
        ithr=0
c
!$OMP PARALLEL private(ithr,iy1,iy2,x,y,z,i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
#ifdef ESSL
        call cft(1,buf1,1,nzpad,buf1,1,nzpad,nzpad,xisz*yjsz/num_thr,
     1       1,factor_yz, 
     1        caux1(1,ithr),cnaux,caux2(1,ithr),cnaux)                 
#endif


c     Perform FFT in z for all x for a given y plane and normalize
c
	do 10 i=1,nv
c
        iy1=iyp1(ithr)
        iy2=iyp1(ithr)+iypi(ithr)-1
  
#ifdef FFTW
      call fftw_execute(plan2_xk(ithr),buf1(1,1,iy1,i),buf1(1,1,iy1,i))
      do y=iy1,iy2
         do x=1,xisz
            do z=1,nzpad
               buf1(z,x,y,i) = buf1(z,x,y,i) * factor_yz
            enddo
         enddo
      enddo
#elif defined ESSL
        call cft(0,buf1(1,1,iy1,i),1,nzpad,buf1(1,1,iy1,i),1,nzpad,nzpad,
     1            xisz*yjsz/num_thr,1,factor_yz, 
     1           caux1(1,ithr),cnaux,caux2(1,ithr),cnaux)                 
#endif
 
 10	continue
c
!$OMP END PARALLEL

c switch to y-lines and take FFT in y
              
         if(num_al_i(0) .gt. 0) then
#ifdef HOMOGENEOUS
            call xkcomm2_trans_cyl(buf1,Dest,nv)
#else
            call xkcomm2_trans_pen(buf1,Dest,1)
#endif
         endif
c
c	call free_work
c
      return
      end

