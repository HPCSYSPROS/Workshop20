          subroutine phase1_om (uy,omy,m)
c
#ifdef VORPHY
            use comp
            implicit none
            include 'intvars'


#include "fft_stuff.f"

            complex(b8) :: uy(nypad,zjsz*xisz,3)
            complex(b8) :: omy(nypad,zjsz*xisz,3)
            integer :: i,j,m
            complex(b8), allocatable :: shift(:)
            integer iy,a

        integer ithr,omp_get_thread_num,ia
	complex bk1i,bk3i

	if (nc.gt.0) return

!$OMP PARALLEL private (ithr,x,xp,z,bk1i,bk2i,bk3i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif


	omy(:,:,:)=0.
c
      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 10 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c

         x = xp + xist-1
	bk1i=imagi*b11(m)*kx(x)
	bk3i=imagi*b33(m)*kz(z)
	
         
c form vorticity in physical space (operate on
c already phase-shifted velocities)
!
            do 30 y=1,nypad
               if(mask(y,a)) then
                  bk2i=imagi*b22(m)*ky(y)
 		omy(y,a,1)=bk2i*uy(y,a,3)-bk3i*uy(y,a,2)
 		omy(y,a,2)=bk3i*uy(y,a,1)-bk1i*uy(y,a,3)
 		omy(y,a,3)=bk1i*uy(y,a,2)-bk2i*uy(y,a,1)
		end if
 30         continue

            call next_xz(xp,z)
 10	continue

!$OMP END PARALLEL

!
#endif
        return 
      end subroutine phase1_om
