          subroutine phase1 (uy,m)
c
            use comp
            implicit none
            include 'intvars'

#include "fft_stuff.f"

            complex(b8) :: uy(nypad,zjsz*xisz,nu)
            integer :: i,j,m
            complex(b8), allocatable :: shift(:)
            complex(b8) :: syz,sxz
            integer iy,a,xz(2)
            real(b8) factor,rk1,rk1sq

        integer ithr,omp_get_thread_num,ia

      allocate (shift(nypad))
c
#ifdef LVSTEP
	if (mod(jstep-1,ivstep).eq.0.or.
     1      (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) then
#endif
c	call phshift (uy(1,1,1),uy(1,1,1),nu,3,kstep)
 	call phshift_inplace (uy(1,1,1),3,kstep)
#ifdef LVSTEP
	end if
#endif

!$OMP PARALLEL private (ithr,sxz,bk2i,x,xp,z)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif


      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 10 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c

         x = xp + xist-1
         
         sxz=sz(z,kstep)*sx(x,kstep)

! shift scalars and form y gradients
!
#ifndef NOSCALAR
         do 30 j=4,3+nc
            do 30 y=1,nypad
               if(mask(y,a)) then
               shift(y) = sxz *sy(y,kstep)
                  bk2i=imagi*b22(m)*ky(y)
                  
                  uy(y,a,j)    = shift(y) * uy(y,a,j)
                  uy(y,a,j+nc) =  bk2i * uy(y,a,j)
               else
                  uy(y,a,j) = 0.0
                  uy(y,a,j+nc) = 0.0
               endif
 30         continue
#endif

            call next_xz(xp,z)
 10	continue

!$OMP END PARALLEL
!
c Note by PKY: in cylindrical truncation version the actual
c y-transform from wavenumber space to physical space is
c now absorbed into kxcomm1_trans_cyl (called from itransform.f)
c
! transform in y-direction to physical space
!
	deallocate(shift)
c
        return 
      end subroutine phase1
