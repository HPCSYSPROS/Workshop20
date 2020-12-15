      subroutine masks
c
c cylindrical truncation version

      use comp
      implicit none
      include 'intvars'
	integer a
      real(b8) :: sqkx,sqkz,sqk,two9,rnxi,rnyi,rnzi
c
        integer ithr,omp_get_thread_num,ia
c
c 12/31/11: note: the following lines are now placed
c in a standalone subroutine, "set_ia_xz" that must
c be called at least once before 'masks'

c
#ifdef TEST
      xp = mystart_x
      z = mystart_z
	allocate (ia_sz(0:num_thr-1),ia_xst(0:num_thr-1),ia_zst(0:num_thr-1))
c
	do ithr=0,num_thr-2
	ia_sz(ithr)=ia_st(ithr+1)-ia_st(ithr)
	end do
	ia_sz(num_thr-1)=num_al+1-ia_st(num_thr-1)
	
	do ithr=0,num_thr-1
	do ia=1,ia_sz(ithr)
	a=ia_st(ithr)+ia-1
	if (ia.eq.1) then
	ia_xst(ithr)=xp
	ia_zst(ithr)=z
	end if
        call next_xz(xp,z)
 	end do
 	end do
#endif
c
c	call set_ia_xz
      
	
! calculate the contents of the 'mask' array
! for 1 z-plane only
c
      rnxi = 1./nx
      rnyi = 1./ny
      rnzi = 1./nz
      two9 = 2./9.

!$OMP PARALLEL private (ithr,sqkx,sqkz,sqk,x,xp,z)

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

         sqkx = ( kx(x)*rnxi )**2 
         sqkz = ( kz(z)*rnzi )**2 + sqkx
c

         do 30 y=1,ny
            sqk = ( ky(y)*rnyi )**2 + sqkz               
            if( sqk .gt. two9 .or. y .eq. nyhp .or. z .eq. nzhp) then
               mask(y,a) = .false.
            else
               mask(y,a) = .true.
            endif
            
 30      continue
c
        call next_xz(xp,z)

 10	continue

!$OMP END PARALLEL
c

      return
      end subroutine masks

