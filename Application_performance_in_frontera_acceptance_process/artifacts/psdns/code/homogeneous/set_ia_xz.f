      subroutine set_ia_xz
c
c routine to set ranges of the "a" index used in managing
c cylindtrical truncation in hybrid MPI-OpenMP code
c
c (this info is used, e.g., by masks.f)

      use comp
      implicit none
      include 'intvars'
	integer a
c
        integer ithr,ia
c

      xp = mystart_x
      z = mystart_z
c
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
      
      return
      end subroutine set_ia_xz

