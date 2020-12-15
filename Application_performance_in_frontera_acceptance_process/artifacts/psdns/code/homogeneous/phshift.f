      subroutine phshift (uy,vy,nv,nps,ev)
c
	use comp
	implicit none
	include 'intvars'
c
        integer nv,ev
      complex(b8) uy(ny,zjsz*xisz,nv)
      complex(b8) vy(ny,zjsz*xisz,nv)

c apply phase-shifting only
c
c	complex sxz,shifty
	complex(b8) sxz,shifty
c
	integer i
        integer a
	integer nps
c
        integer ithr,omp_get_thread_num,ia
c
c
c form shift, and shift the velocities,

!$OMP PARALLEL private (ithr,x,xp,z,sxz,shifty)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 5 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
        x=xp+xist-1
c
      sxz=sx(x,ev)*sz(z,ev)
c
      do 10 y=1,ny

        if (mask(y,a)) then
c
      shifty=sy(y,ev)*sxz
c
	do i=1,nv
	vy(y,a,i)=uy(y,a,i)*shifty
	end do
c
	else
c
	do i=1,nps
	vy(y,a,i)=0.
	end do
c
	end if
c
 10   continue
c
        call next_xz(xp,z)
c
 5      continue
c
!$OMP END PARALLEL
c
      return
      end
!-----------------------------------------------------------------
      subroutine phshift_inplace (uy,nv,ev)
c
	use comp
	implicit none
	include 'intvars'
c
        integer nv,ev
      complex(b8) uy(ny,zjsz*xisz,nv)

c apply phase-shifting only
c
	complex(b8) sxz,shifty
c
	integer i
        integer a
c
        integer ithr,omp_get_thread_num,ia
c
c form shift, and shift the velocities,

!$OMP PARALLEL private (ithr,x,xp,z,sxz,shifty)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 5 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
        x=xp+xist-1
c
      sxz=sx(x,ev)*sz(z,ev)
c
      do 10 y=1,ny

        if (mask(y,a)) then
c
      shifty=sy(y,ev)*sxz

	do i=1,nv
	uy(y,a,i)=uy(y,a,i)*shifty
	end do
c
	else
c
	do i=1,nv
	uy(y,a,i)=0.
	end do
c
	end if
c
 10   continue
c
        call next_xz(xp,z)
c
 5      continue
c
!$OMP END PARALLEL
c
      return
      end
