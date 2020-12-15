       subroutine partsp_4
c
c  routine to perform particle tracking in conjunction
c    with direct numerical simulations of turbulence.
c
#ifdef LAG
#ifdef RKFOUR
c
	use compart
	implicit none
	include 'intvars'
c
c
	integer indexp,j
c
c form splines and interpolate for the velocity
c
#ifdef FROZEN_SP
        call rk_part (0)
#else
        call rk_part (1)
#endif
c
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,(pp(1,j),j=1,6)
  902	format (2i3,3x,'  pp(1,)=',1p,6e12.4)
 	end if
cc
c partial update of particle position

	call rk4_advpos (4, 0., 1./6.)
c
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,(pp(1,j),j=1,6)
	write (10000+taskid,902) istep,kstep,(pp(1,j),j=7,9)
 	end if
c
	iout=0
c
c
#endif
#endif
c
      return
      end
