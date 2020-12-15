       subroutine partsp_1
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
	integer indexp,i,k,j,luwp
c
	do k=1,3
	do i=1,nop/numtasks
	pp(i,ndprop-3+k)=pp(i,k)
	pp(i,6+k)=pp(i,k)
	end do
	end do
#ifdef LAGSC2
	do k=1,3
	do i=1,nop2/numtasks
	pp2(i,ndprop2-3+k)=pp2(i,k)
	pp2(i,6+k)=pp2(i,k)
	end do
	end do
#endif
c
c form splines and interpolate for the velocity
c
#ifdef FROZEN_SP
	if (jstep.eq.1) then
	 call rk_part (1)
	else
	 call rk_part (0)
	end if
#else
	call rk_part (1)
#endif
c
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,(pp(1,j),j=1,6)
  902	format (2i3,3x,'  pp(1,)=',1p,6e12.4)
	write (20000+taskid,912) istep,time,(pp(1,j),j=4,6)
  912	format (i3,1p,4e14.6)
	write (50000+taskid,905) istep,pp(1,1),pp(1,4)
  905	format (i3,'  ppos,pvel=',1p,2e14.6)
 	end if
cc
	iout=0
c
      if ((isflag.eq.1.or.istop.eq.1).and.isave.ge.0) then
      luwp=120
      call wrtpos (luwp)
      end if
c
      if (jstep-1.eq.nsteps.or.(entime.gt.0.and.time.ge.entime-1.e-6).
     1        or.(stend.gt.0.and.shear*time.ge.stend-1.e-6)) then
      iout=2
	go to 90
	end if

	if (ioflag.eq.1) then
c
c check rms particle velocity
c
	call check_pvrms
c
c interpolate for scalars, scalar gradients and velocity gradients  if needed
c
 	indexp=6
  	call part_interp (indexp)
c
c write particle output at current time step
c
	iout=1
c
	end if
c
c
c partial update of particle position

	call rk4_advpos (1, .5, 1./6.)
c
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,(pp(1,j),j=1,6)
	write (10000+taskid,902) istep,kstep,(pp(1,j),j=7,9)
 	end if
c
 90	continue
c
#endif
#endif
c
      return
      end
