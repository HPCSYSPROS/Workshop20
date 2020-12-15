        subroutine partsp (icall)
c
c  routine to perform particle tracking in conjunction
c    with direct numerical simulations of turbulence.
c
#ifdef LAG
c
	use compart
	use lag_timers
	implicit none
	include 'intvars'
c

c the "s(3)" array is renamed as "ps(3)" in this version
c
	real*8 rtime1,rtime2,rtime0
      real um(3)
      real moldif,mvrms,fact1,fact2,factor,hdx,hdt
	integer im1,im2,i,k
      save hdt,idv
	real(4) cof1,cof2
c
	real(p8), allocatable ::  mwork(:)
c
	integer indexp,jj,ii,icall,idv,j,luwp,igm
c


      iout=0
c
c if lpfunc(1).eq.0, no particle tracking is done
c
c
c define scale factors for particle coordinates.
c as far as particle coords are concerned, the length of the box
c is nx / sqrt (b11(2)), etc., with boundaries at gx(1) and
c gx(nx) + 1/sqrt(b11(m))
c
c gx(i) is location of i-th x-plane, etc.
c
      ps(1)=sqrt(b11(2))
      ps(2)=sqrt(b22(2))
      ps(3)=sqrt(b33(2))
c
      do 1 i=1,nx
 1    gx(i)=1+(i-1)/ps(1)+gsh(1,1)
      do 2 i=1,ny
 2    gy(i)=1+(i-1)/ps(2)+gsh(2,1)
      do 3 i=1,nz
 3    gz(i)=1+(i-1)/ps(3)+gsh(3,1)
c
      gxyz(1)=gx(1)
      gxyz(2)=gy(1)
      gxyz(3)=gz(1)
c
c define size of box to pass as argument to interpolation routines
c
	xyzl(1)=nx
	xyzl(2)=ny
	xyzl(3)=nz

c
c	if (taskid.eq.0) write (6,*) 'enter partsp, icall=',icall
c
	if (jstep.gt.1) call ppminmax ('enter partsp, vel',pp(1,4),ndpart)
	if (icall.eq.1) go to 100
	if (icall.eq.2) go to 200
	if (icall.eq.4) go to 400
c
c-------  initialization  ----------------------------------------------
c
 100  continue
c
c if spline derivatives are used to calculate velcocity gradients
c following the fluid particles, the condition
c dudx+dvdy+dwdz=0 is generally not satisfied. to remedy this,
c it is proposed to subtract 1/3 of the dui/dxi from each normal
c component of the velocity gradient tensor
c
c thus, if any of lpgrad(1,2,3) is 1, the other two must be, in order
c that this correction can be applied
c
c if lpmdis is set, velocity gradients must also be obtained to
c calculate the mechanical dissipation. however, if lpgrad(1,2,3)
c were all initialised as zeros, then velocity gradients will not
c be written on file, and are allowed to be overwritten once the
c mechanical dissipation is completely formed and saved
c
      if (lpgrad(1).eq.0.and.lpgrad(2).eq.0.and.lpgrad(3).eq.0) then
      npo(2)=0
      else
      npo(2)=8
      end if
c
c set flags needed for particle velocity and position output
c
      npo(1)=3*iovel
      npo(nplu-1)=3*ioposn
#ifdef LGRAD
      npo(nplu-2)=lplapl(1)+lplapl(2)+lplapl(3)
#endif
c
#ifdef NOSCALAR
#else
c
c set flags for scalars and scalar gradients
      do k=4,3+ncop
      npo(k)=lpfunc(k)+3*lpgrad(k)+lplapl(k)
      end do
#endif
c
c set flags for dissipation, enstrophy and pseudo-dissipation
c
      npo(3)=0
	do i=1,3
	if (lpsdis(i).gt.0) npo(3)=npo(3)+1
	end do
c
c     if (lpmdis.eq.1) lpgrad(1)=1
c
      if (lpgrad(1).eq.lpgrad(2).and.lpgrad(2).eq.lpgrad(3)) go to 111
c
      lpgrad(1)=1
      lpgrad(2)=1
      lpgrad(3)=1
c
      write (6,*) 'warning: lpgrad(1,2,3) are all set equal to 1'
      write (6,*) 'since the use of spline derivatives do not'
      write (6,*) 'satisfy the physical continuity condition'
      write (6,*) 'and a correction is applied based on dui/dxi'
c
 111  continue
c
c validity checks for Lagrangian parameters
c
c     call lvalid
c
c generate/read particle initial positions
c
      call partin

      nxyz(1)=nx
      nxyz(2)=ny
      nxyz(3)=nz
c
c     if (psave.ge.0) call fopen1 (luwp,'savepos','unformatted')
c
c perform preliminary operations for solution of spline equations
c
      do 110 i=1,3
	nd(i)=nxyz(i)
      call solsp1 (p(1,i),q(0,i),t(1,i),nd(i),ss,denom(i))
 110  continue
c
c define order-of-differentiation ordering for variables,
c used in this routine and sub. popsp
c
      call order
c
      go to 90
c
c-------  first call ( predictor step for particles )  ----------------
c
c  on entry,  pp(1-3)  contains the backward positions,
c             pp(4-6)  contains the backward velocities.
c
c  form predictor position x = x + u * dt  in  pp(1-3)
c  store  x + u * dt/2   in  pp(7-9)
c
c  on this call the dns array  uxzr(x,i,z,yp) contains velocity and
c    scalar fields in wavenumber space at the start of the eulerian
c    predictor step.  this corresponds to the forward time for the
c    particle tracking.
c
 200	continue

#ifdef LAGSC2
#ifdef SEP0411
	allocate (bfpos3(3,ndpart2/nrec2,nrec2))
#endif
#endif
c
	if (jstep.eq.1) go to 90
c
	if (ioflag.eq.1) iout=1
      if ((isflag.eq.1.or.istop.eq.1).and.isave.ge.0) then
	iout=2
	end if

       hdt = .5 * dt /2/pi
c
#ifdef NOV15
#else
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,icall,(pp(1,j),j=1,6)
	write (20000+taskid,902) istep,kstep,icall,(pp2(1,j),j=1,6)
	end if
#endif
c
	rtime1=MPI_WTIME()

	cof1=.5
	cof2=.5
#ifdef NOV15
#else
 	call rk_advpos (1,cof1,cof2)
	if (jstep.gt.1) call ppminmax (' partsp, A, vel',pp(1,4),ndpart)
#endif
#ifdef MOL
 	if (nom.gt.0) call rk_advmpos (1)
#endif
c
	jj=1
	if (iout.gt.0) jj=2
c
	if(icall.eq.2) then
	cpu_partsp(2,jj)=cpu_partsp(2,jj) + MPI_WTIME() - rtime1
c	if (taskid.eq.0) write (6,*) 'istep,jj,cpu_partsp=',istep,jj,cpu_partsp(2,jj)
	endif
c

	rtime1=MPI_WTIME()

#ifdef NOV15
#else
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,icall,(pp(1,j),j=1,6)
	write (10000+taskid,902) istep,kstep,icall,(pp(1,j),j=7,9)
	write (20000+taskid,902) istep,kstep,icall,(pp2(1,j),j=1,6)
	write (20000+taskid,902) istep,kstep,icall,(pp2(1,j),j=7,9)
	end if
#endif
c
	if (jstep.gt.1) call ppminmaxv (' partsp, B, vel',pp(1,4),ndpart)
	call rk_part (1)
	if (jstep.gt.1) call ppminmaxv (' partsp, C, vel',pp(1,4),ndpart)
c
	if(icall.eq.2) then
	cpu_partsp(3,jj)=cpu_partsp(3,jj) + MPI_WTIME() - rtime1
	endif

c  determine interpolation quantities based on predictor position,
c  with an appropriate geometric factor to account for the
c  slanted shape of the box (for shear flows)
c
!c	go to 90
c
c-------  second call  ( corrector step for particles )  ---------------
c
c  on this call the dns array  u(x,i,y,z)  contains velocity and scalar
c   fields on the basic grid in physical space at the forward time.
c
c  1.  form the the corrector velocity.
c  2.  perform the corrector step to obtain the new particle positions,
c        stored in pp(1-3).
c  3.  form new velocities pp(4-6).
c  4.  if output is required, interpolate for further particle propertie
c
 300	continue
c
	rtime1 = MPI_WTIME()
c
	cof1=0.
	cof2=0.5
#ifdef NOV15
#else
	call rk_advpos (2, cof1,cof2)
	if (jstep.gt.1) call ppminmax (' partsp, D, vel',pp(1,4),ndpart)
#endif
#ifdef MOL
	call rk_advmpos (2)
#endif

#ifdef NOV15
#else
	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,icall,(pp(1,j),j=1,6)
	write (20000+taskid,902) istep,kstep,icall,(pp2(1,j),j=1,6)
	end if
#endif
c
	if(icall.eq.2) then
	cpu_partsp(4,jj)=cpu_partsp(4,jj) + MPI_WTIME() - rtime1
	endif



	rtime1=MPI_WTIME()
c
c save particle positions on file if desired
c code stops here if particles are tracked
c
#ifdef TEST
 325  if ( (psave.eq.0 .and. jstep-1.eq.nsteps).or.
     1     (psave.gt.0 .and. mod(jstep-1,psave).eq.0).
     1     or.(entime.gt.0..and.time.ge.entime).
     1     or.(stend.gt.0..and.shear*time.ge.entime-1.e-6).
     1     or.(isflag.eq.1).or.(istop.eq.1) ) then
#endif
c
c     if (isflag.eq.1.or.istop.eq.1) then
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
c
	icpvrms=icall
 	call rk_part (0)
	icpvrms=0
  912   format (i3,1p,4e14.6)

	if(icall.eq.2) then
	cpu_partsp(5,jj)=cpu_partsp(5,jj) + MPI_WTIME() - rtime1
	endif
c
	if (ioflag.eq.1) then
c
c check rms particle velocity
c
 	call check_pvrms
c
c write particle output at current time step
c
	iout=1
c
c interpolate for scalars, scalar gradients and velocity gradients  if needed
c
	indexp=6
 	call part_interp (indexp)

 	end if

	go to 90
c
c
 400	continue
c
#ifdef LAGSC2
#ifdef SEP0411
	allocate (bfpos3(3,ndpart2/nrec2,nrec2))
#endif
#endif
c write particle output at current time step
c
	iout=1
c
	call rk_part (1)
c
 	if (mod(taskid,jproc).eq.0) then
	write (10000+taskid,902) istep,kstep,icall,(pp(1,j),j=1,6)
	write (20000+taskid,902) istep,kstep,icall,(pp2(1,j),j=1,6)
  902	format (3i3,'  pp(1,)=',1p,6e12.4)
 	end if
cc
#ifdef MOL
c to ensure correct  checkpointing for molecules, must draw
c some dummy samples for the Brownian-motion random no. generator
c when jstep=1
c

	if (nom.gt.0) then
        if(taskid.eq.0) write(6,*) 'partsp: drawing dummy samples'
	allocate (mwork(nm))
      if (pstart.gt.0.and.jstep.eq.1) then
#ifdef MOLOLD
      do igm=1,ngm
      call ranseq (12+igm,0)
      call rann2 (mwork,1,nm,1,1,nm)
      call rann2 (mwork,1,nm,1,1,nm)
      call rann2 (mwork,1,nm,1,1,nm)
      end do
#else
      call random_seed(put=mseed)
      call random_number(mwork)
      call random_number(mwork)
      call random_number(mwork)
      call random_seed(get=mseed)
#endif
      end if
	deallocate (mwork)
      if (taskid.eq.0) write (6,*) 'partsp: after initial molecule loop'
	end if

c
#endif

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
 90	continue
	
c        if (taskid.eq.0.and.icall.eq.2) then
c        write (6,*) 'partsp: istep,iout=',iout
c        write (6,"('partsp:',i3,1p,6e12.4)") istep,(cpu_partsp(ii,1),ii=0,5)
c        write (6,"('partsp:',i3,1p,6e12.4)") istep,(cpu_partsp(ii,2),ii=0,5)
c        end if

c	if (taskid.eq.0) write (6,*) ' exit partsp, icall=',icall
c
#endif

#ifdef LAGSC2
#ifdef SEP0411
	if(icall.ne.1) deallocate (bfpos3)
#endif
#endif
c
      return
      end
