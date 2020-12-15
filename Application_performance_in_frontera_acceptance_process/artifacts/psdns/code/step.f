      subroutine step
 
! Revised June 2002 to save memory by avoiding xyfac and xydif
! arrays except for shear flow
! ---- some round-off errors versus old results can be expected
! ---- memory saving can be significant for large grids
! 
!
! this version contains corrections of integrating factors
! for pure shear case, by p.k. at penn state, 1/22/90

! generalized to cases where num_thr may not be factor of nxh,nyh,nzh
! (D. Buaria, 8/23/12)
 
	use comp
c	implicit none
#include "intvars"
 
      real, allocatable :: logx(:),logy(:),logz(:),bk12(:)
      real, allocatable :: vmaxz(:)
!
	real tmpvmaxz
 
      real logxy
c
        integer ithr,omp_get_thread_num,ix1,ix2,iy1,iy2,iz1,iz2
	integer, allocatable :: bxp1(:),bxpi(:)
	integer, allocatable :: byp1(:),bypi(:)
	integer, allocatable :: bzp1(:),bzpi(:)
	save bxp1,byp1,bzp1,bxpi,bypi,bzpi

 
      save dt0
      data term1,term2,logxy/3*0./
c     data strof/1.e-5/
      data strof/1.e-6/
c
      allocate (logx(nxhpad))
      allocate (logy(nypad))
      allocate (logz(nzpad))
      allocate (bk12(nxhpad))
      allocate (vmaxz(numtasks))
 
      if (dtout.gt.0) iostep=99999999

	icase=1
 
	
      if (jstep.eq.1.and.cfl.le.0.) dt0=dt
      dt=dt0
!
!  determine time step and advance time  -------------------------------
!
!  calculate dt from cfl condition if appropriate (subject to possible
!  adjustments in shear flows)
!
!
! these lines were found to be slow and causing variability
! in timings (T. Minyard, TACC, 8/12/08)
!     call MPI_ALLGATHER (velmax,1,MPI_REAL,vmaxz,1,MPI_REAL,
!    1                    MPI_COMM_WORLD,mpierr)
!     velmax=maxval(vmaxz)
!
      call MPI_ALLREDUCE (velmax,tmpvmaxz,1,mpireal,MPI_MAX,
     1                    MPI_COMM_WORLD,mpierr)
      velmax=tmpvmaxz

#ifdef LINDBORG
	if (istep.eq.1) then
	dt=dt0
	else
	if (cfl.gt.0.) then
	dt=cfl/velmax
	end if
	end if
#else
      if (cfl.gt.0.) then
      dt=cfl/velmax
      if (velmax.lt..001) then
      dt=.9*dt0
      write (6,*) 'velmax gotten zero and dt=.9*dt0 at istep=',istep,
     1             velmax,dt
      call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
      end if
#endif
!
! for shear-flow, do not allow non-dimensional time to exceed st=stend
! (if stend has been given an non-zero value)
! otherwise, stop at time=entime if entime.gt.0.
!
      if (shear.gt.0..and.stend.gt.0.and.
     1    (time+dt)*shear.ge.stend-strof) then
      dt=stend/shear-time
      write (6,*) 'step, dt=',dt
      if (abs(dt).lt.strof) stop ' exit: stend reached in sub. step'
      else if (icase.eq.1.and.entime.gt.0..and.(time+dt).ge.entime)then
      dt=entime-time
      end if
!
!      advance mesh metric  --------------------------------------------
!      (constant irrotational strain or pure shear)
!
!  b11 = b(1,1)**2, etc. , b12 = b(1,2) / b(2,2)
!
        b11(1)=b11(2)
        b22(1)=b22(2)
        b33(1)=b33(2)
        b12(1)=b12(2)
!
!       b12(2)=b12(1)-s0*dt
!
! in shear flow runs, re-adjust time step size to obtain an orthogonal
! grid if b12 changes sign
!
        b12(2)=b12(1)-s0*b11(1)*dt
      if (b12(1).gt.0.01.and.b12(2).lt.-0.01) then
      dt=b12(1)/s0/b11(1)
      if (taskid.eq.0) write (6,*) ' dt at sign change of b12=',dt
      if (taskid.eq.0) write (6,*) ' b12,s0,b11=',b12(1),s0,b11(1)
      end if
!
      iopflag=ioflag
!
      ioflag=0
!
! also adjust for integral values of st
!
      if (shear.gt.0.) then
      st1=shear*time
      st2=shear*(time+dt)
      ita=st1
      itb=ita+1
      ist1=st1
      if (abs(st1-float(ita)).le.strof) ist1=ita
      if (abs(st1-float(itb)).le.strof) ist1=itb
      ita=st2
      itb=ita+1
      ist2=st2
      if (abs(st2-float(ita)).le.strof) ist2=ita
      if (abs(st2-float(itb)).le.strof) ist2=itb
      if (ist2.gt.ist1) then
      dt=(ist2/shear-time)
      if (taskid.eq.0) write (6,*) ' dt at integral st=',dt
      ioflag=1
      end if
      end if
!
! set tfout and tfsave to time0 if they were initialized with
! negative values
!
      if (jstep.eq.1.and.tfout.lt.0.) tfout=time0
      if (jstep.eq.1.and.tfsave.lt.0.) tfsave=time0
!
! adjust for output at times equal to multiples of dtout since tfout
!
      if (dtout.gt.0.and.shear.eq.0.) then
      t1=time
      t2=time+dt
      it1=(t1-tfout)/dtout
      it2=(t2-tfout)/dtout
      ita=it1
      itb=ita+1
      if (abs(t1-itb*dtout+tfout).le.strof) it1=itb
      ita=it2
      itb=ita+1
      if (abs(t2-itb*dtout+tfout).le.strof) it2=itb
      ioflag=-1
      if (it2.gt.it1) then
      ist2=it1+1
      dt=(it2*dtout+tfout)-time
      ioflag=1
      if (taskid.eq.0) write (6,*) ' dt adjusted in dtout loop, dt=',dt
      end if
      end if
!
! for shear flows, interpret dtout as shear*dtout
! now, output is written only at the specified intervals of St
!
      if (dtout.gt.0.and.shear.gt.0.) then
      st1=shear*time
      st2=shear*(time+dt)
      ist1=(st1-tfout)/dtout
      ist2=(st2-tfout)/dtout
      ita=ist1
      itb=ita+1
      if (abs(st1-itb*dtout+tfout).le.strof) ist1=itb
      ita=ist2
      itb=ita+1
      if (abs(st2-itb*dtout+tfout).le.strof) ist2=itb
      ioflag=-1
      if (ist2.gt.ist1) then
      ist2=ist1+1
      dt=(ist2*dtout+tfout)/shear-time
      ioflag=1
      end if
      end if
!
!
! set the "save" flag at times equal to multiples of dtsave since tfsave
!
      isflag=0
!
      if (dtsave.gt.0..and.shear.eq.0.) then
      t1=time
      t2=time+dt
      if (abs(mod(t2-tfsave,dtsave)).le.strof) isflag=1
      it1=(t1-tfsave)/dtsave
      it2=(t2-tfsave)/dtsave
      ita=it1
      itb=ita+1
      if (abs(t1-itb*dtsave+tfsave).le.strof) it1=itb
      ita=it2
      itb=ita+1
      if (abs(t2-itb*dtsave+tfsave).le.strof) it2=itb
      if (it2.gt.it1) then
      ist2=it1+1
      dt=(it2*dtsave+tfsave)-time
      ioflag=1
      isflag=1
      if (taskid.eq.0)write (6,*)' dt adjusted in dtsave loop, dt=',dt
      end if
      end if
!
! for shear flows, interpret dtsave as shear*dtsave
!
      if (dtsave.gt.0..and.shear.gt.0.) then
      st1=shear*time
      st2=shear*(time+dt)
      ist1=(st1-tfsave)/dtsave
      ist2=(st2-tfsave)/dtsave
      ita=ist1
      itb=ita+1
      if (abs(st1-itb*dtsave+tfsave).le.strof) ist1=itb
      ita=ist2
      itb=ita+1
      if (abs(st2-itb*dtsave+tfsave).le.strof) ist2=itb
      if (ist2.gt.ist1) then
      dt=(ist2*dtsave+tfsave)/shear-time
      if (taskid.eq.0) write (6,*) ' dt adjusted in dtsave loop, dt=',dt
      isflag=1
      end if
      end if
!
 111	continue
!
        b11(2)=b11(1)*exp(-2.*a11*dt)
        b22(2)=b22(1)*exp(-2.*a22*dt)
        b33(2)=b33(1)*exp(-2.*a33*dt)
        b12(2)=b12(1)-s0*b11(1)*dt
!
! increment the time
!
      t1=time
      time=time+dt
      if (dt.eq.0.) then
      write (6,*) ' abort: dt unexpectedly zero in sub. step, istep=',
     1              istep
      stop
      end if
      t2=time
      dt2=dt/2.
!
! print time step or courant number
!
      if( cfl .gt. 0. ) then
!
!  dt from cfl
!     
         dtmin=min(dt,dtmin)
         if (shear.eq.0.) then
            if (taskid.eq.0) write (6,601) istep,dt,time,cfl,velmax
 601        format (' istep,dt,time,cfl,velmax=',i6,1p,4e12.5)
         else
            if (taskid.eq.0) write (6,6010) istep,dt,shear*time,cfl,velmax
 6010       format (' istep,dt,st,cfl,velmax=',i6,1p,4e12.5)
         end if
!     
!     if (dt.gt.10..or.dt.lt.1.e-6) then
         if (dt.gt.10..or.dt.lt.1.e-7) then
            write (6,*) 'dt gotten out of control in sub. step, dt=',dt
            call MPI_ABORT (MPI_COMM_WORLD,ierror)
         end if
!     
      else if( cfl .lt. 0. ) then
!
!  constant dt, monitor cfl
!
         cn=velmax*dt/2./pi
         if (taskid.eq.0) write (6,602) istep,dt,time,cn,velmax
 602     format (' istep,dt,time,courant no=',i6,1p,4e12.5)
         if (cn.gt.1) then 
	call MPI_FINALIZE (mpierr)
	stop 'courant no. limit exceeded'
	end if
         cnmax=max(cn,cnmax)
!     
      endif
!     
! pure shear case calculation of shear-contributions to integrating
! factor
!
      if (icase.eq.2) then
         term1=bet12**2*dt-s0*(t2*t2-t1*t1)*bet12
     1        +s0*s0*(t2**3-t1**3)/3.
         term2=(b12(1)+b12(2))*dt
      end if
!     
!     define physical wave nos. ---------------------------------------
!     
      call waveno
!     
!     define integrating factor  --------------------------------------
!     
!     xyfac*zfac is the inverse of the integrating factor
!     ( similarly xydif*zdif )
!     
#ifdef LINDBORG
      pr_h=visc_h/difs_h
      pr_v=visc_v/difs_v
      term=b11(1)*dt
	do 10 y=1,nypad
        do 10 x=1,nxhpad
	logxy=visc_h*(kx(x)**2+ky(y)**2)**4 * term * b11(1)**3
	xyfac(x,y)=exp(-logxy)
	do i=1,nc
	xydif(x,y,i)=exp(-logxy/pr_h)
	end do
 10	continue
        term=b33(1)*dt
        do 30 z=1,nzpad
        logz(z)=visc_v*kz(z)**8 * term * b11(1)**3
	zfac(z)=exp(-logz(z))
	do i=1,nc
	zdif(z,i)=exp(-logz(z)/pr_v)
	end do
 30	continue
c
c
#else
	if (.not.allocated(bxp1)) then
	allocate (bxp1(0:num_thr-1),bxpi(0:num_thr-1))
	allocate (byp1(0:num_thr-1),bypi(0:num_thr-1))
	allocate (bzp1(0:num_thr-1),bzpi(0:num_thr-1))
	end if

	call divide_thr (nxhpad, bxp1, bxpi, 'step: nxhpad')
	call divide_thr (nypad, byp1, bypi, 'step: nypad')
	call divide_thr (nzpad, bzp1, bzpi, 'step: nzpad')

	ithr=0
c
!$OMP PARALLEL private (ithr,ix1,ix2,iy1,iy2,iz1,iz2,term)
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
c
	ix1=bxp1(ithr)
	ix2=ix1+bxpi(ithr)-1
	iy1=byp1(ithr)
	iy2=iy1+bypi(ithr)-1
	iz1=bzp1(ithr)
	iz2=iz1+bzpi(ithr)-1
c
        term=b11(1)*dt
        if( a11 .ne. 0. ) term=( b11(1) - b11(2) ) / ( 2. * a11 )
        do 10 x=ix1,ix2
        logx(x)=-viscos*kx(x)**2 * term
	  xfac(x)=exp(logx(x))
 10	continue
!
        term=b22(1)*dt
        if( a22 .ne. 0. ) term=( b22(1) - b22(2) ) / ( 2. * a22 )
        do 20 y=iy1,iy2
        logy(y)=-viscos*ky(y)**2 * term
        yfac(y)=exp(logy(y))
 20	continue
!
        term=b33(1)*dt
        if( a33 .ne. 0. ) term=( b33(1) - b33(2) ) / ( 2. * a33 )
        do 30 z=iz1,iz2
        logz(z)=-viscos*kz(z)**2 * term
        zfac(z)=exp(logz(z))
 30	continue
c
#ifndef NOSCALAR
      do i=1,nc
      do z=iz1,iz2
      zdif(z,i)=zfac(z)**(1./pr(i))
      end do
      end do
#endif
!
      if (icase.eq.1) then
!
#ifndef NOSCALAR
      do i=1,nc
      do x=ix1,ix2
      xdif(x,i)=xfac(x)**(1./pr(i))
      end do
 
      do y=iy1,iy2
      ydif(y,i)=yfac(y)**(1./pr(i))
      end do
      end do
#endif

      else if (icase.eq.2) then
!
      vterm=-viscos*beta2**2
      termb1=s0*b11(1)*dt
!
        do 50 y=1,nypad
        do 50 x=1,nxhpad
!
!       logxy=vterm*(term1*kx(x)**2+term2*kx(x)*ky(y))
!       term=logx(x)+logy(y)+logxy
! alternative version (should be equivalent)
#ifdef SHEAR
        term=-viscos*dt*(bkk12(x,y,1)+termb1*b22(1)*kx(x)
     1        *(termb1*kx(x)/3.-k2(x,y,1)))
        xyfac(x,y)=exp(term)
#ifndef NOSCALAR
	do i=1,nc
	xydif(x,y,i)=exp(term/pr(i))
	end do
#endif
#endif
!
 50	continue
!
      end if
!
!$OMP END PARALLEL
c
#endif
	if (ioflag.eq.1) iocount=iocount+1
!
	if (dtout.eq.0..and.mod(istep-istep0,iostep).eq.0) then
	iocount=iocount+1
	end if
c
#ifdef LVSTEP
	if (ivstep.gt.1) then
	xfac(:)=xfac(:)**ivstep
	yfac(:)=yfac(:)**ivstep
	zfac(:)=zfac(:)**ivstep
	end if
#endif
!
      deallocate (logx,logy,logz,bk12,vmaxz)

        return
        end
