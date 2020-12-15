! PSDNS_homo_omp: 
! Code for Petascale Direct Numerical Simulation of Homogeneous Turbulence,
! written in MPI and OpenMP
! Maintained by research group of P.K Yeung at Georgia Tech, USA
! With major contributions below by several individuals (and others)
!
! Fourier pseudo-spectral in space and 2nd or 4th order Runge Kutta in time
! (Rogallo 1981), cubic-spline interpolation for particle tracking
! (Yeung \& Pope 1988)

! --- 2D domain decomposition (D.A. Donzis, D. Pekurovsky)
! --- 1-D FFT using FFTW library (D.A. Donzis) or ESSL on IBMs
! --- cylindrical truncation in wavenumber space (D. Pekurovsky)
! --- re-factoring for stride-1 arithmetic (D. Pekurovsky, N. Malaya)
! --- alltoall communication by Co-Array Fortran (R.A. Fiedler)
! --- changes for particle tracking at 8192^3 on Blue Waters (D. Buaria)
! --- tracking of diffusing molecules (D. Buaria)
! --- axisymmetric spectra for anisotropic flows (M.P. Clay)
! --- parallel HDF5 capability for checkpoint and restart (M.P. Clay)
! --- 1 file per MPI task or 1 file per communicator
! --- allows for uniform solid-body rotation and low-Rm MHD turbulence
! --- options for numerical forcing
! --- passive scalars, including smaller time steps at low Schmidt number

! --- Runs on several HPC platforms: Cray XE, Intel, IBM BGQ
! --- Single and double precision via compiler directives in makefile

      program DNS        
        use rungekutta_module !defines and initializes RK scheme        
        use timestep_module   !contains time discretization scheme
        use wavespace_module  !includes homogeneous, inhomogeneous

       	use comp
	use timers_comm
	use timers_comp
	use timers_io
	use timers_rkstep
#ifdef RANDOM_SEED
	use ranseed
#endif
c
c#ifdef LAG 
#if defined (LAG) || defined (LAGSC2) || defined (MOL)
	use compart, only: nop,nop2,nom,lpgrad,ngm,nompp,nom,iolflag
#endif

        implicit none

      integer dmy(3),hms(3),date_time(8)
!
      integer i,n,j,m,iaux
      real twopi
	integer idim(4),idimc(2)
	real sum1,raux
      real(8) rtime1,rtime2,difftime
 	real totcpu,tcpu,cpumin,cpumax,avcpu,acccpu2
	real(8) tt_comm,gtt_comm(2)
	character*60 tsp
#ifdef IND_STDOUT	
      character*40 fn
#endif	
	logical :: iex 
	integer flag_outcomm
!
	integer itask
	real(8) sum

        character*20 caux
        character*11 form,unform
      data form,unform/'formatted','unformatted'/

        integer num,l,ithr,ia,jthr


c next line added Feb 21, 2013
c 8/10/2016: declaration of clock_chkpt is moved to the 'comp'
c section of module.f
	real(8) clock_start,clock_limit,clock_t,wall_limit,clock0
	integer ichktime, iosteps, extrasteps, stepst

#ifdef OPENMP
      integer prov, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
c
#ifdef OPENMP
      prov = 0
      call MPI_INIT_THREAD (MPI_THREAD_FUNNELED,prov,ierr)
#else
      call MPI_INIT (ierr)
#endif
      call MPI_COMM_SIZE (MPI_COMM_WORLD,numtasks,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,taskid,ierr)

#ifdef CAF
	if (taskid.eq.0) write (6,*) 'This code is compiled with CAF'
#endif
#ifdef ESSL
	if (taskid.eq.0) write (6,*) 'This code uses ESSL'
#endif

#ifdef IND_STDOUT
	write (fn,fmt="('stdout.',i0)") taskid
	open (6,file=fn)
#endif
!
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
	if (taskid.eq.0.or.taskid.eq.numtasks-1) then
       write (6,601) taskid,dmy(2),dmy(1),dmy(3),hms
  601  format (' starting taskid=',i6, ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
	end if
 
	if (taskid.eq.0) then
	inquire (file='clock_chkpt',exist=iex)
	clock_chkpt=0.
	if (iex) then
	open (1113,file='clock_chkpt')
	read (1113,*) clock_chkpt,wall_limit
	end if
	clock_start=MPI_WTIME()
	clock0=clock_start
	end if
	call MPI_BCAST (clock_chkpt,1,MPI_DOUBLE_PRECISION,
     1                  0,MPI_COMM_WORLD,mpierr)
c
c Dec 22, 2015
c lines added by Dhawal to help aim for '75% wallclock accuracy or better'
c (motivated by allocations discount policy on BW)

	ichktime = 0
	if (taskid.eq.0) then
	inquire (file='clock_entime',exist=iex)
	if (iex) then
	open (1113, file='clock_entime')
	read (1113,*) ichktime, wall_limit
	if (wall_limit.lt.100) wall_limit = wall_limit*3600.
	endif
	clock_start = MPI_WTIME()
	clock0 = clock_start
	endif


        if (taskid.eq.0) then
        inquire(file='write_comm',exist=iex)
	flag_outcomm=0
	if (iex) flag_outcomm=1
	write (6,*) 'flag_outcomm=',flag_outcomm
        end if
        call MPI_BCAST (flag_outcomm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	if (taskid.eq.1) write (6,*) 'flag_outcomm=',flag_outcomm

	ichkpt=0
c
      one = 1
      zero = 0
      twopi=2.*pi
 
      call input
 
 
      !#Initializing
      call param_set
      Ntot = nxhp*ny*nz/numtasks
      call rungekutta_initialize(rkmethod)     !initialize RK steps
      call initialize_wavespace !initialize inhomogeneous / homogeneous

	  if (taskid.eq.0) then
	write (6,"('numtasks=',i7,' nx=',i5,' ny=',i5,' nz=', i5,' nc=',i1,' nu=',i2)") numtasks,nx,ny,nz,nc,nu
	write (6,*) 'main: Ntot,nxhp=',Ntot,nxhp
       endif
 
	rtime1=MPI_WTIME()  ! time initialization tasks
	rwall0=rtime1


      call mpisetup
 
      call com_set
	if (taskid.eq.0) write (6,*) 'after com_set'
      call comp_set
	if (taskid.eq.0) write (6,*) 'after comp_set'
      call comsp_set
	if (taskid.eq.0) write (6,*) 'after comsp_set'

ccccccccccccccccccc
c	call compart_set 
c	call pop_file4 (1)
cccccccccccccccccccc

! set some parameters for checkpoint
! not overwrite all previous chkpts	
	iwfld=2
	do i=1,3+nc
 	fnwu(i)='outpenc'
 	luwu(i)=21
	luinit(i)=3
	enddo
c
 	
      if(taskid .eq. 0) then
         print *,'Using processor grid ',iproc,' x ',jproc
      endif
c
        if (mod(taskid,jproc).eq.0) then
        write (6,610) taskid,xist,xien,zjst,zjen
 610    format ('taskid, xist,xien,zjst,zjen=',5i6)
        write (6,620) taskid,zist,zien,yjst,yjen
 620    format ('taskid, zist,zien,yjst,yjen=',5i6)
        end if

 
! initialize some counters 
#ifdef TIMERS_IO	  
	twrite_io = 0.
	iwrite_io = 0
#endif

!ccccccccc for test purposes
	istep0=0
	istop=0
!ccccccccccccccccccccccccccc	

        

! initialize wavenumbers
        if(taskid .eq. 0) then
           print *,'Calling init; nxhpad=',nxhpad
        endif
	call init
        if(taskid .eq. 0) then
           write (6,*) 'Calling waveno ,b11(2)=',b11(2)
        endif
	call waveno

        if(taskid .eq. 0) then
           print *,'Calling mpisetup2'
        endif

      call mpisetup2
c
c lines moved here from epfftw.f since info is needed regardless
c of use of FFTW or ESSL
c
        allocate (ixp1(0:num_thr-1),ixpi(0:num_thr-1))
        allocate (iyp1(0:num_thr-1),iypi(0:num_thr-1))
        allocate (izp1(0:num_thr-1),izpi(0:num_thr-1))
        allocate (num_fft0(0:num_thr-1))
        num=num_al/num_thr
        num_fft0(:)=num
        l=mod(num_al,num_thr)
        if (l.ne.0) then
        do ithr=0,l-1
        num_fft0(ithr)=num+1
        end do
        end if
c       if (ipid.eq.0) then
        if (jpid.eq.0) then
        write (6,"('main: taskid,num_al,num_fft0=',8i6)")
     1       taskid,num_al,num_fft0
        end if


#ifdef FFTW
	call epfftw (u,u,u,u)
#endif
#ifdef ESSL
         call init_work
#endif
c
        allocate (ia_st(0:num_thr-1))
        ia=1
        do jthr=0,num_thr-1
        ia_st(jthr)=ia
        ia=ia+num_fft0(jthr)
        end do
c	write (6,*) 'taskid,ia_st=',taskid,ia_st
c
#ifdef HOMOGENEOUS

        if(taskid .eq. 0) then
           print *,'Calling masks'
        endif
	call set_ia_xz
      call masks
c

        if(taskid .eq. 0) then
           print *,'After masks'
        endif
#endif


	call openf

 
! for testing: initialize a sinusoidal and solenoidal
! velocity field in 3D by calling sub. pfield
! (start with x-z slabs in physical space)
 
#ifdef RANDOM_SEED
        call random_seed (size=seed_size)
        allocate (rseed(seed_size))
#endif

        


        if(taskid .eq. 0) then
           write (6,*) 'main: before call incond'
        endif
 	call incond
        if(taskid .eq. 0) then
           write (6,*) 'main: after call incond'
        endif
c
	if (kinit(1).eq.-2) call sptvar (un)
c

	acccpu=0.
      rtime2 = MPI_WTIME()   ! time initialization task
	difftime= rtime2 - rtime1
      acccpu=acccpu+difftime
      call MPI_REDUCE (acccpu,totcpu,1,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) then
      write (7,711) nx,ny,nz,iproc, jproc
 711  format ('Pencils code, nx,ny,nz=',3i6,'  iproc, jproc=',i4,
     1        ' x ',i6)
#ifdef USE_EVEN
	write (7,"('USE_EVEN: Yes')")
#else
	write (7,"('USE_EVEN: No')")
#endif
#ifndef OPENMP
	num_thr=0
#endif
#ifdef DOUBLE_PREC
	write (7,"('num_thr=',i3,'  double prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
#else
	write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
#endif
#ifndef OPENMP
	num_thr=1
#endif
      write (7,701) numtasks,totcpu/numtasks,'cyl'
 701  format ('total tasks/cpu for initialization=',i6,f12.2,' secs. (',
     1        a4,' data structure)')
      endif


! set u equal to un
 
      do i=1,3+nc
      u(:,:,:,i)=un(:,:,:,i)
      end do

 	if (nsteps.eq.0) then
 	  if (taskid.eq.0) print *,'Writing outpen files before stepping'
	if (flag_outcomm.eq.1) then 
        call outcomm (un)
        else 
        call outpen (un)
        end if
 	  go to 99
 	endif
 
! defined in module timers_comm to time communications	
	i1_comm=0;   i2_comm=0;   i3_comm=0;   i4_comm=0 ; ip_comm=0
	t1_comm=0.;  t2_comm=0.;  t3_comm=0.;  t4_comm=0.; tp1_comm=0.
	t4t_comm=0. ; i4t_comm=0
#ifdef TIMERS	
	t_alltoall=0.
	acccpu2=0.
#endif	
!ccccccccccccccccccccccccccc	
c
	t_rks(:,:)=0.
	t_itrans(:,:)=0.
	t_kxcomm1(:,:)=0.
	ncpusteps=0
	ncpusteps_io=0


	istep=istep0
	ioflag=1
	if (taskid.eq.0) write (6,*) 'before time-stepping loop, kstep=',kstep

 11	istep=istep+1
	jstep=istep-istep0

      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      if (taskid.eq.0) then
         write (6,603) istep,taskid,dmy(2),dmy(1),dmy(3),hms
         write (7,603) istep,taskid,dmy(2),dmy(1),dmy(3),hms
	end if
 603  format ('istep=',i6, ' taskid',i6, ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)

      rtime1=MPI_WTIME()

!	call conchk(un)

 	call shifts
 
#ifdef LAG
c particle initialization call
c
	if (jstep.eq.1) then
	call compart_set 
	if (any(lpgrad.gt.0).and.taskid.eq.1) then
	open (71,file='lag_timings/timing_partvg')
	write (71,"('Detailed timings for operations in sub. partvg (written by task 1)')")
	write (71,"(' istep  kx-xk     udxphy    spxyz_m     nop      nop2     other    overall')")
	end if
#ifdef MOL
	if (nom.gt.0) then
        call init_mol_random
        if(jstep.eq.1.and.istart.gt.0) then
        if (taskid.eq.0) then 
        caux='mranin'
      call fopen1 (817,caux,form)
        end if
      call ranseq (-4,817)
      if (taskid.eq.0) call fclose1 (817)
        endif
	end if
#endif
	call partic (1)
	end if
c
	if (taskid.eq.0.and.jstep.eq.1) then
#ifdef LAGSC2
        write (7,712) nop,nop2
#else
        write (7,712) nop
#endif
 712    format ('Lagrangian run, no. of particles=',2i10)
#ifdef MOL
        write (7,"('Molecules: ngm,nompp,nom=',2i3,i10)")  ngm,nompp,nom
#endif
#ifdef LAG_DOUBLE
	write (7,"('Lagrangian calculation in double precision')")
#else
	write (7,"('Lagrangian calculation in single precision')")
#endif

	end if
#endif

c predictor call for particles
c note that "call partic (2)" is now madefrom sub. realspace



!#Execute single Timestep
!#using Runge Kutta method specified (in input)
      call timestep() 
     
!--------------------------------------------
! At each iostep, check if there is a signal 
! isignal=0: set istop=1 (i.e. write checkpoint and finish execution)
! isignal>0: set nsteps=isignal 
! isignal=-1: read another line and set entime
	if (ioflag.eq.1) then
	  isignal=-2
	  if (taskid.eq.0) then
	    inquire(file='signal',exist=iex)
	    if (iex) then
	      open (999,file='signal')
		read (999,*,end=998,err=998) isignal
		print *, 'file signal exists! isignal = ',isignal
		if (isignal.eq.0) go to 998
		if (isignal.eq.-1) then 
		   read (999,*,end=997,err=997) raux
		   if (raux.le.time) go to 997
		endif   
		go to 998
! 997		print *,'Warning: isignal=-1 but no second line in file signal or entime<time'
 997		print *,'Warning: isignal=-1 but no second line in file'
 	      isignal=0
 998	      close (999)
	    endif
	  endif
	  call MPI_BCAST(isignal,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	  if (isignal.gt.0) then
             if(taskid.eq.0) print *, 'nsteps changed from ',nsteps,' to ',isignal
             nsteps=isignal
	  elseif (isignal.eq.0) then
             if(taskid.eq.0) print *, 'program will write a checkpoint and finish'
             istop=1
	  elseif (isignal.eq.-1) then
             if(taskid.eq.0) print *, 'entime changed to : ',entime
             call MPI_BCAST(entime,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	  endif
       endif
!--------------------------------------------
c force checkpoint after a specified wall clock time

	if (taskid.eq.0) then
	clock_t=MPI_WTIME()
	if (clock_t-clock_start.gt.clock_chkpt.and.clock_chkpt.gt.0.) then
	isflag=1
c	ichkpt=ichkpt+1
c	clock_start=MPI_WTIME()
	if (clock_chkpt.ge.0.5*wall_limit) istop=1
	if (clock_t-clock0+clock_chkpt.ge.wall_limit) istop=1
	write (6,*) 'istep,clock,isflag,istop=',istep,clock_t-clock_start,clock_chkpt,isflag,istop
	clock_start=MPI_WTIME()
	end if
	end if
	call MPI_BCAST (isflag,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (ichkpt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c
        if (isflag.gt.0.and.taskid.eq.0) then
        call chkstf (2)
        if (nc.gt.0.and.scfreq.gt.0) call chkstf (3)
        end if

!----------------------------------------------------
! change entime based on elapsed wall clock time

	if(taskid.eq.0.and.ichktime.eq.1) then
	clock_t = MPI_WTIME()
	if (clock_t-clock0.ge.0.8*wall_limit) then
	if(dtout.ne.0.) then
	stepst = int(time0/dtout)
	iosteps = int((time-time0)/dtout)
	extrasteps = int(iosteps*0.1)
	if(extrasteps.eq.0) extrasteps = 1
	entime = (stepst+iosteps+extrasteps)*dtout
	write(6,*) 'stepst, iosteps, extrasteps = ', stepst,iosteps, extrasteps
	else
	entime = time + 0.1*(time-time0)
	endif
	write(6,*) 'entime changed: time0, entime=', time0, entime
	ichktime = 0
	endif
	endif
        call MPI_BCAST (entime, 1, mpireal, 0, MPI_COMM_WORLD, ierr)


 	
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! checkpointing 
 
      if (isave.gt.0) then
	iaux=mod(istep-istep0,isave)
      else
	iaux=1
      end if
 71   if ((istop.eq.1.and.isave.ge.0) .or.(isflag.gt.0).or.
     1    (isave.gt.0.and.iaux.eq.0)) then
 
	call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c
	isflag=1
 
      if (taskid.eq.0) then
      write (6,*) 'Eulerian checkpointing begins'
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      write (6,621) dmy(2),dmy(1),dmy(3),hms
 621  format ('Eulerian checkpointing: start',
     1        ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
      if (shear.eq.0.) then
      write (6,*) 'velocities and scalars saved at t=',time
      else if (shear.gt.0.) then
      write (6,*) 'velocities and scalars saved at st=',time*shear
      end if
      end if
 
	ichkpt=ichkpt+1
c
      if (taskid.eq.0) call store (2,lustr)
! 7/7/12: Dhawal Buaria : if directive MOL is used than the ranout
! and mranout files are directly written from main.f
! instead of store.f
        caux='ranout'
      if (taskid.eq.0) call fopen1 (817,caux,form)
#ifdef MOL
	if (nom.gt.0) then
        caux='mranout'
      if (taskid.eq.0) call fopen1 (818,caux,form)
	end if
#endif
c
      call ranseq (-2,817)
      if (taskid.eq.0) call fclose1 (817)
#ifdef MOL
	if (nom.gt.0) then
      if (taskid.eq.0) call fclose1 (818)
	end if
#endif
c
	if (taskid.le.1) write (6,*) 'main: before call outpen', taskid
	if (flag_outcomm.eq.1) then 
        call outcomm (un)
        else 
        call outpen (un)
        end if
	if (taskid.eq.0.and.istop.eq.0) call chkstf (1)

      if (taskid.eq.0) then
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      write (6,622) dmy(2),dmy(1),dmy(3),hms
 622  format ('Eulerian checkpointing: ends',
     1        ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
      end if
 
	call time_stamp ('write_timings')
	call write_timings
	call time_stamp ('write_timings')
	call time_stamp1 ('write_iotimings',1)
	call write_iotimings
	call time_stamp1 ('write_iotimings',1)
	call time_stamp1 ('write_iotimings',1)
c
      end if

!  end checkpoint
!cccccccccccccccccccccccccccc

      rtime2=MPI_WTIME()
      tcpu=rtime2-rtime1
      acccpu=acccpu+tcpu
#ifdef TIMERS	
      acccpu2=acccpu2+tcpu
#endif	

      call MPI_REDUCE (tcpu,totcpu,1,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (tcpu,cpumin,1,mpireal,MPI_MIN,0,
     1                 MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (tcpu,cpumax,1,mpireal,MPI_MAX,0,
     1                 MPI_COMM_WORLD,ierr)

      if (taskid.eq.0) then
      avcpu=totcpu/numtasks
      write (7,702) istep,cpumin,avcpu,cpumax
#ifdef IBM
#else
      call flush(7) ! dlw
#endif
 702  format ('cpu at step',i7,
     1            ' :min/ave/max=',3f11.3, '  secs.')
      end if

	! test for unacceptable large t/s/p
	if (taskid.eq.0) then
	if (jstep.ge.3.and.jstep.le.5) then
	  if (cpumax.gt.tsp_max) then
	    write (tsp,fmt="('t/s/p (',f11.3,') greater than specified (',f11.3,')')")
     >            cpumax,tsp_max
	    call abrt(tsp)
	  endif
	endif
	endif

#ifdef LAG
 	if (nop.gt.0.or.nop2.gt.0) then
        iolflag=ioflag
c Fix by D. Buaria, 12/12/2015: the line below is necessary to make
c DNS code count 'npcusteps' properly so that the profiling data
c in *.all files would be written out, when running with fluid
c particles with nonzero dtout
	if( dtout.gt.0. .and. ioflag.eq.-1) iolflag=0
	go to 11
	else
	if (istop.eq.0) go to 11
	end if
#else
	if (istop.eq.0) go to 11
#endif

 99	continue

      call MPI_REDUCE (acccpu,totcpu,1,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)
     
      if (taskid.eq.0) write(6,*)'taskid, cpu/proc used=',
     &   taskid,totcpu/numtasks

#ifdef TIMERS
! timers for communication
      call MPI_Reduce(t1_comm,gt1_comm,1,MPI_REAL8,MPI_MAX,zero,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tp1_comm,gtp1_comm,1,MPI_REAL8,MPI_MAX,zero,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t2_comm,gt2_comm,1,MPI_REAL8,MPI_MAX,zero,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t3_comm,gt3_comm,1,MPI_REAL8,MPI_MAX,zero,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t4_comm,gt4_comm,1,MPI_REAL8,MPI_MAX,zero,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t4t_comm,gt4t_comm,1,MPI_REAL8,MPI_MAX,zero,
     &  MPI_COMM_WORLD,ierr)
	tt_comm = t1_comm + t2_comm + t3_comm + t4_comm + t4t_comm
      call MPI_Reduce(tt_comm,gtt_comm(1),1,MPI_REAL8,MPI_MAX,zero,
     &     MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tt_comm,gtt_comm(2),1,MPI_REAL8,MPI_MIN,zero,
     &     MPI_COMM_WORLD,ierr)

	call time_stamp ('after mpi_reduce')
	if (taskid.eq.0) then
	   write (6,fmt="('xkcomm1: max time(sec), # calls : ',1p,1e13.5,i8)")
     >	   gt1_comm,i1_comm
	   write (6,fmt="('xkcomm2: max time(sec), # calls : ',1p,1e13.5,i8)")
     >	   gt2_comm,i2_comm
	   write (6,fmt="('kxcomm1: max time(sec), # calls : ',1p,1e13.5,i8)")
     >	   gt3_comm,i3_comm
	   write (6,fmt="('kxcomm2: max time(sec), # calls : ',1p,1e13.5,i8)")
     >	   gt4_comm,i4_comm
	   write (6,fmt="('kxcomm2t:max time(sec), # calls : ',1p,1e13.5,i8)")
     >	   gt4t_comm,i4t_comm
	   write (6,fmt="('total in comms: max & min time(sec)   : ',1p,2e13.5)")
     >	   gtt_comm(1),gtt_comm(2)
	endif   

	! write timers/counters to file "timers" 
	call out_timers(acccpu2)
	call time_stamp ('after out_timers')
#endif
c
	if (taskid.eq.0) then
      write (7,604) numtasks,totcpu,totcpu*numtasks
	rwall0 = MPI_WTIME() - rwall0
      write (7,703) rwall0
 604  format ('total cpu over ',i6,' nodes=',f12.3,' secs.',
     1        ' (',f8.2,' s.u.)')
 703  format ('overall wall time for task 0 =',f7.0, ' secs.')
 	end if
c
c write detailed instrumentation timings, added by PKY, 2/3/2012
c moved to a separate routine, 2/23/2012
c	
	if (ichkpt.eq.0) call write_timings
	if (ichkpt.eq.0) call write_iotimings
c
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
 	if (taskid.eq.0.or.taskid.eq.numtasks-1) then
       write (6,602) taskid,dmy(2),dmy(1),dmy(3),hms
 602  format (' DONE, taskid=',i6, ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
 	end if

      call MPI_FINALIZE (mpierr)
c       stop

      end program
!=========================================================
