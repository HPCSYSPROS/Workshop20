      subroutine input
c
#ifdef BOUSS
	use stratification_module
#endif
	use com
#ifdef MODEL_SPECTRUM
      USE spectrum,ONLY: kmaxetaSpec, boxint, kol, cL, intlen, eta, eps
#ifdef SPECTRAL_TRUNCATION
      USE spectrum,ONLY: kLower, kUpper
#endif
#endif
#ifdef ERF_VELOCITY
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      USE ErfVelocity,ONLY: ErfVelocitySetup
#endif
      USE IO

	implicit none
	integer lu,i,j,nchar
	character*256 caux
	real(b8), allocatable :: rsendbuf(:)
	integer, allocatable :: isendbuf(:)
	!integer :: isendbuf(5)
c
	real tl_hat,kmaxeta
!
	integer no_err, no_end
c
	character*1 tail
	integer nchar2,nchar3
	logical exs

#if defined MODEL_SPECTRUM || defined ERF_VELOCITY
      INTEGER :: stat, exsInt
#endif
#ifdef ERF_VELOCITY
      REAL(KIND=REAL64) :: a_erf, b_erf, c_erf, d_erf
      REAL(KIND=REAL64) :: x_start, x_finish
      REAL(KIND=REAL64) :: tol_erf
#endif
      integer :: hdf5_inp,chkp_num,hdf5_out,bwio_flg
c
	no_err=0
	no_end=0
	
	lu=1111

	if (taskid.eq.0) write (6,*) 'enter input'

	allocate (isendbuf(30),rsendbuf(30))
	isendbuf=0
	rsendbuf=0.

!need to add variables to com module 

!--------------------------------------------------------
	if (taskid.eq.0) then ! only task=0 opens/reads input

 	  open(lu,file='input')

	  read (lu,fmt='(1x)',err=10,end=20)
        read (lu,*,err=10,end=20) nx,ny,nz,nc
	write (6,*) 'after read nx...'
#ifdef NOSCALAR
	  nc=0
#endif
	isendbuf(1)=nx; isendbuf(2)=ny; isendbuf(3)=nz; isendbuf(4)=nc
!
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) nsteps,iostep,istart,kstop,isave
	write (6,*) 'after read nsteps...'
	isendbuf(5)=nsteps; isendbuf(6)=iostep; 
	isendbuf(7)=istart; isendbuf(8)=kstop; 
	isendbuf(9)=isave; 
!
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) stend,entime,time0
	write (6,*) 'after read stend...'
	rsendbuf(1)=stend; rsendbuf(2)=entime; rsendbuf(3)=time0;
!
      allocate (kinit(3+nc))
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) (kinit(i),i=1,3)
	isendbuf(10)=kinit(1); isendbuf(11)=kinit(2)
	isendbuf(12)=kinit(3)
!
      allocate (fninit(3+nc))
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) fninit(1)
	write (6,*) 'after read fninit(1)...'
!
      read (lu,fmt='(1x)',err=10,end=20)
      read (lu,fmt='(a256)') caux
	write (6,*) 'after read indir_fn...'
      call blanks(caux,nchar)
      indir_fn=caux(1:nchar)

! lines added 2/24/13 to detect identify last checkpoint
! (this is to enhance robustness of sequence of chained jobs,
! in case a prior run did not reach its own last intended checkpoint)
c
	caux=fninit(1)
	call blanks (caux,nchar3) 
	if (nchar3.ne.7) go to 53

        if (fninit(1)(7:7).ne.'c')  then
c
	if (kinit(1).gt.0.and.istart.gt.0) then
	do ichkpt=9,1,-1
	write (tail,"(i1)") ichkpt
	caux=indir_fn(1:nchar-7)//'/chkptou.'//tail
	call blanks (caux,nchar2)
	inquire (file=caux(1:nchar2),exist=exs) 
	write (6,*) ichkpt,caux(1:nchar2),exs
	if (exs) go to 51
	end do
 51	continue
	if (exs) then
	fninit(1)(7:7)=tail
	write (6,*) 'revised fninit=',fninit(1)
	write (6,*) caux(1:nchar2)
	end if
 52	continue
	end if
c
        end if
c
 53	continue
c
	ichkpt=0
	
	do i=2,3+nc
	fninit(i)=fninit(1)
	enddo
	

	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) irz
	write (6,*) 'after read irz...'
	isendbuf(13)=irz
!
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) viscos
	write (6,*) 'after read viscos...'
	rsendbuf(4)=viscos
!
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) a11,a22,a33,shear
	write (6,*) 'after read a11...'
	rsendbuf(5)=a11
	rsendbuf(6)=a22
	rsendbuf(7)=a33
	rsendbuf(8)=shear
!
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) beta1,beta2,beta3,bet12
	write (6,*) 'after read beta1...'
	rsendbuf(9)=beta1
	rsendbuf(10)=beta2
	rsendbuf(11)=beta3
	rsendbuf(12)=bet12

	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) cfl,dt
	write (6,*) 'after read cfl...'
	rsendbuf(13)=cfl; rsendbuf(14)=dt

	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) dtout,tfout
	write (6,*) 'after read dtout...'
	rsendbuf(15)=dtout; rsendbuf(16)=tfout

	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) dtsave,tfsave
	write (6,*) 'after read dtsave...'
	rsendbuf(17)=dtsave; rsendbuf(18)=tfsave
      
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) kshift
	write (6,*) 'after read kshift...'
	isendbuf(14)=kshift

#ifndef FEK_FORC
	read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) tforce,epsfor,kforce,tfiuo
	write (6,*) 'after read tforce...'
	rsendbuf(19)=tforce; rsendbuf(20)=epsfor
	rsendbuf(21)=kforce; rsendbuf(22)=tfiuo
#else	
	kforce=0.
      read (lu,fmt='(1x)',err=10,end=20)
	allocate (ekinf_input(nx/2))
	ekinf_input(:)=0.
      read (lu,*,err=10,end=20) kf_shell, (ekinf_input(i),i=1,kf_shell+1) 
	write (6,*) 'after read kf_shell...'
#endif
!
	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) rrate
	write (6,*) 'after read rrate...'
	rsendbuf(23)=rrate
!
#ifdef CVC_PRESS
	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) icvc,npout
	isendbuf(15)=icvc
	write (6,*) 'after read icvc...'
#else
	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) iovor
	isendbuf(15)=iovor
	write (6,*) 'after read iovor...'
#endif
!
	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) seed_input
	write (6,*) 'after read seed_input...',seed_input
	isendbuf(16)=seed_input

	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) rkmethod
	write (6,*) 'after read rkmethod...'
	write (6,*) 'rkmethod=',rkmethod
	isendbuf(17)=rkmethod

	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) nxpad,nypad,nzpad
	write (6,*) 'read nxpad...'
	isendbuf(18)=nxpad
	isendbuf(19)=nypad
	isendbuf(20)=nzpad

! The following two parameters are used only if the code
! is run with the EPFOR and SBP_2010 directives but w/o FEK_FORC
	read (lu,fmt='(1x)',err=10,end=20)
	read (lu,*,err=10,end=20) tl_hat,kmaxeta
	write (6,*) 'input: tl_hat,kmaxeta=',tl_hat,kmaxeta
	rsendbuf(24)=tl_hat
	rsendbuf(25)=kmaxeta

! Read HDF5 file IO parameters.
      read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) hdf5_inp,chkp_num
      read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) hdf5_out,bwio_flg
      isendbuf(21)=hdf5_inp
      isendbuf(22)=chkp_num
      isendbuf(23)=hdf5_out
      isendbuf(24)=bwio_flg
      read (lu,fmt='(1x)',err=10,end=20)
      read (lu,*,err=10,end=20) ioaxi(:)
      isendbuf(25:27)=ioaxi(1:3)
	

! For MHD runs, read electrical conductivity
#ifdef MHD
        read (lu,fmt='(1x)')
        read (lu,*) conduc
	rsendbuf(26)=conduc
#endif

 99	close(lu)

 	endif ! if(taskid.eq.0)
!
      luran1=19 ; luran2=20
      kranf1=5  ;     kranf2=5
	ksran=1
!--------------------------------------------------------
c
c PK Yeung, 1/25/10: make all tasks stop if task 0 encounters
c error in reading 'input'
c
	go to 30
 10	no_err=1
	go to 30
 20	no_end=1
 30	continue
	call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (no_end,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (no_err.gt.0) then
	if (taskid.eq.0) write (6,*) 
     1          'Code stops because of error in reading ''input'' file'
	stop
	end if
	if (no_end.gt.0) then
	if (taskid.eq.0) write (6,*) 
     1   'Code stops because of end-of-file in reading ''input'' file'
	stop
	end if

c-----------------------------------------------------------------

! broadcast all values in input file
	call MPI_BCAST(isendbuf,30,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (taskid.ne.0) then
	  nx=isendbuf(1); ny=isendbuf(2)
	  nz=isendbuf(3); nc=isendbuf(4)
	  nsteps=isendbuf(5); iostep=isendbuf(6)
	  istart=isendbuf(7); kstop=isendbuf(8)
	  isave=isendbuf(9)
        allocate (kinit(3+nc))
	  kinit(1)=isendbuf(10); kinit(2)=isendbuf(11); 
	  kinit(3)=isendbuf(12); 
	  irz=isendbuf(13)
	  kshift=isendbuf(14)
#ifdef CVC_PRESS
	  icvc=isendbuf(15)
#else
	  iovor=isendbuf(15)
#endif
	  seed_input=isendbuf(16)
          rkmethod=isendbuf(17)
          nxpad=isendbuf(18)  !x_pad ratio for 3/2 dealiasing
          nypad=isendbuf(19)  !y
          nzpad=isendbuf(20)  !z
	endif
      ! This is put outside of the if block so root also performs these
      ! operations.
      hdf5_inp=isendbuf(21)
      chkp_num=isendbuf(22)
      hdf5_out=isendbuf(23)
      bwio_flg=isendbuf(24)
      ioaxi(1:3)=isendbuf(25:27)
      if (hdf5_inp.gt.0) then
         hdf5_input_flag=.true.
         if (hdf5_inp.eq.1) hdf5_input_type=PENCIL_HDF5_IO
         if (hdf5_inp.eq.2) hdf5_input_type=ROWCOM_HDF5_IO
         hdf5_init=chkp_num
      else
         hdf5_input_flag=.false.
      endif
      if (hdf5_out.gt.0) then
         hdf5_output_flag=.true.
         if (hdf5_out.eq.1) hdf5_output_type=PENCIL_HDF5_IO
         if (hdf5_out.eq.2) hdf5_output_type=ROWCOM_HDF5_IO
         if (bwio_flg.eq.1) then
            hdf5_output_bwio=BWIO_HDF5_ON
         else
            hdf5_output_bwio=BWIO_HDF5_OFF
         endif
      else
         hdf5_output_flag=.false.
      endif
#ifdef CVC_PRESS
	call MPI_BCAST (npout,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	iovor=npout
#endif
	call MPI_BCAST(rsendbuf,30,mpireal,0,MPI_COMM_WORLD,mpierr)
	if (taskid.ne.0) then
	  stend=rsendbuf(1); entime=rsendbuf(2)
	  time0=rsendbuf(3)
	  viscos=rsendbuf(4)
	  a11=rsendbuf(5); a22=rsendbuf(6)
	  a33=rsendbuf(7); shear=rsendbuf(8)
	  beta1=rsendbuf(9)
	  beta2=rsendbuf(10)
	  beta3=rsendbuf(11)
	  bet12=rsendbuf(12)
	  cfl=rsendbuf(13)
	  dt=rsendbuf(14)
	  dtout=rsendbuf(15)
	  tfout=rsendbuf(16)
	  dtsave=rsendbuf(17)
	  tfsave=rsendbuf(18)
	  tforce=rsendbuf(19)
	  epsfor=rsendbuf(20)
	  kforce=rsendbuf(21)
	  tfiuo=rsendbuf(22)
	  rrate=rsendbuf(23)
	tl_hat=rsendbuf(24)
	kmaxeta=rsendbuf(25)
#ifdef MHD
	conduc=rsendbuf(26)
#endif
	endif
c       call MPI_BCAST(indir_fn,120,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(indir_fn,256,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
      if (taskid.ne.0) allocate (fninit(3+nc))
	do i=1,3+nc
	  nchar=len(fninit(i))
	  call MPI_BCAST(fninit(i),nchar,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
	enddo

#ifdef FEK_FORC
	call MPI_BCAST (kf_shell,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
#ifdef SHELL_DK
 	kforce=kf_shell*amin1(beta1,beta2,beta3)
#else
 	kforce=kf_shell
#endif
#endif
 	deallocate (isendbuf,rsendbuf)
!----------------------------------------------
! scalars
#ifndef NOSCALAR

	if (nc.gt.0) then

	allocate (pr(nc),kran(nc),kcps(nc),stat=ierr)
	allocate (grad(3,nc),stat=ierr)

      nceach=nc
	ncps=nc

      kcps(:)=1
!
	luscn=52
	luscd=53

	if (taskid.eq.0) then

	lusc=1113
 	open(lusc,file='input.sc')
!
	write (6,*) 'nc=',nc
!
c 	read (lusc,fmt='(1x)')
 	read (lusc,fmt='(1x)',err=40,end=50)
c      read (lusc,*) (pr(i),i=1,nc)
       read (lusc,*,err=40,end=50) (pr(i),i=1,nc)
 	write (6,*) 'after read pr...'
!
	read (lusc,fmt='(1x)',err=40,end=50)
      read (lusc,*,err=40,end=50) ((grad(j,i),j=1,3),i=1,nc)
	write (6,*) 'after read grad...'

	read (lusc,fmt='(1x)',err=40,end=50)
      read (lusc,*,err=40,end=50) (kinit(i+3),i=1,nc)
	write (6,*) 'after read kinit...'
!
	read (lusc,fmt='(1x)',err=40,end=50)
      read (lusc,*,err=40,end=50) scfreq
	write (6,*) 'after read scfreq...'

	read (lusc,fmt='(1x)',err=40,end=50)
      read (lusc,*,err=40,end=50) gpdflim
	write (6,*) 'after read pdflim...'

#ifdef LVSTEP
c (defaults to 1 if not read)
	ivstep=1
	read (lusc,fmt='(1x)',err=71,end=71)
      read (lusc,*,err=40,end=50) ivstep
	write (6,*) 'after read ivstep..., ivstep=',ivstep
 71	continue
#endif

#ifdef BIF
	read (lusc,fmt='(1x)',err=40,end=50)
      read (lusc,*,err=40,end=50) bvfsq
	write (6,*) 'after read bvfsq...'
#endif

#ifdef BOUSS
        read (lusc,fmt='(1x)',err=40,end=50)
      read (lusc,*,err=40,end=50) froude,noise_ratio,noise_spconst,noise_spslope
	write (6,*) 'after read froude...'
#endif




      do i=1,nc
      kran(i)=i+1
      end do

	close (lusc)

	end if ! taskid.eq.0

c PK Yeung, 1/25/10: make all tasks stop if task 0 encounters
c error in reading 'input'
c
	go to 60
 40	no_err=1
	go to 60
 50	no_end=1
 60	continue
	call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (no_end,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (no_err.gt.0) then
	if (taskid.eq.0) write (6,*) 
     1          'Code stops because of error in reading ''input.sc'' file'
	stop
	end if
	if (no_end.gt.0) then
	if (taskid.eq.0) write (6,*) 
     1   'Code stops because of end-of-file in reading ''input.sc'' file'
	stop
	end if
c
	! broadcasts values in input.sc
	call MPI_BCAST(pr,nc,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(grad,3*nc,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(kinit(4),nc,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(scfreq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(kran,nc,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(kcps,nc,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (gpdflim,1,mpireal,0,MPI_COMM_WORLD,mpierr)

#ifdef LVSTEP
	call MPI_BCAST (ivstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (ivstep.gt.1) then
	if (nc.le.0) then
	if (taskid.eq.0) write (6,*) 
     1      'ivstep > 1 should not be used with no scalars'
	go to 73
	end if
	if (kshift.le.2) then
	if (taskid.eq.0) write (6,*) 
     1      'ivstep > 1 should not be used with kshift < 3'
	go to 73
 	end if
 	end if
	go to 72
 73	call MPI_FINALIZE (mpierr)
	stop 'Code stop in sub. input, due to problems with LVSTEP'
 72	continue
#endif
c
#ifdef BOUSS
        call MPI_BCAST (froude,1,mpireal,0,MPI_COMM_WORLD,mpierr)
        call MPI_BCAST (noise_ratio,1,mpireal,0,MPI_COMM_WORLD,mpierr)
        call MPI_BCAST (noise_spconst,1,mpireal,0,MPI_COMM_WORLD,mpierr)
        call MPI_BCAST (noise_spslope,1,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif


#ifdef BIF
	call MPI_BCAST (bvfsq,1,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif

	end if ! nc.gt.0
#endif
!----------------------------------------------

#ifdef ROTD
	iraxis=3
	if (taskid.eq.0) then
	inquire (file='input.rot',exist=exs)
	if (exs) then
	open (1114,file='input.rot')
	read (1114,*) 
	read (1114,*) iraxis
	end if
	end if
	call MPI_BCAST (iraxis,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.1) write (6,*) 'input: iraxis=',iraxis
#endif 
#ifdef MHD
	imaxis=1
	if (taskid.eq.0) then
	inquire (file='input.mhd',exist=exs)
	if (exs) then
	open (1114,file='input.mhd')
	read (1114,*) 
	read (1114,*) imaxis
	end if
	end if
	call MPI_BCAST (imaxis,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.1) write (6,*) 'input: imaxis=',imaxis
#endif

#ifdef MODEL_SPECTRUM
      ! We only attempt to read and communicate the model spectrum parameters if
      ! kinit(1)=-1, which indicates a spectrum function should be used.
      IF (kinit(1) .EQ. -1) THEN
         stat = 0
         ! The root process opens and reads the spectrum parameters.
         IF (taskid .EQ. 0) THEN
            INQUIRE(FILE='input.spectrum',EXIST=exs)
            IF (exs) THEN
               ! Read the data.
               OPEN(UNIT=1115,FILE='input.spectrum',STATUS='OLD',ACTION='READ')
               READ(1115,*,IOSTAT=stat) kmaxetaSpec, boxint, kol
               CLOSE(UNIT=1115)
               !
               ! The parameter exsInt is zero if the file exists.
               exsInt = 0
            ELSE
               ! Set exsInt to a positive value to indicate error.
               exsInt = 1
            END IF
         END IF
         !
         ! Broadcast logical checks to see if program should terminate.
         CALL MPI_BCAST(exsInt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         CALL MPI_BCAST(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         !
         ! If the input file did not exist, set default parameters.
         IF (exsInt .GT. 0) THEN
            IF (taskid .EQ. 0) THEN
               WRITE(*,*) 'input.spec does not exist, setting default values.'
               kmaxetaSpec = 1.0
               boxint = 12.0
               kol = 1.62
            END IF
         END IF
         !
         ! If there were IO errors, terminate execution.
         IF (stat .GT. 0) THEN
            IF (taskid .EQ. 0) THEN
               WRITE(*,*) 'Error reading data from input.spec.'
            END IF
            CALL MPI_FINALIZE (mpierr)
            STOP
         END IF
         IF (stat .LT. 0) THEN
            IF (taskid .EQ. 0) THEN
               WRITE(*,*) 'End of input.spec file reached.'
            END IF
            CALL MPI_FINALIZE (mpierr)
            STOP
         END IF
         !
         ! Broadcast the parameters to the other processes.
         IF (ALLOCATED(rsendbuf)) THEN
            DEALLOCATE(rsendbuf)
         END IF
         ALLOCATE(rsendbuf(3))
         IF (taskid .EQ. 0) THEN
            rsendbuf(1) = kmaxetaSpec
            rsendbuf(2) = boxint
            rsendbuf(3) = kol
         END IF
         CALL MPI_BCAST(rsendbuf, 3, mpireal, 0, MPI_COMM_WORLD, mpierr)
         !
         ! Unpack the data into the spectrum module variables.
         kmaxetaSpec = rsendbuf(1)
         boxint = rsendbuf(2)
         kol = rsendbuf(3)
         !
         ! Calculate the other quantities required for the model spectrum. NOTE:
         ! right now the calculation of kmax is hardcoded in the calulcation of
         ! eta. This should be changed in the future. NOTE: we determine the size
         ! of the integral scale in the shortest dimension of the box, which
         ! corresponds to the largest value of beta supplied.
         cL = (1.262*kol)**3
         eta = kmaxetaSpec/(SQRT(2.)/3.*MIN(nxpad*beta1,nypad*beta2,nzpad*beta3))
         intlen = 2.0*ACOS(-1.0)/boxint/MAX(beta1, beta2, beta3)
         eps = viscos**3/eta**4
#ifdef SPECTRAL_TRUNCATION
         ! Set the cutoff wavenumbers used for non-cubic grids.
         kLower = MAX(beta1, beta2, beta3)
         kUpper = MIN(0.5*nxpad*beta1, 0.5*nypad*beta2, 0.5*nzpad*beta3)
#endif
         ! Free temporary memory.
         DEALLOCATE(rsendbuf)
      END IF
#endif /* MODEL_SPECTRUM */

#ifdef ERF_VELOCITY
      IF (taskid .EQ. 0) THEN
         INQUIRE(FILE='input.erf',EXIST=exs)
         IF (exs) THEN
            ! Read the data.
            OPEN(UNIT=1115,FILE='input.erf',STATUS='OLD',ACTION='READ')
            READ(1115,*)
            READ(1115,*) a_erf, b_erf, c_erf, d_erf
            READ(1115,*)
            READ(1115,*) x_start, x_finish
            READ(1115,*)
            READ(1115,*) tol_erf
            CLOSE(UNIT=1115)
            !
            ! The parameter exsInt is zero if the file exists.
            exsInt = 0
         ELSE
            ! Set exsInt to a positive value to indicate error.
            exsInt = 1
         END IF
         !
         ! Set up the error function module.
         CALL ErfVelocitySetup(x_start, x_finish, a_erf, b_erf, c_erf, d_erf,
     1                         tol_erf, beta1, beta2, beta3)
      END IF
      !
      ! Broadcast logical checks to see if program should terminate.
      CALL MPI_BCAST(exsInt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      !
      ! If there were IO errors, terminate execution.
      IF (exsInt .GT. 0) THEN
         IF (taskid .EQ. 0) THEN
            WRITE(*,*) 'File input.erf not present. Terminating.'
         END IF
         CALL MPI_FINALIZE (mpierr)
         STOP
      END IF
#endif

	
#ifdef EPFOR
#ifdef SBP_2010
	call epfor_param (tl_hat,kmaxeta)
#endif
#endif
c
#ifdef LVSTEP
	if (nc.eq.0) ivstep=1
#endif

	lustr=2
	if (taskid.eq.0) write (6,*) 'exit input'
 
      return
      end
