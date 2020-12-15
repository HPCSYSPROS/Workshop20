      subroutine input
!
        use bouss_module
	use com
	implicit none
	integer lu,i,j,nchar
	character*120 caux
	real, allocatable :: rsendbuf(:)
	integer, allocatable :: isendbuf(:)
	!integer :: isendbuf(5)
!
	lu=1111

	if (taskid.eq.0) write (6,*) 'enter input'

	allocate (isendbuf(30),rsendbuf(30))
	isendbuf=0
	rsendbuf=0.

!--------------------------------------------------------
	if (taskid.eq.0) then ! only task=0 opens/reads input

 	  open(lu,file='input')

	  read (lu,fmt='(1x)')
        read (lu,*) nx,ny,nz,nc
#ifdef NOSCALAR
	  nc=0
#endif
	isendbuf(1)=nx; isendbuf(2)=ny; isendbuf(3)=nz; isendbuf(4)=nc
!
	read (lu,fmt='(1x)')
      read (lu,*) nsteps,iostep,istart,kstop,isave
	isendbuf(5)=nsteps; isendbuf(6)=iostep; 
	isendbuf(7)=istart; isendbuf(8)=kstop; 
	isendbuf(9)=isave; 
!
	read (lu,fmt='(1x)')
      read (lu,*) stend,entime,time0
	rsendbuf(1)=stend; rsendbuf(2)=entime; rsendbuf(3)=time0;
!
      allocate (kinit(3+nc))
	read (lu,fmt='(1x)')
      read (lu,*) (kinit(i),i=1,3)
	isendbuf(10)=kinit(1); isendbuf(11)=kinit(2)
	isendbuf(12)=kinit(3)
!
      allocate (fninit(3+nc))
	read (lu,fmt='(1x)')
      read (lu,*) fninit(1)
	do i=2,3+nc
	fninit(i)=fninit(1)
	enddo
!
      read (lu,fmt='(1x)')
      read (lu,fmt='(a120)') caux
      call blanks(caux,nchar)
      indir_fn=caux(1:nchar)

	read (lu,fmt='(1x)')
      read (lu,*) irz
	isendbuf(13)=irz
!
	read (lu,fmt='(1x)')
      read (lu,*) viscos
	rsendbuf(4)=viscos
!
	read (lu,fmt='(1x)')
      read (lu,*) a11,a22,a33,shear
	rsendbuf(5)=a11
	rsendbuf(6)=a22
	rsendbuf(7)=a33
	rsendbuf(8)=shear
!
	read (lu,fmt='(1x)')
      read (lu,*) beta1,beta2,beta3,bet12
	rsendbuf(9)=beta1
	rsendbuf(10)=beta2
	rsendbuf(11)=beta3
	rsendbuf(12)=bet12
!
	read (lu,fmt='(1x)')
      read (lu,*) cfl,dt
	rsendbuf(13)=cfl; rsendbuf(14)=dt
!
	read (lu,fmt='(1x)')
      read (lu,*) dtout,tfout
	rsendbuf(15)=dtout; rsendbuf(16)=tfout
!
	read (lu,fmt='(1x)')
      read (lu,*) dtsave,tfsave
	rsendbuf(17)=dtsave; rsendbuf(18)=tfsave
      
	read (lu,fmt='(1x)')
      read (lu,*) kshift
	isendbuf(14)=kshift
!
#ifndef FEK_FORC
	read (lu,fmt='(1x)')
      read (lu,*) tforce,epsfor,kforce,tfiuo
	rsendbuf(19)=tforce; rsendbuf(20)=epsfor
	rsendbuf(21)=kforce; rsendbuf(22)=tfiuo
#else	
	kforce=0.
      read (lu,fmt='(1x)')
!      read (lu,*) kf_shell
	allocate (ekinf_input(nx/2))
	ekinf_input(:)=0.
      read (lu,*) kf_shell, (ekinf_input(i),i=1,kf_shell+1) 
#endif
!
	read (lu,fmt='(1x)')
	read (lu,*) rrate
	rsendbuf(23)=rrate
!
	read (lu,fmt='(1x)')
	read (lu,*) iovor
	isendbuf(15)=iovor
!
	read (lu,fmt='(1x)')
	read (lu,*) seed_input
	write (6,*) 'input: seed_input=',seed_input
	isendbuf(16)=seed_input

	read (lu,fmt='(1x)')
	read (lu,*) rkmethod
	write (6,*) 'rkmethod=',rkmethod
	isendbuf(17)=rkmethod

#ifdef BGW
	read (lu,fmt='(1x)')
	read (lu,*) mrcpin
	isendbuf(18)=mrcpin
	read (lu,fmt='(1x)')
	read (lu,*) nbl_x
	isendbuf(19)=nbl_x
#endif
!	read (lu,*,end=99) dims(1), dims(2)


 99	close(lu)

 	endif ! if(taskid.eq.0)
!
      luran1=19 ; luran2=20
      kranf1=5  ;     kranf2=5
	ksran=1
!--------------------------------------------------------

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
	  iovor=isendbuf(15)
	  seed_input=isendbuf(16)
          rkmethod=isendbuf(17)
#ifdef BGW
	  mrcpin=isendbuf(18)
	  nbl_x=isendbuf(19)
#endif
	endif
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
	endif
	call MPI_BCAST(indir_fn,120,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
      if (taskid.ne.0) allocate (fninit(3+nc))
	do i=1,3+nc
	  nchar=len(fninit(i))
	  call MPI_BCAST(fninit(i),nchar,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
	enddo

#ifdef FEK_FORC
	call MPI_BCAST (kf_shell,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
#endif
 	deallocate (isendbuf,rsendbuf)
!----------------------------------------------
! scalars
#ifndef NOSCALAR

!
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
!
	read (lusc,fmt='(1x)')
      read (lusc,*) (pr(i),i=1,nc)
!
	read (lusc,fmt='(1x)')
      read (lusc,*) ((grad(j,i),j=1,3),i=1,nc)

	read (lusc,fmt='(1x)')
      read (lusc,*) (kinit(i+3),i=1,nc)
!
	read (lusc,fmt='(1x)')
      read (lusc,*) scfreq

	read (lusc,fmt='(1x)')
      read (lusc,*) gpdflim

#ifdef BUOY
	read (lusc,fmt='(1x)')
      read (lusc,*) bvfsq
#endif
#ifdef BUOY2
	read (lusc,fmt='(1x)')
      read (lusc,*) bvfsq
#endif
#ifdef BIF
	read (lusc,fmt='(1x)')
      read (lusc,*) bvfsq
#endif

!# these need to be turned into a subroutine call (hook bouss)
#ifdef BOUSS
	read (lusc,fmt='(1x)')
      read (lusc,*) froude,noise_ratio,noise_spconst,noise_spslope
#ifdef LINDBORG
	read (lusc,fmt='(1x)')
      read (lusc,*) visc_h, visc_v, difs_h,difs_v
	read (lusc,fmt='(1x)')
      read (lusc,*) F_h,ei_rate,ei_waveno
#endif
#endif

      do i=1,nc
      kran(i)=i+1
      end do

	close (lusc)
!
	end if ! taskid.eq.0

	! broadcasts values in input.sc
	call MPI_BCAST(pr,nc,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(grad,3*nc,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(kinit(4),nc,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(scfreq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(kran,nc,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST(kcps,nc,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

#ifdef BUOY
	call MPI_BCAST (bvfsq,1,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif
#ifdef BUOY2
	call MPI_BCAST (bvfsq,1,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif
#ifdef BIF
	call MPI_BCAST (bvfsq,1,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif
#ifdef BOUSS
	call MPI_BCAST (froude,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (noise_ratio,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (noise_spconst,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (noise_spslope,1,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif

#ifdef LINDBORG
	call MPI_BCAST (visc_h,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (visc_v,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (difs_h,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (difs_v,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (F_h,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (ei_rate,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (ei_waveno,1,mpireal,0,MPI_COMM_WORLD,mpierr)
	pi=atan(1.)*4.
	ei_scale=2.*pi/ei_waveno
#endif

	end if ! nc.gt.0
#endif
!----------------------------------------------
	
	lustr=2
	if (taskid.eq.0) write (6,*) 'exit input'

 
      return
      end
!
