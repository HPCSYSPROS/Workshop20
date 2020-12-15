      subroutine store( icall , lu )
!
! PK Yeung, on Jul13 2008: changes to allow for > 9 checkpoints
!
! MPL version: called by task 0 only, and chkptou and ranout are
!              copied directly onto mass storage from SP2 job
!
!  routine to store, write, read and check data for checkpointing.
!
!  icall = 1:   parameters and input data are stored in arrays.
!
!  icall = 2:   the stored parameters and input data, as well as
!               additional data needed for checkpointing, are written
!               ( for use with restart provided by store(icall=3). )
!
!  icall = 3:   restart from checkpointed data.  input data are read
!               from logical  lu  and checked against current values.
!               if there is a discrepancy and kstop=1 then execution is
!               halted.
!  icall = 4:   similar to icall=3, but reset only the forcing sequence
!
	use comp
#ifdef LAG
	use compart, only: nom
#endif
	implicit none
	integer nid,nrd,nk,icall,i,j,nrc,nic,nic2,nrc2,lu
      parameter ( nid=100,  nk=20 )
      integer :: irdata(nid),isdata(nid)
c
	real, allocatable :: rrdata (:),rsdata(:)
	save rrdata, rsdata
c
c	real :: rrdata(nrd),rsdata(nrd)
      save isdata,nic2,nrc2
      data nic2,nrc2/2*0/
      character*90 string
      character*20 savefn
      character*7 savedir
      data savedir/'./save/'/
	integer nn,nchar,nwrong,mwrong
	character*10 caux
!
	character*2 numer
	integer numc
c
	integer no_err
	logical flag

! this is only to be consistent with files from 1d code	
	integer mz,my,nxhvs,nyvs,nzvs
c
	nrd=2*nrproc+100
	if (.not.allocated(rrdata)) then
	allocate (rrdata(nrd),rsdata(nrd))
	end if
!
	call time_stamp ('enter store')

      if (taskid.le.1)write (6,*) 'store, icall,nrd=',icall,nrd
!---------------  icall=1  ---------------------------------------------
!
      if( icall .eq. 1 ) then
!
!  store specified parameters
!
      nic = 1
      call istow (isdata,nic,nid,  nx, 1 )
      call istow (isdata,nic,nid,  ny, 1 )
      call istow (isdata,nic,nid,  nz, 1 )
      call istow (isdata,nic,nid,  my, 1 )
      call istow (isdata,nic,nid,  mz, 1 )
      call istow (isdata,nic,nid,  nxhvs, 1 )
      call istow (isdata,nic,nid,  nyvs, 1 )
      call istow (isdata,nic,nid,  nzvs, 1 )
      call istow (isdata,nic,nid,  nc, 1 )
      call istow (isdata,nic,nid,  ncd, 1 )
!     call istow (isdata,nic,nid,  nceach, 1 )
!     call istow (isdata,nic,nid,  ncop, 1 )
      call istow (isdata,nic,nid,  kfor, 1 )
!
!  store integers from block data
!
      call istow (isdata,nic,nid,  iostep, 1 )
      call istow (isdata,nic,nid,  kshift, 1 )
!     call istow (isdata,nic,nid,  jfor, 1 )
!
!  store reals from block data
!
      nrc = 1
      call rstow (rsdata,nrc,nrd,  viscos, 1 )
#ifndef NOSCALAR	
      call rstow (rsdata,nrc,nrd,  pr , ncd )
#endif
!     else
!     call rstow (rsdata,nrc,nrd,  pr(1), 1 )
      call rstow (rsdata,nrc,nrd,  a11, 1  )
      call rstow (rsdata,nrc,nrd,  a22, 1  )
      call rstow (rsdata,nrc,nrd,  a33, 1 )
      call rstow (rsdata,nrc,nrd,  shear, 1 )
#ifndef NOSCALAR	
      call rstow (rsdata,nrc,nrd,  grad , 3*ncd )
#endif	
      call rstow (rsdata,nrc,nrd,  beta1, 1 )
      call rstow (rsdata,nrc,nrd,  beta2, 1 )
      call rstow (rsdata,nrc,nrd,  beta3, 1 )
      call rstow (rsdata,nrc,nrd,  cfl, 1 )
      call rstow (rsdata,nrc,nrd,  dt, 1 )
      call rstow (rsdata,nrc,nrd,  tforce, 1 )
      call rstow (rsdata,nrc,nrd,  epsfor, 1 )
      call rstow (rsdata,nrc,nrd,  kforce, 1 )
      call rstow (rsdata,nrc,nrd,  tfiuo, 1 )
!
      nic2 = nic -1
      nrc2 = nrc -1
!
      return
!
!-------------  icall=2  --------------------------------------------
!
      elseif( icall .eq. 2 ) then
      write (6,*) 'store: enter icall=2'
!
      nic = nic2 + 1
      nrc = nrc2 + 1
!
!  store additional integers
!
      call istow (isdata,nic,nid,  istep, 1  )
!
!  store additional reals
!
      call rstow (rsdata,nrc,nrd,  time, 1  )
      call rstow(rsdata,nrc,nrd,  gsh(1,2), 3)
      call rstow(rsdata,nrc,nrd,  b11 , 2 )
      call rstow(rsdata,nrc,nrd,  b22 , 2 )
      call rstow(rsdata,nrc,nrd,  b33 , 2 )
      call rstow(rsdata,nrc,nrd,  b12 , 2 )
      call rstow(rsdata,nrc,nrd,  suo , nrproc )
      call rstow(rsdata,nrc,nrd,  svo , nrproc )
!
!  write data
! rewind file first, thus overwriting data previously saved
!
      nic = nic - 1
      nrc = nrc - 1
!
!	ichkpt=ichkpt+1
!
	caux='chkptou'
      call fopen1 (lu,caux,'formatted  ')
      write(lu,900) nic2,nrc2, nic,nrc
      if( nic .gt. 0 ) write(lu,900) ( isdata(j),j=1,nic )
      if( nrc .gt. 0 ) write(lu,910) ( rsdata(j),j=1,nrc )
      call fclose1 (lu)
!
! write current random numbers on logical unit luran2
!
! 7/7/12: Dhawal Buaria : if directive MOL is used than the ranout
! and mranout files are directly written from main.f

	flag=.true.
#ifdef MOL
	if (nom.gt.0) flag=.false. 
#endif
	if (flag) then
	caux='ranout '
      call fopen1 (luran2,caux,'formatted  ')
      call ranseq (-2,luran2)
      call fclose1 (luran2)
	endif
 
c	ichkpt=ichkpt+1
 
	write (numer,600) ichkpt
 600	format (i2)
	numc=1+floor(alog10(1.*ichkpt))
	caux='chkptou.'//numer(3-numc:2)
	call blanks (caux,nn)
      call fopen1 (lu,caux(1:nn),'formatted  ')
 
      write(lu,900) nic2,nrc2, nic,nrc
      if( nic .gt. 0 ) write(lu,900) ( isdata(j),j=1,nic )
      if( nrc .gt. 0 ) write(lu,910) ( rsdata(j),j=1,nrc )
	close (lu)

      if (isave.gt.0.or.dtsave.gt.0) then
 
	flag=.true.
#ifdef MOL
	if (nom.gt.0) flag=.false. 
#endif
	if (flag) then
! write current random numbers on logical unit luran2
!
	write (numer,600) ichkpt
	numc=1+floor(alog10(1.*ichkpt))
	caux='ranout.'//numer(3-numc:2)
	call blanks (caux,nn)
      call fopen1 (luran2,caux(1:nn),'formatted  ')
!
      call ranseq (-2,luran2)
	close (luran2)
	end if
!
      end if
!
      write (6,*) 'store:  exit icall=2'
      return
!
!---------------  icall=3  ---------------------------------------------
!
      elseif ( icall .eq. 3 ) then
!
!  read and check input data
!

	no_err=0
c
	if (taskid.eq.0) then
      read(lu,900,err=315,end=315) nic2,nrc2, nic,nrc
	end if
c
	call MPI_BCAST (nic2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nrc2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nic,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nrc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (nrc.gt.nrd) then
	deallocate (rrdata)
	allocate (rrdata(nrc))
	end if
c
	if (taskid.eq.0) then
	write (6,*) 'store, nic,nrc=',nic,nrc
      if( nic .gt. 0 ) read(lu,900) ( irdata(j),j=1,nic )
	write (6,*) 'store, after read irdata'
!      if( nrc .gt. 0 ) read(lu,910) ( rrdata(j),j=1,nrc )
#ifdef FEK_FORC
c	read (lu,*) (rrdata(j),j=1,12)
	read (lu,"(5e15.8)") (rrdata(j),j=1,12)
	rrdata(13:nrc)=0.
#else
!     if( nrc .gt. 0 ) read(lu,*) ( rrdata(j),j=1,nrc )
      if( nrc .gt. 0 ) read(lu,"(5e15.8)") ( rrdata(j),j=1,nrc )
#endif
	write (6,*) 'store, after read rrdata'
	end if
c
	go to 320
 315	no_err=1
 320	continue
	call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if (no_err.gt.0) then
        if (taskid.eq.0) write (0,*)
     1      'Stops in store.f due to problems with chkptin file'
        if (taskid.eq.0) write (6,*)
     1      'Stops in store.f due to problems with chkptin file'
	call MPI_FINALIZE (mpierr)
        stop
        end if
c
	call MPI_BCAST (irdata,nic,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (rrdata,nrc,mpireal,0,MPI_COMM_WORLD,mpierr)
      nwrong=0
      do 300 j=1,nic2
      if( irdata(j) .ne. isdata(j) ) then
         nwrong = nwrong + 1
         write(6,*)' integer no. ',j,' wrong in store ',isdata(j),
     1       irdata(j)
      endif
300   continue
!
      mwrong=0
      do 310 j=1,nrc2
      if( rrdata(j) .ne. rsdata(j) ) then
         mwrong = mwrong + 1
         write(6,*)' real no. ',j,' wrong in store '
      endif
310   continue
!
      if( kstop .eq. 1 ) then
         if( nwrong .ne. 0  .or.  mwrong .ne. 0 ) then
            write(6,*)' inconsistency in store '
            write(6,*)' execution stopped '
            stop
         endif
      endif
!
!  reset variables for restart
!
      nic = nic2 + 1
      nrc = nrc2 + 1
!
!  reset integers
!
      call irest (irdata,nic,nid,  istep0, 1 )
!
!  reset reals
!
      call rrest (rrdata,nrc,nrd,  time0, 1  )
      call rrest(rrdata,nrc,nrd,  gsh(1,2) , 3 )
      call rrest(rrdata,nrc,nrd,  b11 , 2 )
      call rrest(rrdata,nrc,nrd,  b22 , 2 )
      call rrest(rrdata,nrc,nrd,  b33 , 2 )
      call rrest(rrdata,nrc,nrd,  b12 , 2 )
      call rrest(rrdata,nrc,nrd,  suo , nrproc )
      call rrest(rrdata,nrc,nrd,  svo , nrproc )
!
! added 1/25/99
      if (istart.lt.0) then
!     istep0=0
!     time0=0.
      end if
!
!---------------  icall=4  ---------------------------------------------
!
      elseif ( icall .eq. 4 ) then
!
!  read and check input data
!
	if (taskid.le.1)write (6,*) 'store, icall=4,taskid=',taskid
	if (taskid.eq.0) then
      read(lu,900) nic2,nrc2, nic,nrc
      if( nic .gt. 0 ) read(lu,900) ( irdata(j),j=1,nic )
      if( nrc .gt. 0 ) read(lu,910) ( rrdata(j),j=1,nrc )
	end if
	call MPI_BCAST (nic2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nrc2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nic,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nrc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (irdata,nic,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (rrdata,nrc,mpireal,0,MPI_COMM_WORLD,mpierr)
!
!  reset only suo and svo for forcing. Note that these numbers
!  are stored at the end of the rrdata array
!
      nrc=nrc-2*nrproc
      call rrest(rrdata,nrc,nrd,  suo , nrproc )
      call rrest(rrdata,nrc,nrd,  svo , nrproc )
!
      endif
!
	if (taskid.le.1)write (6,*) 'exit store, icall=',icall,taskid
	call time_stamp (' exit store')
      return
!
900   format( (10i8) )
910   format( (1p,5e15.8) )
!
      end
