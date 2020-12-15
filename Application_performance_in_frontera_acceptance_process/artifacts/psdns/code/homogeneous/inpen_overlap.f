      subroutine inpen_overlap (uny)

c A few changes initated by, and based on input from D. Whitaker
c of Cray, 3/23/2012: allow longer character strings in fn1,fn2
c and generalizing relay mechanism to cases where numtasks
c may not be multiple of 4096 but would be a multiple of 6144
c
	use timers_io
	use comp
	implicit none
	include 'intvars'
        !NPM converted to stride-1 3/6/08
        complex(b8) :: uny(ny,zjsz,xisz,3+nc)

	integer lu,ii,jj,kk,iproc0,jproc0,ip
	integer, allocatable :: ixfst(:),ixlst(:),jzfst(:),jzlst(:)
	integer slabpenc
	integer irz0
	integer iread,mread,nread,relay
	integer mpistatus (MPI_STATUS_SIZE)
c
        integer ind_prec
        logical single_flag
        logical double_flag

	character*256 fn
#ifdef SUBDIR
	character*256 fn1,fn2
	logical ex1,ex2
	integer ndirs,ip2
#endif
	integer jjn,nxin,nyin,nzin,nchar

	integer no_err
c
c 5 lines below added to allow operation with io_scheme=2,
c corresponding to suggestion from Mike Matheson, Aug 2012
c
	integer itask,msgtag,ncount,xi1,zj1,io_scheme,mread0,nvar,imsg
	integer handle1,handle2
	integer, allocatable :: xist_all(:),zjst_all(:)
	logical exs
	complex(b8), allocatable :: unbuf(:)
	integer mpist1(MPI_STATUS_SIZE)
	integer mpist2(MPI_STATUS_SIZE)
c
#ifdef INPUT_PEN
	complex(b8), allocatable :: buf(:,:,:)
	integer i
#endif

#ifdef TIMESTAMP
	character*20 auxiliar
#endif

        logical bwexs
        integer idummy,bwin,lubw

        if (taskid.eq.0) write (6,*) 'enter inpen'

        bwin=0
        if (taskid.eq.0) then
        inquire(file='bwio',exist=bwexs)
	if (bwexs) then
        lubw=819
        open(unit=lubw,file='bwio',action='read')
        read (lubw,'(1x)')
        read (lubw,*) bwin,idummy
        close(lubw)
	endif
        endif
        if (bwin.ne.0) call MPI_BCAST (bwin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


	if (taskid.eq.0) write (6,*) 'enter inpen_overlap'
#ifdef TIMERS_IO	  
	tread_io = 0.
	iread_io = 0
#endif

c	if (kinit(1).lt.0) go to 90
	if (all(kinit.lt.0)) go to 90

	allocate (xist_all(0:numtasks-1),zjst_all(0:numtasks-1))
	call MPI_ALLGATHER (xist,1,MPI_INTEGER,xist_all,1,MPI_INTEGER,
     1                   MPI_COMM_WORLD,mpierr)
	call MPI_ALLGATHER (zjst,1,MPI_INTEGER,zjst_all,1,MPI_INTEGER,
     1                   MPI_COMM_WORLD,mpierr)

	call time_stamp ('enter inpen')
#ifdef DOUBLE_PREC
        single_flag=.false.
        ind_prec=2
        inquire (file='single',exist=single_flag)
        if (single_flag) ind_prec=1
#else
        double_flag=.false.
        ind_prec=1
        inquire (file='double',exist=double_flag)
        if (double_flag) ind_prec=1
#endif

c
#ifdef SUBDIR
	if (taskid.eq.0) then
	    call blanks (indir_fn,nchar)
            fn1=indir_fn(1:nchar)//'/'//fninit(1)//'.p0'
            fn2=indir_fn(1:nchar)//'/0/'//fninit(1)//'.p0'
	    call blanks (fn1,nchar)
	inquire (file=fn1(1:nchar),exist=ex1)
	    call blanks (fn2,nchar)
	inquire (file=fn2(1:nchar),exist=ex2)
	write (6,*) 'inpen: ex1,ex2=',ex1,ex2
	ndirs=0
	if (ex2) ndirs=2**(ip2(irz)/2)
	end if
	call MPI_BCAST (ndirs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
c
	if (taskid.eq.0) then
	    call blanks (indir_fn,nchar)
            fn=indir_fn(1:nchar)//'/'//fninit(1)//'.p0'
            if (ex2) fn=indir_fn(1:nchar)//'/0/'//fninit(1)//'.p0'
	    call blanks (fn,nchar)
	lu=luinit(1)
	open (lu,file=fn(1:nchar),form='unformatted',action='read')
	read (lu) jjn,nxin,nyin,nzin
	if (taskid.eq.0) write (6,*) 'nzin=',nzin
	close (lu)
	end if
	call MPI_BCAST (nxin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nzin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

        if (irz.gt.0) then
           slabpenc=2
           
           allocate (ixfst(irz),ixlst(irz))
           allocate (jzfst(irz),jzlst(irz))
c
	no_err=0

           if (taskid.eq.0) then
              open (300,file='indata.grid',err=310)
              read (300,*,err=310,end=310) irz0,iproc0,jproc0              
              kk=0
              do ii=0,iproc0-1
                 do jj=0,jproc0-1
                    kk=kk+1
      read (300,400,err=310,end=310) ixfst(kk),ixlst(kk),jzfst(kk),jzlst(kk)
                 end do
              end do
 400	format (6x,4i6)
              write (6,*) 'inpen:  after reading indata.grid file'
              close (300)
           end if
c
	go to 320
 310	no_err=1
 320	continue
           call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (no_err.gt.0) then
	if (taskid.eq.0) write (0,*) 
     1      'Stops in inpen.f due to problems with indata.grid file'
	if (taskid.eq.0) write (6,*) 
     1      'Stops in inpen.f due to problems with indata.grid file'
	stop
	end if
c
	if (taskid.eq.0) write (6,*) 'inpen: before bcast iproc0'

c           call MPI_BCAST (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c           call MPI_BCAST (ixfst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c           call MPI_BCAST (ixlst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c           call MPI_BCAST (jzfst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
c           call MPI_BCAST (jzlst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

c
         call MPI_BCAST (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'inpen:  after bcast iproc0'
c
c the following changes  were suggested by Tommy Minyard, TACC, 1/30/10
c
C  Split the following broadcasts into 4K chunk sizes
C       call MPI_BCAST (ixfst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
C       call MPI_BCAST (ixlst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
C       call MPI_BCAST (jzfst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
C       call MPI_BCAST (jzlst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

         jj = irz / 1024
         do ii = 1, jj
            call MPI_BCAST (ixfst((ii-1)*1024+1),1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BCAST (ixlst((ii-1)*1024+1),1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BCAST (jzfst((ii-1)*1024+1),1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BCAST (jzlst((ii-1)*1024+1),1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         enddo

C  Do the remainder of the bcast message if not multiple of 1024 or < 1024

         if ( jj*1024 .lt. irz ) then
            call MPI_BCAST (ixfst(jj*1024+1),irz-jj*1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BCAST (ixlst(jj*1024+1),irz-jj*1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BCAST (jzfst(jj*1024+1),irz-jj*1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BCAST (jzlst(jj*1024+1),irz-jj*1024,
     &          MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         endif
	if (taskid.eq.0) write (6,*) 'inpen:  after all broadcasts'


        else !(irz.lt.0)

           irz=iabs(irz)
           slabpenc=1
           
           allocate (ixfst(irz),ixlst(irz))
           allocate (jzfst(irz),jzlst(irz))
           
           iproc0 = 1 
           do kk=1,irz
              ixfst(kk)=1
              ixlst(kk)=nxin/2
              jzfst(kk)=(kk-1)*nzin/irz+1
              jzlst(kk)=kk*nzin/irz
           end do
        end if

c allow for possibility of nz > nzin
	do kk=1,irz
	if (jzfst(kk).gt.nzin/2) jzfst(kk)=jzfst(kk)+nz-nzin
	if (jzlst(kk).gt.nzin/2) jzlst(kk)=jzlst(kk)+nz-nzin
	end do
c
	if (taskid.eq.0) then
	io_scheme=1
	inquire (file='input.io',exist=exs)
	if (exs) then
	open (4,file='input.io')
	read (4,"(1x)")
	read (4,*) io_scheme
	if (io_scheme.eq.1) then
	read (4,"(1x)")
	read (4,*) mread0
	else 
	mread0=0
	end if
	close (4)
	
	end if
	write (6,*) 'io_scheme=',io_scheme
	end if

	call MPI_BCAST (io_scheme,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (mread0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	
        !
!----------------------------------------------------
	if (io_scheme.eq.1) then

	if (mread0.eq.0) then
c
 	 if (numtasks.le.4096) mread=numtasks
 	if (mod(numtasks,4096).eq.0) mread=4096
c	if (numtasks.le.8192) mread=numtasks
c	if (mod(numtasks,8192).eq.0) mread=8192
	if (mod(numtasks,6144).eq.0) mread=6144
c
	else
	
	mread=mread0
c
	end if

        nread=numtasks/mread
        iread=taskid/mread

	if (iread.eq.0) then
           relay=taskid
#ifdef READR
           if (bwin.eq.1) then
           call inpen_read_bw (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xist,zjst,taskid)
           else
           call inpen_read_r (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xist,zjst,taskid)
           endif
#else
           call inpen_read (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec)
#endif
	else
           call MPI_RECV (relay,1,MPI_INTEGER,taskid-mread,taskid-mread,
     1                 MPI_COMM_WORLD,mpistatus,mpierr)
#ifdef READR
           if (bwin.eq.1) then
           call inpen_read_bw (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xist,zjst,taskid)
           else
           call inpen_read_r (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xist,zjst,taskid)
           endif
#else
           call inpen_read (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec)
#endif
	end if

	if (iread.lt.nread-1) then
           call MPI_SSEND (relay,1,MPI_INTEGER,taskid+mread,taskid,
     1                MPI_COMM_WORLD,mpierr)
	end if
c
!----------------------------------------------------
	else if (io_scheme.eq.2) then

c	mread=iproc
c	nread=jproc
c	iread=taskid/iproc
c	if (taskid.eq.0) write (6,*) 'inpen: numtasks,mread=',
c     1            numtasks,mread
c
c	ncount=ny*zjsz*xisz*3
c	if (nc.gt.0.and.kinit(3+nc).gt.0) ncount=ncount*(3+nc)/3
cc
c	if (jpid.eq.0) then
cc
c	do itask=taskid+numtasks-iproc,taskid,-iproc
c	msgtag=taskid*10000+itask
c	xi1=xist_all(itask)
c	zj1=zjst_all(itask)
c	write (6,*) 'inpen, send : taskid,itask,msg=',taskid,itask,msgtag
c           call inpen_read_r (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xi1,zj1,itask)
c		if (itask.gt.taskid) then
c           call MPI_SSEND (uny,ncount,mpicomplex,itask,msgtag,
c     1                MPI_COMM_WORLD,mpierr)
c		end if
c	end do
c	
c	else
c	msgtag=ipid*10000+taskid
c	write (6,*) 'inpen, recv : taskid,ipid,msg=',taskid,ipid,msgtag
c           call MPI_RECV (uny,ncount,mpicomplex,ipid,msgtag,
c     1                MPI_COMM_WORLD,mpierr)
c
c	end if
c	
	mread=jproc
	nread=iproc
	iread=taskid/jproc
	relay=0
	if (taskid.gt.0) relay=numtasks

	nvar=3
	if (nc.gt.0.and.kinit(3+nc).gt.0) nvar=3+nc
	ncount=ny*zjsz*xisz*nvar
c
	if (ipid.eq.0) then
	allocate (unbuf(ny*zjsz*xisz*nvar))
c
	imsg=0
c
	do 50 itask=taskid+iproc-1,taskid,-1
c
	imsg=imsg+1
	msgtag=taskid*100+itask-taskid
	xi1=xist_all(itask)
	zj1=zjst_all(itask)
c
	if (imsg.gt.2)  then
	if (mod(itask-taskid,2).ne.0) then
	call MPI_WAIT (handle1,mpist1,mpierr)
	else
	call MPI_WAIT (handle2,mpist2,mpierr)
	end if
	end if
c
	if (mod(itask-taskid,2).ne.0) then

           if (bwin.eq.1) then
           call inpen_read_bw (unbuf,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xi1,zj1,itask)
           else
           call inpen_read_r (unbuf,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xi1,zj1,itask)
           endif
 
		else
           call inpen_read_r (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xi1,zj1,itask)
		end if
c
	if (itask.gt.taskid) then
	if (mod(itask-taskid,2).ne.0) then
           call MPI_ISEND (unbuf,ncount,mpicomplex,itask,msgtag,
     1                MPI_COMM_WORLD,handle1,mpierr)
	else
           call MPI_ISEND (uny,ncount,mpicomplex,itask,msgtag,
     1                MPI_COMM_WORLD,handle2,mpierr)
	end if
	end if
c
 50	continue

	
	else
	itask=(taskid/iproc)*iproc
	msgtag=itask*100+taskid-itask
           call MPI_RECV (uny,ncount,mpicomplex,itask,msgtag,
     1                MPI_COMM_WORLD,mpistatus,mpierr)

	end if
	
	if (ipid.eq.0) deallocate (unbuf)
	end if
c
c
!----------------------------------------------
c
	deallocate (ixfst,ixlst,jzfst,jzlst)



#ifdef INPUT_PEN
c
c these lines originally written by Dmitry in inpen_read.f,
c moved here by P.K.Y / D.D. on 12/03/09
c
c Now redistribute the data based on cylindrical data allocation

      allocate(buf(nz,xisz,yjsz))
	if (taskid.eq.0) write (6,*) 'inpen: after allocate'

      do 100 i=1,3+nc
c
        if (kinit(i).le.0) go to 100

c First transpose the data from Y-pencils to Z-pencils

         call kxcomm1_pen (uny(1,1,1,i),buf,1)

c Move from (z,x) space to cyl. variable (num_al) space
c and transpose back to Y cylinders

         if(num_al_i(0) .gt. 0) then
           call xkcomm2_sq2cyl (buf,uny(1,1,1,i),1)
	end if
c
 100	continue

      deallocate(buf)
c
#endif
	call time_stamp ('exit inpen')
c	call MPI_FINALIZE (mpierr)
c	stop
c	write (6,902) taskid,num_al_i(0)
c902	format (' exit inpen, taskid, num_al_i(0)=',i5,i8)

	deallocate (xist_all,zjst_all)

 90	return
	end
