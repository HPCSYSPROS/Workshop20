      subroutine incomm (uny)
	use timers_io
	use comp
	implicit none
	include 'intvars'
        !NPM converted to stride-1 3/6/08
        complex(b8) :: uny(ny,zjsz,xisz,3+nc)

	integer lu,ii,jj,kk,iproc0,jproc0,ip,i
	integer, allocatable :: ixfst(:),ixlst(:),jzfst(:),jzlst(:)
	integer slabpenc
	integer irz0
	integer iread,mread,nread,relay,jsource
	integer mpistatus (MPI_STATUS_SIZE)
c
	integer ind_prec
	logical single_flag
	logical double_flag

	integer no_err
c
        character*100 fn1,fn2,fn
        logical ex1,ex2
        integer ip2,ndirs,nchar
c
        complex(b8), allocatable :: buf(:,:,:)
c
        integer mys
        logical flag_mys
        

#ifdef TIMESTAMP
	character*20 auxiliar
#endif

        logical bwexs
        integer idummy,bwin,lubw

        if (taskid.eq.0) write (6,*) 'enter incomm'

        bwin=0
        if (taskid.eq.0) then
        inquire(file='bwio',exist=bwexs)
        write (6,*) 'bwexs=',bwexs
        if (bwexs) then
        lubw=819
        open(unit=lubw,file='bwio',action='read')
        read (lubw,'(1x)')
        read (lubw,*) bwin,idummy
        close(lubw)
        endif
        endif
        call MPI_BCAST (bwin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (taskid.le.1) write (6,*) 'bwin=',bwin


#ifdef TIMERS_IO	  
	tread_io = 0.
	iread_io = 0
#endif
c
        luinit(1)=3
c
	if (kinit(1).lt.0) go to 90
c
#ifdef DOUBLE_PREC
	single_flag=.false.
	ind_prec=2
	inquire (file='single',exist=single_flag)
	if (single_flag) ind_prec=1
#else
	double_flag=.false.
	ind_prec=1
	inquire (file='double',exist=double_flag)
	if (double_flag) ind_prec=2
#endif

	if (taskid.eq.0) write (6,*) 'incomm,kinit=',taskid,kinit(1),ind_prec

        if (irz.gt.0) then
            slabpenc=2
        if (taskid.eq.0) then
            call blanks (indir_fn,nchar)
            fn1=indir_fn(1:nchar)//'/'//fninit(1)//'.c0'
            fn2=indir_fn(1:nchar)//'/0/'//fninit(1)//'.c0'
            call blanks (fn,nchar)
        inquire (file=fn1(1:nchar),exist=ex1)
            call blanks (fn2,nchar)
        inquire (file=fn2(1:nchar),exist=ex2)
        write (6,*) 'incomm: fn1=',fn1
        write (6,*) 'incomm: fn2=',fn2
        write (6,*) 'incomm: ex1,ex2=',ex1,ex2
        end if

           
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
              write (6,*) 'incomm:  after reading indata.grid file: iproc0,jproc0',iproc0,jproc0
              close (300)
           end if
c
	go to 320
 310	no_err=1
 320	continue
           call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (no_err.gt.0) then
	if (taskid.eq.0) write (0,*) 
     1      'Stops in incomm.f due to problems with indata.grid file'
	if (taskid.eq.0) write (6,*) 
     1      'Stops in incomm.f due to problems with indata.grid file'
	stop
	end if

           call MPI_BCAST (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
           call MPI_BCAST (jproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if (taskid.eq.0) then
        ndirs=0
        if (ex2) ndirs=2**(ip2(jproc0)/2)
        end if
           call MPI_BCAST (ndirs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

           call MPI_BCAST (ixfst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
           call MPI_BCAST (ixlst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
           call MPI_BCAST (jzfst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
           call MPI_BCAST (jzlst,irz,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)


	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
#ifdef TIMESTAMP
           auxiliar='incomm, after bcast of indata.grid info'
           call time_stamp (auxiliar)
#endif
        else !(irz.lt.0)

           irz=iabs(irz)
           slabpenc=1
c
	ndirs=0
           
           allocate (ixfst(irz),ixlst(irz))
           allocate (jzfst(irz),jzlst(irz))
           
           iproc0 = 1 
           do kk=1,irz
              ixfst(kk)=1
              ixlst(kk)=nxh 
              jzfst(kk)=(kk-1)*nz/irz+1
              jzlst(kk)=kk*nz/irz
           end do
        end if

        if (taskid.eq.0) then
        mys=1
        inquire (file='mys',exist=flag_mys)
        if (flag_mys) then
        open (300,file='mys')
        read (300,*) mys
        close (300)
        end if
        write (6,*) 'incomm: value of mys=',mys
        end if
        call MPI_BCAST (mys,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

c relay scheme
c
        

       if (jproc.le.1024) mread=jproc
       if (mod(jproc,1024).eq.0) mread=1024

	nread=jproc/mread
	iread=jpid/mread
        if (taskid.eq.0) write (6,*) 'incomm: batch size=',mread


	if (iread.eq.0) then

           relay=taskid
           call incomm_row (uny,iproc0,jproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,bwin,relay,taskid, mys)
	else
        jsource=taskid-mread*iproc
           call MPI_RECV (relay,1,MPI_INTEGER,jsource,jsource,
     1                 MPI_COMM_WORLD,mpistatus,mpierr)
           call incomm_row (uny,iproc0,jproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,bwin,relay,taskid, mys)
	end if

	if (iread.lt.nread-1) then
           call MPI_SSEND (relay,1,MPI_INTEGER,taskid+mread*iproc,taskid,
     1                MPI_COMM_WORLD,mpierr)
	end if
c
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	deallocate (ixfst,ixlst,jzfst,jzlst)

c these lines originally written by Dmitry in inpen_read.f,
c moved here by P.K.Y / D.D. on 12/03/09
c
c Now redistribute the data based on cylindrical data allocation

      allocate(buf(nz,xisz,yjsz))

      do i=1,3+nc

c First transpose the data from Y-pencils to Z-pencils
c
        if (i.eq.1.and.taskid.eq.0) then
        write (6,*) 'incomm: uny(2,2,2,1)=',uny(2,2,2,1)
        end if
c
         call kxcomm1_pen(uny(1,1,1,i),buf,1)

        if (i.eq.1.and.taskid.eq.0) then
        write (6,*) 'incomm: buf(2,2,2)=',buf(2,2,2)
        end if
c
c Move from (z,x) space to cyl. variable (num_al) space
c and transpose back to Y cylinders

         if(num_al_i(0) .gt. 0) then
         call xkcomm2_sq2cyl(buf,uny(1,1,1,i),1)
        end if
c
        if (i.eq.1.and.taskid.eq.0) then
        write (6,*) 'incomm: uny(2,2,2,1)=',uny(2,2,2,1)
        end if

c
      enddo

      deallocate(buf)
c



 90	return
	end
