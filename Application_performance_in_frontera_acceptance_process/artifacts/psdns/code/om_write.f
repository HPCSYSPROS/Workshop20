      subroutine om_write (omx)
c
#ifdef VORPHY

      use comsp
      implicit none
      include 'intvars'

	real(b8) omx(nx,zisz,yjsz,3)
        real(4), allocatable :: omxsp(:,:)
c
	real(b8) sum,mean,term1,term2,term3

      character*6 numer
      character*6 head
c
      character*40 longname
      character*4 cpid
      character*2 char_icall
      integer nchar_pid,nchar_icall
c
      character*8 fname
      integer nchar,nf,icall
      integer iwgid,ndirs
      integer ip2, ixxx, iyyy, izzz, icomp
      integer iwrite, nwrite, mwrite, relay
      integer mpistatus(MPI_STATUS_SIZE)
c
      save icall
      data icall/0/

c
	if (taskid.eq.0) write (6,*) 'om_write: istep=',istep

	sum=0.
	do yp=1,yjsz
	do zp=1,zisz
	do x=1,nx
	term1=omx(x,zp,yp,1)*omx(x,zp,yp,1)
	term2=omx(x,zp,yp,2)*omx(x,zp,yp,2)
	term3=omx(x,zp,yp,3)*omx(x,zp,yp,3)
	sum=sum+term1+term2+term3
	end do
	end do
	end do
	sum=sum/nx/zisz/yjsz
c
	call MPI_REDUCE (sum,mean,1,mpireal,MPI_SUM,0,
     1                   MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) then
	mean=mean/numtasks
	write (6,*) 'om_write: mean vortsq=', mean,vvij(1,1)+vvij(2,2)+vvij(3,3)
	end if

c
c    write 3d vorticity field
c

      write (numer,"(i6)") taskid
      call blanks (numer,nchar)
c
      ndirs=2**(ip2(numtasks)/2)
      iwgid=mod(taskid,ndirs)
c
      write (cpid,"(i3,'/')") iwgid
      call blanks (cpid,nchar_pid)
c
      icall=icall+1
      head='omega'
      write (char_icall,"(i2)") icall
      call blanks (char_icall,nchar_icall)
        fname=head//char_icall(1:nchar_icall)
c
        call blanks (fname,nf)
      longname='outom/'//cpid(1:nchar_pid)//fname(1:nf)//'.p'//numer(1:nchar)
      call blanks(longname, nchar)

c relay scheme
c
        if (numtasks.le.4096) mwrite=numtasks
        if (mod(numtasks,4096).eq.0) mwrite=4096
        if (mod(numtasks,6144).eq.0) mwrite=6144
        if (mod(numtasks,12288).eq.0) mwrite=12288
        if (mod(numtasks,16384).eq.0) mwrite=16384
        if (mod(numtasks,24576).eq.0) mwrite=24576
c
        nwrite=numtasks/mwrite
        iwrite=taskid/mwrite
c
        if (iwrite.eq.0) then
          relay=taskid
        else
          call MPI_RECV (relay,1,MPI_INTEGER,taskid-mwrite,taskid-mwrite,
     1                   MPI_COMM_WORLD,mpistatus,mpierr)
        end if


        allocate (omxsp(nx, zisz))

        open (201,file=longname(1:nchar),form='unformatted')
c
        write(201) nx,zisz,yjsz
c
        do icomp=1,3
          do iyyy=1,yjsz
            do izzz=1,zisz
              do ixxx=1,nx
                omxsp(ixxx,izzz) = omx(ixxx,izzz,iyyy,icomp)
              end do
            end do
            write(201) ((omxsp(ixxx,izzz),ixxx=1,nx),izzz=1,zisz)
          end do
        end do
c
        close (201)
        deallocate(omxsp)
c
        if (iwrite.lt.nwrite-1) then
      call MPI_SSEND (relay,1,MPI_INTEGER,taskid+mwrite,taskid,
     1                MPI_COMM_WORLD,mpierr)
        end if
c
        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
c
c

#endif
	return
	end
