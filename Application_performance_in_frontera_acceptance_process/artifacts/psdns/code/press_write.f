      subroutine press_write (upx)
c
#ifdef CVC_PRESS

      use comsp
      implicit none
      include 'intvars'

	real(b8) upx(nx,zisz,yjsz)
        real(4), allocatable :: upxsp(:,:)

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
      integer ip2, ixxx, iyyy, izzz
      integer iwrite, nwrite, mwrite, relay
      integer mpistatus(MPI_STATUS_SIZE)
c
      save icall
      data icall/0/
c
	if (taskid.eq.0) write (6,*) 'press_write, istep=',istep
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
      head='press'
      write (char_icall,"(i2)") icall
      call blanks (char_icall,nchar_icall)
        fname=head//char_icall(1:nchar_icall)
c
        call blanks (fname,nf)
      longname='outpress/'//cpid(1:nchar_pid)//fname(1:nf)//'.p'//numer(1:nchar)
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


        allocate (upxsp(nx, zisz))

        open (201,file=longname(1:nchar),form='unformatted')
c
        write(201) nx,zisz,yjsz
c  
        do iyyy=1,yjsz
          do izzz=1,zisz
            do ixxx=1,nx
              upxsp(ixxx,izzz) = upx(ixxx,izzz,iyyy)
            end do
          end do
          write(201) ((upxsp(ixxx,izzz),ixxx=1,nx),izzz=1,zisz)
        end do
c
        close (201)
        deallocate(upxsp)
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
