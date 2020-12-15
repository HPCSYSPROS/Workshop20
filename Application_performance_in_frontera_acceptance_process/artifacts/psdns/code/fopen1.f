      subroutine fopen1 (lu,fn,fmt)
c
c to open file named fn as logical unit lu,
c with fmt as formatted/unformatted option
c
	use mpicom, only: taskid
!      include 'mpif.h'
c
c     character*7 fn,nam
      character*(*) fn
      character*40 nam
      character*11 fmt
      logical opn
      data opn/ .false. /
c
c check if the targetted logical unit is already opened
c
      inquire (lu, opened=opn, name=nam)
c
c if the associated internal filename is the same as the argument fn,
c then just rewind the file
c
c take out any embedded blanks
c
      call blanks (fn,jj)
      call blanks (nam,jjnam)
c
      if (opn.and.nam(1:jjnam).eq.fn(1:jj)) then
      rewind (lu)
      return
      end if
c
c gives a warning if the associated internal filename is not
c the same as the argument fn
c
      if (opn.and.nam(1:jjnam).ne.fn(1:jj)) then
      write (6,*) 'warning: logical unit',lu,'is already connected'
      write (6,*) 'internal filename was   ',nam
      end if
c
      open (lu, file=fn(1:jj), form=fmt, iostat=ios)
      if (taskid.le.1) 
     1write (6,*) 'taskid,lu,fopen1: filename=',taskid,lu,fn(1:jj)
c
      if (ios.gt.0) then
      write (6,*) 'error condition in opening logical unit',lu
      write (6,900) fn(1:jj),fmt
      write (6,*) 'stopped in fopen1'
      write (6,*) 'does the outpen directory exist?'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
 900  format (' fn, fmt=',a10,a14)
      return
      end
