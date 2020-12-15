	subroutine time_stamp (string)
c
	use mpicom, only: taskid
c
      implicit none

      integer dmy(3),hms(3),date_time(8)
c
	character*(*) string
c
	integer num
!
	num=min(len(string),25)
!
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
c
      if (taskid.eq.0) then
      write (6,601) string(1:num),taskid,dmy(2),dmy(1),dmy(3),hms
 601  format (a25,2x,i5, ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
      end if
c
	return
	end
 
	subroutine time_stamp1 (string,itask)
c
	use mpicom, only: taskid
c
      implicit none

      integer dmy(3),hms(3),date_time(8)
c
	character*(*) string
c
	integer num,itask
!
	num=min(len(string),25)
!
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
c
      if (taskid.eq.itask) then
      write (6,601) string(1:num),taskid,dmy(2),dmy(1),dmy(3),hms
 601  format (a25,2x,i5, ' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
      end if
c
	return
	end
 
