	subroutine prni(char)
	use com
	implicit none
!	integer :: in
!	integer :: ia(in)
!	integer :: numc,nchar
	character(*) char

!	nchar=len(char)
!     if (taskid.eq.0) then
!     numc=1
!     else
!     numc=1+floor(alog10(1.*taskid))
!     end if

 	write (6,400) taskid,char
 400	format('[',i5,'] ',A60)

!400	format('[',i<numc>,'] ',A<nchar>)
!400	format(i<numc>)

 	return
	end
