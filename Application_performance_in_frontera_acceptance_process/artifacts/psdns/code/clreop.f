      subroutine clreop (lu,fn,fmt)
c
c close a file, copy its current contents to the savedir directory,
c and re-open it with append option (works for xlf version 3 only)
c
c this version on systems that do not allow system calls (e.g. BG, XT3)
c
      character*(*) fn
      character*11 fmt
      character*40 string,nam
      logical opn
c
c take out any embedded blanks in filename
c
	if (lu.eq.0) then
	write (6,*) 'clreop called with lu=0 for fn=',fn
	return
	end if
c
	write (6,*) 'enter clreop,fn=',fn
	close (lu)

c	call blanks (fn,jj)
c      open (lu, file=fn(1:jj), form=fmt, iostat=ios, position='append')
      open (lu, file=fn, form=fmt, iostat=ios, position='append')
	write (6,*) 'clreop: lu,fn=',lu,fn
c
      if (ios.gt.0) then
      write (6,*) 'error condition in re-opening logical unit',lu
 	call blanks (fn,jj)
      write (6,900) fn(1:jj),fmt
      write (6,*) 'stopped in clreop'
      stop
      end if
	write (6,*) ' exit clreop,fn=',fn
c
 900  format (' fn, fmt=',a10,a14)
      return
      end
