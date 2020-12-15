      subroutine clapf (lu,fn)
c
      character*11 form,unform
      data form,unform/'formatted','unformatted'/
      character*(*) fn
c
c
      close (lu)
      open (lu,file=fn,form=form,iostat=ios,position='append')

        return
        end
