      subroutine fclose1 (lu)
c
c close the logical unit lu and report the closing
c (as an aid to further development of code)
c
c used in MPL, same as fclose before
c
      close (lu)
      if (lu.ne.21) write (6,600) lu
 600  format ('fclose: lu=',i4)
c
      return
      end
