c routine to take out embedded blanks in a character string, and
c keep count of the number of resulting leading non-blank characters
c
      subroutine blanks (char,j)
      character*(*) char
c
      n=len(char)
      j=0
      do 1 i=1,n
      if (char(i:i).ne.' ') then
      j=j+1
      char(j:j)=char(i:i)
      end if
 1    continue
c
      return
      end
