      subroutine istow( iarray , ncount , ndim , ivar , n )
c
c  routine to put the first n elements of the integer array ivar
c    into the integer array iarray, starting at ncount.
c  ncount is incremented, and it is checked that it does not exceed
c    ndim - the dimension of iarray.
c
      dimension iarray(ndim),ivar(n)
c
c  check that there is room
c
      nend = ncount + n - 1
      if( nend .gt. ndim ) then
         write(6,*)' no room in istow '
         write(6,*)' ncount,n,ndim=',ncount,n,ndim
         write(6,*)' stopped in istow '
         stop
      endif
c
      do 10 j = 1 , n
10    iarray( ncount + j - 1) = ivar(j)
c
      ncount = nend + 1
c
      return
      end
