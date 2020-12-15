      subroutine rstow( rarray , ncount , ndim , rvar , n )
c
c  routine to put the first n elements of the real array rvar
c    into the real array rarray, starting at ncount.
c  ncount is incremented, and it is checked that it does not exceed
c    ndim - the dimension of rarray.
c
      dimension rarray(ndim),rvar(n)
c
c  check that there is room
c
      nend = ncount + n - 1
      if( nend .gt. ndim ) then
         write(6,*)' no room in rstow '
         write(6,*)' ncount,n,ndim=',ncount,n,ndim
         write(6,*)' stopped in rstow '
         stop
      endif
c
      do 10 j = 1 , n
10    rarray( ncount + j - 1) = rvar(j)
c
      ncount = nend + 1
c
      return
      end
