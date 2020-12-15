      subroutine irest( iarray , ncount , ndim , ivar , n )
c
c  routine to fill the first n elements of the integer array ivar
c    from the integer array iarray, starting at ncount.
c  ncount is incremented.
c    ndim is the dimension of iarray.
c
      dimension iarray(ndim),ivar(n)
c
      do 10 j = 1 , n
10    ivar(j) = iarray( ncount +j -1 )
c
      ncount = ncount +n
c
      return
      end
