      subroutine rrest( rarray , ncount , ndim , rvar , n )
c
c  routine to fill the first n elements of the real array rvar
c    from the real array rarray, starting at ncount.
c  ncount is incremented.
c    ndim is the dimension of rarray.
c
      dimension rarray(ndim),rvar(n)
c
      do 10 j = 1 , n
10    rvar(j) = rarray( ncount +j -1 )
c
      ncount = ncount +n
c
      return
      end
