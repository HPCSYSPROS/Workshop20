      subroutine order
c
c routine to define ordering of function derivatives
c treated in sub. popsp
c
c should be invoked in 1st call of sub. partsp
c
#ifdef LAG
c
	use compart
	implicit none
c
c iod(i,j)=order of differentiation in ith direction for the jth
c          variable (0 :function, 1: 1st derivative, etc)
c
c e.g., iod(.,j)=0,1,1 gives dydz(function)
c
	iod(:,:)=0
c
c 1st partial derivatives
c
      iod(1,2)=1
      iod(2,3)=1
      iod(3,4)=1
c
c mixed partial derivatives
c
      iod(1,8)=1
      iod(2,8)=1
      iod(1,9)=1
      iod(3,9)=1
      iod(2,10)=1
      iod(3,10)=1
c
c 2nd (normal) partial derivatives
c method: csa, i.e. by differentiating spline fitted to 1st
c         derivative field
c
      iod(1,5)=1
      iod(2,6)=1
      iod(3,7)=1
c
#endif
      return
      end
