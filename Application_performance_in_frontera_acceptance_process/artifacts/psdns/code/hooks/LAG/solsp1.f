c specialsed tridagonal solver for periodic spline:
c operations on coefficient matrix:
c diagonal elements: all equal to 4
c superdiagonal and subdiagonal: all equal to 1
c
c ref. the theory of splines and their applications
c      ahlberg et al. 1967, chap.2
c
      subroutine solsp1 (p,q,t,n,s,d)
c
c "n" is no. of rows/columns of matrix
c p,q,t are arrays to be passed to sub. solsp2
c s is working array
c
	use mpicom
	implicit none
      real(b8)   p(n),q(0:n),s(0:n),t(n)
	real(b8) d
	integer i,n
c
      q(0)=0.
      do 10 i=1,n
      p(i)=q(i-1)+4.
      q(i)=-1./p(i)
 10   continue
c
      s(0)=1.
      do 20 i=1,n
 20   s(i)=-s(i-1)/p(i)
c
      t(n)=1.
      do 30 i=n-1,1,-1
 30   t(i)=q(i)*t(i+1)+s(i)
c
      d=1./(t(1)+t(n-1)+4.)
c
c take reciprocal of p, for later efficiency
c
      do 40 i=1,n
 40   p(i)=1./p(i)
c
      return
      end
