	subroutine logspace(d1,d2,n,out,iopt)
c
c returns a logspaced vector (out). Used for variable intervals PDFs.
c d1: first point 10^d1
c d2: last point 10^d2
c n: number of points
c iopt=1 : one sided. i.e. n points from 10^d1 to 10^d2
c iopt=2 : two sided. i.e. n points from -10^d2 to 10^d2. Last point
c          before zero is 10^d1. (zero is also included).
c          n should be odd.
c
	integer n
	real d1,d2,out(n)
	if (iopt.eq.1) then 
	  do i=0,n-2
	    out(i+1)=10.**(d1+i*(d2-d1)/(n-1))
	  end do
	  out(n)=10.**d2
      else if (iopt.eq.2) then 
	  out(1)=-10.**d2
	  k=2
	  do i=(n-1)/2-2,0,-1
	    out(k)=-10.**(d1+i*(d2-d1)/((n-1)/2-1))
	    k=k+1
	  end do
	  out(k)=0
	  k=k+1
	  do i=0,(n-1)/2-2
	    out(k)=10.**(d1+i*(d2-d1)/((n-1)/2-1))
	    k=k+1
	  end do
	  out(k)=10.**d2
c	  print *,k,n
	end if
	return 
	end
