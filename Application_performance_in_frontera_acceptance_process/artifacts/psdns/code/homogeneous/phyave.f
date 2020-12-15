      subroutine phyave (vx,mean)
c
c compute the mean only
c
	use mpicom
	implicit none
	include 'intvars'
c
c
c beware of possible underfloat errors arising from taking fourth
c powers of numbers of very small absolute value
c
	real(b8) :: vx(nx,zisz,yjsz)
c
	double precision sum,sum0
	real mean,nsamp
c
c initialize
c
	sum=0.
c
	
      do 10 yp=1,yjsz
      do 10 zp=1,zisz
c
c sum by pencils
c
      do 20 x=1,nx
c
      sum=sum+vx(x,zp,yp)
c
 20   continue
c
 10   continue
c
c collect tally over all tasks
c
      call MPI_ALLREDUCE (sum,sum0,1,
     1                    MPI_DOUBLE_PRECISION,
     1                    MPI_SUM,MPI_COMM_WORLD,mpierr)
c
c take average, and convert to normalized central moments
c
      nsamp=float(nx)*ny*nz
	mean=sum0/nsamp
c
c
	return
c
      end
