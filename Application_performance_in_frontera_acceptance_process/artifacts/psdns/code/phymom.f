      subroutine phymom (vx,i1,i2)
c
	use comsp
	implicit none
	include 'intvars'
c
	real(b8) :: vx(nx,zisz,yjsz,nu)
c
	integer i,i1,i2,m,msglen
	real(b8), allocatable :: sumy(:,:),sqterm(:),momt(:,:),momt0(:,:)
	real(b8) nsamp
c
c calculation of mean, variance, skewness coeff. and flatness factor
c for field variables i1 to i2
c
c this version used in the dns code and initialization program
c
c beware of possible underfloat errors arising from taking fourth
c powers of numbers of very small absolute value
c
	allocate (sumy(4,i1:i2),stat=ierr)
	allocate (momt(4,i1:i2),stat=ierr)
	allocate (momt0(4,i1:i2),stat=ierr)
	allocate (sqterm(nx))
c
c initialize
c
	momt(:,:)=0.
c
      do 10 yp=1,yjsz
c
c sum by z-planes
c
	sumy(:,:)=0.
c
      do 20 i=i1,i2
      do 20 zp=1,zisz
      do 20 x=1,nx
c
      sqterm(x)=vx(x,zp,yp,i)*vx(x,zp,yp,i)
      sumy(1,i)=sumy(1,i)+vx(x,zp,yp,i)
      sumy(2,i)=sumy(2,i)+sqterm(x)
      sumy(3,i)=sumy(3,i)+vx(x,zp,yp,i)*sqterm(x)
      sumy(4,i)=sumy(4,i)+sqterm(x)*sqterm(x)
c
 20   continue
c
      do i=i1,i2
      do m=1,4
      momt(m,i)=momt(m,i)+sumy(m,i)
      end do
      end do
c
 10   continue
c
c collect tally over all tasks, and let task 0 finish the job,
c and broadcast the variances back to the other tasks
c
      msglen=4*4*(i2-i1+1)
c
      call MPI_ALLREDUCE (momt(1,i1),momt0(1,i1),msglen/4,mpireal,MPI_SUM,
     1                    MPI_COMM_WORLD,mpierr)
c
c take average, and convert to normalized central moments
c
      nsamp=float(nx)*ny*nz
c
      do 30 i=i1,i2
c
      imean(i)=0 
      umean(i)=0.
c
      momt(1,i)=momt0(1,i)/nsamp
      momt(2,i)=momt0(2,i)/nsamp
      momt(3,i)=momt0(3,i)/nsamp
      momt(4,i)=momt0(4,i)/nsamp
c
c subtract off the mean if zero-mean is not assumed
c
      if (imean(i).ne.0) then
      umean(i)=momt(1,i)
      momt(4,i)=momt(4,i)-4.*umean(i)*momt(3,i)+6.*umean(i)**2*momt(2,i)
     1          -3.*umean(i)**4
      momt(3,i)=momt(3,i)-3.*umean(i)*momt(2,i)+2.*umean(i)**3
      varce(i)=momt(2,i)-umean(i)**2
      else
      umean(i)=0.
      varce(i)=momt(2,i)
      end if
c
      if (varce(i).gt.0.) then
      skcof(i)=momt(3,i)/varce(i)**1.5
      flat(i)=momt(4,i)/varce(i)**2
      else if (varce(i).lt.0.) then
      write (6,*) 'negative variance: i,varce(i)=',i,varce(i)
      write (6,*) 'unexpected error in sub. phymom: aborts'
      stop
      end if
c
      if (taskid.eq.0) write (6,601) i,umean(i),varce(i),skcof(i),flat(i)
 601  format ('moments:',i3,1p,4e14.4)
c
 30   continue
c
	deallocate (sqterm,momt,momt0,sumy)
c
      return
      end
