!NPM converting to stride-1 3/9/09 (not finished)
	subroutine force_eirate (uny,m,rate1)
	
#ifdef LINDBORG
	
	use comsp

	implicit none
	include 'intvars'
	
	integer :: m
	complex(b8) :: uny(ny,zjsz,xisz,3+nc)	  
	real :: term1,term2,kzabs
	real sumh,sumv,sumh0,sumv0,total,rate1
	
	integer :: xst
	
	tfact(:)=2.
	tfact(1)=1.

	sumh=0.
	sumv=0.

	
	if (zjst.eq.1) then
	   zp=1
	   do 30 y=1,ny
	      if (y.eq.nyhp) go to 30
	      xst=1
	      if (y.eq.1) xst=2
	      do 32 xp=xst,xisz
		 x=xp+xist-1
		 term1=uny(y,zp,xp,1)*conjg(ufh(xp,y,1))*b11(m)
		 term2=uny(y,zp,xp,2)*conjg(ufh(xp,y,2))*b22(m)
		 if (istep.le.2) then
		    write (800+taskid,901) xp,y,ufh(xp,y,1),uny(y,zp,xp,1)
		    write (800+taskid,901) xp,y,ufh(xp,y,2),uny(y,zp,xp,2)
		 end if
 901		 format (2i4,1p,4e12.4)
		 sumh=sumh+.5*tfact(x)*(term1+term2)
 32	      continue	
 30	   continue
	end if
       
	if (xist.eq.1) then
	   do zp=1,zjsz
	      z=zp+zjst-1
	      kzabs=abs(kz(z))
	      if ((kzabs.gt.2.999.and.kzabs.lt.5.001)) then
		 term1=uny(1,zp,1,1)*conjg(ufv(zp,1))*b11(m)
		 term2=uny(1,zp,1,2)*conjg(ufv(zp,2))*b22(m)
		 sumv=sumv+.5*(term1+term2)
	      end if
	   end do
	end if
       
	call MPI_ALLREDUCE (sumh,sumh0,1,mpireal,MPI_SUM,
	1    MPI_COMM_WORLD,mpierr)
	call MPI_ALLREDUCE (sumv,sumv0,1,mpireal,MPI_SUM,
	1    MPI_COMM_WORLD,mpierr)

 	rate1=(sumh0+sumv0)
	if (taskid.eq.0) then
c	write (6,*) 'force_eirate: ',istep,sumh0,sumv0,rate1,dt
	end if
#endif       
	return
	end
