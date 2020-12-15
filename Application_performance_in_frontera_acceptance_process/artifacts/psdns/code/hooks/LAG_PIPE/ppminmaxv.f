      subroutine ppminmaxv (string,xp,npdim)
c
#ifdef LAG

c This is a diagnostic routine which can be called in partin.f
c and/or other places to see if the particle coordinates have
c somehow gone far out of reasonable range

	use compart
	implicit none
c
	integer i,j,itask,npdim
	real(p8) pmin,pmax,pmin0,pmax0
	character*(*) string
	real(p8) xp(npdim/numtasks,3)
	real(p8) pmin_all(0:numtasks-1),pmax_all(0:numtasks-1)
c
 	return

c
	pmin=1.e12
	pmax=-1.e12
	do j=1,3
	do i=1,npdim/numtasks
	pmin=min(pmin,xp(i,j))
	pmax=max(pmax,xp(i,j))
	end do
	end do
c
	call MPI_ALLREDUCE (pmin,pmin0,1,pptype,MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)
	call MPI_ALLREDUCE (pmax,pmax0,1,pptype,MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
c	write (6,*) 'taskid,pmin,pmax=',taskid,pmin,pmax
c
	if (taskid.eq.0) then
	write (6,*) 'ppminmaxv:',string,' overall pmin,pmax,npdim=',pmin0,pmax0,npdim
	do i=1,npdim/numtasks
	write (20000,"(i6,1p,3e12.4)") i,xp(i,1),xp(i,2),xp(i,3)
	end do
	end if
	if (taskid.eq.numtasks-1) then
	do i=1,npdim/numtasks
	write (20003,"(i6,1p,3e12.4)") i,xp(i,1),xp(i,2),xp(i,3)
	end do
	end if

	call MPI_GATHER (pmin,1,pptype,pmin_all,1,pptype,0,
     1                   MPI_COMM_WORLD,mpierr)
	call MPI_GATHER (pmax,1,pptype,pmax_all,1,pptype,0,
     1                   MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) then
	pmin0=1.e12
	pmax0=-1.e12
	do itask=0,numtasks-1
	pmin0=min(pmin0,pmin_all(itask))
	pmax0=max(pmax0,pmax_all(itask))
	write (20001,*) 'itask,pmax=',itask,pmax_all(itask),pmax0
	write (20002,*) 'itask,pmin=',itask,pmin_all(itask),pmin0
	end do
	write (20001,*) 'pmin0,pmax0=',pmin0,pmax0
	end if
c
	if (pmin0.lt.-10.or.pmax0.gt.10.) then
	call MPI_FINALIZE (MPI_COMM_WORLD,mpierr)
	stop 'velocities too large: stop in ppinmaxv'
	end if 
#endif
      end
