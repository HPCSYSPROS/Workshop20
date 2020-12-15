	subroutine gltime2 (tcpu,tav,tmin,tmax)
c
c ASSUME USE OF DOUBLE-PRECISION FOR TIMERS
c
	use mpicom
	implicit none
c
	real(8), allocatable :: all(:)
	real(8) tcpu,tav,tmin,tmax
	integer itask
c
	allocate (all(numtasks),stat=ierr)

	call MPI_GATHER (tcpu,1,MPI_DOUBLE_PRECISION,
     1                    all,1,MPI_DOUBLE_PRECISION,
     1		          0,MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) then
c
	tav=0.
	tmin=10000.
	tmax=0.0001
c
	do itask=1,numtasks
	tav=tav+all(itask)
	tmin=dmin1(tmin,all(itask))
	tmax=dmax1(tmax,all(itask))
	end do
c
	tav=tav/numtasks
c
	end if
c
	deallocate (all)
c
	return
	end
