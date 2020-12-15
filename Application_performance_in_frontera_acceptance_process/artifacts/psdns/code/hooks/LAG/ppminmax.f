      subroutine ppminmax (string,xp,npdim)
c
#ifdef LAG

c This is a diagnostic routine which can be called in partin.f
c and/or other places to see if the particle coordinates have
c somehow gone far out of reasonable range

	use compart
	implicit none
c
	integer i,j,npdim
	real(p8) pmin,pmax,pmin0,pmax0
	character*(*) string
	real(p8) xp(npdim/numtasks,3)

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
c
	if (taskid.eq.0) then
	write (6,*) 'ppminmax:',string,' overall pmin,pmax,npdim=',pmin0,pmax0,npdim
	end if

	if (pmin0.lt.-100..or.pmax0.gt.nx*10.) then
	call MPI_FINALIZE (mpierr)
	stop 'ppminmax: positions out of range'
	end if
#endif
      end
