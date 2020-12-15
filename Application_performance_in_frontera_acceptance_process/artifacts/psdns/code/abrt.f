	subroutine abrt(msg)
	use mpicom
	character*(*) msg
	print *, 'ABORT : ', msg
	call MPI_ABORT(MPI_COMM_WORLD,ierr)
	return
	end
