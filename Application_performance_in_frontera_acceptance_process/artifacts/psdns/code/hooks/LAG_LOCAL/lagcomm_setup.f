
        subroutine lagcomm_setup

#ifdef LAG
#ifdef CF_LAG

        use mpilag
	use compart, only : nsubset, nsubsize

	implicit none
	integer icolor, ikey

!	nsubset = jproc
	nsubsize = numtasks/nsubset

        icolor = taskid/nsubsize
        ikey = mod(taskid,nsubsize)

        call MPI_COMM_SPLIT (MPI_COMM_WORLD,icolor,ikey,mpi_comm_lagout,mpierr)
        call MPI_COMM_RANK (mpi_comm_lagout, lagc_id, mpierr)

!        write(6,*) 'taskid,icolor,ikey,lagc_id=',taskid,icolor,ikey,lagc_id
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)


#endif
#endif

        return

        end

