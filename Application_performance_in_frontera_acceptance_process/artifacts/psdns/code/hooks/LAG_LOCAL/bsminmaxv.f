      subroutine bsminmaxv (string,bsxy,nv)

#ifdef LAG
c version for 2D code
c
c
c Fix by D. Buaria 9/1/2012: dimensioning of 'work' array
c is changed to work with cases where nprec is not 
c integer multiple of num_thr
c
        use com, only:istep,nsteps
        use compart, only:nop,nop2,nom,iout
        use mpilag
        use lag_timers
        implicit none
	character*(*) string
c
c
	integer iv,nv,itask
c
	real(p8) bmin,bmax,bmin0,bmax0
        real(p8) bmin_all(0:numtasks-1),bmax_all(0:numtasks-1)

#if defined(CAF) && defined(CF_LAG)
      real(p8) bsxy(xiind,nby,zjind,nv)[0:*]
#else
      real(p8) bsxy(xiind,nby,zjind,nv)
#endif
c
c
      integer xp,zp,y
c
c for function or 1st derivative of spline, just pick the
c correct basis functions to use in the summation, based
c on the order of differentiation
c
	return
c
	bmin=1.e10
	bmax=-1.e10
	do iv=1,nv
	do zp=1,bzjsize
	do y=1,nby
	do xp=1,bxisize
	bmin=min(bmin,bsxy(xp,y,zp,iv))
	bmax=max(bmax,bsxy(xp,y,zp,iv))
	end do
	enddo
	end do
	end do
c
        call MPI_ALLREDUCE (bmin,bmin0,1,pptype,MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (bmax,bmax0,1,pptype,MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
c
        call MPI_GATHER (bmin,1,pptype,bmin_all,1,pptype,0,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_GATHER (bmax,1,pptype,bmax_all,1,pptype,0,
     1                   MPI_COMM_WORLD,mpierr)
c
        if (taskid.eq.0) then
        write (6,*) 'bsminmaxv:',string,' overall bmin,bmax=',bmin0,bmax0
	if (abs(bmin0).gt.20..or.abs(bmax0).gt.20.) then
	do itask=0,numtasks-1
	write (6,*) 'bsminmaxv: taskid:',itask,bmin_all(itask),bmax_all(itask)
	end do
	end if
	end if
c
	if (abs(bmin0).gt.20..or.abs(bmax0).gt.20.) then
	if (taskid.eq.0) write (6,*) 'abort in bsminmaxv'
	call MPI_FINALIZE (mpierr)
	stop
	end if
c
      return


#endif
      end

