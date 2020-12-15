      subroutine bsminmax (string,bsxy,nv)

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
	integer iv,nv
c
	real(p8) bmin,bmax,bmin0,bmax0

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
	do zp=1,bzjsize
	do y=1,nby
	do xp=1,bxisize
	bmin=min(bmin,bsxy(xp,y,zp,nv))
	bmax=max(bmax,bsxy(xp,y,zp,nv))
	end do
	end do
	end do
c
        call MPI_ALLREDUCE (bmin,bmin0,1,pptype,MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (bmax,bmax0,1,pptype,MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
c
        if (taskid.eq.0) then
        write (6,*) 'bsminmax:',string,' overall bmin,bmax=',bmin0,bmax0
        write (6,*) 'bsminmax:',string,' bmin,bmax,bsxy=',bmin,bmax,bsxy(4,4,4,1)
	end if
c
      return


#endif
      end

