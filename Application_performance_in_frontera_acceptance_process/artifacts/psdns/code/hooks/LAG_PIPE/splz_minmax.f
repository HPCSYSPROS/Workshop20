      subroutine splz_minmax (string,bsz,nv)

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
	integer x,yp,zp
	integer mx,my,mz
c
	real(p8) bmin,bmax,bmin0,bmax0

      real(p8) bsz(nv,bxisize,nbz,yjsz)
c
c
	integer xp,z
c
c for function or 1st derivative of spline, just pick the
c correct basis functions to use in the summation, based
c on the order of differentiation
c
c	return
c
	bmin=1.e10
	bmax=-1.e10
	do yp=1,yjsz
	do z=1,nbz
	do xp=1,bxisize
	bmin=min(bmin,bsz(nv,xp,z,yp))
	if (bsz(nv,xp,z,yp).gt.bmax) then
	bmax=bsz(nv,xp,z,yp)
	mx=xp
	my=yp
	mz=z
	end if
	end do
	end do
	end do
c
c	if (taskid.eq.2) write (6,*) 'splz_minmax, task 2:',string,bmin,bmax,mx,mz,my
c
      return


#endif
      end

