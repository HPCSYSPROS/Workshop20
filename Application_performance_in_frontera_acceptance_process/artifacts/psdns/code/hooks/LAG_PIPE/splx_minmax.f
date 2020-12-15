      subroutine splx_minmax (string,bsx,nv)

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
c
	real(p8) bmin,bmax,bmin0,bmax0

      real(p8) bsx(nbx,zisz,yjsz,nv)
c
c
c for function or 1st derivative of spline, just pick the
c correct basis functions to use in the summation, based
c on the order of differentiation
c
 	return
c
	bmin=1.e10
	bmax=-1.e10
	do yp=1,yjsz
	do zp=1,zisz
	do x=1,nbx
	bmin=min(bmin,bsx(x,zp,yp,nv))
	bmax=max(bmax,bsx(x,zp,yp,nv))
	end do
	end do
	end do
c
	if (taskid.eq.2) write (6,*) 'splx_minmax, task 2:',bmin,bmax
c
      return


#endif
      end

