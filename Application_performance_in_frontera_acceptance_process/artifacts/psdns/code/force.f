        subroutine force
c
#ifdef EPFOR
c
c MPL version
c
c  routine to determine velocity increments due to forcing.
c
c  all wavenumber nodes within a sphere of radius kforce are forced.
c  the parameter kfor should be set to at least the integer larger
c    than kforce+1.
c
c  the forcing acceleration is an integrated
c       uhlenbeck-ornstein process (which is continuously
c       differentiable).  the degenerate cases, in which
c       (through the choice of the forcing time scales)
c       the forcing acceleration is a uo process or
c       white noise, are treated correctly.
c
c  forcing parameters:
c       tforce  - integral time scale of forcing acceleration
c                (white noise is obtained by setting tforce=0. )
c       tfiuo   - normalized taylor time scale of forcing
c                 ( tfiuo must be strictly less than 0.5.
c                 uo forcing is obtained by setting tfiuo=0. )
c       epsfor  - the variance of the forcing acceleration is
c                 epsfor/tforce.
c
c forcing is skipped if kforce is 0.
c
c 3/3/04:  If kfor is increased (without changing kforce, etc)
c the results will change slightly because a different number
c of random number samples are drawn into the velinc array
c
c
	use comp
	use force_mod
 	implicit none
	include 'intvars'
c
	integer i,nword,icall,ifwave,iyz,nyz
	complex :: alpha,beta
!     complex velinc(2,kfor,k2fo,k2fo)
!     save velinc
      real w1(nrproc),w2(nrproc)
      real k,k12,k1u,k2u,k3u
	real s1,s2,s3
	data icall/0/
	save icall

	integer i1,i2,i3
        character*15 caux
c
c---------------------------------------------------------------------
c
c  the statement function ifwave returns the integer wavenumber
c    corresponding to the index iyz.
c
c  if  iyz .le.  kfor ,  ifwave = iyz
c  if  iyz .gt.  kfor ,  ifwave = iyz + nyz - k2fo
c
        ifwave( iyz , nyz ) = iyz + ( iyz/(kfor+1) ) * ( nyz-k2fo )
c
       if (kforce.eq.0.) return
       if (taskid.ne.0) go to 111
c
c--------  no forcing       uo forcing; or, iuo forcing  --------------------
c
        if( kforce .eq. 0. ) then
           return
        elseif( tforce * tfiuo .eq. 0. ) then
           call vfuo( tforce, epsfor, dt, nrproc, velinc,
     1          suo, w1, w2, kranf1, kranf2, istart )
        else
           call vfiuo( tforce, tfiuo, epsfor, dt, nrproc, velinc,
     1          suo, svo, w1, w2, kranf1, kranf2, istart )
        endif
c
c  -------------------------------------------------------------------
c
        s1=sqrt(b11(1))
        s2=sqrt(b22(1))
        s3=sqrt(b33(1))
c
      icall=icall+1
        caux='forced_modes'
	if (icall.eq.1) call fopen1 (91,caux,'formatted  ')
c
c  loop over y and z  --------------------------------------------------
c
        do 10 z=1,k2fo
        do 10 y=1,k2fo
c
c    form unnormalised wave number vectors
c
        k2u=s2 * ky( ifwave(y,ny) )
        k3u=s3 * kz( ifwave(z,nz) )
c
        do 10 x=1,kfor
c
c  alpha and beta are components of velocity increment due to
c   forcing over the time step.  alpha is the component normal
c   to k and kz.  beta is the component normal to alpha and k.
c
       alpha = velinc(1,x,y,z)
       beta  = velinc(2,x,y,z)
c
c form unnormalised wave number vectors  k1u
c
        k1u=s1*kx(x)
        k=sqrt(k1u**2+k2u**2+k3u**2)
        k12=sqrt(k1u**2+k2u**2)
c
        if( k .gt. kforce  .or.  k .eq. 0.) then
               for(x,y,z,1)=cmplx(0.,0.)
               for(x,y,z,2)=cmplx(0.,0.)
               for(x,y,z,3)=cmplx(0.,0.)
c
        elseif( k12 .eq. 0. ) then
               for(x,y,z,1)=alpha/s1
               for(x,y,z,2)=beta/s2
               for(x,y,z,3)=0.0
	if (icall.eq.1)
     1write (91,901) x,y,z,k1u,k2u,k3u,k,kforce
 901  format ('x,y,z,k,kforce=',3i4,3f6.1,f10.6,f6.2)
!    1write (91,*) 'x,y,z,k,kforce=',x,y,z,k1u,k2u,k3u,k,kforce
        else
               for(x,y,z,1)=(alpha*k*k2u+beta*k1u*k3u)/(k*k12*s1)
               for(x,y,z,2)=(beta*k2u*k3u-alpha*k*k1u)/(k*k12*s2)
               for(x,y,z,3)=-beta*k12/(k*s3)
	if (icall.eq.1)
     1write (91,901) x,y,z,k1u,k2u,k3u,k,kforce
        end if
c
10       continue
c
	if (icall.eq.1) close (91)
c
c   -------------------------------------------------------------------
c
c        enforce conjugate symmetry on the random amplitudes
c
c  loop over velocities
c
        do 50 i=1,3
c
c kz=0
c
        do 55 y=kfor+1,k2fo
55      for(1,y,1,i)=conjg(for(1,k2fo+2-y,1,i))
c
        do 50 z=kfor+1,k2fo
c
c  ky=0
c
        for(1,1,z,i)=conjg(for(1,1,k2fo+2-z,i))
c
c ky>0 , kz>0
c
           do 50 y=2,k2fo
50         for(1,y,z,i)=conjg(for(1,k2fo+2-y,k2fo+2-z,i))
c
 111    continue
	nword=kfor*3*k2fo*k2fo

	call MPI_BCAST (for,nword,mpicomplex,0,MPI_COMM_WORLD,ierr)
c
c
#endif

        return
        end
