	subroutine shifts
c
	use com
	include 'intvars'

c
c kshift=1: random shifts, in which the aliasing error on the
c             predictor step approximately cancels that of the previous
c             corrector step (rogallo).  over a fixed time period,
c             the residual aliasing error in u is of order dt**2,
c             but over a single step it is of order dt.  not suitable
c             for extraction of multi-time statistics.
c
c kshift=2: random shifts, in which the aliasing error on the
c             corrector step cancels that of the previous
c             predictor step.  over a fixed time period,
c             the residual aliasing error in u is of order dt,
c             but over a single step it is of order dt**2.
c
c kshift=3: fixed shifts.  error on any two adjacent steps cancels.
c             over a fixed time period, and over a single step, the
c             residual aliasing error in u is of order dt**2, but
c             but the residual error is not random.
c
c kshift=4: zero shifts: order one aliasing errors.
c           (mainly for frozen fields)
c
c the index 1 in gsh and sx is for the predictor step.
c the index 2 in gsh and sx is for the corrector step.
c
      complex :: phase,argc
c
c       common/rancom/iseed,rseed(2)
c       real rseed
c
        parameter ( nk=20 )
        common/seqcom/iseq(nk),klast
        dimension iseqdf(nk)
c
      cancel( delta ) = amod( delta + 0.5 , 1. )
      if (taskid.ne.0) go to  111
c
c  loop over directions
c
c initialise gsh(i,2) if needed
c
      if (istart.le.0.and.istep.eq.istep0+1.and.
     1   (kshift.eq.1.or.kshift.eq.3)) then
      call ranseq (ksran,0)
      do 5 i=1,3
 5    gsh(i,2)=ranu(ksran)
#ifdef RKFOUR
#ifndef MAY28
	do i=1,3
	gsh(i,4)=ranu(ksran)
	end do
#endif
#endif
      end if
c
c generate new random numbers for grid shifts if desired
c
      if (kshift.eq.1.or.kshift.eq.2) call ranseq (ksran,0)
c
	do 100 i=1,3

      if( kshift .eq. 1 ) then
         gsh(i,1) = cancel( gsh(i,2) )
         gsh(i,2) = ranu(ksran)
 
      else if( kshift .eq. 2 ) then
         gsh(i,1) = ranu(ksran)
         gsh(i,2) = cancel( gsh(i,1) )
 
      else if( kshift .eq. 3 ) then
         gsh(i,1) = cancel( gsh(i,2) )
         gsh(i,2) = cancel( gsh(i,1) )
 
      else
         gsh(i,1) = 0.
         gsh(i,2) = 0.
 
      endif
 100  continue
c
#ifdef RKFOUR
#ifndef MAY28
	do 300 i=1,3

      if( kshift .eq. 1 ) then
         gsh(i,3) = cancel( gsh(i,4) )
         gsh(i,4) = ranu(ksran)
 
      else if( kshift .eq. 2 ) then
         gsh(i,3) = ranu(ksran)
         gsh(i,4) = cancel( gsh(i,3) )
 
      else if( kshift .eq. 3 ) then
         gsh(i,3) = cancel( gsh(i,4) )
         gsh(i,4) = cancel( gsh(i,3) )
 
      else
         gsh(i,3) = 0.
         gsh(i,4) = 0.
 
      endif
 300  continue
#endif
#endif

#ifdef RKFOUR
#ifdef JUN1411
	do i=1,3
	gsh(i,1)=0.
	gsh(i,2)=1./3.
	gsh(i,3)=2./3.
	gsh(i,4)=0.
	end do
#endif
#endif


 111  continue


#ifdef RKFOUR
c     call MPI_BCAST (gsh,3*rkmethod,mpireal,0,MPI_COMM_WORLD,mpierr)
      call MPI_BCAST (gsh,3*4,mpireal,0,MPI_COMM_WORLD,mpierr)
#else
      call MPI_BCAST (gsh,3*2,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif
c
#ifdef RKFOUR
#ifndef JUN1411
 	gsh(:,3)=gsh(:,1)
 	gsh(:,4)=gsh(:,2)
#endif
#endif
c

c compute phase shift factors from grid shifts
c
c x shifts
!	write (6,*) 'shifts: rkmethod=',rkmethod
c
        phase=imagi*2*pi/nx
#ifdef RKFOUR
        do 1 k=1,rkmethod
#else
        do 1 k=1,2
#endif
        do 1 x=1,nxh
        argc=phase*kx(x)*gsh(1,k)
 1      sx(x,k) = cexp( argc )

c y shifts
c
        phase=imagi*2*pi/ny
#ifdef RKFOUR
        do 2 k=1,rkmethod
#else
        do 2 k=1,2
#endif
        do 2 y=1,ny
        argc=phase*ky(y)*gsh(2,k)
 2      sy(y,k) = cexp( argc )

c z shifts
c
        phase=imagi*2*pi/nz
#ifdef RKFOUR
        do 3 k=1,rkmethod
#else
        do 3 k=1,2
#endif
        do 3 z=1,nz
        argc=phase*kz(z)*gsh(3,k)
 3      sz(z,k) = cexp( argc )
c
	return
	end
