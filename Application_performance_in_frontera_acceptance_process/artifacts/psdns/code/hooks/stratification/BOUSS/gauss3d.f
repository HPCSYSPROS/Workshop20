        subroutine gauss3d (uny,m)
!
c
c routine to generate 3D Gaussian random fields, similar to field.f in DNS code
c
!
      use comp
 	implicit none
	include 'intvars'

	complex(b8) :: uny(ny,zjsz*xisz,3+nc)
      common/rancom/iseed,rseed(2)
	integer iseed
	integer rseed
	real eno

	integer a
c
      integer mpistatus(MPI_STATUS_SIZE)
!
c
      real, allocatable :: k1u(:)
      real :: k,k12,ef,phi,theta1,theta2,kd,phase,angphi
	complex :: alpha,beta
	integer :: yst


      integer :: xst,i,ii,nnn,ir,luran,m,j
      complex :: twopii,argc,ctmp1,ctmp2,ctmp
	real :: k2u,k3u,kyzsq,voljac,s1,s2,s3,twopi,rpid
	real :: s1d,s2d,s3d
	real :: dummy,ekx,frac,rtmp
c
	real ranu
	real e1,e2,e3,emode,tfact

	character*9 fn
	character*30 caux

	integer, allocatable :: ik(:,:)
	allocate (ik(ny,xisz*zjsz))
c
	luran=1111
c
	allocate (k1u(nxh))
!

!
      if (taskid.eq.0) write (6,*) 'enter gauss3d'
! form jacobian
!
      voljac=beta1*beta2*beta3
!
!       form unnormalised wave number vectors ( k1u , k2u , k3u )
!
        twopi=2.*pi
        twopii=twopi*imagi
        rpid=1./sqrt(twopi)
!
        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
!
      s1d=1./s1
      s2d=1./s2
      s3d=1./s3
!
!
! go through the loop if specified spectrum is to be used for
! any of the dependent variables
! otherwise skip all the action
!
c next 5 lines moved here from inside the next IF loop, PK Yeung 12/6/08
c
c	if (taskid.eq.0) iseed=100000
	if (taskid.eq.0) then
	iseed=200000
	open (1112,file='seed') 
	read (1112,*,end=2) iseed
 2	continue
	end if
c
      call MPI_BCAST (iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	iseed=iseed+taskid
      rseed(1)=iseed
      call random_seed (size=ir)
      call random_seed (put=rseed(1:2))
!
      xp = mystart_x
      z = mystart_z

        do 100 a=1,num_al
           x=xp+xist-1
           if (z.eq.nzhp) go to 100
c
        k1u(x)=s1*kx(x)
!
      k3u=s3*kz(z)
!
         if(x .eq. 1 .and. z .eq. 1) then
            yst = 2
         else
            yst = 1
         endif


      do 10 y=yst,ny
      if (.not.mask(y,a)) go to 10
c
      k2u=s2*ky(y)
      kyzsq=k2u**2+k3u**2

        k=sqrt(k1u(x)**2+kyzsq)
!
        kd=1./k
!
        ef=0.0
!
  	ekx = eno(k)
c
        ekx=ekx*voljac
        if (k.gt.kmax) ekx=0.
c
!
        ef=sqrt( ekx ) * rpid*kd
c
        theta1=ranu(1)
        theta2=ranu(1)
        phi=ranu(1)
c
#ifdef TEST (for comparisons without random number generation)
	theta1=float(x)/(nx+1)
	theta2=float(y)/(ny+1)
	phi=float(z)/(nz+1)
#endif
c
        ik(y,a)=k/beta1+1.5

!
        k12=sqrt(k1u(x)**2+k2u**2)
!
         angphi=twopi*phi
        argc=twopii*theta1
        alpha=ef*cexp(argc)*cos(angphi)
        argc=twopii*theta2
        beta =ef*cexp(argc)*sin(angphi)
!
        if( k12 .eq. 0. ) then
                uny(y,a,1)=alpha/s1
                uny(y,a,2)=beta/s2
                uny(y,a,3)=0.
        else
                ctmp1=alpha*k
                ctmp2=beta*k3u
                rtmp=kd/k12
                uny(y,a,1)=(ctmp1*k2u+ctmp2*k1u(x))*rtmp*s1d
                uny(y,a,2)=(ctmp2*k2u-ctmp1*k1u(x))*rtmp*s2d
                uny(y,a,3)=-beta*k12*kd*s3d
        end if
c
 10     continue
c
	call next_xz (xp,z)
!
 100	continue
c
        if (taskid.eq.0) then
        uny(1,1,1)=0.
        uny(1,1,2)=0.
        uny(1,1,3)=0.
	end if
!
!
      call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
!

! enforce conjugate symmetry
!
      call enfsym(uny)

!
 95   continue
!
      deallocate (k1u)

        return
        end
c
	function eno(k)
	real k
	eno=0.
 	eno = .17*k**4/(1.+k**2)**4
	return
	end

