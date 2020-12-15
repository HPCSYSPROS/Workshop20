        subroutine field_pen (uny)
!
c version using a non-cylindrical "pencils" structure
!
      use comp
	
	use ranseed
c
 	implicit none
	include 'intvars'

	complex(b8) :: uny(xisz,zjsz,ny,3+nc)

      integer mpistatus(MPI_STATUS_SIZE)
!
      real, allocatable :: k(:),k12(:),ef(:),phi(:)
      real, allocatable :: theta1(:),theta2(:),k1u(:)
      real, allocatable :: kd(:),phase(:,:),angphi(:)
      complex,allocatable :: alpha(:),beta(:)
      real, allocatable :: ekin(:),kden(:),ekinsc(:)
      integer, allocatable :: xfst(:)

      integer :: xst,i,ii,nnn,ir,ik,luran
      complex :: twopii,argc,ctmp1,ctmp2,ctmp
	real :: k2u,k3u,kyzsq,voljac,s1,s2,s3,twopi,rpid
	real :: s1d,s2d,s3d
	real :: dummy,ekx,frac,rtmp

	character*9 fn
	character*30 caux

	luran=1111
!
      allocate (k(nxh),k12(nxh),ef(nxh),phi(nxh))
      allocate (theta1(nxh),theta2(nxh),k1u(nxh))
      allocate (kd(nxh),phase(nxh,ncd),angphi(nxh))
      allocate (alpha(nxh),beta(nxh))
      allocate (ekin(nxhp),kden(nxhp))
      allocate (xfst(ny))

!
      if (taskid.eq.0) write (6,*) 'enter field'
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
        s1=sqrt(b11(1))
        s2=sqrt(b22(1))
        s3=sqrt(b33(1))
!
      s1d=1./s1
      s2d=1./s2
      s3d=1./s3
!
      do y=1,ny
      xfst(y)=1
      end do
!
! go through the loop if specified spectrum is to be used for
! any of the dependent variables
! otherwise skip all the action
!
	if (all(kinit.ne.-1).and.all(kinit.ne.-3)) go to 95
!
c next 5 lines moved here from inside the next IF loop, PK Yeung 12/6/08
      if (taskid.eq.0)  call ranseq (1,luran)
      call MPI_BCAST (iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!ccc
	iseed=iseed+taskid
!ccc
      if (any(kinit(1:3).lt.0)) then
!
      rseed(1)=iseed
      call random_seed (size=ir)
	write (6,*) 'field: after 1st random_seed call, taskid=',taskid,
     1	iseed
      call random_seed (put=rseed(1:2))
!	write (6,*) 'field: after 2nd random_seed call, taskid=',taskid
!
      if (kinit(1).eq.-3) then
      call fopen1 (3,'ekinp','formatted  ')
      read (3,501)
      nnn=nxh*beta1
      read (3,101) (ekin(ik),ik=1,nnn)
      close (3)

	write (caux,fmt="('den',i0)") nx
      call fopen1 (3,caux,'formatted  ')
!
      if (nx.le.128) then
      do ik=2,nxhp
      read (3,102) kden(ik)
      end do
      else
      do ik=2,nxhp
      read (3,1021) kden(ik)
      end do
      end if
!
	if (taskid.eq.0) 
     1	 write (6,101) (ekin(ii),ii=1,nxh)

      nnn=nx*beta1/2+1
      ekin(nnn)=0.


 201  format (1p,10e13.5)
!
 	do ik=2,kmax+1
 	ekin(ik)=alog10(ekin(ik))
 	end do
!
      end if
!101  format ((3x,10e13.5))
 101  format (10e13.5)
 102  format (28x,f6.4)
 1021 format (29x,f6.4)
 501  format (1x)
!
      call MPI_BCAST (iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!  loop over z and y
!
	do 100 zp=1,zjsz
!
      z=zjst+zp-1

      do 10 y=1,ny
      k2u=s2*ky(y)

      call masks (y)

      xfst(1)=1
      if (z.eq.1) xfst(1)=2
!
      k3u=s3*kz(z)
      kyzsq=k2u**2+k3u**2
      if (z.eq.nzhp) go to 1
!
        xst=1
        if (y.eq.1.and.z.eq.1.and.xist.eq.1) xst=2
!
      if (xst.eq.2) then
      dummy=ranu(1)
      dummy=ranu(1)
      dummy=ranu(1)
      end if
!
      if (kinit(1).eq.-1) then
!
        do 15 xp=xst,xisz
	  x=xist+xp-1
        k1u(x)=s1*kx(x)
!
        k(x)=sqrt(k1u(x)**2+kyzsq)
        kd(x)=1./k(x)
!
        ef(x)=0.0

	  ekx=e( k(x), kmax) * voljac

        ef(x)=mask(xp,zp)*sqrt( ekx ) * rpid*kd(x)
!
!
      if (k(x).gt.kmax.or.mask(xp,zp).eq.0.) go to 16
c
!
        theta1(x)=ranu(1)
        theta2(x)=ranu(1)
        phi(x)=ranu(1)
!
        k12(x)=sqrt(k1u(x)**2+k2u**2)
!
        angphi(x)=twopi*phi(x)
        argc=twopii*theta1(x)
        alpha(x)=ef(x)*cexp(argc)*cos(angphi(x))
        argc=twopii*theta2(x)
        beta(x) =ef(x)*cexp(argc)*sin(angphi(x))
!
!
! **caution**: here we are getting u1/s1, u2/s2 and u3/s3
! in wavenumber space (compare with rogallo 1981, p.53)
!
        if( k12(x) .eq. 0. ) then
                uny(xp,zp,y,1)=alpha(x)/s1
                uny(xp,zp,y,2)=beta(x)/s2
                uny(xp,zp,y,3)=0.
        else
                ctmp1=alpha(x)*k(x)
                ctmp2=beta(x)*k3u
                rtmp=kd(x)/k12(x)
                uny(xp,zp,y,1)=(ctmp1*k2u+ctmp2*k1u(x))*rtmp*s1d
                uny(xp,zp,y,2)=(ctmp2*k2u-ctmp1*k1u(x))*rtmp*s2d
                uny(xp,zp,y,3)=-beta(x)*k12(x)*kd(x)*s3d
        end if
!	  
        go to 15
!
 16     uny(xp,zp,y,1)=0.
        uny(xp,zp,y,2)=0.
        uny(xp,zp,y,3)=0.
!
 15     continue
!
c
      call masks (y)
      else if (kinit(1).eq.-3) then
!
        do 35 xp=xst,xisz
	  x=xist+xp-1
        k1u(x)=s1*kx(x)
!
        k(x)=sqrt(k1u(x)**2+kyzsq)
      if (kx(x).gt.kmax.or.mask(xp,zp).eq.0.) go to 36
!
        kd(x)=1./k(x)
!
        ef(x)=0.0
        ik=k(x)+1.5
        frac=k(x)+1-ik
        ekx=frac*ekin(ik+1)+(1.-frac)*ekin(ik)
!
        ekx=10.**ekx
!
        if (k(x).gt.kmax) ekx=0.
        ekx=ekx*voljac
!
        ef(x)=mask(xp,zp)*sqrt( ekx ) * rpid*kd(x)
!
        theta1(x)=ranu(1)
        theta2(x)=ranu(1)
        phi(x)=ranu(1)
!
      if (ef(x).le.0..or.mask(xp,zp).eq.0.) go to 36
!
        k12(x)=sqrt(k1u(x)**2+k2u**2)
!
         angphi(x)=twopi*phi(x)
        argc=twopii*theta1(x)
        alpha(x)=ef(x)*cexp(argc)*cos(angphi(x))
        argc=twopii*theta2(x)
        beta(x) =ef(x)*cexp(argc)*sin(angphi(x))
!
        if( k12(x) .eq. 0. ) then
                uny(xp,zp,y,1)=alpha(x)/s1
                uny(xp,zp,y,2)=beta(x)/s2
                uny(xp,zp,y,3)=0.
        else
                ctmp1=alpha(x)*k(x)
                ctmp2=beta(x)*k3u
                rtmp=kd(x)/k12(x)
                uny(xp,zp,y,1)=(ctmp1*k2u+ctmp2*k1u(x))*rtmp*s1d
                uny(xp,zp,y,2)=(ctmp2*k2u-ctmp1*k1u(x))*rtmp*s2d
                uny(xp,zp,y,3)=-beta(x)*k12(x)*kd(x)*s3d
        end if
!
        go to 35
!
 36     uny(xp,zp,y,1)=0.
        uny(xp,zp,y,2)=0.
        uny(xp,zp,y,3)=0.
!
 35     continue
      end if
!
 1      continue
 10     continue
 100    continue
#endif
!
      call masks (y)
        if (zjst.eq.1.and.xist.eq.1) then
        uny(1,1,1,1)=0.
        uny(1,1,1,2)=0.
        uny(1,1,1,3)=0.
	end if
!
      end if
!
 1000 continue
!
!
! check of conjugate symmetry
!
!     zg=taskid+1
!     lu=taskid+91
!     do zp=1,mz
!     z=(zg-1)*mz+zp
!     do y=1,ny
!     write (lu,900) y,z,ky(y),kz(z),un(1,2,y,zp)
!900  format ('y,z,ky,kz,un(1,2,y,zp)=',2i4,2f5.0,1p,2e12.4)
!     end do
!     end do
!
      end if
!
!  scalar fields
!  scalars having kinit=-1 are initialised here, with initial
!  spectra same as the 3-d energy spectrum
!  scalars having kinit=-2 are left undone, to be handled by
!  sub. inscf, called from sub. incond after return from this routine
!
#ifdef NOSCALAR
#else
        if (nc.le.0) go to 3000
!
	if (taskid.eq.0) write (6,*) 'field: kinit for scalars=',
     1         kinit(4:3+nc)
        if (any(kinit(4:3+nc).eq.-3)) then
          call fopen1 (3,'ekinp.sc','formatted  ')
	    if (taskid.eq.0) write (6,*) 'field: after opening ekinp.sc'
	  endif
	  allocate (ekinsc(nxhp))

c for scalars: division of \sqrt(2.) here is necessary so 
c that the integral of the specified spectrum is the variance, 
c not half of the mean-square
c
	if (taskid.eq.0) write (6,*) 'rpid:',rpid,rpid/sqrt(2.)
	rpid=rpid/sqrt(2.)
c
        do 2 i=1,nc
        j=i+3
        if (kinit(j).ne.-1 .and. kinit(j).ne.-3) go to 2
        write (6,601) i
 601    format (' scalar no. ',i1,'  initialized as Gaussian')
        call ranseq (kran(i),luran)
c
ccc PK Yeung, 12/6/08
	iseed=iseed+taskid
ccc

      if (kinit(j).eq.-3) then
	  ekinsc(:)=0.
        read (3,501)
        nnn=nxh*beta1
        !read (3,101) (ekinsc(ik),ik=1,nnn)
        read (3,*) (ekinsc(ik),ik=1,nnn)
	ekinsc(kmax+2:nnn)=0.
#ifdef TEST
 	  do ik=2,kmax+1
 	    ekinsc(ik)=alog10(ekinsc(ik))
 	  end do
#endif
	endif
!
      do 2000 zg=1,numtasks
      if (taskid.eq.zg-1) then
!
!     if (taskid.gt.0) call MPI_RECV (iseed,1,MPI_INTEGER,taskid-1,
!    1                 taskid-1,MPI_COMM_WORLD,mpistatus,mpierr)
!     if (taskid.gt.0) call MPI_RECV (rseed,2,MPI_INTEGER,taskid-1,
!    1                 taskid-1,MPI_COMM_WORLD,mpistatus,mpierr)
	  if (kinit(3+i).eq.-1) then

        do 200 zp=1,zjsz
	    z=zjst+zp-1
          k3u=s3*kz(z)
          if (z.eq.nzhp) go to 200
!
	    do 20 y=1,ny
            k2u=s2*ky(y)
            xfst(1)=1
            if (z.eq.1) xfst(1)=2
            call masks (y)
        xst=1
        if (y.eq.1.and.z.eq.1.and.xist.eq.1) xst=2
!
            do 25 xp=xst,xisz
	      x=xist+xp-1
            k1u(x)=s1*kx(x)
            k(x)=sqrt(k1u(x)**2+k2u**2+k3u**2)
      if (k(x).gt.kmax.or.mask(xp,zp).eq.0.) go to 25
            kd(x)=1./k(x)
            ef(x)=0.0

            ekx=e( k(x), kmax ) *voljac
!           if( k(x) .ne. 0.0 ) ef(x)=mask(xp,zp)*sqrt( ekx /(twopi) )*kd(x)
            if( k(x) .ne. 0.0 ) ef(x)=mask(xp,zp)*sqrt(ekx)*rpid*kd(x)
!
            phase(x,i)=ranu(kran(i))
            ctmp=ef(x)*cexp(twopii*phase(x,i))
 	      uny(xp,zp,y,i+3)=ctmp
 25         continue
 20       continue
 200	  continue

	  elseif (kinit(3+i).eq.-3) then

        do 2001 zp=1,zjsz
	    z=zjst+zp-1
          k3u=s3*kz(z)
          if (z.eq.nzhp) go to 2001
!
	    do 21 y=1,ny
            k2u=s2*ky(y)
            xfst(1)=1
            if (z.eq.1) xfst(1)=2
            call masks (y)
        xst=1
        if (y.eq.1.and.z.eq.1.and.xist.eq.1) xst=2

!
            do 26 xp=xst,xisz
	      x=xist+xp-1
            k1u(x)=s1*kx(x)
            k(x)=sqrt(k1u(x)**2+k2u**2+k3u**2)
      if (kx(x).gt.kmax.or.mask(xp,zp).eq.0.) go to 26
            kd(x)=1./k(x)
            ef(x)=0.0

            ik=k(x)+1.5
            frac=k(x)+1-ik
            ekx=frac*ekinsc(ik+1)+(1.-frac)*ekinsc(ik)
c		if (ik.lt.nx/4) ekx=0.
		ekx=abs(ekx)
#ifdef TEST
            ekx=10.**ekx
#endif
            if (k(x).gt.kmax) ekx=0.
            ekx=ekx*voljac
	if (ik.eq.4) write (6,901) x,y,z,k(x),frac,ekinsc(ik),ekinsc(ik+1)
 901	format ('x,y,z,k,frac,ekinsc=',3i4,1p,4e12.4)
	if (ik.eq.4) write (6,902) x,y,z,ekx,rpid,kd(x)
 902	format ('x,y,z,ekx,rpid,kd=',3i4,1p,4e12.4)

            ef(x)=mask(xp,zp)*sqrt( ekx ) * rpid*kd(x)
!
            phase(x,i)=ranu(kran(i))
            ctmp=ef(x)*cexp(twopii*phase(x,i))
 	      uny(xp,zp,y,i+3)=ctmp
 26         continue
 21       continue
 2001	continue
!
!     if (taskid.lt.nzg-1) call MPI_SSEND (iseed,1,MPI_INTEGER,
!    1                     taskid+1,taskid,MPI_COMM_WORLD,mpierr)
!     if (taskid.lt.nzg-1) call MPI_SSEND (rseed,2,MPI_INTEGER,
!    1                     taskid+1,taskid,MPI_COMM_WORLD,mpierr)
!
	endif ! if kinit=-1 else kinit=-3

	
      end if ! if taskid.eq.zg-1
 2000   continue
!
!
 2      continue
	  deallocate (ekinsc)
        if (any(kinit(4:3+nc).eq.-3)) then
		close(3)
	  endif
	  close(1000+taskid)
!
#endif
 3000 continue
!
! pass the current random no. seed back to task 0
!
!     if (taskid.eq.nzg-1) call MPI_SSEND (iseed,1,MPI_INTEGER,
!    1                     0,nzg-1,MPI_COMM_WORLD,mpierr)
!     if (taskid.eq.0) call MPI_RECV (iseed,1,MPI_INTEGER,nzg-1,nzg-1,
!    1                                MPI_COMM_WORLD,mpistatus,mpierr)
      call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
!

! enforce conjugate symmetry
!
      if (taskid.eq.0) write (6,*) 'before call enfsym'
      call enfsym(uny)
!	print *,'field: max = ',maxval(abs(uny)),maxloc(abs(uny))

!	write (fn,fmt='(a5,i4))') 'fort.',8000+taskid
!	open(8000+taskid,file=fn,form='unformatted')
!	write (8000+taskid) xist,zjst,xisz,zjsz
!	do i=1,3
!	write (8000+taskid) (((xist+xp-1,y,zjst+zp-1,uny(xp,zp,y,i),
!    1	xp=1,xisz),y=1,ny),zp=1,zjsz)
!	enddo
!	close(8000+taskid)
! 	do i=1,3
! 	do 2222 zp=1,zjsz
!     z=zjst+zp-1
! 	do 2222 y=1,ny
!	call masks(y)
!	do 2222 xp=1,xisz
!	  x=xist+xp-1
!	  write (9000+taskid,*) x,y,z,i,uny(xp,zp,y,i)
!2222 continue	
! 	enddo
!	close(9000+taskid)
!
      write (6,*) ' exit field, taskid,iseed=',taskid,iseed
!
 95   continue
!
      deallocate (k,k12,ef,phi)
      deallocate (theta1,theta2,k1u)
      deallocate (kd,phase,angphi)
      deallocate (alpha,beta)
      deallocate (ekin,kden)
      deallocate (xfst)
      if (taskid.eq.0) write (6,*) ' exit field'

        return
        end
