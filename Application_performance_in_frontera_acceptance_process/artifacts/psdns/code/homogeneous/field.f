! Routine to initialize velocity and/or scalar fields in wavenumber space
!
! Written in pure MPI
! Each MPI task starts with a different random number seed
! Cylindrical truncation, stride 1
! Spectrum may be set internally or read in from a file
! Conjugate symmetry must be enforced by calling 'enfsym' or equiv.
!
! Major cleanup and extensions by PKY, 7/23/2012
!     
      subroutine field(uny)      
      use comp
#ifdef RANDOM_SEED
	use ranseed
#endif
      implicit none
      include 'intvars'

      complex(b8) :: uny(ny,zjsz*xisz,3+nc)
      complex :: twopii,argc,ctmp1,ctmp2,ctmp,alpha,beta

#ifndef RANDOM_SEED
      common/rancom/iseed,rseed(2)
      integer iseed
      integer rseed
#endif

      integer mpistatus(MPI_STATUS_SIZE)
      integer :: yst,i,ii,ir,ik,luran,a,j

      real ranu, e !npm adding functions here
      real, allocatable :: k1u(:)
      real :: k,k12,ef,phi,theta1,theta2,kd,phase,angphi
      real :: k2u,k3u,kyzsq,voljac,s1,s2,s3,twopi,rpid
      real :: s1d,s2d,s3d
      real :: ekx,frac,rtmp
      real, allocatable :: ekin(:),kden(:),ekinsc(:)

      character*9 caux
	character*8 fnden
	integer nchar
	real(b8) deltak

      luran=1111

      if (taskid.eq.0) write (6,*) 'enter field'

      allocate(k1u(nxh))

!     form jacobian
!     
      voljac=beta1*beta2*beta3
!     
!     form unnormalised wave number vectors ( k1u , k2u , k3u )
!     
      twopi=2.*pi
      twopii=twopi*imagi
      rpid=1./sqrt(twopi)
      
      s1=sqrt(b11(1))
      s2=sqrt(b22(1))
      s3=sqrt(b33(1))

      s1d=1./s1
      s2d=1./s2
      s3d=1./s3

!      do y=1,ny
!         xfst(y)=1
!      end do
!     
!     go through the loop if specified spectrum is to be used for
!     any of the dependent variables
!     otherwise skip all the action
!     
      if (all(kinit.ne.-1).and.all(kinit.ne.-3)) go to 95      

      if (any(kinit.lt.0)) then         
         if (taskid.eq.0)  call ranseq (1,luran)
         call MPI_BCAST (iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         
         iseed=iseed+taskid
        if (taskid.eq.0) write (6,*) 'after broadcast iseed',iseed
         
c        rseed(1)=iseed
c Change by PKY, 1/18/201$ to make sure that user-suppied
c random number seed is fully utilized by the 
c intrinsic Fortran random_number generator (if the RANDOM_SEED
c directive is active)
#ifdef RANDOM_SEED
        rseed(1)=iseed
        rseed(2)=iseed+12345
         rseed(2)=iseed+123456
      call random_seed (put=rseed(1:seed_size))
#else
         call random_seed (size=ir)
         call random_seed (put=rseed(1:2))
#endif

c
c read node-density correction file (if needed)
c
	if (any(kinit.eq.-3)) then
	deltak=min(beta1,beta2,beta3)
	allocate (kden(2:nxhp))
	if (taskid.eq.0) then
	write (fnden,"('den',i5)") nx
	call blanks (fnden,nchar)
	call fopen1 (3,fnden(1:nchar),'formatted  ')
	write (6,*) 'field: nxhp=',nxhp
            if (nx.le.128) then
               do ik=2,nxhp
                  read (3,"(28x,f6.4)") kden(ik)
               end do
            else
               do ik=2,nxhp
                  read (3,"(29x,f6.4)") kden(ik)
               end do
	     end if
	close (3)
	write (6,*) 'field: after reading kden'
	end if
	call MPI_BCAST (kden,nxh,mpireal,0,MPI_COMM_WORLD,ierr)
	end if
c
c read input energy spectrum (if indicated)
c
         if (kinit(1).eq.-3) then
       allocate (ekin(nxhp))
	if (taskid.eq.0) then
	caux='ekinp'
            call fopen1 (3,caux,'formatted  ')
            read (3,501)
c           read (3,101) (ekin(ik),ik=1,nxh)
            read (3,*) (ekin(ik),ik=1,nxh)
	write (6,*) 'ekin(2)=',ekin(2)
		write (6,*) 'sum(ekin)=',sum(ekin)
            close (3)
	end if
	call MPI_BCAST (ekin,nxh,mpireal,0,MPI_COMM_WORLD,ierr)

c       if (taskid.eq.0) write (6,"((3x,1p,10e13.5))") (ekin(ii),ii=1,nxh)
	do ik=2,nxh
	ekin(ik)=ekin(ik)*kden(ik)
	if (taskid.eq.0.and.ik.eq.2) write (6,*) 'ekin(2)=',ekin(2)
	if (ekin(ik).gt.0.) ekin(ik)=alog10(ekin(ik))
	end do
        if (taskid.eq.0) write (6,"((3x,1p,10e13.5))") (ekin(ii),ii=1,nxh)
                        
         end if                 !(kinit(1) .eq. -3)
c
c	if (taskid.eq.0) write (6,*) 'ekin(2)=',ekin(2)
	

 101  format ((3x,10e13.5))
 501     format (1x)
      call MPI_BCAST (iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     
!     loop over z and y
!     
!         do 1000 zg=1,numtasks
!            if (taskid.eq.zg-1) then
!     
!     if (taskid.gt.0) call MPI_RECV (iseed,1,MPI_INTEGER,taskid-1,
!     1                 taskid-1,MPI_COMM_WORLD,mpistatus,ierr)
!     if (taskid.gt.0) call MPI_RECV (rseed,2,MPI_INTEGER,taskid-1,
!     1                 taskid-1,MPI_COMM_WORLD,mpistatus,ierr)

      xp = mystart_x
      z = mystart_z
c
      do a=1,num_al

         x = xp + xist-1
         k1u(x)=s1*kx(x)
         
         k3u=s3*kz(z)               

         if(x .eq. 1 .and. z .eq. 1) then
            yst = 2
         else
            yst = 1
         endif

         do 15 y=yst,ny
            if(.not.mask(y,a)) goto 15
            
            k2u=s2*ky(y)
            kyzsq=k2u**2+k3u**2
               k=sqrt(k1u(x)**2+kyzsq)
               kd=1./k
               
	if (kinit(1).eq.-1) then
               ekx=e( k, kmax) * voljac
               ef=sqrt( ekx ) * rpid*kd
        else if (kinit(1).eq.-3) then
                     ef=0.0
         ik=k/deltak+1.5
         frac=k+1-ik
                     ekx=frac*ekin(ik+1)+(1.-frac)*ekin(ik)
                     ekx=10.**ekx
                     if (k.gt.kmax) ekx=0.
                     ekx=ekx*voljac
                     ef=sqrt( ekx ) * rpid*kd
	end if
                     
               
#ifdef RANDOM_SEED
        call random_number (theta1)
        call random_number (theta2)
        call random_number (phi)
#else
               theta1=ranu(1)
               theta2=ranu(1)
               phi=ranu(1)
#endif
               
	if (taskid.eq.0.and.a.eq.1) then
c	write (6,*) 'field: x,y,z,ef,theta1=',x,y,z,ik,ekx,ef,theta1,ekin(ik),kden(ik)
	end if

               k12=sqrt(k1u(x)**2+k2u**2)
               
               angphi=twopi*phi
               argc=twopii*theta1
                     alpha=ef*cexp(argc)*cos(angphi)
                     argc=twopii*theta2
                     beta =ef*cexp(argc)*sin(angphi)
                     
!     
!     **caution**: here we are getting u1/s1, u2/s2 and u3/s3
!     in wavenumber space (compare with rogallo 1981, p.53)
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
                              
                              
 15            continue

               call next_xz(xp,z)
            enddo
            
            if (taskid .eq. 0) then
               uny(1,1,1)=0.
               uny(1,1,2)=0.
               uny(1,1,3)=0.
            end if
!     x
!     if (taskid.lt.nzg-1) call MPI_SSEND (iseed,1,MPI_INTEGER,
!     1                     taskid+1,taskid,MPI_COMM_WORLD,mpierr)
!     if (taskid.lt.nzg-1) call MPI_SSEND (rseed,2,MPI_INTEGER,
!     1                     taskid+1,taskid,MPI_COMM_WORLD,mpierr)
!     
!               end if
!
! 1000       continue
!     
!     
         end if
!    	 

!     scalar fields
!     scalars having kinit=-1 are initialised here, with initial
!     spectra same as the 3-d energy spectrum
!     scalars having kinit=-2 are left undone, to be handled by
!     sub. inscf, called from sub. incond after return from this routine
!     
#ifdef NOSCALAR
#else
         if (nc.le.0) go to 3000

c read input scalar spectrum (if indicated)
c
         if (any(kinit(4:3+nc).eq.-3)) then
       allocate (ekinsc(nxhp))
	if (taskid.eq.0) then
            call fopen1 (3,'ekinp.sc','formatted  ')
            read (3,501)
c           read (3,101) (ekinsc(ik),ik=1,nxh)
            read (3,*) (ekinsc(ik),ik=1,nxh)
            close (3)
	end if
	call MPI_BCAST (ekinsc,nxh,mpireal,0,MPI_COMM_WORLD,ierr)

	do ik=2,nxh
	ekinsc(ik)=ekinsc(ik)*kden(ik)
	if (ekinsc(ik).gt.0.) ekinsc(ik)=alog10(ekinsc(ik))
	end do
                        
         end if                 !(kinit(1) .eq. -3)
c
	

         do 2 i=1,nc
            j=i+3
            if (kinit(j).ne.-1.and.kinit(j).ne.-3) go to 2
            if (taskid.eq.0) write (6,601) i
 601        format (' scalar no. ',i1,'  initialized as Gaussian')
            call ranseq (kran(i),luran)

	iseed=iseed+taskid
 	if (taskid.le.1) write (6,*) 'scalars: taskid,i,iseed=',taskid,i,iseed
	
	rseed(1)=iseed
#ifdef RANDOM_SEED
      call random_seed (put=rseed(1:seed_size))
#else
         call random_seed (size=ir)
         call random_seed (put=rseed(1:2))
#endif

!            do 2000 zg=1,numtasks
!               if (taskid.eq.zg-1) then
!     
!     if (taskid.gt.0) call MPI_RECV (iseed,1,MPI_INTEGER,taskid-1,
!     1                 taskid-1,MPI_COMM_WORLD,mpistatus,mpierr)
!     if (taskid.gt.0) call MPI_RECV (rseed,2,MPI_INTEGER,taskid-1,
!     1                 taskid-1,MPI_COMM_WORLD,mpistatus,mpierr)

            xp = mystart_x
            z = mystart_z
            do a=1,num_al
               
               x = xp + xist-1
               k1u(x)=s1*kx(x)

               k3u=s3*kz(z)
               
               if (z.eq.1 .and. x .eq. 1) then
                  yst = 2
               else
                  yst = 1
               endif
               

               do 20 y=yst,ny
                  if(.not.mask(y,a)) goto 20
                  
                  k2u=s2*ky(y)
                  k=sqrt(k1u(x)**2+k2u**2+k3u**2)
               kd=1./k
                  ef=0.0
c
		if (kinit(j).eq.-1) then
                  ekx=e( k, kmax ) *voljac
                  if( k .ne. 0.0 ) ef=sqrt( ekx /(twopi) )/k
		else if (kinit(j).eq.-3) then
         	ik=k/deltak+1.5
         	frac=k/deltak+1-ik
                     ekx=frac*ekinsc(ik+1)+(1.-frac)*ekinsc(ik)
                     ekx=10.**ekx
                     if (k.gt.kmax) ekx=0.
c                     ekx=ekx*voljac/2. 
                     ekx=ekx*voljac
                     ef=sqrt( ekx) * rpid * kd 
		end if

#ifdef RANDOM_SEED
        call random_number (phase)
	if (taskid.eq.0.and.a.eq.1.and.y.eq.yst) then
	write (6,*) 'field: i,phase=',i,phase
	end if
#else
                  phase=ranu(kran(i))
#endif
               
c                  phase=ranu(kran(i))
                  ctmp=ef*cexp(twopii*phase)
                  uny(y,a,i+3)=ctmp
                  
                  
 20            continue
               
               call next_xz(xp,z)
            enddo

!     
!     if (taskid.lt.nzg-1) call MPI_SSEND (iseed,1,MPI_INTEGER,
!     1                     taskid+1,taskid,MPI_COMM_WORLD,mpierr)
!     if (taskid.lt.nzg-1) call MPI_SSEND (rseed,2,MPI_INTEGER,
!     1                     taskid+1,taskid,MPI_COMM_WORLD,mpierr)

!               end if
! 2000       continue

 2       continue

 3000    continue

#endif
!     
!     pass the current random no. seed back to task 0
!     
!     if (taskid.eq.nzg-1) call MPI_SSEND (iseed,1,MPI_INTEGER,
!     1                     0,nzg-1,MPI_COMM_WORLD,mpierr)
!     if (taskid.eq.0) call MPI_RECV (iseed,1,MPI_INTEGER,nzg-1,nzg-1,
!     1                                MPI_COMM_WORLD,mpistatus,mpierr)
         call MPI_BARRIER (MPI_COMM_WORLD,mpierr)


!     enforce conjugate symmetry
!     

cc          call enfsym(uny)
	do i=1,3+nc
	if (kinit(i).eq.-1.or.kinit(i).eq.-3) call consym (uny(1,1,i),1)
	end do
	call check_consym (uny(1,1,1),'exit field')

!     
 95      continue

         deallocate(k1u)
c
	if (allocated(kden)) deallocate (kden)
	if (allocated(ekin)) deallocate (ekin)
	if (allocated(ekinsc)) deallocate (ekinsc)

         return
         end subroutine field
