c     
c     Line made redundant upon introduction of new inpen.f 
c     routine, PK Yeung, Dec 2007 have been removed
c     
c     routine to specify initial conditions
c     
      subroutine incond
      use comsp
      USE IO
      implicit none
      include 'intvars'
     
      character*2 str
      character*5 numer
      character*16 fn
      logical pic(ncd)
      integer i,izr,lstart
      real xnorm,ynorm,znorm,tmpvmaxz,vel
      integer m
      character(len=32) :: str1,str2

      m = 2

#ifdef TEST_CYLIO
      allocate(ucomp(ny,zjsz,xisz,3+nc))
      ucomp = 0
#endif

	if (taskid.eq.0) write (6,*) 'enter incond: b11(m)=',b11(m)
c      write (6,*) 'incond, entry, taskid=',taskid, 'irz=',irz
#ifdef NOYET
c     check scalar I.C. parameters
c     
      if (nc.gt.1.and.idfdf.eq.1.and.iscref.gt.nc) then
         write (6,*) ' iscref should be less than no. of scalars'
         write (6,*) ' check the block data'
         stop 'aborts in sub. incond'
      end if
#endif	
c     
c     
c----------------specified initial fields (un)  ----------------------
c     
      time0=0.
      istep0=0
      kstep=2
      if( istart .eq. 0 .or. any(kinit.eq.-1).or.any(kinit.eq.-3) ) then
c     
         un(:,:,:,1:3+nc)=0.
c     
#ifdef NOYET
c     initialize any scalars with kinit=-2 or -7 first:   
c     in MPL version, such scalars must be placed at the beginning
c     of the list
c     
         do i=1,nc
            if (kinit(i+3).eq.-2.or.kinit(i+3).eq.-7) go to 5
         end do
         go to 11
c     error if 1st kinit=-1 scalar not at the beginning of the list
 5       if (i.ne.1) go to 9
         do j=i,nc
            pic(j)=.false.
            if (kinit(j+3).eq.-2.or.kinit(j+3).eq.-7) pic(j)=.true.
         end do
c     the logical flag should not change in value more than once
c     from the 1st kinit=-2 scalar to the end of the list
c     (otherwise this indicates scalars with kinit=-2 or not are
c     illegally interpersed among each other)
         mm=0
         do j=i,nc-1
            if (.not.(pic(j).and.pic(j+1))) mm=mm+1
         end do
         if (mm.gt.1) go to 9
         call inscf
c         write (6,*) 'incond: after inscf'
c     
c     option for initial scalar fields similar to Kronenburg & Bilger
c     
         if (kinit(i+3).eq.-7) call insckb
#endif	
c     
         go to 11
 9       write (6,*) 'all scalars with kinit=-2 must be placed first'
         write (6,*) 'edit ebdata.f and try again'
         stop 'aborts in sub. incond'
 11      continue
c     
c     initialize in Fourier space or read in previous data
c     
         call field(un)
c     
      end if
c     
#ifdef NOYET
#ifndef SCALAR
      if (istart.gt.0.and.kinit(4).eq.-7) then
         call inscf
         call insckb
         write (6,*) 'insckb called'
      end if
#endif
#endif	
c     
#ifdef LINDBORG
c     note that this is an initialization call: contents of un are not changed
      call force_lindborg (un,1)
c      write (6,*) 'after force_lindborg, taskid=',taskid
#endif
      if (taskid.eq.0) write (6,*) 'incond: irz=',irz
      if (any(kinit.gt.0)) then
         if (all(kinit.gt.0)) then
#ifdef AUG1912
            call inpen_box (un)
#else
            if (hdf5_input_flag) then
               CALL RESTART_HDF5_DRIVER(TRIM(ADJUSTL(indir_fn)),hdf5_init,un)
            else
               if (fninit(1)(4:4).eq.'r') then
                  call incomm (un)
               else
                  call inpen (un)
               end if
            end if
#endif
         else
            if (hdf5_input_flag) then
               CALL RESTART_HDF5_DRIVER(TRIM(ADJUSTL(indir_fn)),hdf5_init,un)
            else
               if (fninit(1)(4:4).eq.'r') then
                  call incomm (un)
               else
                  call inpen (un)
               end if
            end if
         end if
         call time_stamp ('incond: after inpen')
      endif
c
#ifdef FEK_FORC
!     if (kinit(1).gt.0) 
      call ek_init (un)
#endif
	call time_stamp ('incond: before do 110')
      do 110 i=1,3
         u(:,:,:,i)=un(:,:,:,i)
 110  continue
	call time_stamp ('incond:  after do 110')

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	call time_stamp ('incond: at barrier')
      
#ifndef NOSCALAR
!     option for differential diffusion: set one scalar to be
!     equal to the other at t=0
      
      do i=4,3+nc
         if (kinit(i).le.-10) then
            iscref=abs(kinit(i))-10
            u(:,:,:,i)=u(:,:,:,3+iscref)
            un(:,:,:,i)=un(:,:,:,3+iscref)
         end if
      end do
#endif
      
#ifdef NOYET
      do j=4,3+nc
         
         if (kinit(j).ne.-2) then
            do zp=1,mz
               do y=1,ny
                  do x=1,nx
                     uxyr(x,y,zp,j)=unr(x,y,zp,j)
                  end do
               end do
            end do
         end if
c     
         if (kinit(j).eq.-2) then
            do zp=1,mz
               do y=1,ny
                  do x=1,nx
                     unr(x,y,zp,j)=uxyr(x,y,zp,j)
                  end do
               end do
            end do
         end if
c     
      end do
c     
#endif
c     call inpvf to specify velocity field in physical space
c     if kinit(1)=-2 (in this case kshift must be equal to 4,
c     i.e., zero-shifts)
c     
      if (kinit(1).eq.-2) then
c     
#ifndef BOUSS
         if (kshift.ne.4.and.nsteps.ge.1) then
            write (6,*) 'initial velocity field to be in physical space'
            write (6,*) 'but zero shifts must be selected,'
            write (6,*) 'i.e., use kshift=4 in ebdata'
            stop 'program aborts in sub. incond'
         end if
#endif
c     

c
         call inpvf(un,u)
c
      end if
c     
#ifdef NOYET
c     
c     if inorm.ne.0, normalize the scalars to make the variance equal
c     to a desired value
c     
      do 300 i=1,nc
         if (kinit(i+3).eq.0) go to 300
         if (inorm.eq.0) go to 300
         if (inorm.eq.-1.and.istart.eq.1) go to 300
         is=i+3
         call sptvar (is,rms(i+3))
         fmult=sqrt(svnorm(i))/rms(i+3)
         do 350 zp=1,mz
            do 350 y=1,ny
               do 350 x=1,nx
                  unr(x,y,zp,is)=unr(x,y,zp,is)*fmult
                  uxyr(x,y,zp,is)=uxyr(x,y,zp,is)*fmult
 350           continue
c     
 300        continue
c     
c     if idfdf=1,
c     make all scalars identical to the iscref-th
c     
            if (idfdf.eq.1.and.nc.gt.1.and.istart.eq.0) then
c     
               do 1020 i=4,3+nc
                  if (i-3.ne.iscref) then
                     j=iscref+3
                     do 1025 zp=1,mz
                        do 1025 y=1,ny
                           do 1025 x=1,nx
                              unr(x,y,zp,i)=unr(x,y,zp,j)
                              uxyr(x,y,zp,i)=uxyr(x,y,zp,j)
 1025                      continue
                        end if
 1020                continue
c     
                  end if
#endif	
c     
c     
c----------------restart from checkpointed data ----------------------
c     
	call time_stamp ('incond: before store')
	
c     if (nx.ge.64) write (6,*) ' before call store'
                  if( istart .ne. 0 ) then
                     lstart = iabs( istart )
                     call store( 3 ,lstart )
c                     if( istart .gt. 0 ) then
                    if( istart .gt. 0.and.all(kinit.ge.0) ) then
                        
                        if (taskid.eq.0) rewind (luran1)
                        call ranseq( -1 , luran1 )
                     endif
                  endif
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
                   if (taskid.eq.0) write (6,*) 'initial data fields in place'    
#ifdef EPFOR
                  if (istart.eq.0.and.tfiuo.lt.0.) then
                     call store (4,1)
                  end if
#endif
c     
c     initialize the time
c     
                  time=time0
                  istep=istep0
c     
c     produce spectra
c     
	if (taskid.eq.0) write (6,*) 'incond: b11(2)=',b11(2)
	if (b11(2).ne.beta1**2) then
	call init
	call waveno
	un(:,:,:,1)=un(:,:,:,1)/beta1
	un(:,:,:,2)=un(:,:,:,2)/beta2
	un(:,:,:,3)=un(:,:,:,3)/beta3
	if (taskid.eq.0) write (6,*) 'incond: b11(2)=',b11(2)
	end if


! spectrum calculated by sptvar routine is at present not correct
! except for the case of a 2pi^3 grid. Also must have nx >= ny & nz in
! order for the arrays in sptvar to not overflow.
      if (beta1.eq.1.0.and.beta2.eq.1.0.and.beta3.eq.1.0
     1    .and.nx.ge.ny.and.nx.ge.nz) then
	      call sptvar (un)
      end if
      str2 = 'AXI_SPTR'
	if (ioaxi(1).gt.0) then
      str1 = 'axi_kx'
      CALL AXI_SPTR(un,1,str1,str2)
	end if
	if (ioaxi(2).gt.0) then
      str1 = 'axi_ky'
      CALL AXI_SPTR(un,2,str1,str2)
	end if
	if (ioaxi(3).gt.0) then
      str1 = 'axi_kz'
      CALL AXI_SPTR(un,3,str1,str2)
	end if

        call sptr(un,1,1)
	call check_consym2 (un,3,'incond',3)
      return
      end
