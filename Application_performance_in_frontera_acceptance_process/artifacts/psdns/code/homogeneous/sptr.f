 	subroutine sptr (uny,m,io)
c
      use comsp
c
      implicit none
	include 'intvars'

	complex(b8) :: uny(ny,zjsz*xisz,3+nc)
	integer :: io,i,j,k,ij,ncount,ik,m,a,get_l
	real(b8) :: two9,rk2,term,dterm,factor,size,s23,s31,s12,ssum,qs
	real(b8) :: grdiss,fact
	real(b8), allocatable :: ssq(:)
	real(b8), allocatable :: kxn2(:),kyn2(:),kzn2(:)
c
	complex(b8), allocatable :: usum(:),uiuj(:,:),bk2y(:,:)
	integer :: iz,iy,nzp2,nyp2,xst
	complex(b8) :: ysum,czero,uiuj1,usum1

	integer iii,jjj,ijc

	complex(b8), allocatable :: xrijy(:,:,:),yrijy(:,:,:),zrijy(:,:,:)

      	integer ithr,OMP_GET_THREAD_NUM
c
	real(b8) beta_min
	
#ifdef SHELL_DK
        beta_min=min(beta1,beta2,beta3)
#endif


	if (taskid.eq.0) write (6,*) 'enter sptr, io,b11(m)=',io,b11(m)
	two9=2./9.
	czero=cmplx(0.,0.)
	nzp2=nz+2
	nyp2=ny+2
	
	allocate (kxn2(nxh),kyn2(ny),kzn2(nz))
	allocate (xrijy(nxhp,ncp,0:num_thr-1))
	allocate (yrijy(nyhp,ncp,0:num_thr-1))
	allocate (zrijy(nzhp,ncp,0:num_thr-1))
	allocate (eijky(mxyz,ncp,0:num_thr-1))
	allocate (dijky(mxyz,ncp,0:num_thr-1))
	allocate (sijky(mxyz,ncp,0:num_thr-1))
	allocate (gijky(mxyz,ncp,0:num_thr-1))
	allocate (grijky(mxyz,ncp,3,0:num_thr-1))
	allocate (ekky(mxyz,0:num_thr-1))
	allocate (der2y(ncp,0:num_thr-1))

	allocate (uiuj(ny,0:num_thr-1))

!	tfact(:)=2.
!	tfact(1)=1.

	allocate(ssq(3+nc))
	
      ssq(1)=b11(m)
      ssq(2)=b22(m)
      ssq(3)=b33(m)
      if (nc.gt.0) then
      do 2 i=4,3+nc
      ssq(i)=1
 2    continue
      end if

	xrijy(:,:,:)=0.
	yrijy(:,:,:)=0.
	zrijy(:,:,:)=0.

 	der2y(:,:)=0.
	eijky(:,:,:)=0.
	dijky(:,:,:)=0.
	sijky(:,:,:)=0.
	gijky(:,:,:)=0.
	grijky(:,:,:,:)=0.
	ekky(:,:)=0.

	laak(:,:)=0.
	taak(:,:)=0.
	raak(:)=0.
	tmre(:)=0.
	rms(:)=0.
	ett(:)=0.

	kxn2(:)=(kx(:)/nx)**2
	kyn2(:)=(ky(:)/ny)**2
	kzn2(:)=(kz(:)/nz)**2

	kx2(:)=b11(m)*kx(:)**2
	ky2(:)=b22(m)*ky(:)**2
	kz2(:)=b33(m)*kz(:)**2

	kyn2(nyhp)=0.
	kzn2(nzhp)=0.
	ky2(nyhp)=0.
	kz2(nzhp)=0.

	ij=0
	do i=1,3+nc
	do j=i,3+nc
	ij=ij+1
 	if (i.eq.j) kt(i)=ij
	end do
	end do

	if (iovor.ge.1) call vorsp2 (1,0,m)
c
	allocate (bk2y(ny,0:num_thr-1))

	ithr=0
!$OMP PARALLEL private (ithr,x,xp,z,ij,tfact_x,ysum,term,rk2,ik,dterm)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif


      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 5 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1

         x = xp + xist-1

	 if(x .eq. 1) then
	    tfact_x=1
	 else
	    tfact_x=2
	 endif


	ij=0
	do 20 i=1,3+nc
	do 20 j=i,3+nc
	ij=ij+1

	 ysum=czero
	 do 10 y=1,ny

	    if(.not.mask(y,a)) goto 10
		 

	    uiuj(y,ithr)=real(uny(y,a,i)*conjg(uny(y,a,j)))
	    ysum=ysum + uiuj(y,ithr)
c
	    if (y.eq.1) then       
	       yrijy(y,ij,ithr)=yrijy(y,ij,ithr) + uiuj(y,ithr)+(tfact_x-1)*conjg(uiuj(y,ithr))
	    else if (y.lt.nyhp) then  
	       yrijy(y,ij,ithr)=yrijy(y,ij,ithr)+uiuj(y,ithr)
	    else 
	       iy=nyp2-y              
	       yrijy(iy,ij,ithr)=yrijy(iy,ij,ithr)+(tfact_x-1)*conjg(uiuj(y,ithr))   
	    end if                 
		 
	    if(.not. mask(y,a)) goto 10
	    term=tfact_x*uiuj(y,ithr)
	    rk2=kx2(x)+ky2(y)+kz2(z)
c
#ifdef SHELL_DK
        ik=sqrt(rk2)/beta_min+1.5
#else
        ik=sqrt(rk2)+1.5
#endif
#ifdef MODEL_SPECTRUM
      IF (ik .GT. mxyz) GOTO 10
#endif
c
#ifdef SHEAR
#else

#ifdef LINDBORG
	    dterm=(visc_h*(kx2(x)+ky2(y))**4+visc_v*kz2(z)**4)*term
#else
	    dterm=rk2*term
#endif

	    eijky(ik,ij,ithr)=eijky(ik,ij,ithr)+term
	    dijky(ik,ij,ithr)=dijky(ik,ij,ithr)+dterm
	    sijky(ik,ij,ithr)=sijky(ik,ij,ithr)+dterm*rk2
	    gijky(ik,ij,ithr)=gijky(ik,ij,ithr)+dterm*rk2**2
	    der2y(ij,ithr)=der2y(ij,ithr)+ky2(y)*term
	    if (i.eq.j.and.i.le.3.and.rk2.ne.0.) then
	       ekky(ik,ithr)=ekky(ik,ithr)+term/sqrt(rk2)*ssq(i)
	    end if

c
c also spectra of individual gradients
c
	    grijky(ik,ij,1,ithr)=grijky(ik,ij,1,ithr)+term*kx2(x)
	    grijky(ik,ij,2,ithr)=grijky(ik,ij,2,ithr)+term*ky2(y)
	    grijky(ik,ij,3,ithr)=grijky(ik,ij,3,ithr)+term*kz2(z)

#endif

 10	 continue		! do y=1,ny

	 if (z.eq.1) then       
	    zrijy(z,ij,ithr)=zrijy(z,ij,ithr)+ysum+ (tfact_x-1)*conjg(ysum)      
	 else if (z.lt.nzhp) then  
	    zrijy(z,ij,ithr)=zrijy(z,ij,ithr)+ysum    
	 else if (z.gt.nzhp) then  
	    iz=nzp2-z              
	    zrijy(iz,ij,ithr)=zrijy(iz,ij,ithr)+(tfact_x-1)*conjg(ysum)   
	 end if                 

	 xrijy(x,ij,ithr)=xrijy(x,ij,ithr)+ysum


 20	continue		! do 20 i,j=1,3+nc

c
	 call next_xz(xp,z)
c
 5	continue

	if (iovor.ge.1)  call vorsp1 (uny,m,ithr,bk2y(1,ithr))
c
c
!$OMP END PARALLEL
!$OMP BARRIER
c
	call sumthr (eijky,mxyz,ncp,1,num_thr)
	call sumthr (dijky,mxyz,ncp,1,num_thr)
	call sumthr (sijky,mxyz,ncp,1,num_thr)
	call sumthr (gijky,mxyz,ncp,1,num_thr)
	call sumthr (grijky,mxyz,ncp,3,num_thr)
	call sumthr (ekky,mxyz,1,1,num_thr)
	call sumthr (der2y(1,0),ncp,1,1,num_thr)
	call csumthr (xrijy(1,1,0),nxhp,ncp,num_thr)
	call csumthr (yrijy(1,1,0),nyhp,ncp,num_thr)
	call csumthr (zrijy(1,1,0),nzhp,ncp,num_thr)
c
	dijky(:,:,0)=dijky(:,:,0)*2*viscos
	sijky(:,:,0)=sijky(:,:,0)*2*viscos
	gijky(:,:,0)=gijky(:,:,0)*2*viscos

	if (iovor.ge.1) then
	 call vorsp2 (0,1,m)
	if (taskid.eq.0) call vorout (35)
	end if
c
	deallocate (bk2y)
	deallocate (kxn2,kyn2,kzn2)
	deallocate(ssq)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c collect contributions from all processors

	allocate (der2ij(ncp),stat=ierr)
	allocate (eijk(mxyz,ncp),stat=ierr)
	allocate (dijk(mxyz,ncp),stat=ierr)
	allocate (sijk(mxyz,ncp),stat=ierr)
	allocate (gijk(mxyz,ncp),stat=ierr)
	allocate (grijk(mxyz,ncp,3),stat=ierr)
	allocate (ekk(mxyz),stat=ierr)

      call MPI_REDUCE (der2y(1,0),der2ij(1),ncp,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)
 
	ncount=mxyz*ncp

 	call MPI_ALLREDUCE (eijky(1,1,0),eijk(1,1),ncount,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)
 	call MPI_ALLREDUCE (dijky(1,1,0),dijk(1,1),ncount,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)
 	call MPI_ALLREDUCE (sijky(1,1,0),sijk(1,1),ncount,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)
 	call MPI_ALLREDUCE (gijky(1,1,0),gijk(1,1),ncount,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)
 	call MPI_ALLREDUCE (grijky(1,1,1,0),grijk(1,1,1),ncount*3,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)
 	call MPI_ALLREDUCE (ekky(1,0),ekk(1),mxyz,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)

! xrij, yrij and zrij are defined in comsp (but allocated and deallocated here)
	allocate (xrij(nxhp,ncp,2))
	allocate (yrij(nyhp,ncp,2))
	allocate (zrij(nzhp,ncp,2))
 	xrij(:,:,io)=0.
 	yrij(:,:,io)=0.
 	zrij(:,:,io)=0.

      ncount=nxhp*ncp
      call MPI_REDUCE (xrijy(1,1,0),xrij(1,1,io),ncount,mpicomplex,
     1                 MPI_SUM,0,MPI_COMM_WORLD,ierr)
      ncount=nyhp*ncp
      call MPI_REDUCE (yrijy(1,1,0),yrij(1,1,io),ncount,mpicomplex,
     1                 MPI_SUM,0,MPI_COMM_WORLD,ierr)
      ncount=nzhp*ncp
      call MPI_REDUCE (zrijy(1,1,0),zrij(1,1,io),ncount,mpicomplex,
     1                 MPI_SUM,0,MPI_COMM_WORLD,ierr)

#ifdef MHD
        call sptr_mhd (uny,m)
#endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c compute other quantities (only taskid=0)
	if (taskid.ne.0) go to 1000

	allocate (s(3+nc))

	s(1)=sqrt(b11(m))
	s(2)=sqrt(b22(m))
	s(3)=sqrt(b33(m))
      s23=s(2)*s(3)
      s31=s(3)*s(1)
      s12=s(1)*s(2)
      size=s(1)*s(2)*s(3) 
c
	do i=4,3+nc
	s(i)=1.
	end do

	ij=0
	do 40 i=1,3+nc
	do 40 j=i,3+nc
	ij=ij+1

	factor=s(i)*s(j)
 	der2ij(ij)=der2ij(ij)*factor

	eijk(:,ij)=eijk(:,ij)*factor
	dijk(:,ij)=dijk(:,ij)*factor
	sijk(:,ij)=sijk(:,ij)*factor
	gijk(:,ij)=gijk(:,ij)*factor

	xrij(:,ij,io)=factor*xrij(:,ij,io)/size
	yrij(:,ij,io)=factor*yrij(:,ij,io)/size
	zrij(:,ij,io)=factor*zrij(:,ij,io)/size

 40	continue

 	do 45 ij=1,ncp

	xrij(:,ij,io)=s23*xrij(:,ij,io)
	yrij(:,ij,io)=s31*yrij(:,ij,io)
	zrij(:,ij,io)=s12*zrij(:,ij,io)

	ssum=real(xrij(1,ij,io))
	do x=2,nxh
	ssum=ssum+2.*real(xrij(x,ij,io))
	enddo
	corr(ij)=s(1)*ssum

 45	continue

 	qs=corr(1)+corr(4+nc)+corr(6+2*nc)

	size=1
	factor=size/2.
      ek(:)=(eijk(:,1)+eijk(:,4+nc)+eijk(:,6+2*nc))*factor
      ek(1)=0.
      dk(:)=(dijk(:,1)+dijk(:,4+nc)+dijk(:,6+2*nc))*factor
	dk(1)=0.
	sk(:)=(sijk(:,1)+sijk(:,4+nc)+sijk(:,6+2*nc))*factor

	tke=sum(ek)
	epslon=sum(dk)
	skew=sum(sk)
c
      if (epslon.gt.0) then
        skew=2./35.*(15.*viscos/epslon)**(1.5)*skew
c compute kolmorgorov length,velocity & time scales
        klen=(viscos**3/epslon)**.25
	  kvel=(viscos*epslon)**.25
	  ktime=(viscos/epslon)**.5
	else
        write (6,*) 'warning: dissipation is non-positive!'
	  write (6,*) 'viscos=',viscos
      end if
c
c scalar dissipation rate
c
#ifndef NOSCALAR
	do 28 i=1,nc
	dijk(:,kt(i+3))=dijk(:,kt(i+3))/pr(i)
	scdiss(i)=sum(dijk(:,kt(i+3)))
 28	continue
#endif
c

c compute integral length scales in all 3 directions for each ij-th products
	lijk(:,:)=0.
	do 47 ij=1,ncp
	if (corr(ij).ne.0.) then
	  lijk(ij,1)=pi*real(xrij(1,ij,io))/corr(ij)
	  lijk(ij,2)=pi*real(yrij(1,ij,io))/corr(ij)
	  lijk(ij,3)=pi*real(zrij(1,ij,io))/corr(ij)
	endif
 47	continue

c compute an integral scale as the integral of e(k)/k
c
      ekl=0.
      do 48 a=1,nxh
 48   ekl=ekl+ekk(a)/2.
      ekl=1.5*pi/(2*tke)*ekl

	do 49 a=1,3
	do 49 i=1,3+nc
 49	laak(i,a)=lijk(kt(i),a)

c take the square roots of the variances to get rms values

 	do 50 i=1,3+nc
 50	rms(i)=sqrt(corr(kt(i)))

c compute integral scale reynolds numbers

      if (viscos.gt.0.) then
      do 51 i=1,3
 51   raak(i)=laak(i,i)*rms(i)/viscos
      end if

c calculate eddy-turnover times
c
      do 34 i=1,3
      if (rms(i).ne.0) ett(i)=laak(i,i)/rms(i)
 34   continue
c
c tmij is to stand for the mean square of the      
c 1st derivatives (ith field variable w.r.t. jth dirn).                 
c each individual mode is multiplied by two because of                  
c conjugate symmetry                               
c                                                  
      do 35 i=1,3+nc                               
c                                                  
      ssum=0.                                       
      do 36 x=2,nxh                                
 36   ssum=ssum+kx(x)**2*real(xrij(x,kt(i),io))      
      tmij(i,1)=2.*s(1)**3*ssum                     
c                                                  
      tmij(i,2)=der2ij(kt(i))                   
c                                                  
      ssum=0. 
      do 37 a=2,nz/2 
 37   ssum=ssum+kz(a)**2*real(zrij(a,kt(i),io))      
      tmij(i,3)=2.*s(3)**3*ssum
c
 35   continue 
c
c use the gradient variances to form dissipation   
c                                                  
      write (6,*) (tmij(i,1),i=1,3)                
      write (6,*) (tmij(i,2),i=1,3)                
      write (6,*) (tmij(i,3),i=1,3)                
      grdiss=tmij(1,1)+tmij(2,2)+tmij(3,3)         
     1      +tmij(1,2)+tmij(2,1)+tmij(1,3)         
     1      +tmij(3,1)+tmij(2,3)+tmij(3,2)         
      grdiss=grdiss*viscos                         
      write (6,*) 'grdiss,size,epslon=',grdiss,size,epslon              

#ifdef MODEL_SPECTRUM
      WRITE(998,400) istep, time, (tmij(i,1),i=1,3),
     1               (tmij(i,2),i=1,3), (tmij(i,3),i=1,3)
 400  FORMAT (I6.6,1P,10ES13.6)
      CLOSE(UNIT=998)
      OPEN(UNIT=998,FILE='tmij',FORM='FORMATTED',POSITION='APPEND')
#endif
c                                                  
c compute taylor microscales based on corr & tmij  
c                                                  
c for scalars, multiply the taylor scales by sqrt(2.) to bring them     
c in line with the conventional definition, which was also used in      
c e&p fda-87-12                                    
! note that this factor is not used in Mike Rogers' thesis
c                                                  
      do 38 i=1,3+nc                               
      do 38 j=1,3+nc                               
      fact=1.                                      
      if (i.gt.3.or.j.gt.3) fact=sqrt(2.)          
      do 38 a=1,3                                  
	taylor(i,j,a)=0.
      if(tmij(j,a).ne.0.) taylor(i,j,a)=sqrt(corr(kt(i))/tmij(j,a))*fact
 38   continue                                     
c                                                  
      do 39 i=1,3+nc                               
      do 39 j=1,3                                  
 39   taak(i,j)=taylor(i,j,j)                      
c                                                  
c compute taylor microscale reynolds numbers       
c                                                  
      if (viscos.gt.0.) then                       
      do 46 i=1,3                                  
 46   tmre(i)=taak(i,i)*rms(i)/viscos              
      end if                                       
c
#ifndef NOSCALAR
c obtain scalar gradient variances in wavenumber space
c
      ij=kt(4)-1
      ijc=0
      do 70 i=1,nc
      do 70 j=i,nc
      ij=ij+1
      ijc=ijc+1
      do 71 a=1,3
      scgvar(ijc,a)=0.
      do 75 k=1,nxh
      scgvar(ijc,a)=scgvar(ijc,a)+grijk(k,ij,a)
 75   continue
 71   continue
 70   continue
c
#endif

c
	deallocate(s)
c
#ifdef OUT1D
	call out1d (13,14,io)
#endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write eulstat
c
#ifdef FEK_FORC
        if (io.eq.1.and.jstep.ge.1) go to 1000
#endif
	call eulout
	call covcof
c
#ifndef NOSCALAR
	call escout (16,io,uny)
      if (nc.eq.2) call scjout (1,2,27,io)
      if (nc.eq.3) then
      call scjout (1,2,27,io)
      call scjout (1,3,28,io)
      call scjout (2,3,29,io)
      end if
#endif
c
#ifdef ROTD
      if (rrate.gt.0.) then
      rossby=epslon/qs/rrate
      write (6,*) 'turbulent Rossby no.=',rossby
      end if
#endif

 1000	continue
c ZHAI: 09/18/2016 anisotropy tensor decomposition
c need tke and ek for taskid 0 in bijdecomp
c see "On the two-dimensionalization of quasistatis MHD turb"
c
        call bijdecomp(uny)
c
#if defined (FEK_FORC) || defined (ROTD)
cc      if (io.eq.1.and.jstep.ge.1) go to 1002
#ifdef BOXPHY
        call dissenst (un,un,u,u)
#endif
 1002   continue
#endif
c
#ifdef FEK_FORC
c       if (taskid.eq.0) write (6,*) 'calling ek_shell from sptr'
        if (jstep.eq.0) call ek_shell (uny,m)
#endif

#ifndef NOSCALAR
      call MPI_BCAST (scgvar,ncgd*3,mpireal,0,MPI_COMM_WORLD,mpierr)
#endif

        !used in BOUSS calls, otherwise empty hook (see null_hooks)
!	call hsptr (uny,m)

        deallocate (eijky,eijk)
     	deallocate (dijky,dijk)
     	deallocate (sijky,sijk)
     	deallocate (gijky,gijk)
     	deallocate (grijky,grijk)
	deallocate (ekky,ekk)
	deallocate (uiuj)
	deallocate (xrijy,yrijy,zrijy)
	deallocate (xrij,yrij,zrij)
	deallocate (der2y,der2ij)

c
 	if (taskid.eq.0) write (6,*) ' exit sptr'

	return
	end
c
	subroutine sumthr (array,n1,n2,n3,num_thr)
	use precision
	real(b8) array(n1,n2,n3,0:num_thr-1)
	if (num_thr.eq.1) return
	do k=1,n3
	do j=1,n2
	do i=1,n1
	sum=array(i,j,k,0)
	do ithr=1,num_thr-1
	sum=sum+array(i,j,k,ithr)
	end do
	array(i,j,k,0)=sum
	end do
	end do
	end do
	return
	end
c
	subroutine csumthr (array,n1,n2,num_thr)
	use precision
	complex(b8) array(n1,n2,0:num_thr-1)
	complex(b8) sum
	if (num_thr.eq.1) return
	do j=1,n2
	do i=1,n1
	sum=array(i,j,0)
	do ithr=1,num_thr-1
	sum=sum+array(i,j,ithr)
	end do
	array(i,j,0)=sum
	end do
	end do
	return
	end
