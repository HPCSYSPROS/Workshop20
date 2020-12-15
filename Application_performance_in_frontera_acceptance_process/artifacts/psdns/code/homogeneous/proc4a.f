! Version with Cylindrical representation 
!
      subroutine proc4a (uy,uny) 
        use comsp        
        implicit none
	include 'intvars'
!new datastructure
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc) 
        complex(b8) untemp
        complex caux
        real(b8) s1,s2,s3,ssq(3),rks
        real(b8) term
        real(b8) factor,bk12
        integer i,j,ipsi,yst,k,iy,iz,j1,m,a,xz(2),ymax,almax
        real xnorm,ynorm,znorm,tmpvmaxz,vel,vm,vm2,tmpvmaxz2

	integer ithr,omp_get_thread_num
	integer modpress
 
	real(b8) rk1,rk2,rk3,rksq

#ifdef ROTD
	complex(b8) si(2),snl(2),cifac,v(3),
     1		                    oldsi(2)
        complex(b8), allocatable :: rk1i(:)
	complex(b8) rk2i
	real(b8) fact12,rk3sq,rkmag,rkmagr,rk13sq
        real(b8), allocatable :: rk12(:),rk12r(:)
        real(b8), allocatable :: rk23(:),rk23r(:)
#elif defined BIF
        complex(b8), allocatable :: phi(:,:),oldphi(:,:),pnl(:,:)
	real fact12,rk3sq,rkmag,rkmagr,rk3r,rk13sq
	real factk,factx,cosvt,sinvt,gfact,gfac2
        real(b8), allocatable :: rk12(:),rk12r(:),rk2r(:)
#else
        real(b8), allocatable :: facty(:,:)
#endif

#ifndef MHD
	integer, parameter:: imaxis = 1
#endif

	logical flag

        m = 2
c
      ipsi=4+nc
!
#ifdef ROTD
	if (iraxis.eq.3) then
        allocate(rk12(ny),rk12r(ny),rk1i(nxh))
	else if (iraxis.eq.1) then
        allocate(rk23(ny),rk23r(ny))
	end if
#elif defined BIF
c        allocate (phi(xisz,2),oldphi(xisz,2),pnl(xisz,2))
        allocate(rk12(ny),rk12r(ny),rk2r(ny))
#else
        allocate(facty(ny,0:num_thr-1))
#endif
!
!c added on 6/26/96: energy input rate from E&P forcing (predictor only)
 
#ifdef EPFOR
      if (kforce.gt.0.) efkz(:,:)=0.
#endif
 
      s1=sqrt(b11(2))
      s2=sqrt(b22(2))
      s3=sqrt(b33(2))
      ssq(1)=b11(2)
      ssq(2)=b22(2)
      ssq(3)=b33(2)
 
#ifdef ROTD
c     rk1i(:)=imagi*s1*kx(:)
#endif

        ithr=0
!$OMP PARALLEL private (ithr,x,xp,z,factor,untemp,yst,bk12,y1,rk1,rk2,rk3,rksq)

c Note more may have to be added to private list for ROTD and BIF cases
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
         x = xp + xist-1


#ifdef ROTD 
         rk1=s1*kx(x)
         rk3=s3*kz(z)

         do 16 y=1,ny
            if(.not.mask(y,a)) goto 16

            rk2=s2*ky(y)
		if (iraxis.eq.3) then
            rk12(y)=sqrt(rk1**2+rk2**2)
            rk12r(y)=1.0/rk12(y)
		else if (iraxis.eq.1) then
            rk23(y)=sqrt(rk3**2+rk2**2)
            rk23r(y)=1.0/rk23(y)
		end if
            
 16     continue
#elif defined BIF

         rk1=s1*kx(x)
         do 16 y=1,ny
            if(.not.mask(y,a)) goto 16

            rk2=s2*ky(y)
            rk12(y)=sqrt(rk1**2+rk2**2)
            rk12r(y)=1.0/rk12(y)
            rk2r(y) = 1.0/rk2
            
 16     continue
#endif
      
!
#ifdef ROTD
      rk3=s3*kz(z)
      rk3sq=rk3**2
      rk13sq=rk1**2+rk3**2
#elif defined BIF
      rk3=s3*kz(z)
      rk3r=1.0/rk3
      rk3sq=rk3**2
      rk13sq=rk1**2+rk3**2
#endif


 
! define product of mask and integrating factors
 
! update the velocities
!
#ifdef ROTD
! use auxiliary variables in diagonized form, in order to integrate
! the rotation term exactly
c
	call helical (uy,uny,x,z,a)

#ifdef TEST
	if (iraxis.eq.3) then
c	call helical_z (uy,uny,x,z,a,rk12,rk12r)
c	call helical_3 (uy,uny,x,z,a)
 	call helical_1 (uy,uny,x,z,a)
	else if (iraxis.eq.1) then
c	call helical_x (uy,uny,x,z,a,rk23,rk23r)
 	call helical_1 (uy,uny,x,z,a)
	else if (iraxis.eq.2) then
c	call helical_y (uy,uny,x,z,a)
 	call helical_1 (uy,uny,x,z,a)
	end if
#endif
!
#elif defined BIF
c
c Buoyancy: use integrating factors for linear terms, re: method
c of Jack Herring based on a Craya decomposition in Fourier space.
c
c	if (bvfsq.eq.0.) go to 239
c
      if (z.ne.1) then
c
      yst=1
      if (x.eq.1) yst=2
      
      do 210 y=yst,ny
         if(.not.mask(y,a)) goto 210

#ifdef LINDBORG
        factor=xyfac(x,y)*zfac(z)
#else
        factor=xfac(x)*yfac(y)*zfac(z)
#endif

c
c define the two auxiliary variables phi_1 and phi_2
c corresponding nonlinear terms
c
c      rk1=s1*kx(x)
      rk2=s2*ky(y)
c      rk12=sqrt(rk1**2+rk2**2)
      rkmag=sqrt(rk13sq+rk2sq)
      rkmagr=1./rkmag
c      rk12r=1./rk12
      phi(1)=(rk2*uny(y,a,1)-rk1*uny(y,a,2))*rk12r(y)
      phi(2)=(rk1*uny(y,a,1)+rk2*uny(y,a,2))*rk12r(y)*rkmag*rk3r
      pnl(1)=(rk2*uy(y,a,1)-rk1*uy(y,a,2))*rk12r(y)
      factk=rk3*rkmagr*rk12r(y)
      pnl(2)=rk1*uy(y,a,1)*factk+rk2*uy(y,a,2)*factk
     1       -rk12(y)*rkmagr*uy(y,a,3)
c
      oldphi(1)=phi(1)
      oldphi(2)=phi(2)
c
c predictor update of phi_1 and phi_2
c
      phi(1)=(oldphi(1)+pnl(1)*2.)*factor
c
      cosvt=rk3*rkmagr
      sinvt=rk12(y)*rkmagr
      gfac2=cos(bvfsq*dt*sinvt)*factor
      gfact=sin(bvfsq*dt*sinvt)*xdif(x,1)*ydif(y,1)*zdif(z,1)
c      
      phi(2)=gfac2*(oldphi(2)+pnl(2)*2.)
     1        +gfact*(uny(y,a,4)+uy(y,a,4)*2.)
c
c recover the corresponding velocities, the third component by
c enforcing continuity
c
c
      uy(y,a,1)=(rk2*phi(1)+rk1*rk3*rkmagr*phi(2))*rk12r(y)
      uy(y,a,2)=(-rk1*phi(1)+rk2*rk3*rkmagr*phi(2))*rk12r(y)
      uy(y,a,3)=-(rk1*uy(y,a,1)+rk2*uy(y,a,2))*rk3r
c
c
c add values at half-dt
c
      phi(1)=(oldphi(1)+pnl(1))*factor
      phi(2)=gfac2*(oldphi(2)+pnl(2))
     1        +gfact*(uny(y,a,4)+uy(y,a,4))
      uny(y,a,1)=(rk2*phi(1)+rk1*rk3*rkmagr*phi(2))*rk12r(y)
      uny(y,a,2)=(-rk1*phi(1)+rk2*rk3*rkmagr*phi(2))*rk12r(y)
      uny(y,a,3)=-(rk1*uny(y,a,1)+rk2*uny(y,a,2))*rk3r
c
c update the temperature
c
      gfac2=sin(bvfsq*dt*sinvt)*factor
      gfact=cos(bvfsq*dt*sinvt)*xdif(x,1)*ydif(y,1)*zdif(z,1)
c
      untemp=uy(y,a,4)
      uy(y,a,4)=-gfac2*(oldphi(2)+pnl(2)*2.)
     1              +gfact*(uny(y,a,4)+untemp*2.)
      uny(y,a,4)=-gfac2*(oldphi(2)+pnl(2))
     1              +gfact*(uny(y,a,4)+untemp)
c
c
 210  continue
c
c
c modes with kx=ky=0
c
      if (x.eq.1) then
         y=1
#ifdef LINDBORG
        factor=xyfac(x,y)*zfac(z)
#else
        factor=xfac(x)*yfac(y)*zfac(z)
#endif
        do i=1,2
           untemp=uny(y,a,i)
           uny(y,a,i)=factor*(untemp+uy(y,a,i))
           uy(y,a,i) =factor*(untemp+uy(y,a,i)*2.)
        end do
        uny(y,a,3)=0.
        uy(y,a,3)=0.
        factx=xdif(x,1)*ydif(y,1)*zdif(z,1)
        untemp=uny(y,a,4)
        uny(y,a,4)=factx*(untemp+uy(y,a,4))
        uy(y,a,4) =factx*(untemp+uy(y,a,4)+uy(y,a,4))
      end if
c
      else 
c
c kz=0...
c
c      if (y.eq.1) then

         y = 1

#ifdef LINDBORG
        factor=xyfac(x,y)*zfac(z)
#else
        factor=xfac(x)*yfac(y)*zfac(z)
#endif

c
c modes with ky=kz=0
c
c      do 220 xp=1,xisz
c	x=xist+xp-1
c
      phi(1)=-uny(y,a,2)
      phi(2)=-uny(y,a,3)
      pnl(1)=-uy(y,a,2)
      pnl(2)=-uy(y,a,3)
c
      oldphi(1)=phi(1)
      oldphi(2)=phi(2)
c
c predictor update of phi_1 and phi_2
c
      phi(1)=(oldphi(1)+pnl(1)*2.)*factor
c
      gfac2=cos(bvfsq*dt)*factor
      gfact=sin(bvfsq*dt)*xdif(x,1)*ydif(y,1)*zdif(z,1)
c      
      phi(2)=gfac2*(oldphi(2)+pnl(2)*2.)
     1        +gfact*(uny(y,a,4)+uy(y,a,4)*2.)
c
c recover the corresponding velocities, the third component by
c enforcing continuity
c
      uy(y,a,1)=0.
      uy(y,a,2)=-phi(1)
      uy(y,a,3)=-phi(2)
c
c add values at half-dt
c
      phi(1)=(oldphi(1)+pnl(1))*factor
      phi(2)=gfac2*(oldphi(2)+pnl(2))
     1        +gfact*(uny(y,a,4)+uy(y,a,4))
      uny(y,a,1)=0.
      uny(y,a,2)=-phi(1)
      uny(y,a,3)=-phi(2)
c
c update the temperature
c
      gfac2=sin(bvfsq*dt)*factor
      gfact=cos(bvfsq*dt)*xdif(x,1)*ydif(y,1)*zdif(z,1)
c
      untemp=uy(y,a,4)
      uy(y,a,4)=-gfac2*(oldphi(2)+pnl(2)*2.)
     1              +gfact*(uny(y,a,4)+untemp*2.)
      uny(y,a,4)=-gfac2*(oldphi(2)+pnl(2))
     1              +gfact*(uny(y,a,4)+untemp)
c
 220  continue
c
c modes kz=0, but ky.ne.0
c
      else 
c
c      do 230 xp=1,xisz
c	x=xist+xp-1

         do 230 y=2,ny

#ifdef LINDBORG
        factor=xyfac(x,y)*zfac(z)
#else
        factor=xfac(x)*yfac(y)*zfac(z)
#endif

c
c      rk1=s1*kx(x)
      rk2=s2*ky(y)
c      rk2r = 1.0/rk2
      rkmag=sqrt(rk1**2+rk2**2)
      rkmagr=1./rkmag
      phi(1)=(rk2*uny(y,a,1)-rk1*uny(y,a,2))*rkmagr
      phi(2)=-uny(y,a,3)
      pnl(1)=(rk2*uy(y,a,1)-rk1*uy(y,a,2))*rkmagr
      pnl(2)=-uy(y,a,3)
c
      oldphi(1)=phi(1)
      oldphi(2)=phi(2)
c
c predictor update of phi_1 and phi_2
c
      phi(1)=(oldphi(1)+pnl(1)*2.)*factor
c
      gfac2=cos(bvfsq*dt)*factor
      gfact=sin(bvfsq*dt)*xdif(x,1)*ydif(y,1)*zdif(z,1)
c      
      phi(2)=gfac2*(oldphi(2)+pnl(2)*2.)
     1        +gfact*(uny(y,a,4)+uy(y,a,4)*2.)
c
c recover the corresponding velocities, the third component by
c enforcing continuity
c
      uy(y,a,1)=(rkmag*rk2r(y))/(1.+(rk1*rk2r(y))**2)*phi(1)
      uy(y,a,2)=-rk1*rk2r(y)*uy(y,a,1)
      uy(y,a,3)=-phi(2)
c
c add values at half-dt
c
      phi(1)=(oldphi(1)+pnl(1))*factor
      phi(2)=gfac2*(oldphi(2)+pnl(2))
     1        +gfact*(uny(y,a,4)+uy(y,a,4))
      uny(y,a,1)=(rkmag*rk2r(y))/(1.+(rk1*rk2r(y))**2)*phi(1)
      uny(y,a,2)=-rk1*rk2r(y)*uny(y,a,1)
      uny(y,a,3)=-phi(2)
c
c update the temperature
c
      gfac2=sin(bvfsq*dt)*factor
      gfact=cos(bvfsq*dt)*xdif(x,1)*ydif(y,1)*zdif(z,1)
c
      untemp=uy(y,a,4)
      uy(y,a,4)=-gfac2*(oldphi(2)+pnl(2)*2.)
     1              +gfact*(uny(y,a,4)+untemp*2.)
      uny(y,a,4)=-gfac2*(oldphi(2)+pnl(2))
     1              +gfact*(uny(y,a,4)+untemp)
c
 230  continue
c
c      end if
c
c
      end if

#else
 
! viscous integrating factor, including low-Rm MHD
c
	y1=1
	if (x.eq.1.and.z.eq.1) y1=2
c
	if (imaxis.eq.1) then
      do 21 y=y1,ny
         if (.not.mask(y,a)) goto 21
#ifdef LINDBORG
        facty(y,ithr)=xyfac(x,y)*zfac(z)
#else
        facty(y,ithr)=xfac(x)*yfac(y)*zfac(z)
#endif
#ifdef MHD
        rk1=s1*kx(x)
        rk2=s2*ky(y)
        rk3=s3*kz(z)
        rksq=rk1*rk1+rk2*rk2+rk3*rk3
        facty(y,ithr)=facty(y,ithr)*exp(-conduc*rk1*rk1/rksq*dt)
#endif
 21       continue

	else if (imaxis.eq.2) then
      do 22 y=y1,ny
         if (.not.mask(y,a)) goto 22
#ifdef LINDBORG
        facty(y,ithr)=xyfac(x,y)*zfac(z)
#else
        facty(y,ithr)=xfac(x)*yfac(y)*zfac(z)
#endif
#ifdef MHD
        rk1=s1*kx(x)
        rk2=s2*ky(y)
        rk3=s3*kz(z)
        rksq=rk1*rk1+rk2*rk2+rk3*rk3
        facty(y,ithr)=facty(y,ithr)*exp(-conduc*rk2*rk2/rksq*dt)
#endif
 22       continue

	else
      do 23 y=y1,ny
         if (.not.mask(y,a)) goto 23
#ifdef LINDBORG
        facty(y,ithr)=xyfac(x,y)*zfac(z)
#else
        facty(y,ithr)=xfac(x)*yfac(y)*zfac(z)
#endif
#ifdef MHD
        rk1=s1*kx(x)
        rk2=s2*ky(y)
        rk3=s3*kz(z)
        rksq=rk1*rk1+rk2*rk2+rk3*rk3
        facty(y,ithr)=facty(y,ithr)*exp(-conduc*rk3*rk3/rksq*dt)
#endif
 23       continue

	end if


#ifdef LVSTEP
  	if (mod(jstep-1,ivstep).eq.0) then
#endif
#ifndef FROZEN
        do i=1,3
           do 20 y=1,ny
              if(.not.mask(y,a)) goto 20
              untemp=uny(y,a,i)
              uny(y,a,i)=facty(y,ithr)*(untemp+uy(y,a,i))
              uy(y,a,i) =facty(y,ithr)*(untemp+uy(y,a,i)*2.)
 20        continue
        enddo
#endif

#ifdef LVSTEP
	end if
#endif
c
#endif
!
! update the scalars
 
#ifndef NOSCALAR
 
! for buoyancy runs, the first scalar (temperature) is
! already taken care of in the lines above
 
c#ifdef BIF
c        j1=5
c#else
c        j1=4
c#endif
c
	do j=4,3+nc

           do 25 y=1,ny
              if(.not.mask(y,a)) goto 25
           
#ifdef LINDBORG
              factor=xydif(x,y,j-3)*zdif(z,j-3)
#else
              factor=xdif(x,j-3)*ydif(y,j-3)*zdif(z,j-3)
#endif
              untemp=uny(y,a,j)
              uny(y,a,j)=factor*(untemp+uy(y,a,j))
              uy(y,a,j) =factor*(untemp+uy(y,a,j)*2.)
 25        continue
        enddo
#endif
	
	
#ifndef ROTD
! loop 30 carries out the functions of sub. press
 
!  special treatment for zero wave number

#ifdef LVSTEP
 	if (mod(jstep-1,ivstep).ne.0) go to 37
#endif
 
        modpress=1
#ifdef CVC_PRESS
        if (flagc) modpress=0
#endif

        if (modpress.eq.1) then

        yst=1
        if (z.eq.1.and.x.eq.1) then
           uy(1,1,ipsi)=cmplx(0.,0.)
           yst=2
        endif
 
!  loop over x
 
c        do 30 xp=xst,xisz
c	  x=xp+xist-1

        do 30 y=yst,ny
           if(mask(y,a)) then
 
! determine modified pressure from poisson equation
!  solve for psi (which is i*pressure)
 
           bk12=bk1(x,2)*kx(x)+b22(2)*ky(y)**2
         uy(y,a,ipsi)=-(bk1(x,2)*uy(y,a,1)
     1                     +b22(2)*ky(y)*uy(y,a,2)
     1   +bk3(z,2)*uy(y,a,3))/(bk12+bkk3(z,2))
 
!  form true velocities, enforcing continuity
 
        uy(y,a,1)=uy(y,a,1)+kx(x)    *uy(y,a,ipsi)
        uy(y,a,2)=uy(y,a,2)+ky(y)    *uy(y,a,ipsi)
        uy(y,a,3)=uy(y,a,3)+kz(z)    *uy(y,a,ipsi)
      endif
 
30     continue

	end if

#ifdef LVSTEP
 37	continue
#endif
#endif

       call next_xz(xp,z)

 100	continue

!$OMP END PARALLEL
 
! forcing -------------------------------------------------------------
 
#ifdef EPFOR
        if( kforce .ne. 0. ) then

           xp = mystart_x
           z = mystart_z
           do a=1,num_al
              
              x = xp + xist-1


c           do xp=1,kfor
c              x=xp+xist-1

c              do 41 zp=1,zjsz
c                 z=zp+zjst-1

              if( z .le. kfor ) then
                 iz=z 
              else if( z .ge. nz+2-kfor ) then
                 iz=z-nz+k2fo
              else
                 iz=0
              endif
              
              if(iz .ne. 0 .and. x .le. kfor) then

                 do i=1,3
                    do 40 y=1,ny
                       if(.not. mask(y,a)) goto 40
                       
                       if( y .le. kfor ) then
                          iy=y 
                       else if( y .ge. ny+2-kfor ) then
                          iy=y-ny+k2fo
                       else
                          iy=0
                       endif
                       
                       if( iy .ne. 0 ) then

c	  do 40 zp=1,zjsz !kfor
c	  z=zp+zjst-1
c 	  do 40 xp=1,xisz !kfor
c 	  x=xp+xist-1
c 	  if (x.le.kfor .and. z.le.kfor) then  
                          rks=kx(x)**2*b11(2)+ky(y)**2*b22(2)+kz(z)**2*b33(2)
                          if (rks.le.kforce**2) then
                             k=sqrt(rks)+1.5
                             term=real(uy(y,a,i)*conjg(for(x,iy,iz,i)))
        if (istep.eq.1.and.i.eq.1) then
        write (6,"('forced: taskid,x,y,z=',4i6)") taskid,x,y,z
        end if
                             efkz(k,i)=efkz(k,i)+2.*term*tfact(x)*ssq(i)
                             uy(y,a,i)=uy(y,a,i)+for(x,iy,iz,i)
                          endif  
                       endif  
 40                 continue
                 enddo
              endif

              call next_xz(xp,z)
           enddo
       endif
#endif


! ---------------------------------------------------------------------
 


	if (kforce.gt.0.) call forcdt

c      deallocate (factor,bk12)

#ifdef ROTD
c	deallocate (si,snl,cifac,v,oldsi,rk1i)
	if (iraxis.eq.3) then
        deallocate(rk12,rk12r,rk1i)
	else if (iraxis.eq.1) then
        deallocate(rk23,rk23r)
	end if
#elif defined BIF
c	deallocate (phi,oldphi,pnl)
        deallocate(rk12,rk12r,rk2r)
#else
        deallocate(facty)
#endif
 
      return
      end
