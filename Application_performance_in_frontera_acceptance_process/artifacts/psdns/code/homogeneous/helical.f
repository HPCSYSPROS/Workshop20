! Version with Cylindrical representation 
!
      subroutine helical (uy,uny,x,z,a) 

! Fixed on 8/21/2015: the 'snl' terms should not be multiplied
! by the metric factors

! 8/22/2015: this version handles all 3 directions of rotation

#ifdef ROTD
        use comsp        
        implicit none
	include 'intvars'

!new datastructure
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc) 
        real(b8) s1,s2,s3,ssq(3)
        real(b8) factor,rkmag,rkmagr,fact23
        integer i,j,yst,k,m,a,ylast,isp

 
	real(b8) rk1,rk2,rk3,rk13sq,rk3sq

	complex(b8) si(2),snl(2),cifac,v(3),oldsi(2)
        complex(b8) rk3i
	complex(b8) rk2i
        real(b8) rk23,vfac
        real(b8) e1(3),e2(3),e3(3),nvec(3),emag
        complex(b8) uhat(3),nkp(3),nkm(3)

!
	s1=sqrt(b11(2))
	s2=sqrt(b22(2))
	s3=sqrt(b33(2))

c note 'iraxis' (1/2/3) is known from input.f
c
         rk1=s1*kx(x)
      rk3=s3*kz(z)
	rk3i=imagi*rk3
      rk3sq=rk3**2
      rk13sq=rk1**2+rk3**2

! use auxiliary variables in diagonized form, in order to integrate
! the rotation term exactly
!
      yst=1
      if (iraxis.eq.1.and.z.eq.1) yst=2
      if (iraxis.eq.3.and.x.eq.1) yst=2
      if (iraxis.eq.2.and.x.eq.1.and.z.eq.1) yst=ny+1

        nvec=0.
        nvec(iraxis)=1.

     
      do 110 y=yst,ny
         if(.not.mask(y,a)) goto 110

#ifdef LINDBORG
        factor=xyfac(x,y)*zfac(z)
#else
        factor=xfac(x)*yfac(y)*zfac(z)
#endif


!
! define the two auxiliary variables (Mansour etal 1994) and their
! corresponding nonlinear terms (projected in the direction normal to k)
!
         rk2=s2*ky(y)
         rk2i=imagi*rk2
         rkmag=sqrt(rk13sq+rk2**2)
         rkmagr=1./rkmag
	rk23=sqrt(rk2*rk2+rk3*rk3)

#ifdef MHD
	if (imaxis.eq.1) then
	factor=factor*exp(-conduc*rk1*rk1/rkmag**2*dt)
	else if (imaxis.eq.2) then
	factor=factor*exp(-conduc*rk2*rk2/rkmag**2*dt)
	else if (imaxis.eq.3) then
	factor=factor*exp(-conduc*rk3*rk3/rkmag**2*dt)
	end if
#endif

c define the new unit vectors: e1 = k/|k| x n, e2 = k/|k| x e1

        e3(1)=rk1/rkmag
        e3(2)=rk2/rkmag
        e3(3)=rk3/rkmag
c
        e1(1)=e3(2)*nvec(3)-e3(3)*nvec(2)
        e1(2)=e3(3)*nvec(1)-e3(1)*nvec(3)
        e1(3)=e3(1)*nvec(2)-e3(2)*nvec(1)
        emag=sqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
        e1(:)=e1(:)/emag
        e2(1)=e3(2)*e1(3)-e3(3)*e1(2)
        e2(2)=e3(3)*e1(1)-e3(1)*e1(3)
        e2(3)=e3(1)*e1(2)-e3(2)*e1(1)
        emag=sqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
        e2(:)=e2(:)/emag
c
c define N+(k) and N-(k)

        do i=1,3
        nkp(i)=e2(i)-imagi*e1(i)
        nkm(i)=e2(i)+imagi*e1(i)
        end do


        uhat(1)=uny(y,a,1)*s1
        uhat(2)=uny(y,a,2)*s2
        uhat(3)=uny(y,a,3)*s3
        si(1)=uhat(1)*nkm(1)+uhat(2)*nkm(2)+uhat(3)*nkm(3)
        si(2)=uhat(1)*nkp(1)+uhat(2)*nkp(2)+uhat(3)*nkp(3)
        uhat(1)=uy(y,a,1)
        uhat(2)=uy(y,a,2)
        uhat(3)=uy(y,a,3)
        snl(1)=uhat(1)*nkm(1)+uhat(2)*nkm(2)+uhat(3)*nkm(3)
        snl(2)=uhat(1)*nkp(1)+uhat(2)*nkp(2)+uhat(3)*nkp(3)


!
! update the si's, using the complex-valued integrating factor
! note rrate is twice the rotation rate
!
!     cifac(x)=cexp(-2.*imagi*rk3*rrate*dt*rkmagr)
!      cifac=cexp(-imagi*rk1*rrate*dt*rkmagr)
      cifac=cexp(-imagi*e3(iraxis)*rrate*dt)
!
      oldsi(1)=si(1)
      oldsi(2)=si(2)
      si(1)=(si(1)+snl(1)*2.)/cifac
      si(2)=(si(2)+snl(2)*2.)*cifac
!
! obtain the v-hats'
      

        do i=1,3
	v(i)=.5*(si(1)*nkp(i)+si(2)*nkm(i))
        end do
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3
      
!
! use viscous integrating factor to recover new u at predictor step

      uy(y,a,1)=factor*v(1)
      uy(y,a,2)=factor*v(2)
      uy(y,a,3)=factor*v(3)
!
! values at "half-dt"
!
      si(1)=(oldsi(1)+snl(1))/cifac
      si(2)=(oldsi(2)+snl(2))*cifac
        do i=1,3
	v(i)=.5*(si(1)*nkp(i)+si(2)*nkm(i))
	end do
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3
!

      uny(y,a,1)=factor*v(1)
      uny(y,a,2)=factor*v(2)
      uny(y,a,3)=factor*v(3)
     
!
 110  continue
!
! special treatment for the mode of zero wavenumber 
! perpendicular to the axis of rotation
!
        isp=0
        if (iraxis.eq.1.and.z.eq.1) then
        i=2 ; j=3 ; ylast=1
        isp=1
        vfac=xfac(x)
        end if
        if (iraxis.eq.3.and.x.eq.1) then
        i=1 ; j=2 ; ylast=1
        isp=3
        vfac=zfac(z)
        end if
        if (iraxis.eq.2.and.x.eq.1.and.z.eq.1) then
        i=3 ; j=1 ; ylast=ny
        isp=2
        end if

        if (isp.gt.0) then
        do 120 y=1,ylast

      si(1)=uny(y,a,i)-imagi*uny(y,a,j)
      si(2)=uny(y,a,i)+imagi*uny(y,a,j)
      snl(1)=uy(y,a,i)-imagi*uy(y,a,j)
      snl(2)=uy(y,a,i)+imagi*uy(y,a,j)
!
      cifac=cexp(-imagi*rrate*dt)
      oldsi(1)=si(1)
      oldsi(2)=si(2)
      si(1)=(si(1)+snl(1)*2.)/cifac
      si(2)=(si(2)+snl(2)*2.)*cifac
!
      v(i)=.5*(si(1)+si(2))
      v(j)=.5*imagi*(si(1)-si(2))
!
! factor(y)=xfac(x) when y=z=1 (ky=kz=0)
!
	if (isp.eq.2) vfac=yfac(y)
      uy(y,a,i)=vfac*v(i)
      uy(y,a,j)=vfac*v(j)
!
! value at "half-dt"
!
      si(1)=(oldsi(1)+snl(1))/cifac
      si(2)=(oldsi(2)+snl(2))*cifac
      v(i)=.5*(si(1)+si(2))
      v(j)=.5*imagi*(si(1)-si(2))
      uny(y,a,i)=vfac*v(i)
      uny(y,a,j)=vfac*v(j)
!
! note that u_1 = 0 at this mode, because of continuity
!
      uy(y,a,iraxis)=0.
      uny(y,a,iraxis)=0.
!
 120    continue

      end if

 
#endif
      return
      end
