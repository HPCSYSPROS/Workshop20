! Version with Cylindrical representation 
!
      subroutine helical_x (uy,uny,x,z,a,rk23,rk23r) 

! Fixed on 8/21/2015: the 'snl' terms should not be multiplied
! by the metric factors

#ifdef ROTD
        use comsp        
        implicit none
	include 'intvars'

!new datastructure
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc) 
        real(b8) s1,s2,s3,ssq(3)
        real(b8) factor,rkmag,rkmagr,fact23
        integer i,j,yst,k,m,a

 
	real(b8) rk1,rk2,rk3,rk13sq,rk3sq

	complex(b8) si(2),snl(2),cifac,v(3),
     1		                    oldsi(2)
        complex(b8) rk3i
	complex(b8) rk2i
	complex(b8) olduny(3)
        real(b8) rk23(ny),rk23r(ny),diff

!
	s1=sqrt(b11(2))
	s2=sqrt(b22(2))
	s3=sqrt(b33(2))

c case of uniform rotation along the x-axis
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
      if (z.eq.1) yst=2
     
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

#ifdef MHD
	if (imaxis.eq.1) then
	factor=factor*exp(-conduc*rk1*rk1/rkmag**2*dt)
	else if (imaxis.eq.2) then
	factor=factor*exp(-conduc*rk2*rk2/rkmag**2*dt)
	else if (imaxis.eq.3) then
	factor=factor*exp(-conduc*rk3*rk3/rkmag**2*dt)
	end if
#endif
#ifdef TEST
        olduny(1)=uny(y,a,1)
        olduny(2)=uny(y,a,2)
        olduny(3)=uny(y,a,3)
#endif


	si(1)=(-rkmag*uny(y,a,1)*s1+rk3i*uny(y,a,2)*s2
     1         -rk2i*uny(y,a,3)*s3)*rk23r(y)
	si(2)=(-rkmag*uny(y,a,1)*s1-rk3i*uny(y,a,2)*s2
     1         +rk2i*uny(y,a,3)*s3)*rk23r(y)
	snl(1)=(-rkmag*uy(y,a,1)+rk3i*uy(y,a,2)
     1         -rk2i*uy(y,a,3))*rk23r(y)
	snl(2)=(-rkmag*uy(y,a,1)-rk3i*uy(y,a,2)
     1         +rk2i*uy(y,a,3))*rk23r(y)

!
! update the si's, using the complex-valued integrating factor
! note rrate is twice the rotation rate
!
!     cifac(x)=cexp(-2.*imagi*rk3*rrate*dt*rkmagr)
      cifac=cexp(-imagi*rk1*rrate*dt*rkmagr)
!
      oldsi(1)=si(1)
      oldsi(2)=si(2)
      si(1)=(si(1)+snl(1)*2.)/cifac
      si(2)=(si(2)+snl(2)*2.)*cifac
!
! obtain the v-hats'
      
      fact23=rk1*rkmagr*rk23r(y)
      v(1)=-.5*rk23(y)*rkmagr*(si(1)+si(2))
      v(2)=-.5*rk3i*rk23r(y)*(si(1)-si(2))+.5*rk2*fact23*(si(1)+si(2))
      v(3)=.5*rk2i*rk23r(y)*(si(1)-si(2))+.5*rk3*fact23*(si(1)+si(2))
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
      v(1)=-.5*rk23(y)*rkmagr*(si(1)+si(2))
      v(2)=-.5*rk3i*rk23r(y)*(si(1)-si(2))+.5*rk2*fact23*(si(1)+si(2))
      v(3)=.5*rk2i*rk23r(y)*(si(1)-si(2))+.5*rk3*fact23*(si(1)+si(2))
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3
!

      uny(y,a,1)=factor*v(1)
      uny(y,a,2)=factor*v(2)
      uny(y,a,3)=factor*v(3)
     
!
#ifdef TEST
      v(1)=-.5*rk23(y)*rkmagr*(oldsi(1)+oldsi(2))
      v(2)=-.5*rk3i*rk23r(y)*(oldsi(1)-oldsi(2))+.5*rk2*fact23*(oldsi(1)+oldsi(2))
      v(3)=.5*rk2i*rk23r(y)*(oldsi(1)-oldsi(2))+.5*rk3*fact23*(oldsi(1)+oldsi(2))
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3

        if (rkmag.lt.2.01) then
        v(1)=olduny(1)-v(1)
        v(2)=olduny(2)-v(2)
        v(3)=olduny(3)-v(3)
        diff=amax1(cabs(v(1)),cabs(v(2)),cabs(v(3)))
        write (6,"('helical_x: x,y,z,diff=',3i4,1p,3e12.4)") x,y,z,
     1    cabs(v(1)),cabs(v(2)),cabs(v(3))
        end if
#endif

!
 110  continue
!
! special treatment for the ky=kz=0 mode
!
      if (z.eq.1) then
!
      y=1
      si(1)=uny(y,a,2)-imagi*uny(y,a,3)
      si(2)=uny(y,a,2)+imagi*uny(y,a,3)
      snl(1)=uy(y,a,2)-imagi*uy(y,a,3)
      snl(2)=uy(y,a,2)+imagi*uy(y,a,3)
!
      cifac=cexp(-imagi*rrate*dt)
      oldsi(1)=si(1)
      oldsi(2)=si(2)
      si(1)=(si(1)+snl(1)*2.)/cifac
      si(2)=(si(2)+snl(2)*2.)*cifac
!
      v(2)=.5*(si(1)+si(2))
      v(3)=.5*imagi*(si(1)-si(2))
!
! factor(y)=xfac(x) when y=z=1 (ky=kz=0)
!
      uy(y,a,2)=xfac(x)*v(2)
      uy(y,a,3)=xfac(x)*v(3)
!
! value at "half-dt"
!
      si(1)=(oldsi(1)+snl(1))/cifac
      si(2)=(oldsi(2)+snl(2))*cifac
      v(2)=.5*(si(1)+si(2))
      v(3)=.5*imagi*(si(1)-si(2))
      uny(y,a,2)=xfac(x)*v(2)
      uny(y,a,3)=xfac(x)*v(3)
!
! note that u_1 = 0 at this mode, because of continuity
!
      uy(y,a,1)=0.
      uny(y,a,1)=0.
!
      end if

 
#endif
      return
      end
