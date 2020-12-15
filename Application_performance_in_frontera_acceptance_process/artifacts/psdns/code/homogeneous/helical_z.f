! Version with Cylindrical representation 
!
      subroutine helical_z (uy,uny,x,z,a,rk12,rk12r) 
#ifdef ROTD
        use comsp        
        implicit none
	include 'intvars'

!new datastructure
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc) 
        real(b8) s1,s2,s3,ssq(3),rks
        real(b8) factor,bk12,rkmag,rkmagr,fact12
        integer i,j,yst,k,m,a

 
	real(b8) rk1,rk2,rk3,rksq,rk13sq,rk3sq

	complex(b8) si(2),snl(2),cifac,v(3),
     1		                    oldsi(2)
        complex(b8) rk1i
	complex(b8) rk2i
        real(b8) rk12(ny),rk12r(ny)
	complex(b8) olduny(3)
	real diff

!
	s1=sqrt(b11(2))
	s2=sqrt(b22(2))
	s3=sqrt(b33(2))
 
         rk1=s1*kx(x)
	rk1i=imagi*rk1
      rk3=s3*kz(z)
      rk3sq=rk3**2
      rk13sq=rk1**2+rk3**2

! use auxiliary variables in diagonized form, in order to integrate
! the rotation term exactly
!
      yst=1
      if (x.eq.1) yst=2
     
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

      si(1)=(rk2i*uny(y,a,1)*s1-rk1i*uny(y,a,2)*s2
     1         -rkmag*uny(y,a,3)*s3)*rk12r(y)
      si(2)=(-rk2i*uny(y,a,1)*s1+rk1i*uny(y,a,2)*s2
     1         -rkmag*uny(y,a,3)*s3)*rk12r(y)
      snl(1)=(rk2i*uy(y,a,1)-rk1i*uy(y,a,2)
     1         -rkmag*uy(y,a,3))*rk12r(y)
      snl(2)=(-rk2i*uy(y,a,1)+rk1i*uy(y,a,2)
     1         -rkmag*uy(y,a,3))*rk12r(y)

!
! update the si's, using the complex-valued integrating factor
! note rrate is twice the rotation rate
!
!     cifac(x)=cexp(-2.*imagi*rk3*rrate*dt*rkmagr)
      cifac=cexp(-imagi*rk3*rrate*dt*rkmagr)
!
      oldsi(1)=si(1)
      oldsi(2)=si(2)
      si(1)=(si(1)+snl(1)*2.)/cifac
      si(2)=(si(2)+snl(2)*2.)*cifac
!
! obtain the v-hats'
      v(3)=-.5*rk12(y)*rkmagr*(si(1)+si(2))
      fact12=rk3*rkmagr*rk12r(y)
      v(1)=-.5*rk2i*rk12r(y)*(si(1)-si(2))+.5*rk1*fact12*(si(1)+si(2))
      v(2)=.5*rk1i*rk12r(y)*(si(1)-si(2))+.5*rk2*fact12*(si(1)+si(2))
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3

! use viscous integrating factor to recover new u at predictor step

      uy(y,a,1)=factor*v(1)
      uy(y,a,2)=factor*v(2)
      uy(y,a,3)=factor*v(3)
!
! values at "half-dt"
!
      si(1)=(oldsi(1)+snl(1))/cifac
      si(2)=(oldsi(2)+snl(2))*cifac
      v(3)=-.5*rk12(y)*rkmagr*(si(1)+si(2))
      v(1)=-.5*rk2i*rk12r(y)*(si(1)-si(2))+.5*rk1*fact12*(si(1)+si(2))
      v(2)=.5*rk1i*rk12r(y)*(si(1)-si(2))+.5*rk2*fact12*(si(1)+si(2))
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3
!

      uny(y,a,1)=factor*v(1)
      uny(y,a,2)=factor*v(2)
      uny(y,a,3)=factor*v(3)
     
#ifdef TEST
      v(3)=-.5*rk12(y)*rkmagr*(oldsi(1)+oldsi(2))
      fact12=rk3*rkmagr*rk12r(y)
      v(1)=-.5*rk2i*rk12r(y)*(oldsi(1)-oldsi(2))+.5*rk1*fact12*(oldsi(1)+oldsi(2))
      v(2)=.5*rk1i*rk12r(y)*(oldsi(1)-oldsi(2))+.5*rk2*fact12*(oldsi(1)+oldsi(2))
      v(1)=v(1)/s1
      v(2)=v(2)/s2
      v(3)=v(3)/s3

	if (rkmag.lt.2.01) then
	v(1)=olduny(1)-v(1)
	v(2)=olduny(2)-v(2)
	v(3)=olduny(3)-v(3)
	diff=amax1(cabs(v(1)),cabs(v(2)),cabs(v(3)))
c	write (6,"('helical_z: x,y,z,diff=',3i4,1p,e12.4)") x,y,z,diff
	write (6,"('helical_z: x,y,z,diff=',3i4,1p,3e12.4)") x,y,z,
     1    cabs(v(1)),cabs(v(2)),cabs(v(3))
	if (x.eq.3.and.y.eq.3.and.z.eq.3) then
	write (6,"('helical_z: rk1,rk2,rk3,rk12,rk12r,rkmag=',1p,6f6.1)") rk1,rk2,rk3,rk12(y),rk12r(y),rkmag
	end if
	end if
#endif

!
 110  continue
!
! special treatment for the kx=ky=0 mode
!
      if (x.eq.1) then
!
      y=1
      si(1)=uny(y,a,1)-imagi*uny(y,a,2)
      si(2)=uny(y,a,1)+imagi*uny(y,a,2)
      snl(1)=uy(y,a,1)-imagi*uy(y,a,2)
      snl(2)=uy(y,a,1)+imagi*uy(y,a,2)
!
      cifac=cexp(-imagi*rrate*dt)
      oldsi(1)=si(1)
      oldsi(2)=si(2)
      si(1)=(si(1)+snl(1)*2.)/cifac
      si(2)=(si(2)+snl(2)*2.)*cifac
!
      v(1)=.5*(si(1)+si(2))
      v(2)=.5*imagi*(si(1)-si(2))
!
! factor(y)=zfac(z) when x=y=1 (kx=ky=0)
!
      uy(y,a,1)=zfac(z)*v(1)
      uy(y,a,2)=zfac(z)*v(2)
!
! value at "half-dt"
!
      si(1)=(oldsi(1)+snl(1))/cifac
      si(2)=(oldsi(2)+snl(2))*cifac
      v(1)=.5*(si(1)+si(2))
      v(2)=.5*imagi*(si(1)-si(2))
      uny(y,a,1)=zfac(z)*v(1)
      uny(y,a,2)=zfac(z)*v(2)
!
! note that u_3 = 0 at this mode, because of continuity
!
      uy(y,a,3)=0.
      uny(y,a,3)=0.
!
      end if

 
#endif
      return
      end
