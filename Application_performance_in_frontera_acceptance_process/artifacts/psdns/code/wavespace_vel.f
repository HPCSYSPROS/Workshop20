!############################################
!#           Wavespace
!#
!#     Program that operates when in wavespace
!#     Therefore, for homogeneous turbulence
!#     Program will transform y to finish forming 
!#     Nonlinear terms              
!#
!#     Built on proc3.f
!############################################
        subroutine wavespace_vel(uy,m)
           use comp
          implicit none
          include 'intvars'
 
          integer :: i,j,m, is1, i2f,a,xz(2)

          complex(b8) :: uy(ny,zjsz*xisz,nu)
          complex(b8), allocatable :: shift(:)

          complex(b8) :: utmp1
          complex(b8) :: syz,sxz,sxyz
	real(b8) norm

        integer ithr,omp_get_thread_num
	real(b8) upy_factor
	integer ii
	
	
        real(b8) rk1,rk1sq
        real s1,s2,s3

#if defined (ROTD) || defined (CVC_PRESS)
        real(b8), allocatable :: ksq(:),pij(:),rk2(:),rk2sq(:),rk3(:),rk3sq(:)
        complex(b8) utmp2
        integer xst,yst
        real p11,p12,p13,p22,p23,p33
        real rksq,denom
#endif

        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
c
#if defined (ROTD) || defined (CVC_PRESS)
#ifdef ROTD
	if (kstep.eq.1) then
#elif defined CVC_PRESS
        if (flagc) then
#endif
        allocate (shift(ny))
        allocate (ksq(nxh),pij(nxh),rk2(ny),rk2sq(ny),rk3(nz),rk3sq(nz))
        rk2(:)=s2*ky(:)
        rk2sq(:)=rk2(:)**2
        rk3(:)=s3*kz(:)
        rk3sq(:)=rk3**2
        end if
#endif


	i2f=4+nc
	upy_factor=.5

#ifdef CVC_PRESS
        if (flagc) i2f=6+nc
	if (rkmethod.eq.4)then
	if (kstep.eq.1) upy_factor=1./6.
	if (kstep.eq.2) upy_factor=1./3.
	if (kstep.eq.3) upy_factor=1./3.
	if (kstep.eq.4) upy_factor=1./6.
	end if
#endif


!#############################################################
!# transpose to y-pencils and take inverse transform 
 

 
      norm=dt2/ny/nz
 
      is1=4+nc

        ithr=0
!$OMP PARALLEL private (ithr,x,xp,z,sxz,sxyz,utmp1)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1


         x = xp + xist-1

         rk1=s1*kx(x)
         rk1sq = rk1**2

            sxz=sz(z,kstep)*sx(x,kstep)


! form shift, shift, and complete formation of convective terms

#ifdef ROTD

	if (kstep.eq.1) then
               do 40 y=1,ny
                  if( mask(y,a)) then
                     shift(y)=-norm*conjg(sxz*sy(y,kstep))
                     bk2i=imagi*s2*ky(y)
                     uy(y,a,1)=shift(y)*(uy(y,a,1)+bk2i*uy(y,a,5+nc))
                     uy(y,a,2)=shift(y)*(uy(y,a,2)+bk2i*uy(y,a,6+nc))
                     uy(y,a,3)=shift(y)*(uy(y,a,3)+bk2i*uy(y,a,4+nc))
                  else
                     uy(y,a,:) = 0.0
                  endif
 40            continue
	else if (kstep.eq.2) then
            do 42 y=1,ny
               if( mask(y,a)) then
                  sxyz=-norm*conjg(sxz*sy(y,kstep))
                  utmp1=imagi*b22(m)*ky(y)*sxyz
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
	end if
#elif defined  CVC_PRESS	
	if (flagc) then
               do 40 y=1,ny
                  if( mask(y,a)) then
                     shift(y)=-norm*conjg(sxz*sy(y,kstep))
                     bk2i=imagi*s2*ky(y)
                     uy(y,a,1)=shift(y)*(uy(y,a,1)+bk2i*uy(y,a,5+nc))
                     uy(y,a,2)=shift(y)*(uy(y,a,2)+bk2i*uy(y,a,6+nc))
                     uy(y,a,3)=shift(y)*(uy(y,a,3)+bk2i*uy(y,a,4+nc))
                  else
                     uy(y,a,:) = 0.0
                  endif
 40            continue
	else
            do 42 y=1,ny
               if( mask(y,a)) then
                  sxyz=-norm*conjg(sxz*sy(y,kstep))
                  utmp1=imagi*b22(m)*ky(y)*sxyz
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
	end if
#else
            do 42 y=1,ny
               if( mask(y,a)) then
                  sxyz=-norm*conjg(sxz*sy(y,kstep))
                  utmp1=imagi*b22(m)*ky(y)*sxyz
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
#endif
 
#if defined (ROTD) || defined (CVC_PRESS)
c project nonlinear terms in the direction perpendicular to wavenumber vector
c
#ifdef ROTD
	if (kstep.eq.1) then
#elif defined CVC_PRESS
        if (flagc) then
#endif

         yst=1
         if (x.eq.1.and.z.eq.1) yst=2
         do 220 y=yst,ny
            if( mask(y,a)) then

               rksq=rk2sq(y)+rk1sq+rk3sq(z)
               denom=1./rksq
               p11=1.-rk1sq*denom
               p22=1.-rk2sq(y)*denom
               p33=1.-rk3sq(z)*denom
               p12=-rk1*rk2(y)*denom
               p13=-rk1*rk3(z)*denom
               p23=-rk2(y)*rk3(z)*denom
	
#ifdef CVC_PRESS
        utmp1=s1*kx(x)*uy(y,a,1)+s2*ky(y)*uy(y,a,2)+s3*kz(z)*uy(y,a,3)
        upy(y,a)=upy(y,a)-upy_factor*imagi*utmp1*denom/dt2
#endif
               utmp1=p11*uy(y,a,1)+p12*uy(y,a,2)+p13*uy(y,a,3)
               utmp2=p12*uy(y,a,1)+p22*uy(y,a,2)+p23*uy(y,a,3)
               uy(y,a,3)=p13*uy(y,a,1)+p23*uy(y,a,2)
     1              +p33*uy(y,a,3)
               uy(y,a,1)=utmp1
               uy(y,a,2)=utmp2

            endif
 220     continue
        end if

#endif

 

      call next_xz(xp,z)

 100	continue

!$OMP END PARALLEL

#if defined (ROTD) || defined (CVC_PRESS)
	if (allocated(shift)) then
        deallocate (ksq,pij,rk2,rk2sq,rk3,rk3sq)
        deallocate (shift)
        end if
#endif

c

      return
	end subroutine wavespace_vel
