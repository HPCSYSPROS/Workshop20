! Version with Cylindrical representation 
!
!################################################
!#      module wavespace_module
!#
!#  contains routines related to operations in wavespace
!#  
!#  NOTE: also contains routines for homogeneous transforms
!#  I.E. the Y inverse and forward transform
!# 
!################################################
      module wavespace_module
      use comp
      implicit none
  

      public  :: wavespace_sub, initialize_wavespace, yitransform_sub
!      private :: 

!############################################
        contains
!###########################################
!#     yitransform_module
!#
!#   transform the y-component 
!#   to physical space
!#
!#   Built on proc1.f
!############################################
          subroutine yitransform_sub(uy,m)
            use comp
            implicit none
            include 'intvars'

#include "fft_stuff.f"

            complex(b8) :: uy(nypad,zjsz*xisz,nu)
            integer :: i,j,m
            complex(b8), allocatable :: shift(:)
            complex(b8) :: syz,sxz
            integer iy,a,xz(2)
            real(b8) factor,rk1,rk1sq

      allocate (shift(nypad))

c

#ifdef ESSL
        call escpft (uy(1,1,1),1,nypad,nypad,num_al,1)
#endif

      xp = mystart_x
      z = mystart_z
      do a=1,num_al

         x = xp + xist-1
         
         sxz=sz(z,kstep)*sx(x,kstep)

#ifdef LVSTEP
	if (mod(jstep-1,ivstep).eq.0.or.
     1      (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) then
#endif
         do 15 y=1,nypad
c      if (y.eq.nyhppad) go to 15
            if(mask(y,a)) then
               shift(y) = sxz *sy(y,kstep)

               uy(y,a,1)=shift(y)*uy(y,a,1)
               uy(y,a,2)=shift(y)*uy(y,a,2)
               uy(y,a,3)=shift(y)*uy(y,a,3)

            else
               uy(y,a,1:3) = 0.0
            endif
 15      continue
#ifdef LVSTEP
	end if
#endif
!
! shift scalars and form y gradients
!
#ifndef NOSCALAR
         do 30 j=4,3+nc
            do 30 y=1,nypad
c              if (y.eq.nyhppad) go to 30
               if(mask(y,a)) then
               shift(y) = sxz *sy(y,kstep)
                  bk2i=imagi*b22(m)*ky(y)
                  
                  uy(y,a,j)    = shift(y) * uy(y,a,j)
                  uy(y,a,j+nc) =  bk2i * uy(y,a,j)
               else
                  uy(y,a,j) = 0.0
                  uy(y,a,j+nc) = 0.0
               endif
 30         continue
#endif

            call next_xz(xp,z)
      enddo

!
c Note by PKY: in cylindrical truncation version the actual
c y-transform from wavenumber space to physical space is
c now absorbed into kxcomm1_trans_cyl (called from itransform.f)
c
! transform in y-direction to physical space
!
!      do 60 i=1,3+2*nc
!
!      do 20 xp=1,xisz
!      do 20 zp=1,zjsz
!      z=zjst+zp-1
!c      if (z.eq.nzhp) go to 20
!	
!      uy(nyhp,a,i)=0.
!
! 20   continue
!
!#ifdef FFTW
!      call fftw_execute(plan2_p1,uy(1,1,1,i),uy(1,1,1,i))
!#else
!      call escfft (uy(1,1,1,i),1,xisz*zjsz,ny,xisz*zjsz,1)
!#endif
!
! 60   continue

	deallocate(shift)
c
        return 
      end subroutine yitransform_sub

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
        subroutine wavespace_sub(uy,m)
          implicit none
!          use comp
          include 'intvars'
 
#include "fft_stuff.f" 
          integer :: i,j,m, is1, i2f,a,xz(2)
          complex(b8) :: uy(ny,zjsz*xisz,nu)
          complex(b8), allocatable :: shift(:)
          complex(b8) :: utmp1
          complex(b8) :: syz,sxz,sxyz
          real(b8) norm,xnorm,ynorm,znorm,vm,vel

        integer ithr,omp_get_thread_num
	integer ii

#if defined (ROTD) || defined (CVC_PRESS)
        real(b8), allocatable :: ksq(:),pij(:),rk2(:),rk2sq(:),rk3(:),rk3sq(:)
        real(b8) rk1,rk1sq
        complex(b8) utmp2
        integer xst,yst
        real s1,s2,s3
        real p11,p12,p13,p22,p23,p33
        real rksq,denom

#ifdef ROTD
        if (kstep.eq.1) then
#elif defined CVC_PRESS
        if (flagc) then
#endif
        allocate (shift(ny))
        allocate (ksq(nxh),pij(nxh),rk2(ny),rk2sq(ny),rk3(nz),rk3sq(nz))
        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
        rk2(:)=s2*ky(:)
        rk2sq(:)=rk2(:)**2
        rk3(:)=s3*kz(:)
        rk3sq(:)=rk3**2
        end if
#endif
        
	i2f=4+nc
#ifdef ROTD
	if (kstep.eq.1) i2f=6+nc
#elif defined CVC_PRESS
        if (flagc) i2f=6+nc
#endif

 
c	do i=1,i2f
c	   call xkcomm2_trans (uy(1,1,1,i),uy(1,1,1,i))
c	end do


!#############################################################
!# transpose to y-pencils and take inverse transform 
 

c	allocate (shift(ny))
 
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

         
#ifdef ROTD
         rk1=s1*kx(x)
         rk1sq = rk1**2
#endif
#ifdef CVC_PRESS
	if (flagc) then
         rk1=s1*kx(x)
         rk1sq = rk1**2
	end if
#endif
 
            sxz=sz(z,kstep)*sx(x,kstep)

#ifdef ROTD
            if (kstep.eq.1) then
               do 40 y=1,ny
                  if( mask(y,a)) then
                     sxyz=-norm*conjg(sxz*sy(y,kstep))
                     bk2i=imagi*s2*ky(y)
                     uy(y,a,1)=sxyz*(uy(y,a,1)+bk2i*uy(y,a,5+nc))
                     uy(y,a,2)=sxyz*(uy(y,a,2)+bk2i*uy(y,a,6+nc))
                     uy(y,a,3)=sxyz*(uy(y,a,3)+bk2i*uy(y,a,4+nc))
                  else
                     uy(y,a,:) = 0.0
                  endif
 40            continue
            else if (kstep.eq.2) then
               do 41 y=1,ny
                  if( mask(y,a)) then
                     sxyz=-norm*conjg(sxz*sy(y,kstep))
                     
                     utmp1=imagi*b22(m)*ky(y)*sxyz
                     uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                     uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                    +bk3(z,m)*uy(y,a,is1))
                     uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
                  else
                     uy(y,a,:) = 0.0
                  endif
 41            continue
            end if
#elif defined CVC_PRESS

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

! form shift, shift, and complete formation of convective terms
 
#ifdef LVSTEP
 	if (mod(jstep-1,ivstep).ne.0) go to 47
#endif
            do 42 y=1,ny
               if( mask(y,a)) then

                  sxyz=-norm*conjg(sxz*sy(y,kstep))
#ifdef LVSTEP
c help get predictor estimate at 'ivstep's forward
	sxyz=sxyz*ivstep
#endif
 
                  utmp1=imagi*b22(m)*ky(y)*sxyz
c
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
#ifdef LVSTEP
 47	continue
#endif
 
	end if
#else

! form shift, shift, and complete formation of convective terms
 
#ifdef LVSTEP
 	if (mod(jstep-1,ivstep).ne.0) go to 47
#endif
            do 42 y=1,ny
               if( mask(y,a)) then

                  sxyz=-norm*conjg(sxz*sy(y,kstep))
#ifdef LVSTEP
c help get predictor estimate at 'ivstep's forward
	sxyz=sxyz*ivstep
#endif
 
                  utmp1=imagi*b22(m)*ky(y)*sxyz
c
                  uy(y,a,1)=sxyz*uy(y,a,1)+utmp1*uy(y,a,2)
                  uy(y,a,2)=sxyz*imagi*(bk1(x,m)*uy(y,a,2)
     1                 +bk3(z,m)*uy(y,a,is1))
                  uy(y,a,3)=sxyz*uy(y,a,3)+utmp1*uy(y,a,is1)
               else
                  uy(y,a,:) = 0.0
               endif
 42         continue
#ifdef LVSTEP
 47	continue
#endif
 

#endif

! shift the (convective terms of) scalars
 
#ifndef NOSCALAR
            do 50 j=4,3+nc
               do 50 y=1,ny
c ?            if (y.eq.nyhp) go to 50
                  if( mask(y,a)) then
                     sxyz=-norm*conjg(sxz*sy(y,kstep))
                     
                     uy(y,a,j)=sxyz*uy(y,a,j)
                  endif
 50            continue
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
        end if

c           endif
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


      return
      end subroutine wavespace_sub
!############################################
!#           Initialize Wavespace
!#
!#     Program that sets up all variables needed
!#     By the wavespace module
!#
!############################################
        subroutine initialize_wavespace
          implicit none
      

        end subroutine initialize_wavespace
!############################################
        end module wavespace_module
