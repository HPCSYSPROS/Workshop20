!#######################################################
!#     Nonlinear
!#
!#     Program to handle all nonlinear terms
!#
!#######################################################
	subroutine realspace_vel (ux,m,rkstep)
          use comp
          implicit none
          include 'intvars'
#include "fft_stuff.f"
 	real(b8)    :: ux(nxpad,zisz,yjsz,nu)

!	real(b8) bv2
c
	integer :: m,i, is1,is2,rkstep,icall
        real(b8) s1,s2,s3,xnorm,ynorm,znorm,vel
	real tmpu,tmpv,tmpw
c
        integer ithr,OMP_GET_THREAD_NUM,yp1,yp2
	real(b8), allocatable :: vmax(:),bv2(:)

      is1=4+nc
      is2=5+nc
!
      s1=sqrt(b11(m))
      s2=sqrt(b22(m))
      s3=sqrt(b33(m))
!
      xnorm=b11(m)*nxpad
      ynorm=b22(m)*nypad
      znorm=b33(m)*nzpad
!
#ifdef LAG

c Lagrangian particle tracking call
c
        call lagpart
c
#endif

	allocate (vmax(0:num_thr-1))
c
	call divide_thr (yjsz,iyp1,iypi,'yjsz: realspace')
c
! Courant number and convective terms in physical space
!
      if (rkstep.eq.1) then
	allocate (bv2(nxpad))

	ithr=0
!$OMP PARALLEL private (ithr,vel,velmax,tmpu,tmpv,tmpw,bv2)
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
	velmax=0.
c
	bv2=0.
c
      do 120 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
      do 120 zp=1,zisz

      do 125 x=1,nxpad
!
#if ROTD
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
      tmpu=s1*ux(x,zp,yp,1)
      tmpv=s2*ux(x,zp,yp,2)
      tmpw=s3*ux(x,zp,yp,3)
!        ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
      ux(x,zp,yp,1)=tmpu**2
      ux(x,zp,yp,2)=tmpu*tmpv
      ux(x,zp,yp,3)=tmpw**2
      ux(x,zp,yp,4+nc)=tmpv*tmpw
      ux(x,zp,yp,5+nc)=tmpu*tmpw
      ux(x,zp,yp,6+nc)=tmpv**2
#elif defined CVC_PRESS
        if (flagc) then
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
      tmpu=s1*ux(x,zp,yp,1)
      tmpv=s2*ux(x,zp,yp,2)
      tmpw=s3*ux(x,zp,yp,3)
      ux(x,zp,yp,1)=tmpu**2
      ux(x,zp,yp,2)=tmpu*tmpv
      ux(x,zp,yp,3)=tmpw**2
      ux(x,zp,yp,4+nc)=tmpv*tmpw
      ux(x,zp,yp,5+nc)=tmpu*tmpw
      ux(x,zp,yp,6+nc)=tmpv**2
        else
         bv2(x)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x)
        end if
#else
         bv2(x)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x)
#endif

 125   continue
 120   continue
	vmax(ithr)=velmax
!

!$OMP END PARALLEL
	deallocate (bv2)
c
	velmax=0.
	do ithr=0,num_thr-1
	velmax=max(velmax,vmax(ithr))
	end do
c
	deallocate (vmax)
c
      else
!
c
	ithr=0
!$OMP PARALLEL private (ithr,bv2)
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
c
	allocate (bv2(nxpad))
	bv2=0.
c
      do 130 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
      do 130 zp=1,zisz
!
      do 135 x=1,nxpad
!
#ifdef CVC_PRESS
	if (flagc) then
      tmpu=s1*ux(x,zp,yp,1)
      tmpv=s2*ux(x,zp,yp,2)
      tmpw=s3*ux(x,zp,yp,3)
      ux(x,zp,yp,1)=tmpu**2
      ux(x,zp,yp,2)=tmpu*tmpv
      ux(x,zp,yp,3)=tmpw**2
      ux(x,zp,yp,4+nc)=tmpv*tmpw
      ux(x,zp,yp,5+nc)=tmpu*tmpw
      ux(x,zp,yp,6+nc)=tmpv**2
	else
         bv2(x)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x)
	end if
#else
         bv2(x)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x)
#endif
!
 135   continue
 130   continue
	deallocate (bv2)
!$OMP END PARALLEL
!
      end if

!################
!#                   Done with realspace, now fft forward
!#################
	return
	end
