!#######################################################
!#     Nonlinear
!#
!#     Program to handle all nonlinear terms
!#
!#######################################################
	subroutine realspace(ux,m,rkstep)
          use comp
          implicit none
          include 'intvars'
#include "fft_stuff.f"
!	real(b8) factor
 	real(b8)    :: ux(nxpad,zisz,yjsz,nu)

        real(b8), allocatable :: bv2(:,:)
	integer :: m,i, is1,is2,novt,rkstep,icall
        real(b8) s1,s2,s3,xnorm,ynorm,znorm,vel,k3z,bk3z
	real tmpu,tmpv,tmpw
c
        integer ithr,OMP_GET_THREAD_NUM,yp1,yp2
	real(b8), allocatable :: vmax(:)
	real sum,sum0

	integer ii
c
	kstep=rkstep

          !1st Step
          if(rkstep.eq.1)then
#ifndef NOSCALAR  
             if (taskid.eq.0.and.schflag)
     1       write (6,*) 'scalar PDFs generated at step=',istep
!
! extract PDFs
!
             do i=1,nc
                if (schflag) then
                   call phymom (u,i+3+nc,i+3+nc)
                   call scnpdf (u(1,1,1,i+3+nc),i)
                end if
             end do 
#endif

#ifdef LVSTEP
        if (mod(jstep-1,ivstep).eq.0.or.
     1    (mod(jstep-2,ivstep).eq.0.and.kstep.eq.1)) go to 100
        return
  100  continue
#endif

#ifdef LAG

c Lagrangian particle tracking call
c
        call lagpart
c
#endif

          endif !done with first step only execution

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
      allocate (bv2(nxpad,0:num_thr-1))
!
	if (rkstep.eq.1) allocate (vmax(0:num_thr-1))

c
	call divide_thr (yjsz,iyp1,iypi,'yjsz: realspace')
c
! Courant number and convective terms in physical space
!
      if (rkstep.eq.1) then

	ithr=0
!$OMP PARALLEL private (ithr,vel,velmax,tmpu,tmpv,tmpw)
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
	velmax=0.
c
c     do 120 yp=1,yjsz
      do 120 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
      do 120 zp=1,zisz

      do 125 x=1,nxpad
!
!
#if ROTD
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
!cc to test 3/22/07: try replacing s1,s2,s3 by 1?
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
         bv2(x,ithr)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x,ithr)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x,ithr)
	end if
#else
         bv2(x,ithr)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x,ithr)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x,ithr)
#endif

 125   continue
 120   continue
	vmax(ithr)=velmax
!$OMP END PARALLEL
c
	velmax=0.
	do ithr=0,num_thr-1
	velmax=max(velmax,vmax(ithr))
	end do
c
	deallocate (vmax)
c  	print *,taskid,': realspace velmax=',velmax
c
#ifdef RKFOUR
      else if (kstep.ge.2) then
#else
      else if (kstep.eq.2) then
#endif
!
c
	ithr=0
!$OMP PARALLEL private (ithr)
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
c
c     do 130 yp=1,yjsz
      do 130 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
      do 130 zp=1,zisz
!
      do 135 x=1,nxpad

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
         bv2(x,ithr)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x,ithr)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x,ithr)
	end if
#else
!
         bv2(x,ithr)=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
         ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
         ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
!
         ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2(x,ithr)
         ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
         ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2(x,ithr)
#endif
!
 135   continue
 130   continue
!$OMP END PARALLEL
!
      end if

!################
!#                   Done with realspace, now fft forward
!#################
      deallocate (bv2)
c
	return
	end
