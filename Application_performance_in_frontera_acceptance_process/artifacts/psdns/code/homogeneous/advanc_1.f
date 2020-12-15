! Version with Cylindrical representation 
!
c PKY: extension of proc4a for 1st substep of RK-4
c ROTD, BIF and LINDBORG directives to be added back later
c
      subroutine advanc_1 (uy,uny,u1y) 
#ifdef RKFOUR
        use comsp        
        implicit none
	include 'intvars'
!new datastructure
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc) 
        complex(b8) :: u1y(ny,zjsz*xisz,3+nc) 
        complex(b8) untemp
        complex caux
        real(b8) s1,s2,s3,ssq(3),rks
        real(b8) term
        real(b8) factor,bk12
        integer i,j,ipsi,yst,k,iy,iz,j1,m,a,xz(2),ymax,almax
!
        integer ithr,omp_get_thread_num
	integer modpress
!
        real(b8), allocatable :: facty(:,:)


        m = 2

        modpress=1
#ifdef CVC_PRESS
        if (flagc) modpress=0
#endif

 
      ipsi=4+nc
!
        allocate(facty(ny,0:num_thr-1))
!
#ifdef EPFOR
      if (kforce.gt.0.) efkz(:,:)=0.
#endif
 
      s1=sqrt(b11(2))
      s2=sqrt(b22(2))
      s3=sqrt(b33(2))
      ssq(1)=b11(2)
      ssq(2)=b22(2)
      ssq(3)=b33(2)

	ithr=0

!$OMP PARALLEL private (ithr,x,xp,z,factor,untemp,yst,bk12)

c Note more may have to be added to private list for ROTD and BIF cases
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)

      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
         x = xp + xist-1



 
! define product of mask and integrating factors
 
! update the velocities
!
! viscous integrating factor
 
      do 21 y=1,ny
         if(.not.mask(y,a)) goto 21

        facty(y,ithr)=xfac(x)*yfac(y)*zfac(z)

 21       continue

c
c note that the nonlinear term is computed with a
c factor of dt/2 attached  to it
c
        do i=1,3
           do 20 y=1,ny
              if(.not.mask(y,a)) goto 20
              untemp=uny(y,a,i)
c
              u1y(y,a,i)=facty(y,ithr)*(uny(y,a,i)+uy(y,a,i)/3.)
              uy(y,a,i) =facty(y,ithr)*(uny(y,a,i)+uy(y,a,i))
c
 20        continue
        enddo
!
! update the scalars
 
#ifndef NOSCALAR
 
! for buoyancy runs, the first scalar (temperature) is
! already taken care of in the lines above
 
	do j=4,3+nc

           do 25 y=1,ny
              if(.not.mask(y,a)) goto 25
           
              factor=xdif(x,j-3)*ydif(y,j-3)*zdif(z,j-3)
              u1y(y,a,j)=factor*(uny(y,a,j)+uy(y,a,j)/3.)
              uy(y,a,j) =factor*(uny(y,a,j)+uy(y,a,j))
 25        continue
        enddo
#endif
	
	
	if (modpress.eq.1) then
c
! loop 30 carries out the functions of sub. press
 
!  special treatment for zero wave number
 
        yst=1
        if (z.eq.1.and.x.eq.1) then
           uy(1,1,ipsi)=cmplx(0.,0.)
           yst=2
        endif
 
!  loop over x
 

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
c
       call next_xz(xp,z)
!
 100	continue

!$OMP END PARALLEL
 
 
! forcing -------------------------------------------------------------
 
#ifdef EPFOR
        if( kforce .ne. 0. ) then

           xp = mystart_x
           z = mystart_z
           do a=1,num_al
              
              x = xp + xist-1


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

c 	  if (x.le.kfor .and. z.le.kfor) then  
                          rks=kx(x)**2*b11(2)+ky(y)**2*b22(2)+kz(z)**2*b33(2)
                          if (rks.le.kforce**2) then
                             k=sqrt(rks)+1.5
                             term=real(uy(y,a,i)*conjg(for(x,iy,iz,i)))
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
 


c	if (kforce.gt.0.) call forcdt


        deallocate(facty)
 

#endif
      return
      end
