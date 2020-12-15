! Version with Cylindrical representation 
c
c PKY: "proc4c" for 1st substep of RK-4
!
      subroutine advanc_3 (uy,uny,u1y)
#ifdef RKFOUR
#ifndef MAY28
        use comsp
        implicit none
	include 'intvars'
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc)
        complex(b8) :: u1y(ny,zjsz*xisz,3+nc)
        real(b8) bk12 
	integer i,ipsi,xst,iy,iz,ix,m,a,xz(2)
!
        integer ithr,omp_get_thread_num
	integer modpress

          m = 2
 
	modpress=1
#ifdef CVC_PRESS
        if (flagc) modpress=0
#endif

	ipsi=4+nc

!        allocate (bk12(nxh))

	ithr=0
 
!$OMP PARALLEL private (ithr,x,xp,z,bk12)

c Note more may have to be added to private list for ROTD and BIF cases
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)

      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1

         x = xp + xist-1

                    ! update the velocities and scalars
         do i=1,3+nc                          
            do y=1,ny
               if(mask(y,a)) then
                  u1y(y,a,i) = u1y(y,a,i)+uy(y,a,i)*2./3.
                  uy(y,a,i) = uny(y,a,i)+uy(y,a,i)*2.
               endif
            enddo
         enddo


	if (modpress.eq.1) then
c
                    !  loop carries out the functions of sub. press 
                    !  loop over x                          
         do 12 y=1,ny
            if(.not.mask(y,a)) goto 12
            
                                !  special treatment for zero wave number
            if (y.eq.1.and.z.eq.1.and.x.eq.1) then
               uy(1,1,ipsi)=cmplx(0.,0.)
            else                          
                                ! determine modified pressure from poisson equation
                                !  solve for psi (which is i*pressure)
               bk12=bk1(x,2)*kx(x)+b22(2)*ky(y)**2
               
               uy(y,a,ipsi)=-(bk1(x,2)*uy(y,a,1)+ b22(2)*ky(y)*uy(y,a,2)
     1              + bk3(z,2)*uy(y,a,3))/(bk12+bkk3(z,2))
               
                                !
                                !  form true velocities, enforcing continuity
                                !
               uy(y,a,1)=uy(y,a,1)+kx(x)*uy(y,a,ipsi)
               uy(y,a,2)=uy(y,a,2)+ky(y)*uy(y,a,ipsi)
               uy(y,a,3)=uy(y,a,3)+kz(z)*uy(y,a,ipsi)
            endif
 12      continue
c
	end if
         

      call next_xz(xp,z)
!
 100    continue

!$OMP END PARALLEL

              !
              ! forcing -------------------------------------------------------------
              !
#ifdef EPFOR

      if( kforce .ne. 0. .and. mystart_x .le. kfor) then

              
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
                       uy(y,a,i)=uy(y,a,i)+for(x,iy,iz,i)
                    endif
 40              continue
           enddo
        endif

        call next_xz(xp,z)
      enddo      
      endif                     !end of (kforce.ne.0) 


#endif
               !done with epfor idef  

#endif
#endif
        return
      end subroutine advanc_3
