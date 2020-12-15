! Version with Cylindrical representation 
!
      subroutine proc4b (uy,uny)
        use comsp
        implicit none
	include 'intvars'
        complex(b8) :: uy(ny,zjsz*xisz,nu)
        complex(b8) :: uny(ny,zjsz*xisz,3+nc)
        real(b8) bk12 
	integer i,ipsi,xst,iy,iz,ix,m,a,xz(2)
	  real xnorm,ynorm,znorm,tmpvmaxz,vel,vm

        integer ithr,omp_get_thread_num
	integer modpress


          m = 2
 
	ithr=0
 
        ipsi=4+nc


!$OMP PARALLEL private (ithr,x,xp,z,bk12)
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 100 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1

         x = xp + xist-1

                    ! update the velocities and scalars
#ifdef LVSTEP
c
 	if (mod(jstep-1,ivstep).eq.0) then
         do i=1,3                          
            do y=1,ny
               if(mask(y,a)) then
                  uy(y,a,i) = uny(y,a,i)+uy(y,a,i)
               endif
            enddo
         enddo
	else
 	uy(:,a,1)=uny(:,a,1)
 	uy(:,a,2)=uny(:,a,2)
 	uy(:,a,3)=uny(:,a,3)
	end if
#else
         do i=1,3                          
            do y=1,ny
#ifndef FROZEN
               if(mask(y,a)) then
                  uy(y,a,i) = uny(y,a,i)+uy(y,a,i)
               endif
#else
               if(mask(y,a)) then
                  uy(y,a,i) = uny(y,a,i))
               endif
#endif
            enddo
         enddo
#endif

         do i=4,3+nc                          
            do y=1,ny
               if(mask(y,a)) then
                  uy(y,a,i) = uny(y,a,i)+uy(y,a,i)
               endif
            enddo
         enddo
                    !  loop carries out the functions of sub. press 
                    !  loop over x                          

	modpress=1
#ifdef CVC_PRESS
 	if (flagc) modpress=0
#endif

	if (modpress.eq.1) then
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
	end if
         
      call next_xz(xp,z)
c
 100	continue

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

              !
              ! finished forcing -------------------------------------------------------------
              !

#ifdef LINDBORG
        ! multiple dt by energy injection rate if not set to unity
        call force_lindborg (uy,2)
#endif

c
! copy new u to uny

	ithr=0
!$OMP PARALLEL private (ithr)
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
	do i=1,3+nc
      do a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
              do y=1,ny
                 if (mask(y,a))then                       
                    uny(y,a,i)=uy(y,a,i)
                 endif
              end do
           end do              
           end do              
!$OMP END PARALLEL
	

c        deallocate (bk12)

        !
        ! enforce zero mean
        !
        if (taskid .eq. 0) then
           do i=1,3+nc
              uy(1,1,i)=0.
              uny(1,1,i)=0.
           enddo
        end if

	 return
      end subroutine proc4b
