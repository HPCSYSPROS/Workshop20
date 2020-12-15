!################################################
!#      module stratification module
!#
!#  This hook module
!#  contains all private variables pertaining to 
!#  use in the BUOY ifdefs
!#
!################################################
      module stratification_module
        implicit none

        !variables for BOUSS
        real(b8), allocatable :: ub(:,:,:,:)
	real froude,noise_ratio,noise_spconst,noise_spslope

        !variables for BUOY
	real bvfsq
        !real(b8), allocatable :: ub(:,:,:,:)

!##########################################################
      contains      
!######################################
!#   Stratification routine
!#   Formerly called buoy
!#   This is the Buoyancy Stratification 
!#   Hook routine
!#   
!#   Proper use of the Runge-kutta module
!#   Would make this routine trivial 
!#   And remove need to have icall
!#
!#####################################
      subroutine stratification(uy,uny,icall) 
      use comp
      implicit none
      include 'intvars'
      
      complex(b8) :: uy(ny,zjsz,xisz,nu)
      complex(b8) :: uny(ny,zjsz,xisz,3+nc)
      complex(b8) :: uby(ny,zjsz,xisz,3:4)
      
      integer icall
      
      go to (1,2,3), icall
      
! ----------------- icall = 1 -----------------------------------------------
     
 1	continue
     
	do 10 xp=1,xisz           
           do 11 zp=1,zjsz
              z=zjst+zp-1
              if (z.eq.nzhp) go to 11              
              do 12 y=1,ny
                 if (.not.mask(y,zp,xp)) go to 12                 
                 !these are for the time step routine -- coefficients
                 uy(y,zp,xp,3)=uy(y,zp,xp,3)+uny(y,zp,xp,4)*dt/2.
                 uy(y,zp,xp,4)=uy(y,zp,xp,4)-uny(y,zp,xp,3)*dt/2.*bvfsq
 12           continue
 11        continue
 10	continue

	return

! ----------------- icall = 2 -----------------------------------------------

 2	continue
        
	do 20 xp=1,xisz
           do 21 zp=1,zjsz
              z=zjst+zp-1
              if (z.eq.nzhp) go to 21
              do 22 y=1,ny
                 if (.not.mask(y,zp,xp)) go to 22
                 uby(y,zp,xp,3)=uy(y,zp,xp,3)
                 uby(y,zp,xp,4)=uy(y,zp,xp,4)
 22           continue
 21        continue
 20	continue
        
	return

! ----------------- icall = 3 -----------------------------------------------

 3	continue
        
	do 30 xp=1,xisz
           do 31 zp=1,zjsz
              z=zjst+zp-1
              if (z.eq.nzhp) go to 31
              do 32 y=1,ny
                 if (.not.mask(y,zp,xp)) go to 32                 
                 uy(y,zp,xp,3)=uy(y,zp,xp,3)+uby(y,zp,xp,4)*dt/2.
                 uy(y,zp,xp,4)=uy(y,zp,xp,4)-uby(y,zp,xp,3)*dt/2.*bvfsq
 32           continue
 31        continue
 30     continue
     
      return
      end
!#######################################################
      end module stratification_module
