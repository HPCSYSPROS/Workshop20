!################################################
!#      module stratification module
!#
!#  This hook module
!#  contains all private variables pertaining to 
!#  use in the BUOY ifdefs
!#
!################################################
      module stratification_module
        use comp
        implicit none
	
        !variables for BOUSS
        real(b8), allocatable :: ub(:,:,:,:)
	real froude,noise_ratio,noise_spconst,noise_spslope

#ifdef LINDBORG
        real visc_h, visc_v, difs_h, difs_v
        real F_h,ei_rate,ei_waveno,ei_scale
        complex,allocatable :: ufh(:,:,:),ufv(:,:)
#endif

!##########################################################
      contains      
!######################################
!#   Stratification routine
!#   Formerly called buoy
!#   This is the BOUSS Stratification 
!#   Hook routine
!#   
!#   Need to write out ub instead of uby
!#
!#####################################
      subroutine stratification(uy,uny,uby,icall) 
        implicit none
	include 'intvars'
 
      complex(b8) :: uy(ny,zjsz*xisz,nu)
      complex(b8) :: uny(ny,zjsz*xisz,3+nc)
      complex(b8) :: uby(ny,zjsz*xisz,3:4)
c

	integer icall
#ifdef NOTYET
	real const
	real bvf

c	if (taskid.eq.0) write (6,*) 'enter bouss, icall=',icall
c
	pi=atan(1.)*4.
	if (froude.gt.1.e6) return
	const=(2.*pi/froude)**2
c
#ifdef LINDBORG
	bvf=ei_rate**(1./3.)/F_h/ei_scale**(2./3.)
#endif
c
	go to (1,2,3), icall
c
 1	continue
c
c	if (taskid.eq.0) write (6,*) 'bouss, A; bvf=',bvf

	do 10 xp=1,xisz

	do 11 zp=1,zjsz
	z=zjst+zp-1
	if (z.eq.nzhp) go to 11

	do 12 y=1,ny
	if (.not. mask(y,zp,xp)) go to 12

#ifdef LINDBORG
 	uy(y,zp,xp,3)=uy(y,zp,xp,3)-uny(y,zp,xp,4)*dt/2.*bvf
   	uy(y,zp,xp,4)=uy(y,zp,xp,4)+uny(y,zp,xp,3)*dt/2.*bvf*beta3
#else
 	uy(y,zp,xp,3)=uy(y,zp,xp,3)-uny(y,zp,xp,4)*dt/2.*const
   	uy(y,zp,xp,4)=uy(y,zp,xp,4)+uny(y,zp,xp,3)*dt/2.*beta3
#endif
 12	continue
 11	continue
 10	continue
c
c	if (taskid.eq.0) write (6,*) ' exit bouss (icall=1)'
	return
c
 2	continue
c
	do 20 xp=1,xisz

	do 21 zp=1,zjsz
	z=zjst+zp-1
	if (z.eq.nzhp) go to 21

	do 22 y=1,ny
	if (.not. mask(y,zp,xp)) go to 22

 	uby(y,zp,xp,3)=uy(y,zp,xp,3)
	uby(y,zp,xp,4)=uy(y,zp,xp,4)
 22	continue
 21	continue
 20	continue
c	if (taskid.eq.0) write (6,*) ' exit bouss (icall=2)'
c
 3	continue

	do 30 xp=1,xisz
	do 31 zp=1,zjsz
	z=zjst+zp-1
	if (z.eq.nzhp) go to 31
	do 32 y=1,ny
	if (.not.mask(y,zp,xp)) go to 32

#ifdef LINDBORG
 	uy(y,zp,xp,3)=uy(y,zp,xp,3)-uby(y,zp,xp,4)*dt/2.*bvf
   	uy(y,zp,xp,4)=uy(y,zp,xp,4)+uby(y,zp,xp,3)*dt/2.*bvf*beta3
#else
 	uy(y,zp,xp,3)=uy(y,zp,xp,3)-uby(y,zp,xp,4)*dt/2.*const
  	uy(y,zp,xp,4)=uy(y,zp,xp,4)+uby(y,zp,xp,3)*dt/2.*beta3
#endif
 32	continue
 31	continue
 30	continue
c
#endif
c	if (taskid.eq.0) write (6,*) ' exit bouss (icall=3)'
c
      return
      end

!#######################################################
      end module stratification_module
