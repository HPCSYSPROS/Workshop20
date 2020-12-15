#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine PNSMapper_Map(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer i,j,k,l,nx,ny,nz
  integer istat
  integer lower_index, upper_index
  real*8 radius
  real*8 tempvel, cosf, sinf, sint, cost, rxy
  real*8 tenthalpy, rdetg,w,vlowx,vlowy,vlowz
! begin EOS Omni vars
  integer :: n,keytemp,anyerr,keyerr(1)
  integer :: eoskey
  real*8  :: eos_rf_prec
  real*8  :: xdummy(1)
! end EOS Omni vars



  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  call CCTK_INFO("Mapping 1D profile to 3D")

  gxx(:,:,:) = 0.0d0
  gxy(:,:,:) = 0.0d0
  gxz(:,:,:) = 0.0d0
  gyy(:,:,:) = 0.0d0
  gyz(:,:,:) = 0.0d0
  gzz(:,:,:) = 0.0d0

  kxx(:,:,:) = 0.0d0
  kxy(:,:,:) = 0.0d0
  kxz(:,:,:) = 0.0d0
  kyy(:,:,:) = 0.0d0
  kyz(:,:,:) = 0.0d0
  kzz(:,:,:) = 0.0d0

  betax(:,:,:) = 0.0d0
  betay(:,:,:) = 0.0d0
  betaz(:,:,:) = 0.0d0

  istat = 0
  call CCTK_QueryGroupStorage(istat,cctkGH,"hydrobase::Y_e")
  if(istat.eq.0) then
     call CCTK_WARN(0,"Must enable storage for hydrobase::Y_e")
  endif

  istat = 0
  call CCTK_QueryGroupStorage(istat,cctkGH,"hydrobase::temperature")
  if(istat.eq.0) then
     call CCTK_WARN(0,"Must enable storage for hydrobase::temperature")
  endif

  ! mapping loop
!$OMP PARALLEL DO PRIVATE(i,j,k,tempvel,radius,rxy,cosf,sinf,sint,lower_index,upper_index)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           radius = sqrt(x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)

           ! find lower index in our 1D radius array:
           call PNSMapper_F_find_index(num_1D,rad1d,radius,upper_index,lower_index)

           ! fluid variables
           if(radius .lt. rad1d(Num_1D)) then
              call PNSMapper_F_Map(rho(i,j,k),radius,rho1d,rad1d,Num_1D,lower_index)
              call PNSMapper_F_Map(temperature(i,j,k),radius,temp1d,rad1d,Num_1D,lower_index)
              call PNSMapper_F_Map(y_e(i,j,k),radius,ye1d,rad1d,Num_1D,lower_index)
              call PNSMapper_F_Map(tempvel,radius,vel1d,rad1d,Num_1D,lower_index)
           else
              rho(i,j,k) = 1.0d-20
              temperature(i,j,k) = GRHydro_hot_atmo_temp
              y_e(i,j,k) = GRHydro_hot_atmo_Y_e
              tempvel = 0.0d0
           end if

           if(rho(i,j,k).le.rho_abs_min) then
              rho(i,j,k) = rho_abs_min
              temperature(i,j,k) = GRHydro_hot_atmo_temp
              y_e(i,j,k) = GRHydro_hot_atmo_Y_e
              tempvel = 0.0d0
           endif

           if(initial_vel_zero.ne.0) then
              tempvel = 0.0d0
           endif

           ! floor Y_e for EOS
           y_e(i,j,k) = max(Y_e_min,y_e(i,j,k))

           ! set up 3D velocity
           rxy = sqrt( x(i,j,k)**2 + y(i,j,k)**2 )

           if(radius .le. 1.0d-5) then
              cosf = 0.0d0
              sinf = 0.0d0
              sint = 0.0d0
              cost = 0.0d0
           else
              if (rxy .le. 1.d-5) then
                 cosf = 0.0d0
                 sinf = 0.0d0
              else
                 cosf = x(i,j,k) / rxy
                 sinf = y(i,j,k) / rxy
              endif

              sint = rxy / radius
              
              if(abs(z(i,j,k)).lt.1.0d-5) then
                 cost = 0.0d0
              else
                 cost = z(i,j,k) / radius
              endif
           
           endif

           ! special treatment: set vel to zero in the center                                            
           if (radius .gt. rad1d(1)) then
              vel(i,j,k,1) = tempvel * cosf * sint
              vel(i,j,k,2) = tempvel * sinf * sint
              vel(i,j,k,3) = tempvel * cost
           else
              !tempvel = vel1D(1)/rad1d(1) * radius
              vel(i,j,k,1) = tempvel * cosf * sint
              vel(i,j,k,2) = tempvel * sinf * sint
              vel(i,j,k,3) = tempvel * cost
           endif

           ! metric vars
           call PNSMapper_F_Map(alp(i,j,k),radius,alp1d,rad1d,Num_1D,lower_index)
           call PNSMapper_F_Map(gxx(i,j,k),radius,psi1d,rad1d,Num_1D,lower_index)

           gxx(i,j,k) = gxx(i,j,k)**4
           gyy(i,j,k) = gxx(i,j,k)
           gzz(i,j,k) = gxx(i,j,k)

        enddo
     enddo
  enddo
!$OMP END PARALLEL DO 

  call CCTK_INFO("Done mapping PNS to 3D grid!!!!")

  ! EOS and conserved vars loop
  n = 1
  keytemp = 1
  keyerr = 0
  anyerr = 0
  eoskey = 4
  eos_rf_prec = 1.0d-10
  !$OMP PARALLEL DO PRIVATE(i,j,k,keyerr,anyerr,xdummy,rdetg,w,vlowx,vlowy,vlowz,tenthalpy)
  do k=1,nz
     do j=1,ny
        do i=1,nx
     
           ! get pressure and eps for each interpolated point
           call EOS_Omni_short(eoskey,keytemp,eos_rf_prec,&
                n,rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                xdummy,entropy(i,j,k),xdummy,xdummy,&
                xdummy,xdummy,xdummy,keyerr,anyerr)           
           if(anyerr.ne.0) then
              call CCTK_WARN(0,"EOS problem in PNSMapper!!!")
           endif

           rdetg = sqrt(gxx(i,j,k)**3)
           w = 1.0d0 / &
                sqrt(1.0d0 - (gxx(i,j,k)*vel(i,j,k,1)*vel(i,j,k,1) &
                + gyy(i,j,k)*vel(i,j,k,2)*vel(i,j,k,2) &
                + gzz(i,j,k)*vel(i,j,k,3)*vel(i,j,k,3) ) )
           
           vlowx = gxx(i,j,k)*vel(i,j,k,1) &
                + gxy(i,j,k)*vel(i,j,k,2)  &
                + gxz(i,j,k)*vel(i,j,k,3)
           vlowy = gxy(i,j,k)*vel(i,j,k,1) &
                + gyy(i,j,k)*vel(i,j,k,2)  &
                + gyz(i,j,k)*vel(i,j,k,3)
           vlowz = gxz(i,j,k)*vel(i,j,k,1) &
                + gyz(i,j,k)*vel(i,j,k,2)  &
                + gzz(i,j,k)*vel(i,j,k,3)

           dens(i,j,k) = rdetg*w*rho(i,j,k)

           tenthalpy = 1.0d0 + eps(i,j,k) + press(i,j,k) / rho(i,j,k)

           tau(i,j,k) = rdetg*( (rho(i,j,k)*(1.0d0+eps(i,j,k))+press(i,j,k))*w*w - press(i,j,k)) &
                - dens(i,j,k)

           w_lorentz(i,j,k) = w

           scon(i,j,k,1) = rdetg*rho(i,j,k)*tenthalpy*(w**2) &
                *vlowx
           scon(i,j,k,2) = rdetg*rho(i,j,k)*tenthalpy*(w**2) &
                *vlowy
           scon(i,j,k,3) = rdetg*rho(i,j,k)*tenthalpy*(w**2) &
                *vlowz

           Y_e_con(i,j,k) = dens(i,j,k)*Y_e(i,j,k)


        enddo
     enddo
  enddo
!$OMP END PARALLEL DO


  call CCTK_INFO("Done initial prim2con!!!!")
  
!  i = 21
!  j = 17
!  k = 8
!  write(6,*) y_e(i,j,k),dens(i,j,k),sqrt(gxx(i,j,k)**3),y_e_con(i,j,k)


end subroutine PNSMapper_Map

subroutine PNSMapper_F_Map(point_value,point_radius,parray,pradius,zones,lower_index)

  implicit none

  CCTK_REAL point_value, point_radius
  CCTK_REAL pradius(*), parray(*)
  integer zones
  integer upper_index, lower_index

  upper_index = lower_index+1

  if (point_radius .ge. pradius(1) .and. point_radius .lt. pradius(zones) ) then

!     call PNSMapper_F_find_index(zones,pradius,point_radius,upper_index,lower_index)

     call PNSMapper_F_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index), point_radius, point_value )

  else if (point_radius .lt. pradius(1)) then

     call PNSMapper_F_linterp(pradius(1),pradius(2), &
          parray(1),parray(2),point_radius,point_value)

  else if (point_radius .gt. pradius(zones)) then
! linear extrapolation                                                                                   
     call PNSMapper_F_linterp(pradius(zones-1),pradius(zones), &
          parray(zones-1),parray(zones),point_radius,point_value)
  endif

end subroutine PNSMapper_F_Map


subroutine PNSMapper_F_linterp(x1,x2,y1,y2,x,y)

  implicit none

  CCTK_REAL slope,x1,x2,y1,y2,x,y

  if (x2.lt.x1) then
     call CCTK_WARN(0,"Error in linterp!")
  endif

  slope = (y2 - y1) / (x2 - x1)


  y = slope*(x-x1) + y1

end subroutine PNSMapper_F_linterp


subroutine PNSMapper_F_find_index(zones,array,goal,upper_index,lower_index)

  implicit none

  integer zones,i
  CCTK_REAL array(*)
  CCTK_REAL goal
  CCTK_INT middle_index,upper_index,lower_index

  lower_index = 1
  upper_index = zones
  
  do i = 2,zones
    if (goal<array(i)) then
      lower_index = i-1
      upper_index = i
      exit
    endif 
  enddo  

  return 

  do while ( (upper_index - lower_index) .gt. 1 )
     middle_index = (lower_index + upper_index) / 2
     if ( (goal .ge. array(lower_index)) &
          .and. (goal .le. array(middle_index)) ) then
        upper_index = middle_index
     else
        if ( (goal .ge. array(middle_index)) &
          .and. (goal .le. array(upper_index)) ) then
           lower_index = middle_index
        endif
     endif
  enddo

  if (upper_index==lower_index) then
    if (goal>=array(lower_index)) then
      upper_index = upper_index+1
    else
      lower_index = lower_index-1
    endif    
  endif 

end subroutine PNSMapper_F_find_index

