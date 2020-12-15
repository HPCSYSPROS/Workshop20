#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine dcm_init_drift (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: n
  
  dc4_calc_error = 1.0
  dc4_calc_error_p = 1.0
  dc4_calc_error_p_p = 1.0
  
end subroutine dcm_init_drift

subroutine dcm_init_drift_recover (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: n
  
  dc4_calc_error = 1.0
  dc4_calc_error_p = 1.0
  dc4_calc_error_p_p = 1.0
  
end subroutine dcm_init_drift_recover



subroutine dcm_calculate_correction (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT,  parameter :: izero = 0
  integer,   parameter :: ik=kind(izero)
  CCTK_REAL, parameter :: one = 1.0d0
  CCTK_REAL :: pi
  
  CCTK_REAL :: size_current, size_desired
  CCTK_REAL :: pos_current(3), pos_desired(3)
  CCTK_REAL :: shape_current(3,3), shape_desired(3,3)
  CCTK_REAL :: delta_size, delta_pos(3), delta_shape(3,3)
  CCTK_REAL :: delta_size_dot, delta_pos_dot(3), delta_shape_dot(3,3)
  CCTK_REAL :: delta_size_dot2, delta_pos_dot2(3), delta_shape_dot2(3,3)
  CCTK_REAL :: delta_omega, delta_theta, delta_radius
  CCTK_REAL :: delta_omega_dot, delta_theta_dot, delta_radius_dot
  CCTK_REAL :: delta_omega_dot2, delta_theta_dot2, delta_radius_dot2
  
  integer :: i, j
  
  CCTK_REAL :: trace
  
  character(len=1000) :: msg
  
  pi = acos(-one)
  
  
  ! Rotate time levels
  dc4_calc_error_p_p = dc4_calc_error_p
  dc4_calc_error_p   = dc4_calc_error
  dc4_time_p_p = dc4_time_p
  dc4_time_p   = dc4_time
  
  dc4_delta_posx_p_p    = dc4_delta_posx_p
  dc4_delta_posx_p      = dc4_delta_posx
  dc4_delta_posy_p_p    = dc4_delta_posy_p
  dc4_delta_posy_p      = dc4_delta_posy
  dc4_delta_posz_p_p    = dc4_delta_posz_p
  dc4_delta_posz_p      = dc4_delta_posz
  
  ! Set up current time level
  dc4_calc_error = 0.0d0
  dc4_time = cctk_time
  
  pos_current(1)     = Mx
  pos_current(2)     = My
  pos_current(3)     = Mz
     
  pos_desired(1) = position_x
  pos_desired(2) = position_y
  pos_desired(3) = position_z
  if(bitant.ne.0) then
     pos_desired(3) = pos_current(3)
  endif

  ! Determine correction
  
  ! The correction has six parts: size, position, shape,
  ! rotational, inclinational, and radial
  
  delta_pos   = pos_current - pos_desired
     
  dc4_delta_posx    = delta_pos(1)
  dc4_delta_posy    = delta_pos(2)
  dc4_delta_posz    = delta_pos(3)

  delta_pos_dot(1)     = timederiv (dc4_delta_posx, &
       &                            dc4_delta_posx_p, &
       &                            dc4_delta_posx_p_p)
  delta_pos_dot(2)     = timederiv (dc4_delta_posy, &
       &                            dc4_delta_posy_p, &
       &                            dc4_delta_posy_p_p)
  delta_pos_dot(3)     = timederiv (dc4_delta_posz, &
       &                            dc4_delta_posz_p, &
       &                            dc4_delta_posz_p_p)
     
  forall (i=1:3)
     delta_pos_dot2(i)    = damped_oscillator (position_timescale, &
          &                                     position_damping_factor, &
          &                                     delta_pos(i), &
          &                                     delta_pos_dot(i))
  end forall

  dc4_delta_posx_dot2    = delta_pos_dot2(1)
  dc4_delta_posy_dot2    = delta_pos_dot2(2)
  dc4_delta_posz_dot2    = delta_pos_dot2(3)

  if(abs(delta_pos(1)).lt.min_dx) then
     dc4_delta_posx_dot2 = 0.0d0
  endif

  if(abs(delta_pos(2)).lt.min_dy) then
     dc4_delta_posy_dot2 = 0.0d0
  endif

  if(abs(delta_pos(3)).lt.min_dz) then
     dc4_delta_posz_dot2 = 0.0d0
  endif

  if(verbose_level.gt.0) then
     write(msg,"(A18,1P10E15.6)") "delta_pos: ",delta_pos(1:3)
     call CCTK_INFO(msg)
     write(msg,"(A18,1P10E15.6)") "delta_pos_dot: ",delta_pos_dot(1:3)
     call CCTK_INFO(msg)
     write(msg,"(A18,1P10E15.6)") "delta_pos_dot2: ",dc4_delta_posx_dot2, dc4_delta_posy_dot2, &
          dc4_delta_posz_dot2
     call CCTK_INFO(msg)
  endif

     
9999 continue

  
contains
  
  pure function timederiv (f0, f1, f2) result (fdot)
    CCTK_REAL, intent(in) :: f0, f1, f2
    CCTK_REAL             :: fdot
    CCTK_REAL :: dt1, dt2
    CCTK_REAL :: fdot1, fdot2
    
    dt1 = dc4_time - dc4_time_p
    dt2 = dc4_time - dc4_time_p_p
    
    if (dc4_calc_error_p/=0) then
       fdot = 0
    else if (dc4_calc_error_p_p/=0) then
       fdot = (f0 - f1) / dt1
    else
       ! f(dt1) = f(0) + dt1 f'(0) + dt1^2/2 f''(0)
       ! f(dt2) = f(0) + dt2 f'(0) + dt2^2/2 f''(0)
       ! f'(0) = [f(dt1) - f(0)] / dt1 - dt1/2 f''(0)
       ! f'(0) = [f(dt2) - f(0)] / dt2 - dt2/2 f''(0)
       fdot1 = (f0 - f1) / dt1
       fdot2 = (f0 - f2) / dt2
       fdot = (fdot1 * dt2 - fdot2 * dt1) / (dt2 - dt1)
    end if
  end function timederiv
  
  pure function damped_oscillator (timescale, damp, f, fdot) result (fdot2)
    CCTK_REAL, intent(in) :: timescale, damp, f, fdot
    CCTK_REAL             :: fdot2
    
    ! T^2 x" + 2Tdx' + x = 0
    !    x(t) = C exp(-t/T)
    ! x" = -1/T^2 (2Tdx' + x)
    fdot2 = - 1/timescale**2 * (2 * timescale * damp * fdot + f)
  end function damped_oscillator
  
end subroutine dcm_calculate_correction



! subroutine dcm_correct_drift (CCTK_ARGUMENTS)
!   implicit none
!   DECLARE_CCTK_ARGUMENTS
!   DECLARE_CCTK_FUNCTIONS
!   DECLARE_CCTK_PARAMETERS
!   
!   CCTK_INT,  parameter :: izero = 0
!   integer,   parameter :: ik=kind(izero)
!   CCTK_REAL, parameter :: local_eps = 1.0e-12
!   
!   CCTK_REAL :: x0(3), dx(3), dt
!   
!   CCTK_REAL :: xpos(3), radius
!   CCTK_REAL :: lpos(3), lradius, lnormal(3)
!   
!   CCTK_REAL :: delta_shape_dot2(3,3), beta_dot2(3)
!   
!   CCTK_REAL :: pos_desired(3), radius_desired
! 
!   CCTK_REAL :: rho0, A2
!   
!   integer :: i, j, k
!   integer :: a, b, c
!   
!   character(len=1000) :: msg
!   
!   integer, parameter :: radial_correction_linear  = 1
!   integer, parameter :: radial_correction_annular = 2
!   integer :: radial_correction
!   
!   
!   
!   x0(:) = CCTK_ORIGIN_SPACE(:)
!   dx(:) = CCTK_DELTA_SPACE(:)
!   dt    = CCTK_DELTA_TIME
!   
!   if (cctk_time < first_driftcorrect_time) goto 9999
!      
!   ! Check for valid horizon data
!   if (dc4_calc_error/=0) then
!      call CCTK_WARN (0, "No valid drift correction data found")
!      goto 9999
!   end if
! 
! 
!   do k = 1, cctk_lsh(3)
!      do j = 1, cctk_lsh(2)
!         do i = 1, cctk_lsh(1)
!               
!            xpos(1) = x0(1) + dx(1) * (cctk_lbnd(1) + i - 1)
!            xpos(2) = x0(2) + dx(2) * (cctk_lbnd(2) + j - 1)
!            xpos(3) = x0(3) + dx(3) * (cctk_lbnd(3) + k - 1)
!            radius = sqrt(sum(xpos(1:2)**2))
!               
!            lpos = xpos - (/ position_x, position_y, position_z /)
!            lradius = sqrt(sum(lpos**2))
!            lnormal = lpos / (lradius + local_eps)
! 
!            ! Position correction
!            if (do_position_correction/=0) then
!               dtbetax(i,j,k) = dtbetax(i,j,k) &
!                    - dc4_delta_posx_dot2* &
!                    exp(-lradius*lradius*position_correction_falloff);
!               dtbetay(i,j,k) = dtbetay(i,j,k) &
!                    - dc4_delta_posy_dot2* &
!                    exp(-lradius*lradius*position_correction_falloff);
!               dtbetaz(i,j,k) = dtbetaz(i,j,k) &
!                    - dc4_delta_posz_dot2* &
!                    exp(-lradius*lradius*position_correction_falloff);
!            end if
!               
!         end do
!      end do
!   end do
!      
! 9999 continue
!   
! end subroutine dcm_correct_drift
