 /*@@
   @file      GRHydro_Particle.F90
   @date      Wed Mar 31 11:10:52 2004
   @author    Ian Hawke
   @desc 
   Track coordinates of particles.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydroParticleRHS
   @date       Wed Mar 31 11:11:26 2004
   @author     Ian Hawke
   @desc 
   Track particles.
   Most of this code is taken from Peter Dieners EHFinder - generator tracking.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydroParticleRHS(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i
  CCTK_INT :: interp_handle, table_handle, status, coord_system_handle

  character(len=200) :: particle_interp
  character(len=128) :: warn_message
  CCTK_INT :: particle_interp_len
  character(len=7) :: particle_order

  CCTK_INT, dimension(1) :: lsh
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(7) :: out_arrays
  CCTK_INT, dimension(7) :: in_arrays
  CCTK_INT, dimension(7), parameter :: op_indices = (/ 0, 1, 2, 3, 4, &
                                                        5, 6 /), &
                                        op_codes = (/ 0, 0, 0, 0, 0, &
                                                      0, 0 /)
  CCTK_INT, dimension(7) :: out_types

  out_types = CCTK_VARIABLE_REAL

! Convert the particle_interpolator string parameter to a Fortran string.
  call CCTK_FortranString ( particle_interp_len, particle_interpolator, &
       particle_interp )

! Get the corresponding interpolator handle.
  call CCTK_InterpHandle ( interp_handle, particle_interp )

  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation.'
    warn_message = trim(warn_message)//' Forgot to activate an implementation'
    warn_message = trim(warn_message)//' providing interpolation operators?'
    call CCTK_ERROR ( trim(warn_message) )
    STOP
  end if

! Convert the interpolation order parameter to a Fortran string to be placed
! in the interpolator table. Note that the order is assumed to contain only
! 1 digit.
  write(particle_order,'(a6,i1)') 'order=', particle_interpolation_order

! Create the table directly from the string.
  call Util_TableCreateFromString ( table_handle, particle_order )
  if ( table_handle .lt. 0 ) then
    call CCTK_ERROR ( 'Cannot create parameter table for interpolator' )
    STOP
  end if

! Get the 3D coordinate system handle.
  call CCTK_CoordSystemHandle ( coord_system_handle, 'cart3d' )
  if ( coord_system_handle .lt. 0) then
    warn_message = 'Cannot get handle for cart3d coordinate system.'
    warn_message = trim(warn_message)//' Forgot to activate an implementation'
    warn_message = trim(warn_message)//' providing coordinates?'
    call CCTK_ERROR ( trim(warn_message) )
    STOP
  endif

! Find out how many interpolation points are located on this processor.
  call CCTK_GrouplshGN ( status, cctkGH, 1, lsh, 'GRHydro::particles' )
  if ( status .lt. 0 ) then
    call CCTK_ERROR( 'cannot get local size for particles')
    STOP
  end if

! Set the pointers to the output arrays.
  out_arrays(1) = CCTK_PointerTo(particle_vx)
  out_arrays(2) = CCTK_PointerTo(particle_vy)
  out_arrays(3) = CCTK_PointerTo(particle_vz)
  out_arrays(4) = CCTK_PointerTo(particle_alp)
  out_arrays(5) = CCTK_PointerTo(particle_betax)
  out_arrays(6) = CCTK_PointerTo(particle_betay)
  out_arrays(7) = CCTK_PointerTo(particle_betaz)

!   Set the pointers to the points to be interpolated to.
  interp_coords(1) = CCTK_PointerTo(particle_x)
  interp_coords(2) = CCTK_PointerTo(particle_y)
  interp_coords(3) = CCTK_PointerTo(particle_z)

!   Set the indices to the input grid functions.
  call CCTK_VarIndex ( in_arrays(1), 'GRHydro::velx' )
  call CCTK_VarIndex ( in_arrays(2), 'GRHydro::vely' )
  call CCTK_VarIndex ( in_arrays(3), 'GRHydro::velz' )
  call CCTK_VarIndex ( in_arrays(4), 'admbase::alp' )
  call CCTK_VarIndex ( in_arrays(5), 'admbase::betax' )
  call CCTK_VarIndex ( in_arrays(6), 'admbase::betay' )
  call CCTK_VarIndex ( in_arrays(7), 'admbase::betaz' )

!     Set the operand indices table entry
  call Util_TableSetIntArray ( status, table_handle, 7, &
       op_indices(1:7), 'operand_indices' )
  if ( status .lt. 0 ) then
    warn_message = 'Cannot set operand indices array in parameter table'
    call CCTK_ERROR ( trim(warn_message) )
    STOP
  endif

!     Set the corresponding table entry for the operation codes.
  call Util_TableSetIntArray ( status, table_handle, 7, &
       op_codes(1:7), 'operation_codes' )
  if ( status .lt. 0 ) then
    warn_message = 'Cannot set operation codes array in parameter table'
    call CCTK_ERROR ( trim(warn_message) )
  endif

!     Call the interpolator.
  call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
       table_handle, coord_system_handle, &
       lsh(1), CCTK_VARIABLE_REAL, &
       interp_coords, 7, in_arrays(1:7), &
       7, out_types(1:7), out_arrays(1:7) )

  if ( status .lt. 0 ) then
    call CCTK_INFO ( 'Interpolation failed.' )
  end if

!     For each point on this processor calculate the right hand side of the
!     characteristic evolution equation.
  do i = 1, lsh(1)

    particle_x_rhs(i) = particle_alp(i) * particle_vx(i) - &
         particle_betax(i)
    particle_y_rhs(i) = particle_alp(i) * particle_vy(i) - &
         particle_betay(i)
    particle_z_rhs(i) = particle_alp(i) * particle_vz(i) - &
         particle_betaz(i)

  end do
  
  return

end subroutine GRHydroParticleRHS

/*@@
   @routine    GRHydroParticleInitial
   @date       Wed Mar 31 11:53:25 2004
   @author     Ian Hawke
   @desc 
   Sets up the initial coordinates of the particles. Can be overwritten
   by something scheduled in HydroBase_Initial.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydroParticleInitial(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, dimension(1) :: lsh, lbnd
  CCTK_INT :: i, status
  CCTK_REAL :: radius, angle, xmin, xmax, ymin, ymax, zmin, zmax, twopi

  twopi = 4.d0 * atan(1.d0)

  if (number_of_particles .gt. 0) then

! Find out how many interpolation points are located on this processor.
    call CCTK_GrouplshGN ( status, cctkGH, 1, lsh, 'GRHydro::particles' )
    if ( status .lt. 0 ) then
      call CCTK_ERROR ( 'cannot get local size for particles' )
      STOP
    end if
    call CCTK_GrouplbndGN ( status, cctkGH, 1, lbnd, 'GRHydro::particles' )
    if ( status .lt. 0 ) then
      call CCTK_ERROR ( 'cannot get lower bounds for particles' )
      STOP
    end if

    call CCTK_CoordRange(status,cctkGH,xmin,xmax,-1,"x","cart3d")
    call CCTK_CoordRange(status,cctkGH,ymin,ymax,-1,"y","cart3d")
    call CCTK_CoordRange(status,cctkGH,zmin,zmax,-1,"z","cart3d")


    radius = 0.1d0 * min(xmax, ymax, zmax)

    do i = 1, lsh(1)

      angle = dble(lbnd(1)+i-1) / dble(number_of_particles) * twopi 

      particle_x(i) = radius * sin(angle)
      particle_y(i) = radius * cos(angle)
      particle_z(i) = 0.d0

      particle_x_p(i) = particle_x(i)
      particle_y_p(i) = particle_y(i)
      particle_z_p(i) = particle_z(i)

      particle_x_p_p(i) = particle_x(i)
      particle_y_p_p(i) = particle_y(i)
      particle_z_p_p(i) = particle_z(i)

      particle_x_rhs(i) = 0.d0
      particle_y_rhs(i) = 0.d0
      particle_z_rhs(i) = 0.d0
      
      particle_vx(i) = 0.d0
      particle_vy(i) = 0.d0
      particle_vz(i) = 0.d0
      particle_alp(i) = 0.d0
      particle_betax(i) = 0.d0
      particle_betay(i) = 0.d0
      particle_betaz(i) = 0.d0

    end do
    
  end if
  
end subroutine GRHydroParticleInitial
