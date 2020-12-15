 /*@@
   @file      GRHydro_TVDReconstruct.F90
   @date      Sat Jan 26 02:11:44 2002
   @author    Luca Baiotti
   @desc 
   The TVD reconstruction routine.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

 /*@@
   @routine    tvdreconstruct
   @date       Sat Jan 26 02:12:12 2002
   @author     Luca Baiotti
   @desc 
   Performs slope limited TVD reconstruction on the given input GF
   @enddesc 
   @calls     
   @calledby   
   @history 
   Follows (in philosophy) old code by Ian Hawke
   @endhistory 

@@*/

subroutine tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
     orig, bextp, bextm, trivial_rp, hydro_excision_mask)
  
  USE GRHydro_Scalars
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, xoffset, yoffset, zoffset, nx, ny, nz
  CCTK_REAL, dimension(nx, ny, nz) :: orig, bextp, bextm
  CCTK_REAL :: dupw, dloc, delta, ratio, hdelta
  logical, dimension(nx,ny,nz) :: trivial_rp
  CCTK_INT, dimension(nx,ny,nz) :: hydro_excision_mask


  bextp = 0.d0
  bextm = 0.d0

!!$ Initially all Riemann problems are NON-trivial

  trivial_rp = .false.
    ! constraint transport needs to be able to average fluxes in the directions
    ! other that flux_direction, which in turn need the primitives on interfaces
    !$OMP PARALLEL DO PRIVATE(i,j,k,dupw,dloc,delta,ratio,hdelta)
    do k = GRHydro_stencil, nz - GRHydro_stencil + 1 + transport_constraints*(1-zoffset)
      do j = GRHydro_stencil, ny - GRHydro_stencil + 1 + transport_constraints*(1-yoffset)
        do i = GRHydro_stencil, nx - GRHydro_stencil + 1 + transport_constraints*(1-xoffset)
          if (GRHydro_enable_internal_excision /= 0 .and. &
              (hydro_excision_mask(i,j,k) .ne. 0)) then
            trivial_rp(i-xoffset, j-yoffset, k-zoffset) = .true.
            trivial_rp(i, j, k) = .true.
            bextm(i, j, k) = orig(i, j, k)
            bextp(i, j, k) = orig(i, j, k)
            if (GRHydro_enable_internal_excision /= 0 .and. &
                (hydro_excision_mask(i+xoffset,j+yoffset,k+zoffset) .eq. 0)) then
              bextm(i, j, k) = orig(i+xoffset, j+yoffset, k+zoffset)
              bextp(i, j, k) = orig(i+xoffset, j+yoffset, k+zoffset)
            end if
          else if (GRHydro_enable_internal_excision /= 0 .and. &
                   ((hydro_excision_mask(i-xoffset,j-yoffset,k-zoffset) .ne. 0) .or. &
                    (hydro_excision_mask(i+xoffset,j+yoffset,k+zoffset) .ne. 0))) then
            bextm(i, j, k) = orig(i, j, k)
            bextp(i, j, k) = orig(i, j, k)
          else
            dupw = orig(i, j, k) - orig(i-xoffset, j-yoffset, k-zoffset)
            dloc = orig(i+xoffset, j+yoffset, k+zoffset) - orig(i, j, k)

            if (MINMOD) then
               delta = minmod_func(dupw,dloc)
            else if (MC2) then
!!$            This is the van Leer MC slope limiter
               if (dupw*dloc < 0.d0) then
                  delta=0.d0
               else 
                  delta=sign(min(2.d0*abs(dupw),2.d0*abs(dloc),&
                       0.5d0*(abs(dupw)+abs(dloc))),dupw+dloc)
               end if
            else if (SUPERBEE) then
               if (dupw*dloc < 0.d0) then
                  delta=0.d0
               else 
                  delta=sign(max(min(2.d0*abs(dupw),abs(dloc)),&
                       min(2.d0*abs(dloc),abs(dupw))),dupw+dloc)
               end if
            else
               call CCTK_ERROR("Type of limiter not recognized")
               STOP
               ! NOTREACHED
               delta = 0d0 
            end if
            hdelta = 0.5d0 * delta 
            bextm(i, j, k) = orig(i, j, k) - hdelta
            bextp(i, j, k) = orig(i, j, k) + hdelta
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

contains
  function minmod_func(a_in,b_in) result(minmod_result)
    implicit none
    CCTK_REAL,intent(IN)::a_in,b_in
    CCTK_REAL::minmod_result
    
    minmod_result=0.5D0*(sign(1.0d0,a_in)+sign(1.0d0,b_in))*min(abs(a_in),abs(b_in))
  end function minmod_func

end subroutine tvdreconstruct

