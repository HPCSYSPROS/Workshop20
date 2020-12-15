#include "cctk.h"
#include "cctk_Functions.h"


  
subroutine lc_get_fortran_type_sizes(type_sizes)
  use loopcontrol_types
  implicit none
  DECLARE_CCTK_FUNCTIONS
  CCTK_POINTER :: type_sizes(3)
  
  type(lc_vec_t) :: vec(2)
  type(lc_space_t) :: space(2)
  type(lc_control_t) :: control(2)
  
  type_sizes(1) = CCTK_PointerTo(vec(2)) - CCTK_PointerTo(vec(1))
  type_sizes(2) = CCTK_PointerTo(space(2)) - CCTK_PointerTo(space(1))
  type_sizes(3) = CCTK_PointerTo(control(2)) - CCTK_PointerTo(control(1))
end subroutine lc_get_fortran_type_sizes
