#include "cctk.h"

#include "loopcontrol.h"



module loopcontrol_types
  
  implicit none
  
  ! Note: These types must correspond to the corresponding C types
  ! declared in loopcontrol.h
  
  type, bind(C) :: lc_vec_t
     CCTK_POINTER :: v(LC_DIM)
  end type lc_vec_t
  
  type, bind(C) :: lc_space_t
     type(lc_vec_t) :: min, max, step, pos
     type(lc_vec_t) :: count, idx
  end type lc_space_t
  
  type, bind(C) :: lc_control_t
     type(lc_vec_t)   :: ash
     type(lc_space_t) :: overall
     type(lc_space_t) :: coarse_thread
     CCTK_POINTER     :: coarse_thread_info_ptr
     integer          :: coarse_thread_done
     type(lc_space_t) :: coarse_loop
     type(lc_space_t) :: fine_loop
     type(lc_space_t) :: fine_thread
     CCTK_POINTER     :: selftest_array
  end type lc_control_t
  
end module loopcontrol_types
