#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine deriv2_mixed ( cctkGH, dir1, dir2, var, ni, nj, nk, &
                                              dvar2, table_handle )

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_POINTER, intent(IN) :: cctkGH
  CCTK_INT, intent(IN) :: dir1, dir2
  CCTK_REAL, dimension(ni,nj,nk), intent(IN) :: var
  CCTK_INT, intent(IN) :: ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(OUT) :: dvar2
  CCTK_INT, intent(IN) :: table_handle
  CCTK_REAL, dimension(:,:,:), allocatable :: tmp

  allocate ( tmp(ni,nj,nk) )

  call Diff_gv ( cctkGH, dir1, var, tmp, table_handle)
  call Diff_gv ( cctkGH, dir2, tmp, dvar2, table_handle )

  deallocate ( tmp )

end subroutine deriv2_mixed
