! Arrays to store dissipation operator coefficients for the full restricted norm case.
! $Header$

#include "cctk.h"

module dissipation_coeff

  type coefficients
    CCTK_REAL, dimension(:,:), pointer :: coeff
  end type coefficients

  type(coefficients), dimension(:), pointer :: xcoeff, ycoeff, zcoeff

  logical :: first = .true.
  logical, dimension(:), allocatable :: savedx, savedy, savedz

end module dissipation_coeff
