#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine CheckFortranParameters (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  if (one /= 1) then
     call CCTK_PARAMWARN ("Fortran parameters do not work -- check the C/Fortran language interface")
  end if
end subroutine CheckFortranParameters
