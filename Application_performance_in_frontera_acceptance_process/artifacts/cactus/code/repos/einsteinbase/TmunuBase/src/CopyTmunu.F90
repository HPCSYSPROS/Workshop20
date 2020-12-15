#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



! Calculate the T2munu copy.
! On input, T2munu contains the contribution to Tmunu that is calculated
! via the CalcTmunu.inc mechanism, and Tmunu contains the complete stress
! energy tensor.
! On output, T2munu should contain the contribution to Tmunu that is
! calculated via the AddToTmunu mechanism.
! We just need to calculate the difference.

subroutine TmunuBase_CopyTmunu (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  if (support_old_CalcTmunu_mechanism == 0) then
     call CCTK_WARN (0, "internal error")
  end if
  
  if (stress_energy_state == 0) then
     call CCTK_WARN (1, "The stress energy tensor does not have storage")
     return
  end if
  
  stress_energy_2_state = stress_energy_state
  
  eT2tt = eTtt - eT2tt
  
  eT2tx = eTtx - eT2tx
  eT2ty = eTty - eT2ty
  eT2tz = eTtz - eT2tz
  
  eT2xx = eTxx - eT2xx
  eT2xy = eTxy - eT2xy
  eT2xz = eTxz - eT2xz
  eT2yy = eTyy - eT2yy
  eT2yz = eTyz - eT2yz
  eT2zz = eTzz - eT2zz
  
end subroutine TmunuBase_CopyTmunu
