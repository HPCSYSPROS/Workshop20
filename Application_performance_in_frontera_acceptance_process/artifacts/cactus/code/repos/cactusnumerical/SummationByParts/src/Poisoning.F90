#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine SBP_Poisoning (ni, nj, nk, bb, offset, dvar)

  implicit none

  DECLARE_CCTK_PARAMETERS

  integer ::  ni, nj, nk
  CCTK_INT, dimension(6), intent(in) :: bb
  CCTK_INT, dimension(6), intent(in) :: offset
  CCTK_REAL, dimension(ni,nj,nk), intent(inout) :: dvar

  CCTK_INT :: il, ir, jl, jr, kl, kr
  
  il = 0; ir = 0
  jl = 0; jr = 0
  kl = 0; kr = 0
  if ( bb(1) /= 0 ) il = offset(1)
  if ( bb(2) /= 0 ) ir = offset(2)
  if ( bb(3) /= 0 ) jl = offset(3)
  if ( bb(4) /= 0 ) jr = offset(4)
  if ( bb(5) /= 0 ) kl = offset(5)
  if ( bb(6) /= 0 ) kr = offset(6)

!$omp parallel workshare
    dvar(1:il,:,:) = poison_value
    dvar(:,1:jl,:) = poison_value
    dvar(:,:,1:kl) = poison_value
    dvar(ni-ir+1:ni,:,:) = poison_value
    dvar(:,nj-jr+1:nj,:) = poison_value
    dvar(:,:,nk-kr+1:nk) = poison_value
!$omp end parallel workshare

end subroutine SBP_Poisoning
