#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)


subroutine ZelmaniAnalysis_UpdateEntropy(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  integer :: keytemp
  integer :: keyerr(1)
  integer :: anyerr
  integer :: i,j,k
  real*8  :: dummy(1),rf_precision
  integer :: n,eoskey
  eoskey = 4
  rf_precision = 1.0d-12

  n=1
  keytemp = 1
  keyerr = 0
  anyerr = 0

  ! this is a bit slow, but will be called only when
  ! output is requested, so should be fine
  !$OMP PARALLEL DO PRIVATE(keyerr,anyerr,dummy)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
          
           call EOS_Omni_short(eoskey,keytemp,rf_precision,n,&
                rho(i,j,k),dummy,temperature(i,j,k),y_e(i,j,k),&
                dummy,entropy(i,j,k),dummy,dummy,dummy,dummy,&
                dummy,keyerr,anyerr)

           if(anyerr.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(0,"EOS error in ZelmaniM1 Prim2Con EOS call. Improve error message. Aborting")
              !$OMP END CRITICAL
           endif
           
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO  



end subroutine ZelmaniAnalysis_UpdateEntropy

