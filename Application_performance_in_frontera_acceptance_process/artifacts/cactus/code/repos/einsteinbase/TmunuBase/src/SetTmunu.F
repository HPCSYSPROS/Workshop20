#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


      
c     Calculate the contribution to the stress energy tensor T_munu
c     which are calcualated via the CalcTmunu.inc mechanism.
c     Then make a copy of that into the T2munu variables.

      subroutine TmunuBase_SetTmunu (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
#define TMUNUBASE_SETTMUNU

#include "EinsteinBase/ADMMacros/src/macro/STRESSENERGY_declare.h"
      
      integer i, j, k
      
      do k = 1, cctk_lsh(3)
         do j = 1, cctk_lsh(2)
            do i = 1, cctk_lsh(1)
               
#include "EinsteinBase/ADMMacros/src/macro/STRESSENERGY_guts.h"

               eTtt(i,j,k) = Ttt
               
               eTtx(i,j,k) = Ttx
               eTty(i,j,k) = Tty
               eTtz(i,j,k) = Ttz
               
               eTxx(i,j,k) = Txx
               eTxy(i,j,k) = Txy
               eTxz(i,j,k) = Txz
               eTyy(i,j,k) = Tyy
               eTyz(i,j,k) = Tyz
               eTzz(i,j,k) = Tzz
               
#include "EinsteinBase/ADMMacros/src/macro/STRESSENERGY_undefine.h"
               
            end do
         end do
      end do
      
#undef TMUNUBASE_SETTMUNU
      
      eT2tt = eTtt
      
      eT2tx = eTtx
      eT2ty = eTty
      eT2tz = eTtz
      
      eT2xx = eTxx
      eT2xy = eTxy
      eT2xz = eTxz
      eT2yy = eTyy
      eT2yz = eTyz
      eT2zz = eTzz
      
      end
