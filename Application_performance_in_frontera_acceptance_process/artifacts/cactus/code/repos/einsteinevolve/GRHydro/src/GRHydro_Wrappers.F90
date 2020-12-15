 /*@@
   @file      GRHydro_Wrappers.F90
   @date      Sat Jun 29 18:57:07 PDT 2013
   @author    Roland Haas
   @desc
      Wrapper routines to provide a consistent interface to the aliased
      functions even if the interface to the underlying implementation changes.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Con2PrimPolyWrapper(handle, dens, sx, sy, sz, tau, rho, &
     velx, vely, velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, detg, x, y, z, r, GRHydro_rho_min, GRHydro_reflevel, GRHydro_C2P_failed)

  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, INTENT(INOUT) :: dens, sx, sy, sz
  CCTK_REAL, INTENT(OUT) :: tau
  CCTK_REAL, INTENT(INOUT) :: rho
  CCTK_REAL, INTENT(INOUT) :: velx, vely, velz, epsilon, press, w_lorentz
  CCTK_REAL, INTENT(IN) :: uxx, uxy, uxz, uyy, uyz, uzz, detg
  CCTK_REAL, INTENT(IN) :: x, y, z, r, GRHydro_rho_min

  CCTK_INT, INTENT(IN) :: handle, GRHydro_reflevel
  CCTK_REAL, INTENT(OUT) :: GRHydro_C2P_failed

  CCTK_REAL :: sdetg

  sdetg = sqrt(detg)

  call Con2Prim_ptPolytype(handle, dens, sx, sy, sz, tau, rho, &
     velx, vely, velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, sdetg, x, y, z, r, GRHydro_rho_min, GRHydro_reflevel, &
     GRHydro_C2P_failed)

end subroutine

subroutine Con2PrimGenWrapper(handle, dens, sx, sy, sz, tau, rho, velx, vely, &
                              velz, epsilon, pressure, w_lorentz, uxx, uxy, &
                              uxz, uyy, uyz, uzz, det, x, y, z, r, epsnegative, &
                              GRHydro_rho_min, pmin, epsmin, GRHydro_reflevel, &
                              retval)

  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle
  CCTK_REAL, INTENT(INOUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(INOUT) :: rho, velx, vely, velz, epsilon, pressure
  CCTK_REAL, INTENT(INOUT) :: w_lorentz
  CCTK_REAL, INTENT(IN) :: uxx, uxy, uxz, uyy, uyz, uzz, det
  CCTK_REAL, INTENT(IN) :: x, y, z, r
  CCTK_INT, INTENT(OUT) :: epsnegative
  CCTK_REAL, INTENT(IN) :: GRHydro_rho_min, pmin, epsmin
  CCTK_INT, INTENT(IN) :: GRHydro_reflevel
  CCTK_REAL, INTENT(OUT) :: retval

  CCTK_REAL :: sdetg
  CCTK_INT, parameter :: cctk_iteration = -1, ii = -1, jj = -1, kk = -1
  logical :: log_epsnegative

  sdetg = sqrt(det)

  call Con2Prim_pt(cctK_iteration, ii, jj, kk, &
                   handle, dens, sx, sy, sz, tau, rho, velx, vely, &
                   velz, epsilon, pressure, w_lorentz, uxx, uxy, &
                   uxz, uyy, uyz, uzz, sdetg, x, y, z, r, log_epsnegative, &
                   GRHydro_rho_min, pmin, epsmin, GRHydro_reflevel, &
                   retval)
  if(log_epsnegative) then
    epsnegative = 1
  else
    epsnegative = 0
  end if

end subroutine

subroutine Con2PrimGenMWrapper(handle, keytemp, prec, gamma_eos, dens, sx, sy, sz, &
                               tau, Bconsx, Bconsy, Bconsz, y_e, temp, rho, velx, &
                               vely, velz, epsilon, pressure, Bvecx, Bvecy, Bvecz, &
                               Bvecsq, w_lorentz, gxx, gxy, gxz, gyy, gyz, gzz, uxx, &
                               uxy, uxz, uyy, uyz, uzz, det, epsnegative, retval)

  use Con2PrimM_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle, keytemp
  CCTK_REAL, INTENT(IN) :: prec, gamma_eos
  CCTK_REAL, INTENT(INOUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(IN) :: Bconsx, Bconsy, Bconsz
  CCTK_REAL, INTENT(INOUT) :: y_e,  temp
  CCTK_REAL, INTENT(INOUT) :: rho, velx, vely, velz, epsilon, pressure
  CCTK_REAL, INTENT(OUT) :: Bvecx, Bvecy, Bvecz, Bvecsq
  CCTK_REAL, INTENT(INOUT) :: w_lorentz
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: uxx, uxy, uxz, uyy, uyz, uzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_INT, INTENT(OUT) :: epsnegative
  CCTK_REAL, INTENT(OUT) :: retval

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  if(handle.eq.1 .or. handle.eq.2) then
    call GRHydro_Con2PrimM_ptold(handle, gamma_eos, dens, sx, sy, sz, &
                                 tau, Bconsx, Bconsy, Bconsz, rho, velx, &
                                 vely, velz, epsilon, pressure, Bvecx, Bvecy, Bvecz, &
                                 Bvecsq, w_lorentz, gxx, gxy, gxz, gyy, gyz, gzz, uxx, &
                                 uxy, uxz, uyy, uyz, uzz, sdetg, epsnegative, retval)
  else
    call GRHydro_Con2PrimM_pt2(handle, keytemp, prec, prec, gamma_eos, dens, sx, sy, sz, &
                               tau, Bconsx, Bconsy, Bconsz, y_e, temp, rho, velx, &
                               vely, velz, epsilon, pressure, Bvecx, Bvecy, Bvecz, &
                               Bvecsq, w_lorentz, gxx, gxy, gxz, gyy, gyz, gzz, uxx, &
                               uxy, uxz, uyy, uyz, uzz, sdetg, epsnegative, retval)
  end if

end subroutine

subroutine Con2PrimPolyMWrapper(handle, gamma_eos, dens, sx, sy, sz, sc,&
                                Bconsx, Bconsy, Bconsz, rho, velx, vely, velz,&
                                epsilon, pressure, Bvecx, Bvecy, Bvecz, Bvecsq,&
                                w_lorentz, gxx, gxy, gxz, gyy, gyz, gzz, uxx,&
                                uxy, uxz, uyy, uyz, uzz, det, epsnegative,&
                                retval)

  use Con2PrimM_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle
  CCTK_REAL, INTENT(IN) :: gamma_eos
  CCTK_REAL, INTENT(INOUT) :: dens, sx, sy, sz, sc
  CCTK_REAL, INTENT(INOUT) :: Bconsx, Bconsy, Bconsz
  CCTK_REAL, INTENT(INOUT) :: rho, velx, vely, velz, epsilon, pressure
  CCTK_REAL, INTENT(OUT) :: Bvecx, Bvecy, Bvecz
  CCTK_REAL, INTENT(OUT) :: Bvecsq
  CCTK_REAL, INTENT(INOUT) :: w_lorentz
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: uxx, uxy, uxz, uyy, uyz, uzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_INT, INTENT(OUT) :: epsnegative
  CCTK_REAL, INTENT(OUT) :: retval

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call GRHydro_Con2PrimM_Polytype_pt(handle, gamma_eos, dens, sx, sy, sz, sc, &
                                     Bconsx, Bconsy, Bconsz, rho, velx, vely, velz, &
                                     epsilon, pressure, Bvecx, Bvecy, Bvecz, Bvecsq, &
                                     w_lorentz, gxx, gxy, gxz, gyy, gyz, gzz, uxx, &
                                     uxy, uxz, uyy, uyz, uzz, sdetg, epsnegative, &
                                     retval)
end subroutine

subroutine Prim2ConGenWrapper(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, dens,&
                              sx, sy, sz, tau, rho, velx, vely, velz, epsilon,&
                              press, w_lorentz)

  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_REAL, INTENT(OUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(IN) :: rho, velx, vely, velz, epsilon
  CCTK_REAL, INTENT(OUT) :: press, w_lorentz

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call prim2con(handle, gxx, gxy, gxz, gyy, gyz, gzz, sdetg, dens,&
                sx, sy, sz, tau, rho, velx, vely, velz, epsilon,&
                press, w_lorentz)
end subroutine

subroutine Prim2ConPolyWrapper(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, dens,&
                               sx, sy, sz, tau, rho, velx, vely, velz, epsilon,&
                               press, w_lorentz)

  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_REAL, INTENT(OUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(IN) :: rho, velx, vely, velz
  CCTK_REAL, INTENT(OUT) :: epsilon, press, w_lorentz

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call prim2conpolytype(handle, gxx, gxy, gxz, gyy, gyz, gzz, sdetg, dens,&
                        sx, sy, sz, tau, rho, velx, vely, velz, epsilon,&
                        press, w_lorentz)
end subroutine

subroutine Prim2ConGenMWrapper(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, dens,&
                               sx, sy, sz, tau, Bconsx, Bconsy, Bconsz, rho,&
                               velx, vely, velz, epsilon, press, Bvecx, Bvecy,&
                               Bvecz, w_lorentz)

  use Con2PrimM_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_REAL, INTENT(OUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(OUT) :: Bconsx, Bconsy, Bconsz
  CCTK_REAL, INTENT(IN) :: rho, velx, vely, velz, epsilon
  CCTK_REAL, INTENT(OUT) :: press
  CCTK_REAL, INTENT(IN) :: Bvecx, Bvecy, Bvecz
  CCTK_REAL, INTENT(OUT) :: w_lorentz

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call prim2conM(handle, gxx, gxy, gxz, gyy, gyz, gzz, sdetg, dens,&
                 sx, sy, sz, tau, Bconsx, Bconsy, Bconsz, rho,&
                 velx, vely, velz, epsilon, press, Bvecx, Bvecy,&
                 Bvecz, w_lorentz)
end subroutine

subroutine Prim2ConGenM_hotWrapper(handle, GRHydro_reflevel, i, j, k, x, y, z,&
                                   gxx, gxy, gxz, gyy, gyz, gzz, det, dens, sx,&
                                   sy, sz, tau, Bconsx, Bconsy, Bconsz, rho,&
                                   velx, vely, velz, epsilon, press, Bvecx,&
                                   Bvecy, Bvecz, w_lorentz, temperature, Y_e)

  use Con2PrimM_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle, GRHydro_reflevel
  CCTK_INT, INTENT(IN) :: i, j, k
  CCTK_REAL, INTENT(IN) :: x, y, z
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_REAL, INTENT(OUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(OUT) :: Bconsx, Bconsy, Bconsz
  CCTK_REAL, INTENT(IN) :: rho, velx, vely, velz, epsilon
  CCTK_REAL, INTENT(OUT) :: press
  CCTK_REAL, INTENT(IN) :: Bvecx, Bvecy, Bvecz
  CCTK_REAL, INTENT(OUT) :: w_lorentz
  CCTK_REAL, INTENT(INOUT) :: temperature
  CCTK_REAL, INTENT(IN) :: Y_e

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call prim2conM_hot(handle, GRHydro_reflevel, i, j, k, x, y, z, gxx,&
                     gxy, gxz, gyy, gyz, gzz, sdetg, dens, sx, sy, sz,&
                     tau, Bconsx, Bconsy, Bconsz, rho, velx, vely, velz,&
                     epsilon, press, Bvecx, Bvecy, Bvecz, w_lorentz,&
                     temperature, Y_e)
end subroutine

subroutine Prim2ConPolyMWrapper(handle, gxx, gxy, gxz, gyy, gyz, gzz, det,&
                                dens, sx, sy, sz, tau, Bconsx, Bconsy, Bconsz,&
                                rho, velx, vely, velz, epsilon, press, Bvecx,&
                                Bvecy, Bvecz, w_lorentz)

  use Con2PrimM_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_REAL, INTENT(OUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(OUT) :: Bconsx, Bconsy, Bconsz
  CCTK_REAL, INTENT(IN) :: rho, velx, vely, velz
  CCTK_REAL, INTENT(OUT) :: epsilon
  CCTK_REAL, INTENT(OUT) :: press
  CCTK_REAL, INTENT(IN) :: Bvecx, Bvecy, Bvecz
  CCTK_REAL, INTENT(OUT) :: w_lorentz

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call prim2conpolytypeM(handle, gxx, gxy, gxz, gyy, gyz, gzz, sdetg,&
                         dens, sx, sy, sz, tau, Bconsx, Bconsy, Bconsz,&
                         rho, velx, vely, velz, epsilon, press, Bvecx,&
                         Bvecy, Bvecz, w_lorentz)
end subroutine

subroutine Con2PrimGenMeeWrapper(handle, keytemp, prec, gamma_eos, dens, sx, sy,&
                                 sz, tau, Bconsx, Bconsy, Bconsz, entropycons,&
                                 y_e, temp, rho, velx, vely, velz, epsilon,&
                                 pressure, Bvecx, Bvecy, Bvecz, Bvecsq,&
                                 w_lorentz, gxx, gxy, gxz, gyy, gyz, gzz, uxx,&
                                 uxy, uxz, uyy, uyz, uzz, det, epsnegative,&
                                 retval)

  use Con2Prim_fortran_interfaces

  implicit none

  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, INTENT(IN) :: handle, keytemp
  CCTK_REAL, INTENT(IN) :: prec, gamma_eos
  CCTK_REAL, INTENT(INOUT) :: dens, sx, sy, sz, tau
  CCTK_REAL, INTENT(IN) :: Bconsx, Bconsy, Bconsz
  CCTK_REAL, INTENT(INOUT) :: entropycons, y_e, temp
  CCTK_REAL, INTENT(INOUT) :: rho, velx, vely, velz, epsilon, pressure
  CCTK_REAL, INTENT(OUT) :: Bvecx, Bvecy, Bvecz
  CCTK_REAL, INTENT(OUT) :: Bvecsq
  CCTK_REAL, INTENT(INOUT) :: w_lorentz
  CCTK_REAL, INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz
  CCTK_REAL, INTENT(IN) :: uxx, uxy, uxz, uyy, uyz, uzz
  CCTK_REAL, INTENT(IN) :: det
  CCTK_INT, INTENT(OUT) :: epsnegative
  CCTK_REAL, INTENT(OUT) :: retval

  CCTK_REAL :: sdetg

  sdetg = sqrt(det)

  call GRHydro_Con2PrimM_ptee(handle, keytemp, prec, gamma_eos, dens, sx, sy, sz,&
                              tau, Bconsx, Bconsy, Bconsz, entropycons, y_e, temp,&
                              rho, velx, vely, velz, epsilon, pressure, Bvecx,&
                              Bvecy, Bvecz, Bvecsq, w_lorentz, gxx, gxy, gxz, gyy,&
                              gyz, gzz, uxx, uxy, uxz, uyy, uyz, uzz, sdetg,&
                              epsnegative, retval)
end subroutine
