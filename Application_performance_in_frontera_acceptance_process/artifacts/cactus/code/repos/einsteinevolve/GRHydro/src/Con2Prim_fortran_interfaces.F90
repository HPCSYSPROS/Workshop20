#include "cctk.h" 

module Con2Prim_fortran_interfaces
  implicit none
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  interface

     subroutine Con2Prim_pt(&
          cctk_iteration,ii,jj,kk,&
          handle, &
          dens, &
          sx, sy, sz, &
          tau, &
          rho, &
          velx, vely, velz, &
          epsilon, press, &
          w_lorentz, &
          uxx, uxy, uxz, &
          uyy, uyz, uzz, &
          det, &
          x, y, z, r, &
          epsnegative, &
          GRHydro_rho_min, pmin, epsmin, &
          GRHydro_reflevel, GRHydro_C2P_failed)
       
       implicit none
       CCTK_INT  cctk_iteration, ii, jj, kk
       CCTK_INT  handle
       CCTK_REAL dens
       CCTK_REAL sx, sy, sz
       CCTK_REAL tau
       CCTK_REAL rho 
       CCTK_REAL velx, vely, velz
       CCTK_REAL epsilon, press
       CCTK_REAL w_lorentz
       CCTK_REAL uxx, uxy, uxz
       CCTK_REAL uyy, uyz, uzz
       CCTK_REAL det
       CCTK_REAL x, y, z, r
       logical  epsnegative
       CCTK_REAL GRHydro_rho_min, pmin, epsmin
       CCTK_INT  GRHydro_reflevel
       CCTK_REAL GRHydro_C2P_failed    
     end subroutine Con2Prim_pt

     subroutine Con2Prim_ptPolytype(GRHydro_polytrope_handle, &
          dens, &
          sx, sy, sz, &
          tau, &
          rho, &
          velx, vely, velz, &
          eps, press, &
          w_lorentz, &
          uxx, uxy, uxz, uyy, uyz, uzz, &
          det, &
          x, y, z, r, &
          GRHydro_rho_min, &
          GRHydro_reflevel, GRHydro_C2P_failed)
          
       implicit none
       CCTK_INT  GRHydro_polytrope_handle
       CCTK_REAL dens
       CCTK_REAL sx, sy, sz
       CCTK_REAL tau
       CCTK_REAL rho 
       CCTK_REAL velx, vely, velz
       CCTK_REAL eps, press
       CCTK_REAL w_lorentz
       CCTK_REAL uxx, uxy, uxz
       CCTK_REAL uyy, uyz, uzz
       CCTK_REAL det
       CCTK_REAL x, y, z, r
       CCTK_REAL GRHydro_rho_min
       CCTK_INT  GRHydro_reflevel
       CCTK_REAL GRHydro_C2P_failed
     end subroutine Con2Prim_ptPolytype

     subroutine Con2Prim_ptTracer(cons_tracer, tracer, dens)
       implicit none  
       CCTK_REAL cons_tracer, tracer, dens
     end subroutine Con2Prim_ptTracer

     subroutine Con2Prim_ptBoundsTracer(cons_tracer, tracer, rho, one_over_w_lorentz, det)
       implicit none
       CCTK_REAL cons_tracer, tracer, rho, one_over_w_lorentz, det
     end subroutine Con2Prim_ptBoundsTracer

  end interface

end module Con2Prim_fortran_interfaces
