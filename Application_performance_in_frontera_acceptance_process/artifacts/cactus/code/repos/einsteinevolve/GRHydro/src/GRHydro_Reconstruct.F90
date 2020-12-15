 /*@@
   @file      GRHydro_Reconstruct.F90
   @date      Sat Jan 26 02:13:25 2002
   @author    Bruno Mundim, Josh Faber, Christian D. Ott
   @desc 
   Wrapper routine to perform the reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "GRHydro_Macros.h"
#include "SpaceMask.h"

subroutine Reconstruction(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT  :: i,j,k
  CCTK_REAL :: local_min_tracer, dummy1, dummy2
  CCTK_INT :: type_bits, not_trivial
  CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz, w

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if
  
  ! Initialize plus and minus states
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           ! must initialize rho and eps plus minus
           ! to cell average in order to have sane values
           ! in the boundary zones for EOS call
           rhoplus(i,j,k) = rho(i,j,k)
           rhominus(i,j,k)= rho(i,j,k)
           epsplus(i,j,k) = eps(i,j,k)
           epsminus(i,j,k) = eps(i,j,k)
           velxplus(i,j,k) = 0.0d0
           velxminus(i,j,k) = 0.0d0
           velyplus(i,j,k) = 0.0d0
           velyminus(i,j,k) = 0.0d0
           velzplus(i,j,k) = 0.0d0
           velzminus(i,j,k) = 0.0d0

           if(evolve_mhd.ne.0) then
              Bvecxplus(i,j,k) = 0.0d0
              Bvecxminus(i,j,k) = 0.0d0
              Bvecyplus(i,j,k) = 0.0d0
              Bvecyminus(i,j,k) = 0.0d0
              Bveczplus(i,j,k) = 0.0d0
              Bveczminus(i,j,k) = 0.0d0
              if(clean_divergence.ne.0) then
                 psidcplus(i,j,k) = 0.0d0
                 psidcminus(i,j,k) = 0.0d0
              endif
           endif

           if (evolve_entropy .ne. 0) then
              entropyplus(i,j,k) = 0.0d0
              entropyminus(i,j,k) = 0.0d0
           endif

           if (evolve_tracer .ne. 0) then
              tracerplus(i,j,k,:) = 0.0d0
              tracerminus(i,j,k,:) = 0.0d0
           endif

           if (evolve_Y_e .ne. 0) then
              ! set this to the cell center values                              
              ! to make sure we have good Y_e even in                           
              ! the boundary region (for full GF EOS calls)                     
              Y_e_plus(i,j,k) = Y_e(i,j,k)
              Y_e_minus(i,j,k) = Y_e(i,j,k)
           endif

          if(evolve_temper .ne. 0) then
              ! set this to cell center value to have                           
              ! good initial guesses at interfaces                              
              ! in case we don't reconstruct temp                               
              tempplus(i,j,k) = temperature(i,j,k)
              tempminus(i,j,k) = temperature(i,j,k)
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  if (CCTK_EQUALS(recon_method,"tvd")) then
    ! this handles MHD and non-MHD
    call GRHydro_TVDReconstruct_drv(CCTK_PASS_FTOF)

  else if (CCTK_EQUALS(recon_method,"ppm")) then

    if(use_optimized_ppm.eq.0) then
      if(evolve_mhd.ne.0) then
        call GRHydro_PPMMReconstruct_drv(CCTK_PASS_FTOF)
      else
        call GRHydro_PPMReconstruct_drv(CCTK_PASS_FTOF)
      end if 
    else 
      call GRHydro_PPMReconstruct_drv_opt(cctkGH)
    end if

  else if (CCTK_EQUALS(recon_method,"eno")) then
    ! this handles MHD and non-MHD
    call GRHydro_ENOReconstruct_drv(CCTK_PASS_FTOF)

  else if (CCTK_EQUALS(recon_method,"weno") .or. CCTK_EQUALS(recon_method,"weno-z")) then
    ! this handles MHD and non-MHD
    call GRHydro_WENOReconstruct_drv(CCTK_PASS_FTOF)
  
  else if (CCTK_EQUALS(recon_method,"mp5")) then
    ! this handles MHD and non-MHD
    call GRHydro_MP5Reconstruct_drv(CCTK_PASS_FTOF)
  
  else
    call CCTK_ERROR("Reconstruction method not recognized!")
    STOP
  
  end if

  if (evolve_tracer .ne. 0) then
    if (use_min_tracer .ne. 0) then
      local_min_tracer = min_tracer
    else
      local_min_tracer = 0.0d0
    end if
  else
    ! pacify compiler warning about unintialized local_min_tracer
    local_min_tracer = 1d42
  end if


  if (flux_direction == 1) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemX")
    call SpaceMask_GetStateBits(not_trivial, "Hydro_RiemannProblemX", &
         &"not_trivial")
  else if (flux_direction == 2) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemY")
    call SpaceMask_GetStateBits(not_trivial, "Hydro_RiemannProblemY", &
         &"not_trivial")
  else if (flux_direction == 3) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemZ")
    call SpaceMask_GetStateBits(not_trivial, "Hydro_RiemannProblemZ", &
         &"not_trivial")
  else
    call CCTK_ERROR("Flux direction not x,y,z")
    STOP
  end if

  !$OMP PARALLEL DO PRIVATE(i,j,k, agxx, agxy, agxz, agyy, agyz, agzz, w, dummy1, dummy2)
  do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 1
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 1
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 1
         SET_ATMO_MIN(dummy2, GRHydro_rho_min, r(i,j,k))
         if(rhoplus(i,j,k).lt.dummy2 .or. &
            rhominus(i,j,k).lt.dummy2) then
            !.or. epsplus(i,j,k) .lt. 0.0d0 .or. epsminus(i,j,k) .lt. 0.0d0) then
                 rhoplus(i,j,k)  = rho(i,j,k)
                 rhominus(i,j,k) = rho(i,j,k)
                 epsplus(i,j,k)  = eps(i,j,k)
                 epsminus(i,j,k) = eps(i,j,k)
                 if (reconstruct_Wv.ne.0) then
                    ! divide out the Loretnz factor for both the
                    ! plus and minus quantities this should by construction ensure
                    ! that any Lorentz factor calculated from them later on is
                    ! physical (ie. > 1.d0)
                    velxplus(i,j,k) = w_lorentz(i,j,k)*vup(i,j,k,1)
                    velyplus(i,j,k) = w_lorentz(i,j,k)*vup(i,j,k,2)
                    velzplus(i,j,k) = w_lorentz(i,j,k)*vup(i,j,k,3)
                    agxx = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
                    agxy = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
                    agxz = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
                    agyy = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
                    agyz = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
                    agzz = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))
                    w = sqrt( 1.d0 + agxx*velxplus(i,j,k)*velxplus(i,j,k) &
                                   + agyy*velyplus(i,j,k)*velyplus(i,j,k) &
                                   + agzz*velzplus(i,j,k)*velzplus(i,j,k) &
                                   + 2.d0*agxy*velxplus(i,j,k)*velyplus(i,j,k) &
                                   + 2.d0*agxz*velxplus(i,j,k)*velzplus(i,j,k) &
                                   + 2.d0*agyz*velyplus(i,j,k)*velzplus(i,j,k) )
                    velxplus(i,j,k) = velxplus(i,j,k)/w
                    velyplus(i,j,k) = velyplus(i,j,k)/w
                    velzplus(i,j,k) = velzplus(i,j,k)/w
 
                    velxminus(i,j,k) = w_lorentz(i,j,k)*vup(i,j,k,1)
                    velyminus(i,j,k) = w_lorentz(i,j,k)*vup(i,j,k,2)
                    velzminus(i,j,k) = w_lorentz(i,j,k)*vup(i,j,k,3)
                    agxx = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
                    agxy = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
                    agxz = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
                    agyy = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
                    agyz = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
                    agzz = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
                    w = sqrt( 1.d0 + agxx*velxminus(i,j,k)*velxminus(i,j,k) &
                                   + agyy*velyminus(i,j,k)*velyminus(i,j,k) &
                                   + agzz*velzminus(i,j,k)*velzminus(i,j,k) &
                                   + 2.d0*agxy*velxminus(i,j,k)*velyminus(i,j,k) &
                                   + 2.d0*agxz*velxminus(i,j,k)*velzminus(i,j,k) &
                                   + 2.d0*agyz*velyminus(i,j,k)*velzminus(i,j,k) )
                    velxminus(i,j,k) = velxminus(i,j,k)/w
                    velyminus(i,j,k) = velyminus(i,j,k)/w
                    velzminus(i,j,k) = velzminus(i,j,k)/w
                 else
                    ! This is the standard way of doing it
                    velxplus(i,j,k) = vup(i,j,k,1)
                    velyplus(i,j,k) = vup(i,j,k,2)
                    velzplus(i,j,k) = vup(i,j,k,3)
                    velxminus(i,j,k) = vup(i,j,k,1)
                    velyminus(i,j,k) = vup(i,j,k,2)
                    velzminus(i,j,k) = vup(i,j,k,3)
                 endif
                 if(evolve_y_e.ne.0) then
                    y_e_plus(i,j,k)  = y_e(i,j,k)
                    y_e_minus(i,j,k) = y_e(i,j,k)
                 endif 
                 if(evolve_tracer.ne.0) then
                    where(tracerplus(i,j,k,:).le.local_min_tracer .or. &
                         tracerminus(i,j,k,:).le.local_min_tracer)
                       tracerplus(i,j,k,:) = tracer(i,j,k,:)
                       tracerminus(i,j,k,:) = tracer(i,j,k,:)
                    end where
                 end if
              end if
              
              ! Check if epsilon became negative for ideal gas EOS.
              if (evolve_temper.eq.0) then
                 if (epsplus(i,j,k) .lt. 0.0d0) then
                    epsplus(i,j,k) = eps(i,j,k)
                 endif
                 if (epsminus(i,j,k) .lt. 0.0d0) then
                    epsminus(i,j,k) = eps(i,j,k)
                 endif
              endif
              
              ! Riemann problem might not be trivial anymore!!!!
              SpaceMask_SetStateBitsF90(space_mask, i-xoffset, j-yoffset, k-zoffset, type_bits, not_trivial)
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, not_trivial)
              
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
  
  if (CCTK_EQUALS(recon_vars,"primitive").or.&
       CCTK_EQUALS(recon_method,"ppm")) then
     if(evolve_mhd.ne.0) then
        call primitive2conservativeM(CCTK_PASS_FTOF)
     else
        call GRHydro_Primitive2Conservative_CforF(cctkGH)
     endif
  else if (CCTK_EQUALS(recon_vars,"conservative")) then
     if(evolve_mhd.ne.0) then
        call Conservative2PrimitiveBoundsM(CCTK_PASS_FTOF)
     else
        call Conservative2PrimitiveBounds(CCTK_PASS_FTOF)
     endif
  else
    call CCTK_ERROR("Variable type to reconstruct not recognized.")
    STOP
  end if

  return

end subroutine Reconstruction


