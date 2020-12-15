 /*@@
   @file      GRHydro_ParamCheck.F90
   @date      Sat Feb  9 23:48:01 2002
   @author    
   @desc 
   Parameter checking routine.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydro_ParamCheck
   @date       Sat Feb  9 23:48:43 2002
   @author     Ian Hawke
   @desc 
   Checks the parameters.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_ParamCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: coordinates_is_active, evolution_mask_idx

  if (GRHydro_stencil > minval(cctk_nghostzones)) then
    call CCTK_PARAMWARN("The stencil is larger than the number of ghost zones. Answer will be dependent on processor number...")
  end if

  if (CCTK_EQUALS(recon_method,"tvd").and.(GRHydro_stencil < 2)) then
    call CCTK_PARAMWARN("The stencil size must be at least 2 to use TVD reconstruction")
  end if
  
  if (CCTK_EQUALS(recon_method,"ppm")) then
    if (CCTK_EQUALS(ppm_flatten,"stencil_3").and.(GRHydro_stencil < 3)) then
      call CCTK_PARAMWARN("The stencil size must be at least 3 to use PPM reconstruction with the stencil-3 variant of the flattening procedure")
    else if (CCTK_EQUALS(ppm_flatten,"stencil_4").and.(GRHydro_stencil < 4)) then
      call CCTK_PARAMWARN("The stencil size must be at least 4 to use PPM reconstruction with the stencil-4 (original) flattening procedure")
    end if
  end if

  if (CCTK_EQUALS(recon_method,"weno").and.(GRHydro_stencil < 3)) then
    call CCTK_PARAMWARN("The stencil size must be at least 3 to use WENO reconstruction.")
  end if

  if (CCTK_EQUALS(recon_method,"mp5").and.(GRHydro_stencil < 3)) then
    call CCTK_PARAMWARN("The stencil size must be at least 3 to use MP5 reconstruction.")
  end if

  if (CCTK_EQUALS(recon_method,"eno").and.(GRHydro_stencil < eno_order)) then
    call CCTK_PARAMWARN("The stencil size must be at least the order of the reconstruction to use ENO reconstruction")
  end if

  if (CCTK_EQUALS(GRHydro_eos_table,"2D_Polytrope").and.&
       (.not.CCTK_EQUALS(GRHydro_eos_type,"Polytype"))) then
    call CCTK_PARAMWARN("When using the 2D_Polytrope EOS you need to set eos_type to Polytype")
  end if

  if (CCTK_EQUALS(GRHydro_eos_table,"Ideal_Fluid").and.&
       (.not.CCTK_EQUALS(GRHydro_eos_type,"General"))) then
    call CCTK_PARAMWARN("When using the Ideal_Fluid EOS you need to set eos_type to General")
  end if

  if (.not.(CCTK_EQUALS(metric_type, "Physical").or.&
       CCTK_EQUALS(metric_type, "Static Conformal"))) then
    call CCTK_ERROR("GRHydro only knows how to deal with physical metric type at the moment. Complain to the maintainers...")
    STOP
  end if
  
  
  if (use_mask .eq. 0) then
    call CCTK_PARAMWARN("GRHydro now requires you to set SpaceMask::use_mask = ""yes""")
  end if

  if (number_of_particles .gt. 0) then
    if (number_of_arrays .ne. 3) then
      call CCTK_PARAMWARN("If tracking particles then you must have number_of_arrays = 3")
    end if
  end if
  ! This check and 'static' as valid keyword option can be removed after the
  ! year 2010
  if (CCTK_EQUALS(bound,"static")) then
    call CCTK_PARAMWARN("GRHydro::bound = 'static' is no longer supported, use 'none' instead");
  end if

  if (CCTK_EQUALS(Bvec_evolution_method,"GRHydro").or.CCTK_EQUALS(Bvec_evolution_method,"GRHydro_Avec")) then
     evolve_MHD = 1
  else
     evolve_MHD = 0
  endif

  if ( CCTK_EQUALS(Bvec_evolution_method,"GRHydro_Avec") ) then 
     if ( CCTK_EQUALS(Avec_gauge,"Lorenz") ) then
        evolve_Lorenz_gge = 1
     else
        evolve_Lorenz_gge = 0
     endif
  endif

  if (CCTK_EQUALS(Y_e_evolution_method,"GRHydro")) then
     evolve_Y_e = 1
  else
     evolve_Y_e = 0
  endif

  if (CCTK_EQUALS(temperature_evolution_method,"GRHydro")) then
     evolve_temper = 1
  else
     evolve_temper = 0
  endif

  if (CCTK_EQUALS(entropy_evolution_method,"GRHydro")) then
     evolve_entropy = 1
  else
     evolve_entropy = 0
  endif

  call CCTK_IsImplementationActive(coordinates_is_active, "Coordinates")
  ! this test is somewhat overzealous but we cannot access
  ! coordinates::general_coordinates yet since it is only set in BaseGrid
  if (CCTK_Equals(riemann_solver,"HLLE").eq.0.and.coordinates_is_active.ne.0) then
    call CCTK_ERROR("There is currently no Riemann solver other than HLLE compatible with multipatch!")
    STOP
  end if
  
  if (constrain_to_1D .eq. 1 .and. coordinates_is_active.eq.1) then
    call CCTK_ERROR("Constrain to 1D option does not work with other than Cartesian coordinates. Ask Christian R. to implement it if required!")
    STOP
  endif

  if (CCTK_EQUALS(riemann_solver,"Roe").and.CCTK_EQUALS(Bvec_evolution_method,"GRHydro")) then
    call CCTK_PARAMWARN("Roe solver is not implemented yet for MHD")
  end if

  if (CCTK_EQUALS(riemann_solver,"Marquina").and.CCTK_EQUALS(Bvec_evolution_method,"GRHydro")) then
    call CCTK_PARAMWARN("Marquina solver is not implemented yet for MHD")
  end if

  
  if(evolve_MHD.ne.0 .and. clean_divergence.ne.0) then
    if (CCTK_EQUALS(psidcspeed,"char speed")) then
      whichpsidcspeed = 0
    else if (CCTK_EQUALS(psidcspeed,"light speed")) then
      whichpsidcspeed = 1
    else if (CCTK_EQUALS(psidcspeed,"set speed")) then
      whichpsidcspeed = 2
    else
      call CCTK_PARAMWARN("Unknown type of speed setting for psidc (psidcspeed)")
    end if
  end if

  if(use_cxx_code.ne.0 .and. evolve_MHD.ne.0) then
    ! C++ code supports only a subset of MHD DC options.
    if(clean_divergence.ne.0) then
      if(whichpsidcspeed.ne.1) then
        call CCTK_PARAMWARN("C++ code only supports psidcspeed == 'light speed'")
      end if
      if(decouple_normal_Bfield.eq.0) then
        call CCTK_PARAMWARN("C++ code requires decouple_normal_Bfield == true")
      end if
    end if
    ! this really should go into RiemannSolverM.F90
    if(.not. CCTK_EQUALS(riemann_solver,"HLLE")) then
        call CCTK_PARAMWARN("C++ code only supports riemnna_solver == 'HLLE'")
    end if
  end if

  if(CCTK_EQUALS(use_evolution_mask, "always")) then
    call CCTK_VarIndex(evolution_mask_idx, "CarpetEvolutionMask::evolution_mask")
    if(evolution_mask_idx .lt.  0) then
      call CCTK_PARAMWARN("You set use_evolution_mask='always' but I cannot find 'CarpetEvolutionMask::evolution_mask'. If you use Carpet, then you should activate thorn 'CarpetEvolutionMask', if using PUGH then evolution_mask makes no sense and you should disable this option.")
    end if
  end if

  if(clean_divergence.ne.0 .and. evolve_MHD.eq.0) then
    call CCTK_WARN(1, "You activated divergence cleaning but do not evolve a magnetic field. At best, nothing should happen, be preapered for a segfault otherwise.")
  end if

  if (Tmunu_damping_radius_min .gt. Tmunu_damping_radius_max) then
     call CCTK_ERROR("Minimum damping radius is greater than maximum damping radius!")
     STOP
  end if

  if(evolve_entropy.ne.0) then
    if(evolve_MHD.eq.0) then
      call CCTK_PARAMWARN("Entropy evolution not properly implemented yet for unmagnetized fluids!")
    end if
    if (CCTK_EQUALS(recon_method,"ppm")) then
      call CCTK_PARAMWARN("Entropy evolution not implemented yet with PPM reconstruction!")
    end if 
  end if

end subroutine GRHydro_ParamCheck

