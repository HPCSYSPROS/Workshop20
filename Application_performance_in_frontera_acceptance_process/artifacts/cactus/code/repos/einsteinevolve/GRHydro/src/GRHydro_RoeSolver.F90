 /*@@
   @file      GRHydro_RoeSolver.F90
   @date      Sat Jan 26 01:55:27 2002
   @author    
   @desc 
   Calculates the Roe flux
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#include "SpaceMask.h"

 /*@@
   @routine    GRHydro_RoeSolve
   @date       Sat Jan 26 01:55:55 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Wrapper routine to calculate the Roe fluxes and hence the update 
   terms.
   @enddesc 
   @calls     
   @calledby   
   @history 
   
   @endhistory 

@@*/

subroutine GRHydro_RoeSolve(CCTK_ARGUMENTS)
  
  USE GRHydro_Eigenproblem
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  CCTK_REAL, dimension(5) :: roeflux,roeave,qdiff,consp,consm_i,&
       fplus,fminus,f_roe,primp,primm_i
  CCTK_REAL :: avg_alp,avg_beta,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
       avg_det,uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,&
       rhoave, velxave, velyave, velzave, epsave, &
       w_lorentzave, usendh
  integer :: m
  integer :: i,j,k
  
  CCTK_INT :: type_bits, trivial

  if (flux_direction == 1) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemX")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemX", &
         &"trivial")
  else if (flux_direction == 2) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemY")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemY", &
         &"trivial")
  else if (flux_direction == 3) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemZ")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemZ", &
         &"trivial")
  else
    call CCTK_ERROR("Flux direction not x,y,z")
    STOP
  end if
  
  f_roe = 0.d0
  
  do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil

!!$        Set the left (p for plus) and right (m_i for minus, i+1) states
        
        consp(1)   = densplus(i,j,k) 
        consp(2)   = sxplus(i,j,k)
        consp(3)   = syplus(i,j,k)
        consp(4)   = szplus(i,j,k)
        consp(5)   = tauplus(i,j,k)
        
        consm_i(1) = densminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(2) = sxminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(3) = syminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(4) = szminus(i+xoffset,j+yoffset,k+zoffset)
        consm_i(5) = tauminus(i+xoffset,j+yoffset,k+zoffset) 
        
        primp(1)   = rhoplus(i,j,k) 
        primp(2)   = velxplus(i,j,k)
        primp(3)   = velyplus(i,j,k)
        primp(4)   = velzplus(i,j,k)
        primp(5)   = epsplus(i,j,k)
        
        primm_i(1) = rhominus(i+xoffset,j+yoffset,k+zoffset)
        primm_i(2) = velxminus(i+xoffset,j+yoffset,k+zoffset)
        primm_i(3) = velyminus(i+xoffset,j+yoffset,k+zoffset)
        primm_i(4) = velzminus(i+xoffset,j+yoffset,k+zoffset)
        primm_i(5) = epsminus(i+xoffset,j+yoffset,k+zoffset) 
        
        roeflux = 0.d0
        qdiff = 0.d0

!!$        Calculate jumps in conserved variables
        
        do m = 1,5
          qdiff(m) = consm_i(m) - consp(m)
        end do

!!$        Set metric terms at interface
        
        if (flux_direction == 1) then
           avg_beta = 0.5d0 * (betax(i+xoffset,j+yoffset,k+zoffset) + &
                betax(i,j,k))
        else if (flux_direction == 2) then
           avg_beta = 0.5d0 * (betay(i+xoffset,j+yoffset,k+zoffset) + &
                 betay(i,j,k))
        else if (flux_direction == 3) then
           avg_beta = 0.5d0 * (betaz(i+xoffset,j+yoffset,k+zoffset) + &
                betaz(i,j,k))
        else
           call CCTK_ERROR("Flux direction not x,y,z")
           STOP
        end if
        
        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gxxh = 0.5d0 * (gxx(i+xoffset,j+yoffset,k+zoffset) + &
             gxx(i,j,k))
        gxyh = 0.5d0 * (gxy(i+xoffset,j+yoffset,k+zoffset) + &
             gxy(i,j,k))
        gxzh = 0.5d0 * (gxz(i+xoffset,j+yoffset,k+zoffset) + &
             gxz(i,j,k))
        gyyh = 0.5d0 * (gyy(i+xoffset,j+yoffset,k+zoffset) + &
             gyy(i,j,k))
        gyzh = 0.5d0 * (gyz(i+xoffset,j+yoffset,k+zoffset) + &
             gyz(i,j,k))
        gzzh = 0.5d0 * (gzz(i+xoffset,j+yoffset,k+zoffset) + &
             gzz(i,j,k))
        
        avg_det = SPATIAL_DETERMINANT(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh)
        
!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell
          
        if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then

          if (flux_direction == 1) then
            call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                 f_roe(1),f_roe(2),f_roe(3),&
                 f_roe(4),f_roe(5),&
                 velxplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 2) then
            call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                 f_roe(1),f_roe(3),f_roe(4),&
                 f_roe(2),f_roe(5),&
                 velyplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 3) then
            call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                 f_roe(1),f_roe(4),f_roe(2),&
                 f_roe(3),f_roe(5),&
                 velzplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if
            
        else !!! The end of this branch is right at the bottom of the routine

          call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
               avg_det,gxxh, gxyh, gxzh, gyyh, gyzh, gzzh)
        
          if (flux_direction == 1) then
            usendh = uxxh
          else if (flux_direction == 2) then
            usendh = uyyh
          else if (flux_direction == 3) then
            usendh = uzzh
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if
          
!!$        Set the Roe average of the fluid variables
        
          call roeaverage(primp, primm_i, roeave)
        
          rhoave = roeave(1)
          velxave = roeave(2)
          velyave = roeave(3)
          velzave = roeave(4)
          epsave = roeave(5)

!!$        Convert to conserved variables and find the part of the Roe
!!$        flux that requires the spectral decomposition.
!!$        The conversion to conserved variables is just to set the 
!!$        pressure at this point (means this routine doesn''t need
!!$        the EOS interface).

!!$        The conversion routine is unnecessary (the pressure is set
!!$        inside the eigenproblem routine) so instead we just have
!!$        to set the average W.

          w_lorentzave = 1.d0 / &
               sqrt(1.d0 - &
                 (gxxh*velxave*velxave + gyyh*velyave*velyave + &
                  gzzh*velzave*velzave + 2*gxyh*velxave*velyave + &
                  2*gxzh*velxave *velzave + 2*gyzh*velyave*velzave))  

          if (flux_direction == 1) then
!!$            call prim2con(GRHydro_eos_handle,gxxh, gxyh, gxzh, gyyh, &
!!$                 gyzh, gzzh, avg_det, &
!!$                 consh(1), consh(2), consh(3), consh(4), consh(5), rhoave, &
!!$                 velxave, velyave, velzave, epsave, pressave, w_lorentzave)
            call eigenproblem(GRHydro_eos_handle,rhoave, velxave, &
                 velyave, velzave, epsave, w_lorentzave, &
                 gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
                 usendh, avg_alp,avg_beta,qdiff(1),qdiff(2), &
                 qdiff(3),qdiff(4),qdiff(5),roeflux(1),roeflux(2),&
                 roeflux(3),roeflux(4),roeflux(5))
          else if (flux_direction == 2) then
!!$            call prim2con(GRHydro_eos_handle,gyyh, gyzh, gxyh, gzzh, &
!!$                 gxzh, gxxh, avg_det, &
!!$                 consh(1), consh(3), consh(4), consh(2), consh(5), rhoave, &
!!$                 velyave, velzave, velxave, epsave, pressave, w_lorentzave)
            call eigenproblem(GRHydro_eos_handle,rhoave, velyave, &
                 velzave, velxave, epsave, w_lorentzave, &
                 gyyh,gyzh,gxyh,gzzh,gxzh,gxxh, &
                 usendh, avg_alp,avg_beta,qdiff(1),qdiff(3), &
                 qdiff(4),qdiff(2),qdiff(5),roeflux(1),roeflux(3),&
                 roeflux(4),roeflux(2),roeflux(5))
          else if (flux_direction == 3) then
!!$            call prim2con(GRHydro_eos_handle,gzzh, gxzh, gyzh, gxxh, &
!!$                 gxyh, gyyh, avg_det, &
!!$                 consh(1), consh(4), consh(2), consh(3), consh(5), rhoave, &
!!$                 velzave, velxave, velyave, epsave, pressave, w_lorentzave)
            call eigenproblem(GRHydro_eos_handle,rhoave, velzave, &
                 velxave, velyave, epsave, w_lorentzave, &
                 gzzh,gxzh,gyzh,gxxh,gxyh,gyyh, &
                 usendh, avg_alp,avg_beta,qdiff(1),qdiff(4), &
                 qdiff(2),qdiff(3),qdiff(5),roeflux(1),roeflux(4),&
                 roeflux(2),roeflux(3),roeflux(5))
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if
          
          fplus = 0.d0
          fminus = 0.d0
        
!!$Calculate the fluxes of the original reconstructed data
        
          if (flux_direction == 1) then
            call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5), &
                 fplus(1),fplus(2),fplus(3),fplus(4), &
                 fplus(5),velxplus(i,j,k),pressplus(i,j,k), &
                 avg_det,avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(2),consm_i(3), &
                 consm_i(4),consm_i(5),fminus(1),fminus(2),fminus(3), &
                 fminus(4), fminus(5), &
                 velxminus(i+xoffset,j+yoffset,k+zoffset), &
                 pressminus(i+xoffset,j+yoffset,k+zoffset), &
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 2) then
            call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5), &
                 fplus(1),fplus(3),fplus(4),fplus(2), &
                 fplus(5),velyplus(i,j,k),pressplus(i,j,k), &
                 avg_det,avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(3),consm_i(4), &
                 consm_i(2),consm_i(5),fminus(1),fminus(3),fminus(4), &
                 fminus(2), fminus(5), &
                 velyminus(i+xoffset,j+yoffset,k+zoffset), &
                 pressminus(i+xoffset,j+yoffset,k+zoffset), &
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 3) then
            call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5), &
                 fplus(1),fplus(4),fplus(2),fplus(3), &
                 fplus(5),velzplus(i,j,k),pressplus(i,j,k),avg_det, &
                 avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(4),consm_i(2), &
                 consm_i(3),consm_i(5),fminus(1),fminus(4),fminus(2), &
                 fminus(3), fminus(5), &
                 velzminus(i+xoffset,j+yoffset,k+zoffset), &
                 pressminus(i+xoffset,j+yoffset,k+zoffset), &
                 avg_det,avg_alp,avg_beta)
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if
        
!!$The combined Roe flux          
        
          do m = 1,5
          
            f_roe(m) = 0.5d0 * (fplus(m) + fminus(m) - roeflux(m))
          
          end do

        end if !!! The end of the SpaceMask check for a trivial RP.
        
        densflux(i,j,k) = f_roe(1)
        sxflux(i,j,k)   = f_roe(2)
        syflux(i,j,k)   = f_roe(3)
        szflux(i,j,k)   = f_roe(4)
        tauflux(i,j,k)  = f_roe(5)

      enddo
    enddo
  enddo
    
end subroutine GRHydro_RoeSolve

