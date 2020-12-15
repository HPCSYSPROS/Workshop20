  /*@@
   @file      GRHydro_Marquina.f90
   @date      Thu Jan  11 11:03:32 2002
   @author    Pedro Montero, Toni Font    
   @desc 
   Routine to obtain the Marquina Fluxes. Note that this is the 
   MODIFIED Marquina formula as given by Aloy et.al. 
   (ApJ Supp 122 (1999) p.151) and not the full Marquina flux 
   of Donat and Marquina.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "GRHydro_Macros.h"
#include "SpaceMask.h"

 /*@@
   @routine    GRHydro_Marquina.f90 
   @date       Wed Feb 13 11:03:32 2002
   @author     Pedro Montero, Toni Font
   @desc 
   Routine to obtain the Marquina Fluxes
   @enddesc 
   @calls     
   @calledby   
   @history 
   Based on routines by Toni Font
   @endhistory 

@@*/


subroutine GRHydro_Marquina(CCTK_ARGUMENTS)
    
    implicit none


    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    CCTK_REAL, dimension(5) :: marquinaflux, &
         consp,consm_i,fplus,fminus,f_marquina,primp,primm_i
    CCTK_REAL :: avg_alp,avg_beta,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
         avg_det,uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,&
         tmp_w_lorentzp, tmp_w_lorentzm_i, w_lorentzp,w_lorentzm_i, usendh
    integer :: m
    integer :: i,j,k
    integer :: keytemp
    
    CCTK_INT :: type_bits, trivial

    if(evolve_temper.eq.1.and.reconstruct_temper.eq.1) then
       keytemp = 1
    else
       keytemp = 0
    endif

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
      !Keep this check in here, it is not checked again later
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if
    
    f_marquina = 0.d0
        
    !$OMP PARALLEL DO PRIVATE(i,j,k,consp,consm_i,primp,primm_i,&
    !$OMP marquinaflux,avg_beta,avg_alp,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
    !$OMP f_marquina,uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,usendh,&
    !$OMP tmp_w_lorentzp, tmp_w_lorentzm_i,w_lorentzp,w_lorentzm_i,&
    !$OMP fplus,fminus,m,avg_det)
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

          marquinaflux = 0.d0
        
!!$        Set metric terms at interface
          
          if (flux_direction == 1) then
             avg_beta = 0.5d0 * (betax(i+xoffset,j+yoffset,k+zoffset) + &
                  betax(i,j,k))
          else if (flux_direction == 2) then
             avg_beta = 0.5d0 * (betay(i+xoffset,j+yoffset,k+zoffset) + &
                  betay(i,j,k))
          else
             avg_beta = 0.5d0 * (betaz(i+xoffset,j+yoffset,k+zoffset) + &
                  betaz(i,j,k))
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

!!$ routine to calculate the determinant of the metric

         avg_det = SPATIAL_DETERMINANT(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh)
          
!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell

          if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then

            if (flux_direction == 1) then
              call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                   f_marquina(1),f_marquina(2),f_marquina(3),&
                   f_marquina(4),f_marquina(5),&
                   velxplus(i,j,k),pressplus(i,j,k),&
                   avg_det,avg_alp,avg_beta)
            else if (flux_direction == 2) then
              call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                   f_marquina(1),f_marquina(3),f_marquina(4),&
                   f_marquina(2),f_marquina(5),&
                   velyplus(i,j,k),pressplus(i,j,k),&
                   avg_det,avg_alp,avg_beta)
            else
              call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                   f_marquina(1),f_marquina(4),f_marquina(2),&
                   f_marquina(3),f_marquina(5),&
                   velzplus(i,j,k),pressplus(i,j,k),&
                   avg_det,avg_alp,avg_beta)
            end if
            
          else !!! The end of this branch is right at the bottom of the routine
            
            call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
                 avg_det,gxxh, gxyh, gxzh, gyyh, gyzh, gzzh)
            
            if (flux_direction == 1) then
              usendh = uxxh
            else if (flux_direction == 2) then
              usendh = uyyh
            else
              usendh = uzzh
            end if

!!$left state

            tmp_w_lorentzp = gxxh*primp(2)*primp(2) + &
                 gyyh*primp(3)*primp(3) + gzzh*primp(4)*primp(4) + &
                 2*gxyh*primp(2)*primp(3) + 2*gxzh*primp(2) *primp(4) + &
                 2*gyzh*primp(3)*primp(4)
            if (tmp_w_lorentzp .ge. 1.d0) then
              w_lorentzp = GRHydro_lorentz_overshoot_cutoff
            else
              w_lorentzp = 1.d0 / sqrt(1.d0 - tmp_w_lorentzp);
            endif


!!$right state

            tmp_w_lorentzm_i = gxxh*primm_i(2)*primm_i(2) + &
                 gyyh*primm_i(3)*primm_i(3) + gzzh*primm_i(4)*primm_i(4) + &
                 2*gxyh*primm_i(2)*primm_i(3) + &
                 2*gxzh*primm_i(2) *primm_i(4)+ &
                 2*gyzh*primm_i(3)*primm_i(4)
            if (tmp_w_lorentzm_i .ge. 1.d0) then
              w_lorentzm_i = GRHydro_lorentz_overshoot_cutoff
            else
              w_lorentzm_i = 1.d0 / sqrt(1.d0 - tmp_w_lorentzm_i);
            endif

            
!!$eigenvalues and right eigenvectors
            
            if (flux_direction == 1) then

               if(evolve_temper.eq.0) then
                  call eigenproblem_marquina(GRHydro_eos_handle,&
                       primm_i(1),primm_i(2), & 
                       primm_i(3),primm_i(4),primm_i(5),primp(1), &
                       primp(2),primp(3),primp(4),primp(5), &
                       gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
                       usendh,avg_det,avg_alp,avg_beta,consp(1),consp(2),&
                       consp(3), consp(4), consp(5),consm_i(1),consm_i(2), &
                       consm_i(3),consm_i(4),consm_i(5),marquinaflux(1), &
                       marquinaflux(2),marquinaflux(3),marquinaflux(4), &
                       marquinaflux(5))
               else
                  call eigenproblem_marquina_hot(GRHydro_eos_handle,keytemp,&
                       primm_i(1),primm_i(2), & 
                       primm_i(3),primm_i(4),primm_i(5),primp(1), &
                       primp(2),primp(3),primp(4),primp(5), &
                       tempminus(i+xoffset,j+yoffset,k+zoffset),&
                       tempplus(i,j,k),&
                       y_e_minus(i+xoffset,j+yoffset,k+zoffset),y_e_plus(i,j,k),&
                       gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
                       usendh,avg_det,avg_alp,avg_beta,consp(1),consp(2),&
                       consp(3), consp(4), consp(5),consm_i(1),consm_i(2), &
                       consm_i(3),consm_i(4),consm_i(5),marquinaflux(1), &
                       marquinaflux(2),marquinaflux(3),marquinaflux(4), &
                       marquinaflux(5))

               endif
              
            else if (flux_direction == 2) then

               if(evolve_temper.eq.0) then              
                  call eigenproblem_marquina(GRHydro_eos_handle,&
                       primm_i(1),primm_i(3), & 
                       primm_i(4),primm_i(2),primm_i(5),primp(1), &
                       primp(3),primp(4),primp(2),primp(5), &
                       gyyh,gyzh,gxyh,gzzh,gxzh,gxxh, &
                       usendh,avg_det,avg_alp,avg_beta,consp(1),consp(3),&
                       consp(4), consp(2), consp(5),consm_i(1),consm_i(3), &
                       consm_i(4),consm_i(2),consm_i(5),marquinaflux(1), &
                       marquinaflux(3),marquinaflux(4),marquinaflux(2), &
                       marquinaflux(5))
               else
                  call eigenproblem_marquina_hot(GRHydro_eos_handle,keytemp,&
                       primm_i(1),primm_i(3), & 
                       primm_i(4),primm_i(2),primm_i(5),primp(1), &
                       primp(3),primp(4),primp(2),primp(5), &
                       tempminus(i+xoffset,j+yoffset,k+zoffset),&
                       tempplus(i,j,k),&
                       y_e_minus(i+xoffset,j+yoffset,k+zoffset),y_e_plus(i,j,k),&
                       gyyh,gyzh,gxyh,gzzh,gxzh,gxxh, &
                       usendh,avg_det,avg_alp,avg_beta,consp(1),consp(3),&
                       consp(4), consp(2), consp(5),consm_i(1),consm_i(3), &
                       consm_i(4),consm_i(2),consm_i(5),marquinaflux(1), &
                       marquinaflux(3),marquinaflux(4),marquinaflux(2), &
                       marquinaflux(5))

               endif
              
            else

               if(evolve_temper.eq.0) then
                  call eigenproblem_marquina(GRHydro_eos_handle,&
                       primm_i(1),primm_i(4), & 
                       primm_i(2),primm_i(3),primm_i(5),primp(1), &
                       primp(4),primp(2),primp(3),primp(5), &
                       gzzh,gxzh,gyzh,gxxh,gxyh,gyyh, &
                       usendh,avg_det,avg_alp,avg_beta,consp(1),consp(4),&
                       consp(2), consp(3), consp(5),consm_i(1),consm_i(4), &
                       consm_i(2),consm_i(3),consm_i(5),marquinaflux(1), &
                       marquinaflux(4),marquinaflux(2),marquinaflux(3), &
                       marquinaflux(5))
               else
                  call eigenproblem_marquina_hot(GRHydro_eos_handle,keytemp,&
                       primm_i(1),primm_i(4), & 
                       primm_i(2),primm_i(3),primm_i(5),primp(1), &
                       primp(4),primp(2),primp(3),primp(5), &
                       tempminus(i+xoffset,j+yoffset,k+zoffset),&
                       tempplus(i,j,k),&
                       y_e_minus(i+xoffset,j+yoffset,k+zoffset),y_e_plus(i,j,k),&
                       gzzh,gxzh,gyzh,gxxh,gxyh,gyyh, &
                       usendh,avg_det,avg_alp,avg_beta,consp(1),consp(4),&
                       consp(2), consp(3), consp(5),consm_i(1),consm_i(4), &
                       consm_i(2),consm_i(3),consm_i(5),marquinaflux(1), &
                       marquinaflux(4),marquinaflux(2),marquinaflux(3), &
                       marquinaflux(5))
               endif
            end if
            
            fplus = 0.d0
            fminus = 0.d0
            
!!$calculate the fluxes
            
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
              
              else
              
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
              
            end if
            
!!$ Marquina flux
            
            do m = 1,5
              
              f_marquina(m) = 0.5d0 * (fplus(m) + fminus(m) - marquinaflux(m))
              
            end do

          end if !!! The end of the SpaceMask check for a trivial RP.

          densflux(i,j,k) = f_marquina(1)
          sxflux(i,j,k)   = f_marquina(2)
          syflux(i,j,k)   = f_marquina(3)
          szflux(i,j,k)   = f_marquina(4)
          tauflux(i,j,k)  = f_marquina(5)

        enddo
      enddo
    enddo
      !$OMP END PARALLEL DO

    if (evolve_tracer .ne. 0) then

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil
        do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil
          do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil

            if (densflux(i, j, k) > 0.d0) then

              cons_tracerflux(i, j, k,:) = &
                   tracerplus(i, j, k,:) * &
                   densflux(i, j, k)

            else

              cons_tracerflux(i, j, k,:) = &
                   tracerminus(i + xoffset, j + yoffset, k + zoffset,:) * &
                   densflux(i, j, k)

            end if

          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end if

    if (evolve_Y_e .ne. 0) then

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil
        do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil
          do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil

            if (densflux(i, j, k) > 0.d0) then

              Y_e_con_flux(i, j, k) = &
                   Y_e_plus(i, j, k) * &
                   densflux(i, j, k)

            else

              Y_e_con_flux(i, j, k) = &
                   Y_e_minus(i + xoffset, j + yoffset, k + zoffset) * &
                   densflux(i, j, k)

            end if

          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end if
    
    return
end subroutine GRHydro_Marquina
























