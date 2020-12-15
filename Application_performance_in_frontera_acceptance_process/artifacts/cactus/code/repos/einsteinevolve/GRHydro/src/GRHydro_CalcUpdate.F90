  /*@@
   @file      GRHydro_CalcUpdate.F90
   @date      Thu Jan  11 11:03:32 2002
   @author    Ian Hawke
   @desc 
   Calculates the update terms given the fluxes. Moved to here so that 
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

 /*@@
   @routine    UpdateCalculation 
   @date       Wed Feb 13 11:03:32 2002
   @author     Ian Hawke
   @desc 
   Calculates the update terms from the fluxes.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Moved out of the Riemann solver routines to make the FishEye /
   weighted flux calculation easier.
   @endhistory 

@@*/


subroutine UpdateCalculation(CCTK_ARGUMENTS)
    
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k,itracer
  CCTK_REAL :: idx, idy, idz, alp_l, alp_r, Bcons_l, Bcons_r, alp_tmp

  idx = 1.d0 / CCTK_DELTA_SPACE(flux_direction)

  if (CCTK_EQUALS(method_type, "RSA FV")) then

    if (use_weighted_fluxes == 0) then

      if(evolve_mhd.ne.0 .and. transport_constraints.ne.0) then
        !$OMP PARALLEL DO PRIVATE(i,j,k,alp_l,alp_r,alp_tmp)
        do k = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(3) - GRHydro_stencil ! we need to compute Evec on all faces/edges where the fluxes are defined
          do j = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(2) - GRHydro_stencil
            do i = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(1) - GRHydro_stencil

              alp_l = 0.5d0 * (alp(i,j,k) + &
                   alp(i-xoffset,j-yoffset,k-zoffset))
              alp_r = 0.5d0 * (alp(i,j,k) + &
                   alp(i+xoffset,j+yoffset,k+zoffset))
             
              ! we have to first compute all components of v\crossB = E and
              ! combine them in the last substep into Bconshs
              ! Evec lives on edges of cell: Evec(i,j,k,1) is at edge i,j+1/2,k+1/2 ie. the lower-front edge of cell (i,j,k)
              if(flux_direction.eq.1) then
                 alp_tmp = 0.5d0 * (alp(i,j,k+1) + alp(i+xoffset,j+yoffset,k+zoffset+1))
                 Evec(i,j,k,2) = Evec(i,j,k,2) + 0.25d0 * (alp_r*Bconszflux(i,j,k) + alp_tmp*Bconszflux(i  ,j  ,k+1))
                 alp_tmp = 0.5d0 * (alp(i,j+1,k) + alp(i+xoffset,j+yoffset+1,k+zoffset))
                 Evec(i,j,k,3) = Evec(i,j,k,3) - 0.25d0 * (alp_r*Bconsyflux(i,j,k) + alp_tmp*Bconsyflux(i  ,j+1,k  ))
              elseif(flux_direction.eq.2) then
                 alp_tmp = 0.5d0 * (alp(i,j,k+1) + alp(i+xoffset,j+yoffset,k+zoffset+1))
                 Evec(i,j,k,1) = Evec(i,j,k,1) - 0.25d0 * (alp_r*Bconszflux(i,j,k) + alp_tmp*Bconszflux(i  ,j  ,k+1))
                 alp_tmp = 0.5d0 * (alp(i+1,j,k) + alp(i+xoffset+1,j+yoffset,k+zoffset))
                 Evec(i,j,k,3) = Evec(i,j,k,3) + 0.25d0 * (alp_r*Bconsxflux(i,j,k) + alp_tmp*Bconsxflux(i+1,j  ,k  ))
              elseif(flux_direction.eq.3) then
                 alp_tmp = 0.5d0 * (alp(i,j+1,k) + alp(i+xoffset,j+yoffset+1,k+zoffset))
                 Evec(i,j,k,1) = Evec(i,j,k,1) + 0.25d0 * (alp_r*Bconsyflux(i,j,k) + alp_tmp*Bconsyflux(i  ,j+1,k  ))
                 alp_tmp = 0.5d0 * (alp(i+1,j,k) + alp(i+xoffset+1,j+yoffset,k+zoffset))
                 Evec(i,j,k,2) = Evec(i,j,k,2) - 0.25d0 * (alp_r*Bconsxflux(i,j,k) + alp_tmp*Bconsxflux(i+1,j  ,k  ))
              end if
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
 
      !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,alp_l,alp_r,alp_tmp,Bcons_l,Bcons_r)
      do k = GRHydro_stencil + 1, cctk_lsh(3) - GRHydro_stencil ! we need to compute Evec on all faces/edges where the fluxes are defined
        do j = GRHydro_stencil + 1, cctk_lsh(2) - GRHydro_stencil
          do i = GRHydro_stencil + 1, cctk_lsh(1) - GRHydro_stencil

            alp_l = 0.5d0 * (alp(i,j,k) + &
                 alp(i-xoffset,j-yoffset,k-zoffset))
            alp_r = 0.5d0 * (alp(i,j,k) + &
                 alp(i+xoffset,j+yoffset,k+zoffset))
            
            densrhs(i,j,k) = densrhs(i,j,k) + &
                 (alp_l * densflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * densflux(i,j,k)) * idx
            srhs(i,j,k,1) = srhs(i,j,k,1) + &
                 (alp_l * sxflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * sxflux(i,j,k)) * idx
            srhs(i,j,k,2) = srhs(i,j,k,2) + &
                 (alp_l * syflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * syflux(i,j,k)) * idx
            srhs(i,j,k,3) = srhs(i,j,k,3) + &
                 (alp_l * szflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * szflux(i,j,k)) * idx
            taurhs(i,j,k) = taurhs(i,j,k) + &
                 (alp_l * tauflux(i-xoffset,j-yoffset,k-zoffset) - &
                 alp_r * tauflux(i,j,k)) * idx
            if(evolve_mhd.ne.0) then
              if(transport_constraints.eq.0) then
                 Bconsrhs(i,j,k,1) = Bconsrhs(i,j,k,1) + &
                      (alp_l * Bconsxflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * Bconsxflux(i,j,k)) * idx
                 Bconsrhs(i,j,k,2) = Bconsrhs(i,j,k,2) + &
                      (alp_l * Bconsyflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * Bconsyflux(i,j,k)) * idx
                 Bconsrhs(i,j,k,3) = Bconsrhs(i,j,k,3) + &
                      (alp_l * Bconszflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * Bconszflux(i,j,k)) * idx
              endif
              if(clean_divergence.ne.0) then
                 psidcrhs(i,j,k) = psidcrhs(i,j,k) + &
                      (alp_l * psidcflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * psidcflux(i,j,k)) * idx
              endif
            endif
    
            if (evolve_tracer .ne. 0) then
               do itracer=1,number_of_tracers
                  cons_tracerrhs(i,j,k,itracer) = cons_tracerrhs(i,j,k,itracer) + &
                       (alp_l * cons_tracerflux(i-xoffset,j-yoffset,k-zoffset,itracer) - &
                       alp_r * cons_tracerflux(i,j,k,itracer)) * idx
               enddo
            end if
            
            if (evolve_Y_e .ne. 0) then
               Y_e_con_rhs(i,j,k) = Y_e_con_rhs(i,j,k) + &
                    (alp_l * Y_e_con_flux(i-xoffset,j-yoffset,k-zoffset) - &
                    alp_r * Y_e_con_flux(i,j,k)) * idx
            end if
            
            if (evolve_entropy .ne. 0) then
               entropyrhs(i,j,k) = entropyrhs(i,j,k) + &
                    (alp_l * entropyflux(i-xoffset,j-yoffset,k-zoffset) - &
                    alp_r * entropyflux(i,j,k)) * idx
            end if

!            densrhs(i,j,k) = 0.0d0
!            taurhs(i,j,k)  = 0.0d0
!            srhs(i,j,k,:) = 0.0d0
!            y_e_con_rhs(i,j,k) = 0.0d0

            
            if (wk_atmosphere .eq. 1) then

              if (atmosphere_mask(i,j,k) .ne. 0) then

!!$                We are in the atmosphere so the momentum flux must vanish

                srhs(i,j,k,:) = 0.d0

                if ( (atmosphere_mask(i-1,j  ,k  ) .ne. 0) .and. &
                     (atmosphere_mask(i+1,j  ,k  ) .ne. 0) .and. &
                     (atmosphere_mask(i  ,j-1,k  ) .ne. 0) .and. &
                     (atmosphere_mask(i  ,j+1,k  ) .ne. 0) .and. &
                     (atmosphere_mask(i  ,j  ,k-1) .ne. 0) .and. &
                     (atmosphere_mask(i  ,j  ,k+1) .ne. 0) &
                   ) then

!!$                    All neighbours are also atmosphere so all rhs vanish

                    densrhs(i,j,k) = 0.d0
                    taurhs(i,j,k)  = 0.d0

                end if
              end if

            end if

          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    else
      
      call CCTK_ERROR("Not supported")
      STOP
      
!!$    do k = GRHydro_stencil + 1, cctk_lsh(3) - GRHydro_stencil
!!$      do j = GRHydro_stencil + 1, cctk_lsh(2) - GRHydro_stencil
!!$        do i = GRHydro_stencil + 1, cctk_lsh(1) - GRHydro_stencil
!!$          
!!$          alp_l = 0.5d0 * (alp(i,j,k) + &
!!$               alp(i-xoffset,j-yoffset,k-zoffset))
!!$          alp_r = 0.5d0 * (alp(i,j,k) + &
!!$               alp(i+xoffset,j+yoffset,k+zoffset))
!!$          
!!$          densrhs(i,j,k) = densrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &densflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &densflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          sxrhs(i,j,k) = sxrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &sxflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &sxflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          syrhs(i,j,k) = syrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &syflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &syflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          szrhs(i,j,k) = szrhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &szflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &szflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          taurhs(i,j,k) = taurhs(i,j,k) + &
!!$               (alp_l * &
!!$               &cell_surface(i-xoffset,j-yoffset,k-zoffset,flux_direction) * &
!!$               &tauflux(i-xoffset,j-yoffset,k-zoffset) - &
!!$               alp_r * &
!!$               &cell_surface(i,j,k,flux_direction) * &
!!$               &tauflux(i,j,k)) * idx / cell_volume(i,j,k)
!!$          
!!$        enddo
!!$      enddo
!!$    enddo

    end if
  
  else if (CCTK_EQUALS(method_type, "Flux split FD")) then

    if (evolve_mhd.ne.0 .and. transport_constraints .ne. 0) then
      call CCTK_ERROR("Not supported")
      STOP
    end if

    do k = GRHydro_stencil + 1, cctk_lsh(3) - GRHydro_stencil
      do j = GRHydro_stencil + 1, cctk_lsh(2) - GRHydro_stencil
        do i = GRHydro_stencil + 1, cctk_lsh(1) - GRHydro_stencil
          
          densrhs(i,j,k) = densrhs(i,j,k) + &
               (densflux(i-xoffset,j-yoffset,k-zoffset) - &
                densflux(i,j,k)) * idx
          srhs(i,j,k,1) = srhs(i,j,k,1) + &
               (sxflux(i-xoffset,j-yoffset,k-zoffset) - &
                sxflux(i,j,k)) * idx
          srhs(i,j,k,2) = srhs(i,j,k,2) + &
               (syflux(i-xoffset,j-yoffset,k-zoffset) - &
                syflux(i,j,k)) * idx
          srhs(i,j,k,3) = srhs(i,j,k,3) + &
               (szflux(i-xoffset,j-yoffset,k-zoffset) - &
                szflux(i,j,k)) * idx
          taurhs(i,j,k) = taurhs(i,j,k) + &
               (tauflux(i-xoffset,j-yoffset,k-zoffset) - &
               tauflux(i,j,k)) * idx
          
          if(evolve_mhd.ne.0) then
             Bconsrhs(i,j,k,1) = Bconsrhs(i,j,k,1) + &
                  (Bconsxflux(i-xoffset,j-yoffset,k-zoffset) - &
                  Bconsxflux(i,j,k)) * idx
             Bconsrhs(i,j,k,2) = Bconsrhs(i,j,k,2) + &
                  (Bconsyflux(i-xoffset,j-yoffset,k-zoffset) - &
                  Bconsyflux(i,j,k)) * idx
             Bconsrhs(i,j,k,3) = Bconsrhs(i,j,k,3) + &
                  (Bconszflux(i-xoffset,j-yoffset,k-zoffset) - &
                  Bconszflux(i,j,k)) * idx
             if(clean_divergence.ne.0) then
                psidcrhs(i,j,k) = psidcrhs(i,j,k) + &
                     (psidcflux(i-xoffset,j-yoffset,k-zoffset) - &
                     psidcflux(i,j,k)) * idx
             endif
             if(track_divB.ne.0) then
               Bcons_l = 0.5d0 * (Bcons(i,j,k,flux_direction) + &
                        Bcons(i-xoffset,j-yoffset,k-zoffset,flux_direction))
               Bcons_r = 0.5d0 * (Bcons(i,j,k,flux_direction) + &
                        Bcons(i+xoffset,j+yoffset,k+zoffset,flux_direction))
               divB(i,j,k) = divB(i,j,k) + ( Bcons_l - Bcons_r ) * idx 
             endif
          endif

        enddo
      enddo
    enddo
      
  end if

  if (evolve_mhd.ne.0 .and. transport_constraints.ne.0 .and. &
      flux_direction.eq.1) then ! HACK: x direction is last
    ! FIXME: I think one could wrap all of this into a single do loop and remove the
    !        Evec storage
    ! idx differs from idx which was 1d0/CCTK_DELTA_SPACE(flux_direction)
    idx = 1d0 / CCTK_DELTA_SPACE(1)
    idy = 1d0 / CCTK_DELTA_SPACE(2)
    idz = 1d0 / CCTK_DELTA_SPACE(3)
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = GRHydro_stencil + 1, cctk_lsh(3) - GRHydro_stencil
      do j = GRHydro_stencil + 1, cctk_lsh(2) - GRHydro_stencil
        do i = GRHydro_stencil + 1, cctk_lsh(1) - GRHydro_stencil
          Bconsrhs(i,j,k,1) = - 0.5d0 * ((Evec(i-1,j  ,k-1,2)-Evec(i-1,j  ,k  ,2)) &
                                       + (Evec(i  ,j  ,k-1,2)-Evec(i  ,j  ,k  ,2))) * idz &
                              - 0.5d0 * ((Evec(i-1,j  ,k  ,3)-Evec(i-1,j-1,k  ,3)) &
                                       + (Evec(i  ,j  ,k  ,3)-Evec(i  ,j-1,k  ,3))) * idy
          Bconsrhs(i,j,k,2) = - 0.5d0 * ((Evec(i-1,j-1,k  ,3)-Evec(i  ,j-1,k  ,3)) &
                                       + (Evec(i-1,j  ,k  ,3)-Evec(i  ,j  ,k  ,3))) * idx &
                              - 0.5d0 * ((Evec(i  ,j-1,k  ,1)-Evec(i  ,j-1,k-1,1)) &
                                       + (Evec(i  ,j  ,k  ,1)-Evec(i  ,j  ,k-1,1))) * idz
          Bconsrhs(i,j,k,3) = - 0.5d0 * ((Evec(i  ,j-1,k-1,1)-Evec(i  ,j  ,k-1,1)) &
                                       + (Evec(i  ,j-1,k  ,1)-Evec(i  ,j  ,k  ,1))) * idy &
                              - 0.5d0 * ((Evec(i  ,j  ,k-1,2)-Evec(i-1,j  ,k-1,2)) &
                                       + (Evec(i  ,j  ,k  ,2)-Evec(i-1,j  ,k  ,2))) * idx
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end if
  
  return

end subroutine UpdateCalculation
  




 /*@@
   @routine    ConstrainSconTo1D
   @date       Tue  24 14:12 2012
   @author     Christian Reisswig
   @desc 
   Constrains the conserved fluid velocity to radial direction
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 

@@*/


subroutine ConstrainSconTo1D(CCTK_ARGUMENTS)
   implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_REAL :: rnorm, rnormI, scon_tmp1, scon_tmp2, scon_tmp3

      !$OMP PARALLEL DO PRIVATE(i,j,k,rnorm,rnormI,scon_tmp1,scon_tmp2,scon_tmp3)
      do k = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(3) - GRHydro_stencil ! we need to compute Evec on all faces/edges where the fluxes are defined
        do j = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(2) - GRHydro_stencil
          do i = GRHydro_stencil + 1 - transport_constraints, cctk_lsh(1) - GRHydro_stencil

            
            ! Eliminate non-radial fluid velocities to obtain pseudo 1D scheme
               rnorm = (x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2)
               if (rnorm.lt.1.0d-10) then
                  rnormI = 0.0d0
               else
                  rnormI = 1.0d0/rnorm
               endif
               
               scon_tmp1 = (x(i,j,k)*scon(i,j,k,1) &
                          + y(i,j,k)*scon(i,j,k,2) &
                          + z(i,j,k)*scon(i,j,k,3)) * rnormI * x(i,j,k)
               scon_tmp2 = (x(i,j,k)*scon(i,j,k,1) &
                          + y(i,j,k)*scon(i,j,k,2) &
                          + z(i,j,k)*scon(i,j,k,3)) * rnormI * y(i,j,k)
               scon_tmp3 = (x(i,j,k)*scon(i,j,k,1) &
                          + y(i,j,k)*scon(i,j,k,2) &
                          + z(i,j,k)*scon(i,j,k,3)) * rnormI * z(i,j,k)
               
               scon(i,j,k,1) = scon_tmp1
               scon(i,j,k,2) = scon_tmp2
               scon(i,j,k,3) = scon_tmp3
               
        end do
      end do
    end do
end subroutine ConstrainSconTo1D

