 /*@@
   @file      GRHydro_HLLE.F90
   @date      Sat Jan 26 01:40:14 2002
   @author    Ian Hawke, Pedro Montero, Toni Font
   @desc 
   The HLLE solver. Called from the wrapper function, so works in 
   all directions.
   @enddesc 
 @@*/
   
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_Macros.h"
#include "SpaceMask.h"

 /*@@
   @routine    GRHydro_HLLE
   @date       Sat Jan 26 01:41:02 2002
   @author     Ian Hawke, Pedro Montero, Toni Font
   @desc 
   The HLLE solver. Sufficiently simple that its just one big routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Altered from Cactus 3 routines originally written by Toni Font.
   @endhistory 

@@*/

! eta across face in x-direction
#define etaX(i,j,k, v, c) \
  (0.5d0*(abs(v((i)+1,j,k,1)-v(i,j,k,1)) + abs(c((i)+1,j,k)-c(i,j,k))))

! eta across face in x-direction
#define etaY(i,j,k, v, c) \
  (0.5d0*(abs(v(i,(j)+1,k,2)-v(i,j,k,2)) + abs(c(i,(j)+1,k)-c(i,j,k))))
  
! eta across face in x-direction
#define etaZ(i,j,k, v, c) \
  (0.5d0*(abs(v(i,j,(k)+1,3)-v(i,j,k,3)) + abs(c(i,j,(k)+1)-c(i,j,k))))

  

subroutine H_viscosity(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, nx, ny, nz
  CCTK_REAL dpdrho,dpdeps, cs2
  character*256 :: warnline
  
  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)
  
! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1)
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars

  ! save memory when MP is not used
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

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)   

  !$OMP PARALLEL DO PRIVATE(i,j,k,&
  !$OMP anyerr, keyerr, keytemp,&
  !$OMP warnline, xpress,xeps,xtemp,xye, dpdrho, dpdeps, cs2)
  do k = 1, nz 
    do j = 1, ny 
      do i = 1, nx

        ! set to zero initially
        eos_c(i,j,k) = 0.0d0
      
        !do not compute if in atmosphere or in excised region
        if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
            (GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .gt. 0))) cycle

!!$  Set required fluid quantities
        if (evolve_temper.ne.1) then

           call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
              rho(i,j,k),eps(i,j,k),xtemp,xye,xpress,keyerr,anyerr)

            call EOS_Omni_DPressByDEps(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
              rho(i,j,k),eps(i,j,k),xtemp,xye,dpdeps,keyerr,anyerr)

           call EOS_Omni_DPressByDRho(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
              rho(i,j,k),eps(i,j,k),xtemp,xye,dpdrho,keyerr,anyerr)

           cs2 = (dpdrho + xpress(1) * dpdeps / (rho(i,j,k)**2))/ &
              (1.0d0 + eps(i,j,k) + xpress(1)/rho(i,j,k))

           eos_c(i,j,k) = sqrt(cs2)
         
         else
           
           call EOS_Omni_cs2(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
              rho(i,j,k),eps(i,j,k),temperature(i,j,k),Y_e(i,j,k),cs2,keyerr,anyerr)
  
           eos_c(i,j,k) = sqrt(cs2)
         
         end if
  
       end do
    end do
  end do
  
end subroutine H_viscosity
  
  
subroutine GRHydro_HLLE(CCTK_ARGUMENTS)
  USE GRHydro_Eigenproblem

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k, m
  CCTK_INT :: keytemp
  CCTK_REAL, dimension(5) :: consp,consm_i,fplus,fminus,lamplus
  CCTK_REAL, dimension(5) :: f1,lamminus
  CCTK_REAL, dimension(5) :: qdiff
  CCTK_REAL ::  charmin, charmax, charpm,avg_alp,avg_det, etabar
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_beta, usendh
    
  CCTK_INT :: type_bits, trivial

  ! sign requires its arguments to be of identical KIND
  CCTK_REAL, parameter :: one = 1d0
  
  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: beta1, beta2, beta3
  pointer (pbeta1,beta1), (pbeta2,beta2), (pbeta3,beta3)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)
  
  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pbeta1 = loc(betaa)
    pbeta2 = loc(betab)
    pbeta3 = loc(betac)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pbeta1 = loc(betax)
    pbeta2 = loc(betay)
    pbeta3 = loc(betaz)
    pvup = loc(vel)
  end if

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

  if(evolve_temper.eq.1.and.reconstruct_temper.eq.1) then
     keytemp = 1
  else
     keytemp = 0
  endif

  !$OMP PARALLEL DO PRIVATE(k,j,i,f1,lamminus,lamplus,consp,consm_i, &
  !$OMP                     fplus,fminus,qdiff,avg_beta,avg_alp,   &
  !$OMP                     avg_det,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
  !$OMP                     uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,         &
  !$OMP                     usendh, charmin, charmax, charpm, &
  !$OMP                     m,etabar)
  do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil
        
        f1 = 0.d0
        lamminus = 0.d0
        lamplus = 0.d0
        consp = 0.d0
        consm_i = 0.d0
        fplus = 0.d0
        fminus = 0.d0
        qdiff = 0.d0
        
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
          
!!$        Calculate various metric terms here.
!!$        Note also need the average of the lapse at the 
!!$        left and right points.

        if (flux_direction == 1) then
           avg_beta = 0.5d0 * (beta1(i+xoffset,j+yoffset,k+zoffset) + &
               beta1(i,j,k))
        else if (flux_direction == 2) then
           avg_beta = 0.5d0 * (beta2(i+xoffset,j+yoffset,k+zoffset) + &
               beta2(i,j,k))
        else if (flux_direction == 3) then
           avg_beta = 0.5d0 * (beta3(i+xoffset,j+yoffset,k+zoffset) + &
               beta3(i,j,k))
        else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
         end if

        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gxxh = 0.5d0 * (g11(i+xoffset,j+yoffset,k+zoffset) + g11(i,j,k))
        gxyh = 0.5d0 * (g12(i+xoffset,j+yoffset,k+zoffset) + g12(i,j,k))
        gxzh = 0.5d0 * (g13(i+xoffset,j+yoffset,k+zoffset) + g13(i,j,k))
        gyyh = 0.5d0 * (g22(i+xoffset,j+yoffset,k+zoffset) + g22(i,j,k))
        gyzh = 0.5d0 * (g23(i+xoffset,j+yoffset,k+zoffset) + g23(i,j,k))
        gzzh = 0.5d0 * (g33(i+xoffset,j+yoffset,k+zoffset) + g33(i,j,k))

        avg_det =  SPATIAL_DETERMINANT(gxxh,gxyh,gxzh,\
             gyyh,gyzh,gzzh)

!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell
          
        if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then

          if (flux_direction == 1) then
            call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                 f1(1),f1(2),f1(3),&
                 f1(4),f1(5),&
                 velxplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 2) then
            call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                 f1(1),f1(3),f1(4),&
                 f1(2),f1(5),&
                 velyplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 3) then
            call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                 f1(1),f1(4),f1(2),&
                 f1(3),f1(5),&
                 velzplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if
          
        else !!! The end of this branch is right at the bottom of the routine
            
          call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
               avg_det,gxxh, gxyh, gxzh, &
               gyyh, gyzh, gzzh)
          
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
          
!!$        Calculate the jumps in the conserved variables
          
          qdiff(1) = consm_i(1) - consp(1)
          qdiff(2) = consm_i(2) - consp(2)
          qdiff(3) = consm_i(3) - consp(3)
          qdiff(4) = consm_i(4) - consp(4)
          qdiff(5) = consm_i(5) - consp(5)
          
!!$        Eigenvalues and fluxes either side of the cell interface
          
          if (flux_direction == 1) then
             if(evolve_temper.ne.1) then
                call eigenvalues(GRHydro_eos_handle,&
                     rhominus(i+xoffset,j+yoffset,k+zoffset),&
                     velxminus(i+xoffset,j+yoffset,k+zoffset),&
                     velyminus(i+xoffset,j+yoffset,k+zoffset),&
                     velzminus(i+xoffset,j+yoffset,k+zoffset),&
                     epsminus(i+xoffset,j+yoffset,k+zoffset),&
                     w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                     lamminus,gxxh,gxyh,gxzh,gyyh,&
                     gyzh,gzzh,&
                     usendh,avg_alp,avg_beta)
                call eigenvalues(GRHydro_eos_handle,rhoplus(i,j,k),&
                     velxplus(i,j,k),velyplus(i,j,k),&
                     velzplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
                     lamplus,gxxh,gxyh,gxzh,gyyh,&
                     gyzh,gzzh,&
                     usendh,avg_alp,avg_beta)
             else
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,&
                     int(i,ik),int(j,ik),int(k,ik),&
                     rhominus(i+xoffset,j+yoffset,k+zoffset),&
                     velxminus(i+xoffset,j+yoffset,k+zoffset),&
                     velyminus(i+xoffset,j+yoffset,k+zoffset),&
                     velzminus(i+xoffset,j+yoffset,k+zoffset),&
                     epsminus(i+xoffset,j+yoffset,k+zoffset),&
                     tempminus(i+xoffset,j+yoffset,k+zoffset),&
                     y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                     w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                     lamminus,gxxh,gxyh,gxzh,gyyh,&
                     gyzh,gzzh,&
                     usendh,avg_alp,avg_beta)
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,&
                     int(i,ik),int(j,ik),int(k,ik),&
                     rhoplus(i,j,k),&
                     velxplus(i,j,k),velyplus(i,j,k),&
                     velzplus(i,j,k),epsplus(i,j,k), &
                     tempplus(i,j,k),&
                     y_e_plus(i,j,k),&
                     w_lorentzplus(i,j,k),&
                     lamplus,gxxh,gxyh,gxzh,gyyh,&
                     gyzh,gzzh,&
                     usendh,avg_alp,avg_beta)
             endif
            call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                 fplus(1),fplus(2),fplus(3),fplus(4),&
                 fplus(5),velxplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(2),consm_i(3),&
                 consm_i(4),consm_i(5),fminus(1),fminus(2),fminus(3),&
                 fminus(4), fminus(5),&
                 velxminus(i+xoffset,j+yoffset,k+zoffset),&
                 pressminus(i+xoffset,j+yoffset,k+zoffset),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 2) then
             if(evolve_temper.ne.1) then
                call eigenvalues(GRHydro_eos_handle,&
                     rhominus(i+xoffset,j+yoffset,k+zoffset),&
                     velyminus(i+xoffset,j+yoffset,k+zoffset),&
                     velzminus(i+xoffset,j+yoffset,k+zoffset),&
                     velxminus(i+xoffset,j+yoffset,k+zoffset),&
                     epsminus(i+xoffset,j+yoffset,k+zoffset),&
                     w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                     lamminus,gyyh,gyzh,gxyh,gzzh,&
                     gxzh,gxxh,&
                     usendh,avg_alp,avg_beta)
                call eigenvalues(GRHydro_eos_handle,rhoplus(i,j,k),&
                     velyplus(i,j,k),velzplus(i,j,k),&
                     velxplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
                     lamplus,gyyh,gyzh,gxyh,gzzh,&
                     gxzh,gxxh,&
                     usendh,avg_alp,avg_beta)
             else
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,int(i,ik),int(j,ik),int(k,ik),&
                     rhominus(i+xoffset,j+yoffset,k+zoffset),&
                     velyminus(i+xoffset,j+yoffset,k+zoffset),&
                     velzminus(i+xoffset,j+yoffset,k+zoffset),&
                     velxminus(i+xoffset,j+yoffset,k+zoffset),&
                     epsminus(i+xoffset,j+yoffset,k+zoffset),&
                     tempminus(i+xoffset,j+yoffset,k+zoffset),&
                     y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                     w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                     lamminus,gyyh,gyzh,gxyh,gzzh,&
                     gxzh,gxxh,&
                     usendh,avg_alp,avg_beta)
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,int(i,ik),int(j,ik),int(k,ik),&
                     rhoplus(i,j,k),&
                     velyplus(i,j,k),velzplus(i,j,k),&
                     velxplus(i,j,k),epsplus(i,j,k),&
                     tempplus(i,j,k),y_e_plus(i,j,k),&
                     w_lorentzplus(i,j,k),&
                     lamplus,gyyh,gyzh,gxyh,gzzh,&
                     gxzh,gxxh,&
                     usendh,avg_alp,avg_beta)
             endif
            call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                 fplus(1),fplus(3),fplus(4),fplus(2),&
                 fplus(5),velyplus(i,j,k),pressplus(i,j,k),&
                 avg_det,avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(3),consm_i(4),&
                 consm_i(2),consm_i(5),fminus(1),fminus(3),fminus(4),&
                 fminus(2), fminus(5),&
                 velyminus(i+xoffset,j+yoffset,k+zoffset),&
                 pressminus(i+xoffset,j+yoffset,k+zoffset),&
                 avg_det,avg_alp,avg_beta)
          else if (flux_direction == 3) then
             if(evolve_temper.ne.1) then
                call eigenvalues(GRHydro_eos_handle,&
                     rhominus(i+xoffset,j+yoffset,k+zoffset),&
                     velzminus(i+xoffset,j+yoffset,k+zoffset),&
                     velxminus(i+xoffset,j+yoffset,k+zoffset),&
                     velyminus(i+xoffset,j+yoffset,k+zoffset),&
                     epsminus(i+xoffset,j+yoffset,k+zoffset),&
                     w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                     lamminus,gzzh,gxzh,gyzh,&
                     gxxh,gxyh,gyyh,&
                     usendh,avg_alp,avg_beta)
                call eigenvalues(GRHydro_eos_handle,rhoplus(i,j,k),&
                     velzplus(i,j,k),velxplus(i,j,k),&
                     velyplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
                     lamplus,gzzh,gxzh,gyzh,&
                     gxxh,gxyh,gyyh,&
                     usendh,avg_alp,avg_beta)
             else
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,int(i,ik),int(j,ik),int(k,ik),&
                     rhominus(i+xoffset,j+yoffset,k+zoffset),&
                     velzminus(i+xoffset,j+yoffset,k+zoffset),&
                     velxminus(i+xoffset,j+yoffset,k+zoffset),&
                     velyminus(i+xoffset,j+yoffset,k+zoffset),&
                     epsminus(i+xoffset,j+yoffset,k+zoffset),&
                     tempminus(i+xoffset,j+yoffset,k+zoffset),&
                     y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                     w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
                     lamminus,gzzh,gxzh,gyzh,&
                     gxxh,gxyh,gyyh,&
                     usendh,avg_alp,avg_beta)
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,&
                     int(i,ik),int(j,ik),int(k,ik),&
                     rhoplus(i,j,k),&
                     velzplus(i,j,k),velxplus(i,j,k),&
                     velyplus(i,j,k),epsplus(i,j,k),&
                     tempplus(i,j,k),y_e_plus(i,j,k),&
                     w_lorentzplus(i,j,k),&
                     lamplus,gzzh,gxzh,gyzh,&
                     gxxh,gxyh,gyyh,&
                     usendh,avg_alp,avg_beta)
             endif
            call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                 fplus(1),fplus(4),fplus(2),fplus(3),&
                 fplus(5),velzplus(i,j,k),pressplus(i,j,k),avg_det,&
                 avg_alp,avg_beta)
            call num_x_flux(consm_i(1),consm_i(4),consm_i(2),&
                 consm_i(3),consm_i(5),fminus(1),fminus(4),fminus(2),&
                 fminus(3), fminus(5),&
                 velzminus(i+xoffset,j+yoffset,k+zoffset),&
                 pressminus(i+xoffset,j+yoffset,k+zoffset),&
                 avg_det,avg_alp,avg_beta)
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if

!!$        Compute H viscosity if requested
          
          if (apply_H_viscosity .ne. 0) then
          
            if (flux_direction == 1) then
              etabar = max(etaX(i,j,k, vup, eos_c),&
                           etaY(i,j,k, vup, eos_c),&
                           etaY(i+1,j,k, vup, eos_c),&
                           etaY(i,j-1,k, vup, eos_c),&
                           etaY(i+1,j-1,k, vup, eos_c),&
                           etaZ(i,j,k, vup, eos_c),&
                           etaZ(i+1,j,k, vup, eos_c),&
                           etaZ(i,j,k-1, vup, eos_c),&
                           etaZ(i+1,j,k-1, vup, eos_c))
            else if (flux_direction == 2) then
              etabar = max(etaY(i,j,k, vup, eos_c),&
                           etaX(i,j,k, vup, eos_c),&
                           etaX(i,j+1,k, vup, eos_c),&
                           etaX(i-1,j,k, vup, eos_c),&
                           etaX(i-1,j+1,k, vup, eos_c),&
                           etaZ(i,j,k, vup, eos_c),&
                           etaZ(i,j+1,k, vup, eos_c),&
                           etaZ(i,j,k-1, vup, eos_c),&
                           etaZ(i,j+1,k-1, vup, eos_c))
            else if (flux_direction == 3) then
              etabar = max(etaZ(i,j,k, vup, eos_c),&
                           etaX(i,j,k, vup, eos_c),&
                           etaX(i,j,k+1, vup, eos_c),&
                           etaX(i-1,j,k, vup, eos_c),&
                           etaX(i-1,j,k+1, vup, eos_c),&
                           etaY(i,j,k, vup, eos_c),&
                           etaY(i,j,k+1, vup, eos_c),&
                           etaY(i,j-1,k, vup, eos_c),&
                           etaY(i,j-1,k+1, vup, eos_c))
            else
              call CCTK_ERROR("Flux direction not x,y,z")
              STOP
            end if
            
            ! modify eigenvalues of Roe's matrix by computed H viscosity
            do m = 1,5
               lamplus(m) = sign(one, lamplus(m))*max(abs(lamplus(m)), etabar)
               lamminus(m) = sign(one, lamminus(m))*max(abs(lamminus(m)), etabar)
            end do
          endif
          
          
!!$        Find minimum and maximum wavespeeds
      
          charmin = min(0.d0, lamplus(1), lamplus(2), lamplus(3), &
               lamplus(4),lamplus(5),  lamminus(1),lamminus(2),lamminus(3),&
               lamminus(4),lamminus(5))  
          
          charmax = max(0.d0, lamplus(1), lamplus(2), lamplus(3), &
               lamplus(4),lamplus(5),  lamminus(1),lamminus(2),lamminus(3),&
               lamminus(4),lamminus(5))
          
          charpm = charmax - charmin
          
!!$        Calculate flux by standard formula
          
          do m = 1,5 
            
            qdiff(m) = consm_i(m) - consp(m)
            
            f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
                 charmax * charmin * qdiff(m)) / charpm
            
          end do

        end if !!! The end of the SpaceMask check for a trivial RP.
            
        densflux(i, j, k) = f1(1)
        sxflux(i, j, k) = f1(2)
        syflux(i, j, k) = f1(3)
        szflux(i, j, k) = f1(4)
        tauflux(i, j, k) = f1(5)

        if(evolve_Y_e.ne.0) then
           if (densflux(i, j, k) > 0.d0) then
              Y_e_con_flux(i, j, k) = &
                   Y_e_plus(i, j, k) * &
                   densflux(i, j, k)
           else
              Y_e_con_flux(i, j, k) = &
                   Y_e_minus(i + xoffset, j + yoffset, k + zoffset) * &
                   densflux(i, j, k)
           endif
        endif

      end do
    end do
  end do
  !$OMP END PARALLEL DO
      
end subroutine GRHydro_HLLE

 /*@@
   @routine    GRHydro_HLLE_Tracer
   @date       Mon Mar  8 13:47:13 2004
   @author     Ian Hawke
   @desc 
   HLLE just for the tracer.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_HLLE_Tracer(CCTK_ARGUMENTS)

  USE GRHydro_Eigenproblem

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k,  m
  CCTK_REAL, dimension(number_of_tracers) :: consp,consm_i,fplus,fminus,f1
  CCTK_REAL, dimension(5) :: lamminus,lamplus
  CCTK_REAL, dimension(number_of_tracers) :: qdiff
  CCTK_REAL ::  charmin, charmax, charpm,avg_alp,avg_det
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_beta, usendh
    
  CCTK_INT :: type_bits, trivial

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: beta1, beta2, beta3
  pointer (pbeta1,beta1), (pbeta2,beta2), (pbeta3,beta3)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pbeta1 = loc(betaa)
    pbeta2 = loc(betab)
    pbeta3 = loc(betac)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pbeta1 = loc(betax)
    pbeta2 = loc(betay)
    pbeta3 = loc(betaz)
  end if

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

  do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil
        
        f1 = 0.d0
        lamminus = 0.d0
        lamplus = 0.d0
        consp = 0.d0
        consm_i = 0.d0
        fplus = 0.d0
        fminus = 0.d0
        qdiff = 0.d0
        
!!$        Set the left (p for plus) and right (m_i for minus, i+1) states
        
        do m=1,number_of_tracers
           consp(m)   = cons_tracerplus(i,j,k,m)         
           consm_i(m) = cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,m)
        enddo
!!$        Calculate various metric terms here.
!!$        Note also need the average of the lapse at the 
!!$        left and right points.

        if (flux_direction == 1) then
           avg_beta = 0.5d0 * (beta1(i+xoffset,j+yoffset,k+zoffset) + &
               beta1(i,j,k))
        else if (flux_direction == 2) then
           avg_beta = 0.5d0 * (beta2(i+xoffset,j+yoffset,k+zoffset) + &
               beta2(i,j,k))
        else if (flux_direction == 3) then
           avg_beta = 0.5d0 * (beta3(i+xoffset,j+yoffset,k+zoffset) + &
               beta3(i,j,k))
        else
           call CCTK_ERROR("Flux direction not x,y,z")
           STOP
        end if

        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gxxh = 0.5d0 * (g11(i+xoffset,j+yoffset,k+zoffset) + g11(i,j,k))
        gxyh = 0.5d0 * (g12(i+xoffset,j+yoffset,k+zoffset) + g12(i,j,k))
        gxzh = 0.5d0 * (g13(i+xoffset,j+yoffset,k+zoffset) + g13(i,j,k))
        gyyh = 0.5d0 * (g22(i+xoffset,j+yoffset,k+zoffset) + g22(i,j,k))
        gyzh = 0.5d0 * (g23(i+xoffset,j+yoffset,k+zoffset) + g23(i,j,k))
        gzzh = 0.5d0 * (g33(i+xoffset,j+yoffset,k+zoffset) + g33(i,j,k))

        avg_det =  SPATIAL_DETERMINANT(gxxh,gxyh,gxzh,\
             gyyh,gyzh,gzzh)

        call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
             avg_det,gxxh, gxyh, gxzh, &
             gyyh, gyzh, gzzh)
        
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
          
!!$        Calculate the jumps in the conserved variables
          
        qdiff = consm_i - consp
          
!!$        Eigenvalues and fluxes either side of the cell interface
          
        if (flux_direction == 1) then
          call eigenvalues(GRHydro_eos_handle,&
               rhominus(i+xoffset,j+yoffset,k+zoffset),&
               velxminus(i+xoffset,j+yoffset,k+zoffset),&
               velyminus(i+xoffset,j+yoffset,k+zoffset),&
               velzminus(i+xoffset,j+yoffset,k+zoffset),&
               epsminus(i+xoffset,j+yoffset,k+zoffset),&
               w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
               lamminus,gxxh,gxyh,gxzh,gyyh,&
               gyzh,gzzh,&
               usendh,avg_alp,avg_beta)
          call eigenvalues(GRHydro_eos_handle,rhoplus(i,j,k),&
               velxplus(i,j,k),velyplus(i,j,k),&
               velzplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
               lamplus,gxxh,gxyh,gxzh,gyyh,&
               gyzh,gzzh,&
               usendh,avg_alp,avg_beta)
          fplus(:)  = (velxplus(i,j,k)  - avg_beta / avg_alp) * &
               cons_tracerplus(i,j,k,:)
          fminus(:) = (velxminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
               cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else if (flux_direction == 2) then
          call eigenvalues(GRHydro_eos_handle,&
               rhominus(i+xoffset,j+yoffset,k+zoffset),&
               velyminus(i+xoffset,j+yoffset,k+zoffset),&
               velzminus(i+xoffset,j+yoffset,k+zoffset),&
               velxminus(i+xoffset,j+yoffset,k+zoffset),&
               epsminus(i+xoffset,j+yoffset,k+zoffset),&
               w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
               lamminus,gyyh,gyzh,gxyh,gzzh,&
               gxzh,gxxh,&
               usendh,avg_alp,avg_beta)
          call eigenvalues(GRHydro_eos_handle,rhoplus(i,j,k),&
               velyplus(i,j,k),velzplus(i,j,k),&
               velxplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
               lamplus,gyyh,gyzh,gxyh,gzzh,&
               gxzh,gxxh,&
               usendh,avg_alp,avg_beta)
          fplus(:)  = (velyplus(i,j,k)  - avg_beta / avg_alp) * &
               cons_tracerplus(i,j,k,:)
          fminus(:) = (velyminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
               cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else if (flux_direction == 3) then
          call eigenvalues(GRHydro_eos_handle,&
               rhominus(i+xoffset,j+yoffset,k+zoffset),&
               velzminus(i+xoffset,j+yoffset,k+zoffset),&
               velxminus(i+xoffset,j+yoffset,k+zoffset),&
               velyminus(i+xoffset,j+yoffset,k+zoffset),&
               epsminus(i+xoffset,j+yoffset,k+zoffset),&
               w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),&
               lamminus,gzzh,gxzh,gyzh,&
               gxxh,gxyh,gyyh,&
               usendh,avg_alp,avg_beta)
          call eigenvalues(GRHydro_eos_handle,rhoplus(i,j,k),&
               velzplus(i,j,k),velxplus(i,j,k),&
               velyplus(i,j,k),epsplus(i,j,k),w_lorentzplus(i,j,k),&
               lamplus,gzzh,gxzh,gyzh,&
               gxxh,gxyh,gyyh,&
               usendh,avg_alp,avg_beta)
          fplus(:)  = (velzplus(i,j,k) - avg_beta / avg_alp) * &
               cons_tracerplus(i,j,k,:)
          fminus(:) = (velzminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
               cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else
          call CCTK_ERROR("Flux direction not x,y,z")
          STOP
        end if
        
!!$        Find minimum and maximum wavespeeds
      
        charmin = min(0.d0, lamplus(1), lamplus(2), lamplus(3), &
             lamplus(4),lamplus(5),  lamminus(1),lamminus(2),lamminus(3),&
             lamminus(4),lamminus(5))  
          
        charmax = max(0.d0, lamplus(1), lamplus(2), lamplus(3), &
             lamplus(4),lamplus(5),  lamminus(1),lamminus(2),lamminus(3),&
             lamminus(4),lamminus(5))
        
        charpm = charmax - charmin
          
!!$        Calculate flux by standard formula
          
        do m = 1,number_of_tracers
            
          qdiff(m) = consm_i(m) - consp(m)
          
          f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
               charmax * charmin * qdiff(m)) / charpm
          
        end do
            
        cons_tracerflux(i, j, k,:) = f1(:)
!!$
!!$        if ( ((flux_direction.eq.3).and.(i.eq.4).and.(j.eq.4)).or.&
!!$             ((flux_direction.eq.2).and.(i.eq.4).and.(k.eq.4)).or.&
!!$             ((flux_direction.eq.1).and.(j.eq.4).and.(k.eq.4))&
!!$             ) then
!!$          write(*,*) flux_direction, i, j, k, f1(1), consm_i(1), consp(1)
!!$        end if
        
      end do
    end do
  end do

      
end subroutine GRHydro_HLLE_Tracer

