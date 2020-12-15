 /*@@
   @file      GRHydro_HLLC.F90
   @date      Sat Jan 26 01:40:14 2002
   @author    Ian Hawke, Pedro Montero, Toni Font, Roland Haas
   @desc 
   The HLLC solver. Called from the wrapper function, so works in 
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
   @routine    GRHydro_HLLC
   @date       Sat Jan 26 01:41:02 2002
   @author     Ian Hawke, Pedro Montero, Toni Font, Roland Haas
   @desc 
   The HLLC solver. Sufficiently simple that its just one big routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Based entirely on Matt Duez notes on HLL and HLLC
   @endhistory 

@@*/

module GRHydro_hllc_helpers

  implicit none

  contains
  ! using nested subroutines seems to not work (either because gfortran has a
  ! bug or because they are not allowed with OpenMP). So we have to pass each
  ! argument explictly. ok for now since I plan to re-implement this in C
  ! anyway.
  subroutine compute_stardata(charstar,pressstar,consm_i,consp,sqrt_avg_det,pressminus,pressplus,&
                              uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,usendh,avg_beta,avg_alp,vxR,vxL,charmin,charmax,i,j,k,&
                              cctk_ash1,cctk_ash2,cctk_ash3,xoffset,yoffset,zoffset,flux_direction)

    implicit none

    CCTK_REAL, intent(out) :: charstar, pressstar

    integer, intent(in) :: cctk_ash1, cctk_ash2, cctk_ash3
    CCTK_REAL, intent(in), DIMENSION(5) :: consm_i, consp
    CCTK_REAL, intent(in) :: sqrt_avg_det, avg_beta, avg_alp
    CCTK_REAL, intent(in), DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: pressminus, pressplus
    CCTK_REAL, intent(in) :: uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,usendh
    CCTK_REAL, intent(in) :: charmin, charmax
    CCTK_REAL, intent(in) :: vxR, vxL
    integer, intent(in) :: i,j,k
    CCTK_INT, intent(in) :: xoffset, yoffset, zoffset, flux_direction

    ! naming convention: in this routine we use the convention of Duez (and by
    ! extenstion Mignone et al.MNRAS 364 126 (2005) which means that L<->plus,
    ! R<->minus
    ! for convenience, use s^k = D h u^k = D h g^{k \mu} u_\mu
    ! ie. this is NOT s^k = gamma^{kj} s_j
    ! doing the math, I find: 
    ! s^k = gamma^{kj} s_j - beta^k/alpha sqrt(detgamma) rho h W**2
    ! where the last part is rho h W**2 = tau+D+press
    CCTK_REAL :: sXL, sXR
    ! the coefficients of the quadratic equaion determining lambdastar
    CCTK_REAL :: AL,AR,BL,BR,CL,CR,EL,ER

    ! the combination D*h*W which appear more than once
    CCTK_REAL :: DhWL,DhWR
    CCTK_REAL :: lamstar_tilde ! eq (22)


    ! compute sX
    DhWR = consm_i(1)+consm_i(5)+ &
           sqrt_avg_det*pressminus(i+xoffset,j+yoffset,k+zoffset)
    DhWL = consp  (1)+consp  (5)+ &
           sqrt_avg_det*pressplus (i        ,j        ,k        )
    if (flux_direction == 1) then
      sXR = uxxh*consm_i(2)+uxyh*consm_i(3)+uxzh*consm_i(4) - &
            avg_beta/avg_alp*DhWR
      sXL = uxxh*consp  (2)+uxyh*consp  (3)+uxzh*consp  (4) - &
            avg_beta/avg_alp*DhWL
    else if (flux_direction == 2) then
      sXR = uxyh*consm_i(2)+uyyh*consm_i(3)+uyzh*consm_i(4) - &
            avg_beta/avg_alp*DhWR
      sXL = uxyh*consp  (2)+uyyh*consp  (3)+uyzh*consp  (4) - &
            avg_beta/avg_alp*DhWL
    else if (flux_direction == 3) then
      sXR = uxzh*consm_i(2)+uyzh*consm_i(3)+uzzh*consm_i(4) - &
            avg_beta/avg_alp*DhWR
      sXL = uxzh*consp  (2)+uyzh*consp  (3)+uzzh*consp  (4) - &
            avg_beta/avg_alp*DhWL
    else
      call CCTK_ERROR("Flux direction not x,y,z")
      STOP
    end if

    ! compute the coefficients A-D of Duez (23)-(26)
    AR = (consm_i(5)+consm_i(1))*(charmax+avg_beta) - (avg_alp*sXR+DhWR*avg_beta)
    AL = (consp  (5)+consp  (1))*(charmin+avg_beta) - (avg_alp*sXL+DhWL*avg_beta)
    CR = avg_alp*usendh
    CL = CR
    BR = CR*sqrt_avg_det*pressminus(i+xoffset,j+yoffset,k+zoffset) - &
         (sXR+DhWR*avg_beta/avg_alp)*(charmax-vXR)
    BL = CL*sqrt_avg_det*pressplus (i        ,j        ,k        ) - &
         (sXL+DhWL*avg_beta/avg_alp)*(charmin-vXL)
    ER = -(avg_beta + charmax) ! sing ERROR in Duez (26)
    EL = -(avg_beta + charmin)

    ! we require the negative root only (see Mignone et al), this is (27) of
    ! Duez
    call solve_quadratic(EL*AR-AL*ER, -(AL*CR+BL*ER-CL*AR-EL*BR), BR*CL-BL*CR, &
                         lamstar_tilde)

    charstar = avg_alp * lamstar_tilde - avg_beta ! eq (22)
    ! pressstar is equal left and right of contact discontinuity
    pressstar = (AL*lamstar_tilde + BL)/((CL+EL*lamstar_tilde)*sqrt_avg_det)

  end subroutine

  ! this function is a straightforward transcription of the GSL solve_quadratic
  ! functions. Both codes are GPL. 
  ! Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
  subroutine solve_quadratic(a,b,c, res)

    implicit none

    CCTK_REAL, intent(in) :: a,b,c
    CCTK_REAL, intent(out) :: res

    CCTK_REAL, parameter :: one = 1
    CCTK_REAL :: disc, r, sgnb, temp, r1, r2
    character(len=256) :: warnline

    disc = b * b - 4 * a * c
  
    if (a == 0) then! Handle linear case
        if (b == 0) then
            write (warnline, '("Could not solve quadratic equation with coefficients ",3g22.15)') a,b,c
            call CCTK_ERROR(warnline)
            STOP
        else
            res = -c / b
        end if
        return
    end if
  
    if (disc > 0) then
        if (b == 0) then
            r = abs (0.5d0 * sqrt (disc) / a)
            res = -r
            !*x1 =  r ! only want negative root
        else
            sgnb = sign(one, b)
            temp = -0.5d0 * (b + sgnb * sqrt (disc))
            r1 = temp / a
            r2 = c / temp

            res = r2 ! only want negative (second) root
  
            !if (r1 < r2) then
            !    res = r1
            !    !*x1 = r2 ! only want negative root
            !else 
            !    res = r2
            !    !*x1 = r1 ! only want negative root
            !end if
        end if
    else if (disc == 0) then
        res = -0.5d0 * b / a
        !*x1 = -0.5 * b / a ! only want one root
    else
        write (warnline, '("Could not solve quadratic equation with coefficients ",3g22.15)') a,b,c
        call CCTK_ERROR(warnline)
        STOP
    end if

  end subroutine
end module GRHydro_hllc_helpers


subroutine GRHydro_HLLC(CCTK_ARGUMENTS)
  USE GRHydro_Eigenproblem
  USE GRHydro_hllc_helpers

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k
  CCTK_INT :: keytemp
  CCTK_REAL, dimension(5) :: consp,consm_i,fplus,fminus,lamplus
  CCTK_REAL, dimension(5) :: f1,lamminus, consstar
  CCTK_REAL ::  charmin, charmax, charstar, avg_alp,avg_det
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_beta, usendh
  CCTK_REAL :: pressstar ! NOTE: this name is unfortuante for MHD where it
                         ! conflicts
  CCTK_REAL :: sqrt_avg_det
  ! Duez und GRHydro define different quantities v^i, we have 
  ! v^i_Duez = alp*v^i_GRHydro - beta^i
  CCTK_REAL :: vXL, vXR
  CCTK_REAL :: velstar ! GRhydro's velocity in the star region
    
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

  !$OMP PARALLEL DO PRIVATE(k,j,i,f1,lamminus,lamplus,consp,consm_i,consstar, &
  !$OMP                     fplus,fminus,avg_beta,avg_alp,   &
  !$OMP                     avg_det,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
  !$OMP                     uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,         &
  !$OMP                     usendh, charmin, charmax, charstar, &
  !$OMP                     pressstar, sqrt_avg_det, vXL, vXR, velstar)
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
        sqrt_avg_det = sqrt(avg_det)

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
            vXR = avg_alp*velxminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta
            vXL = avg_alp*velxplus (i        ,j        ,k        ) - avg_beta
          else if (flux_direction == 2) then
            usendh = uyyh
            vXR = avg_alp*velyminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta
            vXL = avg_alp*velyplus (i        ,j        ,k        ) - avg_beta
          else if (flux_direction == 3) then
            usendh = uzzh
            vXR = avg_alp*velzminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta
            vXL = avg_alp*velzplus (i        ,j        ,k        ) - avg_beta
          else
            call CCTK_ERROR("Flux direction not x,y,z")
            STOP
          end if
          
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
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,&
                     int(i,ik),int(j,ik),int(k,ik),&
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
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,&
                     int(i,ik),int(j,ik),int(k,ik),&
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
                call eigenvalues_hot(GRHydro_eos_handle,keytemp,&
                     int(i,ik),int(j,ik),int(k,ik),&
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

          
!!$        we need characteristic speeds for leftmost and rightmost state only.
      
          ! Duez and GRHydro differ by a factor of lapse that we pull out of the
          ! fluxes. To aid the overall confusion, we are going to adhere to
          ! Duez' convetion for a while until we compute the fluxes. Note that
          ! the contained function compute_lamstar has access to all of our
          ! local variables and does indeed look at charmin and charmax.
          charmin = min(lamminus(1),lamplus(1))*avg_alp
          charmax = max(lamminus(5),lamplus(5))*avg_alp
          call compute_stardata(charstar, pressstar,consm_i,consp,sqrt_avg_det,pressminus,pressplus,&
                              uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,usendh,avg_beta,avg_alp,vxR,vxL,charmin,charmax,i,j,k,&
                              cctk_ash1,cctk_ash2,cctk_ash3,xoffset,yoffset,zoffset,flux_direction)
          
!!$        Calculate flux by standard formula

          ! TODO: possibly one can safe a bit of time but oly computing charstar
          ! when needed? No need to compute it when 0.lt.charmin or 0.gt.charmax
          ! TODO: only compute fplus and fminus when we need them
          ! TODO: clean this up. This is horrible. Consider reviving the GT C
          ! HLL code for this.
          if (0.lt.charmin) then
            f1 = fplus
          else if (0.lt.charstar) then
            ! solve for state vector ULstar using (12)-(16)
            ! note the change in sign in both numerator and denom
            consstar = (vXL - charmin)/(charstar-charmin)*consp
            ! sx, sy, sz
            consstar(1+flux_direction) = consstar(1+flux_direction) + &
              avg_alp*sqrt_avg_det* &
                (pressplus (i        ,j        ,k        )-pressstar) / &
                (charstar-charmin)
            ! tau
            consstar(5) = consstar(5) + sqrt_avg_det * &
              ((vXL+avg_beta)*pressplus (i        ,j        ,k        ) - &
               (charstar+avg_beta)*pressstar) / &
               (charstar-charmin)
            ! need to pass in GRHydro velocity computed from Duez velocity 
            ! which is charstar in the star region
            velstar = (charstar+avg_beta)/avg_alp
            if (flux_direction == 1) then
              call num_x_flux(consstar(1),consstar(2),consstar(3),&
                   consstar(4),consstar(5),&
                   f1(1),f1(2),f1(3),f1(4),f1(5),&
                   velstar,pressstar,avg_det,avg_alp,avg_beta)
            else if (flux_direction == 2) then
              call num_x_flux(consstar(1),consstar(3),consstar(4),&
                   consstar(2),consstar(5),&
                   f1(1),f1(3),f1(4),f1(2),f1(5),&
                   velstar,pressstar,avg_det,avg_alp,avg_beta)
            else if (flux_direction == 3) then
              call num_x_flux(consstar(1),consstar(4),consstar(2),&
                   consstar(3),consstar(5),&
                   f1(1),f1(4),f1(2),f1(3),f1(5),&
                   velstar,pressstar,avg_det,avg_alp,avg_beta)
            else
              call CCTK_ERROR("Flux direction not x,y,z")
              STOP
            end if
          else if (0.lt.charmax) then
            ! solve for state vector URstar using (12)-(16)
            ! note the change in sign in both numerator and denom
            consstar = (vXR - charmax)/(charstar-charmax)*consm_i
            ! sx, sy, sz
            consstar(1+flux_direction) = consstar(1+flux_direction) + &
              avg_alp*sqrt_avg_det* &
                (pressminus(i+xoffset,j+yoffset,k+zoffset)-pressstar) / &
                (charstar-charmax)
            ! tau
            consstar(5) = consstar(5) + sqrt_avg_det * &
              ((vXR+avg_beta)*pressminus(i+xoffset,j+yoffset,k+zoffset) - &
               (charstar+avg_beta)*pressstar) / &
               (charstar-charmax)
            ! need to pass in GRHydro velocity computed from Duez velocity 
            ! which is charstar in the star region
            velstar = (charstar+avg_beta)/avg_alp
            if (flux_direction == 1) then
              call num_x_flux(consstar(1),consstar(2),consstar(3),&
                   consstar(4),consstar(5),&
                   f1(1),f1(2),f1(3),f1(4),f1(5),&
                   velstar,pressstar,avg_det,avg_alp,avg_beta)
            else if (flux_direction == 2) then
              call num_x_flux(consstar(1),consstar(3),consstar(4),&
                   consstar(2),consstar(5),&
                   f1(1),f1(3),f1(4),f1(2),f1(5),&
                   velstar,pressstar,avg_det,avg_alp,avg_beta)
            else if (flux_direction == 3) then
              call num_x_flux(consstar(1),consstar(4),consstar(2),&
                   consstar(3),consstar(5),&
                   f1(1),f1(4),f1(2),f1(3),f1(5),&
                   velstar,pressstar,avg_det,avg_alp,avg_beta)
            else
              call CCTK_ERROR("Flux direction not x,y,z")
              STOP
            end if
          else
            f1 = fminus
          end if
          
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


end subroutine GRHydro_HLLC

