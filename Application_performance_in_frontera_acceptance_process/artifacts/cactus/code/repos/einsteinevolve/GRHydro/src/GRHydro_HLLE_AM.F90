 /*@@
   @file      GRHydro_HLLEPolyM.F90
   @date      Aug 30, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke, Pedro Montero, Toni Font
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
   @routine    GRHydro_HLLE_AM
   @date       Aug 30, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke, Pedro Montero, Toni Font
   @desc 
   The HLLE solver. Sufficiently simple that its just one big routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Altered from Cactus 3 routines originally written by Toni Font.
   @endhistory 

@@*/

subroutine GRHydro_HLLE_AM(CCTK_ARGUMENTS)
  USE GRHydro_EigenproblemM
  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k, m
  CCTK_REAL, dimension(8) :: cons_p,cons_m,fplus,fminus,f1,qdiff
  CCTK_REAL, dimension(10) :: prim_p, prim_m
  CCTK_REAL, dimension(5) :: lamminus,lamplus
  CCTK_REAL ::  charmin, charmax, charpm,chartop,avg_alp,avg_det, sdet
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_beta, usendh
  CCTK_REAL :: rhoenth_p, rhoenth_m, avg_betax, avg_betay, avg_betaz
  CCTK_REAL :: vxtp,vytp,vztp,vxtm,vytm,vztm,ab0p,ab0m,b2p,b2m,bdotvp,bdotvm
  CCTK_REAL :: wp,wm,v2p,v2m,bxlowp,bxlowm,bylowp,bylowm,bzlowp,bzlowm,vA2m,vA2p
  CCTK_REAL :: Bvecxlowp,Bvecxlowm,Bvecylowp,Bvecylowm,Bveczlowp,Bveczlowm
  CCTK_REAL :: pressstarp,pressstarm,velxlowp,velxlowm,velylowp,velylowm,velzlowp,velzlowm

  CCTK_INT :: type_bits, trivial
  CCTK_REAL :: xtemp

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
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define vel faulty_vel
#define Bvec faulty_Bvec

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
  end if

  ! constraint transport needs to be able to average fluxes in the directions
  ! other that flux_direction

  !$OMP PARALLEL DO PRIVATE(k,j,i,f1,lamminus,lamplus,cons_p,cons_m,fplus,fminus,qdiff,prim_p,prim_m,&
  !$OMP avg_betax,avg_betay,avg_betaz,avg_beta,avg_alp,&
  !$OMP gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,avg_det,sdet,uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,&
  !$OMP vxtp,vxtm,vytp,vytm,vztp,vztm,&
  !$OMP velxlowp,velxlowm,Bvecxlowp,Bvecxlowm,&
  !$OMP velylowp,velylowm,Bvecylowp,Bvecylowm,&
  !$OMP velzlowp,velzlowm,Bveczlowp,Bveczlowm,&
  !$OMP bdotvp,bdotvm,b2p,b2m,v2p,v2m,wp,wm,&
  !$OMP bxlowp,bxlowm,bylowp,bylowm,bzlowp,bzlowm,&
  !$OMP rhoenth_p,rhoenth_m,ab0p,ab0m,vA2p,vA2m,pressstarp,pressstarm,usendh,charmin,charmax,chartop,charpm,m,xtemp)
  do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + transport_constraints*(1-zoffset)
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + transport_constraints*(1-yoffset)
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + transport_constraints*(1-xoffset)
        
        f1 = 0.d0
        lamminus = 0.d0
        lamplus = 0.d0
        cons_p = 0.d0
        cons_m = 0.d0
        fplus = 0.d0
        fminus = 0.d0
        qdiff = 0.d0

!!$        Set the left (p for plus) and right (m_i for minus, i+1) states
        
        cons_p(1)   = densplus(i,j,k) 
        cons_p(2)   = sxplus(i,j,k)
        cons_p(3)   = syplus(i,j,k)
        cons_p(4)   = szplus(i,j,k)
        cons_p(5)   = tauplus(i,j,k)
        cons_p(6)   = Avecxplus(i,j,k)
        cons_p(7)   = Avecyplus(i,j,k)
        cons_p(8)   = Aveczplus(i,j,k)
        
        cons_m(1) = densminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(2) = sxminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(3) = syminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(4) = szminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(5) = tauminus(i+xoffset,j+yoffset,k+zoffset) 
        cons_m(6) = Avecxminus(i+xoffset,j+yoffset,k+zoffset) 
        cons_m(7) = Avecyminus(i+xoffset,j+yoffset,k+zoffset) 
        cons_m(8) = Aveczminus(i+xoffset,j+yoffset,k+zoffset) 

        prim_p(1)   = rhoplus(i,j,k) 
        prim_p(2)   = velxplus(i,j,k)
        prim_p(3)   = velyplus(i,j,k) 
        prim_p(4)   = velzplus(i,j,k)
        prim_p(5)   = epsplus(i,j,k)
        prim_p(6)   = pressplus(i,j,k)
        prim_p(7)   = w_lorentzplus(i,j,k)
        prim_p(8)   = Bvecxplus(i,j,k)
        prim_p(9)   = Bvecyplus(i,j,k) 
        prim_p(10)  = Bveczplus(i,j,k)
        
        prim_m(1) = rhominus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(2) = velxminus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(3) = velyminus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(4) = velzminus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(5) = epsminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(6) = pressminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(7) = w_lorentzminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(8) = Bvecxminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(9) = Bvecyminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(10)= Bveczminus(i+xoffset,j+yoffset,k+zoffset) 
          
!!$        Calculate various metric terms here.
!!$        Note also need the average of the lapse at the 
!!$        left and right points.
!!$        
!!$        In MHD, we need all three shift components regardless of the flux direction

        avg_betax = 0.5d0 * (beta1(i+xoffset,j+yoffset,k+zoffset) + &
             beta1(i,j,k))
        avg_betay = 0.5d0 * (beta2(i+xoffset,j+yoffset,k+zoffset) + &
             beta2(i,j,k))
        avg_betaz = 0.5d0 * (beta3(i+xoffset,j+yoffset,k+zoffset) + &
             beta3(i,j,k))
        if (flux_direction == 1) then
           avg_beta=avg_betax
        else if (flux_direction == 2) then
           avg_beta=avg_betay
        else if (flux_direction == 3) then
           avg_beta=avg_betaz
        else
           call CCTK_ERROR("Flux direction not x,y,z")
        end if

        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gxxh = 0.5d0 * (g11(i+xoffset,j+yoffset,k+zoffset) + g11(i,j,k))
        gxyh = 0.5d0 * (g12(i+xoffset,j+yoffset,k+zoffset) + g12(i,j,k))
        gxzh = 0.5d0 * (g13(i+xoffset,j+yoffset,k+zoffset) + g13(i,j,k))
        gyyh = 0.5d0 * (g22(i+xoffset,j+yoffset,k+zoffset) + g22(i,j,k))
        gyzh = 0.5d0 * (g23(i+xoffset,j+yoffset,k+zoffset) + g23(i,j,k))
        gzzh = 0.5d0 * (g33(i+xoffset,j+yoffset,k+zoffset) + g33(i,j,k))

        avg_det =  SPATIAL_DETERMINANT(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh)
        sdet = sqrt(avg_det)

        call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
             avg_det,gxxh, gxyh, gxzh, &
             gyyh, gyzh, gzzh)
          

        vxtp = prim_p(2)-avg_betax/avg_alp
        vytp = prim_p(3)-avg_betay/avg_alp
        vztp = prim_p(4)-avg_betaz/avg_alp
        vxtm = prim_m(2)-avg_betax/avg_alp
        vytm = prim_m(3)-avg_betay/avg_alp
        vztm = prim_m(4)-avg_betaz/avg_alp

        call calc_vlow_blow(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
             prim_p(2),prim_p(3),prim_p(4),prim_p(8),prim_p(9),prim_p(10), &
             velxlowp,velylowp,velzlowp,Bvecxlowp,Bvecylowp,Bveczlowp, &
             bdotvp,b2p,v2p,wp,bxlowp,bylowp,bzlowp)
        call calc_vlow_blow(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
             prim_m(2),prim_m(3),prim_m(4),prim_m(8),prim_m(9),prim_m(10), &
             velxlowm,velylowm,velzlowm,Bvecxlowm,Bvecylowm,Bveczlowm, &
             bdotvm,b2m,v2m,wm,bxlowm,bylowm,bzlowm)

        rhoenth_p = prim_p(1)*(1.d0+prim_p(5))+prim_p(6)
        rhoenth_m = prim_m(1)*(1.d0+prim_m(5))+prim_m(6)
        
        ab0p = wp*bdotvp
        ab0m = wm*bdotvm

        vA2p = b2p/(rhoenth_p+b2p)
        vA2m = b2m/(rhoenth_m+b2m)

!!$ p^*  = p+b^2/2  in Anton et al.
        pressstarp = prim_p(6)+0.5d0*b2p
        pressstarm = prim_m(6)+0.5d0*b2m        


!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell
          
        if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then
          
!!$ we now pass in the B-field as conserved and flux, the vtilde's instead of v's,
!!$ pressstar instead of P, b_i, alp b^0, w, metric determinant, 
!!$ alp, and beta in the flux dir

          if (flux_direction == 1) then
            call num_x_fluxM(cons_p(1),cons_p(2),cons_p(3),cons_p(4),cons_p(5),&
                 cons_p(6),cons_p(7),cons_p(8),&
                 f1(1),f1(2),f1(3),f1(4),f1(5),f1(6),f1(7),f1(8),&
                 vxtp,vytp,vztp,pressstarp,bxlowp,bylowp,bzlowp,ab0p,wp, &
                 avg_det,avg_alp,avg_beta)

          else if (flux_direction == 2) then
            call num_x_fluxM(cons_p(1),cons_p(3),cons_p(4),cons_p(2),cons_p(5),&
                 cons_p(7),cons_p(8),cons_p(6),&
                 f1(1),f1(3),f1(4),f1(2),f1(5),f1(7),f1(8),f1(6),&
                 vytp,vztp,vxtp,pressstarp,bylowp,bzlowp,bxlowp,ab0p,wp, &
                 avg_det,avg_alp,avg_beta)

          else if (flux_direction == 3) then
            call num_x_fluxM(cons_p(1),cons_p(4),cons_p(2),cons_p(3),cons_p(5),&
                 cons_p(8),cons_p(6),cons_p(7),&
                 f1(1),f1(4),f1(2),f1(3),f1(5),f1(8),f1(6),f1(7), &
                 vztp,vxtp,vytp,pressstarp,bzlowp,bxlowp,bylowp,ab0p,wp, &
                 avg_det,avg_alp,avg_beta)
            
          else
            call CCTK_ERROR("Flux direction not x,y,z")
          end if
          
        else !!! The end of this branch is right at the bottom of the routine
            
          if (flux_direction == 1) then
            usendh = uxxh
          else if (flux_direction == 2) then
            usendh = uyyh
          else if (flux_direction == 3) then
            usendh = uzzh
          else
            call CCTK_ERROR("Flux direction not x,y,z")
          end if
          
!!$        Calculate the jumps in the conserved variables
          
          qdiff(1) = cons_m(1) - cons_p(1)
          qdiff(2) = cons_m(2) - cons_p(2)
          qdiff(3) = cons_m(3) - cons_p(3)
          qdiff(4) = cons_m(4) - cons_p(4)
          qdiff(5) = cons_m(5) - cons_p(5)
          qdiff(6) = cons_m(6) - cons_p(6)
          qdiff(7) = cons_m(7) - cons_p(7)
          qdiff(8) = cons_m(8) - cons_p(8)

!!$        Eigenvalues and fluxes either side of the cell interface
          
          if (flux_direction == 1) then
            if(evolve_temper.ne.1) then
              call eigenvaluesM(GRHydro_eos_handle,&
                   prim_m(1),prim_m(2),prim_m(3),prim_m(4),prim_m(5),prim_m(6),prim_m(7), &
                   prim_m(8),prim_m(9),prim_m(10),&
                   lamminus,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
                   usendh,avg_alp,avg_beta)                                    
              call eigenvaluesM(GRHydro_eos_handle, &
                   prim_p(1),prim_p(2),prim_p(3),prim_p(4),prim_p(5),prim_p(6),prim_p(7), &
                   prim_p(8),prim_p(9),prim_p(10),&
                   lamplus,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
                   usendh,avg_alp,avg_beta)
            else
              xtemp = temperature(i,j,k)
              call eigenvaluesM_hot(GRHydro_eos_handle,&
                   int(i,ik),int(j,ik),int(k,ik),&
                   prim_m(1),prim_m(2),prim_m(3),prim_m(4),prim_m(5),prim_m(6),prim_m(7), &
                   prim_m(8),prim_m(9),prim_m(10),&
                   xtemp,y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                   lamminus,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
                   usendh,avg_alp,avg_beta)
              xtemp = temperature(i,j,k)
              call eigenvaluesM_hot(GRHydro_eos_handle,&
                   int(i,ik),int(j,ik),int(k,ik),&
                   prim_p(1),prim_p(2),prim_p(3),prim_p(4),prim_p(5),prim_p(6),prim_p(7), &
                   prim_p(8),prim_p(9),prim_p(10),&
                   xtemp,y_e_plus(i,j,k),&
                   lamplus,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
                   usendh,avg_alp,avg_beta)
            endif 

            call num_x_fluxM(cons_p(1),cons_p(2),cons_p(3),cons_p(4),cons_p(5),&
                 cons_p(6),cons_p(7),cons_p(8),&
                 fplus(1),fplus(2),fplus(3),fplus(4),fplus(5),fplus(6),fplus(7),fplus(8),&
                 vxtp,vytp,vztp,pressstarp,bxlowp,bylowp,bzlowp,ab0p,wp, &
                 avg_det,avg_alp,avg_beta)
            call num_x_fluxM(cons_m(1),cons_m(2),cons_m(3),cons_m(4),cons_m(5),&
                 cons_m(6),cons_m(7),cons_m(8),&
                 fminus(1),fminus(2),fminus(3),fminus(4),fminus(5),&
                 fminus(6),fminus(7),fminus(8),&
                 vxtm,vytm,vztm,pressstarm,bxlowm,bylowm,bzlowm,ab0m,wm, &
                 avg_det,avg_alp,avg_beta)

          else if (flux_direction == 2) then
            if(evolve_temper.ne.1) then
              call eigenvaluesM(GRHydro_eos_handle,&
                   prim_m(1),prim_m(3),prim_m(4),prim_m(2),prim_m(5),prim_m(6),prim_m(7), &
                   prim_m(9),prim_m(10),prim_m(8),&
                   lamminus,gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,&
                   usendh,avg_alp,avg_beta)                                    
              call eigenvaluesM(GRHydro_eos_handle, &
                   prim_p(1),prim_p(3),prim_p(4),prim_p(2),prim_p(5),prim_p(6),prim_p(7), &
                   prim_p(9),prim_p(10),prim_p(8),&
                   lamplus,gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,&
                   usendh,avg_alp,avg_beta)
            else
              xtemp = temperature(i,j,k)
              call eigenvaluesM_hot(GRHydro_eos_handle,&
                   int(i,ik),int(j,ik),int(k,ik),&
                   prim_m(1),prim_m(3),prim_m(4),prim_m(2),prim_m(5),prim_m(6),prim_m(7), &
                   prim_m(9),prim_m(10),prim_m(8),&
                   xtemp,y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                   lamminus,gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,&
                   usendh,avg_alp,avg_beta)                                    
              xtemp = temperature(i,j,k)                                       
              call eigenvaluesM_hot(GRHydro_eos_handle,&
                   int(i,ik),int(j,ik),int(k,ik),&
                   prim_p(1),prim_p(3),prim_p(4),prim_p(2),prim_p(5),prim_p(6),prim_p(7), &
                   prim_p(9),prim_p(10),prim_p(8),&
                   xtemp,y_e_plus(i,j,k),&
                   lamplus,gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,&
                   usendh,avg_alp,avg_beta)
            endif 

            call num_x_fluxM(cons_p(1),cons_p(3),cons_p(4),cons_p(2),cons_p(5),&
                 cons_p(7),cons_p(8),cons_p(6),&
                 fplus(1),fplus(3),fplus(4),fplus(2),fplus(5),fplus(7),fplus(8),fplus(6),&
                 vytp,vztp,vxtp,pressstarp,bylowp,bzlowp,bxlowp,ab0p,wp, &
                 avg_det,avg_alp,avg_beta)
            call num_x_fluxM(cons_m(1),cons_m(3),cons_m(4),cons_m(2),cons_m(5),&
                 cons_m(7),cons_m(8),cons_m(6),&
                 fminus(1),fminus(3),fminus(4),fminus(2),fminus(5),&
                 fminus(7),fminus(8),fminus(6),&
                 vytm,vztm,vxtm,pressstarm,bylowm,bzlowm,bxlowm,ab0m,wm, &
                 avg_det,avg_alp,avg_beta)

          else if (flux_direction == 3) then
            if(evolve_temper.ne.1) then
              call eigenvaluesM(GRHydro_eos_handle,&
                   prim_m(1),prim_m(4),prim_m(2),prim_m(3),prim_m(5),prim_m(6),prim_m(7), &
                   prim_m(10),prim_m(8),prim_m(9),&
                   lamminus,gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,&
                   usendh,avg_alp,avg_beta)                                    
              call eigenvaluesM(GRHydro_eos_handle, &
                   prim_p(1),prim_p(4),prim_p(2),prim_p(3),prim_p(5),prim_p(6),prim_p(7), &
                   prim_p(10),prim_p(8),prim_p(9),&
                   lamplus,gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,&
                   usendh,avg_alp,avg_beta)
            else
              xtemp = temperature(i,j,k)
              call eigenvaluesM_hot(GRHydro_eos_handle,&
                   int(i,ik),int(j,ik),int(k,ik),&
                   prim_m(1),prim_m(4),prim_m(2),prim_m(3),prim_m(5),prim_m(6),prim_m(7), &
                   prim_m(10),prim_m(8),prim_m(9),&
                   xtemp,y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                   lamminus,gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,&
                   usendh,avg_alp,avg_beta)                                    
              xtemp = temperature(i,j,k)                                       
              call eigenvaluesM_hot(GRHydro_eos_handle,&
                   int(i,ik),int(j,ik),int(k,ik),&
                   prim_p(1),prim_p(4),prim_p(2),prim_p(3),prim_p(5),prim_p(6),prim_p(7), &
                   prim_p(10),prim_p(8),prim_p(9),&
                   xtemp,y_e_plus(i,j,k),&
                   lamplus,gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,&
                   usendh,avg_alp,avg_beta)
            endif 

            call num_x_fluxM(cons_p(1),cons_p(4),cons_p(2),cons_p(3),cons_p(5),&
                 cons_p(8),cons_p(6),cons_p(7),&
                 fplus(1),fplus(4),fplus(2),fplus(3),fplus(5),fplus(8),fplus(6),fplus(7), &
                 vztp,vxtp,vytp,pressstarp,bzlowp,bxlowp,bylowp,ab0p,wp, &
                 avg_det,avg_alp,avg_beta)
            call num_x_fluxM(cons_m(1),cons_m(4),cons_m(2),cons_m(3),cons_m(5),&
                 cons_m(8),cons_m(6),cons_m(7),&
                 fminus(1),fminus(4),fminus(2),fminus(3),fminus(5), &
                 fminus(8),fminus(6),fminus(7), &
                 vztm,vxtm,vytm,pressstarm,bzlowm,bxlowm,bylowm,ab0m,wm, &
                 avg_det,avg_alp,avg_beta)

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

          chartop = max(-charmin,charmax)
          
          charpm = charmax - charmin

!!$        Calculate flux by standard formula
          
          do m = 1,8 
             
             qdiff(m) = cons_m(m) - cons_p(m)
             
            if (HLLE) then
              f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
                       charmax * charmin * qdiff(m)) / charpm
            else if (LLF) then
              f1(m) = 0.5d0 * (fplus(m) + fminus(m) - chartop * qdiff(m)) 
            end if 
             
          end do

        end if !!! The end of the SpaceMask check for a trivial RP.
            
        densflux(i, j, k) = f1(1)
        sxflux(i, j, k) = f1(2)
        syflux(i, j, k) = f1(3)
        szflux(i, j, k) = f1(4)
        tauflux(i, j, k) = f1(5)

        if ( evolve_Lorenz_gge.gt.0 ) then
           !! FIX: These aren't zero
           Avecxflux(i,j,k) = 0.d0
           Avecyflux(i,j,k) = 0.d0
           Aveczflux(i,j,k) = 0.d0
           Aphiflux(i,j,k) = 0.d0
        else
           Avecxflux(i,j,k) = 0.d0
           Avecyflux(i,j,k) = 0.d0
           Aveczflux(i,j,k) = 0.d0
        end if

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
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_betax
#undef faulty_betay
#undef faulty_betaz
#undef faulty_vel
#undef faulty_Bvec
      
end subroutine GRHydro_HLLE_AM

 /*@@
   @routine    GRHydro_HLLE_TracerAM
   @date       Aug 30, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   HLLE just for the tracer.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_HLLE_TracerAM(CCTK_ARGUMENTS)

  USE GRHydro_EigenproblemM

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k,  m
  CCTK_REAL, dimension(number_of_tracers) :: cons_p,cons_m,fplus,fminus,f1
  CCTK_REAL, dimension(5) :: lamminus,lamplus
  CCTK_REAL, dimension(number_of_tracers) :: qdiff
  CCTK_REAL, dimension(7) :: prim_p, prim_m
  CCTK_REAL, dimension(3) :: mag_p, mag_m
  CCTK_REAL ::  charmin, charmax, charpm,avg_alp,avg_det
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_beta, usendh
  CCTK_REAL :: b2p,b2m,vA2m,vA2p
    
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
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define vel faulty_vel
#define Bvec faulty_Bvec

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
  end if

  do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil
        
        f1 = 0.d0
        lamminus = 0.d0
        lamplus = 0.d0
        cons_p = 0.d0
        cons_m = 0.d0
        mag_p = 0.d0
        mag_m = 0.d0
        fplus = 0.d0
        fminus = 0.d0
        qdiff = 0.d0
        
!!$        Set the left (p for plus) and right (m_i for minus, i+1) states
        
        cons_p(:) = cons_tracerplus(i,j,k,:) 
        cons_m(:) = cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        
        mag_p(1) = Bvecxplus(i,j,k)
        mag_p(2) = Bvecyplus(i,j,k)
        mag_p(3) = Bveczplus(i,j,k)

        mag_m(1) = Bvecxminus(i+xoffset,j+yoffset,k+zoffset)
        mag_m(2) = Bvecyminus(i+xoffset,j+yoffset,k+zoffset)
        mag_m(3) = Bveczminus(i+xoffset,j+yoffset,k+zoffset)

        prim_p(1)   = rhoplus(i,j,k) 
        prim_p(2)   = velxplus(i,j,k)
        prim_p(3)   = velyplus(i,j,k) 
        prim_p(4)   = velzplus(i,j,k)
        prim_p(5)   = epsplus(i,j,k)
        prim_p(6)   = pressplus(i,j,k)
        prim_p(7)   = w_lorentzplus(i,j,k)
        
        prim_m(1) = rhominus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(2) = velxminus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(3) = velyminus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(4) = velzminus(i+xoffset,j+yoffset,k+zoffset)
        prim_m(5) = epsminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(6) = pressminus(i+xoffset,j+yoffset,k+zoffset) 
        prim_m(7) = w_lorentzminus(i+xoffset,j+yoffset,k+zoffset) 
        
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
        end if

        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gxxh = 0.5d0 * (g11(i+xoffset,j+yoffset,k+zoffset) + g11(i,j,k))
        gxyh = 0.5d0 * (g12(i+xoffset,j+yoffset,k+zoffset) + g12(i,j,k))
        gxzh = 0.5d0 * (g13(i+xoffset,j+yoffset,k+zoffset) + g13(i,j,k))
        gyyh = 0.5d0 * (g22(i+xoffset,j+yoffset,k+zoffset) + g22(i,j,k))
        gyzh = 0.5d0 * (g23(i+xoffset,j+yoffset,k+zoffset) + g23(i,j,k))
        gzzh = 0.5d0 * (g33(i+xoffset,j+yoffset,k+zoffset) + g33(i,j,k))

        avg_det =  SPATIAL_DETERMINANT(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh)

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
        end if

!!$        b^2 = (1-v^2)B^2+(B dot v)^2 
        b2p=DOTP2(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,mag_p(1),mag_p(2),mag_p(3))/prim_p(7)**2 + &
             (DOTP(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,prim_p(2),prim_p(3),prim_p(4),mag_p(1),mag_p(2),mag_p(3)))**2
        b2m=DOTP2(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,mag_m(1),mag_m(2),mag_m(3))/prim_m(7)**2 + &
             (DOTP(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,prim_m(2),prim_m(3),prim_m(4),mag_m(1),mag_m(2),mag_m(3)))**2
        
        vA2p = b2p/(prim_p(1)*(1.0d0+prim_p(5))+prim_p(6)+b2p)
        vA2m = b2m/(prim_m(1)*(1.0d0+prim_m(5))+prim_m(6)+b2m)
          
!!$        Calculate the jumps in the conserved variables
          
        qdiff = cons_m - cons_p
          
!!$        Eigenvalues and fluxes either side of the cell interface
          
        if (flux_direction == 1) then
           call eigenvaluesM(GRHydro_eos_handle,&
                prim_m(1),prim_m(2),prim_m(3),prim_m(4),prim_m(5),prim_m(6),prim_m(7), &
                mag_m(1),mag_m(2),mag_m(3),&
                lamminus,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
                usendh,avg_alp,avg_beta)                                    
           call eigenvaluesM(GRHydro_eos_handle, &
                prim_p(1),prim_p(2),prim_p(3),prim_p(4),prim_p(5),prim_p(6),prim_p(7), &
                mag_p(1),mag_p(2),mag_p(3),&
                lamplus,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,&
                usendh,avg_alp,avg_beta)
           fplus(:)  = (velxplus(i,j,k)  - avg_beta / avg_alp) * &
                cons_tracerplus(i,j,k,:)
           fminus(:) = (velxminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
                cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else if (flux_direction == 2) then
           call eigenvaluesM(GRHydro_eos_handle,&
                prim_m(1),prim_m(3),prim_m(4),prim_m(2),prim_m(5),prim_m(6),prim_m(7), &
                mag_m(2),mag_m(3),mag_m(1),&
                lamminus,gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,&
                usendh,avg_alp,avg_beta)                                    
           call eigenvaluesM(GRHydro_eos_handle, &
                prim_p(1),prim_p(3),prim_p(4),prim_p(2),prim_p(5),prim_p(6),prim_p(7), &
                mag_p(2),mag_p(3),mag_p(1),&
                lamplus,gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,&
                usendh,avg_alp,avg_beta)
           fplus(:)  = (velyplus(i,j,k)  - avg_beta / avg_alp) * &
                cons_tracerplus(i,j,k,:)
           fminus(:) = (velyminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
                cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else if (flux_direction == 3) then
           call eigenvaluesM(GRHydro_eos_handle,&
                prim_m(1),prim_m(4),prim_m(2),prim_m(3),prim_m(5),prim_m(6),prim_m(7), &
                mag_m(3),mag_m(1),mag_m(2),&
                lamminus,gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,&
                usendh,avg_alp,avg_beta)                                    
           call eigenvaluesM(GRHydro_eos_handle,&
                prim_p(1),prim_p(4),prim_p(2),prim_p(3),prim_p(5),prim_p(6),prim_p(7), &
                mag_p(3),mag_p(1),mag_p(2),&
                lamplus,gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,&
                usendh,avg_alp,avg_beta)
           fplus(:)  = (velzplus(i,j,k) - avg_beta / avg_alp) * &
                cons_tracerplus(i,j,k,:)
           fminus(:) = (velzminus(i+xoffset,j+yoffset,k+zoffset) - avg_beta / avg_alp) * &
                cons_tracerminus(i+xoffset,j+yoffset,k+zoffset,:)
        else
           call CCTK_ERROR("Flux direction not x,y,z")
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
           
           qdiff(m) = cons_m(m) - cons_p(m)
           
           f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
                charmax * charmin * qdiff(m)) / charpm
           
        end do
        
        cons_tracerflux(i, j, k,:) = f1(:)
!!$
!!$        if ( ((flux_direction.eq.3).and.(i.eq.4).and.(j.eq.4)).or.&
!!$             ((flux_direction.eq.2).and.(i.eq.4).and.(k.eq.4)).or.&
!!$             ((flux_direction.eq.1).and.(j.eq.4).and.(k.eq.4))&
!!$             ) then
!!$          write(*,*) flux_direction, i, j, k, f1(1), cons_m(1), cons_p(1)
!!$        end if
        
     end do
  end do
end do
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_betax
#undef faulty_betay
#undef faulty_betaz
#undef faulty_vel
#undef faulty_Bvec


end subroutine GRHydro_HLLE_TracerAM

