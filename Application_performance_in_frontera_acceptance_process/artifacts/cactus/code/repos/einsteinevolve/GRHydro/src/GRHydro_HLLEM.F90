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
   @routine    GRHydro_HLLEM
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

subroutine GRHydro_HLLEM(CCTK_ARGUMENTS)
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
  CCTK_REAL, dimension(3,3) :: gh, uh
  CCTK_REAL :: usendh
  CCTK_REAL :: rhoenth_p, rhoenth_m
  CCTK_REAL, dimension(3) :: avg_beta
  CCTK_REAL, dimension(3) :: vtp,vtm,blowp,blowm,Bveclowp,Bveclowm
  CCTK_REAL, dimension(3) :: vellowp,vellowm
  CCTK_REAL :: ab0p,ab0m,b2p,b2m,bdotvp,bdotvm
  CCTK_REAL :: wp,wm,v2p,v2m
  CCTK_REAL :: pressstarp,pressstarm

  CCTK_REAL :: entropyconsp,entropyconsm,entropyp,entropym,entropyf,entropydiff,entropyfp,entropyfm

  CCTK_REAL :: psidcp, psidcm, psidcf, psidcdiff, psidcfp, psidcfm
  CCTK_REAL :: charmax_dc, charmin_dc, charpm_dc   
 
  CCTK_INT :: type_bits, trivial
  CCTK_REAL :: xtemp

  integer :: ix,iy,iz

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
    ix = 1 ; iy = 2 ; iz = 3
  else if (flux_direction == 2) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemY")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemY", &
         &"trivial")
    ix = 2 ; iy = 3 ; iz = 1
  else if (flux_direction == 3) then
    call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemZ")
    call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemZ", &
         &"trivial")
    ix = 3 ; iy = 1 ; iz = 2
  else
    call CCTK_ERROR("Flux direction not x,y,z")
    STOP
  end if

  ! constraint transport needs to be able to average fluxes in the directions
  ! other that flux_direction

  !$OMP PARALLEL PRIVATE(k,j,i,f1,lamminus,lamplus,cons_p,cons_m,fplus,fminus,qdiff,psidcf,psidcp,psidcm,prim_p,prim_m,&
  !$OMP avg_beta,avg_alp,&
  !$OMP gh,avg_det,sdet,uh,&
  !$OMP vtp,vtm,vellowp,vellowm,Bveclowp,Bveclowm,&
  !$OMP bdotvp,bdotvm,b2p,b2m,v2p,v2m,wp,wm,&
  !$OMP blowp,blowm,&
  !$OMP rhoenth_p,rhoenth_m,ab0p,ab0m,pressstarp,pressstarm,&
  !$OMP usendh,psidcdiff,psidcfp,psidcfm,charmin,charmax,chartop,charpm,&
  !$OMP charmin_dc,charmax_dc,charpm_dc,m,xtemp,&
  !$OMP entropyconsp,entropyconsm,entropyp,entropym,entropyf,entropydiff,entropyfp,entropyfm)
  ! avoid compiler warnings
  psidcf = 0d0
  psidcp = 0.d0
  psidcm = 0d0
  psidcfp = 0d0
  psidcfm = 0d0

  entropyf = 0d0
  entropyfp = 0d0
  entropyfm = 0d0
  entropyconsm = 0d0
  entropyconsp = 0d0
  !$OMP DO
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
        cons_p(6)   = Bconsxplus(i,j,k)
        cons_p(7)   = Bconsyplus(i,j,k)
        cons_p(8)   = Bconszplus(i,j,k)
        
        cons_m(1) = densminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(2) = sxminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(3) = syminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(4) = szminus(i+xoffset,j+yoffset,k+zoffset)
        cons_m(5) = tauminus(i+xoffset,j+yoffset,k+zoffset) 
        cons_m(6) = Bconsxminus(i+xoffset,j+yoffset,k+zoffset) 
        cons_m(7) = Bconsyminus(i+xoffset,j+yoffset,k+zoffset) 
        cons_m(8) = Bconszminus(i+xoffset,j+yoffset,k+zoffset) 

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
          
        if(clean_divergence.ne.0) then
           psidcp = psidcplus(i,j,k)
           psidcm = psidcminus(i+xoffset,j+yoffset,k+zoffset)
        endif

        if(evolve_entropy.ne.0) then
           entropyp = entropyplus(i,j,k)
           entropym = entropyminus(i+xoffset,j+yoffset,k+zoffset)
           entropyconsp = entropyconsplus(i,j,k)
           entropyconsm = entropyconsminus(i+xoffset,j+yoffset,k+zoffset)
        endif

!!$        Calculate various metric terms here.
!!$        Note also need the average of the lapse at the 
!!$        left and right points.
!!$        
!!$        In MHD, we need all three shift components regardless of the flux direction

        avg_beta(1) = 0.5d0 * (beta1(i+xoffset,j+yoffset,k+zoffset) + &
             beta1(i,j,k))
        avg_beta(2) = 0.5d0 * (beta2(i+xoffset,j+yoffset,k+zoffset) + &
             beta2(i,j,k))
        avg_beta(3) = 0.5d0 * (beta3(i+xoffset,j+yoffset,k+zoffset) + &
             beta3(i,j,k))

        avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        
        gh(1,1) = 0.5d0 * (g11(i+xoffset,j+yoffset,k+zoffset) + g11(i,j,k))
        gh(1,2) = 0.5d0 * (g12(i+xoffset,j+yoffset,k+zoffset) + g12(i,j,k))
        gh(1,3) = 0.5d0 * (g13(i+xoffset,j+yoffset,k+zoffset) + g13(i,j,k))
        gh(2,2) = 0.5d0 * (g22(i+xoffset,j+yoffset,k+zoffset) + g22(i,j,k))
        gh(2,3) = 0.5d0 * (g23(i+xoffset,j+yoffset,k+zoffset) + g23(i,j,k))
        gh(3,3) = 0.5d0 * (g33(i+xoffset,j+yoffset,k+zoffset) + g33(i,j,k))
        gh(2,1) = gh(1,2) ; gh(3,1) = gh(1,3) ; gh(3,2) = gh(2,3)

        avg_det =  SPATIAL_DETERMINANT(gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3))
        sdet = sqrt(avg_det)

        call UpperMetric(uh(1,1), uh(1,2), uh(1,3), uh(2,2), uh(2,3), uh(3,3), &
             avg_det,gh(1,1), gh(1,2), gh(1,3), &
             gh(2,2), gh(2,3), gh(3,3))
        uh(2,1) = uh(1,2) ; uh(3,1) = uh(1,3) ; uh(3,2) = uh(2,3)
          

        vtp(1) = prim_p(2)-avg_beta(1)/avg_alp
        vtp(2) = prim_p(3)-avg_beta(2)/avg_alp
        vtp(3) = prim_p(4)-avg_beta(3)/avg_alp
        vtm(1) = prim_m(2)-avg_beta(1)/avg_alp
        vtm(2) = prim_m(3)-avg_beta(2)/avg_alp
        vtm(3) = prim_m(4)-avg_beta(3)/avg_alp

        call calc_vlow_blow(gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3), &
             prim_p(2),prim_p(3),prim_p(4),prim_p(8),prim_p(9),prim_p(10), &
             vellowp(1),vellowp(2),vellowp(3),Bveclowp(1),Bveclowp(2),Bveclowp(3), &
             bdotvp,b2p,v2p,wp,blowp(1),blowp(2),blowp(3))
        call calc_vlow_blow(gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3), &
             prim_m(2),prim_m(3),prim_m(4),prim_m(8),prim_m(9),prim_m(10), &
             vellowm(1),vellowm(2),vellowm(3),Bveclowm(1),Bveclowm(2),Bveclowm(3), &
             bdotvm,b2m,v2m,wm,blowm(1),blowm(2),blowm(3))

        rhoenth_p = prim_p(1)*(1.d0+prim_p(5))+prim_p(6)
        rhoenth_m = prim_m(1)*(1.d0+prim_m(5))+prim_m(6)
        
        ab0p = wp*bdotvp
        ab0m = wm*bdotvm

!!$ p^*  = p+b^2/2  in Anton et al.
        pressstarp = prim_p(6)+0.5d0*b2p
        pressstarm = prim_m(6)+0.5d0*b2m        


!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell
          
        if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then
          
!!$ we now pass in the B-field as conserved and flux, the vtilde's instead of v's,
!!$ pressstar instead of P, b_i, alp b^0, w, metric determinant, 
!!$ alp, and beta in the flux dir

          call num_x_fluxM(cons_p(1),cons_p(1+ix),cons_p(1+iy),cons_p(1+iz),&
               cons_p(5),cons_p(5+ix),cons_p(5+iy),cons_p(5+iz),&
               f1(1),f1(1+ix),f1(1+iy),f1(1+iz),f1(5),f1(5+ix),f1(5+iy),&
               f1(5+iz),&
               vtp(ix),vtp(iy),vtp(iz),pressstarp,blowp(ix),blowp(iy),&
               blowp(iz),ab0p,wp, &
               avg_det,avg_alp,avg_beta(flux_direction))
          if(clean_divergence.ne.0) then
             f1(6)=f1(6) + 1.0d0*sdet*uh(flux_direction,1)*psidcp - cons_p(5+flux_direction)*avg_beta(1)/avg_alp
             f1(7)=f1(7) + 1.0d0*sdet*uh(flux_direction,2)*psidcp - cons_p(5+flux_direction)*avg_beta(2)/avg_alp
             f1(8)=f1(8) + 1.0d0*sdet*uh(flux_direction,3)*psidcp - cons_p(5+flux_direction)*avg_beta(3)/avg_alp
             psidcf = cons_p(5+flux_direction)/sdet-psidcp*avg_beta(flux_direction)/avg_alp
          endif
          if(evolve_entropy.ne.0) then
             entropyf = entropyconsp*vtp(flux_direction)
          endif
          
        else !!! The end of this branch is right at the bottom of the routine
            
          usendh = uh(flux_direction,flux_direction)
          
!!$        Calculate the jumps in the conserved variables
          
          qdiff = cons_m - cons_p

          if (clean_divergence.ne.0) then
             psidcdiff = psidcm - psidcp
          endif
          if(evolve_entropy.ne.0) then
             entropydiff = entropyconsm - entropyconsp
          endif

!!$       Eigenvalues and fluxes either side of the cell interface
          
          if(evolve_temper.ne.1) then
            call eigenvaluesM(GRHydro_eos_handle,&
                 prim_m(1),prim_m(1+ix),prim_m(1+iy),prim_m(1+iz),prim_m(5),prim_m(6),prim_m(7), &
                 prim_m(7+ix),prim_m(7+iy),prim_m(7+iz),&
                 lamminus,gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3),&
                 usendh,avg_alp,avg_beta(flux_direction))
            call eigenvaluesM(GRHydro_eos_handle,&
                 prim_p(1),prim_p(1+ix),prim_p(1+iy),prim_p(1+iz),prim_p(5),prim_p(6),prim_p(7), &
                 prim_p(7+ix),prim_p(7+iy),prim_p(7+iz),&
                 lamplus,gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3),&
                 usendh,avg_alp,avg_beta(flux_direction))
          else
            xtemp = temperature(i,j,k)
            call eigenvaluesM_hot(GRHydro_eos_handle,&
                 int(i,ik),int(j,ik),int(k,ik),&
                 prim_m(1),prim_m(1+ix),prim_m(1+iy),prim_m(1+iz),prim_m(5),prim_m(6),prim_m(7), &
                 prim_m(7+ix),prim_m(7+iy),prim_m(7+iz),&
                 xtemp,y_e_minus(i+xoffset,j+yoffset,k+zoffset),&
                 lamminus,gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3),&
                 usendh,avg_alp,avg_beta(flux_direction))
            xtemp = temperature(i,j,k)
            call eigenvaluesM_hot(GRHydro_eos_handle,&
                 int(i,ik),int(j,ik),int(k,ik),&
                 prim_p(1),prim_p(1+ix),prim_p(1+iy),prim_p(1+iz),prim_p(5),prim_p(6),prim_p(7), &
                 prim_p(7+ix),prim_p(7+iy),prim_p(7+iz),&
                 xtemp,y_e_plus(i,j,k),&
                 lamplus,gh(1,1),gh(1,2),gh(1,3),gh(2,2),gh(2,3),gh(3,3),&
                 usendh,avg_alp,avg_beta(flux_direction))
          endif 

          call num_x_fluxM(cons_p(1),cons_p(1+ix),cons_p(1+iy),cons_p(1+iz),&
               cons_p(5),&
               cons_p(5+ix),cons_p(5+iy),cons_p(5+iz),&
               fplus(1),fplus(1+ix),fplus(1+iy),fplus(1+iz),fplus(5),&
               fplus(5+ix),fplus(5+iy),fplus(5+iz),&
               vtp(ix),vtp(iy),vtp(iz),pressstarp,blowp(ix),blowp(iy),&
               blowp(iz),ab0p,wp, &
               avg_det,avg_alp,avg_beta(flux_direction))
          call num_x_fluxM(cons_m(1),cons_m(1+ix),cons_m(1+iy),cons_m(1+iz),&
               cons_m(5),&
               cons_m(5+ix),cons_m(5+iy),cons_m(5+iz),&
               fminus(1),fminus(1+ix),fminus(1+iy),fminus(1+iz),fminus(5),&
               fminus(5+ix),fminus(5+iy),fminus(5+iz),&
               vtm(ix),vtm(iy),vtm(iz),pressstarm,blowm(ix),blowm(iy),&
               blowm(iz),ab0m,wm, &
               avg_det,avg_alp,avg_beta(flux_direction))

          if(clean_divergence.ne.0) then
             fminus(6)=fminus(6) + 1.0d0*sdet*uh(flux_direction,1)*psidcm - cons_m(5+flux_direction)*avg_beta(1)/avg_alp
             fminus(7)=fminus(7) + 1.0d0*sdet*uh(flux_direction,2)*psidcm - cons_m(5+flux_direction)*avg_beta(2)/avg_alp
             fminus(8)=fminus(8) + 1.0d0*sdet*uh(flux_direction,3)*psidcm - cons_m(5+flux_direction)*avg_beta(3)/avg_alp
             fplus(6)=fplus(6) + 1.0d0*sdet*uh(flux_direction,1)*psidcp - cons_p(5+flux_direction)*avg_beta(1)/avg_alp
             fplus(7)=fplus(7) + 1.0d0*sdet*uh(flux_direction,2)*psidcp - cons_p(5+flux_direction)*avg_beta(2)/avg_alp
             fplus(8)=fplus(8) + 1.0d0*sdet*uh(flux_direction,3)*psidcp - cons_p(5+flux_direction)*avg_beta(3)/avg_alp
             psidcfp = cons_p(5+flux_direction)/sdet-avg_beta(flux_direction)*psidcp/avg_alp
             psidcfm = cons_m(5+flux_direction)/sdet-avg_beta(flux_direction)*psidcm/avg_alp
          endif
          if(evolve_entropy.ne.0) then
             entropyfp = entropyconsp*vtp(flux_direction)
             entropyfm = entropyconsm*vtm(flux_direction)
          endif

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

          if(clean_divergence.ne.0) then

           psidcdiff = psidcm - psidcp
           select case(whichpsidcspeed)
             case(0)
               if (HLLE) then
                 psidcf = (charmax * psidcfp - charmin * psidcfm + &
                           charmax * charmin * psidcdiff) / charpm
               else if (LLF) then
                 psidcf = 0.5d0 * (psidcfp + psidcfm - chartop * psidcdiff)
               end if
             case(1)
             !!$ Wavespeeds for psidc are +/-c, not Fast Magnetosonic?
             !!$ psidcf = 0.5d0 * (1.d0 * psidcfp - (-1.d0) * psidcfm + &
             !!$      1.d0 * (-1.d0) * psidcdiff) 

             !!$ The fastest speed for psidc comes from the condition
             !!$ that the normal vector to the characteristic hypersurface
             !!$ be spacelike (Eq. 60 of Anton et al.)

             charmax_dc = sqrt(usendh) - avg_beta(flux_direction)/avg_alp
             charmin_dc = -1.d0*sqrt(usendh) - avg_beta(flux_direction)/avg_alp
             charpm_dc = charmax_dc - charmin_dc

             psidcf = (charmax_dc * psidcfp - charmin_dc * psidcfm + &
                       charmax_dc * charmin_dc * psidcdiff) / charpm_dc

             if(decouple_normal_Bfield .ne. 0) then ! same expression for HLLE and LLF
               !!$ B^i field decouples from the others and has same propagation
               !!$ speed as divergence -null direction,  
               !!$ \pm sqrt(g^{xx}} - beta^x/alpha
               f1(5+flux_direction) = (charmax_dc * fplus(5+flux_direction) &
                                     - charmin_dc * fminus(5+flux_direction) + &
                  charmax_dc * charmin_dc * qdiff(5+flux_direction)) / charpm_dc
             end if
             case(2)
             charmax = setcharmax
             charmin = setcharmin
               if (HLLE) then
                 psidcf = (charmax * psidcfp - charmin * psidcfm + &
                           charmax * charmin * psidcdiff) / charpm
               else if (LLF) then
                 chartop = max(-charmin,charmax)
                 psidcf = 0.5d0 * (psidcfp + psidcfm - chartop * psidcdiff)
               end if
           end select
          endif
          
          if(evolve_entropy.ne.0) then
            entropydiff = entropyconsm - entropyconsp
            if (HLLE) then
              entropyf = (charmax * entropyfp - charmin * entropyfm + &
                       charmax * charmin * entropydiff) / charpm
            else if (LLF) then
              entropyf = 0.5d0 * (entropyfp + entropyfm - chartop * entropydiff) 
            end if 
          endif
             
             

        end if !!! The end of the SpaceMask check for a trivial RP.
            
        densflux(i, j, k) = f1(1)
        sxflux(i, j, k) = f1(2)
        syflux(i, j, k) = f1(3)
        szflux(i, j, k) = f1(4)
        tauflux(i, j, k) = f1(5)

        Bconsxflux(i, j, k) = f1(6)
        Bconsyflux(i, j, k) = f1(7)
        Bconszflux(i, j, k) = f1(8)

        if(clean_divergence.ne.0) then
           psidcflux(i,j,k) = psidcf
        endif

        if(evolve_entropy.ne.0) then
           entropyflux(i,j,k) = entropyf
        endif

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
  !$OMP END DO
  !$OMP END PARALLEL
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
      
end subroutine GRHydro_HLLEM

 /*@@
   @routine    GRHydro_HLLE_TracerM
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

subroutine GRHydro_HLLE_TracerM(CCTK_ARGUMENTS)

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
  CCTK_REAL :: b2p,b2m
    
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
    STOP
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
           STOP
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
          STOP
        end if

!!$        b^2 = (1-v^2)B^2+(B dot v)^2 
        b2p=DOTP2(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,mag_p(1),mag_p(2),mag_p(3))/prim_p(7)**2 + &
             (DOTP(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,prim_p(2),prim_p(3),prim_p(4),mag_p(1),mag_p(2),mag_p(3)))**2
        b2m=DOTP2(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,mag_m(1),mag_m(2),mag_m(3))/prim_m(7)**2 + &
             (DOTP(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,prim_m(2),prim_m(3),prim_m(4),mag_m(1),mag_m(2),mag_m(3)))**2
        
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


end subroutine GRHydro_HLLE_TracerM

