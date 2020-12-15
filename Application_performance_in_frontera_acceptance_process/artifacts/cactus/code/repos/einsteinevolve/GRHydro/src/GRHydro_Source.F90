 /*@@
   @file      GRHydro_Source.F90
   @date      Sat Jan 26 02:03:56 2002
   @author    Ian Hawke
   @desc 
   The geometric source terms for the matter evolution
   @enddesc 
 @@*/

! Second order f.d.

#define DIFF_X_2(q) (0.5d0 * (q(i+1,j,k) - q(i-1,j,k)) * ida)
#define DIFF_Y_2(q) (0.5d0 * (q(i,j+1,k) - q(i,j-1,k)) * idb)
#define DIFF_Z_2(q) (0.5d0 * (q(i,j,k+1) - q(i,j,k-1)) * idc)

! Fourth order f.d.

#if(1)
#define DIFF_X_4(q) (((q(i-2,j,k)-q(i+2,j,k)) + 8.d0 * (q(i+1,j,k) - q(i-1,j,k))) / 12.d0 * ida)
#define DIFF_Y_4(q) (((q(i,j-2,k)-q(i,j+2,k)) + 8.d0 * (q(i,j+1,k) - q(i,j-1,k))) / 12.d0 * idb)
#define DIFF_Z_4(q) (((q(i,j,k-2)-q(i,j,k+2)) + 8.d0 * (q(i,j,k+1) - q(i,j,k-1))) / 12.d0 * idc)
#else
#define DIFF_X_4(q) ((-q(i+2,j,k) + 8.d0 * q(i+1,j,k) - 8.d0 * q(i-1,j,k) + q(i-2,j,k)) / 12.d0 * ida)
#define DIFF_Y_4(q) ((-q(i,j+2,k) + 8.d0 * q(i,j+1,k) - 8.d0 * q(i,j-1,k) + q(i,j-2,k)) / 12.d0 * idb)
#define DIFF_Z_4(q) ((-q(i,j,k+2) + 8.d0 * q(i,j,k+1) - 8.d0 * q(i,j,k-1) + q(i,j,k-2)) / 12.d0 * idc)
#endif

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#define vela(i,j,k) vup(i,j,k,1)
#define velb(i,j,k) vup(i,j,k,2)
#define velc(i,j,k) vup(i,j,k,3)
 /*@@
   @routine    SourceTerms
   @date       Sat Jan 26 02:04:21 2002
   @author     Ian Hawke
   @desc 
   Calculate the geometric source terms and add to the update GFs
   @enddesc 
   @calls     
   @calledby   
   @history 
   Minor alterations of routine from GR3D.
   @endhistory 

@@*/

subroutine SourceTerms(CCTK_ARGUMENTS)
      
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, na, nb, nc
  CCTK_REAL :: one, two, half
  CCTK_REAL :: t00, t0a, t0b, t0c, taa, tab, tac, tbb, tbc, tcc
  CCTK_REAL :: uaa, uab, uac, ubb, ubc, ucc, rhoenthalpyW2
  CCTK_REAL :: shifta, shiftb, shiftc, velashift, velbshift, velcshift 
  CCTK_REAL :: vlowa, vlowb, vlowc
  CCTK_REAL :: da_betaa, da_betab, da_betac, db_betaa, db_betab,&
       db_betac, dc_betaa, dc_betab, dc_betac
  CCTK_REAL :: da_alp, db_alp, dc_alp
  CCTK_REAL :: tau_source, sa_source, sb_source, sc_source
  CCTK_REAL :: localgaa,localgab,localgac,localgbb,localgbc,localgcc
  CCTK_REAL :: da_gaa, da_gab, da_gac, da_gbb, da_gbc, da_gcc
  CCTK_REAL :: db_gaa, db_gab, db_gac, db_gbb, db_gbc, db_gcc
  CCTK_REAL :: dc_gaa, dc_gab, dc_gac, dc_gbb, dc_gbc, dc_gcc
  CCTK_REAL :: da, db, dc, ida, idb, idc
  CCTK_REAL :: shiftshiftk, shiftka, shiftkb, shiftkc
  CCTK_REAL :: sumTK
  CCTK_REAL :: halfshiftdga, halfshiftdgb, halfshiftdgc
  CCTK_REAL :: halfTdga, halfTdgb, halfTdgc
  CCTK_REAL :: invalp, invalp2
  CCTK_INT  :: local_spatial_order

  logical, allocatable, dimension (:,:,:) :: force_spatial_second_order

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: k11, k12, k13, k22, k23, k33
  pointer (pk11,k11), (pk12,k12), (pk13,k13), (pk22,k22), (pk23,k23), (pk33,k33)
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
    pk11 = loc(kaa)
    pk12 = loc(kab)
    pk13 = loc(kac)
    pk22 = loc(kbb)
    pk23 = loc(kbc)
    pk33 = loc(kcc)
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
    pk11 = loc(kxx)
    pk12 = loc(kxy)
    pk13 = loc(kxz)
    pk22 = loc(kyy)
    pk23 = loc(kyz)
    pk33 = loc(kzz)
    pbeta1 = loc(betax)
    pbeta2 = loc(betay)
    pbeta3 = loc(betaz)
    pvup = loc(vel)
  end if

  one = 1.0d0
  two = 2.0d0
  half = 0.5d0
  na = cctk_lsh(1)
  nb = cctk_lsh(2)
  nc = cctk_lsh(3)
  da = CCTK_DELTA_SPACE(1)
  db = CCTK_DELTA_SPACE(2)
  dc = CCTK_DELTA_SPACE(3)
  ida = 1.d0/da
  idb = 1.d0/db
  idc = 1.d0/dc
 
!!$  Initialize the update terms to be zero.
!!$  This will guarantee that no garbage in the boundaries is updated.

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0

  if (evolve_tracer .ne. 0) then
    cons_tracerrhs = 0.0d0
  end if

  if (evolve_Y_e .ne. 0) then
     y_e_con_rhs = 0.0d0
  endif

!!$  Set up the array for checking the order. We switch to second order
!!$  differencing at boundaries and near excision regions.
!!$  Copied straight from BSSN.

  allocate (force_spatial_second_order(na,nb,nc))
  force_spatial_second_order = .FALSE.
  
  if (sources_spatial_order > 2) then
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1 + GRHydro_stencil, nc - GRHydro_stencil
      do j = 1 + GRHydro_stencil, nb - GRHydro_stencil
        do i = 1 + GRHydro_stencil, na - GRHydro_stencil
          if ((i < 3).or.(i > cctk_lsh(1) - 2).or. &
               (j < 3).or.(j > cctk_lsh(2) - 2).or. &
               (k < 3).or.(k > cctk_lsh(3) - 2) ) then
            force_spatial_second_order(i,j,k) = .TRUE.
          else if ( use_mask > 0 ) then
            if (minval(emask(i-2:i+2,j-2:j+2,k-2:k+2)) < 0.75d0) then
              force_spatial_second_order(i,j,k) = .TRUE.
            end if
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if
  
  !$OMP PARALLEL DO PRIVATE(i, j, k, local_spatial_order,&
  !$OMP localgaa,localgab,localgac,localgbb,localgbc,localgcc,&
  !$OMP rhoenthalpyW2,shifta,shiftb,shiftc,&
  !$OMP da_betaa,da_betab,da_betac,db_betaa,db_betab,db_betac,&
  !$OMP dc_betaa,dc_betab,dc_betac,velashift,velbshift,velcshift,&
  !$OMP vlowa,vlowb,vlowc,t00,t0a,t0b,t0c,taa,tab,tac,tbb,tbc,tcc,&
  !$OMP da_alp,db_alp,dc_alp,tau_source,sa_source,sb_source,sc_source,&
  !$OMP uaa, uab, uac, ubb, ubc, ucc,&
  !$OMP da_gaa, da_gab, da_gac, da_gbb, da_gbc, da_gcc,&
  !$OMP db_gaa, db_gab, db_gac, db_gbb, db_gbc, db_gcc,&
  !$OMP dc_gaa, dc_gab, dc_gac, dc_gbb, dc_gbc, dc_gcc,&
  !$OMP shiftshiftk,shiftka,shiftkb,shiftkc,&
  !$OMP sumTK,halfshiftdga,halfshiftdgb,halfshiftdgc,&
  !$OMP halfTdga,halfTdgb,halfTdgc,invalp,invalp2)
  do k=1 + GRHydro_stencil,nc - GRHydro_stencil
    do j=1 + GRHydro_stencil,nb - GRHydro_stencil
      do i=1 + GRHydro_stencil,na - GRHydro_stencil

        local_spatial_order = sources_spatial_order
        if (force_spatial_second_order(i,j,k)) then
          local_spatial_order = 2
        end if
        
!!$        Set the metric terms.

        localgaa = g11(i,j,k)
        localgab = g12(i,j,k)
        localgac = g13(i,j,k)
        localgbb = g22(i,j,k)
        localgbc = g23(i,j,k)
        localgcc = g33(i,j,k)

        call UpperMetric(uaa, uab, uac, ubb, ubc, ucc, &
             sdetg(i,j,k)*sdetg(i,j,k), localgaa,&
             localgab, localgac, localgbb, localgbc, localgcc)
        
!!$        All the matter variables (except velocity) always appear
!!$        together in this form

        rhoenthalpyW2 = (rho(i,j,k)*(one + eps(i,j,k)) + press(i,j,k))*&
             w_lorentz(i,j,k)**2
        
        shifta = beta1(i,j,k)
        shiftb = beta2(i,j,k)
        shiftc = beta3(i,j,k)

        if (local_spatial_order .eq. 2) then
           
           da_betaa = DIFF_X_2(beta1)
           da_betab = DIFF_X_2(beta2)
           da_betac = DIFF_X_2(beta3)
           
           db_betaa = DIFF_Y_2(beta1)
           db_betab = DIFF_Y_2(beta2)
           db_betac = DIFF_Y_2(beta3)
           
           dc_betaa = DIFF_Z_2(beta1)
           dc_betab = DIFF_Z_2(beta2)
           dc_betac = DIFF_Z_2(beta3)
           
        else

           da_betaa = DIFF_X_4(beta1)
           da_betab = DIFF_X_4(beta2)
           da_betac = DIFF_X_4(beta3)
           
           db_betaa = DIFF_Y_4(beta1)
           db_betab = DIFF_Y_4(beta2)
           db_betac = DIFF_Y_4(beta3)
            
           dc_betaa = DIFF_Z_4(beta1)
           dc_betab = DIFF_Z_4(beta2)
           dc_betac = DIFF_Z_4(beta3)
           
        end if
          
        invalp = 1.0d0 / alp(i,j,k)
        invalp2 = invalp**2
        velashift = vela(i,j,k) - shifta*invalp
        velbshift = velb(i,j,k) - shiftb*invalp
        velcshift = velc(i,j,k) - shiftc*invalp
        vlowa = vela(i,j,k)*localgaa + velb(i,j,k)*localgab +&
             velc(i,j,k)*localgac
        vlowb = vela(i,j,k)*localgab + velb(i,j,k)*localgbb +&
             velc(i,j,k)*localgbc
        vlowc = vela(i,j,k)*localgac + velb(i,j,k)*localgbc +&
             velc(i,j,k)*localgcc

!!$        For a change, these are T^{ij}

        t00 = (rhoenthalpyW2 - press(i,j,k))*invalp2
        t0a = rhoenthalpyW2*velashift/alp(i,j,k) +&
             press(i,j,k)*shifta*invalp2
        t0b = rhoenthalpyW2*velbshift/alp(i,j,k) +&
             press(i,j,k)*shiftb*invalp2
        t0c = rhoenthalpyW2*velcshift/alp(i,j,k) +&
             press(i,j,k)*shiftc*invalp2
        taa = rhoenthalpyW2*velashift*velashift +&
             press(i,j,k)*(uaa - shifta*shifta*invalp2)
        tab = rhoenthalpyW2*velashift*velbshift +&
             press(i,j,k)*(uab - shifta*shiftb*invalp2)
        tac = rhoenthalpyW2*velashift*velcshift +&
             press(i,j,k)*(uac - shifta*shiftc*invalp2)
        tbb = rhoenthalpyW2*velbshift*velbshift +&
             press(i,j,k)*(ubb - shiftb*shiftb*invalp2)
        tbc = rhoenthalpyW2*velbshift*velcshift +&
             press(i,j,k)*(ubc - shiftb*shiftc*invalp2)
        tcc = rhoenthalpyW2*velcshift*velcshift +&
             press(i,j,k)*(ucc - shiftc*shiftc*invalp2)

!!$        Derivatives of the lapse, metric and shift

        if (local_spatial_order .eq. 2) then

          da_alp = DIFF_X_2(alp)
          db_alp = DIFF_Y_2(alp)
          dc_alp = DIFF_Z_2(alp)

        else

          da_alp = DIFF_X_4(alp)
          db_alp = DIFF_Y_4(alp)
          dc_alp = DIFF_Z_4(alp)

        end if
        
        if (local_spatial_order .eq. 2) then

           da_gaa = DIFF_X_2(g11)
           da_gab = DIFF_X_2(g12)
           da_gac = DIFF_X_2(g13)
           da_gbb = DIFF_X_2(g22)
           da_gbc = DIFF_X_2(g23)
           da_gcc = DIFF_X_2(g33)
           db_gaa = DIFF_Y_2(g11)
           db_gab = DIFF_Y_2(g12)
           db_gac = DIFF_Y_2(g13)
           db_gbb = DIFF_Y_2(g22)
           db_gbc = DIFF_Y_2(g23)
           db_gcc = DIFF_Y_2(g33)
           dc_gaa = DIFF_Z_2(g11)
           dc_gab = DIFF_Z_2(g12)
           dc_gac = DIFF_Z_2(g13)
           dc_gbb = DIFF_Z_2(g22)
           dc_gbc = DIFF_Z_2(g23)
           dc_gcc = DIFF_Z_2(g33)
           
        else

           da_gaa = DIFF_X_4(g11)
           da_gab = DIFF_X_4(g12)
           da_gac = DIFF_X_4(g13)
           da_gbb = DIFF_X_4(g22)
           da_gbc = DIFF_X_4(g23)
           da_gcc = DIFF_X_4(g33)
           db_gaa = DIFF_Y_4(g11)
           db_gab = DIFF_Y_4(g12)
           db_gac = DIFF_Y_4(g13)
           db_gbb = DIFF_Y_4(g22)
           db_gbc = DIFF_Y_4(g23)
           db_gcc = DIFF_Y_4(g33)
           dc_gaa = DIFF_Z_4(g11)
           dc_gab = DIFF_Z_4(g12)
           dc_gac = DIFF_Z_4(g13)
           dc_gbb = DIFF_Z_4(g22)
           dc_gbc = DIFF_Z_4(g23)
           dc_gcc = DIFF_Z_4(g33)

        end if
          
!!$        Contract the shift with the eatrinsic curvature

        shiftshiftk = shifta*shifta*k11(i,j,k) + &
                      shiftb*shiftb*k22(i,j,k) + &
                      shiftc*shiftc*k33(i,j,k) + &
             two*(shifta*shiftb*k12(i,j,k) + &
                  shifta*shiftc*k13(i,j,k) + &
                  shiftb*shiftc*k23(i,j,k))

        shiftka = shifta*k11(i,j,k) + shiftb*k12(i,j,k) + shiftc*k13(i,j,k)
        shiftkb = shifta*k12(i,j,k) + shiftb*k22(i,j,k) + shiftc*k23(i,j,k)
        shiftkc = shifta*k13(i,j,k) + shiftb*k23(i,j,k) + shiftc*k33(i,j,k)

!!$        Contract the matter terms with the extrinsic curvature

        sumTK = taa*k11(i,j,k) + tbb*k22(i,j,k) + tcc*k33(i,j,k) &
             + two*(tab*k12(i,j,k) + tac*k13(i,j,k) + tbc*k23(i,j,k))

!!$        Update term for tau
        
        tau_source = t00* &
             (shiftshiftk - (shifta*da_alp + shiftb*db_alp + shiftc*dc_alp) )&
             + t0a*(-da_alp + two*shiftka) &
             + t0b*(-db_alp + two*shiftkb) &
             + t0c*(-dc_alp + two*shiftkc) &
             + sumTK

!!$        The following looks verb little like the terms in the
!!$        standard papers. Take a look in the ThornGuide to see why
!!$        it is really the same thing.

!!$        Contract the shift with derivatives of the metric

        halfshiftdga = half*(shifta*shifta*da_gaa + &
             shiftb*shiftb*da_gbb + shiftc*shiftc*da_gcc) + &
             shifta*shiftb*da_gab + shifta*shiftc*da_gac + &
             shiftb*shiftc*da_gbc
        halfshiftdgb = half*(shifta*shifta*db_gaa + &
             shiftb*shiftb*db_gbb + shiftc*shiftc*db_gcc) + &
             shifta*shiftb*db_gab + shifta*shiftc*db_gac + &
             shiftb*shiftc*db_gbc
        halfshiftdgc = half*(shifta*shifta*dc_gaa + &
             shiftb*shiftb*dc_gbb + shiftc*shiftc*dc_gcc) + &
             shifta*shiftb*dc_gab + shifta*shiftc*dc_gac + &
             shiftb*shiftc*dc_gbc

!!$        Contract the matter with derivatives of the metric

        halfTdga = half*(taa*da_gaa + tbb*da_gbb + tcc*da_gcc) +&
             tab*da_gab + tac*da_gac + tbc*da_gbc
        halfTdgb = half*(taa*db_gaa + tbb*db_gbb + tcc*db_gcc) +&
             tab*db_gab + tac*db_gac + tbc*db_gbc
        halfTdgc = half*(taa*dc_gaa + tbb*dc_gbb + tcc*dc_gcc) +&
             tab*dc_gab + tac*dc_gac + tbc*dc_gbc

        sa_source = t00*&
             (halfshiftdga - alp(i,j,k)*da_alp) +&
             t0a*(shifta*da_gaa + shiftb*da_gab + shiftc*da_gac) +&
             t0b*(shifta*da_gab + shiftb*da_gbb + shiftc*da_gbc) +&
             t0c*(shifta*da_gac + shiftb*da_gbc + shiftc*da_gcc) +&
             halfTdga + rhoenthalpyW2*&
             (vlowa*da_betaa + vlowb*da_betab + vlowc*da_betac)*&
             invalp
        sb_source = t00*&
             (halfshiftdgb - alp(i,j,k)*db_alp) +&
             t0a*(shifta*db_gaa + shiftb*db_gab + shiftc*db_gac) +&
             t0b*(shifta*db_gab + shiftb*db_gbb + shiftc*db_gbc) +&
             t0c*(shifta*db_gac + shiftb*db_gbc + shiftc*db_gcc) +&
             halfTdgb + rhoenthalpyW2*&
             (vlowa*db_betaa + vlowb*db_betab + vlowc*db_betac)*&
             invalp
        sc_source = t00*&
             (halfshiftdgc - alp(i,j,k)*dc_alp) +&
             t0a*(shifta*dc_gaa + shiftb*dc_gab + shiftc*dc_gac) +&
             t0b*(shifta*dc_gab + shiftb*dc_gbb + shiftc*dc_gbc) +&
             t0c*(shifta*dc_gac + shiftb*dc_gbc + shiftc*dc_gcc) +&
             halfTdgc + rhoenthalpyW2*&
             (vlowa*dc_betaa + vlowb*dc_betab + vlowc*dc_betac)*&
             invalp

        densrhs(i,j,k) = 0.d0
        srhs(i,j,k,1)  = alp(i,j,k)*sdetg(i,j,k)*sa_source
        srhs(i,j,k,2)  = alp(i,j,k)*sdetg(i,j,k)*sb_source
        srhs(i,j,k,3)  = alp(i,j,k)*sdetg(i,j,k)*sc_source
        taurhs(i,j,k) = alp(i,j,k)*sdetg(i,j,k)*tau_source
        
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(force_spatial_second_order)

  

end subroutine SourceTerms



