 /*@@
   @file      GRHydro_SourceM.F90
   @date      Oct 11, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   The geometric source terms for the matter evolution
   @enddesc 
 @@*/

! Second order f.d.

#define DIFF_X_2(q) (0.5d0 * (q(i+1,j,k) - q(i-1,j,k)) * idx)
#define DIFF_Y_2(q) (0.5d0 * (q(i,j+1,k) - q(i,j-1,k)) * idy)
#define DIFF_Z_2(q) (0.5d0 * (q(i,j,k+1) - q(i,j,k-1)) * idz)

! Fourth order f.d.

#define DIFF_X_4(q) ((-q(i+2,j,k) + 8.d0 * q(i+1,j,k) - 8.d0 * q(i-1,j,k) + \
                      q(i-2,j,k)) / 12.d0 * idx)
#define DIFF_Y_4(q) ((-q(i,j+2,k) + 8.d0 * q(i,j+1,k) - 8.d0 * q(i,j-1,k) + \
                      q(i,j-2,k)) / 12.d0 * idy)
#define DIFF_Z_4(q) ((-q(i,j,k+2) + 8.d0 * q(i,j,k+1) - 8.d0 * q(i,j,k-1) + \
                      q(i,j,k-2)) / 12.d0 * idz)

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)
#define Bvecx(i,j,k) Bprim(i,j,k,1)
#define Bvecy(i,j,k) Bprim(i,j,k,2)
#define Bvecz(i,j,k) Bprim(i,j,k,3)
#define Avecx(i,j,k) Avec(i,j,k,1)
#define Avecy(i,j,k) Avec(i,j,k,2)
#define Avecz(i,j,k) Avec(i,j,k,3)
#define Avecrhsx(i,j,k) Avecrhs(i,j,k,1)
#define Avecrhsy(i,j,k) Avecrhs(i,j,k,2)
#define Avecrhsz(i,j,k) Avecrhs(i,j,k,3)

 /*@@
   @routine    SourceTermsAM
   @date       Aug 30, 2010
   @author     Tanja Bode, Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   Calculate the geometric source terms and add to the update GFs
   @enddesc 
   @calls     
   @calledby   
   @history 
   Minor alterations of routine from GR3D.
   @endhistory 
   
@@*/

subroutine SourceTermsAM(CCTK_ARGUMENTS)
      
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: one, two, half
  CCTK_REAL :: t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz
  CCTK_REAL :: sqrtdet, uxx, uxy, uxz, uyy, uyz, uzz
  CCTK_REAL :: shiftx, shifty, shiftz, velxshift, velyshift, velzshift 
  CCTK_REAL :: vlowx, vlowy, vlowz
  CCTK_REAL :: dx_betax, dx_betay, dx_betaz, dy_betax, dy_betay,&
       dy_betaz, dz_betax, dz_betay, dz_betaz
  CCTK_REAL :: dx_alp, dy_alp, dz_alp
  CCTK_REAL :: tau_source, sx_source, sy_source, sz_source
  CCTK_REAL :: localgxx,localgxy,localgxz,localgyy,localgyz,localgzz
  CCTK_REAL :: dx_gxx, dx_gxy, dx_gxz, dx_gyy, dx_gyz, dx_gzz
  CCTK_REAL :: dy_gxx, dy_gxy, dy_gxz, dy_gyy, dy_gyz, dy_gzz
  CCTK_REAL :: dz_gxx, dz_gxy, dz_gxz, dz_gyy, dz_gyz, dz_gzz
  CCTK_REAL :: dx, dy, dz, idx, idy, idz
  CCTK_REAL :: shiftshiftk, shiftkx, shiftky, shiftkz
  CCTK_REAL :: sumTK
  CCTK_REAL :: halfshiftdgx, halfshiftdgy, halfshiftdgz
  CCTK_REAL :: halfTdgx, halfTdgy, halfTdgz
  CCTK_REAL :: invalp, invalp2
  CCTK_REAL :: Avcx_source, Avcy_source, Avcz_source
  CCTK_REAL :: dx_det_bydet, dy_det_bydet, dz_det_bydet
  CCTK_REAL :: gdg_x, gdg_y, gdg_z !! g^{ik} d_k g_{ij}
  CCTK_INT  :: local_spatial_order

  CCTK_REAL :: Bvecxlow,Bvecylow,Bveczlow,bdotv,b2,dum1,dum2,bxlow,bylow,bzlow
  CCTK_REAL :: bt,bx,by,bz,rhohstarW2,pstar

  logical, allocatable, dimension (:,:,:) :: force_spatial_second_order

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: k11, k12, k13, k22, k23, k33
  pointer (pk11,k11), (pk12,k12), (pk13,k13), (pk22,k22), (pk23,k23), (pk33,k33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: beta1, beta2, beta3
  pointer (pbeta1,beta1), (pbeta2,beta2), (pbeta3,beta3)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)

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
    pBprim = loc(lBvec)
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
    pBprim = loc(Bvec)
  end if
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define kxx faulty_kxx
#define kxy faulty_kxy
#define kxz faulty_kxz
#define kyy faulty_kyy
#define kyz faulty_kyz
#define kzz faulty_kzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define gaa faulty_gaa
#define gab faulty_gab
#define gac faulty_gac
#define gbb faulty_gbb
#define gbc faulty_gbc
#define gcc faulty_gcc
#define kaa faulty_kaa
#define kab faulty_kab
#define kac faulty_kac
#define kbb faulty_kbb
#define kbc faulty_kbc
#define kcc faulty_kcc
#define betaa faulty_betaa
#define betab faulty_betab
#define betac faulty_betac
#define vel faulty_vel
#define Bvec faulty_Bvec

  one = 1.0d0
  two = 2.0d0
  half = 0.5d0
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  idx = 1.d0/dx
  idy = 1.d0/dy
  idz = 1.d0/dz
 
!!$  Initialize the update terms to be zero.
!!$  This will guarantee that no garbage in the boundaries is updated.

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0
  Avecrhs = 0.d0

  if (evolve_tracer .ne. 0) then
    cons_tracerrhs = 0.d0
  end if

  if (evolve_Y_e .ne. 0) then
     Y_e_con_rhs = 0.0d0
  endif
  
  if (clean_divergence .ne. 0) then
     psidcrhs=0.d0
  endif

  if (track_divB .ne. 0) then
     divB=0.d0
  endif

  if (transport_constraints .ne. 0) then
     Evec = 0.d0
  endif
  

!!$  Set up the array for checking the order. We switch to second order
!!$  differencing at boundaries and near excision regions.
!!$  Copied straight from BSSN.

  allocate (force_spatial_second_order(nx,ny,nz))
  force_spatial_second_order = .FALSE.
  
  if (sources_spatial_order > 2) then
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1 + GRHydro_stencil, nz - GRHydro_stencil
      do j = 1 + GRHydro_stencil, ny - GRHydro_stencil
        do i = 1 + GRHydro_stencil, nx - GRHydro_stencil
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
  !$OMP localgxx,localgxy,localgxz,localgyy,localgyz,localgzz,&
  !$OMP sqrtdet,shiftx,shifty,shiftz,&
  !$OMP dx_betax,dx_betay,dx_betaz,dy_betax,dy_betay,dy_betaz,&
  !$OMP dz_betax,dz_betay,dz_betaz,velxshift,velyshift,velzshift,&
  !$OMP vlowx,vlowy,vlowz,Bvecxlow,Bvecylow,Bveczlow, &
  !$OMP  bdotv,b2,dum1,dum2,bxlow,bylow,bzlow,bt,bx,by,bz,rhohstarW2,pstar,&
  !$OMP t00,t0x,t0y,t0z,txx,txy,txz,tyy,tyz,tzz,&
  !$OMP dx_alp,dy_alp,dz_alp,&
  !$OMP tau_source,sx_source,sy_source,sz_source,&
  !$OMP dx_det_bydet,dy_det_bydet,dz_det_bydet,&
  !$OMP gdg_x,gdg_y,gdg_z,&
  !$OMP Avcx_source,Avcy_source,Avcz_source,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz,&
  !$OMP dx_gxx, dx_gxy, dx_gxz, dx_gyy, dx_gyz, dx_gzz,&
  !$OMP dy_gxx, dy_gxy, dy_gxz, dy_gyy, dy_gyz, dy_gzz,&
  !$OMP dz_gxx, dz_gxy, dz_gxz, dz_gyy, dz_gyz, dz_gzz,&
  !$OMP shiftshiftk,shiftkx,shiftky,shiftkz,&
  !$OMP sumTK,halfshiftdgx,halfshiftdgy,halfshiftdgz,&
  !$OMP halfTdgx,halfTdgy,halfTdgz,invalp,invalp2)
  do k=1 + GRHydro_stencil,nz - GRHydro_stencil
    do j=1 + GRHydro_stencil,ny - GRHydro_stencil
      do i=1 + GRHydro_stencil,nx - GRHydro_stencil

        local_spatial_order = sources_spatial_order
        if (force_spatial_second_order(i,j,k)) then
          local_spatial_order = 2
        end if
        
!!$        Set the metric terms.
        
        localgxx = g11(i,j,k)
        localgxy = g12(i,j,k)
        localgxz = g13(i,j,k)
        localgyy = g22(i,j,k)
        localgyz = g23(i,j,k)
        localgzz = g33(i,j,k)

        sqrtdet = sdetg(i,j,k)
        call UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, sqrtdet*sqrtdet, localgxx,&
             localgxy, localgxz, localgyy, localgyz, localgzz)
        

        shiftx = beta1(i,j,k)
        shifty = beta2(i,j,k)
        shiftz = beta3(i,j,k)
        
        if (local_spatial_order .eq. 2) then
           
           dx_betax = DIFF_X_2(beta1)
           dx_betay = DIFF_X_2(beta2)
           dx_betaz = DIFF_X_2(beta3)
           
           dy_betax = DIFF_Y_2(beta1)
           dy_betay = DIFF_Y_2(beta2)
           dy_betaz = DIFF_Y_2(beta3)
            
           dz_betax = DIFF_Z_2(beta1)
           dz_betay = DIFF_Z_2(beta2)
           dz_betaz = DIFF_Z_2(beta3)

        else

           dx_betax = DIFF_X_4(beta1)
           dx_betay = DIFF_X_4(beta2)
           dx_betaz = DIFF_X_4(beta3)
           
           dy_betax = DIFF_Y_4(beta1)
           dy_betay = DIFF_Y_4(beta2)
           dy_betaz = DIFF_Y_4(beta3)
           
           dz_betax = DIFF_Z_4(beta1)
           dz_betay = DIFF_Z_4(beta2)
           dz_betaz = DIFF_Z_4(beta3)
           
        end if
          
        invalp = 1.0d0 / alp(i,j,k)
        invalp2 = invalp**2
        velxshift = velx(i,j,k) - shiftx*invalp
        velyshift = vely(i,j,k) - shifty*invalp
        velzshift = velz(i,j,k) - shiftz*invalp

        call calc_vlow_blow(localgxx,localgxy,localgxz,localgyy,localgyz,localgzz, &
             velx(i,j,k),vely(i,j,k),velz(i,j,k),Bvecx(i,j,k),Bvecy(i,j,k),Bvecz(i,j,k), &
             vlowx,vlowy,vlowz,Bvecxlow,Bvecylow,Bveczlow, &
             bdotv,b2,dum1,dum2,bxlow,bylow,bzlow)

!!$ These are the contravariant components
        bt = w_lorentz(i,j,k)/alp(i,j,k)*bdotv
        bx = Bvecx(i,j,k)/w_lorentz(i,j,k)+w_lorentz(i,j,k)*bdotv*velxshift
        by = Bvecy(i,j,k)/w_lorentz(i,j,k)+w_lorentz(i,j,k)*bdotv*velyshift
        bz = Bvecz(i,j,k)/w_lorentz(i,j,k)+w_lorentz(i,j,k)*bdotv*velzshift

        rhohstarW2 = (rho(i,j,k)*(one + eps(i,j,k)) + press(i,j,k)+ b2)*&
             w_lorentz(i,j,k)**2
        pstar = press(i,j,k)+0.5d0*b2

!!$        For a change, these are T^{ij}

        t00 = (rhohstarW2 - pstar)*invalp2-bt**2
        t0x = rhohstarW2*velxshift*invalp +&
             pstar*shiftx*invalp2-bt*bx
        t0y = rhohstarW2*velyshift*invalp +&
             pstar*shifty*invalp2-bt*by
        t0z = rhohstarW2*velzshift*invalp +&
             pstar*shiftz*invalp2-bt*bz
        txx = rhohstarW2*velxshift*velxshift +&
             pstar*(uxx - shiftx*shiftx*invalp2)-bx**2
        txy = rhohstarW2*velxshift*velyshift +&
             pstar*(uxy - shiftx*shifty*invalp2)-bx*by
        txz = rhohstarW2*velxshift*velzshift +&
             pstar*(uxz - shiftx*shiftz*invalp2)-bx*bz
        tyy = rhohstarW2*velyshift*velyshift +&
             pstar*(uyy - shifty*shifty*invalp2)-by**2
        tyz = rhohstarW2*velyshift*velzshift +&
             pstar*(uyz - shifty*shiftz*invalp2)-by*bz
        tzz = rhohstarW2*velzshift*velzshift +&
             pstar*(uzz - shiftz*shiftz*invalp2)-bz**2

!!$        Derivatives of the lapse, metric and shift

        if (local_spatial_order .eq. 2) then

          dx_alp = DIFF_X_2(alp)
          dy_alp = DIFF_Y_2(alp)
          dz_alp = DIFF_Z_2(alp)

        else

          dx_alp = DIFF_X_4(alp)
          dy_alp = DIFF_Y_4(alp)
          dz_alp = DIFF_Z_4(alp)

        end if
        
        if (local_spatial_order .eq. 2) then

           dx_gxx = DIFF_X_2(g11)
           dx_gxy = DIFF_X_2(g12)
           dx_gxz = DIFF_X_2(g13)
           dx_gyy = DIFF_X_2(g22)
           dx_gyz = DIFF_X_2(g23)
           dx_gzz = DIFF_X_2(g33)
           dy_gxx = DIFF_Y_2(g11)
           dy_gxy = DIFF_Y_2(g12)
           dy_gxz = DIFF_Y_2(g13)
           dy_gyy = DIFF_Y_2(g22)
           dy_gyz = DIFF_Y_2(g23)
           dy_gzz = DIFF_Y_2(g33)
           dz_gxx = DIFF_Z_2(g11)
           dz_gxy = DIFF_Z_2(g12)
           dz_gxz = DIFF_Z_2(g13)
           dz_gyy = DIFF_Z_2(g22)
           dz_gyz = DIFF_Z_2(g23)
           dz_gzz = DIFF_Z_2(g33)
           
        else

           dx_gxx = DIFF_X_4(g11)
           dx_gxy = DIFF_X_4(g12)
           dx_gxz = DIFF_X_4(g13)
           dx_gyy = DIFF_X_4(g22)
           dx_gyz = DIFF_X_4(g23)
           dx_gzz = DIFF_X_4(g33)
           dy_gxx = DIFF_Y_4(g11)
           dy_gxy = DIFF_Y_4(g12)
           dy_gxz = DIFF_Y_4(g13)
           dy_gyy = DIFF_Y_4(g22)
           dy_gyz = DIFF_Y_4(g23)
           dy_gzz = DIFF_Y_4(g33)
           dz_gxx = DIFF_Z_4(g11)
           dz_gxy = DIFF_Z_4(g12)
           dz_gxz = DIFF_Z_4(g13)
           dz_gyy = DIFF_Z_4(g22)
           dz_gyz = DIFF_Z_4(g23)
           dz_gzz = DIFF_Z_4(g33)

        end if
          
!!$        Contract the shift with the extrinsic curvature

        shiftshiftk = shiftx*shiftx*k11(i,j,k) + &
                      shifty*shifty*k22(i,j,k) + &
                      shiftz*shiftz*k33(i,j,k) + &
             two*(shiftx*shifty*k12(i,j,k) + &
                  shiftx*shiftz*k13(i,j,k) + &
                  shifty*shiftz*k23(i,j,k))

        shiftkx = shiftx*k11(i,j,k) + shifty*k12(i,j,k) + shiftz*k13(i,j,k)
        shiftky = shiftx*k12(i,j,k) + shifty*k22(i,j,k) + shiftz*k23(i,j,k)
        shiftkz = shiftx*k13(i,j,k) + shifty*k23(i,j,k) + shiftz*k33(i,j,k)

!!$        Contract the matter terms with the extrinsic curvature

        sumTK = txx*k11(i,j,k) + tyy*k22(i,j,k) + tzz*k33(i,j,k) &
             + two*(txy*k12(i,j,k) + txz*k13(i,j,k) + tyz*k23(i,j,k))

!!$        Update term for tau
        
        tau_source = t00* &
             (shiftshiftk - (shiftx*dx_alp + shifty*dy_alp + shiftz*dz_alp) )&
             + t0x*(-dx_alp + two*shiftkx) &
             + t0y*(-dy_alp + two*shiftky) &
             + t0z*(-dz_alp + two*shiftkz) &
             + sumTK

!!$        The following looks very little like the terms in the
!!$        standard papers. Take a look in the ThornGuide to see why
!!$        it is really the same thing.

!!$        Contract the shift with derivatives of the metric

        halfshiftdgx = half*(shiftx*shiftx*dx_gxx + &
             shifty*shifty*dx_gyy + shiftz*shiftz*dx_gzz) + &
             shiftx*shifty*dx_gxy + shiftx*shiftz*dx_gxz + &
             shifty*shiftz*dx_gyz
        halfshiftdgy = half*(shiftx*shiftx*dy_gxx + &
             shifty*shifty*dy_gyy + shiftz*shiftz*dy_gzz) + &
             shiftx*shifty*dy_gxy + shiftx*shiftz*dy_gxz + &
             shifty*shiftz*dy_gyz
        halfshiftdgz = half*(shiftx*shiftx*dz_gxx + &
             shifty*shifty*dz_gyy + shiftz*shiftz*dz_gzz) + &
             shiftx*shifty*dz_gxy + shiftx*shiftz*dz_gxz + &
             shifty*shiftz*dz_gyz

!!$        Contract the matter with derivatives of the metric

        halfTdgx = half*(txx*dx_gxx + tyy*dx_gyy + tzz*dx_gzz) +&
             txy*dx_gxy + txz*dx_gxz + tyz*dx_gyz
        halfTdgy = half*(txx*dy_gxx + tyy*dy_gyy + tzz*dy_gzz) +&
             txy*dy_gxy + txz*dy_gxz + tyz*dy_gyz
        halfTdgz = half*(txx*dz_gxx + tyy*dz_gyy + tzz*dz_gzz) +&
             txy*dz_gxy + txz*dz_gxz + tyz*dz_gyz

     


       sx_source = t00*&
             (halfshiftdgx - alp(i,j,k)*dx_alp) + halfTdgx + &
             t0x*(shiftx*dx_gxx + shifty*dx_gxy + shiftz*dx_gxz) +&
             t0y*(shiftx*dx_gxy + shifty*dx_gyy + shiftz*dx_gyz) +&
             t0z*(shiftx*dx_gxz + shifty*dx_gyz + shiftz*dx_gzz) +&
             rhohstarW2*invalp*(vlowx*dx_betax + vlowy*dx_betay + vlowz*dx_betaz) -&
             bt*(bxlow*dx_betax + bylow*dx_betay + bzlow*dx_betaz)
        
       sy_source = t00*&
             (halfshiftdgy - alp(i,j,k)*dy_alp) + halfTdgy + &
             t0x*(shiftx*dy_gxx + shifty*dy_gxy + shiftz*dy_gxz) +&
             t0y*(shiftx*dy_gxy + shifty*dy_gyy + shiftz*dy_gyz) +&
             t0z*(shiftx*dy_gxz + shifty*dy_gyz + shiftz*dy_gzz) +&
             rhohstarW2*invalp*(vlowx*dy_betax + vlowy*dy_betay + vlowz*dy_betaz) -&
             bt*(bxlow*dy_betax + bylow*dy_betay + bzlow*dy_betaz)

       sz_source = t00*&
             (halfshiftdgz - alp(i,j,k)*dz_alp) + halfTdgz + &
             t0x*(shiftx*dz_gxx + shifty*dz_gxy + shiftz*dz_gxz) +&
             t0y*(shiftx*dz_gxy + shifty*dz_gyy + shiftz*dz_gyz) +&
             t0z*(shiftx*dz_gxz + shifty*dz_gyz + shiftz*dz_gzz) +&
             rhohstarW2*invalp*(vlowx*dz_betax + vlowy*dz_betay + vlowz*dz_betaz) -&
             bt*(bxlow*dz_betax + bylow*dz_betay + bzlow*dz_betaz)

        !! B^i and A^i both live in cell centers currently
        Avcx_source = sqrtdet*(vely(i,j,k)*Bvecz(i,j,k) - velz(i,j,k)*Bvecy(i,j,k))
        Avcy_source = sqrtdet*(velz(i,j,k)*Bvecx(i,j,k) - velx(i,j,k)*Bvecz(i,j,k))
        Avcz_source = sqrtdet*(velx(i,j,k)*Bvecy(i,j,k) - vely(i,j,k)*Bvecx(i,j,k))

        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = 0.d0 
        end if

        densrhs(i,j,k) = 0.d0
        srhs(i,j,k,1)  = alp(i,j,k)*sqrtdet*sx_source
        srhs(i,j,k,2)  = alp(i,j,k)*sqrtdet*sy_source
        srhs(i,j,k,3)  = alp(i,j,k)*sqrtdet*sz_source
        taurhs(i,j,k)  = alp(i,j,k)*sqrtdet*tau_source
        Avecrhsx(i,j,k) = Avcx_source
        Avecrhsy(i,j,k) = Avcy_source
        Avecrhsz(i,j,k) = Avcz_source

      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(force_spatial_second_order)

#if(0) /* poison edges of domain */
  if(last_iteration_seen .ne. cctk_iteration .or. reflevel .ne. grhydro_reflevel) then
    last_iteration_seen = cctk_iteration
    reflevel = grhydro_reflevel
    mol_substep = 0
  else
    mol_substep = mol_substep + 1
  end if
  do k = 1, GRHydro_stencil*mol_substep
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        dens(i,j,k) = -1d100
        Scon(i,j,k,1) = -1d100
        Scon(i,j,k,2) = -1d100
        Scon(i,j,k,3) = -1d100
        tau(i,j,k) = -1d100
        Avecx(i,j,k) = -1d100
        Avecy(i,j,k) = -1d100
        Avecz(i,j,k) = -1d100
        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = -1d100
        end if
      end do
    end do
  end do
  do k = cctk_lsh(3)-GRHydro_stencil*mol_substep+1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        dens(i,j,k) = -1d100
        Scon(i,j,k,1) = -1d100
        Scon(i,j,k,2) = -1d100
        Scon(i,j,k,3) = -1d100
        tau(i,j,k) = -1d100
        Avecx(i,j,k) = -1d100
        Avecy(i,j,k) = -1d100
        Avecz(i,j,k) = -1d100
        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = -1d100
        end if
      end do
    end do
  end do
  do i = 1, GRHydro_stencil*mol_substep
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        dens(i,j,k) = -1d100
        Scon(i,j,k,1) = -1d100
        Scon(i,j,k,2) = -1d100
        Scon(i,j,k,3) = -1d100
        tau(i,j,k) = -1d100
        Avecx(i,j,k) = -1d100
        Avecy(i,j,k) = -1d100
        Avecz(i,j,k) = -1d100
        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = -1d100
        end if
      end do
    end do
  end do
  do i = cctk_lsh(1)-GRHydro_stencil*mol_substep+1, cctk_lsh(1)
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        dens(i,j,k) = -1d100
        Scon(i,j,k,1) = -1d100
        Scon(i,j,k,2) = -1d100
        Scon(i,j,k,3) = -1d100
        tau(i,j,k) = -1d100
        Avecx(i,j,k) = -1d100
        Avecy(i,j,k) = -1d100
        Avecz(i,j,k) = -1d100
        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = -1d100
        end if
      end do
    end do
  end do
  do j = 1, GRHydro_stencil*mol_substep
    do i = 1, cctk_lsh(1)
      do k = 1, cctk_lsh(3)
        dens(i,j,k) = -1d100
        Scon(i,j,k,1) = -1d100
        Scon(i,j,k,2) = -1d100
        Scon(i,j,k,3) = -1d100
        tau(i,j,k) = -1d100
        Avecx(i,j,k) = -1d100
        Avecy(i,j,k) = -1d100
        Avecz(i,j,k) = -1d100
        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = -1d100
        end if
      end do
    end do
  end do
  do j = cctk_lsh(2)-GRHydro_stencil*mol_substep+1, cctk_lsh(2)
    do i = 1, cctk_lsh(1)
      do k = 1, cctk_lsh(3)
        dens(i,j,k) = -1d100
        Scon(i,j,k,1) = -1d100
        Scon(i,j,k,2) = -1d100
        Scon(i,j,k,3) = -1d100
        tau(i,j,k) = -1d100
        Avecx(i,j,k) = -1d100
        Avecy(i,j,k) = -1d100
        Avecz(i,j,k) = -1d100
        if ( evolve_Lorenz_gge.gt.0 ) then
           Aphi(i,j,k) = -1d100
        end if
      end do
    end do
  end do
#endif

#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_vel
#undef faulty_Bvec
#undef faulty_gxx_p
#undef faulty_gxy_p
#undef faulty_gxz_p
#undef faulty_gyy_p
#undef faulty_gyz_p
#undef faulty_gzz_p
#undef faulty_vel_p
#undef faulty_Bvec_p
#undef faulty_gxx_p_p
#undef faulty_gxy_p_p
#undef faulty_gxz_p_p
#undef faulty_gyy_p_p
#undef faulty_gyz_p_p
#undef faulty_gzz_p_p
#undef faulty_vel_p_p
#undef faulty_Bvec_p_p
end subroutine SourceTermsAM



