#include <cassert>
#include <cmath>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


// some prototypes
extern "C" CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);
static inline double SpatialDeterminantC(double gxx, double gxy,
                                        double gxz, double gyy,
                                        double gyz, double gzz);

static inline __attribute__((always_inline)) void prim2conC(double *w, double *dens, double *sx,
                             double *sy, double *sz, double *tau,
                             const double gxx, const double gxy,
                             const double gxz, const double gyy, const double gyz,
                             const double gzz, const double sdet,
                             const double rho, const double vx,
                             const double vy, const double vz,
                             const double eps, const double press);

static inline __attribute__((always_inline)) void prim2conMC(double *w, double *dens, double *sx,
                             double *sy, double *sz, double *tau,
			     double *Bconsx, double *Bconsy, double *Bconsz,
                             const double gxx, const double gxy,
                             const double gxz, const double gyy, const double gyz,
                             const double gzz, const double sdet,
                             const double rho, const double vx,
                             const double vy, const double vz,
                             const double eps, const double press,
			     double Bx, double By, double Bz);

static inline void GRHydro_SpeedOfSound(CCTK_ARGUMENTS);


extern "C" void GRHydro_Primitive2Conservative_CC(CCTK_ARGUMENTS);

extern "C"
CCTK_FCALL void CCTK_FNAME(GRHydro_Primitive2Conservative_CforF)(
  cGH * const * const p_cctkGH) {
  GRHydro_Primitive2Conservative_CC(*p_cctkGH);
}

extern "C"
CCTK_FCALL void CCTK_FNAME(GRHydro_SpeedOfSound)(cGH * const * const p_cctkGH) {
  GRHydro_SpeedOfSound(*p_cctkGH);
}

template<const bool do_mhd>
void GRHydro_Primitive2Conservative_CC_LL(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   // save memory when multipatch is not used
   CCTK_REAL * restrict g11, * restrict g12, * restrict g13;
   CCTK_REAL * restrict g22, * restrict g23, * restrict g33;

   if(GRHydro_UseGeneralCoordinates(cctkGH)) {
     g11 = gaa;
     g12 = gab;
     g13 = gac;
     g22 = gbb;
     g23 = gbc;
     g33 = gcc;
   } else {
     g11 = gxx;
     g12 = gxy;
     g13 = gxz;
     g22 = gyy;
     g23 = gyz;
     g33 = gzz;
   }

   GRHydro_SpeedOfSound(CCTK_PASS_CTOC);

   // if padding is used the the "extra" vector elements could contain junk
   // which we do not know how to handle
   assert(cctk_lsh[0] == cctk_ash[0]);
   assert(cctk_lsh[1] == cctk_ash[1]);
   assert(cctk_lsh[2] == cctk_ash[2]);

#pragma omp parallel for
     for(int k = GRHydro_stencil-1; k < cctk_lsh[2]-GRHydro_stencil+1; k++)
       for(int j = GRHydro_stencil-1; j < cctk_lsh[1]-GRHydro_stencil+1; j++)
#pragma ivdep // force compiler to vectorize the loop
         for(int i = GRHydro_stencil-1; i < cctk_lsh[0]-GRHydro_stencil+1; i++) {

           const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

           const int idxl = CCTK_GFINDEX3D(cctkGH,i-*xoffset,j-*yoffset,k-*zoffset);
           const int idxr = CCTK_GFINDEX3D(cctkGH,i+*xoffset,j+*yoffset,k+*zoffset);

           const double g11l = 0.5 * (g11[idx] + g11[idxl]);
           const double g12l = 0.5 * (g12[idx] + g12[idxl]);
           const double g13l = 0.5 * (g13[idx] + g13[idxl]);
           const double g22l = 0.5 * (g22[idx] + g22[idxl]);
           const double g23l = 0.5 * (g23[idx] + g23[idxl]);
           const double g33l = 0.5 * (g33[idx] + g33[idxl]);

           const double g11r = 0.5 * (g11[idx] + g11[idxr]);
           const double g12r = 0.5 * (g12[idx] + g12[idxr]);
           const double g13r = 0.5 * (g13[idx] + g13[idxr]);
           const double g22r = 0.5 * (g22[idx] + g22[idxr]);
           const double g23r = 0.5 * (g23[idx] + g23[idxr]);
           const double g33r = 0.5 * (g33[idx] + g33[idxr]);

           const double savg_detl =
             sqrt(SpatialDeterminantC(g11l,g12l,g13l,g22l,g23l,g33l));
           const double savg_detr =
             sqrt(SpatialDeterminantC(g11r,g12r,g13r,g22r,g23r,g33r));

           if (do_mhd) {
             // minus call to p2c
             prim2conMC(&w_lorentzminus[idx], &densminus[idx], &sxminus[idx],
                       &syminus[idx], &szminus[idx], &tauminus[idx],
                       &Bconsxminus[idx], &Bconsyminus[idx], &Bconszminus[idx],
                       g11l,g12l,g13l,g22l,g23l,g33l,
                       savg_detl,rhominus[idx], velxminus[idx], velyminus[idx],
                       velzminus[idx], epsminus[idx], pressminus[idx],
                       Bvecxminus[idx], Bvecyminus[idx], Bveczminus[idx]);


             // plus call to p2c
             prim2conMC(&w_lorentzplus[idx], &densplus[idx], &sxplus[idx],
                       &syplus[idx], &szplus[idx], &tauplus[idx],
                       &Bconsxplus[idx], &Bconsyplus[idx], &Bconszplus[idx],
                       g11r,g12r,g13r,g22r,g23r,g33r,
                       savg_detr,rhoplus[idx], velxplus[idx], velyplus[idx],
                       velzplus[idx], epsplus[idx], pressplus[idx],
                       Bvecxplus[idx], Bvecyplus[idx], Bveczplus[idx]);
           } else {
             // minus call to p2c
             prim2conC(&w_lorentzminus[idx], &densminus[idx], &sxminus[idx],
                       &syminus[idx], &szminus[idx], &tauminus[idx],
                       g11l,g12l,g13l,g22l,g23l,g33l,
                       savg_detl,rhominus[idx], velxminus[idx], velyminus[idx],
                       velzminus[idx], epsminus[idx], pressminus[idx]);


             // plus call to p2c
             prim2conC(&w_lorentzplus[idx], &densplus[idx], &sxplus[idx],
                       &syplus[idx], &szplus[idx], &tauplus[idx],
                       g11r,g12r,g13r,g22r,g23r,g33r,
                       savg_detr,rhoplus[idx], velxplus[idx], velyplus[idx],
                       velzplus[idx], epsplus[idx], pressplus[idx]);
	   }
         }

} // end function Conservative2PrimitiveC

static inline void GRHydro_SpeedOfSound(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_PARAMETERS;
   DECLARE_CCTK_ARGUMENTS;

   // if padding is used the the "extra" vector elements could contain junk
   // which we do not know how to handle
   assert(cctk_lsh[0] == cctk_ash[0]);
   assert(cctk_lsh[1] == cctk_ash[1]);
   assert(cctk_lsh[2] == cctk_ash[2]);

   // EOS calls (now GF-wide)
   if(!*evolve_temper) {
     // n needs to be computed using ash since ash is used when computing the
     // linear index in CCTK_GFINDEX3D
     const size_t n = size_t(cctk_ash[0]*cctk_ash[1]*cctk_ash[2]);
     std::vector<CCTK_INT> keyerr(n);
     CCTK_INT anyerr = 0;
     CCTK_INT keytemp = 0;

     // don't need special error handling for analytic EOS
#pragma omp parallel for
       for(int k=0;k<cctk_lsh[2];k++)
         for(int j=0;j<cctk_lsh[1];j++) {
           int i = CCTK_GFINDEX3D(cctkGH,0,j,k);
           EOS_Omni_press_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,cctk_lsh[0],
                          &(rhominus[i]),&(epsminus[i]),NULL,NULL,&(pressminus[i]),&(cs2minus[i]),
                          &(keyerr[i]),&anyerr);
           EOS_Omni_press_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,cctk_lsh[0],
                          &(rhoplus[i]),&(epsplus[i]),NULL,NULL,&(pressplus[i]),&(cs2plus[i]),
                          &(keyerr[i]),&anyerr);
         }
   } else {
     if(reconstruct_temper) {
       const size_t n = size_t(cctk_ash[0]*cctk_ash[1]*cctk_ash[2]);
       int nx = cctk_lsh[0] - (GRHydro_stencil - 1) - (GRHydro_stencil) + 1 + (transport_constraints && !*xoffset);
       std::vector<CCTK_INT> keyerr(n);
       int keytemp = 1;
       // ensure Y_e and temperature within bounds
#pragma omp parallel for
       for(int i=0;i<n;i++) { // walks over slightly too many elements but cannot fail
         Y_e_minus[i] = MAX(MIN(Y_e_minus[i],GRHydro_Y_e_max),GRHydro_Y_e_min);
         Y_e_plus[i] = MAX(MIN(Y_e_plus[i],GRHydro_Y_e_max),GRHydro_Y_e_min);
         tempminus[i] = MIN(MAX(tempminus[i],GRHydro_hot_atmo_temp),GRHydro_max_temp);
         tempplus[i] = MIN(MAX(tempplus[i],GRHydro_hot_atmo_temp),GRHydro_max_temp);
         keyerr[i] = 0;
        }

       // call the EOS with slices
#pragma omp parallel for
       for(int k=GRHydro_stencil-1;k<cctk_lsh[2]-GRHydro_stencil+1 + (transport_constraints && !*zoffset);k++)
         for(int j=GRHydro_stencil-1;j<cctk_lsh[1]-GRHydro_stencil+1 + (transport_constraints && !*yoffset);j++) {
           CCTK_INT anyerr = 0;
           int i = CCTK_GFINDEX3D(cctkGH,GRHydro_stencil-1,j,k);
           EOS_Omni_press_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,nx,
                          &rhominus[i],&epsminus[i],&tempminus[i],&Y_e_minus[i],
                          &pressminus[i],&cs2minus[i],
                          &keyerr[i],&anyerr);
           if(anyerr) {
             for(int ii=GRHydro_stencil-1;ii<cctk_lsh[0]-GRHydro_stencil+1;ii++) {
               int idx = CCTK_GFINDEX3D(cctkGH,ii,j,k);
               if(keyerr[idx]!=0) {
#pragma omp critical
                 {
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d i,j,k,x,y,z: %d,%d,%d %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), ii, j, k, x[idx], y[idx], z[idx], int(keyerr[idx]));
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d r,t,ye: %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), rhominus[idx], tempminus[idx], Y_e_minus[idx], int(keyerr[idx]));
                 }
               }
             }
             CCTK_ERROR("Aborting!");
           } // if (anyerr)
         } // loop


#pragma omp parallel for
       for(int k=GRHydro_stencil-1;k<cctk_lsh[2]-GRHydro_stencil+1 + (transport_constraints && !*zoffset);k++)
         for(int j=GRHydro_stencil-1;j<cctk_lsh[1]-GRHydro_stencil+1 + (transport_constraints && !*yoffset);j++) {
           CCTK_INT anyerr = 0;
           int i = CCTK_GFINDEX3D(cctkGH,GRHydro_stencil-1,j,k);
           EOS_Omni_press_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,nx,
                          &rhoplus[i],&epsplus[i],&tempplus[i],&Y_e_plus[i],
                          &pressplus[i],&cs2plus[i],
                          &keyerr[i],&anyerr);
           if(anyerr) {
             for(int ii=GRHydro_stencil-1;ii<cctk_lsh[0]-GRHydro_stencil+1;ii++) {
               int idx = CCTK_GFINDEX3D(cctkGH,ii,j,k);
               if(keyerr[idx]!=0) {
#pragma omp critical
                 {
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d i,j,k,x,y,z: %d,%d,%d %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), ii, j, k, x[idx], y[idx], z[idx], int(keyerr[idx]));
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d r,t,ye: %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), rhoplus[idx], tempplus[idx], Y_e_plus[idx], int(keyerr[idx]));
                 }
               }
             }
             CCTK_ERROR("Aborting!");
           } // if (anyerr)
         } // loop
     } else {
       // ******************** EPS RECONSTRUCTION BRANCH ******************
       int nx = cctk_lsh[0] - (GRHydro_stencil - 1) - (GRHydro_stencil) + 1 + (transport_constraints && !*xoffset);
       const int n = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
       std::vector<CCTK_INT> keyerr(n);
       int keytemp = 0;

      // ensure Y_e and temperature within bounds
#pragma omp parallel for
       for(int i=0;i<n;i++) {
         Y_e_minus[i] = MAX(MIN(Y_e_minus[i],GRHydro_Y_e_max),GRHydro_Y_e_min);
         Y_e_plus[i] = MAX(MIN(Y_e_plus[i],GRHydro_Y_e_max),GRHydro_Y_e_min);
         tempminus[i] = MIN(MAX(tempminus[i],GRHydro_hot_atmo_temp),GRHydro_max_temp);
         tempplus[i] = MIN(MAX(tempplus[i],GRHydro_hot_atmo_temp),GRHydro_max_temp);
         temperature[i] = MIN(MAX(temperature[i],GRHydro_hot_atmo_temp),GRHydro_max_temp);
         keyerr[i] = 0;
        }

       // call the EOS with slices for minus states
#pragma omp parallel for
       for(int k=GRHydro_stencil-1;k<cctk_lsh[2]-GRHydro_stencil+1 + (transport_constraints && !*zoffset);k++)
         for(int j=GRHydro_stencil-1;j<cctk_lsh[1]-GRHydro_stencil+1 + (transport_constraints && !*yoffset);j++) {
           int i = CCTK_GFINDEX3D(cctkGH,GRHydro_stencil-1,j,k);
           CCTK_INT anyerr = 0;
           EOS_Omni_press_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,nx,
                          &rhominus[i],&epsminus[i],&tempminus[i],&Y_e_minus[i],
                          &pressminus[i],&cs2minus[i],
                          &keyerr[i],&anyerr);
           if(anyerr) {
             for(int ii=GRHydro_stencil-1;ii<cctk_lsh[0]-GRHydro_stencil+1;ii++) {
               int idx = CCTK_GFINDEX3D(cctkGH,ii,j,k);
               if(keyerr[idx]!=0) {
#pragma omp critical
                 {
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d i,j,k,x,y,z: %d,%d,%d %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), ii, j, k, x[idx], y[idx], z[idx], int(keyerr[idx]));
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d r,t,ye: %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), rhominus[idx], tempminus[idx], Y_e_minus[idx], int(keyerr[idx]));

                   if(keyerr[idx] == 668) {
                     // This means the temperature came back negative.
                     // We'll try using piecewise constant for the temperature
                     tempminus[idx] = temperature[idx];
                     const int ln=1;
                     CCTK_INT lkeyerr[1];
                     CCTK_INT lanyerr = 0;
                     CCTK_INT lkeytemp = 1;
                     EOS_Omni_press_cs2(*GRHydro_eos_handle,lkeytemp,GRHydro_eos_rf_prec,ln,
                                        &rhominus[idx],&epsminus[idx],&tempminus[idx],
                                        &Y_e_minus[idx],&pressminus[idx],&cs2minus[idx],lkeyerr,&lanyerr);
                     if(lanyerr !=0) {
                       CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                                  "Aborting! keyerr=%d, r=%15.6E, t=%15.6E, ye=%15.6E",
                                  int(lkeyerr[0]),rhominus[idx],tempminus[idx],Y_e_minus[idx]);
                     }
                   } else {
                     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                                "Aborting! keyerr=%d, r=%15.6E, t=%15.6E, ye=%15.6E",
                                int(keyerr[idx]),rhominus[idx],tempminus[idx],Y_e_minus[idx]);
                   }
                 } // omp critical pragma
               } // if keyerr
             } // loop ii
           } // if (anyerr)
         } // big loop

       // call the EOS with slices for plus states
#pragma omp parallel for
       for(int k=GRHydro_stencil-1;k<cctk_lsh[2]-GRHydro_stencil+1 + (transport_constraints && !*zoffset);k++)
         for(int j=GRHydro_stencil-1;j<cctk_lsh[1]-GRHydro_stencil+1 + (transport_constraints && !*yoffset);j++) {
           int i = CCTK_GFINDEX3D(cctkGH,GRHydro_stencil-1,j,k);
           CCTK_INT anyerr = 0;
           EOS_Omni_press_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,nx,
                          &rhoplus[i],&epsplus[i],&tempplus[i],&Y_e_plus[i],
                          &pressplus[i],&cs2plus[i],
                          &keyerr[i],&anyerr);
           if(anyerr) {
             for(int ii=GRHydro_stencil-1;ii<cctk_lsh[0]-GRHydro_stencil+1;ii++) {
               int idx = CCTK_GFINDEX3D(cctkGH,ii,j,k);
               if(keyerr[idx]!=0) {
#pragma omp critical
                 {
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d i,j,k,x,y,z: %d,%d,%d %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), ii, j, k, x[idx], y[idx], z[idx], int(keyerr[idx]));
                   CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "rl: %d r,t,ye: %15.6E %15.6E %15.6E, keyerr: %d",
                              int(*GRHydro_reflevel), rhoplus[idx], tempplus[idx], Y_e_plus[idx], int(keyerr[idx]));

                   if(keyerr[idx] == 668) {
                     // This means the temperature came back negative.
                     // We'll try using piecewise constant for the temperature
                     tempplus[idx] = temperature[idx];
                     const int ln=1;
                     CCTK_INT lkeyerr[1];
                     CCTK_INT lanyerr = 0;
                     CCTK_INT lkeytemp = 1;
                     EOS_Omni_press_cs2(*GRHydro_eos_handle,lkeytemp,GRHydro_eos_rf_prec,ln,
                                        &rhoplus[idx],&epsplus[idx],&tempplus[idx],
                                        &Y_e_plus[idx],&pressplus[idx],&cs2plus[idx],lkeyerr,&lanyerr);
                     if(lanyerr !=0) {
                       CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                                  "Aborting! keyerr=%d, r=%15.6E, t=%15.6E, ye=%15.6E",
                                  int(lkeyerr[0]),rhoplus[idx],tempplus[idx],Y_e_plus[idx]);
                     }
                   } else {
                     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                                "Aborting! keyerr=%d, r=%15.6E, t=%15.6E, ye=%15.6E",
                                int(keyerr[idx]),rhoplus[idx],tempplus[idx],Y_e_plus[idx]);
                   }
                 } // omp critical pragma
               } // if keyerr
             } // loop ii
           } // if (anyerr)
         } // big loop over plus states

     } // end branch for no temp reconsturction
   } // end of evolve temper branch
}

static inline double SpatialDeterminantC(double gxx, double gxy,
                                         double gxz, double gyy,
                                         double gyz, double gzz)
{
  return -gxz*gxz*gyy + 2.0*gxy*gxz*gyz - gxx*gyz*gyz
    - gxy*gxy*gzz + gxx*gyy*gzz;
}


static inline __attribute__((always_inline)) void prim2conC(double *w, double *dens, double *sx,
                             double *sy, double *sz, double *tau,
                             const double gxx, const double gxy,
                             const double gxz, const double gyy,
                             const double gyz,
                             const double gzz, const double sdet,
                             const double rho, const double vx,
                             const double vy, const double vz,
                             const double eps, const double press)
{

  // local helpers
  const double wtemp = 1.0 / sqrt(1.0 - (gxx*vx*vx + gyy*vy*vy +
                                    gzz*vz*vz + 2.0*gxy*vx*vy
                                    + 2.0*gxz*vx*vz + 2.0*gyz*vy*vz));

  const double vlowx = gxx*vx + gxy*vy + gxz*vz;
  const double vlowy = gxy*vx + gyy*vy + gyz*vz;
  const double vlowz = gxz*vx + gyz*vy + gzz*vz;

  const double hrhow2 = (rho*(1.0+eps)+press)*(wtemp)*(wtemp);
  const double denstemp = sdet*rho*(wtemp);

  *w = wtemp;
  *dens = denstemp;
  *sx = sdet*hrhow2 * vlowx;
  *sy = sdet*hrhow2 * vlowy;
  *sz = sdet*hrhow2 * vlowz;
  *tau = sdet*( hrhow2 - press) - denstemp;

}

static inline __attribute__((always_inline)) void prim2conMC(double *w, double *dens, double *sx,
                             double *sy, double *sz, double *tau, 
			     double *Bconsx, double *Bconsy, double *Bconsz,
                             const double gxx, const double gxy,
                             const double gxz, const double gyy,
                             const double gyz,
                             const double gzz, const double sdet,
                             const double rho, const double vx,
                             const double vy, const double vz,
                             const double eps, const double press,
			     double Bx, double By, double Bz)
{

  // local helpers
  const double wtemp = 1.0 / sqrt(1.0 - (gxx*vx*vx + gyy*vy*vy +
                                    gzz*vz*vz + 2.0*gxy*vx*vy
                                    + 2.0*gxz*vx*vz + 2.0*gyz*vy*vz));

  const double vlowx = gxx*vx + gxy*vy + gxz*vz;
  const double vlowy = gxy*vx + gyy*vy + gyz*vz;
  const double vlowz = gxz*vx + gyz*vy + gzz*vz;

  const CCTK_REAL Bxlow = gxx*Bx + gxy*By + gxz*Bz;
  const CCTK_REAL Bylow = gxy*Bx + gyy*By + gyz*Bz;
  const CCTK_REAL Bzlow = gxz*Bx + gyz*By + gzz*Bz;
  
  const CCTK_REAL B2 =Bxlow*Bx+Bylow*By+Bzlow*Bz;

  const CCTK_REAL Bdotv = Bxlow*vx+Bylow*vy+Bzlow*vz;
  const CCTK_REAL Bdotv2 = Bdotv*Bdotv;
  const CCTK_REAL wtemp2 = wtemp*wtemp;
  const CCTK_REAL b2 = B2/wtemp2+Bdotv2;
  const CCTK_REAL ab0 = wtemp*Bdotv;
  const CCTK_REAL blowx = (gxx*Bx + gxy*By + gxz*Bz)/wtemp + wtemp*Bdotv*vlowx;
  const CCTK_REAL blowy = (gxy*Bx + gyy*By + gyz*Bz)/wtemp + wtemp*Bdotv*vlowy;
  const CCTK_REAL blowz = (gxz*Bx + gyz*By + gzz*Bz)/wtemp + wtemp*Bdotv*vlowz;

  const double hrhow2 = (rho*(1.0+eps)+press+b2)*(wtemp)*(wtemp);
  const double denstemp = sdet*rho*(wtemp);
  
  *w = wtemp;
  *dens = denstemp;
  
  *sx = sdet*hrhow2 * vlowx - ab0*blowx;
  *sy = sdet*hrhow2 * vlowy - ab0*blowy;
  *sz = sdet*hrhow2 * vlowz - ab0*blowz;

  *tau = sdet*( hrhow2 - press - b2/2.0 - ab0*ab0) - denstemp;

  *Bconsx = sdet*Bx;
  *Bconsy = sdet*By;
  *Bconsz = sdet*Bz;

}

extern "C" void GRHydro_Primitive2Conservative_CC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(*evolve_MHD)
    GRHydro_Primitive2Conservative_CC_LL<true>(CCTK_PASS_CTOC);
  else
    GRHydro_Primitive2Conservative_CC_LL<false>(CCTK_PASS_CTOC);
}
