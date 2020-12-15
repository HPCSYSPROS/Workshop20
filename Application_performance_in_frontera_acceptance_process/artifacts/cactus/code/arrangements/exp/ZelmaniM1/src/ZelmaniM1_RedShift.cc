#include <cassert>
#include <algorithm>
#include <iostream>
#include <time.h>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "ZelmaniM1_Metric.hh"
#include "ZelmaniM1_Closure.hh"
#include "carpet.hh"

#define INV_RHO_GF 6.18735016707159e17
#define INV_EPS_GF 8.98755175549085e20
#define EPS_GF 1.11265006e-21
#define C_L 2.99792458e10
#define G_GRAV 6.673e-8
#define M_SUN 1.98892e33
#define PI 3.14159265359
#define YE_CONV_FAC 931.4944016646

using namespace std;


namespace ZelmaniM1 {

  extern "C" 
  void zm1_Redshift(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];

    const double dt = CCTK_DELTA_TIME;
    const double idx = 1.0/CCTK_DELTA_SPACE(0); // 1/delta x of zone  
    const double idy = 1.0/CCTK_DELTA_SPACE(1); // 1/delta x of zone  
    const double idz = 1.0/CCTK_DELTA_SPACE(2); // 1/delta x of zone  
    
    double *enusr, *fnuxsr, *fnuysr, *fnuzsr;
    
    const int flux_ndirs = 3;   // x y z
    const int flux_nelts = ngroups*nspecies;
    const int flux_nvars = 4;   // enu fnux fnuy fnuz
    
    // These pointers should never actually be accessed if do_m1_reflux = .false. 
    const double * flux;
    double * register_fine;
    double * register_coarse;
    if (do_m1_reflux){  
      //const double * flux =
      flux =
        (const CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "Refluxing::flux[0]");
      //double * register_fine =
      register_fine = 
        (CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "Refluxing::register_fine[0]");
      //double * register_coarse =
      register_coarse = 
        (CCTK_REAL*)CCTK_VarDataPtr(cctkGH, 0, "Refluxing::register_coarse[0]");
      assert(flux);
      assert(register_fine);
      assert(register_coarse);
    }
    
     
    if (zm1_verbose)
      CCTK_VInfo(CCTK_THORNSTRING,"Redhshift : %d "
        "  cctk_time: %8.5f",*zm1_RKstep,cctk_time);
    
    if (*zm1_RKstep == 2) {
      enusr  = enu;
      fnuxsr = fnux;
      fnuysr = fnuy;
      fnuzsr = fnuz;
#pragma omp parallel for
      // This is currently implemented in an exceptionally stupid way
      for(int ig=0;ig<ngroups*nspecies;ig++)
        for(int k=0;k<nz;k++) 
          for(int j=0;j<ny;j++)
            for(int i=0;i<nx;i++) {
              // Shouldn't have to do this if we are setting boundary conditions
              //int index4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
              //enuh[index4D]  = enu[index4D];
              //fnuxh[index4D] = fnux[index4D];
              //fnuyh[index4D] = fnuy[index4D];
              //fnuzh[index4D] = fnuz[index4D];
              if (do_m1_reflux) {
                for (int var=0; var<flux_nvars; ++var) {
                  for (int dir=0; dir<flux_ndirs; ++dir) {
                    const int n = (var * flux_nelts + ig) * flux_ndirs + dir;
                    const int flux_idx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, n);
                    register_fine_h[flux_idx] = register_fine[flux_idx];
                    register_coarse_h[flux_idx] = register_coarse[flux_idx];
                  }
                }
              }
      } // End grid loop
    } else { 
      if (do_m1_RK2) {
        enusr  = enuh;
        fnuxsr = fnuxh;
        fnuysr = fnuyh;
        fnuzsr = fnuzh;
      } else {
        enusr  = enu;
        fnuxsr = fnux;
        fnuysr = fnuy;
        fnuzsr = fnuz;
      }
    }

    #pragma omp parallel
    {
    #pragma omp for 
    for(int k=zm1_ghost;k<nz-zm1_ghost;k++) 
      for(int j=zm1_ghost;j<ny-zm1_ghost;j++)
        for(int i=zm1_ghost;i<nx-zm1_ghost;i++) {
          
          int i3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
          int i3Dxp = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
          int i3Dxm = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
          int i3Dyp = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
          int i3Dym = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
          int i3Dzp = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
          int i3Dzm = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
         
          ThreeTensor::Metric gammac(gxx[i3D],gxy[i3D],gxz[i3D],
              gyy[i3D],gyz[i3D],gzz[i3D]);
          ThreeTensor::Metric gammaco(gxx_p[i3D],gxy_p[i3D],gxz_p[i3D],
              gyy_p[i3D],gyz_p[i3D],gzz_p[i3D]);
          gammac.loadExtCurv(kxx[i3D],kxy[i3D],kxz[i3D],
              kyy[i3D],kyz[i3D],kzz[i3D]);

          // Quantities required for velocity depenedent red-shifting       
          double alph   = alp[i3D];
          double dlnadx = (log(alp[i3Dxp]/alp[i3Dxm]))*idx*0.5;
          double dlnady = (log(alp[i3Dyp]/alp[i3Dym]))*idy*0.5;
          double dlnadz = (log(alp[i3Dzp]/alp[i3Dzm]))*idz*0.5;
           
          Closure cl(0.0,0.0,0.0,gammac);
          cl.setDebugInfo("red");
          
          for(int isp=0;isp<nspecies;isp++) { 
            // Calculate momentum space fluxes
            // this loop modifies enurhs at ieg and ieg-1 and ieg+1
            // this needs to be taken into account if OpenMP parallelizing it
            for(int ieg=0;ieg<ngroups;ieg++) {
              int ig = isp*ngroups + ieg;
              int index4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
              int igl = isp*ngroups + ieg-1;
              int index4Dl = ieg>0 ? CCTK_VectGFIndex3D(cctkGH,i,j,k,igl) : -1;
              int igu = isp*ngroups + ieg+1;
              int index4Du = ieg<ngroups-1 ? 
                  CCTK_VectGFIndex3D(cctkGH,i,j,k,igu) : -1;
              
              // Need un-energy integrated quantities, divide through by group width 
              double E  = enusr[index4D]/bin_widths[ieg];
              double Fx = fnuxsr[index4D]/bin_widths[ieg]; 
              double Fy = fnuysr[index4D]/bin_widths[ieg]; 
              double Fz = fnuzsr[index4D]/bin_widths[ieg];
    
              cl.setClosure(E,Fx,Fy,Fz,0.5*(xioldzL[index4D]+xioldzR[index4D]));
          
              double fE,fN,fFx,fFy,fFz;
              cl.getMomentumFluxes(dlnadx,dlnady,dlnadz,&fE,&fN,&fFx,&fFy,&fFz); 
                    
              // Calculate the fluxes using upwinding
              if (fE>0 && ieg<ngroups-1) {
                enurhs[index4D]  -= alph*fE*bin_top[ieg];
                enurhs[index4Du] += alph*fE*bin_top[ieg];
              } else if (ieg>0) {
                enurhs[index4D]  += alph*fE*bin_bottom[ieg];
                enurhs[index4Dl] -= alph*fE*bin_bottom[ieg];
              }
              if (fFx*Fx>0 && ieg<ngroups-1) {
                fnuxrhs[index4D]  -= alph*fFx*bin_top[ieg];
                fnuxrhs[index4Du] += alph*fFx*bin_top[ieg];
              } else if (ieg>0) {
                fnuxrhs[index4D]  += alph*fFx*bin_bottom[ieg];
                fnuxrhs[index4Dl] -= alph*fFx*bin_bottom[ieg];
              }
              if (fFy*Fy>0 && ieg<ngroups-1) {
                fnuyrhs[index4D]  -= alph*fFy*bin_top[ieg];
                fnuyrhs[index4Du] += alph*fFy*bin_top[ieg];
              } else if (ieg>0) {
                fnuyrhs[index4D]  += alph*fFy*bin_bottom[ieg];
                fnuyrhs[index4Dl] -= alph*fFy*bin_bottom[ieg];
              }
              if (fFz*Fz>0 && ieg<ngroups-1) {
                fnuzrhs[index4D]  -= alph*fFz*bin_top[ieg];
                fnuzrhs[index4Du] += alph*fFz*bin_top[ieg];
              } else if (ieg>0) {
                fnuzrhs[index4D]  += alph*fFz*bin_bottom[ieg];
                fnuzrhs[index4Dl] -= alph*fFz*bin_bottom[ieg];
              }
            } // Loop over energies
          } // Loop over species
       } // Loop over zones
    } // OpenMP parallel region
  }
}

