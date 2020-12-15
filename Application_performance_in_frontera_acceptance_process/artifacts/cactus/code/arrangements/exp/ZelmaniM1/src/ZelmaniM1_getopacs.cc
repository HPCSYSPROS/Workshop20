#include <iostream>
#include <algorithm>
#include <time.h>
#include <string.h>
#include <math.h>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "ZelmaniM1_Metric.hh"
#include "ZelmaniM1_Closure.hh"

#define SIGMA0 1.705e-44
#define C_L 2.99792458e10
#define G_GRAV 6.673e-8
#define M_SUN 1.98892e33
#define PI 3.14159265359

using namespace std;


namespace ZelmaniM1 {

  extern "C" 
  void zm1_setEquil(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const int nx = cctk_lsh[0];
    const int ny = cctk_lsh[1];
    const int nz = cctk_lsh[2];
    
    // loop over 3D grid
    #pragma omp parallel
    {
      // Allocate opac_tmp inside the parallel region, so that we
      // don't have to declare it as private. Variable-length arrays
      // are an extension to the C++ standard, and not all compilers
      // can handle private VLAs.
      const int nvars = ngroups*nspecies*3;
      double opac_tmp[nvars];
      #pragma omp for
      for(int k=0;k<nz;k++) 
        for(int j=0;j<ny;j++) 
          for(int i=0;i<nx;i++) {
            int i3D = CCTK_GFINDEX3D(cctkGH,i,j,k);
            ThreeTensor::Metric gamma(gxx[i3D],gxy[i3D],gxz[i3D],
                              gyy[i3D],gyz[i3D],gzz[i3D],
                              alp[i3D],betax[i3D],betay[i3D],betaz[i3D]);
                  
            double xrho  = log10(INV_RHO_GF * rho[i3D]);
            double xtemp = log10(temperature[i3D]);
                  
            double vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
            double vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
            double vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
            
            if (use_zerov) {
              vx = 0.0;
              vy = 0.0;
              vz = 0.0;
            }
            
            Closure cl(vx,vy,vz,gamma);
            double sdetg = sqrt(gamma.getDetg());
            double ye = Y_e[i3D];
            
            if (ye<GRHydro_Y_e_min) ye = GRHydro_Y_e_min;
            if (ye>GRHydro_Y_e_max) ye = GRHydro_Y_e_max;

            linterp_many(xrho, xtemp,ye,
                opac_tmp, alltables,
                nrho, ntemp, nye, nvars,
                rho_points, temp_points, ye_points);
            
             
            for (int ig=0; ig<ngroups*nspecies;ig++) { 
              int i4D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
              double ab = pow(10.0,opac_tmp[ig]);
              double em = pow(10.0,opac_tmp[ngroups*nspecies + ig])*sdetg;
              double sc = pow(10.0,opac_tmp[2*ngroups*nspecies + ig]);
          
              if (isnan(em))
                CCTK_VInfo(CCTK_THORNSTRING, 
                    "bad emis: %d %d %d %d %14.5e %14.5e %14.5e %14.5e",
                    i,j,k,ig,pow(10.0,xrho),pow(10.0,xtemp),Y_e[i3D],em);
              if (isnan(ab))
                CCTK_VInfo(CCTK_THORNSTRING, 
                    "bad abs : %d %d %d %d %14.5e %14.5e %14.5e %14.5e",
                    i,j,k,ig,pow(10.0,xrho),pow(10.0,xtemp),Y_e[i3D],ab);  
              double dtau = 1.e15;
              double E,N,Fx,Fy,Fz;
              double J,dJ,dN,dE,dSx,dSy,dSz;

              cl.setDebugInfo("ini");
              cl.setClosure(enu[i4D],fnux[i4D],fnuy[i4D],fnuz[i4D],
                  0.5*(xioldzL[i4D]+xioldzR[i4D]));
          
              cl.getSourceUpdate(0.0,0.0,0.0,0.0,0.0,ab,sc,em,
                  dtau,&E,&N,&Fx,&Fy,&Fz,&J,&dJ,&dN,&dE,&dSx,&dSy,&dSz);
    
              // Assume zero neutrino content at low density
              double free_stream_fac = tanh(rho[i3D]*INV_RHO_GF/1.e10); 
              enu[i4D]  = max(E*free_stream_fac,1.e-40);
              fnux[i4D] = Fx*free_stream_fac;
              fnuy[i4D] = Fy*free_stream_fac;
              fnuz[i4D] = Fz*free_stream_fac;
              enuh[i4D]  = enu[i4D];
              fnuxh[i4D] = fnux[i4D];
              fnuyh[i4D] = fnuy[i4D];
              fnuzh[i4D] = fnuz[i4D];
            }
          } // Energy and species loop
    } // Grid loop

    CCTK_Info(CCTK_THORNSTRING, "Radiation set to equilibrium");

    return;
  }

  extern "C" 
  void zm1_getopacs(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // We are treating the opacity GFs in a non-standard way.
    // This means that they will not be accessible from Fortran,
    // because we use a non-Cactus way of storing the data
    // the fastest varying index will be the last index,
    // running from 0 to ngroups*nspecies.

    const int nx = cctk_lsh[0];
    const int ny = cctk_lsh[1];
    const int nz = cctk_lsh[2];
    if (zm1_verbose) 
      CCTK_Info(CCTK_THORNSTRING, "Start interpolating opacs");
    // timer
    clock_t start, end;
    double elapsed;
    start = clock();
    // loop over 3D grid
#pragma omp parallel
    {
      // Allocate opac_tmp inside the parallel region, so that we
      // don't have to declare it as private. Variable-length arrays
      // are an extension to the C++ standard, and not all compilers
      // can handle private VLAs.
      const int nvars = ngroups*nspecies*3;
      double opac_tmp[nvars];
#pragma omp for
      for(int k=0;k<nz;k++) 
        for(int j=0;j<ny;j++) 
          for(int i=0;i<nx;i++) {
            int ind3D = CCTK_GFINDEX3D(cctkGH,i,j,k);
            double xrho = log10(INV_RHO_GF * rho[ind3D]);
            double xtemp = log10(temperature[ind3D]);
            
            linterp_many(xrho, xtemp, Y_e[ind3D],
                         opac_tmp, alltables,
                         nrho, ntemp, nye, nvars,
                         rho_points, temp_points, ye_points);
        
            for (int ig=0; ig<ngroups*nspecies;ig++) { 
              int ind4D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
              if (pow(10.0,xrho)>rho_opac_min) {
                absorb[ind4D] = pow(10.0,opac_tmp[ig]);
                emis[ind4D]   = pow(10.0,opac_tmp[ngroups*nspecies + ig]);
                scat[ind4D]   = pow(10.0,opac_tmp[2*ngroups*nspecies + ig]);
              } else {
                absorb[ind4D] = 0.0; 
                emis[ind4D]   = 0.0; 
                scat[ind4D]   = 0.0; 
              }   

#if 0 //For debugging
    if (absorb[ind4D] <= 0.0 || isnan(absorb[ind4D])) 
      CCTK_VInfo(CCTK_THORNSTRING, 
          "bad abs : %d %d %d %d %14.5e %14.5e %14.5e %14.5e",
          i,j,k,ig,pow(10.0,xrho),pow(10.0,xtemp),Y_e[ind3D],absorb[ind4D]);
    if (emis[ind4D] <= 0.0 || isnan(emis[ind4D]))     
      CCTK_VInfo(CCTK_THORNSTRING, 
          "bad emis: %d %d %d %d %14.5e %14.5e %14.5e %14.5e",
          i,j,k,ig,pow(10.0,xrho),pow(10.0,xtemp),Y_e[ind3D],emis[ind4D]);
    if (scat[ind4D] <= 0.0 || isnan(scat[ind4D]))     
      CCTK_VInfo(CCTK_THORNSTRING, 
          "bad scat: %d %d %d %d %14.5e %14.5e %14.5e %14.5e",
          i,j,k,ig,pow(10.0,xrho),pow(10.0,xtemp),Y_e[ind3D],scat[ind4D]);
#endif
            }

            //// now we need to distribute these into the opac arrays 
            //int ind4D = ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k);
            //// copy into absorb
            //memcpy( &absorb[ind4D], &opac_tmp[0], sizeof(double)*ngroups*nspecies);
            //// copy into emis
            //memcpy( &emis[ind4D], &opac_tmp[ngroups*nspecies], sizeof(double)*ngroups*nspecies);
            //// copy into scat
            //memcpy( &scat[ind4D], &opac_tmp[2*ngroups*nspecies], sizeof(double)*ngroups*nspecies);
      }
    }

    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    if (zm1_verbose) {
      CCTK_VInfo (CCTK_THORNSTRING,"time used in proc 0: %5.2f",elapsed);
      CCTK_Info(CCTK_THORNSTRING, "Done interpolating opacs");
    }

    return;
  }

  extern "C" 
  void zm1_getopacs_gray(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // We are treating the opacity GFs in a non-standard way.
    // This means that they will not be accessible from Fortran,
    // because we use a non-Cactus way of storing the data
    // the fastest varying index will be the last index,
    // running from 0 to ngroups*nspecies.

    const int nx = cctk_lsh[0];
    const int ny = cctk_lsh[1];
    const int nz = cctk_lsh[2];

    if (zm1_verbose) 
      CCTK_Info(CCTK_THORNSTRING, "Start interpolating opacs");
    // timer
    clock_t start, end;
    double elapsed;
    start = clock();
    // loop over 3D grid
#pragma omp parallel
    {
      // Allocate opac_tmp inside the parallel region, so that we
      // don't have to declare it as private. Variable-length arrays
      // are an extension to the C++ standard, and not all compilers
      // can handle private VLAs.
      const int nvars = ngtab*nspecies*3;
      double opac_tmp[nvars];
#pragma omp for
      for(int k=0;k<nz;k++) 
        for(int j=0;j<ny;j++) 
          for(int i=0;i<nx;i++) {
              int ind3D = CCTK_GFINDEX3D(cctkGH,i,j,k);
              double xrho = log10(INV_RHO_GF * rho[ind3D]);
              double xtemp = log10(temperature[ind3D]);
              
              linterp_many(xrho, xtemp, Y_e[ind3D],
                           opac_tmp, alltables,
                           nrho, ntemp, nye, nvars,
                           rho_points, temp_points, ye_points);

              for (int is=0; is<nspecies; is++) {
                // Calculate local planck opacities 
                double eta   = 0.0;
                double bbtot = 0.0;
                double ksbb  = 0.0;
                
                for (int ig=0;ig<ngtab;ig++)
                {
                  int idx = ig + ngtab*is;
                  double bb = opac_tmp[ngtab*nspecies + idx]/opac_tmp[idx];
                  bbtot  += bb;
                  eta    += opac_tmp[ngtab*nspecies + idx];
                  ksbb   += opac_tmp[2*ngtab*nspecies + idx]*bb;
                }

                for (int ig=0;ig<ngroups;ig++) 
                { 
                  int idx = ig + ngroups*is;
                        int ind4D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,idx);
                  absorb[ind4D] = eta/bbtot;
                  emis[ind4D]   = eta;
                  scat[ind4D]   = ksbb/bbtot;

                }
              }
      }
    }

    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    if (zm1_verbose) {
      CCTK_VInfo (CCTK_THORNSTRING,"time used in proc 0: %5.2f",elapsed);
      CCTK_Info(CCTK_THORNSTRING, "Done interpolating gray opacs");
    }

    return;
  }
}

