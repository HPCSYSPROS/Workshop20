#include <algorithm>
#include <cassert>
#include <iostream>
#include <time.h>
#include "ZelmaniM1.hh"
//#include "ZelmaniM1_AsymptoticF.hh"
#include "ZelmaniM1_Closure.hh"
#include "ZelmaniM1_Metric.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "carpet.hh"

#define C_L 2.99792458e10
#define G_GRAV 6.673e-8
#define M_SUN 1.98892e33

using namespace std;


namespace ZelmaniM1 {

extern "C" 
void zm1_RiemannHLLE(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];
  
  if (zm1_verbose) 
    CCTK_VInfo(CCTK_THORNSTRING,"RiemannHLLE: %d %d cctk_time: %8.5f",
      *zm1_RKstep,*zm1_flux_direction,cctk_time);

  //    if (Carpet::reflevel <= 2) return;

  int direction;
  if      (*zm1_xoffset==1) direction = 0;  
  else if (*zm1_yoffset==1) direction = 1;  
  else if (*zm1_zoffset==1) direction = 2;  
  else direction = -1;
  
  double idx = 1.0;
  //double idy = 1.0;
  //double idz = 1.0;

  if (direction == 0) {
    idx = 1.0/CCTK_DELTA_SPACE(0);
    //idy = 1.0/CCTK_DELTA_SPACE(1);
    //idz = 1.0/CCTK_DELTA_SPACE(2);
  } else if (direction == 1) {
    idx = 1.0/CCTK_DELTA_SPACE(1);
    //idy = 1.0/CCTK_DELTA_SPACE(2);
    //idz = 1.0/CCTK_DELTA_SPACE(0);
  } else if (direction == 2) { 
    idx = 1.0/CCTK_DELTA_SPACE(2);
    //idy = 1.0/CCTK_DELTA_SPACE(0);
    //idz = 1.0/CCTK_DELTA_SPACE(1);
  }

  // Point at the correct arrays for the current RK step
  double *enupt;
  double *fnuxpt;
  double *fnuypt;
  double *fnuzpt;
  if (do_m1_RK2 && (*zm1_RKstep == 1)) {
    enupt  = enuh;
    fnuxpt = fnuxh;
    fnuypt = fnuyh;
    fnuzpt = fnuzh;
  
  } else {
    enupt  = enu;
    fnuxpt = fnux;
    fnuypt = fnuy;
    fnuzpt = fnuz;
  }
  
  // Refluxing: capture the fluxes
  const int flux_ndirs = 3;   // x y z
  const int flux_nelts = ngroups*nspecies;
  //const int flux_nvars = 4;   // enu fnux fnuy fnuz
  int flux_var_offset;

  CCTK_REAL * restrict flux;
  if (do_m1_reflux) {
    const int flux_idx = CCTK_VarIndex("Refluxing::flux[0]");
    //CCTK_REAL *restrict const flux =
    flux = 
      static_cast<CCTK_REAL*>(CCTK_VarDataPtrI(cctkGH, 0, flux_idx));
    assert(flux);
    // Indexing order: flux[var][elt][dir]
    //const int flux_var_offset = 
    flux_var_offset = 
      CCTK_VECTGFINDEX3D(cctkGH, 0,0,0, flux_nelts * flux_ndirs) -
      CCTK_VECTGFINDEX3D(cctkGH, 0,0,0, 0);
  }
  
  #pragma omp parallel shared(idx,direction)
  {     
  #pragma omp for
  // Iteration over zone edges
  // Start on the last ghost zone and end on the last real zone
  for(int k=zm1_ghost-1;k<nz-zm1_ghost;k++) 
    for(int j=zm1_ghost-1;j<ny-zm1_ghost;j++)
      for(int i=zm1_ghost-1;i<nx-zm1_ghost;i++) {
        
        int index3D  = CCTK_GFINDEX3D(cctkGH,i,j,k);
        int index3Dm = CCTK_GFINDEX3D(cctkGH,i + *zm1_xoffset,
                          j + *zm1_yoffset,k + *zm1_zoffset);
        
        // Calculate some stuff at the zone edge
        const double xedge = 0.5*(x[index3D] + x[index3Dm]);
        const double yedge = 0.5*(y[index3D] + y[index3Dm]);
        const double zedge = 0.5*(z[index3D] + z[index3Dm]);
        const double redge = sqrt(xedge*xedge + yedge*yedge + zedge*zedge);
        const double rhx = xedge/redge; 
        const double rhy = yedge/redge; 
        const double rhz = zedge/redge; 
        
        // Closure factor for fixed closure 
        const double xi_guess = 0.5 + 0.5*tanh(((redge 
            - closure_transition_radius)/closure_transition_scale));
         
        // Interpolate lapse, shift, and spatial metric to both 
        // edges of zone (i,j,k) in given offset direction
        // beta has a raised index 
        // gamma has two lowered indices
        ThreeTensor::Metric gammaLc(gxx[index3D],gxy[index3D],gxz[index3D],
                            gyy[index3D],gyz[index3D],gzz[index3D],
                            alp[index3D],betax[index3D],
                            betay[index3D],betaz[index3D]);
        ThreeTensor::Metric gammaRc(gxx[index3Dm],gxy[index3Dm],gxz[index3Dm],
                            gyy[index3Dm],gyz[index3Dm],gzz[index3Dm],
                            alp[index3Dm],betax[index3Dm],
                            betay[index3Dm],betaz[index3Dm]);
        ThreeTensor::Metric gammac(0.5 * (gxx[index3Dm]   + gxx[index3D]),
                                   0.5 * (gxy[index3Dm]   + gxy[index3D]),
                                   0.5 * (gxz[index3Dm]   + gxz[index3D]),
                                   0.5 * (gyy[index3Dm]   + gyy[index3D]),
                                   0.5 * (gyz[index3Dm]   + gyz[index3D]),
                                   0.5 * (gzz[index3Dm]   + gzz[index3D]),
                                   0.5 * (alp[index3Dm]   + alp[index3D]),
                                   0.5 * (betax[index3Dm] + betax[index3D]),
                                   0.5 * (betay[index3Dm] + betay[index3D]),
                                   0.5 * (betaz[index3Dm] + betaz[index3D]));
 
        // This loading needs to depend on the direction we are solving in
        double vxL,vyL,vzL;
        double vxR,vyR,vzR;
        
        vxL = velxp[index3D];
        vxR = velxm[index3Dm];
        vyL = velyp[index3D];
        vyR = velym[index3Dm];
        vzL = velzp[index3D];
        vzR = velzm[index3Dm];
        
        vxL = velp[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]; 
        vyL = velp[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]; 
        vzL = velp[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]; 
        vxR = velm[CCTK_VECTGFINDEX3D(cctkGH,i + *zm1_xoffset,
                                             j + *zm1_yoffset,
                                             k + *zm1_zoffset,0)];
        vyR = velm[CCTK_VECTGFINDEX3D(cctkGH,i + *zm1_xoffset,
                                             j + *zm1_yoffset,
                                             k + *zm1_zoffset,1)];
        vzR = velm[CCTK_VECTGFINDEX3D(cctkGH,i + *zm1_xoffset,
                                             j + *zm1_yoffset,
                                             k + *zm1_zoffset,2)];
         
         
        double vxLc = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]; 
        double vyLc = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]; 
        double vzLc = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]; 
        double vxRc = vel[CCTK_VECTGFINDEX3D(cctkGH,i + *zm1_xoffset,
                                                    j + *zm1_yoffset,
                                                    k + *zm1_zoffset,0)];
        double vyRc = vel[CCTK_VECTGFINDEX3D(cctkGH,i + *zm1_xoffset,
                                                    j + *zm1_yoffset,
                                                    k + *zm1_zoffset,1)];
        double vzRc = vel[CCTK_VECTGFINDEX3D(cctkGH,i + *zm1_xoffset,
                                                    j + *zm1_yoffset,
                                                    k + *zm1_zoffset,2)];
         
        if (use_zerov) {
          vxL  = 0.0; vyL  = 0.0; vzL  = 0.0;
          vxR  = 0.0; vyR  = 0.0; vzR  = 0.0;
          vxLc = 0.0; vyLc = 0.0; vzLc = 0.0;
          vxRc = 0.0; vyRc = 0.0; vzRc = 0.0;
        }

        double lambda = 1.0; 
        if (do_fixed_closure || force_radial) lambda = 0.0; 
        Closure clL(vxL,vyL,vzL,gammac,lambda,rhx,rhy,rhz);
        Closure clR(vxR,vyR,vzR,gammac,lambda,rhx,rhy,rhz);
        Closure clLc(vxLc,vyLc,vzLc,gammaLc,lambda,rhx,rhy,rhz);
        Closure clRc(vxRc,vyRc,vzRc,gammaRc,lambda,rhx,rhy,rhz);
        
        clL.setDebugInfo("L"); 
        clR.setDebugInfo("R"); 
        clLc.setDebugInfo("Lc"); 
        clRc.setDebugInfo("Rc"); 
       
        for(int ig=0;ig<ngroups*nspecies;ig++) {
  
          int index4D  = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
          int index4Dm = CCTK_VectGFIndex3D(cctkGH,i + *zm1_xoffset,
                            j + *zm1_yoffset,k + *zm1_zoffset,ig);
          
          // Load the reconstructed left and right hand states 
          double EL   = enuplus[index4D];
          double FxL  = fnuxp[index4D];
          double FyL  = fnuyp[index4D];
          double FzL  = fnuzp[index4D]; 
          
          double ER   = enuminus[index4Dm];
          double FxR  = fnuxm[index4Dm];
          double FyR  = fnuym[index4Dm];
          double FzR  = fnuzm[index4Dm]; 
  
#if 0
          if(isnan(EL) || isnan(ER) || isnan(FxL) || isnan(FxR)) {
            if(CCTK_MyProc(cctkGH)==608) 
              CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
              "Riemann rl: %d ig: %d ijk: %d,%d,%d, %15.6E %15.6E %15.6E %15.6E",
               Carpet::reflevel,ig,i,j,k,EL,ER,FxL,FxR);
          }
#endif      
                 
          // Calculate the closure 
          if (do_fixed_closure) {
            xioldxL[index4D]  = clL.setClosure(EL,FxL,FyL,FzL,xi_guess,true);
            xioldxR[index4Dm] = clR.setClosure(ER,FxR,FyR,FzR,xi_guess,true);
          } else {
            if (r[index3D]>freestreaming_radius) {
              xioldxL[index4D]  = clL.setClosure(EL,FxL,FyL,FzL,1.0,true);
              xioldxR[index4Dm] = clR.setClosure(ER,FxR,FyR,FzR,1.0,true);
            } else {  
              xioldxL[index4D]  = clL.setClosure(EL,FxL,FyL,FzL,xioldxL[index4D]);
              xioldxR[index4Dm] = clR.setClosure(ER,FxR,FyR,FzR,xioldxR[index4Dm]);
            }    
          }                                                                  
          // Calculate the zone centered radiation quantities assuming 
          // diffusive regime (they will only be used at high optical depth anyway)
          clLc.setClosure(enupt[index4D],fnuxpt[index4D],fnuypt[index4D],
              fnuzpt[index4D],0.0,true);
          clRc.setClosure(enupt[index4Dm],fnuxpt[index4Dm],fnuypt[index4Dm],
              fnuzpt[index4Dm],0.0,true);
          
          // Calculate the wave speeds using the Eddington factor 
          // found above for either side of the interface
          std::vector<double> lamp = clL.getWaveSpeeds(direction);
          std::vector<double> lamm = clR.getWaveSpeeds(direction);
          lamm.insert(lamm.end(),lamp.begin(),lamp.end()); 

          // Find minimum and maximum wave speeds 
          double psip = max(*std::max_element(lamm.begin(),lamm.end()), 1.0e-10); 
          double psim = min(*std::min_element(lamm.begin(),lamm.end()),-1.0e-10); 
          
          // Calculate the fluxes on either side of the interface
          double fEL,fFxL,fFyL,fFzL;
          clL.getFluxes(direction,&fEL,&fFxL,&fFyL,&fFzL);
          double fER,fFxR,fFyR,fFzR;
          clR.getFluxes(direction,&fER,&fFxR,&fFyR,&fFzR);
  
          double mindiff = 0.e-2;
          double diff = (psip*psim)/(psip-psim);
          if (abs(diff)<mindiff) {
            if (diff > 0.0) diff = mindiff;
            else            diff =-mindiff;
          }
          double eflx  = (psip*fEL  - psim*fER )/(psip-psim) + diff*(ER  - EL ); 
          double Fxflx = (psip*fFxL - psim*fFxR)/(psip-psim) + diff*(FxR - FxL); 
          double Fyflx = (psip*fFyL - psim*fFyR)/(psip-psim) + diff*(FyR - FyL); 
          double Fzflx = (psip*fFzL - psim*fFzR)/(psip-psim) + diff*(FzR - FzL); 
           
          // Mix with the asymptotic flux for high Peclet number
          // The current form is only appropriate for the zero velocity flat space case
          
          double kappa = sqrt((absorb[index4D ] + scat[index4D ]) 
                             *(absorb[index4Dm] + scat[index4Dm]));
          double A = min(1.0,idx/kappa);
          double aa = tanh(idx/kappa); 
          
          // Correct the asymptotic flux ala Jin 
          if (kappa<=0.0) aa = 1.e20;
          double eflxasym;
          if (aa>=1.0) {
              aa = 1.0;
              eflxasym = 0.0;
          } else {
            
            eflxasym = getAsymptoticDiffusiveFlux(kappa,idx,clLc,clRc,gammac);  
            // Calculate the E flux from advection
            if (do_m1_advect) {
              psip = clL.getAsymptoticWaveSpeeds(direction)[0];
              psim = clR.getAsymptoticWaveSpeeds(direction)[0];
              psip  = lamp[0];
              psim  = lamm[0];
                    
              double a,b,c;
              clL.getAsymptoticAdvectiveFluxes(direction,&fEL,&a,&b,&c);
              clR.getAsymptoticAdvectiveFluxes(direction,&fER,&a,&b,&c);
              if (psip>0.0&&psim>0.0) {
                eflxasym  = fEL + eflxasym; 
                //nflxasym  = fNL + nflxasym; 
              } else if (psip<0.0&&psim<0.0) {
                eflxasym  = fER + eflxasym; 
                //nflxasym  = fNR + nflxasym; 
              }
            }         
          }
          
          eflx  = aa*eflx   + (1.0 - aa)*eflxasym;
          Fxflx = A*A*Fxflx + (1.0 - A*A)*(fFxL + fFxR)/2.0;
          Fyflx = A*A*Fyflx + (1.0 - A*A)*(fFyL + fFyR)/2.0;
          Fzflx = A*A*Fzflx + (1.0 - A*A)*(fFzL + fFzR)/2.0;
  
          // Add this flux contribution to the total flux variables
          enuflux[index4D]   = eflx*idx;
          fnuxflux[index4D]  = Fxflx*idx;
          fnuyflux[index4D]  = Fyflx*idx;
          fnuzflux[index4D]  = Fzflx*idx;
  
          // Refluxing: capture the fluxes
          if (do_m1_reflux) { // TODO: re-use fluxes used for evolution ?
          // Indexing order: flux[var][elt][dir]
          const int flux_idx =
            CCTK_VECTGFINDEX3D(cctkGH, i,j,k, flux_ndirs * ig + direction);
          // TODO: Check factors, in particular alpha (lapse),
          // dx,dy,dz (resolution) etc.
          flux[flux_idx + 0 * flux_var_offset] = eflx *idx;
          flux[flux_idx + 1 * flux_var_offset] = Fxflx*idx;
          flux[flux_idx + 2 * flux_var_offset] = Fyflx*idx;
          flux[flux_idx + 3 * flux_var_offset] = Fzflx*idx;
          if (i==40 && j==3 && k==3) {
            if (zm1_verbose) 
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Refluxing: calculating flux %g",
                          (double)flux[flux_idx + 0 * flux_var_offset]);
            }
          }
        } // Loop over energy and species groups 
  } // Loop over zones

  // Symmetrization (by Roland)
  // #pragma omp flush // implied by #pragma omp for
  // TODO: check if it is faster to keep the loop sizes the same as before
  // and skip the first k iteration (so that the caches match)
  #pragma omp for
  for(int k=zm1_ghost;k<nz-zm1_ghost;k++)
    for(int j=zm1_ghost;j<ny-zm1_ghost;j++)
      for(int i=zm1_ghost;i<nx-zm1_ghost;i++)
        for(int ig=0;ig<ngroups*nspecies;ig++) {
          int index4D  = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
          // naming convention differs from GRHydro. index4DM is + *xm1_xoffset.
          int index4Dp = CCTK_VectGFIndex3D(cctkGH,i - *zm1_xoffset,
                           j - *zm1_yoffset,k - *zm1_zoffset,ig);
          // Add this flux contribution to the RHS variables
          enurhs[index4D]  -= (enuflux[index4D] - enuflux[index4Dp]);
          fnuxrhs[index4D] -= (fnuxflux[index4D] - fnuxflux[index4Dp]);
          fnuyrhs[index4D] -= (fnuyflux[index4D] - fnuyflux[index4Dp]);
          fnuzrhs[index4D] -= (fnuzflux[index4D] - fnuzflux[index4Dp]);
  } // Loop over energy and species groups
  } // OpenMP parallel region
} // End of zm1_RiemannHLLE
} // End of namespace ZelmaniM1

