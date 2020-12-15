#include <algorithm>
#include <cassert>
#include <iostream>
#include <time.h>
#include "ZelmaniM1.hh"
#include "ZelmaniM1_Closure.hh"
#include "ZelmaniM1_Metric.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define C_L 2.99792458e10
#define G_GRAV 6.673e-8
#define M_SUN 1.98892e33

using namespace std;


namespace ZelmaniM1 {

  extern "C" 
  void zm1_FluxesDG(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
    
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"RiemannDG: %d %d cctk_time: %8.5f",*zm1_RKstep,*zm1_flux_direction,cctk_time);

    int direction;
    if      (*zm1_xoffset==1) direction = 0;  
    else if (*zm1_yoffset==1) direction = 1;  
    else if (*zm1_zoffset==1) direction = 2;  
    else direction = -1;
    
    double idx = 1.0;
    double idy = 1.0;
    double idz = 1.0;

    if (direction == 0) {
	idx = 1.0/CCTK_DELTA_SPACE(0);
	idy = 1.0/CCTK_DELTA_SPACE(1);
	idz = 1.0/CCTK_DELTA_SPACE(2);
    } else if (direction == 1) {
	idx = 1.0/CCTK_DELTA_SPACE(1);
	idy = 1.0/CCTK_DELTA_SPACE(2);
	idz = 1.0/CCTK_DELTA_SPACE(0);
    } else if (direction == 2) { 
	idx = 1.0/CCTK_DELTA_SPACE(2);
	idy = 1.0/CCTK_DELTA_SPACE(0);
	idz = 1.0/CCTK_DELTA_SPACE(1);
    }

    // Point at the correct arrays for the current RK step
    double *enupt;
    double *nnupt;
    double *fnuxpt;
    double *fnuypt;
    double *fnuzpt;
    if (do_m1_RK2 && (*zm1_RKstep == 1)) {
	    if (do_m1_nevolve) 
	      nnupt  = nnuh;
	    enupt  = enuh;
	    fnuxpt = fnuxh;
	    fnuypt = fnuyh;
	    fnuzpt = fnuzh;
    
    } else {
	    if (do_m1_nevolve) 
	      nnupt  = nnu;
	    enupt  = enu;
	    fnuxpt = fnux;
	    fnuypt = fnuy;
	    fnuzpt = fnuz;
    }
    
    // Refluxing: capture the fluxes
    const int flux_ndirs = 3;   // x y z
    const int flux_nelts = ngroups*nspecies;
    const int flux_nvars = 4;   // enu fnux fnuy fnuz
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
	    
	    // Check which DG flux to calculate 
	    bool do_interior = true;
	    if (direction == 0 && (i-zm1_ghost+1)%2==0 ) do_interior = false;
	    if (direction == 1 && (j-zm1_ghost+1)%2==0 ) do_interior = false;
	    if (direction == 2 && (k-zm1_ghost+1)%2==0 ) do_interior = false;

	    int index3D  = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index3Dm = CCTK_GFINDEX3D(cctkGH,i + *zm1_xoffset,
					  j + *zm1_yoffset,k + *zm1_zoffset);
	   
	    // Interpolate lapse, shift, and spatial metric to both 
	    // edges of zone (i,j,k) in given offset direction
	    // beta has a raised index 
	    // gamma has two lowered indices
	    ThreeMetric gammaLc,gammaRc,gammac;
	   
	    double alph = 0.5 * (alp[index3Dm] + alp[index3D]);
						
	    double betaxh = 0.5 * (betax[index3Dm] + betax[index3D]);
	    double betayh = 0.5 * (betay[index3Dm] + betay[index3D]);
	    double betazh = 0.5 * (betaz[index3Dm] + betaz[index3D]);
		
	    double gxxh = 0.5 * (gxx[index3Dm] + gxx[index3D]); 
	    double gxyh = 0.5 * (gxy[index3Dm] + gxy[index3D]); 
	    double gxzh = 0.5 * (gxz[index3Dm] + gxz[index3D]); 
	    double gyyh = 0.5 * (gyy[index3Dm] + gyy[index3D]); 
	    double gyzh = 0.5 * (gyz[index3Dm] + gyz[index3D]); 
	    double gzzh = 0.5 * (gzz[index3Dm] + gzz[index3D]); 
	     
	    gammaLc.loadLapseShiftAll(alp[index3D],betax[index3D ],betay[index3D ],betaz[index3D ]);
	    gammaRc.loadLapseShiftAll(alp[index3D],betax[index3Dm],betay[index3Dm],betaz[index3Dm]);
	    gammac. loadLapseShiftAll(alph,betaxh,betayh,betazh);
	    
	    gammaLc.loadG(gxx[index3D],gxy[index3D],gxz[index3D],
	    	      gyy[index3D],gyz[index3D],gzz[index3D]);
	    gammaRc.loadG(gxx[index3Dm],gxy[index3Dm],gxz[index3Dm],
	    	      gyy[index3Dm],gyz[index3Dm],gzz[index3Dm]);
	    gammac.loadG(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh);
	    
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

	    Closure clL,clR;
	    Closure clLc,clRc;
	    
	    clL.loadBackground(vxL,vyL,vzL,gammac);
	    clR.loadBackground(vxR,vyR,vzR,gammac);
	    clLc.loadBackground(vxLc,vyLc,vzLc,gammaLc);
	    clRc.loadBackground(vxRc,vyRc,vzRc,gammaRc);
	    
	    for(int ig=0;ig<ngroups*nspecies;ig++) {

	      int index4Dp  = CCTK_VectGFIndex3D(cctkGH,i - *zm1_xoffset,
	          		j - *zm1_yoffset,k - *zm1_zoffset,ig);
	      int index4D   = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	      int index4Dm  = CCTK_VectGFIndex3D(cctkGH,i + *zm1_xoffset,
	          		j + *zm1_yoffset,k + *zm1_zoffset,ig);
	      int index4Dmm = CCTK_VectGFIndex3D(cctkGH,i + 2*(*zm1_xoffset),
	          		j + 2*(*zm1_yoffset),k + 2*(*zm1_zoffset),ig);
	      
	      // Load the reconstructed left and right hand states 
	      double EL, FxL, FyL, FzL;
	      double ER, FxR, FyR, FzR;
	      double NL = 0.0;
	      double NR = 0.0;

              if (do_interior) {
	        EL   = enupt[index4D];
	        FxL  = fnuxpt[index4D];
	        FyL  = fnuypt[index4D];
	        FzL  = fnuzpt[index4D]; 
	        
	        ER   = enupt[index4Dm];
	        FxR  = fnuxpt[index4Dm];
	        FyR  = fnuypt[index4Dm];
	        FzR  = fnuzpt[index4Dm]; 
	      } else {
	        EL   = 1.5*enupt[index4D]  - 0.5*enupt[index4Dp];
	        FxL  = 1.5*fnuxpt[index4D] - 0.5*fnuxpt[index4Dp];
	        FyL  = 1.5*fnuypt[index4D] - 0.5*fnuypt[index4Dp];
	        FzL  = 1.5*fnuzpt[index4D] - 0.5*fnuzpt[index4Dp];
	        
	        ER   = 1.5*enupt[index4Dm]  - 0.5*enupt[index4Dmm];
	        FxR  = 1.5*fnuxpt[index4Dm] - 0.5*fnuxpt[index4Dmm];
	        FyR  = 1.5*fnuypt[index4Dm] - 0.5*fnuypt[index4Dmm];
	        FzR  = 1.5*fnuzpt[index4Dm] - 0.5*fnuzpt[index4Dmm];
	      }
	      if (do_m1_nevolve) {
	        NL   = nnuplus[index4D];
	        NR   = nnuminus[index4Dm];
	      }	
	             
	      // Calculate the closure and pressure tensor on either side 
	      // of the interface and at cell centers 
	      
	      clL.loadRad(EL,FxL,FyL,FzL);
	      clR.loadRad(ER,FxR,FyR,FzR);
	      clL.setFhat(x[index3D]+0.5/idx,y[index3D],z[index3D]);
	      clR.setFhat(x[index3D]+0.5/idx,y[index3D],z[index3D]);
	      xioldxL[index4D]  = clL.setClosure(xioldxL[index4D]);
	      xioldxR[index4Dm] = clR.setClosure(xioldxR[index4Dm]);
	      
	      // Calculate the wave speeds using the Eddington factor 
	      // found above for either side of the interface
	      double lamp[4];
	      clL.getWaveSpeeds(direction,lamp);
	      double lamm[4];
	      clR.getWaveSpeeds(direction,lamm);
	      
	      // Find minimum and maximum wave speeds 
	      double lminp = min(min(lamp[0],lamp[1]),min(lamp[2],lamp[3]));
	      double lmaxp = max(max(lamp[0],lamp[1]),max(lamp[2],lamp[3]));
	      double lminm = min(min(lamm[0],lamm[1]),min(lamm[2],lamm[3]));
	      double lmaxm = max(max(lamm[0],lamm[1]),max(lamm[2],lamm[3]));
	      double psip  = max(max(lmaxp,lmaxm), 1.0e-20);
	      double psim  = min(min(lminp,lminm),-1.0e-20);
	       
	      // Calculate the left fluxes
	      double fEL,fNL,fFxL,fFyL,fFzL;
	      clL.getFluxes(direction,NL,&fEL,&fNL,&fFxL,&fFyL,&fFzL);

	      // Calculate the right fluxes 
	      double fER,fNR,fFxR,fFyR,fFzR;
	      clR.getFluxes(direction,NR,&fER,&fNR,&fFxR,&fFyR,&fFzR);
	      
	      double eflx,nflx,Fxflx,Fyflx,Fzflx;
	      if (do_interior) {
	        eflx  = 0.5*fEL  + 0.5*fER ; 
	        nflx  = 0.5*fNL  + 0.5*fNR ;
	        Fxflx = 0.5*fFxL + 0.5*fFxR;
	        Fyflx = 0.5*fFyL + 0.5*fFyR;
	        Fzflx = 0.5*fFzL + 0.5*fFzR;
	       
	        // Add this flux contribution to the RHS variables
		double ac = 0.5; // David coefficient (include mass lumping)
		//double ac = 1.5; // re-derived coefficient
	        #pragma omp atomic      
                enurhs[index4D]   -= ac*eflx*idx;
	        #pragma omp atomic      
                enurhs[index4Dm]  += ac*eflx*idx;
	        if (do_m1_nevolve) {
	          #pragma omp atomic      
                  nnurhs[index4D]   -= ac*nflx*idx;
	          #pragma omp atomic     
                  nnurhs[index4Dm]  += ac*nflx*idx;
                }
	        #pragma omp atomic      
	        fnuxrhs[index4D]  -= ac*Fxflx*idx;
	        #pragma omp atomic      
	        fnuxrhs[index4Dm] += ac*Fxflx*idx;
	        #pragma omp atomic      
	        fnuyrhs[index4D]  -= ac*Fyflx*idx;
	        #pragma omp atomic      
	        fnuyrhs[index4Dm] += ac*Fyflx*idx;
	        #pragma omp atomic      
	        fnuzrhs[index4D]  -= ac*Fzflx*idx;
	        #pragma omp atomic      
	        fnuzrhs[index4Dm] += ac*Fzflx*idx;
	      
	      } else {
	        // Calculate the flux at the interface using the HLLE solver 
	        eflx  = (psip*fEL  - psim*fER  + psip*psim*(ER  - EL ))/(psip-psim); 
	        nflx  = (psip*fNL  - psim*fNR  + psip*psim*(NR  - NL ))/(psip-psim); 
	        Fxflx = (psip*fFxL - psim*fFxR + psip*psim*(FxR - FxL))/(psip-psim); 
	        Fyflx = (psip*fFyL - psim*fFyR + psip*psim*(FyR - FyL))/(psip-psim); 
	        Fzflx = (psip*fFzL - psim*fFzR + psip*psim*(FzR - FzL))/(psip-psim); 
	        
		double ac = 0.75; // David coefficients (include mass lumping)
		double af = 0.25;	
		//double ac =  1.25; // re-derived coefficients
		//double af = -0.25;	
		#pragma omp atomic      
                enurhs[index4D]   -= ac*eflx*idx;
		#pragma omp atomic      
                enurhs[index4Dm]  += ac*eflx*idx;
		#pragma omp atomic      
                enurhs[index4Dp]  -= af*eflx*idx;
		#pragma omp atomic      
                enurhs[index4Dmm] += af*eflx*idx;

		#pragma omp atomic      
                fnuxrhs[index4D]   -= ac*Fxflx*idx;
		#pragma omp atomic      
                fnuxrhs[index4Dm]  += ac*Fxflx*idx;
		#pragma omp atomic      
                fnuxrhs[index4Dp]  -= af*Fxflx*idx;
		#pragma omp atomic      
                fnuxrhs[index4Dmm] += af*Fxflx*idx;
		#pragma omp atomic      
                fnuyrhs[index4D]   -= ac*Fyflx*idx;
		#pragma omp atomic      
                fnuyrhs[index4Dm]  += ac*Fyflx*idx;
		#pragma omp atomic      
                fnuyrhs[index4Dp]  -= af*Fyflx*idx;
		#pragma omp atomic      
                fnuyrhs[index4Dmm] += af*Fyflx*idx;
		#pragma omp atomic      
                fnuzrhs[index4D]   -= ac*Fzflx*idx;
		#pragma omp atomic      
                fnuzrhs[index4Dm]  += ac*Fzflx*idx;
		#pragma omp atomic      
                fnuzrhs[index4Dp]  -= af*Fzflx*idx;
		#pragma omp atomic      
                fnuzrhs[index4Dmm] += af*Fzflx*idx;

	      }
	      
              
              // Refluxing: capture the fluxes
              if (do_m1_reflux) {
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
} // OpenMP parallel region
    return;

  }

}

