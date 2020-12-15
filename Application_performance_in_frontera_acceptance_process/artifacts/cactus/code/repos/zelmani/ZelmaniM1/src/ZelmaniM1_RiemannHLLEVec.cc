#include <algorithm>
#include <cassert>
#include <iostream>
#include <time.h>
#include "ZelmaniM1.hh"
#include "ZelmaniM1_AsymptoticF.hh"
#include "ZelmaniM1_ClosureVec.hh"
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
  void zm1_RiemannHLLE(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
    
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"RiemannHLLE: %d %d cctk_time: %8.5f",*zm1_RKstep,*zm1_flux_direction,cctk_time);

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
      double * restrict alph;
      double * restrict betaxh, * restrict betayh, * restrict betazh;
      double * restrict gxxh, * restrict gxyh, * restrict gxzh;
      double * restrict gyyh, * restrict gyzh, * restrict gzzh;
      double * restrict vxR, * restrict vyR, * restrict vzR;
      double * restrict lamLmin, * restrict lamLmax;
      double * restrict lamRmin, * restrict lamRmax;
      double * restrict NL, * restrict NR;
      double * restrict fNL, * restrict fNR;
      double * restrict fEL, * restrict fER;
      double * restrict fFxL, * restrict fFyL, * restrict fFzL;
      double * restrict fFxR, * restrict fFyR, * restrict fFzR;
      
      int xsize = nx - 2*zm1_ghost + 1;

      alph = new double[xsize];
      betaxh = new double[xsize];
      betayh = new double[xsize];
      betazh = new double[xsize];
      gxxh = new double[xsize];
      gxyh = new double[xsize];
      gxzh = new double[xsize];
      gyyh = new double[xsize];
      gyzh = new double[xsize];
      gzzh = new double[xsize];
//      vxR = new double[xsize];
//      vyR = new double[xsize];
//      vzR = new double[xsize];
      lamLmin = new double[xsize];
      lamLmax = new double[xsize];
      lamRmin = new double[xsize];
      lamRmax = new double[xsize];
      NL = new double[xsize];
      NR = new double[xsize];
      fNL = new double[xsize];
      fNR = new double[xsize];
      fEL = new double[xsize];
      fER = new double[xsize];
      fFxL = new double[xsize];
      fFyL = new double[xsize];
      fFzL = new double[xsize];
      fFxR = new double[xsize];
      fFyR = new double[xsize];
      fFzR = new double[xsize];

      ClosureVec clL(xsize), clR(xsize);

#pragma omp for
      // Iteration over zone edges
      // Start on the last ghost zone and end on the last real zone
      for(int k=zm1_ghost-1;k<nz-zm1_ghost;k++) 
	for(int j=zm1_ghost-1;j<ny-zm1_ghost;j++) {

	  for(int i=zm1_ghost-1;i<nx-zm1_ghost;i++) {
	    
	    int index3D  = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index3Dm = CCTK_GFINDEX3D(cctkGH,i + *zm1_xoffset,
					  j + *zm1_yoffset,k + *zm1_zoffset);
	   
	    // Interpolate lapse, shift, and spatial metric to both 
	    // edges of zone (i,j,k) in given offset direction
	    // beta has a raised index 
	    // gamma has two lowered indices
	   
	    int ii = i-zm1_ghost+1;
	    alph[ii] = 0.5 * (alp[index3Dm] + alp[index3D]);
						
	    betaxh[ii] = 0.5 * (betax[index3Dm] + betax[index3D]);
	    betayh[ii] = 0.5 * (betay[index3Dm] + betay[index3D]);
	    betazh[ii] = 0.5 * (betaz[index3Dm] + betaz[index3D]);
		
	    gxxh[ii] = 0.5 * (gxx[index3Dm] + gxx[index3D]); 
	    gxyh[ii] = 0.5 * (gxy[index3Dm] + gxy[index3D]); 
	    gxzh[ii] = 0.5 * (gxz[index3Dm] + gxz[index3D]); 
	    gyyh[ii] = 0.5 * (gyy[index3Dm] + gyy[index3D]); 
	    gyzh[ii] = 0.5 * (gyz[index3Dm] + gyz[index3D]); 
	    gzzh[ii] = 0.5 * (gzz[index3Dm] + gzz[index3D]); 
	  }
	  
          int ii  = zm1_ghost-1;
          int iim = ii + *zm1_xoffset;
          int jjm = j  + *zm1_yoffset;
          int kkm = k  + *zm1_zoffset;
	  
	  clL.loadBackground(&velp[CCTK_VECTGFINDEX3D(cctkGH,ii,j,k,0)],
                             &velp[CCTK_VECTGFINDEX3D(cctkGH,ii,j,k,1)],
                             &velp[CCTK_VECTGFINDEX3D(cctkGH,ii,j,k,2)],
                             alph, betaxh, betayh, betazh,
                             gxxh, gxyh, gxzh, gyyh, gyzh, gzzh);
          clR.loadBackground(&velm[CCTK_VECTGFINDEX3D(cctkGH,iim,jjm,kkm,0)],
                             &velm[CCTK_VECTGFINDEX3D(cctkGH,iim,jjm,kkm,1)],
                             &velm[CCTK_VECTGFINDEX3D(cctkGH,iim,jjm,kkm,2)],
                             alph, betaxh, betayh, betazh,
                             gxxh, gxyh, gxzh, gyyh, gyzh, gzzh);

// To be implemented for the vector class.
//	    if (use_zerov) {
//	      vxL  = 0.0; vyL  = 0.0; vzL  = 0.0;
//	      vxR  = 0.0; vyR  = 0.0; vzR  = 0.0;
//	      vxLc = 0.0; vyLc = 0.0; vzLc = 0.0;
//	      vxRc = 0.0; vyRc = 0.0; vzRc = 0.0;
//	    }
	    
	  for(int ig=0;ig<ngroups*nspecies;ig++) {

	    int index4D  = CCTK_VectGFIndex3D(cctkGH,ii ,j  ,k  ,ig);
	    int index4Dm = CCTK_VectGFIndex3D(cctkGH,iim,jjm,kkm,ig);

	    // Load the reconstructed left and right hand states 
            clL.loadRad(&enuplus[index4D],
                        &fnuxp[index4D],
                        &fnuyp[index4D],
                        &fnuzp[index4D]);
            clR.loadRad(&enuminus[index4Dm],
                        &fnuxm[index4Dm],
                        &fnuym[index4Dm],
                        &fnuzm[index4Dm]);
            
#pragma ivdep
            for(int i=zm1_ghost-1;i<nx-zm1_ghost;i++) {
	      int ii = i-zm1_ghost+1;
	      NL[ii] = 0.0;
	      NR[ii] = 0.0;
	      if (do_m1_nevolve) {
	        NL[ii] = nnuplus[index4D];
	        NR[ii] = nnuminus[index4Dm];
	      }	
	    }
     
	    // Calculate the closure and pressure tensor on either side 
	    // of the interface and at cell centers 

// Fix passing in the right grid function for different directions.
            if (direction==0) {
              clL.setClosureVecZV(&xioldxL[index4D]);
              clR.setClosureVecZV(&xioldxR[index4Dm]);
            } else if (direction==1) {
              clL.setClosureVecZV(&xioldyL[index4D]);
              clR.setClosureVecZV(&xioldyR[index4Dm]);
            } else if (direction==2) {
              clL.setClosureVecZV(&xioldzL[index4D]);
              clR.setClosureVecZV(&xioldzR[index4Dm]);
	    }
	    
	    clL.setStressTensor();
	    clR.setStressTensor();

	    // Calculate the wave speeds using the Eddington factor 
	    // found above for either side of the interface
	    clL.getWaveSpeeds(direction,lamLmin,lamLmax);
	    clR.getWaveSpeeds(direction,lamRmin,lamRmax);
	       
	    // Calculate the left fluxes
	    clL.getFluxes(direction,NL,fEL,fNL,fFxL,fFyL,fFzL);

	    // Calculate the right fluxes 
	    clR.getFluxes(direction,NR,fER,fNR,fFxR,fFyR,fFzR);
	      
	    // Calculate the flux at the interface using the HLLE solver 
#pragma ivdep
            for(int i=zm1_ghost-1;i<nx-zm1_ghost;i++) {

              int index4D  = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
              int index4Dm = CCTK_VectGFIndex3D(cctkGH,i + *zm1_xoffset,
                                         j + *zm1_yoffset,k + *zm1_zoffset,ig);
              int index3D  = CCTK_GFINDEX3D(cctkGH,i,j,k);
              int index3Dm = CCTK_GFINDEX3D(cctkGH,i + *zm1_xoffset,
                                            j + *zm1_yoffset,k + *zm1_zoffset);
	      int ii = i-zm1_ghost+1;

	      // Find minimum and maximum wave speeds 
	      double psip  = max(max(lamLmax[ii],lamRmax[ii]), 1.0e-20);
	      double psim  = min(min(lamLmin[ii],lamRmin[ii]),-1.0e-20);

              double fELL = fEL[ii];
              double fERL = fER[ii];
              double EL = enuplus[index4D];
              double ER = enuminus[index4Dm];
              double fNLL = fNL[ii];
              double fNRL = fNR[ii];
              double NRL = NR[ii];
              double NLL = NL[ii];
              double fFxLL = fFxL[ii];
              double fFxRL = fFxR[ii];
              double FxL  = fnuxp[index4D];
              double FxR  = fnuxm[index4Dm];
              double fFyLL = fFyL[ii];
              double fFyRL = fFyR[ii];
              double FyL  = fnuyp[index4D];
              double FyR  = fnuym[index4Dm];
              double fFzLL = fFzL[ii];
              double fFzRL = fFzR[ii];
              double FzL  = fnuzp[index4D];
              double FzR  = fnuzm[index4Dm];

	      double eflx,nflx,Fxflx,Fyflx,Fzflx;
              double psiprod = psip*psim;
              double invpsidiff = 1.0/(psip-psim);

	//      if (i==15 && j==15 && k==15) cout << i  << " " << index3D << " " << index4D   << " " 
	//	                                << sqrt(FxL*FxL + FyL*FyL + FzL*FzL)/EL << " " << xioldxL[index4D] << " " 
	//					<< psim << "\n";

	      // Mix with the asymptotic flux for high Peclet number
	      // The current form is only appropriate for the zero velocity flat space case
	      int idxopac  = ig + ngroups*nspecies*index3D;
	      int idxopacm = ig + ngroups*nspecies*index3Dm;
	      idxopac  = index4D;
	      idxopacm = index4Dm;
	      
              double absorbL = absorb[idxopac];
              double absorbR = absorb[idxopacm];
              double scatL = scat[idxopac];
              double scatR = scat[idxopacm];

	      double kappa = sqrt((absorbL + scatL) 
	      	                 *(absorbR + scatR));
	      double kappas = sqrt(scatL*scatR);
	      double kappaa = sqrt(absorbL*absorbR);
	      

              double idxdivkappa = ( (kappa==0.0)?1.0:idx/kappa);
	      double A = min(1.0,idxdivkappa);
              double Asq = A*A;
	      
	      // The Audit et al. prescription 
	      eflx  = (psip*fELL - psim*fERL + A*psiprod*(ER - EL))*invpsidiff; 
	      nflx  = (psip*fNLL - psim*fNRL + A*psiprod*(NRL - NLL))*invpsidiff; 
	      Fxflx = A*(A*psip*fFxLL - A*psim*fFxRL + psiprod*(FxR - FxL))*invpsidiff; 
	      Fyflx = A*(A*psip*fFyLL - A*psim*fFyRL + psiprod*(FyR - FyL))*invpsidiff; 
	      Fzflx = A*(A*psip*fFzLL - A*psim*fFzRL + psiprod*(FzR - FzL))*invpsidiff; 
	      Fxflx = Fxflx + (1.0 - Asq)*(fFxLL + fFxRL)*0.5;
	      Fyflx = Fyflx + (1.0 - Asq)*(fFyLL + fFyRL)*0.5;
	      Fzflx = Fzflx + (1.0 - Asq)*(fFzLL + fFzRL)*0.5;
              
	      if (eflx!=eflx) {
#pragma omp critical
		cout << "Eflx bad : " << eflx << " " << psip << " " << psim << " " << fELL << " " <<  fERL << "\n";
	        CCTK_WARN(0,"aborting!");
              }
	      if (Fxflx!=Fxflx){
#pragma omp critical
		cout << "Fxflx bad : " << Fxflx << " " << fFxLL << " " << fFxRL << " " << invpsidiff << "\n";
	        CCTK_WARN(0,"aborting!");
              }
#pragma omp critical
              if (Fyflx!=Fyflx) cout << "Fyflx bad : " << Fyflx << "\n";
#pragma omp critical
              if (Fzflx!=Fzflx) cout << "Fzflx bad : " << Fzflx << "\n";
	      
	      //eflx  = (psip*fELL  - psim*fERL  + psiprod*(ER  - EL ))*invpsidiff; 
	      //nflx  = (psip*fNLL  - psim*fNRL  + psiprod*(NRL  - NLL ))*invpsidiff; 
	      //Fxflx = (psip*fFxLL - psim*fFxRL + psiprod*(FxR - FxL))*invpsidiff; 
	      //Fyflx = (psip*fFyLL - psim*fFyRL + psiprod*(FyR - FyL))*invpsidiff; 
	      //Fzflx = (psip*fFzLL - psim*fFzRL + psiprod*(FzR - FzL))*invpsidiff; 
	      
	      //double aa = tanh(idxdivkappa); 
	      // Correct the asymptotic flux ala Jin 
	      //if (kappa<=0.0) aa = 1.e20;
	      //double eflxasym,nflxasym;
	      //if (aa>=1.0) {
	      //    aa = 1.0;
	      //    eflxasym = 0.0;
	      //} else {
	      //  AsymptoticF asymF;
	      //  asymF.loadAllL(enupt[index4D],fnuxpt[index4D],fnuypt[index4D],fnuzpt[index4D],
	      //  vxL,vyL,vzL,&gammaLc);
	      //  asymF.loadAllR(enupt[index4Dm],fnuxpt[index4Dm],fnuypt[index4Dm],fnuzpt[index4Dm],
	      //  vxR,vyR,vzR,&gammaRc);
	      //  asymF.loadMetricC(&gammac);
	      //     
	      //  eflxasym = asymF.getAsymptoticFa(direction,kappa,idx,&nflxasym);
	      //     
	      //  // Calculate the E flux from advection
	      //  if (do_m1_advect) {
	      //    clL.getAsymptoticWaveSpeeds(direction,lamp);
	      //    clR.getAsymptoticWaveSpeeds(direction,lamm);
	      //    psip  = lamp[0];
	      //    psim  = lamm[0];
                   
	      //    double a,b,c;
	      //    clL.getAsymptoticFluxes(direction,NL,&fEL,&fNL,&a,&b,&c);
	      //    clR.getAsymptoticFluxes(direction,NR,&fER,&fNR,&a,&b,&c);
              //    if (psip>0.0&&psim>0.0) {
	      //      eflxasym  = fEL + eflxasym; 
	      //      nflxasym  = fNL + nflxasym; 
              //    } else if (psip<0.0&&psim<0.0) {
	      //      eflxasym  = fER + eflxasym; 
	      //      nflxasym  = fNR + nflxasym; 
              //    }
	      //  }         
	      //}
	      //
	      //eflx  = aa*eflx   + (1.0 - aa)*eflxasym;
	      //nflx  = aa*nflx   + (1.0 - aa)*nflxasym;
	      //Fxflx = A*A*Fxflx + (1.0 - A*A)*(fFxL + fFxR)/2.0;
	      //Fyflx = A*A*Fyflx + (1.0 - A*A)*(fFyL + fFyR)/2.0;
	      //Fzflx = A*A*Fzflx + (1.0 - A*A)*(fFzL + fFzR)/2.0;

	      // Add this flux contribution to the RHS variables
	      #pragma omp atomic      
              enurhs[index4D]   -= eflx*idx;
	      #pragma omp atomic      
              enurhs[index4Dm]  += eflx*idx;
	      if (do_m1_nevolve) {
	        #pragma omp atomic      
                nnurhs[index4D]   -= nflx*idx;
	        #pragma omp atomic     
                nnurhs[index4Dm]  += nflx*idx;
              }
	      #pragma omp atomic      
	      fnuxrhs[index4D]  -= Fxflx*idx;
	      #pragma omp atomic      
	      fnuxrhs[index4Dm] += Fxflx*idx;
	      #pragma omp atomic      
	      fnuyrhs[index4D]  -= Fyflx*idx;
	      #pragma omp atomic      
	      fnuyrhs[index4Dm] += Fyflx*idx;
	      #pragma omp atomic      
	      fnuzrhs[index4D]  -= Fzflx*idx;
	      #pragma omp atomic      
	      fnuzrhs[index4Dm] += Fzflx*idx;
              
// This doesn't work according to Luke, so commented out for now.
//              // Refluxing: capture the fluxes
//              if (do_m1_reflux) {
//	        // Indexing order: flux[var][elt][dir]
//                const int flux_idx =
//                  CCTK_VECTGFINDEX3D(cctkGH, i,j,k, flux_ndirs * ig + direction);
//                // TODO: Check factors, in particular alpha (lapse),
//                // dx,dy,dz (resolution) etc.
//                flux[flux_idx + 0 * flux_var_offset] = eflx *idx;
//                flux[flux_idx + 1 * flux_var_offset] = Fxflx*idx;
//                flux[flux_idx + 2 * flux_var_offset] = Fyflx*idx;
//                flux[flux_idx + 3 * flux_var_offset] = Fzflx*idx;
//                if (i==40 && j==3 && k==3) {
//                  if (zm1_verbose) 
//		    CCTK_VInfo(CCTK_THORNSTRING,
//                             "Refluxing: calculating flux %g",
//                             (double)flux[flux_idx + 0 * flux_var_offset]);
//                }
//              }

	    } // Loop over x-direction
	  } // Loop over energy and species groups 
        } // Loop over y- and z-direction.
} // OpenMP parallel region
    return;

  }

}
