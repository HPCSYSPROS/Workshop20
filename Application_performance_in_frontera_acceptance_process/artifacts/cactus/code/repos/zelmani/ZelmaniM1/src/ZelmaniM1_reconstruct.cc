#include <iostream>
#include <algorithm>
#include <time.h>
#include <math.h>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


using namespace std;

// define signum
double inline signof(double a) { 
  return (a == 0.0) ? 0.0 : (a<0.0 ? -1.0 : 1.0); 
}

// minmod limiter
double inline minmodfunc(double a_in, double b_in) {

  double sig;

  if (a_in<0.0 && b_in<0.0) {
    sig = max(a_in,b_in);
  } else if (a_in>0.0 && b_in>0.0) {
    sig = min(a_in,b_in);
  } else {
    sig = 0.0;
  }

  return sig;
  
  //return 0.5 * ( signof(a_in) + signof(b_in) ) *
  //  MIN(abs(a_in),abs(b_in));

}

inline double minmod3(double sig1, double sig2, double sig3){
  double sig;
  if (sig1<0.0 && sig2<0.0 && sig3<0.0) sig = max(max(sig1,sig2),sig3);
  else if (sig1>0.0 && sig2>0.0 && sig3>0.0) sig = min(min(sig1,sig2),sig3);
  else sig = 0.0;
  return sig;
}

void TVD_minmod_3D(cGH* cctkGH, int nx, int ny, int nz,
		   int xoffset, int yoffset, int zoffset,
		   double* orig, double* bextp, double* bextm) {
  
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for shared(nx,ny,nz)
  for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
    for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
      for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {

        int index3 = CCTK_GFIndex3D(cctkGH,i,j,k);
        int lindex = CCTK_GFIndex3D(cctkGH,i-xoffset,j-yoffset,k-zoffset);
        int rindex = CCTK_GFIndex3D(cctkGH,i+xoffset,j+yoffset,k+zoffset);
   
        double dupw = orig[index3] - orig[lindex];
        double dloc = orig[rindex] - orig[index3];

        bextm[index3] = orig[index3] - 0.5 * minmodfunc(dupw,dloc);
        bextp[index3] = orig[index3] + 0.5 * minmodfunc(dupw,dloc);

      }

  return;
}

void TVD_minmod(cGH* cctkGH, int nx, int ny, int nz, int ng,
		int xoffset, int yoffset, int zoffset,
		double* orig, double* bextp, double* bextm) {
  
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for shared(nx,ny,nz)
  for(int ig=0;ig<ng;ig++) 
    for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
      for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
	for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {

	  int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	  int lindex = CCTK_VectGFIndex3D(cctkGH,i-xoffset,j-yoffset,k-zoffset,ig);
	  int rindex = CCTK_VectGFIndex3D(cctkGH,i+xoffset,j+yoffset,k+zoffset,ig);
     
	  double dupw = orig[index] - orig[lindex];
	  double dloc = orig[rindex] - orig[index];

	  bextm[index] = orig[index] - 0.5 * minmodfunc(dupw,dloc);
	  bextp[index] = orig[index] + 0.5 * minmodfunc(dupw,dloc);

	}

  return;
}

void TVD_Piecewise(cGH* cctkGH, int nx, int ny, int nz,
		int xoffset, int yoffset, int zoffset,
		double* orig, double* bextp, double* bextm) {
  
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for shared(nx,ny,nz)
  for(int ig=0;ig<ngroups*nspecies;ig++) 
    for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
      for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
	for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {

	  int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
          
	  bextm[index] = orig[index];
	  bextp[index] = orig[index];
	   
	}

  return;
}

void TVD_MCLimiter(cGH* cctkGH, int nx, int ny, int nz,
		int xoffset, int yoffset, int zoffset,
		double* orig, double* bextp, double* bextm) {
  
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for shared(nx,ny,nz)
  for(int ig=0;ig<ngroups*nspecies;ig++) 
    for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
      for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
	for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {

	  int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	  int lindex = CCTK_VectGFIndex3D(cctkGH,i-xoffset,j-yoffset,k-zoffset,ig);
	  int rindex = CCTK_VectGFIndex3D(cctkGH,i+xoffset,j+yoffset,k+zoffset,ig);
     
          double sig1 = orig[index]  - orig[lindex];
	  double sig2 = orig[rindex] - orig[index];
	  double sig3 = orig[rindex] - orig[lindex];
	  bextm[index] = orig[index] - 0.5 * minmod3(2.0*sig1,2.0*sig2,0.5*sig3);
	  bextp[index] = orig[index] + 0.5 * minmod3(2.0*sig1,2.0*sig2,0.5*sig3);

	}

  return;
}

void TVD_MC_f_over_e(cGH* cctkGH, int nx, int ny, int nz,
		int xoffset, int yoffset, int zoffset,
		double* orig, double* bextp, double* bextm, 
		double* norm, double* normp, double* normm) {
  
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for shared(nx,ny,nz)
  for(int ig=0;ig<ngroups*nspecies;ig++) 
    for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
      for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
	for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {
	  
	  int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	  int lindex = CCTK_VectGFIndex3D(cctkGH,i-xoffset,j-yoffset,k-zoffset,ig);
	  int rindex = CCTK_VectGFIndex3D(cctkGH,i+xoffset,j+yoffset,k+zoffset,ig);
     
	  double cent = orig[index]/norm[index];
	  double left = orig[lindex]/norm[lindex];
	  double right = orig[rindex]/norm[rindex];
	    

	  double dupw = cent - left;
	  double dloc = right - cent;

          double sig1 = cent  - left;
	  double sig2 = right - cent;
	  double sig3 = right - left;
	  
	  bextm[index] = (cent - 0.5 * minmod3(2.0*sig1,2.0*sig2,0.5*sig3)) * normm[index];
	  bextp[index] = (cent + 0.5 * minmod3(2.0*sig1,2.0*sig2,0.5*sig3)) * normp[index];
	  
	}

  return;
}

void TVD_minmod_f_over_e(cGH* cctkGH, int nx, int ny, int nz, int ng,
		int xoffset, int yoffset, int zoffset,
		double* orig, double* bextp, double* bextm, 
		double* norm, double* normp, double* normm) {
  
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for shared(nx,ny,nz)
  for(int ig=0;ig<ng;ig++) 
    for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
      for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
	for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {
	  
	  int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	  int lindex = CCTK_VectGFIndex3D(cctkGH,i-xoffset,j-yoffset,k-zoffset,ig);
	  int rindex = CCTK_VectGFIndex3D(cctkGH,i+xoffset,j+yoffset,k+zoffset,ig);
     
	  double cent = orig[index]/norm[index];
	  double left = orig[lindex]/norm[lindex];
	  double right = orig[rindex]/norm[rindex];
	    

	  double dupw = cent - left;
	  double dloc = right - cent;

	  bextm[index] = (cent - 0.5 * minmodfunc(dupw,dloc)) * normm[index];
	  bextp[index] = (cent + 0.5 * minmodfunc(dupw,dloc)) * normp[index];

	}

  return;
}



namespace ZelmaniM1 {

  extern "C" 
  void zm1_Reconstruct(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
    
    double *enupt,*nnupt,*fnuxpt,*fnuypt,*fnuzpt; 
    if (do_m1_RK2 && (*zm1_RKstep == 1)) {
      if (do_m1_nevolve)
        nnupt  = enuh;
      enupt  = enuh;
      fnuxpt = fnuxh;
      fnuypt = fnuyh;
      fnuzpt = fnuzh;
    } else { 
      if (do_m1_nevolve)
        nnupt  = enu;
      enupt  = enu;
      fnuxpt = fnux;
      fnuypt = fnuy;
      fnuzpt = fnuz;
    }

    double *myvelx = &vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,0)];
    double *myvely = &vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,1)];
    double *myvelz = &vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,2)];

    //double *myvelx = vel;
    //double *myvely = vel + nx*ny*nz;
    //double *myvelz = myvely + nx*ny*nz;
    
    // x velocity
    //TVD_minmod_3D(cctkGH, nx, ny, nz,
    //           *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    //           myvelx, velxp, velxm);

    //// y velocity
    //TVD_minmod_3D(cctkGH, nx, ny, nz,
    //           *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    //           myvely, velyp, velym);

    //// z velocity
    //TVD_minmod_3D(cctkGH, nx, ny, nz,
    //	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    //	       myvelz, velzp, velzm);

    TVD_minmod(cctkGH, nx, ny, nz, 3,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	        vel, velp, velm);
    
    if (use_zerov) { 
      for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
        for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++) 
	  for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {
	    int ix = CCTK_VectGFIndex3D(cctkGH,i,j,k,0);
	    int iy = CCTK_VectGFIndex3D(cctkGH,i,j,k,1);
	    int iz = CCTK_VectGFIndex3D(cctkGH,i,j,k,2);
	    velp[ix] = 0.0;
	    velm[ix] = 0.0;
	    velp[iy] = 0.0;
	    velm[iy] = 0.0;
	    velp[iz] = 0.0;
	    velm[iz] = 0.0;
	  }
    } 
     
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"Reconstruct: %d %d cctk_time: %8.5f",*zm1_RKstep,*zm1_flux_direction,cctk_time);

    int ng = ngroups*nspecies;
    if (CCTK_EQUALS(reconstruction_type,"TVD-minmod-FoE")) {
      // Minmod reconstruction, E, F^i/E
      if (do_m1_nevolve)
        TVD_minmod(cctkGH, nx, ny, nz, ng,
    	          *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	           nnupt, nnuplus, nnuminus);
      TVD_minmod(cctkGH, nx,  ny, nz, ng,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       enupt, enuplus, enuminus);
      TVD_minmod_f_over_e(cctkGH, nx,  ny, nz, ng,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    		fnuxpt, fnuxp, fnuxm, enupt, enuplus, enuminus);
      TVD_minmod_f_over_e(cctkGH, nx,  ny, nz, ng,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    		fnuypt, fnuyp, fnuym, enupt, enuplus, enuminus);
      TVD_minmod_f_over_e(cctkGH, nx,  ny, nz, ng,
    			*zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    		fnuzpt, fnuzp, fnuzm, enupt, enuplus, enuminus);
    } else if (CCTK_EQUALS(reconstruction_type,"TVD-minmod")) { 
      // Minmod Reconstruction, E, F^i 
      if (do_m1_nevolve)
        TVD_minmod(cctkGH, nx,  ny, nz, ng,
    	          *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	           nnupt, nnuplus, nnuminus);
      TVD_minmod(cctkGH, nx,  ny, nz, ng,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       enupt, enuplus, enuminus);
      TVD_minmod(cctkGH, nx,  ny, nz, ng,
     	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuxpt, fnuxp, fnuxm);
      TVD_minmod(cctkGH, nx,  ny, nz, ng,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuypt, fnuyp, fnuym);
      TVD_minmod(cctkGH, nx,  ny, nz, ng,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuzpt, fnuzp, fnuzm);
    } else if (CCTK_EQUALS(reconstruction_type,"TVD-MC")) { 
      // MC Reconstruction, E, F^i 
      if (do_m1_nevolve)
        TVD_MCLimiter(cctkGH, nx,  ny, nz,
    	          *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	           nnupt, nnuplus, nnuminus);
      TVD_MCLimiter(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       enupt, enuplus, enuminus);
      TVD_MCLimiter(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuxpt, fnuxp, fnuxm);
      TVD_MCLimiter(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuypt, fnuyp, fnuym);
      TVD_MCLimiter(cctkGH, nx,  ny, nz,
     	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuzpt, fnuzp, fnuzm);
    } else if (CCTK_EQUALS(reconstruction_type,"TVD-MC-FoE")) {
      // Minmod reconstruction, E, F^i/E
      if (do_m1_nevolve)
        TVD_MCLimiter(cctkGH, nx,  ny, nz,
    	          *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	           nnupt, nnuplus, nnuminus);
      TVD_MCLimiter(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       enupt, enuplus, enuminus);
      TVD_MC_f_over_e(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    		fnuxpt, fnuxp, fnuxm, enupt, enuplus, enuminus);
      TVD_MC_f_over_e(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    		fnuypt, fnuyp, fnuym, enupt, enuplus, enuminus);
      TVD_MC_f_over_e(cctkGH, nx,  ny, nz,
    			*zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    		fnuzpt, fnuzp, fnuzm, enupt, enuplus, enuminus);
    } else if (CCTK_EQUALS(reconstruction_type, "PC")) { 
      // do Piecewise reconstruction by default
      if (do_m1_nevolve)
        TVD_Piecewise(cctkGH, nx,  ny, nz,
    	          *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	           nnupt, nnuplus, nnuminus);
      TVD_Piecewise(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       enupt, enuplus, enuminus);
      TVD_Piecewise(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuxpt, fnuxp, fnuxm);
      TVD_Piecewise(cctkGH, nx,  ny, nz,
    	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuypt, fnuyp, fnuym);
      TVD_Piecewise(cctkGH, nx,  ny, nz,
     	       *zm1_xoffset, *zm1_yoffset, *zm1_zoffset,
    	       fnuzpt, fnuzp, fnuzm);
    } else {
      CCTK_WARN(0,"Selected reconstruction type not available!");
    }

    return;
  }

}

