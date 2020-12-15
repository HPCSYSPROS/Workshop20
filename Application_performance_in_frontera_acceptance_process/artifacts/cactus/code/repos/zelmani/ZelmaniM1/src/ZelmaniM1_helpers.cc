#include <iostream>
#include <algorithm>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <util_Table.h>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Symmetry.h"
#include "ZelmaniM1_Closure.hh"
#include "ZelmaniM1_Metric.hh"

#define INV_RHO_GF 6.18735016707159e17
#define PRESS_NU_CONSTANT 3.52127727e24
#define PI 3.1415926535897ed0
#define PI4 97.4090910340024e0
#define PI2 9.86960440108936e0
#define MEV_TO_ERG 1.60217733e-6
#define PRESS_GF 1.80123683248503e-39

using namespace std;


namespace ZelmaniM1 {

  extern "C" 
  void zm1_Initialize(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // set stuff to zero initially
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
#pragma omp parallel for
    for(int ig=0;ig<ngroups*nspecies;ig++) 
      for(int k=0;k<nz;k++)
        for(int j=0;j<ny;j++)
          for(int i=0;i<nx;i++) {

            int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
            int index3D = CCTK_GFIndex3D(cctkGH,i,j,k);
            
	    bound_table = -1;
	    
	    if (do_m1_nevolve) {
	      nnu[index]  = 1.e-40;
	      nnuh[index]  = 1.e-40;
	      nnuplus[index] = 0.0;
	      nnuminus[index] = 0.0;
	    }
	    
	    if (do_m1_pnu)
	      pnu[index3D] = 0.0;
	    
	    enu[index]  = 1.e-40;
            fnux[index] = 0.0;
            fnuy[index] = 0.0;
            fnuz[index] = 0.0;
	    
	    enuh[index]  = 1.e-40;
            fnuxh[index] = 0.0;
            fnuyh[index] = 0.0;
            fnuzh[index] = 0.0;
	    
	    enuplus[index] = 0.0;
	    enuminus[index] = 0.0;
	    fnuxp[index] = 0.0;
	    fnuxm[index] = 0.0;
	    fnuyp[index] = 0.0;
	    fnuym[index] = 0.0;
	    fnuzp[index] = 0.0;
	    fnuzm[index] = 0.0;
	    
	    enurhs[index]  = 0.0;
	    fnuxrhs[index] = 0.0;
	    fnuyrhs[index] = 0.0;
	    fnuzrhs[index] = 0.0;
            
	    absorb[index] = 0.0;
	    emis[index]   = 0.0;
	    scat[index]   = 0.0;
	    
	    xioldxL[index] = -1.0;
	    xioldxR[index] = -1.0;
	    xioldyL[index] = -1.0;
	    xioldyR[index] = -1.0;
	    xioldzL[index] = -1.0;
	    xioldzR[index] = -1.0;
	    
	    velxp[index3D] = 0.0;
	    velxm[index3D] = 0.0;
	    velyp[index3D] = 0.0;
	    velym[index3D] = 0.0;
	    velzp[index3D] = 0.0;
	    velzm[index3D] = 0.0;

	    heatcool[index3D] = 0.0;
	    netheat[index3D] = 0.0;
	    netcool[index3D] = 0.0;

          }

    if(use_nuTmunu) {
#pragma omp parallel for
      for(int k=0;k<nz;k++)
        for(int j=0;j<ny;j++)
          for(int i=0;i<nx;i++) {
            int index3D = CCTK_GFIndex3D(cctkGH,i,j,k);
            eTnutt[index3D] = 0.;
            eTnutx[index3D] = 0.;
            eTnuty[index3D] = 0.;
            eTnutz[index3D] = 0.;
            eTnuxx[index3D] = 0.;
            eTnuxy[index3D] = 0.;
            eTnuxz[index3D] = 0.;
            eTnuyy[index3D] = 0.;
            eTnuyz[index3D] = 0.;
            eTnuzz[index3D] = 0.;
          }
    }
    
    CCTK_VInfo(CCTK_THORNSTRING,"Initialized Radiation.");
    
    return;
  }

  extern "C" 
  void zm1_ZeroRHS(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"ZeroRHS    : %d   cctk_time: %8.5f",*zm1_RKstep,cctk_time);
    
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
#pragma omp parallel for
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++)
	  for(int i=0;i<nx;i++) {
	    
	    int index  = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
    	    if (do_m1_nevolve) nnurhs[index] = 0.0;
	    enurhs[index]  = 0.0;
	    fnuxrhs[index] = 0.0;
	    fnuyrhs[index] = 0.0;
	    fnuzrhs[index] = 0.0;
    
    }
    return;
  }

  extern "C" 
  void zm1_UpdateEoS(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    //DECLARE_CCTK_FUNCTIONS;

    if (zm1_verbose)
      CCTK_VInfo(CCTK_THORNSTRING,"UpdateEoS  : %d   cctk_time: %8.5f",*zm1_RKstep,cctk_time);
    
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
  
    CCTK_INT eoskey = 4;
    CCTK_INT havetemp = 0;
    CCTK_REAL rf_precision = 1.e-10;
    CCTK_INT npts = 1;
    
    CCTK_REAL ss,pp; 
    CCTK_REAL cs2temp, dedttemp, dpderhotemp;
    CCTK_REAL dpdrhoetemp, munutemp;
    CCTK_INT  keyerr,anyerr;


#pragma omp parallel for private(ss,pp,cs2temp,dedttemp,dpderhotemp,dpdrhoetemp,munutemp,keyerr,anyerr)
    for(int k=0;k<nz;k++) 
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++) {
          
          int index = CCTK_GFIndex3D(cctkGH,i,j,k);
	  EOS_Omni_short(eoskey,
	  	         havetemp,
		         rf_precision,
		         npts,
		         &rho[index],
		         &eps[index],
		         &temperature[index],
		         &Y_e[index],
			 &pp,
			 &ss,
		         &cs2temp,
			 &dedttemp,
			 &dpderhotemp,
			 &dpdrhoetemp,
			 &munutemp,
		         &keyerr,&anyerr);
    
    }
    return;
  }

  extern "C" 
  void zm1_StartRKLoop(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    if (do_m1_RK2) *zm1_RKstep  = 2;
    else *zm1_RKstep = 1;
    *zm1_RKstep_half = 1;
    *zm1_RKstep_full = 0;
    return;
  }
  
  extern "C" 
  void zm1_AdvanceRKLoop(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    (*zm1_RKstep)--;
    if (*zm1_RKstep==1) {
      *zm1_RKstep_half = 0;
      *zm1_RKstep_full = 1;
    } else {
      *zm1_RKstep_half = 1;
      *zm1_RKstep_full = 0;
    }
    return;

  }

  extern "C" 
  void zm1_StartLoop(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Note that we can't use c-style loop counting,
    // because scheduler does not like it.
    *zm1_flux_direction = 3;
    *zm1_xoffset = 0;
    *zm1_yoffset = 0;
    *zm1_zoffset = 1;

    return;
  }


  extern "C" 
  void zm1_AdvanceLoop(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    (*zm1_flux_direction)--;

    if(*zm1_flux_direction == 2) {
      *zm1_xoffset = 0;
      *zm1_yoffset = 1;
      *zm1_zoffset = 0;
    } else if(*zm1_flux_direction == 1) {
      *zm1_xoffset = 1;
      *zm1_yoffset = 0;
      *zm1_zoffset = 0;
    } else {
      *zm1_xoffset = -10000;
      *zm1_yoffset = -10000;
      *zm1_zoffset = -10000;
    }
    return;
  }

  extern "C"
  void zm1_ZeroVelocityCurv(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"Zero velocity and curvature (cctk_time: %8.5f)",cctk_time);
    
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
      for(int k=1;k<nz;k++) 
	for(int j=1;j<ny;j++)
	  for(int i=1;i<nx;i++) {
	    int idx2 = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
	    ////double vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
	    //double vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
	    //double vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
            //CCTK_VInfo(CCTK_THORNSTRING,"vel[%d] = %e => 0.0 (i=%d,j=%d,k=%d) ",idx2,vx,i,j,k);
	    vel[idx2             ] = 0.0;
	    vel[idx2 +   nx*ny*nz] = 0.0;
	    vel[idx2 + 2*nx*ny*nz] = 0.0;
    
            kxx[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kxy[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kxz[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kyy[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kyz[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kzz[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
    }
    return;
  }
  
  extern "C"
  void zm1_ComputeNuTmunu(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"Compute Tmunu (cctk_time: %8.5f)",cctk_time);
    
    const double CONV =  2.88114021e-6; //Convert from MeV/fm^3 to M_SUN^-2

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];

#pragma omp parallel for 
      for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++)
	  for(int i=0;i<nx;i++) {
	    
	    int i3D = CCTK_GFIndex3D(cctkGH,i,j,k);
	    ThreeMetric gamma;
	    double alpha = alp[i3D];
	    gamma.loadG(gxx[i3D],gxy[i3D],gxz[i3D],gyy[i3D],gyz[i3D],gzz[i3D]);
	    gamma.loadLapseShift(alpha,betax[i3D]);
	    double sdetg = sqrt(gamma.getDetg()); 
	    double vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
	    double vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
	    double vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
	    
	    if (use_zerov) {
	      vx = 0.0;
	      vy = 0.0;
	      vz = 0.0;
	    }

	    Closure cl;
	    cl.loadBackground(vx,vy,vz,gamma);
	    
	    for (int ig=0;ig<ngroups*nspecies;ig++) {
	      
	      int i4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);

	      // Calculate the neutrino pressure 3-tensor in the lab frame
	      double Pxx,Pxy,Pxz,Pyy,Pyz,Pzz;
	      cl.loadRad(enu[i4D],fnux[i4D],fnuy[i4D],fnuz[i4D]);
	      cl.setClosure(xioldzL[i4D]);
	      cl.getPll(&Pxx,&Pxy,&Pxz,&Pyy,&Pyz,&Pzz); 
	      
	      eTnutt[i3D] += alpha*alpha*enu[i4D]/sdetg;
	      eTnutx[i3D] -= alpha*fnux[i4D]/sdetg;
	      eTnuty[i3D] -= alpha*fnuy[i4D]/sdetg;
	      eTnutz[i3D] -= alpha*fnuz[i4D]/sdetg;
	      eTnuxx[i3D] += Pxx/sdetg;
	      eTnuxy[i3D] += Pxy/sdetg;
	      eTnuxz[i3D] += Pxz/sdetg;
	      eTnuyy[i3D] += Pyy/sdetg;
	      eTnuyz[i3D] += Pyz/sdetg;
	      eTnuzz[i3D] += Pzz/sdetg;
	    }
    }

    // need to call boundary conditions for Tmunu. Hard coded Dirichlet with
    // value = 0 are just fine.
    int ierr = 0;
    const int faces = CCTK_ALL_FACES;
    assert(cctk_nghostzones[0] == cctk_nghostzones[1] &&
           cctk_nghostzones[0] == cctk_nghostzones[2]);
    ierr += Boundary_SelectGroupForBC(cctkGH, faces, cctk_nghostzones[0]-zm1_ghost, -1,
				      "ZelmaniM1::stress_energy_scalar_nu","scalar");
    ierr += Boundary_SelectGroupForBC(cctkGH, faces, cctk_nghostzones[0]-zm1_ghost, -1,
				      "ZelmaniM1::stress_energy_vector_nu","scalar");
    ierr += Boundary_SelectGroupForBC(cctkGH, faces, cctk_nghostzones[0]-zm1_ghost, -1,
				      "ZelmaniM1::stress_energy_tensor_nu","scalar");
    if (ierr)
      CCTK_ERROR("Problems registering boundary condition 'scalar' for evolution vars ZelmaniM1::stress_energy_X_nu");

    return;
  }


  extern "C"
  void zm1_UpdateTmunu(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    if (zm1_verbose)
      CCTK_VInfo(CCTK_THORNSTRING,"Update Tmunu (cctk_time: %8.5f)",cctk_time);

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];

#pragma omp parallel for
      for(int k=0;k<nz;k++)
	for(int j=0;j<ny;j++)
	  for(int i=0;i<nx;i++) {

	    int i3D = CCTK_GFIndex3D(cctkGH,i,j,k);

	    eTtt[i3D] += eTnutt[i3D];
	    eTtx[i3D] += eTnutx[i3D];
	    eTty[i3D] += eTnuty[i3D];
	    eTtz[i3D] += eTnutz[i3D];
	    eTxx[i3D] += eTnuxx[i3D];
	    eTxy[i3D] += eTnuxy[i3D];
	    eTxz[i3D] += eTnuxz[i3D];
	    eTyy[i3D] += eTnuyy[i3D];
	    eTyz[i3D] += eTnuyz[i3D];
	    eTzz[i3D] += eTnuzz[i3D];
    }
    return;
  }

  extern "C"
  void zm1_CopyNuTmunu(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const size_t nbytes = size_t(sizeof(*eTnutt) * cctk_ash[0] *
                                 cctk_ash[1] * cctk_ash[2]);

    // OpenMP parallelization does not help here since this is purely memory
    // bound
    memcpy(eTnutt, eTnutt_p, nbytes);
    memcpy(eTnutx, eTnutx_p, nbytes);
    memcpy(eTnuty, eTnuty_p, nbytes);
    memcpy(eTnutz, eTnutz_p, nbytes);
    memcpy(eTnuxx, eTnuxx_p, nbytes);
    memcpy(eTnuxy, eTnuxy_p, nbytes);
    memcpy(eTnuxz, eTnuxz_p, nbytes);
    memcpy(eTnuyy, eTnuyy_p, nbytes);
    memcpy(eTnuyz, eTnuyz_p, nbytes);
    memcpy(eTnuzz, eTnuzz_p, nbytes);
  }
  
  extern "C" 
  void zm1_pnu(CCTK_ARGUMENTS) { 
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
   
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];

    double i2dx = 0.5/CCTK_DELTA_SPACE(0); 
    double i2dy = 0.5/CCTK_DELTA_SPACE(1); 
    double i2dz = 0.5/CCTK_DELTA_SPACE(2); 
    
    
    if (pnu_equil) { 
#pragma omp parallel for   
      for(int k=zm1_ghost-1;k<nz-zm1_ghost+1;k++) 
	for(int j=zm1_ghost-1;j<ny-zm1_ghost+1;j++)
	  for(int i=zm1_ghost-1;i<nx-zm1_ghost+1;i++) {
    
            CCTK_INT eoskey  = 4;
            CCTK_INT keytemp = 1;
            CCTK_REAL rf_precision = 1.e-10;
            CCTK_INT npts = 1;
            
            CCTK_REAL epstemp,Tout,ss,pp; 
            CCTK_REAL cs2temp, dedttemp, dpderhotemp;
            CCTK_REAL dpdrhoetemp, munu;
            CCTK_INT  keyerr,anyerr;
            
            double F3const  = 7.0*PI4/60.0;
            double pnuconst = PRESS_NU_CONSTANT*PRESS_GF;
	    
	    int i3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    pnu[i3D] = 0.0;		    
	    
	    if (rho[i3D]*INV_RHO_GF>0.01*pnu_dens) {

	      EOS_Omni_short(eoskey,keytemp,rf_precision,npts,
	                   &rho[i3D],&epstemp,&temperature[i3D],
	                   &Y_e[i3D],&pp,&ss,&cs2temp,&dedttemp,
	          	   &dpderhotemp,&dpdrhoetemp,&munu,
	                   &keyerr,&anyerr);
              
	      if (anyerr != 0) 
	        CCTK_WARN(0,"Bad EoS call in pnu equilibrium");
	      
	      double eta = munu/temperature[i3D];
	      if (abs(eta)>10.0) eta = 10.0;
	      double F3  = F3const + 0.5*eta*eta*(PI2 + 0.5*eta*eta);
	      pnu[i3D]   = F3 * pnuconst * pow(temperature[i3D],4.0);
	    }
      }    
    }

#pragma omp parallel for   
    for(int k=zm1_ghost;k<nz-zm1_ghost;k++) 
      for(int j=zm1_ghost;j<ny-zm1_ghost;j++)
        for(int i=zm1_ghost;i<nx-zm1_ghost;i++) {
          
          int i4Dx  = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
          int i4Dy  = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
          int i4Dz  = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);
          
          int i3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
          int i3Dxm = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
          int i3Dxp = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
          int i3Dym = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
          int i3Dyp = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
          int i3Dzm = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
          int i3Dzp = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
	      
	  double sdetg = 1.0;
	   
	  if (pnu_equil) {
	         sdetg = sqrt(2.0*gxy[i3D]*gxz[i3D]*gyz[i3D] + gxx[i3D]*gyy[i3D]*gzz[i3D] 
	                        - gyy[i3D]*gxz[i3D]*gxz[i3D] - gxx[i3D]*gyz[i3D]*gyz[i3D] 
	                     	  - gzz[i3D]*gxy[i3D]*gxy[i3D]);
          }
	   
          double pnu_rolloff = tanh(rho[i3D]*INV_RHO_GF/pnu_dens);
          
	  if (rho[i3D]*INV_RHO_GF>10.0*pnu_dens) {
		  pnu_rolloff = 1.0;
	  } else if (rho[i3D]*INV_RHO_GF<0.1*pnu_dens) {
	          pnu_rolloff = 0.0;
	  }
          
	  srhs[i4Dx] -= (pnu[i3Dxp] - pnu[i3Dxm])*i2dx*sdetg*alp[i3D]*pnu_rolloff; 
          srhs[i4Dy] -= (pnu[i3Dyp] - pnu[i3Dym])*i2dy*sdetg*alp[i3D]*pnu_rolloff; 
          srhs[i4Dz] -= (pnu[i3Dzp] - pnu[i3Dzm])*i2dz*sdetg*alp[i3D]*pnu_rolloff; 
          
    }   
  }
 
  extern "C"
  void zm1_RegisterVars(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int ierr=0;

    CCTK_INFO("Registering evolution variables.");
    if (do_m1_nevolve) 
      ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ZelmaniM1::nnu"));
    if (do_m1_pnu) 
      ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ZelmaniM1::pnu"));
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ZelmaniM1::enu"));
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ZelmaniM1::fnu"));
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ZelmaniM1::heatcoolanalysis"));
    if (do_m1_reflux) {
      const int flux_ndirs = 3;   // x y z
      const int flux_nelts = ngroups*nspecies;
      const int flux_nvars = 4;   // enu fnux fnuy fnuz
      const int vi_register_fine =
        CCTK_VarIndex("Refluxing::register_fine[0]");
      const int vi_register_coarse =
        CCTK_VarIndex("Refluxing::register_coarse[0]");
      for (int n=0; n < flux_nvars * flux_nelts * flux_ndirs; ++n) {
        ierr += MoLRegisterConstrained(vi_register_fine + n);
        ierr += MoLRegisterConstrained(vi_register_coarse + n);
      }
    }
    if (ierr) CCTK_WARN(0,"Problems registering variables with MoL");

    return;
  }

  extern "C"
  void zm1_Boundaries(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int ierr=0;

    const int faces = CCTK_ALL_FACES;

    // boundary condition
    if (bound_table<0 && CCTK_EQUALS("Robin",zm1_bound_type)) 
      bound_table = Util_TableCreateFromString ("DECAY_POWER = 2");
    
    if (*zm1_RKstep == 1) {
      ierr += EnableProlongating(1);
      if (do_m1_nevolve)
        ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
               			      "ZelmaniM1::nnu",zm1_bound_type);
      if (do_m1_pnu)
        ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
               			      "ZelmaniM1::pnu",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::enu",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::fnu",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::heatcoolanalysis","flat");
      if (zm1_verbose)
        CCTK_INFO("RK2 Full step boundary condition");
    } else if (*zm1_RKstep == 2) { 
      ierr += EnableProlongating(0);
      if (do_m1_nevolve) 
        ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::nnuh",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::enuh",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::fnuh",zm1_bound_type);
      if (zm1_verbose)
        CCTK_INFO("RK2 Half step boundary condition");
    } else { 
      // This should be called during POSTREGRID and POSTRESTRICT 
      if (do_m1_nevolve) {
        ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
               			      "ZelmaniM1::nnu",zm1_bound_type);
      }
      
      if (do_m1_pnu)
        ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
               			      "ZelmaniM1::pnu",zm1_bound_type);
      
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::enu",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::fnu",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
           			      "ZelmaniM1::heatcoolanalysis","flat");
      if (zm1_verbose)
        CCTK_INFO("All boundary conditions");

    }
    
    if (ierr) {
      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                "Problems registering boundary condition '%s' for evolution vars ZelmaniM1::enu and ZelmaniM1::fnu",
                zm1_bound_type);
    }

    return;
  }

}

