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
  void zm1_CalcUpdate(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];

    const double dt = CCTK_DELTA_TIME;
    const double idx = 1.0/CCTK_DELTA_SPACE(0); // 1/delta x of zone  
    const double idy = 1.0/CCTK_DELTA_SPACE(1); // 1/delta x of zone  
    const double idz = 1.0/CCTK_DELTA_SPACE(2); // 1/delta x of zone  
    
    double *enupt, *nnupt, *fnuxpt, *fnuypt, *fnuzpt;
    double *enusr, *nnusr, *fnuxsr, *fnuysr, *fnuzsr;
    double *register_fine_pt, *register_coarse_pt;
    
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
    
     
    double dts = dt; // delta t for partial step 
    if (zm1_verbose)
      CCTK_VInfo(CCTK_THORNSTRING,"CalcUpdate : %d   cctk_time: %8.5f",*zm1_RKstep,cctk_time);
    
    if (*zm1_RKstep == 2) {
      dts = 0.5*dt;
      if (do_m1_nevolve) 
        nnupt  = nnuh;
      enupt  = enuh;
      fnuxpt = fnuxh;
      fnuypt = fnuyh;
      fnuzpt = fnuzh;
      if (do_m1_nevolve) 
        nnusr  = nnu;
      enusr  = enu;
      fnuxsr = fnux;
      fnuysr = fnuy;
      fnuzsr = fnuz;
      if (do_m1_reflux) {
        register_fine_pt = register_fine_h;
        register_coarse_pt = register_coarse_h;
      } 
#pragma omp parallel for
      // This is currently implemented in an exceptionally stupid way
      for(int ig=0;ig<ngroups*nspecies;ig++)
        for(int k=0;k<nz;k++) 
	  for(int j=0;j<ny;j++)
	    for(int i=0;i<nx;i++) {
	      int index4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
              // Shouldn't have to do this if we are setting boundary conditions
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
      }
    } else { 
      
      dts = dt;
      if (do_m1_nevolve) 
        nnupt  = nnu;
      enupt  = enu;
      fnuxpt = fnux;
      fnuypt = fnuy;
      fnuzpt = fnuz;
      
      if (do_m1_RK2) {
        if (do_m1_nevolve) 
          nnusr  = nnuh;
        enusr  = enuh;
        fnuxsr = fnuxh;
        fnuysr = fnuyh;
        fnuzsr = fnuzh;
      } else {
        if (do_m1_nevolve) 
          nnusr  = nnu;
        enusr  = enu;
        fnuxsr = fnux;
        fnuysr = fnuy;
        fnuzsr = fnuz;
      }
      
      if (do_m1_reflux) {
        register_fine_pt = register_fine;
        register_coarse_pt = register_coarse;
      }
    }

    double Etotnew  = 0.0;
    double Etotold  = 0.0;
    double Fxtotnew = 0.0;
    double Fxtotold = 0.0;

    double sconfac = 1.0e0;
    if(scon_backreact_taper_on_time > 1.0e-10) {
      sconfac = (1.0e0 - exp(- cctk_time / scon_backreact_taper_on_time));
    }

  const double idt = 1.0/dt; 

#pragma omp parallel
{

#pragma omp for reduction(+:Etotnew,Etotold,Fxtotold,Fxtotnew)  
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
	   
            CCTK_INT eoskey = 4;
            CCTK_INT havetemp = 1;
            CCTK_REAL rf_precision = 1.e-10;
            CCTK_INT npts = 1;
            
            CCTK_REAL ss,pp; 
            CCTK_REAL cs2temp, dedttemp, dpderhotemp;
            CCTK_REAL dpdrhoetemp, epstemp;
            CCTK_INT  keyerr,anyerr;
	    CCTK_REAL munu = 0.0;
	    if (do_gray && *zm1_RKstep<2) 
              EOS_Omni_short(eoskey,
	      	             havetemp,
	    	             rf_precision,
	    	             npts,
	    	             &rho[i3D],
	    	             &epstemp,
	    	             &temperature[i3D],
	    	             &Y_e[i3D],
	    		     &pp,
	    		     &ss,
	    	             &cs2temp,
	    		     &dedttemp,
	    		     &dpderhotemp,
	    		     &dpdrhoetemp,
	    		     &munu,
	    	             &keyerr,&anyerr);

	    ThreeMetric gammac,gammaco;
	    Closure cl;
            
	    gammac.loadG(gxx[i3D],gxy[i3D],gxz[i3D],gyy[i3D],gyz[i3D],gzz[i3D]);
	    gammac.loadExtCurv(kxx[i3D],kxy[i3D],kxz[i3D],kyy[i3D],kyz[i3D],kzz[i3D]);
	    gammaco.loadG(gxx_p[i3D],gxy_p[i3D],gxz_p[i3D],gyy_p[i3D],gyz_p[i3D],gzz_p[i3D]);
            
	    double vx,vy,vz;
	    double dWdx,dWdy,dWdz,dWdt;
	    double dWvxdx,dWvydx,dWvzdx;
	    double dWvxdy,dWvydy,dWvzdy;
	    double dWvxdz,dWvydz,dWvzdz;
	    double dWvxdt,dWvydt,dWvzdt;
            double W,ndW,ndWvx,ndWvy,ndWvz;

	    if (use_zerov) {
	      vx = 0.0;
	      vy = 0.0;
	      vz = 0.0;
	      dWvxdx = 0; dWvxdy = 0; dWvxdz = 0;
	      dWvydx = 0; dWvydy = 0; dWvydz = 0;
	      dWvzdx = 0; dWvzdy = 0; dWvzdz = 0;
	      dWvxdt = 0; dWvydt = 0; dWvzdt = 0;
	      dWdx = 0; dWdy = 0; dWdz = 0; dWdt = 0;
	      W=1; ndW = 0; ndWvx = 0; ndWvy = 0; ndWvz = 0;
	    }  else {
	      vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
	      vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
	      vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
	      double vx_p = vel_p[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
	      double vy_p = vel_p[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
	      double vz_p = vel_p[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
	      W = 1/sqrt(1 -    gxx[i3D]*vx*vx - gyy[i3D]*vy*vy - gzz[i3D]*vz*vz 
	                   - 2*(gxy[i3D]*vx*vy + gxz[i3D]*vx*vz + gyz[i3D]*vy*vz));
	      double W_p = 1/sqrt(1 -    gxx_p[i3D]*vx_p*vx_p - gyy_p[i3D]*vy_p*vy_p - gzz_p[i3D]*vz_p*vz_p 
	                            - 2*(gxy_p[i3D]*vx_p*vy_p + gxz_p[i3D]*vx_p*vz_p + gyz_p[i3D]*vy_p*vz_p));

              if (do_m1_redshift) {
	        double vxp = vel[CCTK_VECTGFINDEX3D(cctkGH,i+1,j,k,0)];
	        double vyp = vel[CCTK_VECTGFINDEX3D(cctkGH,i+1,j,k,1)];
	        double vzp = vel[CCTK_VECTGFINDEX3D(cctkGH,i+1,j,k,2)];
	        double vxm = vel[CCTK_VECTGFINDEX3D(cctkGH,i-1,j,k,0)];
	        double vym = vel[CCTK_VECTGFINDEX3D(cctkGH,i-1,j,k,1)];
	        double vzm = vel[CCTK_VECTGFINDEX3D(cctkGH,i-1,j,k,2)]; 
	        double Wp = 1/sqrt(1 -    gxx[i3Dxp]*vxp*vxp - gyy[i3Dxp]*vyp*vyp - gzz[i3Dxp]*vzp*vzp 
	                             - 2*(gxy[i3Dxp]*vxp*vyp + gxz[i3Dxp]*vxp*vzp + gyz[i3Dxp]*vyp*vzp));
	        double Wm = 1/sqrt(1 -    gxx[i3Dxm]*vxm*vxm - gyy[i3Dxm]*vym*vym - gzz[i3Dxm]*vzm*vzm 
	                             - 2*(gxy[i3Dxm]*vxm*vym + gxz[i3Dxm]*vxm*vzm + gyz[i3Dxm]*vym*vzm));
	        dWvxdx = (Wp*vxp - Wm*vxm)*idx*0.5; 
	        dWvydx = (Wp*vyp - Wm*vym)*idx*0.5; 
	        dWvzdx = (Wp*vzp - Wm*vzm)*idx*0.5; 
                dWdx   = (Wp - Wm)*idx*0.5;
                
	        vxp = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j+1,k,0)];
	        vyp = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j+1,k,1)];
	        vzp = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j+1,k,2)];
	        vxm = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j-1,k,0)];
	        vym = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j-1,k,1)];
	        vzm = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j-1,k,2)]; 
	        Wp = 1/sqrt(1 -    gxx[i3Dyp]*vxp*vxp - gyy[i3Dyp]*vyp*vyp - gzz[i3Dyp]*vzp*vzp 
	                      - 2*(gxy[i3Dyp]*vxp*vyp + gxz[i3Dyp]*vxp*vzp + gyz[i3Dyp]*vyp*vzp));
	        Wm = 1/sqrt(1 -    gxx[i3Dym]*vxm*vxm - gyy[i3Dym]*vym*vym - gzz[i3Dym]*vzm*vzm 
	                      - 2*(gxy[i3Dym]*vxm*vym + gxz[i3Dym]*vxm*vzm + gyz[i3Dym]*vym*vzm));
	        dWvxdy = (Wp*vxp - Wm*vxm)*idy*0.5; 
	        dWvydy = (Wp*vyp - Wm*vym)*idy*0.5; 
	        dWvzdy = (Wp*vzp - Wm*vzm)*idy*0.5; 
                dWdy   = (Wp - Wm)*idy*0.5;
                
	        vxp = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k+1,0)];
	        vyp = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k+1,1)];
	        vzp = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k+1,2)];
	        vxm = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k-1,0)];
	        vym = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k-1,1)];
	        vzm = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k-1,2)]; 
	        Wp = 1/sqrt(1 -    gxx[i3Dzp]*vxp*vxp - gyy[i3Dzp]*vyp*vyp - gzz[i3Dzp]*vzp*vzp 
	                      - 2*(gxy[i3Dzp]*vxp*vyp + gxz[i3Dzp]*vxp*vzp + gyz[i3Dzp]*vyp*vzp));
	        Wm = 1/sqrt(1 -    gxx[i3Dzm]*vxm*vxm - gyy[i3Dzm]*vym*vym - gzz[i3Dzm]*vzm*vzm 
	                      - 2*(gxy[i3Dzm]*vxm*vym + gxz[i3Dzm]*vxm*vzm + gyz[i3Dzm]*vym*vzm));
	        dWvxdz = (Wp*vxp - Wm*vxm)*idz*0.5; 
	        dWvydz = (Wp*vyp - Wm*vym)*idz*0.5; 
	        dWvzdz = (Wp*vzp - Wm*vzm)*idz*0.5; 
                dWdz   = (Wp - Wm)*idz*0.5;
	         
                dWdt = (W-W_p)/dt;
                dWvxdt = (W*vx-W_p*vx_p)/dt;
                dWvydt = (W*vy-W_p*vy_p)/dt;
                dWvzdt = (W*vz-W_p*vz_p)/dt;

	        ndW = dWdt/alp[i3D] - (betax[i3D]*dWdx + betay[i3D]*dWdy + betaz[i3D]*dWdz)/alp[i3D];
	        ndWvx = dWvxdt/alp[i3D] - (betax[i3D]*dWvxdx + betay[i3D]*dWvxdy + betaz[i3D]*dWvxdz)/alp[i3D];
	        ndWvy = dWvydt/alp[i3D] - (betax[i3D]*dWvydx + betay[i3D]*dWvydy + betaz[i3D]*dWvydz)/alp[i3D];
	        ndWvz = dWvzdt/alp[i3D] - (betax[i3D]*dWvzdx + betay[i3D]*dWvzdy + betaz[i3D]*dWvzdz)/alp[i3D];
	      } else {
	        dWvxdx = 0; dWvxdy = 0; dWvxdz = 0;
	        dWvydx = 0; dWvydy = 0; dWvydz = 0;
	        dWvzdx = 0; dWvzdy = 0; dWvzdz = 0;
	        dWvxdt = 0; dWvydt = 0; dWvzdt = 0;
	        dWdx = 0; dWdy = 0; dWdz = 0; dWdt = 0;
	        W=1; ndW = 0; ndWvx = 0; ndWvy = 0; ndWvz = 0;
	      }
	    }

	    cl.loadBackground(vx,vy,vz,gammac);
	     
	    double alph   = alp[i3D];
	    double sdetg  = sqrt(gammac.getDetg());
	    double sdetgo = sqrt(gammaco.getDetg());
            
            // TODO: combine log() via log(a)-log(b) = log(a/b)
	    double dlnadx = (log(alp[i3Dxp]) - log(alp[i3Dxm]))*idx*0.5;
	    double dlnady = (log(alp[i3Dyp]) - log(alp[i3Dym]))*idy*0.5;
	    double dlnadz = (log(alp[i3Dzp]) - log(alp[i3Dzm]))*idz*0.5;
	     
	    double dbxdx = (betax[i3Dxp] - betax[i3Dxm])*idx*0.5;
	    double dbxdy = (betax[i3Dyp] - betax[i3Dym])*idy*0.5;
	    double dbxdz = (betax[i3Dzp] - betax[i3Dzm])*idz*0.5;
	    double dbydx = (betay[i3Dxp] - betay[i3Dxm])*idx*0.5;
	    double dbydy = (betay[i3Dyp] - betay[i3Dym])*idy*0.5;
	    double dbydz = (betay[i3Dzp] - betay[i3Dzm])*idz*0.5;
	    double dbzdx = (betaz[i3Dxp] - betaz[i3Dxm])*idx*0.5;
	    double dbzdy = (betaz[i3Dyp] - betaz[i3Dym])*idy*0.5;
	    double dbzdz = (betaz[i3Dzp] - betaz[i3Dzm])*idz*0.5;
	    
	    double dgxxdx = (gxx[i3Dxp] - gxx[i3Dxm])*idx*0.5;
	    double dgxxdy = (gxx[i3Dyp] - gxx[i3Dym])*idy*0.5;
	    double dgxxdz = (gxx[i3Dzp] - gxx[i3Dzm])*idz*0.5;
	    double dgxydx = (gxy[i3Dxp] - gxy[i3Dxm])*idx*0.5;
	    double dgxydy = (gxy[i3Dyp] - gxy[i3Dym])*idy*0.5;
	    double dgxydz = (gxy[i3Dzp] - gxy[i3Dzm])*idz*0.5;
	    double dgxzdx = (gxz[i3Dxp] - gxz[i3Dxm])*idx*0.5;
	    double dgxzdy = (gxz[i3Dyp] - gxz[i3Dym])*idy*0.5;
	    double dgxzdz = (gxz[i3Dzp] - gxz[i3Dzm])*idz*0.5;
	    double dgyydx = (gyy[i3Dxp] - gyy[i3Dxm])*idx*0.5;
	    double dgyydy = (gyy[i3Dyp] - gyy[i3Dym])*idy*0.5;
	    double dgyydz = (gyy[i3Dzp] - gyy[i3Dzm])*idz*0.5;
	    double dgyzdx = (gyz[i3Dxp] - gyz[i3Dxm])*idx*0.5;
	    double dgyzdy = (gyz[i3Dyp] - gyz[i3Dym])*idy*0.5;
	    double dgyzdz = (gyz[i3Dzp] - gyz[i3Dzm])*idz*0.5;
	    double dgzzdx = (gzz[i3Dxp] - gzz[i3Dxm])*idx*0.5;
	    double dgzzdy = (gzz[i3Dyp] - gzz[i3Dym])*idy*0.5;
	    double dgzzdz = (gzz[i3Dzp] - gzz[i3Dzm])*idz*0.5;
	   
	    
            double aa, bbx, bby, bbz;
	    double cxx, cxy, cxz, cyy, cyz, czz; 
	    if (do_m1_redshift) { 
	      // Calculat some stuff to contract with the rest frame moments
	      aa   = -W*(vx*dlnadx + vy*dlnady + vz*dlnadz);
              double bbxu = ndWvx + W/alp[i3D]*(vx*dbxdx + vy*dbxdy + vz*dbxdz); 
              double bbyu = ndWvy + W/alp[i3D]*(vx*dbydx + vy*dbydy + vz*dbydz); 
              double bbzu = ndWvz + W/alp[i3D]*(vx*dbzdx + vy*dbzdy + vz*dbzdz); 
	      gammac.lowerAu(bbxu,bbyu,bbzu,&bbx,&bby,&bbz);
	      
	      cxx = 0.5*(gxx[i3D]*dWvxdx + gxy[i3D]*dWvydx + gxz[i3D]*dWvzdx)  
	           +0.5*(gxx[i3D]*dWvxdx + gxy[i3D]*dWvydx + gxz[i3D]*dWvzdx);   	    
	      cxy = 0.5*(gxx[i3D]*dWvxdy + gxy[i3D]*dWvydy + gxz[i3D]*dWvzdy)  
	           +0.5*(gxy[i3D]*dWvxdx + gyy[i3D]*dWvydx + gyz[i3D]*dWvzdx);   	    
	      cxz = 0.5*(gxx[i3D]*dWvxdz + gxy[i3D]*dWvydz + gxz[i3D]*dWvzdz)  
	           +0.5*(gxz[i3D]*dWvxdx + gyz[i3D]*dWvydx + gzz[i3D]*dWvzdx);   
	      cyy = 0.5*(gxy[i3D]*dWvxdy + gyy[i3D]*dWvydy + gyz[i3D]*dWvzdy)  
	           +0.5*(gxy[i3D]*dWvxdy + gyy[i3D]*dWvydy + gyz[i3D]*dWvzdy);   	    
	      cyz = 0.5*(gxy[i3D]*dWvxdz + gyy[i3D]*dWvydz + gyz[i3D]*dWvzdz)  
	           +0.5*(gxz[i3D]*dWvxdy + gyz[i3D]*dWvydy + gzz[i3D]*dWvzdy);   	    
	      czz = 0.5*(gxz[i3D]*dWvxdz + gyz[i3D]*dWvydz + gzz[i3D]*dWvzdz)  
	           +0.5*(gxz[i3D]*dWvxdz + gyz[i3D]*dWvydz + gzz[i3D]*dWvzdz);   	  
	          	   
	      cxx += W*0.5*(vx*dgxxdx + vy*dgxxdy + vz*dgxxdz);		 
	      cxy += W*0.5*(vx*dgxydx + vy*dgxydy + vz*dgxydz);		 
	      cxz += W*0.5*(vx*dgxzdx + vy*dgxzdy + vz*dgxzdz);		 
	      cyy += W*0.5*(vx*dgyydx + vy*dgyydy + vz*dgyydz);		 
	      cxz += W*0.5*(vx*dgyzdx + vy*dgyzdy + vz*dgyzdz);		 
	      czz += W*0.5*(vx*dgzzdx + vy*dgzzdy + vz*dgzzdz);		 
            } 
	     			 	    
	    double dtau = dts*alph;
	    dtautmp[i3D] = 0.0;
	    if (do_m1_pnu) pnu[i3D] = 0.0;
	     
            for(int isp=0;isp<nspecies;isp++) { 
	      
	      // Calculate momentum space fluxes
	      if (do_m1_redshift) {
                // this loop modifies enurhs at ieg and ieg-1 and ieg+1
                // this needs to be taken into account if OpenMP parallelizing it
	        for(int ieg=0;ieg<ngroups;ieg++) {
	          
		  int ig = isp*ngroups + ieg;
	          int index4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
		  int igl = isp*ngroups + ieg-1;
		  int index4Dl =
                    ieg>0 ? CCTK_VectGFIndex3D(cctkGH,i,j,k,igl) : -1;
		  int igu = isp*ngroups + ieg+1;
	          int index4Du =
                    ieg<ngroups-1 ? CCTK_VectGFIndex3D(cctkGH,i,j,k,igu) : -1;
	          
		  // Need un-energy integrated quantities, divide through by group width 
	          double E  = enusr[index4D]/bin_widths[ieg];
	          double Fx = fnuxsr[index4D]/bin_widths[ieg]; 
		  double Fy = fnuysr[index4D]/bin_widths[ieg]; 
		  double Fz = fnuzsr[index4D]/bin_widths[ieg];
		  
	          cl.loadRad(E,Fx,Fy,Fz);
	          cl.setClosure(0.5*(xioldzL[index4D]+xioldzR[index4D]));
		  
	          double fE,fN,fFx,fFy,fFz;
	          cl.getMomentumFluxes(dlnadx,dlnady,dlnadz,ndW,dWdx,dWdy,dWdz,aa,bbx,bby,bbz,
	                               cxx,cxy,cxz,cyy,cyz,czz,&fE,&fN,&fFx,&fFy,&fFz); 
                  
		  // Calculate the fluxes using upwinding
		  if (fE*E>0 && ieg<ngroups-1) {
			  enurhs[index4D]  -= alph*fE*bin_top[ieg];
			  enurhs[index4Du] += alph*fE*bin_top[ieg];
		  } else if (ieg>0) {
			  enurhs[index4D]  += alph*fE*bin_bottom[ieg];
			  enurhs[index4Dl] -= alph*fE*bin_bottom[ieg];
		  }
	          
		  if (do_m1_nevolve) {	
		    if (fN>0 && ieg<ngroups-1) {
		            nnurhs[index4D]  -= alph*fN*bin_top[ieg];
		            nnurhs[index4Du] += alph*fN*bin_top[ieg];
		    } else if (ieg>0) {
		            nnurhs[index4D]  += alph*fN*bin_bottom[ieg];
		            nnurhs[index4Dl] -= alph*fN*bin_bottom[ieg];
		    }
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
		
		}
	      } // do_m1_redshift
              

	      for(int ieg=0;ieg<ngroups;ieg++) {
	    
	        // Calculate closure 
	        int ig = isp*ngroups + ieg;
	        int index4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	        double E  = enusr[index4D];
	        double Fx = fnuxsr[index4D]; double Fy = fnuysr[index4D]; double Fz = fnuzsr[index4D];
	        double E0  = enu[index4D];
	        double Fx0 = fnux[index4D]; double Fy0 = fnuy[index4D]; double Fz0 = fnuz[index4D];
                double N0 = 0.0;
		double N  = 0.0;
	        
		cl.loadRad(E,Fx,Fy,Fz);
	        cl.setClosure(0.5*(xioldzL[index4D]+xioldzR[index4D]));
                
	        // Calculate explicit portion of updates
	        double fexplicitx = Fx0 + dts*fnuxrhs[index4D];
	        double fexplicity = Fy0 + dts*fnuyrhs[index4D];
	        double fexplicitz = Fz0 + dts*fnuzrhs[index4D];
	        double eexplicit  = E0  + dts*enurhs[index4D];
                
		if (do_m1_GRsource) {
	        
		  double Fdotdlna  = gammac.contractGuuAlBl(Fx,Fy,Fz,dlnadx,dlnady,dlnadz); 
	          double Fdotdbdx  = gammac.contractGuuAlBl(Fx,Fy,Fz,dbxdx,dbydx,dbzdx); 
	          double Fdotdbdy  = gammac.contractGuuAlBl(Fx,Fy,Fz,dbxdy,dbydy,dbzdy); 
	          double Fdotdbdz  = gammac.contractGuuAlBl(Fx,Fy,Fz,dbxdz,dbydz,dbzdz); 
	          
	          double PK    = cl.contractPuuTllSymmetric(kxx[i3D],kxy[i3D],kxz[i3D],kyy[i3D],kyz[i3D],kzz[i3D]); 
	          double Pdgdx = cl.contractPuuTllSymmetric(dgxxdx,dgxydx,dgxzdx,dgyydx,dgyzdx,dgzzdx); 
	          double Pdgdy = cl.contractPuuTllSymmetric(dgxxdy,dgxydy,dgxzdy,dgyydy,dgyzdy,dgzzdy); 
	          double Pdgdz = cl.contractPuuTllSymmetric(dgxxdz,dgxydz,dgxzdz,dgyydz,dgyzdz,dgzzdz); 
	        
	          fexplicitx = fexplicitx + dtau*(- E*dlnadx + Fdotdbdx/alph + Pdgdx/2.0);
	          fexplicity = fexplicity + dtau*(- E*dlnady + Fdotdbdy/alph + Pdgdy/2.0);
	          fexplicitz = fexplicitz + dtau*(- E*dlnadz + Fdotdbdz/alph + Pdgdz/2.0);
	          eexplicit  = eexplicit  + dtau*(+ PK       - Fdotdlna                 );
		}
		 
	        double nexplicit  = 0.0;
		if (do_m1_nevolve) {
		  N  = nnusr[index4D];
		  N0 = nnu[index4D];
	          nexplicit  = N0 + dts*nnurhs[index4D];
		}

    	        if (*zm1_RKstep<2) {
	                Etotold  += enupt[index4D]; 
	                Fxtotold += fnuxpt[index4D];
	        }

	        // Do implicit updates of radiation quantities 
	        
	        double J,dJ,dN,dE,dSx,dSy,dSz;
	        int idxopac = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k);
	        idxopac = index4D;
	         
	        double ab = absorb[idxopac]; 
	        double sc = scat[idxopac]; 
	        double em = emis[idxopac]*sdetg; // The densitized emissivity
                
	        // Use the source term defined in the fluid frame, which is correct 
	        dE = enupt[index4D]; 
	        cl.getSourceUpdate(eexplicit,nexplicit,fexplicitx,fexplicity,fexplicitz,ab,sc,em,
	            	         dtau,&E,&N,&Fx,&Fy,&Fz,&J,&dJ,&dN,&dE,&dSx,&dSy,&dSz);
	        if (do_m1_nevolve) {
		  nnupt[index4D]  = N;
		  dJ = dN; // Set the correct Ye source term
		}
		enupt[index4D]  = E;
	        fnuxpt[index4D] = Fx;
	        fnuypt[index4D] = Fy;
	        fnuzpt[index4D] = Fz;
	     
	        // Check for updates to fluid quantities and other things 
	        if (*zm1_RKstep<2){
	          
	          double v2 = gammac.contractGllAuBu(vx,vy,vz,vx,vy,vz);
	          double W  = 1.0/sqrt(1.0-v2);
	          
	          // Calculate the magnitude of the flux 
	          fnumag[index4D] = sqrt(gammac.contractGuuAlBl(fnuxpt[index4D],fnuypt[index4D],fnuzpt[index4D],
	                                                        fnuxpt[index4D],fnuypt[index4D],fnuzpt[index4D])); 
	          // Calculate the contribution to the neutrino pressure
		  if (do_m1_pnu) {
		    pnu[i3D] += E/3.0;
		  }
		     	
	          // Internal energy update (Conversion from MeV/fm^3 to geometrized eps)
	          if (do_m1_eps_backreact) {
	          	//eps[i3D] -= dJ/(sdetg*rho[i3D]*W); 
		    //	          	tau[i3D] -= dE;
		    dtautmp[i3D] -= dE;
	          }

		  if (do_m1_scon_backreact) {
                    double pnu_rolloff = 1.0;
		    if (do_m1_pnu) {
		       if (rho[i3D]*INV_RHO_GF>10.0*pnu_dens) {
		         pnu_rolloff = 0.0;
		       } else if (rho[i3D]*INV_RHO_GF<0.1*pnu_dens) {	  
		         pnu_rolloff = 1.0;
		       } else {	 	       
		         pnu_rolloff = 1.0 - tanh(rho[i3D]*INV_RHO_GF/pnu_dens);
	               } 
		    }		    
		    
   		    scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] -= dSx*sconfac*pnu_rolloff;
		    scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] -= dSy*sconfac*pnu_rolloff;
		    scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] -= dSz*sconfac*pnu_rolloff;
		  }

	          // Y_e update	
	          // Spectral transport case
	          if (do_m1_ye_backreact && do_opac) {
	          	if ((isp==0 && spec_idx[isp]<1) || spec_idx[isp]==1 ) { 
	          		//Y_e[i3D] -= YE_CONV_FAC*dJ/(W*neutrino_energies[sgroup+ieg]*rho[i3D]*sdetg);
	          		Y_e_con[i3D] -= YE_CONV_FAC*dJ/neutrino_energies[sgroup+ieg];
	          	} else if ((isp==1 && spec_idx[isp]<1) || spec_idx[isp]==2 ){ 
	          		//Y_e[i3D] += YE_CONV_FAC*dJ/(W*neutrino_energies[sgroup+ieg]*rho[i3D]*sdetg);
	          		Y_e_con[i3D] += YE_CONV_FAC*dJ/neutrino_energies[sgroup+ieg];
	          	}
	          }
	          
	          // Y_e update	
	          // Gray opacity case 	
	          if (do_m1_ye_backreact && rho[i3D]>1.e-20 && do_gray) {
	          	if ((isp==0 && spec_idx[isp]<1) || spec_idx[isp]==1 ) { 
	          	    double eta  = munu/temperature[i3D];
	          	    double eavg = temperature[i3D]*Fermi_3(eta)/Fermi_2(eta);  
	          	    if (eavg<8.0) eavg = 8.0;  // Put a floor on the neutrino average energy
	          	    Y_e[i3D] -= YE_CONV_FAC*dJ/(eavg*rho[i3D]*sdetg);
	          	} else if ((isp==1 && spec_idx[isp]<1) || spec_idx[isp]==2 ){ 
	          	    double eta  = -munu/temperature[i3D];
	          	    double eavg = temperature[i3D]*Fermi_3(eta)/Fermi_2(eta);  
	          	    if (eavg<8.0) eavg = 8.0;  // Put a floor on the neutrino average energy
	          	    Y_e[i3D] += YE_CONV_FAC*dJ/(eavg*rho[i3D]*sdetg);
	          	}
	          }

	        	Etotnew  += enupt[index4D];
	          Fxtotnew += fnuxpt[index4D]; 
	        }

	        //Refluxing
                if (do_m1_reflux) {             
                  // Refluxing: Integrate flux registers
                  for (int var=0; var< flux_nvars; ++var) {
                    for (int dir=0; dir<flux_ndirs; ++dir) {
                      const int n = (var * flux_nelts + ig) * flux_ndirs + dir;
                      const int flux_idx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, n);
                      register_fine_pt[flux_idx] += dts * flux[flux_idx];
                      register_coarse_pt[flux_idx] += dts * flux[flux_idx];
                      if (var==0 && i==40 && j==3 && k==3) {
                        if (zm1_verbose)
	                  CCTK_VInfo(CCTK_THORNSTRING,
                                   "Refluxing: integrating flux %g %g",
                                   (double)register_fine_pt[flux_idx],
                                   (double)register_coarse_pt[flux_idx]);
                      }
                    }
                  }
                }
	      } // Loop over energies	
	    } // Loop over species
	    
	    // apply actual update to tau (dtautmp will be 0.0 if do_m1_backreact_eps is false)
	    tau[i3D] += dtautmp[i3D];
	    // set analysis GFs
	    heatcool[i3D] = dtautmp[i3D]*idt*volume_form[i3D];
	    if(dtautmp[i3D] > 0.0) {
	      netheat[i3D] = dtautmp[i3D]*idt*volume_form[i3D];
	      netcool[i3D] = 0.0;
	    } else {
	      netheat[i3D] = 0.0;
	      netcool[i3D] = dtautmp[i3D]*idt*volume_form[i3D];
	    }


    } // Loop over zones

} // OpenMP parallel region
    
    if (*zm1_RKstep<2 && do_M1_testing) 
      
      CCTK_VInfo(CCTK_THORNSTRING,"[Conservation] Enu => %13.5e (%13.5e) Fnux => %13.5e (%13.5e)",(Etotold-Etotnew)/Etotold,Etotnew,(Fxtotold-Fxtotnew)/Fxtotold,Fxtotnew);
    
    return;

  }

}

