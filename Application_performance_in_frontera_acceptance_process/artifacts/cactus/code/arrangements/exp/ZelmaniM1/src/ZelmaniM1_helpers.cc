#include <iostream>
#include <algorithm>
#include <time.h>
#include <cstdio>
#include <util_Table.h>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Symmetry.h"
#include "ZelmaniM1_Closure.hh"
#include "ZelmaniM1_Metric.hh"
#include "carpet.hh"

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

            enue_tot[index3D] = 0.0;
            enub_tot[index3D] = 0.0;
            enux_tot[index3D] = 0.0;
    }
    
    CCTK_VInfo(CCTK_THORNSTRING,"Initialized Radiation.");
  }

  extern "C" 
  void zm1_ZeroRHS(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    if (zm1_verbose) 
      CCTK_VInfo(CCTK_THORNSTRING,"ZeroRHS    : %d   cctk_time: %8.5f",
          *zm1_RKstep,cctk_time);
    
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
#pragma omp parallel for
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=0;k<nz;k++) 
        for(int j=0;j<ny;j++)
          for(int i=0;i<nx;i++) {
            int index  = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
            enurhs[index]  = 0.0;
            fnuxrhs[index] = 0.0;
            fnuyrhs[index] = 0.0;
            fnuzrhs[index] = 0.0;
    }
  }

  extern "C" 
  void zm1_UpdateEoS(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    //DECLARE_CCTK_FUNCTIONS;

    if (zm1_verbose)
      CCTK_VInfo(CCTK_THORNSTRING,"UpdateEoS  : %d   cctk_time: %8.5f",
          *zm1_RKstep,cctk_time);
    
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
  }

  extern "C" 
  void zm1_StartRKLoop(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    if (do_m1_RK2) *zm1_RKstep  = 2;
    else *zm1_RKstep = 1;
    *zm1_RKstep_half = 1;
    *zm1_RKstep_full = 0;
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
      CCTK_VInfo(CCTK_THORNSTRING,
          "Zero velocity and curvature (cctk_time: %8.5f)",cctk_time);
    
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
      //      CCTK_VInfo(CCTK_THORNSTRING,"vel[%d] = %e => 0.0 (i=%d,j=%d,k=%d) ",idx2,vel[idx2],i,j,k);
      //vel[idx2             ] = 0.0;
      //vel[idx2 +   nx*ny*nz] = 0.0;
      //vel[idx2 + 2*nx*ny*nz] = 0.0;
    
            kxx[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kxy[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kxz[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kyy[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kyz[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
            kzz[CCTK_GFIndex3D(cctkGH,i,j,k)] = 0.0;
    }
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
          double alpha = alp[i3D];
          ThreeTensor::Metric gamma(gxx[i3D],gxy[i3D],gxz[i3D],
                                    gyy[i3D],gyz[i3D],gzz[i3D],
                                    alpha,betax[i3D],betay[i3D],betaz[i3D]);
          double sdetg = sqrt(gamma.getDetg()); 
          double vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
          double vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
          double vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
          
          if (use_zerov) {
            vx = 0.0;
            vy = 0.0;
            vz = 0.0;
          }

          Closure cl(vx,vy,vz,gamma);
      
          for (int ig=0;ig<ngroups*nspecies;ig++) {
            int i4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
#if 0
            if(isnan(enu[i4D]) || isnan(fnux[i4D])) {
              CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "rl: %d ig: %d ijk: %d,%d,%d, %15.6E %15.6E "
                "%15.6E %15.6E %15.6E", Carpet::reflevel,ig,i,j,k,
                r[i3D],enu[i4D],fnux[i4D],fnuy[i4D],fnuz[i4D]);
            } 
#endif

            // Calculate the neutrino pressure 3-tensor in the lab frame
            double Pxx,Pxy,Pxz,Pyy,Pyz,Pzz;
            cl.setDebugInfo("tmunu");
            cl.setClosure(enu[i4D],fnux[i4D],fnuy[i4D],fnuz[i4D],xioldzL[i4D]);
            cl.getPll(&Pxx,&Pxy,&Pxz,&Pyy,&Pyz,&Pzz); 
            
            eTtt[i3D] += alpha*alpha*enu[i4D]/sdetg;    
            eTtx[i3D] -= alpha*fnux[i4D]/sdetg;   
            eTty[i3D] -= alpha*fnuy[i4D]/sdetg;   
            eTtz[i3D] -= alpha*fnuz[i4D]/sdetg;   
            eTxx[i3D] += Pxx/sdetg;   
            eTxy[i3D] += Pxy/sdetg;   
            eTxz[i3D] += Pxz/sdetg;   
            eTyy[i3D] += Pyy/sdetg;   
            eTyz[i3D] += Pyz/sdetg;   
            eTzz[i3D] += Pzz/sdetg; 
          }
    }
    
    return;
  
  }

  extern "C"
  void zm1_ProbeShit(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];

#pragma omp parallel for 
    for(int k=0;k<nz;k++) 
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++) {
          int i3D = CCTK_GFIndex3D(cctkGH,i,j,k);
          for (int ig=0;ig<ngroups*nspecies;ig++) {
            int i4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
            if(isnan(enu[i4D]) || isnan(fnux[i4D]) || isnan(fnuxh[i4D]) ) {
                CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "rl: %d ig: %d ijk: %d,%d,%d, %15.6E %15.6E"
                " %15.6E %15.6E %15.6E %15.6E",Carpet::reflevel,ig,i,j,k,
                r[i3D],enu[i4D],fnux[i4D],fnuxh[i4D],fnuy[i4D],fnuz[i4D]);
            } 
          }
    }
  }
  
  extern "C"
  void zm1_RegisterVars(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int ierr=0;

    CCTK_INFO("Registering evolution variables.");
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
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
                      "ZelmaniM1::enuh",zm1_bound_type);
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, zm1_ghost, bound_table, 
                      "ZelmaniM1::fnuh",zm1_bound_type);
      if (zm1_verbose)
        CCTK_INFO("RK2 Half step boundary condition");
    } else { 
      // This should be called during POSTREGRID and POSTRESTRICT 
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
          "Problems registering boundary condition '%s' for "
          "evolution vars ZelmaniM1::enu and ZelmaniM1::fnu",zm1_bound_type);
    }
  }
}

