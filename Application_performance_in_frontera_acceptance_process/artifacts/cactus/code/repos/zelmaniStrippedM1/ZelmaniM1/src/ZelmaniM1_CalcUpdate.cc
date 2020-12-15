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
    
    double *enupt, *fnuxpt, *fnuypt, *fnuzpt;
    double *enusr, *fnuxsr, *fnuysr, *fnuzsr;
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
      CCTK_VInfo(CCTK_THORNSTRING,"CalcUpdate : %d "
        "  cctk_time: %8.5f",*zm1_RKstep,cctk_time);
    
    if (*zm1_RKstep == 2) {
      dts = 0.5*dt;
      enupt  = enuh;
      fnuxpt = fnuxh;
      fnuypt = fnuyh;
      fnuzpt = fnuzh;
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
              //int i4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
              // Shouldn't have to do this if we are setting boundary conditions
              //enuh[i4D]  = enu[i4D];
              //fnuxh[i4D] = fnux[i4D];
              //fnuyh[i4D] = fnuy[i4D];
              //fnuzh[i4D] = fnuz[i4D];
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
      dts = dt;
      enupt  = enu;
      fnuxpt = fnux;
      fnuypt = fnuy;
      fnuzpt = fnuz;
      
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
          
          const int i3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
          const int i3Dxp = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
          const int i3Dxm = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
          const int i3Dyp = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
          const int i3Dym = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
          const int i3Dzp = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
          const int i3Dzm = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
         
                
          ThreeTensor::Metric gammac(gxx[i3D],gxy[i3D],gxz[i3D],gyy[i3D],gyz[i3D],gzz[i3D]);
          ThreeTensor::Metric gammaco(gxx_p[i3D],gxy_p[i3D],gxz_p[i3D],gyy_p[i3D],gyz_p[i3D],gzz_p[i3D]);
          gammac.loadExtCurv(kxx[i3D],kxy[i3D],kxz[i3D],kyy[i3D],kyz[i3D],kzz[i3D]);
                
          double vx=0.0,vy=0.0,vz=0.0;
          if (!use_zerov) {
            vx = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
            vy = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
            vz = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
          }
    
          Closure cl(vx,vy,vz,gammac);
          cl.setDebugInfo("guenni1");
           
          const double alph   = alp[i3D];
          const double sdetg  = sqrt(gammac.getDetg());
                 
          const double dlnadx = (log(alp[i3Dxp]/alp[i3Dxm]))*idx*0.5;
          const double dlnady = (log(alp[i3Dyp]/alp[i3Dym]))*idy*0.5;
          const double dlnadz = (log(alp[i3Dzp]/alp[i3Dzm]))*idz*0.5;
            
          const double dbxdx = (betax[i3Dxp] - betax[i3Dxm])*idx*0.5;
          const double dbxdy = (betax[i3Dyp] - betax[i3Dym])*idy*0.5;
          const double dbxdz = (betax[i3Dzp] - betax[i3Dzm])*idz*0.5;
          const double dbydx = (betay[i3Dxp] - betay[i3Dxm])*idx*0.5;
          const double dbydy = (betay[i3Dyp] - betay[i3Dym])*idy*0.5;
          const double dbydz = (betay[i3Dzp] - betay[i3Dzm])*idz*0.5;
          const double dbzdx = (betaz[i3Dxp] - betaz[i3Dxm])*idx*0.5;
          const double dbzdy = (betaz[i3Dyp] - betaz[i3Dym])*idy*0.5;
          const double dbzdz = (betaz[i3Dzp] - betaz[i3Dzm])*idz*0.5;
           
          const double dgxxdx = (gxx[i3Dxp] - gxx[i3Dxm])*idx*0.5;
          const double dgxxdy = (gxx[i3Dyp] - gxx[i3Dym])*idy*0.5;
          const double dgxxdz = (gxx[i3Dzp] - gxx[i3Dzm])*idz*0.5;
          const double dgxydx = (gxy[i3Dxp] - gxy[i3Dxm])*idx*0.5;
          const double dgxydy = (gxy[i3Dyp] - gxy[i3Dym])*idy*0.5;
          const double dgxydz = (gxy[i3Dzp] - gxy[i3Dzm])*idz*0.5;
          const double dgxzdx = (gxz[i3Dxp] - gxz[i3Dxm])*idx*0.5;
          const double dgxzdy = (gxz[i3Dyp] - gxz[i3Dym])*idy*0.5;
          const double dgxzdz = (gxz[i3Dzp] - gxz[i3Dzm])*idz*0.5;
          const double dgyydx = (gyy[i3Dxp] - gyy[i3Dxm])*idx*0.5;
          const double dgyydy = (gyy[i3Dyp] - gyy[i3Dym])*idy*0.5;
          const double dgyydz = (gyy[i3Dzp] - gyy[i3Dzm])*idz*0.5;
          const double dgyzdx = (gyz[i3Dxp] - gyz[i3Dxm])*idx*0.5;
          const double dgyzdy = (gyz[i3Dyp] - gyz[i3Dym])*idy*0.5;
          const double dgyzdz = (gyz[i3Dzp] - gyz[i3Dzm])*idz*0.5;
          const double dgzzdx = (gzz[i3Dxp] - gzz[i3Dxm])*idx*0.5;
          const double dgzzdy = (gzz[i3Dyp] - gzz[i3Dym])*idy*0.5;
          const double dgzzdz = (gzz[i3Dzp] - gzz[i3Dzm])*idz*0.5;
         
          const double dtau = dts*alph;
          dtautmp[i3D] = 0.0;
          enue_tot[i3D] = 0.0; 
          enub_tot[i3D] = 0.0; 
          enux_tot[i3D] = 0.0; 
          for(int isp=0;isp<nspecies;isp++) { 
            for(int ieg=0;ieg<ngroups;ieg++) {
          
              // Calculate closure 
              const int ig = isp*ngroups + ieg;
              const int i4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);

              double E  = enusr[i4D];
              double Fx = fnuxsr[i4D]; 
              double Fy = fnuysr[i4D]; 
              double Fz = fnuzsr[i4D];
              
              const double E0  = enu[i4D];
              const double Fx0 = fnux[i4D]; 
              const double Fy0 = fnuy[i4D]; 
              const double Fz0 = fnuz[i4D];
              double N   = 0.0;
              
              cl.setClosure(E,Fx,Fy,Fz,0.5*(xioldzL[i4D]+xioldzR[i4D]));
                    
              // Calculate explicit portion of updates 
              double fexplicitx = Fx0 + dts*fnuxrhs[i4D];
              double fexplicity = Fy0 + dts*fnuyrhs[i4D];
              double fexplicitz = Fz0 + dts*fnuzrhs[i4D];
              double eexplicit  = E0  + dts*enurhs[i4D];
                    
              if (do_m1_GRsource) {
                const double Fdotdlna  = gammac.contractGuuAlBl(Fx,Fy,Fz,
                    dlnadx,dlnady,dlnadz); 
                const double Fdotdbdx  = gammac.contractGuuAlBl(Fx,Fy,Fz,
                    dbxdx,dbydx,dbzdx); 
                const double Fdotdbdy  = gammac.contractGuuAlBl(Fx,Fy,Fz,
                    dbxdy,dbydy,dbzdy); 
                const double Fdotdbdz  = gammac.contractGuuAlBl(Fx,Fy,Fz,
                    dbxdz,dbydz,dbzdz); 
                
                const double PK    = cl.contractPuuTllSymmetric(kxx[i3D],
                    kxy[i3D],kxz[i3D],kyy[i3D],kyz[i3D],kzz[i3D]); 
                const double Pdgdx = cl.contractPuuTllSymmetric(dgxxdx,dgxydx,
                    dgxzdx,dgyydx,dgyzdx,dgzzdx); 
                const double Pdgdy = cl.contractPuuTllSymmetric(dgxxdy,dgxydy,
                    dgxzdy,dgyydy,dgyzdy,dgzzdy); 
                const double Pdgdz = cl.contractPuuTllSymmetric(dgxxdz,dgxydz,
                    dgxzdz,dgyydz,dgyzdz,dgzzdz); 
                
                fexplicitx += dtau*(- E*dlnadx + Fdotdbdx/alph + Pdgdx/2.0);
                fexplicity += dtau*(- E*dlnady + Fdotdbdy/alph + Pdgdy/2.0);
                fexplicitz += dtau*(- E*dlnadz + Fdotdbdz/alph + Pdgdz/2.0);
                eexplicit  += dtau*(+ PK       - Fdotdlna                 );
              }
               
    
              if (*zm1_RKstep<2) {
                Etotold  += enupt[i4D]; 
                Fxtotold += fnuxpt[i4D];
              }
    
              // Do implicit updates of radiation quantities 
              double J,dJ,dN,dE,dSx,dSy,dSz;
              double ab = absorb[i4D]; 
              double sc = scat[i4D]; 
              double em = emis[i4D]*sdetg; // The densitized emissivity
                    
              // Use the source term defined in the fluid frame, which is correct 
              dE = enupt[i4D]; 
              double nexplicit  = 0.0;
              cl.getSourceUpdate(eexplicit,nexplicit,fexplicitx,fexplicity,fexplicitz,
                  ab,sc,em,dtau,&E,&N,&Fx,&Fy,&Fz,&J,&dJ,&dN,&dE,&dSx,&dSy,&dSz);
              //#pragma omp critical
              //cout << emis[i4D] << " "  << em/ab << " " << sqrt(x[i3D]*x[i3D] 
              //    + y[i3D]*y[i3D] + z[i3D]*z[i3D]) << "\n";
              
              //E = (eexplicit + em*dtau)/(1.0 + ab*dtau);
              //Fx = (fexplicitx)/(1.0 + (ab+sc)*dtau);
              //Fy = (fexplicity)/(1.0 + (ab+sc)*dtau);
              //Fz = (fexplicitz)/(1.0 + (ab+sc)*dtau);
              enupt[i4D]  = E;
              fnuxpt[i4D] = Fx;
              fnuypt[i4D] = Fy;
              fnuzpt[i4D] = Fz;
           
              // Check for updates to fluid quantities and other things 
              if (*zm1_RKstep < 2){
                if (isp==0) enue_tot[i3D] += E;
                if (isp==1) enub_tot[i3D] += E;
                if (isp==2) enux_tot[i3D] += E;
                if (isnan(dE)) dE = 0.0; 
                if (isnan(dJ)) dJ = 0.0; 
                if (isnan(dSx)) dSx = 0.0; 
                if (isnan(dSy)) dSy = 0.0; 
                if (isnan(dSy)) dSz = 0.0; 
    
                // Calculate the magnitude of the flux 
                fnumag[i4D] = sqrt(gammac.contractGuuAlBl(fnuxpt[i4D],fnuypt[i4D],fnuzpt[i4D],
                                                              fnuxpt[i4D],fnuypt[i4D],fnuzpt[i4D])); 
              
                // Internal energy update (Conversion from MeV/fm^3 to geometrized eps)
                if (do_m1_eps_backreact)
                  dtautmp[i3D] -= dE;
              
                if (do_m1_scon_backreact) {
                  scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] -= dSx*sconfac;
                  scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] -= dSy*sconfac;
                  scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] -= dSz*sconfac;
                }
    
                // Y_e update 
                // Spectral transport case
                if (do_m1_ye_backreact && do_opac) {
                  if ((isp==0 && spec_idx[isp]<1) || spec_idx[isp]==1 ) { 
                    Y_e_con[i3D] -= YE_CONV_FAC*dJ/neutrino_energies[sgroup+ieg];
                  } else if ((isp==1 && spec_idx[isp]<1) || spec_idx[isp]==2 ){ 
                    Y_e_con[i3D] += YE_CONV_FAC*dJ/neutrino_energies[sgroup+ieg];
                  }
                }
                
                Etotnew  += enupt[i4D];
                Fxtotnew += fnuxpt[i4D]; 
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
                      if (zm1_verbose) CCTK_VInfo(CCTK_THORNSTRING,
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
      CCTK_VInfo(CCTK_THORNSTRING,"[Conservation] Enu => %13.5e (%13.5e)"
        " Fnux => %13.5e (%13.5e)",(Etotold-Etotnew)/Etotold,Etotnew,
        (Fxtotold-Fxtotnew)/Fxtotold,Fxtotnew);
  }
}

