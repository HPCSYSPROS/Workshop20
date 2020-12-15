#include <stdio.h>
#include <cassert>
#include <cmath>
#include <mpi.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Functions.h>
#include <cctk_Faces.h>
#include <cctk_Schedule.h>
#include <Symmetry.h>
#include <time.h>

#include "carpet.hh"

extern "C" {
  void PNSHelper_MapSpacetime(CCTK_ARGUMENTS);
  void PNSHelper_ChangeStatus(CCTK_ARGUMENTS);
}

void PNSHelper_ChangeStatus(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *have_interp_data = 0;
  
}

void PNSHelper_LinExt(int npoints, double times, double* data, double* out);
void PNSHelper_QuadExt(int npoints, double times, double* data, double* out);

void PNSHelper_LinExt(int npoints, double* times, double* data, double time, 
		      double* out) {
  double lam0 = (time-times[1])/(times[0]-times[1]);
  double lam1 = (time-times[0])/(times[1]-times[0]);
  for(int i=0;i<npoints;i++) {
    int iind2D_0 = i + (npoints*0);
    int iind2D_1 = i + (npoints*1);
    out[i] = data[iind2D_0]*lam0 + data[iind2D_1]*lam1;
  }
  return;
}

void PNSHelper_QuadExt(int npoints, double* times, double* data, double time, 
		      double* out) {
  double lam0 = (time-times[2])/(times[0]-times[2])*(time-times[1])/(times[0]-times[1]);
  double lam1 = (time-times[2])/(times[1]-times[2])*(time-times[0])/(times[1]-times[0]);
  double lam2 = (time-times[1])/(times[2]-times[1])*(time-times[0])/(times[2]-times[0]);
  for(int i=0;i<npoints;i++) {
    int iind2D_0 = i + (npoints*0);
    int iind2D_1 = i + (npoints*1);
    int iind2D_2 = i + (npoints*2);
    out[i] = data[iind2D_0]*lam0 + data[iind2D_1]*lam1 + data[iind2D_2]*lam2;
  }
  return;
}

void PNSHelper_MapSpacetime(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace Carpet;

  int Num_1D = nrad+nrad_outer;
  int nouter = 500;
  int np = nrad + 3 + nouter;
  double dr = *drad;

  // metric and lapse extrapolation
  double* ex_alp = (double*) malloc(sizeof(double)*Num_1D);
  double* ex_psi = (double*) malloc(sizeof(double)*Num_1D);
  double ex_mgrav;

  int rl = trunc(log10((1.0*cctk_levfac[0]))/log10(2.0));

  int const do_every
    = ipow(mgfact, mglevel) * (maxtimereflevelfact / timereffacts.AT(rl));

  const cFunctionData * funcdata =
    CCTK_ScheduleQueryCurrentFunction(cctkGH);

  //  if(CCTK_Equals(funcdata->where,"PNSHelper_Group")) {
  //  //    if( (cctk_iteration-1) % do_every != 0) return;
  //  if( (cctk_iteration-1) % do_every != 0) return;
  // }

    if(CCTK_Equals(funcdata->where,"MoL_PostStep") && cctk_iteration < 1) {
    CCTK_Info(CCTK_THORNSTRING,"Bailing out -- don't want to be called here at Initial Data");
    return;
  }

  if(update_GR_every < 0) return;

  //fprintf(stdout,"I am %s::%s being called from %s; reflevel: %d \n",
  //        funcdata->thorn, funcdata->routine,
  //        funcdata->where, rl);

  //CCTK_VInfo(CCTK_THORNSTRING,"MapMetric: cctk_time: %8.5f",cctk_time);
  
  double timec = cctk_time - CCTK_DELTA_TIME;

  if (!*have_interp_data) {

    
    //void PNSHelper_LinExt(int npoints, double* times, double* data, double time, 
    //	      double* out) {

    if(extrapolation_order==0 || cctk_iteration < 16) {
      if (verbose)
      fprintf(stdout,"Using stale data data! time: %8.5f rl: %d\n",cctk_time,rl);

      ex_mgrav = *pns_mgrav;
      for(int i=0;i<Num_1D;i++) {
	ex_alp[i] = pns_alp[i];
	ex_psi[i] = pns_psi[i];
      }    

    } else if(extrapolation_order == 1 && cctk_iteration >= 16) {

      if (verbose)
      fprintf(stdout,
	      "extrapolating metric: rl: %d extrap_time: %6.2f  metric_time: %6.2f  metric_time_old: %6.2f\n",
	      rl, timec, metric_times[0], metric_times[1]);

      PNSHelper_LinExt(1,      metric_times, pns_mgrav_store, timec, &ex_mgrav); 
      PNSHelper_LinExt(Num_1D, metric_times, pns_psi_store,   timec, ex_psi); 
      PNSHelper_LinExt(Num_1D, metric_times, pns_alp_store,   timec, ex_alp); 

    } else if (extrapolation_order==2 && cctk_iteration >= 16) {
      
      if (verbose)
      fprintf(stdout,
	      "Quad extrapolating metric: rl: %d extrap_time: %6.2f  metric_times[0]: %6.2f  metric_times[1]: %6.2f metric_times[2]: %6.2f\n",
	      rl, timec, metric_times[0], metric_times[1], metric_times[2]);
      
      PNSHelper_QuadExt(1,      metric_times, pns_mgrav_store, timec, &ex_mgrav); 
      PNSHelper_QuadExt(Num_1D, metric_times, pns_psi_store,   timec, ex_psi); 
      PNSHelper_QuadExt(Num_1D, metric_times, pns_alp_store,   timec, ex_alp); 

      //CCTK_WARN(0,"This metric extrapolation order is not implemented!");

    } else {

      CCTK_WARN(0,"This metric extrapolation order is not implemented!");

    }

#if 0
    ex_mgrav = (*pns_mgrav - *pns_mgrav_old) / (*metric_time - *metric_time_old) 
      * (cctk_time - *metric_time) + *pns_mgrav;
    for(int i=0;i<Num_1D;i++) {

      ex_alp[i] = (pns_alp[i] - pns_alp_old[i]) / (*metric_time - *metric_time_old) 
	* (cctk_time - *metric_time) + pns_alp[i];
      ex_psi[i] = (pns_psi[i] - pns_psi_old[i]) / (*metric_time - *metric_time_old) 
	* (cctk_time - *metric_time) + pns_psi[i];
    }
#endif
  } else {
    if (verbose)
    fprintf(stdout,"We have fresh data! time: %8.5f rl: %d\n",cctk_time,rl);
    ex_mgrav = *pns_mgrav;
    for(int i=0;i<Num_1D;i++) {
      ex_alp[i] = pns_alp[i];
      ex_psi[i] = pns_psi[i];
    }    
  }

  if (verbose)
    fprintf(stdout,"ex_alp[0] = %15.6E    ex_psi[0] = %15.6E\n",ex_alp[0],ex_psi[0]);

  //  fprintf(stdout,"rl: %d ex_psi[0] = %15.6E  ex_alp[0] = %15.6E  ex_mgrav = %15.6E\n",
  //	  rl, ex_psi[0],ex_alp[0],ex_mgrav);


  // first deal with inner and outer boundary; we just slap on
  // some zones left and right
  double* myrad = (double*) malloc(sizeof(double)*np);
  double* mypsi4 = (double*) malloc(sizeof(double)*np);
  double* myalp = (double*) malloc(sizeof(double)*np);

  // radius
  myrad[0] = -rad[2];
  myrad[1] = -rad[1];
  myrad[2] = -rad[0];
  for(int i=3;i<nrad+3;i++) {
    myrad[i] = rad[i-3];
  }
  
  for(int i=nrad+3;i<np;i++) {
    myrad[i] = myrad[i-1] + dr;
  }

  // psi
  mypsi4[0] = pow(ex_psi[2],4);
  mypsi4[1] = pow(ex_psi[1],4);
  mypsi4[2] = pow(ex_psi[0],4);
  for(int i=3;i<nrad+3;i++) {
    mypsi4[i] = pow(ex_psi[i-3],4);
  }
  
  for(int i=nrad+3;i<np;i++) {
    mypsi4[i] = pow(1.0 + ex_mgrav / 2.0 / (myrad[i]),4.0);
  }

  // alp
  myalp[0] = (ex_alp[2]);
  myalp[1] = (ex_alp[1]);
  myalp[2] = (ex_alp[0]);
  for(int i=3;i<nrad+3;i++) {
    myalp[i] = (ex_alp[i-3]);
  }
  
  for(int i=nrad+3;i<np;i++) {
    myalp[i] = (1.0 - ex_mgrav / 2.0 / (myrad[i])) / (1.0 + ex_mgrav / 2.0 / (myrad[i]));
  }

  // preparing for interpolation
  int n = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  double rmin = myrad[0];
  double rmax = myrad[np-1];

#pragma omp parallel for
  for(int i=0;i<n;i++) {
    
    if(r[i] >= rmax) {
      #pragma omp critical
      CCTK_VWarn(1,__LINE__, __FILE__, CCTK_THORNSTRING,
		 "r=%15.6E > r_max=%15.6E in map_spacetime",r[i],rmax);
    }
    
    int ip = trunc ((r[i] - rmin - 1.0e-10) / dr) + 1;
    int im = ip - 1;

    //    #pragma omp critical
    //fprintf(stdout,"%8d %d %d %15.6E %15.6E %15.6E\n",i,im,ip,myrad[im],r[i],myrad[ip]);

    gxx[i] = (mypsi4[ip] - mypsi4[im]) / (myrad[ip] - myrad[im]) 
      * (r[i]-myrad[im])  + mypsi4[im];
    gyy[i] = gxx[i];
    gzz[i] = gxx[i];
    gxy[i] = 0.0e0;
    gxz[i] = 0.0e0;
    gyz[i] = 0.0e0;

    alp[i] = (myalp[ip] - myalp[im]) / (myrad[ip] - myrad[im])	
      * (r[i]-myrad[im])  + myalp[im];

    betax[i] = 0.0e0;
    betay[i] = 0.0e0;
    betaz[i] = 0.0e0;

  }

  free(ex_alp);
  free(ex_psi);
  free(myrad);
  free(mypsi4);
  free(myalp);
}

