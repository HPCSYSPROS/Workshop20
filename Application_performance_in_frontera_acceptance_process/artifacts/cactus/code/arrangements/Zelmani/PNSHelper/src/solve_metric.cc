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
#include <Symmetry.h>
#include <time.h>

#include <gsl/gsl_linalg.h>

extern "C"
CCTK_FCALL void CCTK_FNAME(dgtsv)(int *n, int *nrhs, double *dl,                              
        double *d__, double *du, double *b, int *ldb, int *info);


#define IND2D(n,i,j) ((n)*(i) + (j))

extern "C" {
  void PNSHelper_SolveMetric(CCTK_ARGUMENTS);
}

double pns_eval_mgrav(int n, double* rad, double* psi, double* rho);

void pns_eval_F(int n, double* psi, double*rho, double* rad,
            double mout, double dr, double* F);

void pns_eval_F_alp_psi(int n, double* alp, double* psi, double*rho, double* rad,
            double mout, double dr, double* F);

void pns_eval_jacobian_tridiag(int n, int rho_scale, double* psi, double* rho, double* rad, double dr, 
                           double* ld, double* d, double* ud);

void pns_eval_jacobian_alp_psi_tridiag(int n, double* alp, double* psi, double* rho, double* rad, double dr, 
                                   double* ld, double* d, double* ud);

void PNSHelper_SolveMetric(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!*have_interp_data) return;
  
  if (verbose) 
    CCTK_VInfo(CCTK_THORNSTRING,"SolveMetric: %d  cctk_time: %8.5f",*have_interp_data,cctk_time);

  int Num_1D = nrad+nrad_outer;

  for(int i=0;i<Num_1D;i++) {
    if (CCTK_EQUALS(collect,"interpolate")) { 
      pns_rho_star[i] = (pns_av_rho[i] + pns_av_rho[i]*pns_av_eps[i] 
        	      + pns_av_press[i])*pns_av_w_lorentz[i]*pns_av_w_lorentz[i]
                      - pns_av_press[i];
    } else {  
      pns_rho_star[i] = pns_rho_psi[i];
    }
    pns_real_rho_star[i] = pns_rho_star[i];
    pns_rho_star[i] = pns_real_rho_star[i]*pow(pns_psi[i],6);
  }


  double err = 1.0;
  double tol = 1.0e-10;
  int i = 0;
  int imax = 50;
  double*  F = (double*)(malloc(sizeof(double)*Num_1D));
  double* ld = (double*)(malloc(sizeof(double)*Num_1D));
  double*  d = (double*)(malloc(sizeof(double)*Num_1D));
  double* ud = (double*)(malloc(sizeof(double)*Num_1D));

  // find psi via NR iteration
  while(err > tol && i < imax) {
# if 1
    if (verbose)
    printf("step: %d\n",i);
#endif
    i++;

    // PNS gravitational mass for outer boundary condition
    // This is incorrect because it assumes the PNS gravitational 
    // from the last step.  We can do it correctly by just 
    // reformulating the equations in terms of the enclosed mass
    // and changing the boundary condition.  Probably does not 
    // matter.
    *pns_mgrav = pns_eval_mgrav(Num_1D,rad,pns_psi,pns_rho_star);

    // Calculate the function to be zeroed
    pns_eval_F(Num_1D,pns_psi,pns_rho_star,rad,*pns_mgrav,*drad,F);
    
    // Calculate the Jacobian of F (only non-zero elements)
    // First formulation does not garantee convergence, but 
    // works better in practice.  Second one is guaranteed to 
    // converge.
    if (i<10)
      pns_eval_jacobian_tridiag(Num_1D,0,pns_psi,pns_rho_star,rad,*drad,ld,d,ud);
    else
      pns_eval_jacobian_tridiag(Num_1D,1,pns_psi,pns_rho_star,rad,*drad,ld,d,ud);
    
    // Solve linear equation using LAPACK tri-diagonal solver
    int nrhs[1]; nrhs[0] = 1;
    int nn[1];   nn[0]   = Num_1D;
    int ldF[1];  ldF[0]  = Num_1D;
    int info[1];
    CCTK_FNAME(dgtsv)( nn, nrhs, ld, d, ud, F, ldF, info );
    
    // Apply corrections to solution vector and check dx error
    err = 0.0;
    for (int j=0;j<Num_1D;j++){ 
      if( fabs(F[j]/pns_psi[j]) > err ) err = fabs(F[j]/pns_psi[j]);
      pns_psi[j] = pns_psi[j] - F[j];
      pns_rho_star[j] = pns_real_rho_star[j]*pow(pns_psi[j],6);
    }
    	     
#if 1
    if (verbose)
    printf("err: %15.6E\n",err);
#endif

  }

  // find alp from linear constraint equation
  
  // setup new rho for alp*psi calculation:
  for(int j=0; j<Num_1D;j++) {
    if (CCTK_EQUALS(collect,"interpolate")) { 
      pns_rho_star[j] = ( pns_av_rho[j] + pns_av_rho[j]*pns_av_eps[j] + pns_av_press[j] ) *
                        (3.0*pns_av_w_lorentz[j]*pns_av_w_lorentz[j] - 2.0) + 5.0*pns_av_press[j];
    } else { 
      pns_rho_star[j] = pns_rho_alp[j];
    }
  } 

  // Evaluate alpha constraint function => F
  pns_eval_F_alp_psi(Num_1D,pns_alp,pns_psi,pns_rho_star,rad,*pns_mgrav,*drad,F);
    
  // Evaluate jacobian of constraint function => ld,d,ud (diagonals of Jacobian)
  pns_eval_jacobian_alp_psi_tridiag(Num_1D,pns_alp,pns_psi,pns_rho_star,rad,*drad,ld,d,ud);

  // Solve linear equation using LAPACK tri-diagonal solver
  int nrhs[1]; nrhs[0] = 1;
  int nn[1];     nn[0] = Num_1D;
  int ldF[1];   ldF[0] = Num_1D;
  int info[1];
  CCTK_FNAME(dgtsv)( nn, nrhs, ld, d, ud, F, ldF, info );
    
  // Apply corrections to solution vector
  // This equation is linear, so machine precision
  // should be obtained on the first step
  for (int j=0;j<Num_1D;j++) 
      pns_alp[j] = pns_alp[j] - F[j]/pns_psi[j];
  
  // Clean up memory
  free(F);
  free(ud);
  free(ld);
  free(d);

  // copy old spacetime data: shuffle timelevels
  metric_times[2] = metric_times[1];
  metric_times[1] = metric_times[0];
  metric_times[0] = cctk_time;
  pns_mgrav_store[2] = pns_mgrav_store[1];
  pns_mgrav_store[1] = pns_mgrav_store[0];
  pns_mgrav_store[0] = *pns_mgrav;
  for(int i=0;i<nrad+nrad_outer;i++) {
    int j_0 = 0;
    int j_1 = 1;
    int j_2 = 2;
    int iind2D_0 = i + (nrad+nrad_outer)*j_0; 
    int iind2D_1 = i + (nrad+nrad_outer)*j_1; 
    int iind2D_2 = i + (nrad+nrad_outer)*j_2; 

    pns_psi_store[iind2D_2] = pns_psi_store[iind2D_1];  
    pns_psi_store[iind2D_1] = pns_psi_store[iind2D_0];  
    pns_psi_store[iind2D_0] = pns_psi[i];

    pns_alp_store[iind2D_2] = pns_alp_store[iind2D_1];  
    pns_alp_store[iind2D_1] = pns_alp_store[iind2D_0];  
    pns_alp_store[iind2D_0] = pns_alp[i];

  }

}

double pns_eval_mgrav(int n, double* rad, double* psi, double* rho) {
  double m=0.0;
  double vol;
  double dr = rad[1]-rad[0];
  vol = 4.0/3.0 * M_PI * pow((rad[1] - dr/2.0),3);
  m = rho[0] * vol / psi[0];
  for(int i=0;i<n;i++) {
    vol = 4.0/3.0 * M_PI *
      ( pow(rad[i]+dr/2.0,3) - pow(rad[i]-dr/2.0,3) );
    m = m + vol*rho[i]/psi[i];
  }

#if 0
  printf("mass: %15.6E\n",m);
#endif

  return m;
}

double pns_linterp(double x1,double x2,double f1,double f2, double x) {

  return (f2-f1)/(x2-x1) * (x-x1) + f1;

}

double pns_map(int n, double* prad, double* pq, double rad) {

  double q = 0.0;
  int i = 0;

  if(rad < prad[0]) {
    q = pns_linterp(prad[0],prad[1],pq[0],pq[1],rad);
  } else {
    while(rad > prad[i] && i < n) {
      i++;
    }
    q = pns_linterp(prad[i],prad[i-1],pq[i],pq[i-1],rad);
  }

  return q;
}


void pns_eval_F(int n, double* psi, double*rho, double* rad,
            double mout, double dr, double* F) {
  double dri = 1.0/dr;
  double dri2 = dri*dri;
  double psiout;
  int il = n-1;
  // inner boundary
  F[0] = (psi[1]-psi[0])*dri2 +
    (psi[1]-psi[0]) * dri / rad[0] +
    2.0*M_PI * rho[0] / psi[0];

  for(int i=1;i<(n-1);i++) {
    F[i] = (psi[i-1]-2.0*psi[i]+psi[i+1])*dri2 +
      (psi[i+1]-psi[i-1])*dri/rad[i] +
      2.0*M_PI*rho[i] / psi[i];
  }
  psiout = 1.0 + mout / 2.0 / (rad[il]+dr);
  F[il] = (psi[il-1]- 2*psi[il] + psiout) * dri2 +
    (psiout - psi[il-1])/rad[il] * dri +
    2.0*M_PI*rho[il] / psi[il];

#if 0
  printf("psiout: %15.6E %15.6E\n",psiout,psi[il]);
#endif

}

void pns_eval_jacobian_tridiag(int n, int rho_scale, double* psi, double* rho, double* rad, double dr, 
                           double* ld, double* d, double* ud) {
  // Fill the lower, middle, and upper diagonals of the Jacobian
  // if rho_scale == 1, then the rho_adm*psi^6 is assumed constant, otherwise 
  // rho_adm is assumed constant over the iteration
  // ld and ud are vectors of length n-1
  // d is a vector of length n

  double dri = 1.0/dr;
  double dri2 = dri*dri;
  
  //Boundary
  if (rho_scale==1) d[0] = -1.0 * dri2 - dri / rad[0] - 2.0 *M_PI / (psi[0]*psi[0]) *rho[0];
  else              d[0] = -1.0 * dri2 - dri / rad[0] + 10.0*M_PI / (psi[0]*psi[0]) *rho[0];

  // Non-Boundary
  for(int k=0;k<n-1;k++) {
    ld[k  ] = dri2 - dri / rad[k+1];	  
    ud[k  ] = dri2 + dri / rad[k];	  
    if (rho_scale==1) d[k+1] = - 2.0 * dri2 - 2.0 *M_PI / (psi[k+1]*psi[k+1]) * rho[k+1];
    else              d[k+1] = - 2.0 * dri2 + 10.0*M_PI / (psi[k+1]*psi[k+1]) * rho[k+1];
  }

} // pns_eval_jacobian_tridiag

void pns_eval_F_alp_psi(int n, double* alp, double* psi, double*rho, double* rad,
            double mout, double dr, double* F) {
  double dri = 1.0/dr;
  double dri2 = dri*dri;
  double psiout, alpout;
  int il = n-1;
  // inner boundary
  F[0] = (alp[1]*psi[1]-alp[0]*psi[0])*dri2 +
    (alp[1]*psi[1]-alp[0]*psi[0]) * dri / rad[0] -
    2.0*M_PI * rho[0] * alp[0]*pow(psi[0],5);

  for(int i=1;i<(n-1);i++) {
    F[i] = (alp[i-1]*psi[i-1]-2.0*alp[i]*psi[i]+alp[i+1]*psi[i+1])*dri2 +
      (alp[i+1]*psi[i+1]-alp[i-1]*psi[i-1])*dri/rad[i] -
      2.0*M_PI*rho[i]*alp[i]*pow(psi[i],5);
  }
  psiout = 1.0 + mout / 2.0 / (rad[il]+dr);
  alpout = (1.0 - mout / 2.0 / (rad[il]+dr)) / (1.0 + mout / 2.0 / (rad[il]+dr));;
  F[il] = (alp[il-1]*psi[il-1]- 2*alp[il]*psi[il] + alpout*psiout) * dri2 +
    (alpout*psiout - alp[il-1]*psi[il-1])/rad[il] * dri -
    2.0*M_PI*rho[il]*alp[il]*pow(psi[il],5);

} //pns_eval_F_alp_psi

void pns_eval_jacobian_alp_psi_tridiag(int n, double* alp, double* psi, double* rho, double* rad, double dr, 
                                   double* ld, double* d, double* ud) {
  // Fill the lower, middle, and upper diagonals of the Jacobian
  // ld and ud are vectors of length n-1
  // d is a vector of length n

  double dri = 1.0/dr;
  double dri2 = dri*dri;
  
  //Boundary
  d[0] = - 1.0 * dri2 - dri / rad[0] - 2.0*M_PI * pow(psi[0],4) * rho[0];
  
  // Non-Boundary
  for(int k=0;k<n-1;k++) {
    ld[k  ] = dri2 - dri / rad[k+1];	  
    ud[k  ] = dri2 + dri / rad[k];	  
    d[k+1] = - 2.0 * dri2 - 2.0 * M_PI * pow(psi[k+1],4)  * rho[k+1];
  }

} // pns_eval_jacobian_alp_psi_tridiag
