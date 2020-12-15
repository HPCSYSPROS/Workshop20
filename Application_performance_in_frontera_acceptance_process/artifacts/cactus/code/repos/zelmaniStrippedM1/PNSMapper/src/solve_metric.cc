#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "ZelmaniM1.hh"

#define IND2D(n,i,j) (n*i + j)

#define PRESS_NU_CONSTANT 3.52127727e24
#define PRESS_GF 1.80123683248503e-39
#define PI 3.14159265358979e0
#define PI4 97.4090910340024e0
#define PI2 9.86960440108936e0
#define MEV_TO_ERG 1.60217733e-6


extern "C"
CCTK_FCALL void CCTK_FNAME(dgtsv)(int *n, int *nrhs, double *dl,
				  double *d__, double *du, double *b, int *ldb, int *info);
extern "C" {
  void PNSMapper_SolveMetric(CCTK_ARGUMENTS);
}

namespace PNSMapper {
  
  double get_edens(const double rho, const double temp, const double ye, const double munu);
  double pns_linterp(double x1,double x2,double f1,double f2,double x);
  double pns_map(int n, double* prad, double* pq, double rad);
  double pns_eval_mgrav(int n, double* rad, double* psi, double* rho);
  void pns_eval_F(int n, double* psi, double*rho, double* rad,
		double mout, double dr, double* F);
  void pns_eval_F_alp_psi(int n, double* alp, double* psi, double*rho, double* rad,
			  double mout, double dr, double* F);
  
  void pns_eval_jacobian_tridiag(int n, int rho_scale, double* psi, double* rho, double* rad, double dr,
				 double* ld, double* d, double* ud);

  void pns_eval_jacobian_alp_psi_tridiag(int n, double* alp, double* psi, double* rho, double* rad, double dr,
					 double* ld, double* d, double* ud);


}


using namespace std;

void PNSMapper_SolveMetric(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  using namespace PNSMapper;

  // setup grid
  double dr = Rad_1D / Num_1D;
  for(int i=0;i<Num_1D;i++) {
    rad1d[i] = dr*i + 0.5*dr;
  }


  // map profile
  for(int i=0;i<Num_1D;i++) {
    rho1d[i] = pns_map(*pnszones,pradius,prho,rad1d[i]);
    temp1d[i] = pns_map(*pnszones,pradius,ptemp,rad1d[i]);
    ye1d[i] = pns_map(*pnszones,pradius,pye,rad1d[i]);
    psi1d[i] = 1.0e0;
    vel1d[i] = pns_map(*pnszones,pradius,pvel,rad1d[i]);
    alp1d[i] = 1.0e0;
  
    int keytemp = 1;
    int keyerr = 0;
    int anyerr = 0;
    double xdummy;
    double prec = 1.0e-10;
    int eoskey = 4;
    int n = 1;

    // get pressure and internal energy
    EOS_Omni_short(eoskey,
		   keytemp,
		   prec,
		   n,
		   (const double*)(&rho1d[i]),
		   &eps1d[i],
		   &temp1d[i],
		   (const double*)&ye1d[i],
		   &press1d[i],
		   &xdummy, // entropy
		   &xdummy, // cs2
		   &xdummy, // dedt
		   &xdummy, // dpderho
		   &xdummy, // dpdrhoe
		   &munu1d[i],
		   &keyerr,
		   &anyerr);

    if(anyerr) {
      fprintf(stderr,"keyerr: %d\n",keyerr);
      CCTK_WARN(0,"EOS error in PNS Mapper!!! Is the hot EOS active in EOS_Omni?");
    }

    wlorentz1d[i] = 1.0 / sqrt((1.0 - vel1d[i]*vel1d[i]));
    //    fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,rho1d[i],temp1d[i],eps1d[i],press1d[i]);

    if(isinf(rho1d[i]) || isnan(rho1d[i]) || isnan(eps1d[i]) || isnan(vel1d[i]) || isnan(wlorentz1d[i]) || isnan(press1d[i]) ) {
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E %15.6E %15.6E %15.6E\n",i,rho1d[i],temp1d[i],ye1d[i],eps1d[i],press1d[i],wlorentz1d[i],vel1d[i]);
      CCTK_WARN(0,"NAN in profile!");      
    }

    double nuedens = get_edens(rho1d[i],temp1d[i],ye1d[i],munu1d[i]);
    double nupress = nuedens/3.0;
    rhostar1d[i] = (rho1d[i] + rho1d[i]*eps1d[i] + press1d[i] + nupress + nuedens)*wlorentz1d[i]*wlorentz1d[i] 
      - press1d[i] - nupress;
    realrhostar1d[i] = rhostar1d[i];
    // ugly hack to test some stuff 
    // *(1.0+0.01*(exp(-rad1d[i]/20.0e0)));
    rhostar1d[i] = realrhostar1d[i]*pow(psi1d[i],6);
  }

  double err = 1.0;
  double tol = 1.0e-10;
  int i = 0;
  int imax = 50;
  double mgrav1d = 0.0;
  double*  F = (double*)(malloc(sizeof(double)*Num_1D));
  double* ld = (double*)(malloc(sizeof(double)*Num_1D));
  double*  d = (double*)(malloc(sizeof(double)*Num_1D));
  double* ud = (double*)(malloc(sizeof(double)*Num_1D));



  // find psi
  while(err > tol && i < imax) {
    printf("step: %d\n",i);
    i++;

    mgrav1d = pns_eval_mgrav(Num_1D,rad1d,psi1d,rhostar1d);

    pns_eval_F(Num_1D,psi1d,rhostar1d,rad1d,mgrav1d,dr,F);
    if (i<10)
      pns_eval_jacobian_tridiag(Num_1D,0,psi1d,rhostar1d,rad1d,dr,ld,d,ud);
    else
      pns_eval_jacobian_tridiag(Num_1D,1,psi1d,rhostar1d,rad1d,dr,ld,d,ud);

    // Solve linear equation using LAPACK tri-diagonal solver
    int nrhs[1]; nrhs[0] = 1;
    int nn[1];   nn[0]   = Num_1D;
    int ldF[1];  ldF[0]  = Num_1D;
    int info[1];
    CCTK_FNAME(dgtsv)( nn, nrhs, ld, d, ud, F, ldF, info );


    err = 0.0;
    for(int j=0;j<Num_1D;j++) {
      if( fabs(F[j]/psi1d[j]) > err )
        err = fabs(F[j]/psi1d[j]);

      psi1d[j] = psi1d[j] - F[j];
      rhostar1d[j] = realrhostar1d[j]*pow(psi1d[j],6);
    }


    printf("err: %15.6E\n",err);

  }

  // find alp
  err = 1.0;
  tol = 1.0e-10;
  i = 0;
  imax = 50;
  // setup new rho for alp*psi calculation:
  for(int j=0; j<Num_1D;j++)
    rhostar1d[j] = ( rho1d[j] + rho1d[j]*eps1d[j] + press1d[j] ) *
      (3.0*wlorentz1d[j]*wlorentz1d[j] - 2.0) + 5.0*press1d[j];

  pns_eval_F_alp_psi(Num_1D,alp1d,psi1d,rhostar1d,rad1d,mgrav1d,dr,F);
  pns_eval_jacobian_alp_psi_tridiag(Num_1D,alp1d,psi1d,rhostar1d,rad1d,dr,ld,d,ud);

  // Solve linear equation using LAPACK tri-diagonal solver
  int nrhs[1]; nrhs[0] = 1;
  int nn[1];     nn[0] = Num_1D;
  int ldF[1];   ldF[0] = Num_1D;
  int info[1];
  CCTK_FNAME(dgtsv)( nn, nrhs, ld, d, ud, F, ldF, info );
  for (int j=0;j<Num_1D;j++)
    alp1d[j] = alp1d[j] - F[j]/psi1d[j];
    
  // Clean up memory
  free(F);
  free(ud);
  free(ld);
  free(d);
}



namespace PNSMapper {
  double get_edens(const double rho, const double temp, const double ye, const double munu) {

  DECLARE_CCTK_PARAMETERS;
  using namespace ZelmaniM1;

  double edens = 0.0;

  if(CCTK_IsThornActive("ZelmaniM1")) {

       const int nvars = ngroups*nspecies*3;
       double opac_tmp[nvars];
       
       const double xrho = log10(INV_RHO_GF*rho);
       const double xtemp = log10(temp);
       
       linterp_many(xrho,xtemp,ye,
		    opac_tmp, alltables, nrho, ntemp, nye, nvars,
		    rho_points, temp_points, ye_points);

       for (int ig=0;ig<ngroups*nspecies;ig++) {
	 double ab = pow(10.0,opac_tmp[ig]);
	 double em = pow(10.0,opac_tmp[ngroups*nspecies + ig]);
	 edens += em / (ab + 1.0e-30);
       }
  } else if(CCTK_IsThornActive("ZelmaniLeak")) {

    const double pnuconst = PRESS_NU_CONSTANT * PRESS_GF;
    const double F3const = 7.0e0*PI4/60.0e0;
    const double eta = munu/temp;
    const double F3 = F3const + 0.5e0*eta*eta*(PI2 + 0.5e0*eta*eta);
    const double pnu = F3*pnuconst*temp*temp*temp*temp;
    edens = 3.0e0 * pnu * exp(- 2.0e12*RHO_GF / rho);
  }

  return edens;
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

  printf("mass: %15.6E\n",m);

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

  printf("psi: in %15.6E out %15.6E\n",psi[0],psiout);

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

}
