#include "nuc_eos.h"
#include "nuc_eos_c.h"

extern double P_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *struct_ptr);

void prim2con_MHD_(int EOS, double * EOS_const, int numprims, double * prim, double * con, double g_con[4][4],
double g_cov[4][4], double g_det,int nrho_in, int ntemp_in, int nye_in,
  double *fourtables_in, double *logrho_in, double *logtemp_in, double *yes_in,double energy_shift_in, double dtemp_in, double dtempi_in,
  double drho_in, double drhoi_in, 
  double dye_in, double dyei_in, double eos_rhomax_in, double eos_rhomin_in, double eos_tempmin_in, double eos_tempmax_in, 
  double eos_yemin_in, double eos_yemax_in, int *ivs_short_in){
  
  struct nuc_eos_vars myvars_prim2con;
  if(EOS==3){
    myvars_prim2con.nrho = nrho_in;
    myvars_prim2con.ntemp = ntemp_in;
    myvars_prim2con.nye = nye_in;
    myvars_prim2con.fourtables = fourtables_in;
    myvars_prim2con.logrho = logrho_in;
    myvars_prim2con.logtemp = logtemp_in;
    myvars_prim2con.yes = yes_in;
    myvars_prim2con.energy_shift = energy_shift_in;
    myvars_prim2con.dtemp = dtemp_in;
    myvars_prim2con.dtempi = dtempi_in;
    myvars_prim2con.drho = drho_in;
    myvars_prim2con.drhoi = drhoi_in;
    myvars_prim2con.dye = dye_in;
    myvars_prim2con.dyei = dyei_in;
    myvars_prim2con.eos_rhomax = eos_rhomax_in;
    myvars_prim2con.eos_rhomin = eos_rhomin_in;
    myvars_prim2con.eos_tempmin = eos_tempmin_in;
    myvars_prim2con.eos_tempmax = eos_tempmax_in;
    myvars_prim2con.eos_yemin = eos_yemin_in;
    myvars_prim2con.eos_yemax = eos_yemax_in;
    myvars_prim2con.ivs_short = ivs_short_in;
  }
  
  
  // FORMAT
  // prim - rho, v_j, eps, B^k
  // con - D, S_j, tau, B^k
  
  int i;
  
  // Raise indices
  double v1_con = g_con[1][1]*prim[v1_cov]+g_con[1][2]*prim[v2_cov]+g_con[1][3]*prim[v3_cov];
  double v2_con = g_con[1][2]*prim[v1_cov]+g_con[2][2]*prim[v2_cov]+g_con[2][3]*prim[v3_cov];
  double v3_con = g_con[1][3]*prim[v1_cov]+g_con[2][3]*prim[v2_cov]+g_con[3][3]*prim[v3_cov];
  
  // v^2 = v^i * v_i
  double v_squared = prim[v1_cov]*v1_con + prim[v2_cov]*v2_con + prim[v3_cov]*v3_con;
  
  // Lorentz factor
  double W = 1/sqrt(1-v_squared);

 // printf("W %g\n",W);
  
  // Anton et al. Equation (44)
  double alpha_b0 = W*(prim[B1_con]*prim[v1_cov]+prim[B2_con]*prim[v2_cov]+prim[B3_con]*prim[v3_cov]);
  
  // Lower indices - covariant
  double B1_cov = g_cov[1][1]*prim[B1_con]+g_cov[1][2]*prim[B2_con]+g_cov[1][3]*prim[B3_con];
  double B2_cov = g_cov[1][2]*prim[B1_con]+g_cov[2][2]*prim[B2_con]+g_cov[2][3]*prim[B3_con];
  double B3_cov = g_cov[1][3]*prim[B1_con]+g_cov[2][3]*prim[B2_con]+g_cov[3][3]*prim[B3_con];
  
  // Anton et al. Equations (44) & (45)
  double b1_cov = B1_cov/W + alpha_b0*prim[v1_cov];
  double b2_cov = B2_cov/W + alpha_b0*prim[v2_cov];
  double b3_cov = B3_cov/W + alpha_b0*prim[v3_cov];
  
  // B^2 = B^i * B_i
  double B_squared = prim[B1_con]*B1_cov + prim[B2_con]*B2_cov + prim[B3_con]*B3_cov;
  
  // Anton et al. Equation (46)
  double b_squared = (B_squared + pow(alpha_b0,2.0))/pow(W,2.0);
 
//  double PRESS = P_iter_(EOS,EOS_const,prim,&myvars_prim2con);
  double PRESS = prim[10];
  double P_star = PRESS+b_squared/2.0;
  double eps_star = prim[EPS] + b_squared/(2.0*prim[RHO]);
  double hstar;
  if(EOS==1) // h = h(rho)
    hstar = 1.0 + b_squared/(2.0*prim[RHO]) + P_star/prim[RHO];
  else // h = h(rho, eps)
    hstar = 1.0 + eps_star + P_star/prim[RHO];
  
  // Anton et al. Equation (37)
  // Need factor sqrt(g_det)?
//  printf("prim[0] %g\n",prim[0]);
  con[D] = prim[RHO]*W; // factor or g_det
  
  // Anton et al. Equation (38)
  // S_j = rho*h*W^2*v_j - a*b^0*b_j
  // rho*h = rho + rho*eps + press
  con[S1_cov] = prim[RHO]*hstar*prim[v1_cov]*pow(W,2.0) - alpha_b0*b1_cov;
  con[S2_cov] = prim[RHO]*hstar*prim[v2_cov]*pow(W,2.0) - alpha_b0*b2_cov;
  con[S3_cov] = prim[RHO]*hstar*prim[v3_cov]*pow(W,2.0) - alpha_b0*b3_cov;

//  double S1_con = g_con[1][1]*con[S1_cov]+g_con[1][2]*con[S2_cov]+g_con[1][3]*con[S3_cov];
//  double S2_con = g_con[1][2]*con[S1_cov]+g_con[2][2]*con[S2_cov]+g_con[2][3]*con[S3_cov];
//  double S3_con = g_con[1][3]*con[S1_cov]+g_con[2][3]*con[S2_cov]+g_con[3][3]*con[S3_cov];
 
//  con[S1_cov] = S1_con;
//  con[S2_cov] = S2_con;
//  con[S3_cov] = S3_con;
 
  // Anton et al. Equation (39)
  // tau = rho*h*W^2 - p - (a*b^0)^2 - D
  con[TAU] = prim[RHO]*hstar*pow(W,2.0) - P_star - pow(alpha_b0,2.0) - con[D]; 
  // + con[D] to recover Noble & Gammie
  
  // Need factor sqrt(g_det)?
  con[B1_con] = prim[B1_con];
  con[B2_con] = prim[B2_con];
  con[B3_con] = prim[B3_con];
  
  if(EOS==3){
    con[YE] = prim[YE]*con[D]*sqrt(g_det);
    con[TEMP] = prim[TEMP];
  }
  
  for(i=0;i<numprims;i++){
//    con[i] = sqrt(-g_det)*con[i];
  }
}
