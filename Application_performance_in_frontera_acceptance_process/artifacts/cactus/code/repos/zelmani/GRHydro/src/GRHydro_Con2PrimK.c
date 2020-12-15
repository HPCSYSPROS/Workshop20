// GRHydro_LastMoLPostStep.c
//
// Compute is this is the last MoL PostStep call. Code taken from Christian
// Reisswig's rejected MoL changes.
//
// Roland Haas
// Sun Jun  3 17:35:53 PDT 2012

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
//#include "GRHydro_Macros.h"
//#include "GRHydro_Macros.h"
#include "nuc_eos_c2p/nuc_eos.h"

#define DEBUG 0
#define OPENCL 0
#define GPU 0

// TO FIX: define max iterations here instead

#include "nuc_eos_c2p/nuc_eos_c.h"
//const double rho_gf   =  1.61930347e-18;
//const double press_gf =  1.80171810e-39;
//const double eps_gf   =  1.11265006e-21;
//const double time_gf  =  2.03001708e05;
//const double mass_gf  =  5.02765209e-34;
//const double length_gf = 6.77140812e-06;
const double rho_gf     =  1.61887093132742e-18;
const double press_gf =  1.80123683248503e-39;
const double eps_gf   =  1.1126500560536e-21;
const double time_gf = 2.03040204956746e05;
const double mass_gf =  5.02916918125126e-34;
const double length_gf = 6.77269222552442e-06;
  

const double min_press = 1.0e-30;
const double clite    = 2.99792458e10;
//const double pi       = 3.14159265358979e0;
const double ggrav    = 6.673e-8;
const double msun     = 1.98899999860571e33;
const double b_gf     = 4.24439887e-20; // CHECK ACCURACY

extern void initialize_host(int);

extern double con2prim_opencl(int, int, double *, int, double *, double *, double *, double *, double,int *,int,int,
		      int *ints,double *doubles,double *fourtables,double *logrho,double *logtemp,double *yes,int *ivs_short);

extern void cpp_createbuffers_(int, int, int, int, int *);
extern void cpp_setkernelargs_(int EOS, double g_det, int CODETEST);
extern void cpp_writeonce_(int EOS, double * EOS_const,double g_con[4][4], double g_cov[4][4],
		    int *ints,double *doubles,double *fourtables,double *logrho,double *logtemp,double *yes,int *ivs_short);

extern double P_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *myvars);

void Conservative2PrimitiveK(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
//  DECLARE_CCTK_FUNCTIONS;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];
  int ax = cctk_ash[0];
  int ay = cctk_ash[1];
  int az = cctk_ash[2];

  // Choose EOS
  // 0: Ideal gas, 1: Polytropic, 2: Hybrid, 3: table
  int o; 
  int EOS = 3; 
  int l,numprims;

  if(EOS==3)
    numprims=12; // add ye and temp
  else
    numprims=8;

  // Need to create struct because function
  // P_prim_ requires a struct be passed in
  int ints[3];
  double doubles[13];

  double rho_cgs = 3.0e14; 
  double epsilon = 0.0;
  double ye = 0.0;
  double temp = 0.0;
  struct nuc_eos_vars myvars_driver;

  // CODETEST = 0: 10% perturbation
  // CODETEST = 1: use "safe guess" values from Cerdá-Durán
  int CODETEST = 0;

  // TO FIX: use a structure to make this clearer?
  double EOS_const[5];
  for(o=0;o<5;o++)
    EOS_const[o] = 0.0;
  EOS_const[2] = 2.5; //gamma2
  EOS_const[3] = 1.5; //gammath
  EOS_const[4] = ((2.0e14)*rho_gf); //rho_nuc
  
//  char * input = argv[1];

  if(EOS==3){
//    nuc_eos_C_ReadTable("/panfs/ds06/sxs/pmoesta/simulations/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5");
//    printf("nrho: %d\n",nrho);
//    printf("ntemp: %d\n",ntemp);
//    printf("nye: %d\n\n",nye);

   ints[0] = nrho;
   ints[1] = ntemp;
   ints[2] = nye;
   doubles[0] = energy_shift;
   doubles[1] = dtemp;
   doubles[2] = dtempi;
   doubles[3] = drho;
   doubles[4] = drhoi;
   doubles[5] = dye;
   doubles[6] = dyei;
   doubles[7] = eos_rhomax;
   doubles[8] = eos_rhomin;
   doubles[9] = eos_tempmin;
   doubles[10] = eos_tempmax;
   doubles[11] = eos_yemin;
   doubles[12] = eos_yemax;
   myvars_driver.nrho = nrho;
   myvars_driver.ntemp = ntemp;
   myvars_driver.nye = nye;
   myvars_driver.fourtables = fourtables;
   myvars_driver.logrho = logrho;
   myvars_driver.logtemp = logtemp;
   myvars_driver.yes = yes;
   myvars_driver.energy_shift = energy_shift;
   myvars_driver.dtemp = dtemp;
   myvars_driver.dtempi = dtempi;
   myvars_driver.drho = drho;
   myvars_driver.drhoi = drhoi;
   myvars_driver.dye = dye;
   myvars_driver.dyei = dyei;
   myvars_driver.eos_rhomax = eos_rhomax;
   myvars_driver.eos_rhomin = eos_rhomin;
   myvars_driver.eos_tempmin = eos_tempmin;
   myvars_driver.eos_tempmax = eos_tempmax;
   myvars_driver.eos_yemin = eos_yemin;
   myvars_driver.eos_yemax = eos_yemax;
   myvars_driver.ivs_short = ivs_short;
  }


#pragma omp parallel 
{
  int i,j,k;
  int index;
  int count;  
  int m;
  double radius; 
  double* original_prim = malloc(sizeof(double) * numprims);
  double g[4][4];
  double glin[4][4];
  double* prim = malloc(sizeof(double) * numprims);
  double* con = malloc(sizeof(double) * numprims);
  double* con2 = malloc(sizeof(double) * numprims);
  if (prim == NULL)
  {
	printf("Could not allocate that much memory\n");
	exit(1);
  }
  // Set up the metric
  double g_det = -1;
  double inv_g_det = -1;
  double sqrt_g_det = -1;
  for(i=0;i<4;i++){
    for(m=0;m<4;m++){
      g[i][m]=0.0;
      glin[i][m] = 0.0;
    }
  }
  g[0][0] = -1.0;
  glin[0][0] = -1.0;
  g[1][1] = 1.0;
  glin[1][1] = 1.0;
  g[2][2] = 1.0;
  glin[2][2] = 1.0;
  g[3][3] = 1.0;
  glin[3][3] = 1.0;

  int n,keytemp,keyerr,anyerr;
  double lrho;
  double leps;
  double ltemp;
  double ly_e;
  double lpress;
  double Bdotv,Bdotv2,b2,bc2; 
  double dBvcx,dBvcy,dBvcz; 
  double ab0;
  double blowx,blowy,blowz; 
 
  #pragma omp for 
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {

          index = CCTK_GFINDEX3D(cctkGH,i,j,k);
          for(l=0;l<(numprims);l++){
            con[l]=0.0;
            con2[l]=0.0;
            prim[l]=0.0;
          }

  g[1][1] = gxx[index];
  g[1][2] = gxy[index];
  g[1][3] = gxz[index];
  g[2][2] = gyy[index];
  g[2][3] = gyz[index];
  g[3][3] = gzz[index];

  g_det = -pow(gxz[index],2)*gyy[index] + 2.0*gxy[index]*gxz[index]*gyz[index] - gxx[index]*pow(gyz[index],2) - pow(gxy[index],2)*gzz[index] + gxx[index]*gyy[index]*gzz[index];
  sqrt_g_det = sqrt(g_det);
  inv_g_det = 1.0/g_det;

//g_det = 1.0;
//sqrt_g_det = 1.0;
//inv_g_det = 1.0;

  
  glin[1][1] = (-pow(gyz[index],2) + gyy[index]*gzz[index])*inv_g_det;
  glin[1][2] = (gxz[index]*gyz[index] - gxy[index]*gzz[index])*inv_g_det;
  glin[2][2] = (-pow(gxz[index],2) + gxx[index]*gzz[index])*inv_g_det;
  glin[1][3] = (-gxz[index]*gyy[index] + gxy[index]*gyz[index])*inv_g_det;
  glin[2][3] = (gxy[index]*gxz[index] - gxx[index]*gyz[index])*inv_g_det;
  glin[3][3] = (-pow(gxy[index],2) + gxx[index]*gyy[index])*inv_g_det;
// Set primitive variables
  //con[0] = con[D];
  con[0] = dens[index]/sqrt_g_det;
  con[1] = scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]/sqrt_g_det;
  con[2] = scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]/sqrt_g_det;
  con[3] = scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]/sqrt_g_det;
  con[4] = tau[index]/sqrt_g_det;
  con[5] = Bcons[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]/sqrt_g_det;
  con[6] = Bcons[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]/sqrt_g_det;
  con[7] = Bcons[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]/sqrt_g_det;
  //  if(EOS==3){
	con[8] = Y_e_con[index];
	con[9] = temperature[index];
        prim[YE]=con[YE]/(con[0]*sqrt_g_det);
        prim[TEMP]=con[TEMP];
//  printf("Inside C2P: sconx = %g\n",  scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]);
//  }

//  printf("Inside C2P: Success \n");
//  printf("Inside C2P: Temp %g\n", temperature[index]);
//  printf("Inside C2P: Ye %g\n", Y_e_con[index]);
//  printf("Inside C2P: con[8] %g\n",con[8]);
//  printf("Inside C2P: con[9] %g\n",con[9]);

//  struct timeval tim;
//  gettimeofday(&tim, NULL);
//  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  //double t1_v2 = clock();
  prim[1] = g[1][1]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + g[1][2]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + g[1][3]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
  prim[2] = g[1][2]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + g[2][2]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + g[2][3]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
  prim[3] = g[1][3]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + g[2][3]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + g[3][3]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];

  prim[EPS] = eps[index];
//  prim[TEMP] = temperature[index];
//  prim[YE] = Y_e[index];
//  prim[1] = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
//  prim[2] = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
//  prim[3] = vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];

  count = con2prim_MHD_(EOS,EOS_const,prim,con,glin,g,g_det,CODETEST,
                             nrho, ntemp, nye, fourtables, logrho,logtemp,yes,
                             energy_shift, dtemp,dtempi,drho,drhoi,dye,dyei,eos_rhomax,eos_rhomin,
                             eos_tempmin,eos_tempmax,eos_yemin,eos_yemax,ivs_short);

//  rho[index] = 0.000001; 

//  printf("Inside C2P2: Y_e,Y_e = %g,%g\n",prim[YE],Y_e[index]);
  // Convert prim to con
//  prim2con_MHD_(EOS, EOS_const,numprims,prim,con2,glin,g,g_det,
//  nrho, ntemp, nye, fourtables, logrho,logtemp,yes,
//  energy_shift, dtemp,dtempi,drho,drhoi,dye,dyei,eos_rhomax,eos_rhomin,eos_tempmin,eos_tempmax,eos_yemin,eos_yemax,ivs_short
//  );
//////
//  int print = 0;
//  for (i=0;i,10;i++) {
//    double diff=fabs(con[4]-con2[4]);
//    if (diff > 1.0e-10)  {
//  //  if (x[index]==52 && y[index] ==28 && z[index] ==60) {
//    printf("Inside C2P: con[0] %g\n", con[0]);
//    printf("Inside C2P: con[1] %g\n", con[1]);
//    printf("Inside C2P: con[2] %g\n", con[2]);
//    printf("Inside C2P: con[3] %g\n", con[3]);
//    printf("Inside C2P: con[4] %g\n", con[4]);
//    printf("Inside C2P: con[5] %g\n", con[5]);
//    printf("Inside C2P: con[6] %g\n", con[6]);
//    printf("Inside C2P: con[7] %g\n", con[7]);
//    printf("Inside C2P: con[8] %g\n", con[8]);
//    printf("Inside C2P: con[9] %g\n", con[9]);
//    printf("Inside C2P: con2[0] %g\n",con2[0]);
//    printf("Inside C2P: con2[1] %g\n",con2[1]);
//    printf("Inside C2P: con2[2] %g\n",con2[2]);
//    printf("Inside C2P: con2[3] %g\n",con2[3]);
//    printf("Inside C2P: con2[4] %g\n",con2[4]);
//    printf("Inside C2P: con2[5] %g\n",con2[5]);
//    printf("Inside C2P: con2[6] %g\n",con2[6]);
//    printf("Inside C2P: con2[7] %g\n",con2[7]);
//    printf("Inside C2P: con2[8] %g\n", con[8]);
//    printf("Inside C2P: con2[9] %g\n", con[9]);
//    }
//  }

  if (count > 126) {
//  if (x[index]==52 && y[index] ==28 && z[index] ==60) {
//  if (print ==1) {
//  if (fabs(con[4]-con2[4]) > 1.0e-10) {
  printf("Inside C2P: con[0] %g\n", con[0]);
  printf("Inside C2P: con[1] %g\n", con[1]);
  printf("Inside C2P: con[2] %g\n", con[2]);
  printf("Inside C2P: con[3] %g\n", con[3]);
  printf("Inside C2P: con[4] %g\n", con[4]);
  printf("Inside C2P: con[5] %g\n", con[5]);
  printf("Inside C2P: con[6] %g\n", con[6]);
  printf("Inside C2P: con[7] %g\n", con[7]);
  printf("Inside C2P: con[8] %g\n", con[8]);
  printf("Inside C2P: con[9] %g\n", con[9]);
  printf("Inside C2P: con2[0] %g\n",con2[0]);
  printf("Inside C2P: con2[1] %g\n",con2[1]);
  printf("Inside C2P: con2[2] %g\n",con2[2]);
  printf("Inside C2P: con2[3] %g\n",con2[3]);
  printf("Inside C2P: con2[4] %g\n",con2[4]);
  printf("Inside C2P: con2[5] %g\n",con2[5]);
  printf("Inside C2P: con2[6] %g\n",con2[6]);
  printf("Inside C2P: con2[7] %g\n",con2[7]);
  printf("Inside C2P: con2[8] %g\n", con[8]);
  printf("Inside C2P: con2[9] %g\n", con[9]);
  }

//  if(GRHydro_eos_hot_eps_fix != 0 && rho[index] > 5.0e-7 && temperature[index] < 1.2) {
//
//              int keytemp = 1;
//              int anyerr = 0;
//              int keyerr = 0;
//
//              double b2 = bcom_sq[index];
//              temperature[index] = 1.5;
//
//              EOS_Omni_press(3,keytemp,1e-10,1,&rho[index],&eps[index],&temperature[index],&Y_e[index],&press[index],&keyerr,&anyerr);
//
//              tau[index]  = sqrt_g_det * (rho[index]*(1.0+eps[index]+b2/2.0) + press[index])* pow(w_lorentz[index],2) - dens[index] - press[index];        
//
//              keytemp = 0;
//  }

//  printf("Before mask\n");
//  printf("Inside C2P: count %d\n",count);
  if (count > 127) {
    GRHydro_C2P_failed[index] = 21.0;
    printf("Setting mask to 21 at %g %g %g\n", x[index],y[index],z[index]);
    printf("Inside C2P: con[1] %g\n", con[0]);
    printf("Inside C2P: con[1] %g\n", con[1]);
  }

  if (prim[0] != prim[0]) {
    printf("Inside C2P NaN: prim[0] %g\n", prim[0]);
  }
  if (prim[3] != prim[3]) {
    printf("Inside C2P NaN: prim[0] %g\n", prim[0]);
  }

  if (count > 100) {
  printf("Inside C2P: Fail \n");
  printf("Inside C2P: count %d\n",count);
  printf("Inside C2P: x: %g\n", x[index]);
  printf("Inside C2P: y: %g\n", y[index]);
  printf("Inside C2P: z: %g\n", z[index]);
  printf("Inside C2P: dens %g\n", dens[index]);
  printf("Inside C2P: sconx %g\n", scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]);
  printf("Inside C2P: scony %g\n", scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]);
  printf("Inside C2P: sconz %g\n", scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]);
  printf("Inside C2P: Bconx %g\n", Bcons[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]);
  printf("Inside C2P: Bcony %g\n", Bcons[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]);
  printf("Inside C2P: Bconz %g\n", Bcons[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]);
  printf("Inside C2P: temp %g\n", temperature[index]);
  printf("Inside C2P: Y_e %g\n", Y_e[index]);
  printf("Inside C2P: con[0] %g\n",con[0]);
  printf("Inside C2P: con[1] %g\n",con[1]);
  printf("Inside C2P: con[2] %g\n",con[2]);
  printf("Inside C2P: con[3] %g\n",con[3]);
  printf("Inside C2P: con[4] %g\n",con[4]);
  printf("Inside C2P: con[5] %g\n",con[5]);
  printf("Inside C2P: con[6] %g\n",con[6]);
  printf("Inside C2P: con[7] %g\n",con[7]);
  printf("Inside C2P: con[8] %g\n",con[8]);
  printf("Inside C2P: con[9] %g\n",con[9]);
  printf("Inside C2P: prim[0] %g\n",prim[0]);
  printf("Inside C2P: prim[1] %g\n",prim[1]);
  printf("Inside C2P: prim[2] %g\n",prim[2]);
  printf("Inside C2P: prim[3] %g\n",prim[3]);
  printf("Inside C2P: vel[1] %g\n", vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]);
  printf("Inside C2P: vel[2] %g\n", vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]);
  printf("Inside C2P: vel[3] %g\n", vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]);
  printf("Inside C2P: prim[4] %g\n",prim[4]);
  printf("Inside C2P: prim[5] %g\n",prim[5]);
  printf("Inside C2P: prim[6] %g\n",prim[6]);
  printf("Inside C2P: prim[7] %g\n",prim[7]);
  printf("Inside C2P: prim[8] %g\n",prim[8]);
  printf("Inside C2P: prim[9] %g\n",prim[9]);
  }

  if (count < 128) {
  
  rho[index] = prim[0];
  vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = glin[1][1]*prim[1] + glin[1][2]*prim[2] + glin[1][3]*prim[3];
  vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] = glin[1][2]*prim[1] + glin[2][2]*prim[2] + glin[2][3]*prim[3];
  vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] = glin[1][3]*prim[1] + glin[2][3]*prim[2] + glin[3][3]*prim[3];
  Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = prim[5];
  Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] = prim[6];
  Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] = prim[7];
  eps[index] = prim[4];
  press[index] = prim[10];
  w_lorentz[index] = prim[11];
  Y_e[index] = prim[8];
  temperature[index] = prim[9];

  }
  
//  n =1;
//  keytemp = 1;
//  keyerr =0;
//  anyerr =0;
//  lrho = rho[index];
//  leps = 0.0;
//  ltemp = temperature[index];
//  ly_e = Y_e[index];
//  lpress = press[index];

  //EOS_Omni_press(4,keytemp,1.0e-10,n,&lrho,&leps,&ltemp,&ly_e,&lpress,&keyerr,&anyerr);

  //double diff=fabs(eps[index]-leps);
  //if (diff > 1.0e-10)  {
  //  printf("Eps mismatch at x,y,z: %g,%g,%g Temp = %g Y_e = %g, EpsGF = %g, lEps = %g\n", x[index],y[index],z[index],temperature[index],Y_e[index],eps[index],leps);
  //}
//  radius = sqrt(x[index]*x[index] + y[index]*y[index]+z[index]*z[index]);

//  if (radius < 350.0 && temperature[index] < 0.4) {
//    temperature[index] = 0.5;
//    keytemp = 1;
//    n =1;
//    keyerr =0;
//    anyerr =0;
//
//    lrho = rho[index];
//    leps = eps[index];
//    ltemp = temperature[index];
//    ly_e = Y_e[index];
//    lpress = press[index];
//
//    GRHydro_C2P_failed[index] = 0.0;
//    printf("Resetting temperature1 on level %d at x,y,z: %g,%g,%g Temp: %g Y_e = %g\n",GRHydro_reflevel,x[index],y[index],z[index],temperature[index],Y_e[index]);
//
//    EOS_Omni_press(4,keytemp,1.0e-10,n,&lrho,&leps,&ltemp,&ly_e,&lpress,&keyerr,&anyerr);
//
//    rho[index] = lrho;
//    eps[index] = leps;
//    temperature[index] = ltemp;
//    Y_e[index] = ly_e;
//    press[index] = lpress;
//
//    printf("Resetting temperature1 at x,y,z: %g,%g,%g Temp: %g Y_e = %g\n",x[index],y[index],z[index],temperature[index],Y_e[index]);
//
//    double b2 = gxx[index]*pow(Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)],2)+gyy[index]*pow(Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)],2)+gzz[index]*pow(Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)],2) + 
//                   2.0*(gxy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]+gxz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] + 
//                   gyz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]);
//
//    double Bdotv = gxx[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + 
//                   gxy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + 
//                   gxz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] + 
//                   gxy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + 
//                   gyy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + 
//                   gyz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] + 
//                   gxz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + 
//                   gyz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + 
//                   gzz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]; 
//
//    double Bdotv2 = Bdotv*Bdotv;
//    double bc2 = b2/(w_lorentz[index]*w_lorentz[index]) + Bdotv2;
//    tau[index]  = sqrt_g_det * ((rho[index]*(1.0+eps[index]) + press[index]) * pow(w_lorentz[index],2) + b2 - (press[index] + bc2/2.0)) - dens[index];       
//    keytemp = 0;
//    double dBvcx = Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)];
//    double dBvcy = Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)];
//    double dBvcz = Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)];
////    Bdotv=DOT(dvelx,dvely,dvelz,dBvcx,dBvcy,dBvcz)
//    double ab0=w_lorentz[index]*Bdotv;
////    b2 = DOT2(dBvcx,dBvcy,dBvcz)/w**2+Bdotv**2
//    double blowx = (gxx[index]*dBvcx + gxy[index]*dBvcy + gxz[index]*dBvcz)/w_lorentz[index] +
//         w_lorentz[index]*Bdotv*prim[1];
//    double blowy = (gxy[index]*dBvcx + gyy[index]*dBvcy + gyz[index]*dBvcz)/w_lorentz[index] +
//         w_lorentz[index]*Bdotv*prim[2];
//    double blowz = (gxz[index]*dBvcx + gyz[index]*dBvcy + gzz[index]*dBvcz)/w_lorentz[index] + 
//         w_lorentz[index]*Bdotv*prim[3];
//  
//    dens[index] = sqrt_g_det * rho[index] * w_lorentz[index]; 
//
//    scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] = sqrt_g_det * ((rho[index]*(1.0+eps[index])+press[index]+bc2)*w_lorentz[index]*w_lorentz[index] * prim[1] - 
//         ab0*blowx);
//    scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] = sqrt_g_det * ((rho[index]*(1.0+eps[index])+press[index]+bc2)*w_lorentz[index]*w_lorentz[index] * prim[2] - 
//         ab0*blowy);
//    scon[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] = sqrt_g_det * ((rho[index]*(1.0+eps[index])+press[index]+bc2)*w_lorentz[index]*w_lorentz[index] * prim[3] - 
//         ab0*blowz);
////    dtau = sqrt_g_det * ((drho*(1+deps)+dpress+b2)*w*w - dpress-b2/2.0-ab0**2) - ddens 
//  
////    dBconsx = sqrt(det)*dBvcx
////    dBconsy = sqrt(det)*dBvcy
////    dBconsz = sqrt(det)*dBvcz
//  }

//  if (temperature[index] < 0.75 && radius > 70.0 && radius < 100.0) {
//    temperature[index] = 0.8;
//    int keytemp = 1;
//    int n =1;
//    int keyerr =0;
//    int anyerr =0;
//
//    double lrho = rho[index];
//    double leps = eps[index];
//    double ltemp = temperature[index];
//    double ly_e = Y_e[index];
//    double lpress = press[index];
//
//    printf("Resetting temperature2 on level %d at x,y,z: %g,%g,%g Temp: %g Y_e = %g\n",GRHydro_reflevel,x[index],y[index],z[index],temperature[index],Y_e[index]);
//
//    EOS_Omni_press(4,keytemp,1.0e-10,n,&lrho,&leps,&ltemp,&ly_e,&lpress,&keyerr,&anyerr);
//
//    rho[index] = lrho;
//    eps[index] = leps;
//    temperature[index] = ltemp;
//    Y_e[index] = ly_e;
//    press[index] = lpress;
//
//    printf("Resetting temperature2 at x,y,z: %g,%g,%g Temp: %g Y_e = %g\n",x[index],y[index],z[index],temperature[index],Y_e[index]);
//
//    double b2 = gxx[index]*pow(Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)],2)+gyy[index]*pow(Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)],2)+gzz[index]*pow(Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)],2) + 
//                   2.0*(gxy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]+gxz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] + 
//                   gyz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]);
//
//    double Bdotv = gxx[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + 
//                   gxy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + 
//                   gxz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] + 
//                   gxy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + 
//                   gyy[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + 
//                   gyz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)] + 
//                   gxz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0)] + 
//                   gyz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1)] + 
//                   gzz[index]*Bvec[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]*vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2)]; 
//
//    double Bdotv2 = Bdotv*Bdotv;
//    double bc2 = b2/(w_lorentz[index]*w_lorentz[index]) + Bdotv2;
//    tau[index]  = sqrt_g_det * ((rho[index]*(1.0+eps[index]) + press[index]) * pow(w_lorentz[index],2) + b2 - (press[index] + bc2/2.0)) - dens[index];       
//    keytemp = 0;
//  }
}
}


//  gettimeofday(&tim, NULL);
//  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  //double t2_v2 = clock();

//  kernelExecTime=(double)(t2-t1);

}


