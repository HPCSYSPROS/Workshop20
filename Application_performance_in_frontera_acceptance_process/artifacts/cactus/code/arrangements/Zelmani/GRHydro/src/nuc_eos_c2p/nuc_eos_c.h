#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define RHO 0
#define v1_cov 1
#define v2_cov 2
#define v3_cov 3
#define EPS 4

#define D 0
#define S1_cov 1
#define S2_cov 2
#define S3_cov 3
#define TAU 4

#define B1_con 5
#define B2_con 6
#define B3_con 7

#define YE 8
#define TEMP 9

#define NSUBTABLES 4

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern const double rho_gf;
extern const double press_gf;
extern const double eps_gf;
extern const double time_gf;
extern const double mass_gf;
extern const double length_gf;
extern const double min_press;
extern const double clite;
extern const double pi;
extern const double ggrav;
extern const double msun;
extern const double b_gf;

struct nuc_eos_vars {
   int nrho;
   int ntemp;
   int nye;
   double *fourtables;
   double *logrho; 
   double *logtemp;
   double *yes;
   double energy_shift;
   double dtemp, dtempi;
   double drho, drhoi;
   double dye, dyei;
   // min and max values
   double eos_rhomax, eos_rhomin;
   double eos_tempmin, eos_tempmax;
   double eos_yemin, eos_yemax;
   int *ivs_short;
};

void nuc_eos_C_linterp_some(double x, double y, double z,
			    double* f, double* ft, 
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt,struct nuc_eos_vars *struct_ptr);

void nuc_eos_C_linterp_some2(double x, double y, double z,
			    double* f, double* ft, 
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt,struct nuc_eos_vars *struct_ptr,double* dlepsdlrho, double* dlepsdlt, double* dlPdlrho, double* dlPdlt);

void nuc_eos_C_linterp_some3(double x, double y, double z,
			    double* f, double* ft, 
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt,struct nuc_eos_vars *struct_ptr,double* dlepsdlrho, double* dlepsdlt, double* dlPdlrho, double* dlPdlt);

void nuc_eos_C_linterp_for_temp(double x, double y, double z,
				double* f, double* ft, 
				int nx, int ny, int nz, 
				double* xt, double*yt, double* zt,
				double* linterp_for_temp);

void nuc_eos_C_findtemp(double lr, double lt0, double ye, 
			double leps, double prec, double *lt,
			int *keyerr,struct nuc_eos_vars *struct_ptr);

double P_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *struct_ptr);

double T_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *struct_ptr);

double P_iter_(int EOS, double * EOS_const, double * prim, double * con, double * x,struct nuc_eos_vars *struct_ptr);

void P_dpdr_dpde(int EOS, double * EOS_const, double * x, double * con, double * prim, double * Pprim, double * dPdrho,double * dPde,struct nuc_eos_vars *struct_ptr);

void P_dPdW_dPdZ(double * Pprim, double * dPdW, double * dPdZ, int EOS, double * EOS_const, double * x, double * con, double * prim, struct nuc_eos_vars *struct_ptr);

void nuc_eos_P_dpdr_dpde(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     double* xdpderho,
		     double *xdpdrhoe, int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr);

void nuc_eos_P_only(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr);

void nuc_eos_E_only(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr);

void nuc_eos_dpdrho_dpdt_dedrho_dedt(double xrho, double *xtemp, double xye,
		     double *xenr, double *xenr2, double* xprs, double* xdedt2,
		     double *dlepsdlrho,
		     double * dlepsdlt, double *dlPdlrho, double *dlPdlrho2, double *dlPdlt, int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr, int stepsize); 

void E_dEdr_dEdt_dPdr_dPdt(int EOS, double * EOS_const, double * x, double * con, double * prim, double * Eprim, double * dEdrho, double * dEdt, double * dPdrho, double * dPdt, struct nuc_eos_vars *struct_ptr, int stepsize);

void E_dEdW_dEdZ(double * Eprim, double * dEdW, double * dEdZ, double * dEdT, double * dPdT, int EOS, double * EOS_const, double * x, double * con, double * prim, struct nuc_eos_vars *struct_ptr, int stepsize);
