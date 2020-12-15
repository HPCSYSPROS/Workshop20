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


struct nuc_eos_vars {
   int nrho;
   int ntemp;
   int nye;
   __global double *fourtables;
   __global double *logrho; 
   __global double *logtemp;
   __global double *yes;
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
			    double* f, __global double* ft, 
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    __global double* xt,__global double*yt, __global double* zt,struct nuc_eos_vars *struct_ptr);

void nuc_eos_C_linterp_one(double x, double y, double z,
			    double* f, __global double* ft, 
			    int* ivs,
			    int nx, int ny, int nz,
			    __global double* xt,__global double*yt, __global double* zt,struct nuc_eos_vars *struct_ptr);

void nuc_eos_C_linterp_for_temp(double x, double y, double z,
				double* f, __global double* ft, 
				int nx, int ny, int nz, 
				__global double* xt, __global double*yt, __global double* zt,
				double* linterp_for_temp);



void nuc_eos_C_findtemp(double lr, double lt0, double ye, 
			double leps, double prec, double *lt,
			int *keyerr,struct nuc_eos_vars *struct_ptr);

double P_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *struct_ptr);

void P_dpdr_dpde(int EOS, double * EOS_const, double * x, double * con,double * Pprim, double * dPdrho,double * dPde,struct nuc_eos_vars *struct_ptr);

void P_dPdW_dPdZ(double * Pprim, double * dPdW, double * dPdZ, int EOS, double * EOS_const, double * x, double * con, struct nuc_eos_vars *struct_ptr);

void nuc_eos_P_dpdr_dpde(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     double* xdpderho,
		     double *xdpdrhoe, int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr);

void nuc_eos_P_only(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr);