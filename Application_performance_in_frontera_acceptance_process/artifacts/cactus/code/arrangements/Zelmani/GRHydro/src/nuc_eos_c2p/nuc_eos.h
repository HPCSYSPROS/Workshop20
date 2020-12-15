#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NSUBTABLES 4
#define NTABLES 19
//#define DEBUG 1

int nrho;
int ntemp;
int nye;

double *alltables;
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

// table key
// 0 logpress 
// 1 logenergy
// 2 entropy
// 3 munu
// 4 cs2
// 5 dedt
// 6 dpdrhoe
// 7 dpderho
// 8 muhat
// 9 mu_e
// 10 mu_p
// 11 mu_n
// 12 Xa
// 13 Xh
// 14 Xn
// 15 Xp
// 16 Abar
// 17 Zbar
// 18 Gamma

// some vectors for selecting variables for more
// efficient interpolation
int ivs_short[4];

// frontend function declarations

//void nuc_eos_C_ReadTable(char* nuceos_table_name);
