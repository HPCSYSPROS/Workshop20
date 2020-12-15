#ifndef NUC_EOS_HH
#define NUC_EOS_HH

#include "cctk.h"

// TODO: remove hard coded constants
// TODO: introduce defines for table index of variables
#define HAVEGR 1
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 19
#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38

namespace nuc_eos {

  extern double temp0, temp1;
  extern double energy_shift;

// min and max values

  extern double eos_rhomax, eos_rhomin;
  extern double eos_tempmin, eos_tempmax;
  extern double eos_yemin, eos_yemax;
  
  extern double c2p_tempmin;
  extern double c2p_tempmax;

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
  enum eos_var {i_logpress=0, i_logenergy, i_entropy, i_munu, i_cs2, i_dedt,
                i_dpdrhoe, i_dpderho, i_muhat, i_mu_e, i_mu_p, i_mu_n, i_Xa,
                i_Xh, i_Xn, i_Xp, i_Abar, i_Zbar, i_Gamma};
}

namespace nuc_eos_private {

// table data

  extern int nrho;
  extern int ntemp;
  extern int nye;

  extern double * restrict alltables;
  extern double * restrict epstable;
  extern double * restrict logrho;
  extern double * restrict logtemp;
  extern double dlintemp,dlintempi;
  extern double drholintempi;
  extern double dlintempyei;
  extern double drholintempyei;
  extern double * restrict yes;
  extern double dtemp, dtempi;
  extern double drho, drhoi;
  extern double dye, dyei;
  extern double drhotempi;
  extern double drhoyei;
  extern double dtempyei;
  extern double drhotempyei;

}

#endif // NUC_EOS_HH
