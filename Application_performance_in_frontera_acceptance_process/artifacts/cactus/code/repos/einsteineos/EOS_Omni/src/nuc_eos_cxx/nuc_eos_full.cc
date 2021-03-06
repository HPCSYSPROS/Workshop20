#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.hh"
#include "helpers.hh"

namespace nuc_eos {
extern "C"
void CCTK_FNAME(nuc_eos_m_kt1_full)(const int *restrict n_in,
			  const double *restrict rho, 
			  const double *restrict temp,
			  const double *restrict ye,
			  double *restrict eps,
			  double *restrict prs,
			  double *restrict ent,
			  double *restrict cs2,
			  double *restrict dedt,
			  double *restrict dpderho,
			  double *restrict dpdrhoe,
			  double *restrict xa,
			  double *restrict xh,
			  double *restrict xn,
			  double *restrict xp,
			  double *restrict abar,
			  double *restrict zbar,
			  double *restrict mue,
			  double *restrict mun,
			  double *restrict mup,
			  double *restrict muhat,
			  int *restrict keyerr,
			  int *restrict anyerr)
{

  using namespace nuc_eos;

  const int n = *n_in;

  *anyerr = 0;
  for(int i=0;i<n;i++) {
    // check if we are fine
    keyerr[i] = checkbounds(rho[i], temp[i], ye[i]);
    if(CCTK_BUILTIN_EXPECT(keyerr[i] != 0, false)) {
      *anyerr = 1;
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
  
    int idx[8];
    double delx,dely,delz;
    const double xrho = log(rho[i]);
    const double xtemp = log(temp[i]);

    get_interp_spots(xrho,xtemp,ye[i],&delx,&dely,&delz,idx);
    // get prs
    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }
    // get eps
    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }
    // get entropy
    {
      const int iv = 2;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(ent[i]),iv);
    }
    // get cs2
    {
      const int iv = 4;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
    }
    // get dedT
    {
      const int iv = 5;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dedt[i]),iv);
    }
    // get dpdrhoe
    {
      const int iv = 6;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
    }
    // get dpderho
    {
      const int iv = 7;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
    }
    // get muhat
    {
      const int iv = 8;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(muhat[i]),iv);
    }
    // get mue
    {
      const int iv = 9;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mue[i]),iv);
    }
    // get mup
    {
      const int iv = 10;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mup[i]),iv);
    }
    // get mun
    {
      const int iv = 11;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mun[i]),iv);
    }
    // get xa
    {
      const int iv = 12;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xa[i]),iv);
    }
    // get xh
    {
      const int iv = 13;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xh[i]),iv);
    }
    // get xn
    {
      const int iv = 14;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xn[i]),iv);
    }
    // get xp
    {
      const int iv = 15;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xp[i]),iv);
    }
    // get abar
    {
      const int iv = 16;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(abar[i]),iv);
    }
    // get zbar
    {
      const int iv = 17;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(zbar[i]),iv);
    }
  }
  
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
#if HAVEGR
    cs2[i] = rho[i] * cs2[i] / (rho[i] + rho[i] * eps[i] + prs[i]);
#endif
  }

  return;
}

extern "C"
void CCTK_FNAME(nuc_eos_m_kt0_full)(const int *restrict n_in,
			  const double *restrict rho, 
			  double *restrict temp,
			  const double *restrict ye,
			  const double *restrict eps,
			  double *restrict prs,
			  double *restrict ent,
			  double *restrict cs2,
			  double *restrict dedt,
			  double *restrict dpderho,
			  double *restrict dpdrhoe,
			  double *restrict xa,
			  double *restrict xh,
			  double *restrict xn,
			  double *restrict xp,
			  double *restrict abar,
			  double *restrict zbar,
			  double *restrict mue,
			  double *restrict mun,
			  double *restrict mup,
			  double *restrict muhat,
			  const double *restrict prec,
			  int *restrict keyerr,
			  int *restrict anyerr)
{

  using namespace nuc_eos;

  const int n = *n_in;

  *anyerr = 0;

  for(int i=0;i<n;i++) {
    
    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    keyerr[i] = checkbounds_kt0_noTcheck(rho[i], ye[i]);
    if(CCTK_BUILTIN_EXPECT(keyerr[i] != 0, false)) {
      *anyerr = 1;
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
    double ltout;
    const double epstot = eps[i]+energy_shift;
    if(CCTK_BUILTIN_EXPECT(epstot>0.0e0, true)) {
      // this is the standard scenario; eps is larger than zero
      // and we can operate with logarithmic tables
      const double lxeps = log(epstot);
#if DEBUG
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,lr,lt,ye[i],lxeps);
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,
	      exp(lr),exp(lt),ye[i],exp(lxeps));
#endif
      nuc_eos_findtemp(lr,lt,ye[i],lxeps,*prec,
		       (double *restrict)(&ltout),&keyerr[i]);
    } else {
      keyerr[i] = 667;
    } // epstot > 0.0


    if(CCTK_BUILTIN_EXPECT(keyerr[i] != 0, false)) {
      *anyerr=1;
    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      // get prs
      {
	const int iv = 0;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
      // get entropy
      {
	const int iv = 2;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(ent[i]),iv);
      }
      // get cs2
      {
	const int iv = 4;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
      }
      // get dedT
      {
	const int iv = 5;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dedt[i]),iv);
      }
      // get dpdrhoe
      {
	const int iv = 6;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
      }
      // get dpderho
      {
	const int iv = 7;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
      }
      // get muhat
      {
	const int iv = 8;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(muhat[i]),iv);
      }
      // get mue
      {
	const int iv = 9;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mue[i]),iv);
      }
      // get mup
      {
	const int iv = 10;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mup[i]),iv);
      }
      // get mun
      {
	const int iv = 11;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mun[i]),iv);
      }
      // get xa
      {
	const int iv = 12;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xa[i]),iv);
      }
      // get xh
      {
	const int iv = 13;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xh[i]),iv);
      }
      // get xn
      {
	const int iv = 14;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xn[i]),iv);
      }
      // get xp
      {
	const int iv = 15;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(xp[i]),iv);
      }
      // get abar
      {
	const int iv = 16;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(abar[i]),iv);
      }
      // get zbar
      {
	const int iv = 17;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(zbar[i]),iv);
      }
    }
  }

  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
#if HAVEGR
    cs2[i] = rho[i]*cs2[i] / (rho[i] + rho[i]*eps[i] + prs[i]);
#endif
  }

  return;
}



} // namespace nuc_eos
