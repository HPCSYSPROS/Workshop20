/***********************************************************************************
    Copyright 2006 Scott C. Noble, Charles F. Gammie, Jonathan C. McKinney, 
                   and Luca Del Zanna.

                        PVS_GRMHD

    This file was derived from PVS_GRMHD.  The authors of PVS_GRMHD include 
    Scott C. Noble, Charles F. Gammie, Jonathan C. McKinney, and Luca Del Zanna.
    PVS_GRMHD is available under the GPL from:
    http://rainman.astro.uiuc.edu/codelib/  

    You are morally obligated to cite the following  paper in his/her 
    scientific literature that results from use of this file:

    [1] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

    PVS_GRMHD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PVS_GRMHD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PVS_GRMHD; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    If the user has any questions, please direct them to Scott C. Noble at 
    scn@astro.rit.edu  . 

***********************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include "cctk.h"
#include "cctk_Parameters.h"

/* Set this to be 1 if you want debug output */
#define DEBUG_CON2PRIMM (0)


/* Adiabatic index used for the state equation */

#define MAX_NEWT_ITER   (30)       /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL        (1.0e-10)  /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL    (1.0e-10)  /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER (2)

#define NEWT_TOL2       (1.0e-15)  /* TOL of new 1D^*_{v^2} gnr2 method */
#define MIN_NEWT_TOL2   (1.0e-10)  /* TOL of new 1D^*_{v^2} gnr2 method */

#define W_TOO_BIG       (1.e20)    /* \gamma^2 (\rho_0 + u + p) is assumed
                                      to always be smaller than this.  This
                                      is used to detect solver failures */

#define FAIL_VAL        (1.e30)    /* Generic value to which we set variables when a problem arises */

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/* Local Globals */
struct LocGlob {
  CCTK_REAL Bsq, QdotBsq, Qtsq, Qdotn, D, half_Bsq, Sc, g_o_gm1,
           W_for_gnr2, rho_for_gnr2, W_for_gnr2_old, rho_for_gnr2_old, drho_dW ;
} ;

// Declarations: 


void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_ptee)  ( 
    CCTK_INT  *handle, CCTK_INT *keytemp, CCTK_REAL *prec,
    CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *tau_in,     
    CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
    CCTK_REAL *entropycons_in,
    CCTK_REAL *y_e_in, CCTK_REAL* temp_in,
    CCTK_REAL *rho, 
    CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz,
    CCTK_REAL *epsilon, CCTK_REAL *pressure,
    CCTK_REAL *Bx, CCTK_REAL *By, CCTK_REAL *Bz, 
    CCTK_REAL *bsq,
    CCTK_REAL *w_lorentz, 
    CCTK_REAL *gxx, CCTK_REAL *gxy, CCTK_REAL *gxz, 
    CCTK_REAL *gyy, CCTK_REAL *gyz, CCTK_REAL *gzz, 
    CCTK_REAL *uxx, CCTK_REAL *uxy, CCTK_REAL *uxz,
    CCTK_REAL *uyy, CCTK_REAL *uyz, CCTK_REAL *uzz,
    CCTK_REAL *sdet,
    CCTK_INT  *epsnegative,
    CCTK_REAL *retval) ;

static CCTK_INT general_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos,
               struct LocGlob *lgp,
               void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                    CCTK_REAL [][1], CCTK_REAL *,
                    CCTK_REAL *, CCTK_REAL, struct LocGlob *) );

static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
          CCTK_REAL jac[][1], CCTK_REAL *f, CCTK_REAL *df, 
          CCTK_REAL gammaeos, struct LocGlob *lgp);

/**********************************************************************/
/**********************************************************************************

  Con2PrimM_pt():
  -----------------------------

   -- Attempts an inversion from GRMHD conserved variables to primitive variables assuming a guess.

   -- Uses the 2D method of Noble et al. (2006): 
       -- Solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method 
       -- Can be used (in principle) with a general equation of state. 

   -- Minimizes two residual functions using a homemade Newton-Raphson routine. 
       -- It is homemade so that it can catch exceptions and handle them correctly, plus it is 
          optimized for this particular  problem. 

  -- Note that the notation used herein is that of Noble et al. (2006) except for the argument 
     list. 


INPUT:  (using GRHydro variable defintions)

   s[x,y,z]   =  scons[0,1,2]  = \alpha \sqrt(\gamma) T^0_i 
   dens, tau  =  as defined in GRHydro and are assumed to be densitized (i.e. with sqrt(\gamma))   
   dens       =  D = \sqrt(\gamma) W \rho    
   tau        =  \alpha^2 \sqrt(\gamma) T^{00} - D 
   g[x,y,z][x,y,x] = spatial metric corresponding to \gamma 
   u[x,y,z][x,y,z] = inverse of the spatial metric, g[x,y,z][x,y,x]
   det        =  sqrt(\gamma)
   B[x,y,z]   =  Bvec[0,1,2] 
   bsq        = b^\mu b_\mu  

   epsnegative = (integer)
               = 0  if rho and epsilon are positive
              != 0  otherwise 


  --  (users should set B[x,y,z] = 0  for hydrodynamic runs) 


OUTPUT:  (using GRHydro variable defintions)
   rho, eps   =  as defined in GRHydro, primitive variables
   vel[x,y,z] =  as defined in GRHydro, primitive variables


RETURN VALUE: of retval = (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
             
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: unphysical vsq = v^2  value at initial guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             ( used to be  5 -> failure: rho,uu <= 0   but now sets epsnegative to non-zero )

**********************************************************************************/

void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_ptee)  ( 
    CCTK_INT  *handle, CCTK_INT *keytemp, CCTK_REAL *prec,
    CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *tau_in,     
    CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
    CCTK_REAL *entropycons_in,
    CCTK_REAL *y_e_in, CCTK_REAL* temp_in,
    CCTK_REAL *rho, 
    CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz,
    CCTK_REAL *epsilon, CCTK_REAL *pressure,
    CCTK_REAL *Bx, CCTK_REAL *By, CCTK_REAL *Bz, 
    CCTK_REAL *bsq,
    CCTK_REAL *w_lorentz, 
    CCTK_REAL *gxx, CCTK_REAL *gxy, CCTK_REAL *gxz, 
    CCTK_REAL *gyy, CCTK_REAL *gyz, CCTK_REAL *gzz, 
    CCTK_REAL *uxx, CCTK_REAL *uxy, CCTK_REAL *uxz,
    CCTK_REAL *uyy, CCTK_REAL *uyz, CCTK_REAL *uzz,
    CCTK_REAL *sdet,
    CCTK_INT  *epsnegative,
    CCTK_REAL *retval)

{
  CCTK_REAL x_1d[1];
  CCTK_REAL sx, sy, sz;
  CCTK_REAL usx, usy, usz;
  CCTK_REAL tau, dens, gammaeos;
  CCTK_REAL QdotB;
  CCTK_REAL rho0,u,p,w,gammasq,gamma,W,vsq;
  CCTK_REAL g_o_WBsq, QdB_o_W;
  CCTK_REAL sqrt_detg = *sdet;
  CCTK_REAL inv_sqrt_detg = 1./sqrt_detg; 
  CCTK_REAL rho_gm1;
  CCTK_REAL utsq;
 
  DECLARE_CCTK_PARAMETERS;
 
  struct LocGlob lg; 

  gammaeos = *gamma_eos;

  /* Assume ok initially: */
  // BCM: But let the driver function take care of its initialization
  //*retval = 0.; 
  *epsnegative = 0; 

#if(DEBUG_CON2PRIMM)
  fprintf(stdout," *dens        = %26.16e \n", *dens_in         );
  fprintf(stdout," *sx          = %26.16e \n", *sx_in     );    
  fprintf(stdout," *sy          = %26.16e \n", *sy_in     );    
  fprintf(stdout," *sz          = %26.16e \n", *sz_in     );    
  fprintf(stdout," *tau         = %26.16e \n", *tau_in    );    
  fprintf(stdout," *Bconsx      = %26.16e \n", *Bconsx_in );    
  fprintf(stdout," *Bconsy      = %26.16e \n", *Bconsy_in );    
  fprintf(stdout," *Bconsz      = %26.16e \n", *Bconsz_in );    
  fprintf(stdout," *entropycons      = %26.16e \n", *entropycons_in );    
  fprintf(stdout," *rho         = %26.16e \n", *rho       );    
  fprintf(stdout," *velx        = %26.16e \n", *velx      );    
  fprintf(stdout," *vely        = %26.16e \n", *vely      );    
  fprintf(stdout," *velz        = %26.16e \n", *velz      );    
  fprintf(stdout," *epsilon     = %26.16e \n", *epsilon     );
  fprintf(stdout," *pressure    = %26.16e \n", *pressure    );
  fprintf(stdout," *Bx          = %26.16e \n", *Bx        );    
  fprintf(stdout," *By          = %26.16e \n", *By        );    
  fprintf(stdout," *Bz          = %26.16e \n", *Bz        );    
  fprintf(stdout," *bsq         = %26.16e \n", *bsq       );    
  fprintf(stdout," *w_lorentz   = %26.16e \n", *w_lorentz   );
  fprintf(stdout," *gxx         = %26.16e \n", *gxx       );    
  fprintf(stdout," *gxy         = %26.16e \n", *gxy       );    
  fprintf(stdout," *gxz         = %26.16e \n", *gxz       );    
  fprintf(stdout," *gyy         = %26.16e \n", *gyy       );    
  fprintf(stdout," *gyz         = %26.16e \n", *gyz       );    
  fprintf(stdout," *gzz         = %26.16e \n", *gzz       );    
  fprintf(stdout," *uxx         = %26.16e \n", *uxx       );    
  fprintf(stdout," *uxy         = %26.16e \n", *uxy       );    
  fprintf(stdout," *uxz         = %26.16e \n", *uxz       );    
  fprintf(stdout," *uyy         = %26.16e \n", *uyy       );    
  fprintf(stdout," *uyz         = %26.16e \n", *uyz       );    
  fprintf(stdout," *uzz         = %26.16e \n", *uzz       );    
  fprintf(stdout," *sdet        = %26.16e \n", *sdet       );    
  fprintf(stdout," *epsnegative = %10d    \n", *epsnegative );
  fprintf(stdout," *retval      = %26.16e \n", *retval      );
  fflush(stdout);
#endif

  /* First undensitize all conserved variables : */
  sx   = (  *sx_in)   * inv_sqrt_detg;
  sy   = (  *sy_in)   * inv_sqrt_detg;
  sz   = (  *sz_in)   * inv_sqrt_detg;
  tau  = ( *tau_in)   * inv_sqrt_detg;
  dens = (*dens_in)   * inv_sqrt_detg;

  usx  =  (*uxx)*sx + (*uxy)*sy + (*uxz)*sz;
  usy  =  (*uxy)*sx + (*uyy)*sy + (*uyz)*sz;
  usz  =  (*uxz)*sx + (*uyz)*sy + (*uzz)*sz;

  *Bx = (*Bconsx_in) * inv_sqrt_detg;
  *By = (*Bconsy_in) * inv_sqrt_detg;
  *Bz = (*Bconsz_in) * inv_sqrt_detg;

  lg.Sc = (*entropycons_in) * inv_sqrt_detg;
  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:

  lg.g_o_gm1 = gammaeos/(gammaeos-1.0);
  
  lg.Bsq = 
    (*gxx) * (*Bx) * (*Bx) + 
    (*gyy) * (*By) * (*By) + 
    (*gzz) * (*Bz) * (*Bz) + 
    2*( 
       (*gxy) * (*Bx) * (*By) +
       (*gxz) * (*Bx) * (*Bz) +
       (*gyz) * (*By) * (*Bz) );
    
  QdotB = (sx * (*Bx) + sy * (*By) + sz * (*Bz)) ;
  lg.QdotBsq = QdotB*QdotB ;

  lg.Qdotn = -(tau + dens) ;

  lg.Qtsq = (usx * sx  +  usy * sy  +  usz * sz) ;

  lg.D = dens;

  lg.half_Bsq = 0.5*lg.Bsq;

  /* calculate W from last timestep and use for guess */
  vsq = 
    (*gxx) * (*velx) * (*velx) + 
    (*gyy) * (*vely) * (*vely) + 
    (*gzz) * (*velz) * (*velz) + 
    2*( 
       (*gxy) * (*velx) * (*vely) +
       (*gxz) * (*velx) * (*velz) +
       (*gyz) * (*vely) * (*velz) );

  if( (vsq < 0.) && (fabs(vsq) < 1.0e-13) ) { 
    vsq = fabs(vsq);
  }
  if(vsq < 0. || vsq > 1. ) {
    *retval = 2.;
    return;
  }

  gammasq = 1. / (1. - vsq);
  gamma  = sqrt(gammasq);
        
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = lg.D / gamma ;
  rho_gm1 = pow(rho0,(gammaeos-1.));
  p = lg.Sc * rho_gm1 / gamma;
  u = p / (gammaeos-1.);
  w = rho0 + u + p ;

//  W_last = w*gammasq ;

  // Calculate W and vsq: 
  x_1d[0] = rho0;
//  *retval = 1.0*twod_newton_raphson( x_2d, gammaeos, &lg, func_vsq ) ;  
  *retval = general_newton_raphson( x_1d, gammaeos, &lg, func_rho ) ;  
  rho0 = x_1d[0];
        
  /* Problem with solver, so return denoting error before doing anything further */
  if( ((*retval) != 0.) || (rho0 == FAIL_VAL) ) {
    *retval = *retval*100.+1.;
    return;
  }
  else{
    if( rho0 > W_TOO_BIG) {
      *retval = 3.;
      return;
    }
  }

  // Calculate v^2:
  utsq = (lg.D-rho0)*(lg.D+rho0)/(rho0*rho0);
  gammasq = 1.+utsq;
  gamma = sqrt(gammasq);

  if( utsq < 0. ) {
    *retval = 4.;
    return;
  }


  // Recover the primitive variables from the scalars and conserved variables:
  rho0 = lg.D / gamma;
  rho_gm1 = pow(rho0,(gammaeos-1.));
  p = lg.Sc * rho_gm1 / gamma;
  u = p / (gammaeos-1.);
  w = rho0 + u + p ;
  W = w * gammasq;

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( (rho0 <= 0.) || (u <= 0.) ) { 
    *epsnegative = 1; 
    return;
  }

  *rho = rho0;
  *epsilon = u / rho0;
  *w_lorentz = gamma; 
  *pressure = p ; 

  g_o_WBsq = 1./(W+lg.Bsq);
  QdB_o_W  = QdotB / W; 
  *bsq = lg.Bsq / gammasq + QdB_o_W*QdB_o_W;

  *velx = g_o_WBsq * ( usx + QdB_o_W*(*Bx) ) ;
  *vely = g_o_WBsq * ( usy + QdB_o_W*(*By) ) ;
  *velz = g_o_WBsq * ( usz + QdB_o_W*(*Bz) ) ;

  if (*rho <= rho_abs_min*(1.0+GRHydro_atmo_tolerance) ) {
    *rho = rho_abs_min;
    *velx = 0.0;
    *vely = 0.0;
    *velz = 0.0;
    *w_lorentz = 1.0;
  }


#if(DEBUG_CON2PRIMM)
  fprintf(stdout,"rho          = %26.16e \n",*rho      );
  fprintf(stdout,"epsilon      = %26.16e \n",*epsilon  );
  fprintf(stdout,"pressure     = %26.16e \n",*pressure );
  fprintf(stdout,"w_lorentz    = %26.16e \n",*w_lorentz);
  fprintf(stdout,"bsq          = %26.16e \n",*bsq      );
  fprintf(stdout,"velx         = %26.16e \n",*velx     );
  fprintf(stdout,"vely         = %26.16e \n",*vely     );
  fprintf(stdout,"velz         = %26.16e \n",*velz     );
  fprintf(stdout,"gammaeos     = %26.16e \n",gammaeos  );
  fflush(stdout);
#endif

  /* done! */
  return;

}

/**********************************************************************/
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static CCTK_INT general_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos,
               struct LocGlob *lgp,
               void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                    CCTK_REAL [][1], CCTK_REAL *,
                    CCTK_REAL *, CCTK_REAL, struct LocGlob *) )
{
  CCTK_REAL f, df, dx[1], x_old[1], resid[1],
    jac[1][1];
  CCTK_REAL errx, x_orig[1];
  CCTK_INT  n_iter, i_extra, doing_extra;
  CCTK_REAL W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ;
  df = f = 1.;
  i_extra = doing_extra = 0;
  x_old[0] = x_orig[0] = x[0] ;

  W = W_old = 0.;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) {

    (*funcd) (x, dx, resid, jac, &f, &df, gammaeos, lgp);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    x_old[0] = x[0] ;
    /* don't use line search : */
    x[0] += dx[0]  ;

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = fabs(x[0]);


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/

    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) ||
   (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (!finite(f)) || (!finite(df)) || (!finite(x[0]))  ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}

/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho Qtsq equation with 
            the definition of W 
        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho
              substituted in. 

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
          CCTK_REAL jac[][1], CCTK_REAL *f, CCTK_REAL *df, 
          CCTK_REAL gammaeos, struct LocGlob *lgp)
{

  CCTK_REAL t1,t2;
  CCTK_REAL t12;
  CCTK_REAL t17;
  CCTK_REAL t3, t100,rhosq,t200,rho,W,dWdrho,dvsqdrho,vsq;

  rho = x[0];
  rhosq = rho*rho;
  t200 = 1./(lgp->D*lgp->D);
  t100 = lgp->g_o_gm1*lgp->Sc*pow(rho,(gammaeos-1.));
  W = lgp->D*( lgp->D + t100 ) / rho ;
  dWdrho = lgp->D * ( -lgp->D + t100*(gammaeos-2.) ) / rhosq;
  t1 = W*W;
  t2 = lgp->Bsq+W;
  //    t3 = pow(Bsq+W,2.0);
  t3 = t2*t2;
  vsq = (lgp->D-rho)*(lgp->D+rho)*t200;
  dvsqdrho = -2*rho*t200;
  resid[0] = t1*(lgp->Qtsq-vsq*t3)+lgp->QdotBsq*(t2+W);
  t12 = lgp->Bsq*lgp->Bsq;
  t17 = dWdrho*vsq;
  jac[0][0] = 2*lgp->QdotBsq*dWdrho
    +((lgp->Qtsq-vsq*t12)*2*dWdrho+(-6*t17*lgp->Bsq-dvsqdrho*t12
                +(-2*dvsqdrho*lgp->Bsq-4*t17-dvsqdrho*W)*W)*W)*W;

  dx[0] = -resid[0]/jac[0][0];
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

  return;

}


/****************************************************************************** 
             END   
 ******************************************************************************/


#undef DEBUG_CON2PRIMM
