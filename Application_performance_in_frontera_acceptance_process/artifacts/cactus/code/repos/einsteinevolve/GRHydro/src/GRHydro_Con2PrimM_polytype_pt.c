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

/* Set this to be 1 if you want debug output */
#define DEBUG_CON2PRIMM (0)


/* Adiabatic index used for the state equation */

#define MAX_NEWT_ITER   (30)       /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL        (1.0e-10)  /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL    (1.0e-10)  /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER (2)

#define NEWT_TOL2       (1.0e-15)  /* TOL of new 1D^*_{v^2} gnr2 method */
#define MIN_NEWT_TOL2   (1.0e-10)  /* TOL of new 1D^*_{v^2} gnr2 method */

#define W_TOO_BIG        (1.e20)    /* \gamma^2 (\rho_0 + u + p) is assumed
                                      to always be smaller than this.  This
                                      is used to detect solver failures */

#define FAIL_VAL        (1.e30)    /* Generic value to which we set variables when a problem arises */

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/* Local Globals */
struct LocGlob {
  CCTK_REAL Bsq,QdotBsq,Qtsq,D;
  CCTK_REAL W_for_gnr2, rho_for_gnr2, W_for_gnr2_old, rho_for_gnr2_old;
  CCTK_REAL t4,t7,t24,two_Bsq,t300,t400,s200,s100;
} ;


// Declarations: 
static CCTK_INT oned_newton_raphson( CCTK_REAL x[], CCTK_INT n, CCTK_REAL gammaeos, struct LocGlob *lgp,
                                     void (*funcd) (CCTK_REAL [], CCTK_REAL [], 
                                                    CCTK_REAL [], CCTK_REAL [][1], 
                                                    CCTK_REAL *, CCTK_REAL *, CCTK_INT, CCTK_REAL, 
                                                    struct LocGlob *lgp) );

static  void func_W( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][1], 
                     CCTK_REAL *f, CCTK_REAL *df, CCTK_INT n, CCTK_REAL gammaeos, 
                     struct LocGlob * lgp);

static CCTK_INT gnr2( CCTK_REAL x[], CCTK_INT n, CCTK_REAL gammaeos, struct LocGlob *lgp,
                            void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
                                           CCTK_REAL [][1],CCTK_REAL *,CCTK_REAL *,CCTK_INT, CCTK_REAL,
                                           struct LocGlob *lgp) );

static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
                     CCTK_REAL jac[][1], CCTK_REAL *f, CCTK_REAL *df, CCTK_INT n, CCTK_REAL gammaeos, 
                     struct LocGlob *lgp);

void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_Polytype_pt) (
    CCTK_INT  *handle, CCTK_REAL *gamma_eos, 
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *sc_in, 
    CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
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
    CCTK_REAL *det,
    CCTK_INT  *epsnegative,
    CCTK_REAL *retval);

/**********************************************************************/
/**********************************************************************************

  Con2PrimM_Polytype_pt():
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
   dens       =  as defined in GRHydro and are assumed to be densitized (i.e. with sqrt(\gamma))   
   dens       =  D = \sqrt(\gamma) W \rho    
   sc_in      =  K D, where K is the polytropic constant
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
void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_Polytype_pt)  ( 
    CCTK_INT  *handle, CCTK_REAL *gamma_eos, 
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *sc_in,
    CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
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
  CCTK_REAL dens, gammaeos;
  CCTK_REAL QdotB,utsq,gamma_sq;
  CCTK_REAL rho0,u,p,w,gammasq,gamma,W_last,W,vsq;
  CCTK_REAL g_o_WBsq, QdB_o_W;
  CCTK_REAL rho_g, x_rho[1];
  CCTK_REAL sqrt_detg = *sdet;
  CCTK_REAL inv_sqrt_detg = 1./sqrt_detg; 
  CCTK_REAL t2, rho_gm1;
  CCTK_INT  i_increase,ntries ;

  struct LocGlob lg;


  /* Assume ok initially: */
  *retval = 0.;
  *epsnegative = 0; 

  gammaeos = *gamma_eos;

#if(DEBUG_CON2PRIMM)
  fprintf(stdout," *dens        = %26.16e \n", *dens_in         );
  fprintf(stdout," *sx                 = %26.16e \n", *sx_in           );    
  fprintf(stdout," *sy          = %26.16e \n", *sy_in            );    
  fprintf(stdout," *sz                 = %26.16e \n", *sz_in           );    
  fprintf(stdout," *Sc                 = %26.16e \n", *sc_in           );    
  fprintf(stdout," *Bconsx         = %26.16e \n", *Bconsx_in );    
  fprintf(stdout," *Bconsy         = %26.16e \n", *Bconsy_in );    
  fprintf(stdout," *Bconsz         = %26.16e \n", *Bconsz_in );    
  fprintf(stdout," *rho         = %26.16e \n", *rho           );    
  fprintf(stdout," *velx        = %26.16e \n", *velx           );    
  fprintf(stdout," *vely        = %26.16e \n", *vely          );    
  fprintf(stdout," *velz        = %26.16e \n", *velz          );    
  fprintf(stdout," *epsilon     = %26.16e \n", *epsilon     );
  fprintf(stdout," *pressure    = %26.16e \n", *pressure    );
  fprintf(stdout," *Bx                 = %26.16e \n", *Bx           );    
  fprintf(stdout," *By                 = %26.16e \n", *By           );    
  fprintf(stdout," *Bz                 = %26.16e \n", *Bz           );    
  fprintf(stdout," *bsq                = %26.16e \n", *bsq          );    
  fprintf(stdout," *w_lorentz   = %26.16e \n", *w_lorentz   );
  fprintf(stdout," *gxx         = %26.16e \n", *gxx           );    
  fprintf(stdout," *gxy         = %26.16e \n", *gxy           );    
  fprintf(stdout," *gxz         = %26.16e \n", *gxz           );    
  fprintf(stdout," *gyy         = %26.16e \n", *gyy           );    
  fprintf(stdout," *gyz         = %26.16e \n", *gyz           );    
  fprintf(stdout," *gzz         = %26.16e \n", *gzz           );    
  fprintf(stdout," *uxx         = %26.16e \n", *uxx           );    
  fprintf(stdout," *uxy         = %26.16e \n", *uxy           );    
  fprintf(stdout," *uxz                = %26.16e \n", *uxz          );    
  fprintf(stdout," *uyy         = %26.16e \n", *uyy           );    
  fprintf(stdout," *uyz         = %26.16e \n", *uyz           );    
  fprintf(stdout," *uzz                = %26.16e \n", *uzz          );    
  fprintf(stdout," *sdet                = %26.16e \n", *sdet          );    
  fprintf(stdout," *epsnegative = %10d    \n", *epsnegative );
  fprintf(stdout," *retval      = %26.16e \n", *retval      );
  fflush(stdout);
#endif

  /* First undensitize all conserved variables : */
  sx   = (  *sx_in)   * inv_sqrt_detg;
  sy   = (  *sy_in)   * inv_sqrt_detg;
  sz   = (  *sz_in)   * inv_sqrt_detg;
  dens = (*dens_in)   * inv_sqrt_detg;

  usx  =  (*uxx)*sx + (*uxy)*sy + (*uxz)*sz;
  usy  =  (*uxy)*sx + (*uyy)*sy + (*uyz)*sz;
  usz  =  (*uxz)*sx + (*uyz)*sy + (*uzz)*sz;

  *Bx = (*Bconsx_in) * inv_sqrt_detg; 
  *By = (*Bconsy_in) * inv_sqrt_detg; 
  *Bz = (*Bconsz_in) * inv_sqrt_detg; 

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  
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

  lg.Qtsq = (usx * sx  +  usy * sy  +  usz * sz) ;

  lg.D = dens;

  t2 = lg.D*lg.D;
  lg.t4 = lg.QdotBsq*t2;
  lg.t7 = lg.Bsq*lg.Bsq;
  lg.t24 = 1/t2;
  lg.two_Bsq = lg.Bsq + lg.Bsq;
  lg.t300 = lg.QdotBsq*lg.Bsq*t2;
  lg.t400 = lg.Qtsq*t2;

  lg.s200 = lg.D*gammaeos*(*sc_in);
  CCTK_REAL g_o_gm1 = (gammaeos/(gammaeos-1.));
  lg.s100 = g_o_gm1*(*sc_in);

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
  CCTK_REAL gm1 = gammaeos-1.;
  rho_gm1 = pow(rho0,gm1);
  p = (*sc_in) * rho_gm1 / gamma; 
  u = p / gm1;
  w = rho0 + u + p ;

  W_last = w*gammasq ;


  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*lg.Bsq ) 
            - lg.QdotBsq*(2.*W_last + lg.Bsq) ) <= W_last*W_last*(lg.Qtsq-lg.Bsq*lg.Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  lg.W_for_gnr2 = lg.W_for_gnr2_old = W_last;
  lg.rho_for_gnr2 = lg.rho_for_gnr2_old = rho0;

  // Calculate W: 
  x_1d[0] = W_last;

  *retval = 1.0*oned_newton_raphson( x_1d, 1, gammaeos, &lg, func_W ) ;  

  W = x_1d[0];
        
  /* Problem with solver, so return denoting error before doing anything further */
  if( ((*retval) != 0.) || (W == FAIL_VAL) ) {
    *retval = *retval*100.+1.;
    return;
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      *retval = 3.;
      return;
    }
  }

  rho_g = x_rho[0] = lg.rho_for_gnr2; 

  ntries = 0;
  while (  (*retval = gnr2( x_rho, 1, gammaeos, &lg, func_rho)) &&  ( ntries++ < 10 )  ) { 
    rho_g *= 10.;
    x_rho[0] = rho_g;
  }

  lg.rho_for_gnr2 = x_rho[0];

  if( (*retval != 0) ) {
    *retval = 10;
    return;
  }

  // Calculate v^2 : 
  rho0 = lg.rho_for_gnr2;
  rho_gm1 = pow(rho0,gm1);

  utsq = (lg.D-rho0)*(lg.D+rho0)/(rho0*rho0);

  gamma_sq = 1.+utsq;
  gamma = sqrt(gamma_sq);

  // Calculate v^2:
  if( vsq >= 1. ) {
    *retval = 4.;
    return;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  
  w = W / gamma_sq;
  
  //  printf("doublecheck - S, rho, gamma: %e %e %e\n",*sc_in, rho_gm1,gamma);

  p = (*sc_in) * inv_sqrt_detg * rho_gm1 / gamma;

  u = p / gm1;

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
  *bsq = lg.Bsq * (1.-vsq) + QdB_o_W*QdB_o_W;

  *velx = g_o_WBsq * ( usx + QdB_o_W*(*Bx) ) ;
  *vely = g_o_WBsq * ( usy + QdB_o_W*(*By) ) ;
  *velz = g_o_WBsq * ( usz + QdB_o_W*(*Bz) ) ;


#if(DEBUG_CON2PRIMM)
  fprintf(stdout,"rho          = %26.16e \n",*rho      );
  fprintf(stdout,"epsilon      = %26.16e \n",*epsilon  );
  fprintf(stdout,"pressure     = %26.16e \n",*pressure );
  fprintf(stdout,"w_lorentz    = %26.16e \n",*w_lorentz);
  fprintf(stdout,"bsq          = %26.16e \n",*bsq      );
  fprintf(stdout,"velx         = %26.16e \n",*velx     );
  fprintf(stdout,"vely         = %26.16e \n",*vely     );
  fprintf(stdout,"velz         = %26.16e \n",*velz     );
  fprintf(stdout,"gam          = %26.16e \n",gammaeos       );
  fflush(stdout);
#endif

  /* done! */
  return;

}


#if(DEBUG_CON2PRIMM)
/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon 
       their definitions:
    
*********************************************************************/

static void validate_x(CCTK_REAL x[2], CCTK_REAL x0[2] ) 
{
  
  const CCTK_REAL dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */ 
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];
  
  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}
#endif

/************************************************************

  oned_newton_raphson(): 

    -- performs Newton-Rapshon method on an 2d system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static CCTK_INT oned_newton_raphson( CCTK_REAL x[], CCTK_INT n, CCTK_REAL gammaeos, struct LocGlob *lgp,
                                   void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
                                                  CCTK_REAL [][1], CCTK_REAL *, 
                                                  CCTK_REAL *, CCTK_INT, CCTK_REAL, struct LocGlob *) )
{
  CCTK_REAL f, df, dx[1], x_old[1], resid[1], 
    jac[1][1];
  CCTK_REAL errx, x_orig[1];
  CCTK_INT    n_iter, i_extra, doing_extra;

  CCTK_INT   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  //-fast  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n, gammaeos, lgp);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;

    //-fast    for( id = 0; id < n ; id++) { x_old[id] = x[id] ;  }
    x_old[0] = x[0] ;

    /* don't use line search : */
    //-fast    for( id = 0; id < n ; id++) { x[id] += dx[id]  ;  }
    x[0] += dx[0]  ;

//    //METHOD specific:
//    i_increase = 0;
//    while( (( x[0]*x[0]*x[0] * ( x[0] + 2.*lgp->Bsq ) - 
//              lgp->QdotBsq*(2.*x[0] + lgp->Bsq) ) <= x[0]*x[0]*(lgp->Qtsq-lgp->Bsq*lgp->Bsq))
//           && (i_increase < 10) ) {
//      x[0] -= (1.*i_increase) * dx[0] / 10. ;
//      i_increase++;
//    }

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
#if(DEBUG_CON2PRIMM)
    fprintf(stderr,"\ngnr not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
                f,df,x[0],x_old[0],lgp->W_for_gnr2_old,lgp->W_for_gnr2,lgp->rho_for_gnr2_old,lgp->rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }

  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  } 

  

}


/**********************************************************************/
/*********************************************************************************
   func_W()

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void func_W(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
                   CCTK_REAL jac[][1], CCTK_REAL *f, CCTK_REAL *df, CCTK_INT n, CCTK_REAL gammaeos, 
                   struct LocGlob *lgp)
{
  CCTK_INT retval, ntries;
  CCTK_REAL  W, x_rho[1], rho, rho_g ;
  CCTK_REAL  t15,   t200,t1000, drho_dW, rho_gm1 ;
  
  W  = x[0];
  lgp->W_for_gnr2_old = lgp->W_for_gnr2;
  lgp->W_for_gnr2 = W;

  // get rho from NR:
  rho_g = x_rho[0] = lgp->rho_for_gnr2;
  
  ntries = 0;
  while (  (retval = gnr2( x_rho, 1, gammaeos, lgp, func_rho)) &&  ( ntries++ < 10 )  ) { 
    rho_g *= 10.;
    x_rho[0] = rho_g;
  }

#if(DEBUG_CON2PRIMM)
  if( x_rho[0] <= 0. ) { 
    fprintf(stderr,"gnr2 neg rho = %d ,rho_n,rho,rho_o,W,W_o = %26.20e %26.20e %26.20e %26.20e %26.20e \n", retval, x_rho[0], lgp->rho_for_gnr2, lgp->rho_for_gnr2_old, x[0], lgp->W_for_gnr2_old);
    fflush(stderr);
  }
  
  if( retval ) { 
    fprintf(stderr,"gnr2 retval = %d ,rho_n,rho,rho_o,W,W_o = %26.20e %26.20e %26.20e %26.20e %26.20e \n", retval, x_rho[0], lgp->rho_for_gnr2, lgp->rho_for_gnr2_old, x[0], lgp->W_for_gnr2_old);
    fflush(stderr);
  }
#endif 

  lgp->rho_for_gnr2_old = lgp->rho_for_gnr2; 
  rho = lgp->rho_for_gnr2 = x_rho[0];

  CCTK_REAL gm1 = gammaeos-1.;

  rho_gm1 = pow(rho,gm1);
  drho_dW = -rho*rho/( -rho_gm1*lgp->s200 + W*rho);

  t15 = -(lgp->D-rho)*(lgp->D+rho);  // t6-t2
  t200 = W + lgp->two_Bsq;
  t1000 = rho*drho_dW;
  resid[0] = (lgp->t300+(lgp->t4+lgp->t4+(lgp->t400+t15*(lgp->t7+(t200)*W))*W)*W)*lgp->t24;
  jac[0][0] = 2*(lgp->t4+(lgp->t400+t15*lgp->t7+(3.0*t15*lgp->Bsq+lgp->t7*t1000+(t15+t15+t1000*(t200))*W)*W)*W)*lgp->t24;

  dx[0] = -resid[0]/jac[0][0];

  *df = - resid[0]*resid[0];
  *f = -0.5*(*df);


//  fprintf(stdout,"QdotBsq = %28.18e ; \n",lgp->QdotBsq );
//  fprintf(stdout,"Sc      = %28.18e ; \n",Sc     );
//  fprintf(stdout,"Bsq     = %28.18e ; \n",lgp->Bsq     );
//  fprintf(stdout,"Qtsq    = %28.18e ; \n",lgp->Qtsq    );
//  fprintf(stdout,"Dc      = %28.18e ; \n",lgp->D       );
//  fprintf(stdout,"drhodW  = %28.18e ; \n",drho_dW  );
//  fprintf(stdout,"W       = %28.18e ; \n",W       );
//  fprintf(stdout,"rho     = %28.18e ; \n",rho     );
//  fprintf(stdout,"resid_W = %28.18e ; \n",resid[0] );
//  fprintf(stdout,"jac_W   = %28.18e ; \n",jac[0][0]);
//  fprintf(stdout,"deriv1 %g %g %g %g \n",W,resid[0],jac[0][0],dx[0]);

  return;

}


/***********************************************************/
/********************************************************************** 

  gnr2()

    -- used to calculate rho from W

*****************************************************************/
static CCTK_INT gnr2( CCTK_REAL x[], CCTK_INT n, CCTK_REAL gammaeos, struct LocGlob *lgp,
                            void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
                                           CCTK_REAL [][1],CCTK_REAL *,CCTK_REAL *,CCTK_INT,CCTK_REAL, struct LocGlob *lgp) )
{
  CCTK_REAL f, df, dx[1], x_old[1], resid[1], 
    jac[1][1];
  CCTK_REAL errx, x_orig[1];
  CCTK_INT    n_iter, i_extra, doing_extra;
  

  CCTK_INT   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  //-fast  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n,gammaeos, lgp);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    //-fast  errx = 0.;
    //-fast    for( id = 0; id < n ; id++) { x_old[id] = x[id] ;  }
    x_old[0] = x[0] ;  

    /* Make the newton step: */
    //-fast    for( id = 0; id < n ; id++) { x[id] += dx[id] ;  }
    x[0] += dx[0] ;

    /* Calculate the convergence criterion */
    //-fast    for( id = 0; id < n ; id++) { errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]); }
    //-fast    errx /= 1.*n;
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]); 

    /* Make sure that the new x[] is physical : */
    // METHOD specific:
    x[0] = fabs(x[0]);


    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*   before stopping                                                         */
    if( (fabs(errx) <= NEWT_TOL2) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 0;
    }

    if( doing_extra == 1 ) i_extra++ ;

    // See if we've done the extra iterations, or have done too many iterations:
    if( ((fabs(errx) <= NEWT_TOL2)&&(doing_extra == 0)) || 
        (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }  


  /*  Check for bad untrapped divergences : */
  if( (!finite(f)) || (!finite(df)) || (!finite(x[0]))  ) {
#if(DEBUG_CON2PRIMM)
    fprintf(stderr,"\ngnr2 not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
                f,df,x[0],x_old[0],lgp->W_for_gnr2_old,lgp->W_for_gnr2,lgp->rho_for_gnr2_old,lgp->rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }
  
  // Return in different ways depending on whether a solution was found:
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  } 

  

}

/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho from W via the polytrope:

        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
// for the isentropic version:   eq.  (27)
static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
                     CCTK_REAL jac[][1], CCTK_REAL *f, CCTK_REAL *df, CCTK_INT n, CCTK_REAL gammaeos, 
                     struct LocGlob *lgp)
{

  CCTK_REAL rho, W;
  CCTK_REAL t40,t14;
  
  CCTK_REAL gm1 = gammaeos-1.;

  rho = x[0];
  W = lgp->W_for_gnr2;

  t40 = pow(rho,gm1);

  resid[0] = (rho*W+(-t40*lgp->s100-lgp->D)*lgp->D);
  t14 = t40/rho;  // rho^(g-2)
  jac[0][0] = -t14*lgp->s200 + W;
  //  drho_dW = -rho/jac[0][0];

  dx[0] = -resid[0]/jac[0][0];
  *df = - resid[0]*resid[0];
  *f = -0.5*(*df);

  //  fprintf(stdout,"deriv3 %g %g %g %g %g \n",rho,W,resid[0],jac[0][0],dx[0]);
//  fprintf(stdout,"Dc        := %28.18e ; \n",lgp->D);
//  fprintf(stdout,"Sc        := %28.18e ; \n",Sc);
//  fprintf(stdout,"rho       := %28.18e ; \n",rho);
//  fprintf(stdout,"W         := %28.18e ; \n",W);
//  fprintf(stdout,"resid_rho := %28.18e ; \n",resid[0] );
//  fprintf(stdout,"jac_rho   := %28.18e ; \n",jac[0][0] );

  return;

}


/****************************************************************************** 
             END   
 ******************************************************************************/


#undef DEBUG_CON2PRIMM
