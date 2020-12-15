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
  CCTK_REAL Bsq,QdotBsq,Qtsq,Qdotn,D,half_Bsq ;
} ;

// Declarations: 
static CCTK_REAL vsq_calc(CCTK_REAL W, struct LocGlob *lgp);

static CCTK_INT twod_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos, struct LocGlob *lgp,
                                     void (*funcd) (CCTK_REAL [], CCTK_REAL [], 
                                                    CCTK_REAL [], CCTK_REAL [][2], 
                                                    CCTK_REAL *, CCTK_REAL *, CCTK_REAL, struct LocGlob *) );

static  void func_vsq( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][2], 
                       CCTK_REAL *f, CCTK_REAL *df, CCTK_REAL gammaeos, struct LocGlob *lgp);

static CCTK_REAL x1_of_x0(CCTK_REAL x0, struct LocGlob *lgp ) ;


// EOS STUFF:
static CCTK_REAL eos_info(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL *dpdw, CCTK_REAL *dpdvsq, CCTK_REAL gammaeos, struct LocGlob *lgp);
/* pressure as a function of rho0 and u */
static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u, CCTK_REAL gammaeos)
{
  return((gammaeos - 1.)*u) ;
}

/* Pressure as a function of rho0 and w = rho0 + u + p */
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w,CCTK_REAL gammaeos)
{
  return((gammaeos-1.)*(w - rho0)/gammaeos) ;
}


void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_ptold) (
    CCTK_INT  *handle, CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *tau_in, CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
    CCTK_REAL *rho, 
    CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz,
    CCTK_REAL *epsilon, CCTK_REAL *pressure,
    CCTK_REAL *Bx, CCTK_REAL *By, CCTK_REAL *Bz, CCTK_REAL *bsq,
    CCTK_REAL *w_lorentz, 
    CCTK_REAL *gxx, CCTK_REAL *gxy, CCTK_REAL *gxz, 
    CCTK_REAL *gyy, CCTK_REAL *gyz, CCTK_REAL *gzz, 
    CCTK_REAL *uxx, CCTK_REAL *uxy, CCTK_REAL *uxz,
    CCTK_REAL *uyy, CCTK_REAL *uyz, CCTK_REAL *uzz,
    CCTK_REAL *sdet,
    CCTK_INT  *epsnegative,
    CCTK_REAL *retval);

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
void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_ptold)  ( 
    CCTK_INT  *handle, CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *tau_in,     
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
  CCTK_REAL x_2d[2];
  CCTK_REAL sx, sy, sz;
  CCTK_REAL usx, usy, usz;
  CCTK_REAL tau, dens, gammaeos;
  CCTK_REAL QdotB;
  CCTK_REAL rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,vsq;
  CCTK_REAL g_o_WBsq, QdB_o_W;
  CCTK_REAL sqrt_detg = *sdet;
  CCTK_REAL inv_sqrt_detg = 1./sqrt_detg; 
  CCTK_INT  i_increase;
 
  DECLARE_CCTK_PARAMETERS;
 
  struct LocGlob lg; 

  gammaeos = *gamma_eos;

  /* Assume ok initially: */
  *retval = 0.;
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
  u = (*epsilon) * rho0;
  p = pressure_rho0_u(rho0,u,gammaeos) ;  // EOS
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
  
  // Calculate W and vsq: 
  x_2d[0] =  fabs( W_last );
  x_2d[1] = x1_of_x0( W_last, &lg ) ;
  *retval = 1.0*twod_newton_raphson( x_2d, gammaeos, &lg, func_vsq ) ;  

  W = x_2d[0];
  vsq = x_2d[1];
        
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

  // Calculate v^2:
  if( vsq >= 1. ) {
    *retval = 4.;
    return;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = lg.D * gtmp;

  w = W * (1. - vsq) ;
  p = pressure_rho0_w(rho0,w,gammaeos) ;  // EOS
  u = w - (rho0 + p) ;

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
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static CCTK_REAL vsq_calc(CCTK_REAL W, struct LocGlob *lgp)
{
        CCTK_REAL Wsq,Xsq,Bsq_W;
        
        Wsq = W*W ;
        Bsq_W = (lgp->Bsq + W);
        Xsq = Bsq_W * Bsq_W;

        return(  ( Wsq * lgp->Qtsq  + lgp->QdotBsq * (Bsq_W + W)) / (Wsq*Xsq) );
}


/********************************************************************

  x1_of_x0(): 
           
    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

static CCTK_REAL x1_of_x0(CCTK_REAL x0, struct LocGlob *lgp ) 
{
  CCTK_REAL vsq;
  CCTK_REAL dv = 1.e-15;

  vsq = fabs(vsq_calc(x0,lgp)) ; // guaranteed to be positive 

  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   ); 

}

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

/************************************************************

  twod_newton_raphson(): 

    -- performs Newton-Rapshon method on an 2d system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static CCTK_INT twod_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos, struct LocGlob *lgp, 
                            void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
                                           CCTK_REAL [][2], CCTK_REAL *, 
                                           CCTK_REAL *, CCTK_REAL, struct LocGlob *) )
{
  CCTK_REAL f, df, dx[2], x_old[2];
  CCTK_REAL resid[2], jac[2][2];
  CCTK_REAL errx, x_orig[2];
  CCTK_INT    n_iter, i_extra, doing_extra;
  CCTK_REAL vsq,W,W_old;
  const CCTK_REAL dv = (1.-1.e-15);

  CCTK_INT   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;

  vsq = W = W_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, gammaeos, lgp);  /* returns with new dx, f, df */
      

    /* Save old values before calculating the new: */
    errx = 0.;
    x_old[0] = x[0] ; 
    x_old[1] = x[1] ; 

    /* Make the newton step: */
    x[0] += dx[0]  ; 
    x[1] += dx[1]  ; 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    if( x[0] < 0. ) {  x[0] = fabs(x[0]);  } 
    else { 
     if(x[0] > W_TOO_BIG)  { x[0] = x_old[0] ; }
    }

    if( x[1] < 0. ) {  x[1] = 0.; } 
    else { 
      if( x[1] > 1. ) { x[1] = dv; }
    }
      
    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


    /*  Check for bad untrapped divergences : */
  if( (!finite(f)) ||  (!finite(df)) ) {
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
   func_vsq(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_vsq(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
                     CCTK_REAL jac[][2], CCTK_REAL *f, CCTK_REAL *df, CCTK_REAL gammaeos, struct LocGlob *lgp)
{

  
  CCTK_REAL  W, vsq, p_tmp, dPdvsq, dPdW;
  CCTK_REAL t11, t16,t18,t2,t21,   t23,   t24,   t25,   t3,   t35,  t4,   t40;


  W = x[0];
  vsq = x[1];
  
  p_tmp = eos_info(W, vsq, &dPdW, &dPdvsq, gammaeos, lgp);

  // These expressions were calculated using Mathematica, but made into efficient 
  // code using Maple.  Since we know the analytic form of the equations, we can 
  // explicitly calculate the Newton-Raphson step: 

  t2 = -lgp->half_Bsq+dPdvsq;
  t3 = lgp->Bsq+W;
  t4 = t3*t3;
  t23 = 1/W;
  t16 = lgp->QdotBsq*t23*t23;
  t11 = lgp->Qtsq-vsq*t4+t16*(lgp->Bsq+W+W);
  t18 = -lgp->Qdotn-lgp->half_Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(lgp->Bsq-2.0*dPdvsq)*(t16+vsq*W)*t23;
  //  t21 = 1/t3;
  //  t36 = 1/t35;
  t21 = 1/(t3*t35);
  dx[0] = -(t2*t11+t4*t18)*t21;
  t40 = -2*(vsq+t24)*t3;
  dx[1] = -(-t25*t11+t40*t18)*t21;
  //  detJ = t3*t35;
  jac[0][0] = t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}


/********************************************************************** 
 ********************************************************************** 
   
 The following routines specify the equation of state.  All routines 
  above here should be indpendent of EOS.  If the user wishes 
  to use another equation of state, the below functions must be replaced 
  by equivalent routines based upon the new EOS. 

 **********************************************************************
**********************************************************************/

/**********************************************************************/
/********************************************************************** 
  eos_info():
 
      -- returns with all the EOS-related values needed;
 **********************************************************************/
static CCTK_REAL eos_info(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL *dpdw, CCTK_REAL *dpdvsq, CCTK_REAL gammaeos, struct LocGlob *lgp)
{
  register CCTK_REAL ftmp,gtmp;

  ftmp = 1. - vsq;
  gtmp = sqrt(ftmp);

  CCTK_REAL gam_m1_o_gam = ((gammaeos-1.)/gammaeos);

  *dpdw =  gam_m1_o_gam * ftmp ;
  *dpdvsq =  gam_m1_o_gam * ( 0.5 * lgp->D/gtmp  -  W ) ;

  return( gam_m1_o_gam * ( W * ftmp  -  lgp->D * gtmp )  );  // p 

}


/****************************************************************************** 
             END   
 ******************************************************************************/


#undef DEBUG_CON2PRIMM
