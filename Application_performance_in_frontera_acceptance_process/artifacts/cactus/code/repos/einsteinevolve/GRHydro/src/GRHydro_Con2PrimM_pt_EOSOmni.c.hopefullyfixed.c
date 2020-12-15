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

#define W_TOO_BIG        (1.e20)    /* \gamma^2 (\rho_0 + u + p) is assumed
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

struct eosomnivars {
  CCTK_INT eoshandle,eoskeytemp;
  CCTK_REAL eosprec,eos_y_e[1],eos_temp[1];
  CCTK_INT eoskeyerr[1],eosanyerr[1];
} ;

// Declarations: 
static CCTK_REAL vsq_calc(CCTK_REAL W, struct LocGlob *lgp);

static CCTK_INT twod_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos, struct LocGlob *lgp,
                                     void (*funcd) (CCTK_REAL [], CCTK_REAL [], 
                                                    CCTK_REAL [], CCTK_REAL [][2], 
                                                    CCTK_REAL *, CCTK_REAL *, CCTK_REAL, struct LocGlob *) );

static CCTK_INT threed_newton_raphson_omni( CCTK_REAL x[], struct eosomnivars *eosvars, struct LocGlob *lgp,
                                     void (*funcd) (CCTK_REAL [], CCTK_REAL [], 
                                                    CCTK_REAL [], CCTK_REAL [][3], 
                                                    CCTK_REAL *, CCTK_REAL *, 
                                                    struct eosomnivars *eosvars, struct LocGlob *) );

static  void func_vsq( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][2], 
                       CCTK_REAL *f, CCTK_REAL *df, CCTK_REAL gammaeos, struct LocGlob *lgp);

static  void func_vsq_eosomni( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][3], 
                       CCTK_REAL *f, CCTK_REAL *df, struct eosomnivars *eosvars, struct LocGlob *lgp);

static CCTK_REAL x1_of_x0(CCTK_REAL x0, struct LocGlob *lgp ) ;


// EOS STUFF:
static CCTK_REAL eos_info(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL *dpdw, CCTK_REAL *dpdvsq, CCTK_REAL gammaeos, struct LocGlob *lgp);
static CCTK_REAL eos_info_eosomni(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL eps0, struct LocGlob *lgp);

/* pressure as a function of rho0 and u */
static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u, CCTK_REAL gammaeos)
{
  return((gammaeos - 1.)*u) ;
}

static CCTK_REAL pressure_rho0_eps_eosomni(CCTK_REAL rho0, CCTK_REAL eps, CCTK_REAL* dpdrho, CCTK_REAL* dpdeps, struct eosomnivars *eosvars);

/* Pressure as a function of rho0 and w = rho0 + u + p */
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w,CCTK_REAL gammaeos)
{
  return((gammaeos-1.)*(w - rho0)/gammaeos) ;
}

void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_pt) (
    CCTK_INT  *handle, CCTK_INT *keytemp, CCTK_REAL *prec, 
    CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *tau_in, CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
    CCTK_REAL *y_e_in, CCTK_REAL *temp_in,
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
   det        =  \gamma
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


void CCTK_FCALL CCTK_FNAME(GRHydro_Con2PrimM_pt)  ( 
    CCTK_INT  *handle, CCTK_INT *keytemp, CCTK_REAL *prec,
    CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in, 
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in, 
    CCTK_REAL *tau_in,     
    CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in, 
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
  CCTK_REAL x_3d[3];
  CCTK_REAL sx, sy, sz;
  CCTK_REAL usx, usy, usz;
  CCTK_REAL tau, dens, gammaeos;
  CCTK_REAL QdotB;
  CCTK_REAL rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,vsq;
  CCTK_REAL g_o_WBsq, QdB_o_W;
  CCTK_REAL sqrt_detg = *sdetg;
  CCTK_REAL inv_sqrt_detg = 1./sqrt_detg; 
  CCTK_INT i,j, i_increase ;

  DECLARE_CCTK_PARAMETERS;
  
  struct LocGlob lg; 
  struct eosomnivars eosvars;

  eosvars.eoshandle = *handle;
  //printf("handle = %i\n",*handle);
  eosvars.eoskeytemp = *keytemp;
  eosvars.eosprec = *prec;
  eosvars.eos_y_e[0] = *y_e_in;
  eosvars.eos_temp[0] = *temp_in;
  eosvars.eoskeyerr[0] = 0;
  eosvars.eosanyerr[0] = 0;

  gammaeos = *gamma_eos;

  /* Assume ok initially: */
  *retval = 0.;
  *epsnegative = 0; 

#if(DEBUG_CON2PRIMM)
  fprintf(stdout," *dens        = %26.16e \n", *dens_in         );
  fprintf(stdout," *sx                 = %26.16e \n", *sx_in           );    
  fprintf(stdout," *sy          = %26.16e \n", *sy_in            );    
  fprintf(stdout," *sz                 = %26.16e \n", *sz_in           );    
  fprintf(stdout," *tau         = %26.16e \n", *tau_in           );    
  fprintf(stdout," *Bconsx         = %26.16e \n", *Bconsx_in );    
  fprintf(stdout," *Bconsy         = %26.16e \n", *Bconsy_in );    
  fprintf(stdout," *Bconsz      = %26.16e \n", *Bconsz_in );    
  fprintf(stdout," *rho         = %26.16e \n", *rho           );    
  fprintf(stdout," *velx        = %26.16e \n", *velx           );    
  fprintf(stdout," *vely        = %26.16e \n", *vely          );    
  fprintf(stdout," *velz        = %26.16e \n", *velz          );    
  fprintf(stdout," *epsilon     = %26.16e \n", *epsilon     );
  fprintf(stdout," *temp_in     = %26.16e \n", *temp_in     );
  fprintf(stdout," *y_e_in     = %26.16e \n", *y_e_in     );
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
    fprintf(stdout," *retval      = %26.16e \n", *retval      );
    return;
  }

  gammasq = 1. / (1. - vsq);
  gamma  = sqrt(gammasq);
        
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = lg.D / gamma ;
  u = (*epsilon) * rho0;
  CCTK_REAL uold = u; 

    CCTK_REAL dum1,dum2;

  if (*handle==1 || *handle==2) {     
    p = pressure_rho0_u(rho0,u,gammaeos) ;  // EOS
  } else {
    p = pressure_rho0_eps_eosomni(rho0,*epsilon,&dum1,&dum2,&eosvars) ;  // EOSOMNI
  }
  
  w = rho0 + u + p ;

  W_last = w*gammasq ;

  //fprintf(stdout," p                      = %26.16e \n", p                );

  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*lg.Bsq ) 
            - lg.QdotBsq*(2.*W_last + lg.Bsq) ) <= W_last*W_last*(lg.Qtsq-lg.Bsq*lg.Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  // Calculate W and vsq: 
  x_3d[0] =  fabs( W_last );
  x_3d[1] = x1_of_x0( W_last, &lg ) ;

  //Use 2d NR for polytropes!
  if (*handle==1 || *handle==2) {    
    *retval = 1.0*twod_newton_raphson( x_3d, gammaeos, &lg, func_vsq ) ;  

  } else {
    //USE 3d NR for non-polytropes!
    x_3d[2] = u;
    *retval = 1.0*threed_newton_raphson_omni( x_3d, &eosvars, &lg, func_vsq_eosomni ) ;  
  }
  
  W = x_3d[0];
  vsq = x_3d[1];
        
  /* Problem with solver, so return denoting error before doing anything further */
  if( ((*retval) != 0.) || (W == FAIL_VAL) ) {
    *retval = *retval*100.+1.;
    fprintf(stdout," *retval      = %26.16e \n", *retval      );
    return;
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      *retval = 3.;
      fprintf(stdout," *retval      = %26.16e \n", *retval      );
      return;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    *retval = 4.;
    fprintf(stdout," *retval      = %26.16e \n", *retval      );
    return;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = lg.D * gtmp;

  w = W * (1. - vsq) ;
 
  if (*handle==1 || *handle==2) {     
    p = pressure_rho0_w(rho0,w,gammaeos) ;  // EOS
    u = w - (rho0 + p) ;
  } else {
    u=x_3d[2];
    *epsilon = u/rho0;
    CCTK_REAL dum1,dum2;
    p = pressure_rho0_eps_eosomni(rho0,*epsilon,&dum1,&dum2,&eosvars) ;  // EOSOMNI
   // printf("%g %g %g %g\n",rho0,u,*epsilon,p);
  }

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:

  if((*handle==1 || *handle==2) &&  (rho0 <= 0.) || (u <= 0.)) {
    *epsnegative = 1; 
    fprintf(stdout," *epsnegative = %10d    \n", *epsnegative );
    fprintf(stdout," rho0              = %26.16e \n", rho0                ); 
    fprintf(stdout," u                      = %26.16e \n", u                );
    fprintf(stdout," W                      = %26.16e \n", W                );
    fprintf(stdout," vsq                      = %26.16e \n", vsq                );
    fprintf(stdout," uold                      = %26.16e \n", uold                );
    return;
  } else  if(rho0 <= 0.) { 
//    *epsnegative = 1; 
    fprintf(stdout," *epsnegative = %10d    \n", *epsnegative );
    fprintf(stdout," rho0              = %26.16e \n", rho0                ); 
    fprintf(stdout," u                      = %26.16e \n", u                );
    fprintf(stdout," W                      = %26.16e \n", W                );
    fprintf(stdout," vsq                      = %26.16e \n", vsq                );
    fprintf(stdout," uold                      = %26.16e \n", uold                );
    return;
  }

  *rho = rho0;

  if(*handle==1 || *handle==2) {     
    *epsilon = u / rho0;
  } 

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

    -- performs Newton-Rapshon method on an 2d system for polytropes.

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

  return(0);

}


/************************************************************

  threed_newton_raphson_omni(): 

    -- performs Newton-Rapshon method on an 2d system for polytropes.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static CCTK_INT threed_newton_raphson_omni( CCTK_REAL x[], struct eosomnivars *eosvars, struct LocGlob *lgp, 
                            void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
                                           CCTK_REAL [][3], CCTK_REAL *, 
                                           CCTK_REAL *, struct eosomnivars *, struct LocGlob *) )
{
  CCTK_REAL f, df, dx[3], x_old[3];
  CCTK_REAL resid[3], jac[3][3];
  CCTK_REAL errx, x_orig[3];
  CCTK_INT    n_iter, id, jd, i_extra, doing_extra;
  CCTK_REAL dW,dvsq,du,vsq_old,vsq,W,W_old,u,u_old;
  const CCTK_REAL dv = (1.-1.e-15);

  CCTK_INT   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;
  x_old[2] = x_orig[2] = x[2] ;

  vsq_old = vsq = W = W_old = u = u_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, eosvars, lgp);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    x_old[0] = x[0] ; 
    x_old[1] = x[1] ; 
    x_old[2] = x[2] ; 

    /* Make the newton step: */
    x[0] += dx[0]  ; 
    x[1] += dx[1]  ; 
    x[2] += dx[2]  ; 

    //printf("Updating vars: %g %g %g %g %g %g\n",x[0],dx[0],x[1],dx[1],x[2],dx[2]);

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

    // HOT: if( x[2] < 0. ) {  x[2] = 0.; } 
      
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

  return(0);

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

  
  CCTK_REAL  W, vsq, Wsq, p_tmp, dPdvsq, dPdW;
  CCTK_REAL res0, QB2Winv2,res1,j11,detJinv,   Winv,   QB2Winv3,   j10,   B2plusW,   detJ,   t36,   mj01,   j00;


  W = x[0];
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp = eos_info(W, vsq, &dPdW, &dPdvsq, gammaeos, lgp);

  // These expressions were calculated using Mathematica, but made into efficient 
  // code using Maple.  Since we know the analytic form of the equations, we can 
  // explicitly calculate the Newton-Raphson step: 

  //j11 = dP/dv^2-B^2/2
  j11 = -lgp->half_Bsq+dPdvsq;

  B2plusW = lgp->Bsq+W;

  //mj01 is B2plusW squared = - (partial Eq. 4 / partial v^2)
  mj01 = B2plusW*B2plusW;

  Winv = 1/W;

  QB2Winv2 = lgp->QdotBsq*Winv*Winv;

  //Eq. 4 - Residual 0: Qtsq - v^2(B^2+W)^2 -QdotBsq(B^2+2W)/W^2
  res0 = lgp->Qtsq-vsq*mj01+QB2Winv2*(lgp->Bsq+W+W);

  //Eq. 5 - Residual 1: -Qdotn - B^2/2(1+v^2)+1/2 QdotBsq/W^2 - W+p
  res1 = -lgp->Qdotn-lgp->half_Bsq*(1.0+vsq)+0.5*QB2Winv2-W+p_tmp;

  QB2Winv3 = QB2Winv2*Winv;

  //j10 is -QB2Winv3 - 1 + dp/dW - (partial Eq. 5 / partial W)
  j10 = -1.0+dPdW-QB2Winv3;

  //This is detJ: j10*mj01/B2W +  -2 j11            *   -j00/2
  detJ = B2plusW*(j10*B2plusW+(lgp->Bsq-2.0*dPdvsq)*(QB2Winv2+vsq*W)*Winv);

  detJinv = 1/detJ;

  // - (Jinv00 * res0 + Jinv01 * res 1)/detJ
  dx[0] = -(j11*res0+mj01*res1)*detJinv;

  //j00 is -2v^2(B^2+W)-2QB2 (B2+W)/W^3 - (partial Eq. 4 / partial W)
  j00 = -2*(vsq+QB2Winv3)*B2plusW;

  // (-Jinv10 * res0 -Jinv11 * res1) / DetJ
  dx[1] = (j10*res0-j00*res1)*detJinv;
  //  detJ = B2plusW*detJ_gcf;
  jac[0][0] = j00;
  jac[0][1] = -mj01;
  jac[1][0] = j10;
  jac[1][1] = j11;
  resid[0] = res0;
  resid[1] = res1;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}

/**********************************************************************/
/*********************************************************************************
   func_vsq_eosomni(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq,u here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_vsq_eosomni(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
                             CCTK_REAL jac[][3], CCTK_REAL *f, CCTK_REAL *df, 
                             struct eosomnivars *eosvars, struct LocGlob *lgp)
{

  CCTK_REAL  W, vsq, u, p_tmp,epsilon,Wsq, dPdvsq, dPdW;
  CCTK_REAL res0, LorInv, rho0, Winv, QB2Winv2, drho0_dv, res1, j11, detJ,detJinv;
  CCTK_REAL QB2Winv3, j10,B2plusW,detJ_gcf,t36,B2plusW_sq, j00,j01,j12,j20,j21,j22;
  CCTK_REAL dpress_dv,dpress_du,dpdrho,dpdeps,c00,c01,c02,c10,c11,c12,c20,c21,c22;


  W = x[0];
  vsq = x[1];
  u = x[2];
  
  Wsq = W*W;
  
  LorInv = sqrt(1.0-vsq);
  rho0 = lgp->D * LorInv;
  epsilon = u/rho0;
  p_tmp = pressure_rho0_eps_eosomni(rho0,epsilon,&dpdrho,&dpdeps,eosvars);

  B2plusW = lgp->Bsq+W;
  B2plusW_sq = B2plusW*B2plusW;
  Winv = 1/W;
  QB2Winv2 = lgp->QdotBsq*Winv*Winv;
  QB2Winv3=QB2Winv2*Winv;

  //Eq. 4: Qtsq-v^2(B^2+W)^2-QdotBsq(B^2+2W)/W^2=0 <-No u or p dependence

  resid[0] = lgp->Qtsq - vsq*B2plusW_sq+QB2Winv2*(lgp->Bsq+W+W);

  j00 = -2*(vsq+QB2Winv3)*B2plusW;
  j01 = -1.0*B2plusW_sq;


  //Eq. 5: -Qdotn - B^2(1+vsq)/2 + QdotBsq/2/W^2 - vsq W - rho_0(vsq) - u = 0
  //rho0 = D * sqrt(1-vsq)

  resid[1] = -lgp->Qdotn - lgp->half_Bsq*(1.0+vsq) + 0.5*QB2Winv2 - vsq*W - rho0 - u;

  drho0_dv = -0.5*lgp->D / LorInv;

  j10 = -vsq-QB2Winv3;
  j11 = -lgp->half_Bsq - W - drho0_dv;

  //Eq. 6: u+ p - W(1-vsq) + rho0 = 0  =>  p-W = -W vsq - rho0 - u

  resid[2] = u + p_tmp - W*(1.0-vsq) + rho0;

  //dp/dv = (dp/drho)u *  drho/dv
  //dp/drho_u = dp/drho_eps - eps/rho0 dpdeps 
  dpress_dv = drho0_dv*dpdrho+
    u/2.0/lgp->D*pow(1.0-vsq,-1.5)*dpdeps;

  //dp/du = 1/rho dp deps, since rho0 is function of v only
  dpress_du=dpdeps/rho0;

  jac[0][0] = j00;
  jac[0][1] = j01;
  jac[0][2] = 0.0;
  jac[1][0] = j10;
  jac[1][1] = j11;
  jac[1][2] = -1.0;
  jac[2][0] = vsq-1.0; 
  jac[2][1] = dpress_dv + W + drho0_dv;
  jac[2][2] = dpress_du + 1.0;

  c00 = j11*jac[2][2]+jac[2][1];
  c01 = -1*j01*jac[2][2];
  c02 = j01*jac[1][2];

  c10 = -jac[2][0]-j10*jac[2][2];
  c11 = j00*jac[2][2];
  c12 = j00;

  c20 = j10*jac[2][1]-j11*jac[2][0];
  c21 = j01*jac[2][0]-j00*jac[2][1];
  c22 = j00*j11-j01*j10;

  detJ=j00*c00+j01*c10;
  detJinv = 1/detJ;
  
  dx[0] = -(c00*resid[0]+c01*resid[1]+c02*resid[2])*detJinv;
  dx[1] = -(c10*resid[0]+c11*resid[1]+c12*resid[2])*detJinv;
  dx[2] = -(c20*resid[0]+c21*resid[1]+c22*resid[2])*detJinv;

  *df = -resid[0]*resid[0] - resid[1]*resid[1] - resid[2]*resid[2];

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

static CCTK_REAL pressure_rho0_eps_eosomni(CCTK_REAL rho,CCTK_REAL epsilon, CCTK_REAL* dpdrho, CCTK_REAL* dpdeps, struct eosomnivars *eosvars)
{

  CCTK_REAL rhopt[1],epspt[1],press[1];
  rhopt[0]=rho;
  epspt[0]=epsilon;

  EOS_Omni_press(eosvars->eoshandle,eosvars->eoskeytemp,eosvars->eosprec,1,
                 &rho,&epsilon,eosvars->eos_temp,
                 eosvars->eos_y_e,press,eosvars->eoskeyerr,eosvars->eosanyerr);

  EOS_Omni_DPressByDRho(eosvars->eoshandle,eosvars->eoskeytemp,eosvars->eosprec,1,
                 &rho,&epsilon,eosvars->eos_temp,
                 eosvars->eos_y_e,dpdrho,eosvars->eoskeyerr,eosvars->eosanyerr);

  EOS_Omni_DPressByDEps(eosvars->eoshandle,eosvars->eoskeytemp,eosvars->eosprec,1,
                 &rho,&epsilon,eosvars->eos_temp,
                 eosvars->eos_y_e,dpdeps,eosvars->eoskeyerr,eosvars->eosanyerr);

  return press[0];

}


/****************************************************************************** 
             END   
 ******************************************************************************/


#undef DEBUG_CON2PRIMM
