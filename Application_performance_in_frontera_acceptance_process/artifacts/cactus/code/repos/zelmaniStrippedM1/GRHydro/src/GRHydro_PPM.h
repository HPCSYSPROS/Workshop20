#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if 1
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static inline void steep(double *x, double *dx, double* dmx, const int i) {
  if ( (x[i+1] - x[i]) * (x[i]-x[i-1]) > 0.0 ) {
    dmx[i] = copysign(1.0,dx[i]) * MIN(fabs(dx[i]),
					MIN(2.0*fabs(x[i]-x[i-1]),
					     2.0*fabs(x[i+1]-x[i])));
  } else {
    dmx[i] = 0.0;
  }
}

static inline double approx_at_cell_interface(double* a, const int i) {
  return 7.0/12.0*(a[i]+a[i+1]) - 1.0/12.0*(a[i-1]+a[i+2]);
}

static inline void monotonize(double* restrict xminus,
			      double* restrict x,
			      double* restrict xplus,
			      const int i) {

  if (  !(xplus[i]==x[i] && x[i]==xminus[i]) 
	&& ( (xplus[i]-x[i])*(x[i]-xminus[i]) <= 0.0 ) ) 
    {
      xminus[i] = x[i];
      xplus[i] = x[i];
    }  else if( 6.0 * (xplus[i]-xminus[i]) * 
		(x[i]-0.5*(xplus[i]+xminus[i])) >
		(xplus[i]-xminus[i])*(xplus[i]-xminus[i]) )
    {
      xminus[i] = 3.0*x[i]-2.0*xplus[i]; 
    } else if( 6.0 * (xplus[i]-xminus[i]) * 
	       (x[i]-0.5*(xplus[i]+xminus[i])) <
	       -(xplus[i]-xminus[i])*(xplus[i]-xminus[i]) ) 
    {
      xplus[i] = 3.0*x[i]-2.0*xminus[i]; 
    }
  
  return;
}	     	    



