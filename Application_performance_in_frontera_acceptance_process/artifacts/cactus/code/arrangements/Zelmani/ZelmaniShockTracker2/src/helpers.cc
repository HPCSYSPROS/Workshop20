#include <stdlib.h>
#include <stdio.h>
#include <cmath>


double alpha_grid_zst(double rmin, double rmax, int nzones,
		  double drmin, double prec) {

  double alpha,alpha_p;
  double f,fp,dfda;
  double rad,rad_p,dr;
  int it;

  rad = rmin;
  alpha = 1.0e0;
  dr = drmin;
  for(int i=0;i<nzones;i++) {
    rad = rad + dr;
    dr = alpha * dr;
  }
  rad_p = rad;
  alpha_p = alpha;

  rad = rmin;
  alpha = 1.01e0;
  dr = drmin;
  for(int i=0;i<nzones;i++) {
    rad = rad + dr;
    dr = alpha * dr;
  }

  it = 0;
  f = rmax - rad;
  while( (fabs(f/rmax) > prec) && (it < 100) ) {
    dfda =  ( (rmax - rad) - (rmax - rad_p) ) / (alpha - alpha_p);
    rad_p = rad;
    alpha_p = alpha;
    alpha = alpha - f / dfda;
    rad = rmin;
    dr = drmin;
    for(int i=0;i<nzones;i++) {
      rad = rad + dr;
      dr = alpha * dr;
    }
    f = rmax - rad;
    //    printf("%d %15.6E %15.6E %15.6E %15.6E\n",it,alpha,f,fabs(f/rmax),prec);
    it++;
  }
  if (it >= 200) {
    alpha = -1.0e0;
  }

  return alpha;
}
