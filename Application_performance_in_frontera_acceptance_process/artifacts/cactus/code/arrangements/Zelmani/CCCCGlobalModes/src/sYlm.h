#ifndef _CCCC_sYlm_
#define _CCCC_sYlm_


#include <math.h>
#include <assert.h>


/*
   Calculates the real and imaginary spin-weighted spherical harmonics 
   with given spin s for l,m in spherical coordinates.
   For s=0, we obtain the standard spherical harmonics Ylm.
*/
extern void sYlm(const int ss, const int ll, const int mm, const double theta, const double phi, double* const res_re, double* const res_im);



#endif
