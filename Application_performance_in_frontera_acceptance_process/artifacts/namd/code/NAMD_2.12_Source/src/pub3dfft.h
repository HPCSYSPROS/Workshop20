/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PUB3DFFT_H
#define PUB3DFFT_H

typedef struct { double r, i; } doublecomplex;

/* ntable should be 4*max(n1,n2,n3) +15 */
/* size of table should be 3*ntable doubles */
/* size of work should be 2*max(n1,n2,n3) doubles */

int pubz3di(int *n1, int *n2, int *n3, double *table, int *ntable);

int pubz3d(int *isign, int *n1, int *n2,
   int *n3, doublecomplex *w, int *ld1, int *ld2, double
   *table, int *ntable, doublecomplex *work);

/* for real to complex n1 and ld1 must be even */

int pubd3di(int n1, int n2, int n3, double *table, int ntable);

int pubdz3d(int isign, int n1, int n2,
   int n3, double *w, int ld1, int ld2, double
   *table, int ntable, double *work);

int pubzd3d(int isign, int n1, int n2,
   int n3, double *w, int ld1, int ld2, double
   *table, int ntable, double *work);

#endif

