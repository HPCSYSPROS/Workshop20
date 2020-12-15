/*

Code in this file was taken from PyMol v0.90 and user by permissing under
the following license agreement contained in the PyMol distribution.  
Trivial modifications have been made to permit incorporation into VMD.



PyMOL Copyright Notice
======================

The PyMOL source code is copyrighted, but you can freely use and copy
it as long as you don't change or remove any of the copyright notices.

----------------------------------------------------------------------
PyMOL is Copyright 1998-2003 by Warren L. DeLano of 
DeLano Scientific LLC, San Carlos, CA, USA (www.delanoscientific.com).

                        All Rights Reserved

Permission to use, copy, modify, distribute, and distribute modified 
versions of this software and its documentation for any purpose and 
without fee is hereby granted, provided that the above copyright 
notice appear in all copies and that both the copyright notice and 
this permission notice appear in supporting documentation, and that 
the names of Warren L. DeLano or DeLano Scientific LLC not be used in 
advertising or publicity pertaining to distribution of the software 
without specific, written prior permission.

WARREN LYFORD DELANO AND DELANO SCIENTIFIC LLC DISCLAIM ALL WARRANTIES 
WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL WARREN LYFORD DELANO
OR DELANO SCIENTIFIC LLC BE LIABLE FOR ANY SPECIAL, INDIRECT OR 
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE 
USE OR PERFORMANCE OF THIS SOFTWARE.
----------------------------------------------------------------------

Where indicated, portions of the PyMOL system are instead protected
under the copyrights of the respective authors.  However, all code in
the PyMOL system is released as non-restrictive open-source software
under the above license or an equivalent license.  

PyMOL Trademark Notice
======================

PyMOL(TM) is a trademark of DeLano Scientific LLC.  Derivate software
which contains PyMOL source code must be plainly distinguished from
the PyMOL package distributed by DeLano Scientific LLC in all publicity,
advertising, and documentation.

The slogans, "Includes PyMOL(TM).", "Based on PyMOL(TM) technology.",
"Contains PyMOL(TM) source code.", and "Built using PyMOL(TM).", may
be used in advertising, publicity, and documentation of derivate
software provided that the notice, "PyMOL is a trademark of DeLano
Scientific LLC.", is included in a footnote or at the end of the document.

All other endorsements employing the PyMOL trademark require specific,
written prior permission.

--Warren L. DeLano (warren@delanoscientific.com)

*/

#include <stdio.h>
#include <math.h>

#include "common.h"

#ifdef R_SMALL
#undef R_SMALL
#endif
#define R_SMALL 0.000000001

static void normalize3d(BigReal *v) {
  BigReal vlen;
  vlen = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (vlen > R_SMALL) {
    v[0] /= vlen;
    v[1] /= vlen;
    v[2] /= vlen;
  } else {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }
}

/*========================================================================*/
BigReal MatrixFitRMS(int n, BigReal *v1, BigReal *v2, const BigReal *wt, BigReal *ttt)
{
  /*
	Subroutine to do the actual RMS fitting of two sets of vector coordinates
	This routine does not rotate the actual coordinates, but instead returns 
	the RMS fitting value, along with the center-of-mass translation vectors 
	T1 and T2 and the rotation vector M, which rotates the translated 
	coordinates of molecule 2 onto the translated coordinates of molecule 1.
  */

  BigReal *vv1,*vv2;
  BigReal m[3][3],aa[3][3],x[3],xx[3];
  BigReal sumwt, tol, sig, gam;
  BigReal sg, bb, cc, err, etmp, tmp;
  int a, b, c, maxiter, iters, ix, iy, iz;
  BigReal t1[3],t2[3];

  /* Initialize arrays. */

  for(a=0;a<3;a++) {
	for(b=0;b<3;b++) {
	  m[a][b] = 0.0F;
	  aa[a][b] = 0.0F;
	}
	m[a][a] = 1.0F;
	t1[a]=0.0F;
	t2[a]=0.0F;
  }

  sumwt = 0.0F;
#if 0
  tol = SettingGet(cSetting_fit_tolerance);
#else
  tol = 0.00001F;
#endif
#if 0
  maxiter = (int)SettingGet(cSetting_fit_iterations);
#else
  maxiter = 10000;
#endif

  /* Calculate center-of-mass vectors */

  vv1=v1;
  vv2=v2;

  if(wt) {
	for(c=0;c<n;c++)
	  {
		for(a=0;a<3;a++) {
		  t1[a] += wt[c]*vv1[a];
		  t2[a] += wt[c]*vv2[a];
		}
		if (wt[c]!=0.0F) {
		  sumwt = sumwt + wt[c];
		} else {
		  sumwt = sumwt + 1.0F; /* WHAT IS THIS? */
		}
		vv1+=3;
		vv2+=3;
	  }
  } else {
	for(c=0;c<n;c++)
	  {
		for(a=0;a<3;a++) {
		  t1[a] += vv1[a];
		  t2[a] += vv2[a];
		}
		sumwt+=1.0F;
		vv1+=3;
		vv2+=3;
	  }
  }
  if(sumwt==0.0F) sumwt = 1.0F;
  for(a=0;a<3;a++) {
	t1[a] /= sumwt;
	t2[a] /= sumwt;
  }
  /* Calculate correlation matrix */
  vv1=v1;
  vv2=v2;
  for(c=0;c<n;c++)
	{
	  if(wt) {
		for(a=0;a<3;a++) {
		  x[a] = wt[c]*(vv1[a] - t1[a]);
		  xx[a] = wt[c]*(vv2[a] - t2[a]);
		}
	  } else {
		for(a=0;a<3;a++) {
		  x[a] = vv1[a] - t1[a];
		  xx[a] = vv2[a] - t2[a];
		}
	  }
	  for(a=0;a<3;a++)
		for(b=0;b<3;b++)
		  aa[a][b] = aa[a][b] + xx[a]*x[b];
	  vv1+=3;
	  vv2+=3;
	}
  if(n>1) {
    /* Primary iteration scheme to determine rotation matrix for molecule 2 */
    iters = 0;
    while(1) {
      /*	for(a=0;a<3;a++)
         {
         for(b=0;b<3;b++) 
         printf("%8.3f ",m[a][b]);
         printf("\n");
         }
         for(a=0;a<3;a++)
         {
         for(b=0;b<3;b++) 
         printf("%8.3f ",aa[a][b]);
         printf("\n");
         }
         printf("\n");
      */
      
      /* IX, IY, and IZ rotate 1-2-3, 2-3-1, 3-1-2, etc.*/
      iz = (iters+1) % 3;
      iy = (iz+1) % 3;
      ix = (iy+1) % 3;
      sig = aa[iz][iy] - aa[iy][iz];
      gam = aa[iy][iy] + aa[iz][iz];

      if(iters>=maxiter) 
        {
#if 0
          PRINTFB(FB_Matrix,FB_Details)
#else
            fprintf(stderr,
#endif
            " Matrix: Warning: no convergence (%1.8f<%1.8f after %d iterations).\n",(BigReal)tol,(BigReal)gam,iters
#if 0
            ENDFB;
#else
                );
#endif
          break;
        }

      /* Determine size of off-diagonal element.  If off-diagonals exceed the
         diagonal elements * tolerance, perform Jacobi rotation. */
      tmp = sig*sig + gam*gam;
      sg = sqrt(tmp);
      if((sg!=0.0F) &&(fabs(sig)>(tol*fabs(gam)))) {
        sg = 1.0F / sg;
        for(a=0;a<3;a++)
          {
            bb = gam*aa[iy][a] + sig*aa[iz][a];
            cc = gam*aa[iz][a] - sig*aa[iy][a];
            aa[iy][a] = bb*sg;
            aa[iz][a] = cc*sg;
            
            bb = gam*m[iy][a] + sig*m[iz][a];
            cc = gam*m[iz][a] - sig*m[iy][a];
            m[iy][a] = bb*sg;
            m[iz][a] = cc*sg;
          }
      } else {
        break;
      }
      iters++;
    }
  }
  /* At this point, we should have a converged rotation matrix (M).  Calculate
	 the weighted RMS error. */
  err = 0.0F;
  vv1=v1;
  vv2=v2;

  normalize3d(m[0]);
  normalize3d(m[1]);
  normalize3d(m[2]);
  for(c=0;c<n;c++) {
	etmp = 0.0F;
	for(a=0;a<3;a++) {
	  tmp = m[a][0]*(vv2[0]-t2[0])
		+ m[a][1]*(vv2[1]-t2[1])
		+ m[a][2]*(vv2[2]-t2[2]);
	  tmp = (vv1[a]-t1[a])-tmp;
	  etmp += tmp*tmp;
	}
	if(wt)
	  err += wt[c] * etmp;
	else 
	  err += etmp;
	vv1+=3;
	vv2+=3;
  }

  err=err/sumwt;
  err=sqrt(err);

  ttt[0]=(BigReal)m[0][0];
  ttt[1]=(BigReal)m[0][1];
  ttt[2]=(BigReal)m[0][2];
  ttt[3]=(BigReal)-t1[0];
  ttt[4]=(BigReal)m[1][0];
  ttt[5]=(BigReal)m[1][1];
  ttt[6]=(BigReal)m[1][2];
  ttt[7]=(BigReal)-t1[1];
  ttt[8]=(BigReal)m[2][0];
  ttt[9]=(BigReal)m[2][1];
  ttt[10]=(BigReal)m[2][2];
  ttt[11]=(BigReal)-t1[2];
  ttt[12]=(BigReal)t2[0];
  ttt[13]=(BigReal)t2[1];
  ttt[14]=(BigReal)t2[2];
  ttt[15]=1.0F; /* for compatibility with normal 4x4 matrices */

#if 0
  if(fabs(err)<R_SMALL4)
    err=0.0F;
#endif

  return((BigReal)err);
}


