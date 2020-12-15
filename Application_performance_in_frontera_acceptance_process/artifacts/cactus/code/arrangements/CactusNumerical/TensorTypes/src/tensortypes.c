/* $Header$ */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "tensortypes_imp.h"



/* Calculate the number of elements in an array */
#define ARRSIZE(T) (sizeof(T) / sizeof(*(T)))



/******************************************************************************/



/* A scalar */
static int const scalar_vars[] = {
  0
};
static int const scalar_parity[] = {
  +1
};
static int const scalar_comps[] = {
  0
};
struct tensor const TT_scalar = {
  "scalar", 3, 0,
  ARRSIZE(scalar_vars), scalar_vars, scalar_parity,
  ARRSIZE(scalar_comps), scalar_comps
};



/* A vector */
static int const vector_vars[] = {
  0,1,2
};
static int const vector_parity[] = {
  +1,+1,+1
};
static int const vector_comps[] = {
  0,1,2
};
struct tensor const TT_vector = {
  "vector", 3, 1,
  ARRSIZE(vector_vars), vector_vars, vector_parity,
  ARRSIZE(vector_comps), vector_comps
};



/* A second rank tensor without symmetries */
static int const tensor_vars[] = {
  0,1,2,   3,4,5,   6,7,8
};
static int const tensor_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int const tensor_comps[] = {
  0,1,2,   3,4,5,   6,7,8
};
struct tensor const TT_tensor = {
  "tensor", 3, 2,
  ARRSIZE(tensor_vars), tensor_vars, tensor_parity,
  ARRSIZE(tensor_comps), tensor_comps
};



/* A symmetric second rank tensor */
static int const symmtensor_vars[] = {
  0,1,2,   1,3,4,   2,4,5
};
static int const symmtensor_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int const symmtensor_comps[] = {
  0,1,2,   4,5,   8
};
struct tensor const TT_symmtensor = {
  "symmetric tensor T_(ij)", 3, 2,
  ARRSIZE(symmtensor_vars), symmtensor_vars, symmtensor_parity,
  ARRSIZE(symmtensor_comps), symmtensor_comps
};



/* A third rank tensor with symmetry in the first two indices */
static int const symmtensor3a_vars[] = {
#if 1
   0, 1, 2,    1, 3, 4,    2, 4, 5,
   6, 7, 8,    7, 9,10,    8,10,11,
  12,13,14,   13,15,16,   14,16,17
#else
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
   3, 4, 5,
   9,10,11,
  12,13,14,
  
   6, 7, 8,
  12,13,14,
  15,16,17
#endif
};
static int const symmtensor3a_parity[] = {
  +1,+1,+1,
  +1,+1,+1,
  +1,+1,+1,
  
  +1,+1,+1,
  +1,+1,+1,
  +1,+1,+1,
  
  +1,+1,+1,
  +1,+1,+1,
  +1,+1,+1,
};
static int const symmtensor3a_comps[] = {
#if 1
   0, 1, 2,    4, 5,    8,
   9,10,11,   13,14,   17,
  18,19,20,   22,23,   26
#else
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
  12,13,14,
  15,16,17,
  
  24,25,26
#endif
};
struct tensor const TT_symmtensor3a = {
  "symmetric tensor T_(ij)k", 3, 3,
  ARRSIZE(symmtensor3a_vars), symmtensor3a_vars, symmtensor3a_parity,
  ARRSIZE(symmtensor3a_comps), symmtensor3a_comps
};



/* A third rank tensor with symmetry in the last two indices */
static int const symmtensor3b_vars[] = {
#if 1
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
   3, 4, 5,
   9,10,11,
  12,13,14,
  
   6, 7, 8,
  12,13,14,
  15,16,17
#else
   0, 1, 2,    1, 3, 4,    2, 4, 5,
   6, 7, 8,    7, 9,10,    8,10,11,
  12,13,14,   13,15,16,   14,16,17
#endif
};
static int const symmtensor3b_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int const symmtensor3b_comps[] = {
#if 1
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
  12,13,14,
  15,16,17,
  
  24,25,26
#else
   0, 1, 2,    4, 5,    8,
   9,10,11,   13,14,   17,
  18,19,20,   22,23,   26
#endif
};
struct tensor const TT_symmtensor3b = {
  "symmetric tensor T^i_(jk)", 3, 3,
  ARRSIZE(symmtensor3b_vars), symmtensor3b_vars, symmtensor3b_parity,
  ARRSIZE(symmtensor3b_comps), symmtensor3b_comps
};



/* A third rank tensor with symmetry in the last two indices */
static int const symmtensor3c_vars[] = {
#if 1
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
   3, 4, 5,
   9,10,11,
  12,13,14,
  
   6, 7, 8,
  12,13,14,
  15,16,17
#else
   0, 1, 2,    1, 3, 4,    2, 4, 5,
   6, 7, 8,    7, 9,10,    8,10,11,
  12,13,14,   13,15,16,   14,16,17
#endif
};
static int const symmtensor3c_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int const symmtensor3c_comps[] = {
#if 1
   0, 1, 2,
   3, 4, 5,
   6, 7, 8,
  
  12,13,14,
  15,16,17,
  
  24,25,26
#else
   0, 1, 2,    4, 5,    8,
   9,10,11,   13,14,   17,
  18,19,20,   22,23,   26
#endif
};
struct tensor const TT_symmtensor3c = {
  "symmetric tensor T_i(jk)", 3, 3,
  ARRSIZE(symmtensor3c_vars), symmtensor3c_vars, symmtensor3c_parity,
  ARRSIZE(symmtensor3c_comps), symmtensor3c_comps
};



/* A fourth rank tensor with symmetries both in its first and last two
   indices */
static int const symmtensor4_vars[] = {
   0, 1, 2,    1, 3, 4,    2, 4, 5,
   6, 7, 8,    7, 9,10,    8,10,11,
  12,13,14,   13,15,16,   14,16,17,
  
   6, 7, 8,    7, 9,10,    8,10,11,
  18,19,20,   19,21,22,   20,22,23,
  24,25,26,   25,27,28,   26,28,29,
  
  12,13,14,   13,15,16,   14,16,17,
  24,25,26,   25,27,28,   26,28,29,
  30,31,32,   31,33,34,   32,34,35
};
static int const symmtensor4_parity[] = {
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1,
  +1,+1,+1,   +1,+1,+1,   +1,+1,+1
};
static int const symmtensor4_comps[] = {
   0, 1, 2,    4, 5,    8,
   9,10,11,   13,14,   17,
  18,19,20,   22,23,   26,
  
  36,37,38,   40,41,   44,
  45,46,47,   49,50,   53,
  
  72,73,74,   76,77,   80
};
struct tensor const TT_symmtensor4 = {
  "symmetric tensor T_(ij)(kl)", 3, 4,
  ARRSIZE(symmtensor4_vars), symmtensor4_vars, symmtensor4_parity,
  ARRSIZE(symmtensor4_comps), symmtensor4_comps
};



/******************************************************************************/



/* A 4-scalar */
static int const fourscalar_vars[] = {
  0
};
static int const fourscalar_parity[] = {
  +1
};
static int const fourscalar_comps[] = {
  0
};
struct tensor const TT_4scalar = {
  "4-scalar", 4, 0,
  ARRSIZE(fourscalar_vars), fourscalar_vars, fourscalar_parity,
  ARRSIZE(fourscalar_comps), fourscalar_comps
};



/* A vector */
static int const fourvector_vars[] = {
  0,1,2,3
};
static int const fourvector_parity[] = {
  +1,+1,+1,+1
};
static int const fourvector_comps[] = {
  0,1,2,3
};
struct tensor const TT_4vector = {
  "4-vector", 4, 1,
  ARRSIZE(fourvector_vars), fourvector_vars, fourvector_parity,
  ARRSIZE(fourvector_comps), fourvector_comps
};



/* A symmetric second rank tensor */
static int const foursymmtensor_vars[] = {
  0,1,2,3,   1,4,5,6,   2,5,7,8,   3,6,8,9
};
static int const foursymmtensor_parity[] = {
  +1,+1,+1,+1,   +1,+1,+1,+1,   +1,+1,+1,+1,   +1,+1,+1,+1
};
static int const foursymmtensor_comps[] = {
  0,1,2,3,   5,6,7,   10,11,   15
};
struct tensor const TT_4symmtensor = {
  "symmetric 4-tensor T_(ij)", 4, 2,
  ARRSIZE(foursymmtensor_vars), foursymmtensor_vars, foursymmtensor_parity,
  ARRSIZE(foursymmtensor_comps), foursymmtensor_comps
};



/******************************************************************************/



/* The Weyl scalars are treated as special case.  */

/* This expects a tetrad as specified in gr-qc/0104063 (Baker,
   Campanelli, Lousto: The Lazarus project: A pragmatic approach to
   binary black hole evolutions) eqns. (5.6).  In particular: l^a and
   n^a point into the future, in the sense that their temporal
   components are positive; l^a points outwards, n^a points inwards,
   Re(m) is aligned with e_theta, and Im(m) is aligned with e_phi.
   The orthonormalisation conserves the direction of Im(m).  */

/* The Weyl scalars Psi_n, stored as 5 complex numbers */
static int const weylscalars_vars[] = {
  0,1,2,3,4
};
static int const weylscalars_parity[] = {
  +1,+1,+1,+1,+1
};
static int const weylscalars_comps[] = {
  0,1,2,3,4
};
struct tensor const TT_weylscalars = {
  "Weyl scalars Psi_n", 4, 0,
  ARRSIZE(weylscalars_vars), weylscalars_vars, weylscalars_parity,
  ARRSIZE(weylscalars_comps), weylscalars_comps
};



/* The Weyl scalars Psi_n, stored as 10 real numbers */
static int const weylscalars_real_vars[] = {
  0,1,2,3,4,5,6,7,8,9
};
static int const weylscalars_real_parity[] = {
  +1,+1,+1,+1,+1,+1,+1,+1,+1,+1
};
static int const weylscalars_real_comps[] = {
  0,1,2,3,4,5,6,7,8,9
};
struct tensor const TT_weylscalars_real = {
  "Weyl scalars Psi_n", 4, 0,
  ARRSIZE(weylscalars_real_vars), weylscalars_real_vars,
  weylscalars_real_parity,
  ARRSIZE(weylscalars_real_comps), weylscalars_real_comps
};



/******************************************************************************/



struct tensor const * const TT_alltensors[] = {
  &TT_scalar,                   /* []              T          */
  &TT_vector,                   /* u               T^i        */
  &TT_tensor,                   /* ud              T^i_j      */
  &TT_symmtensor,               /* dd_sym          T_(ij)     */
  &TT_symmtensor3a,             /* dd_sym_d        T_(ij)k    */
  &TT_symmtensor3b,             /* u_dd_sym        T^i_(jk)   */
  &TT_symmtensor3c,             /* d_dd_sym        T_i(jk)    */
  &TT_symmtensor4,              /* dd_sym_dd_sym   T_(ij)(kl) */
  &TT_weylscalars,              /*                 Psi_n      */
  &TT_weylscalars_real,         /*                 Psi_n      */
};

int const TT_numtensors = ARRSIZE(TT_alltensors);



/* Ensure that all tensor declarations are internally consistent */
void
CheckTensorType (struct tensor const * restrict const atensor)
{
  int i, n;
  assert (atensor->name);
  assert (atensor->dim>=0);
  assert (atensor->rank>=0);
  assert (atensor->ncomps>=0);
  assert (atensor->ncomps == floor(pow(atensor->dim, atensor->rank) + 0.5));
  assert (atensor->nvars>=0 && atensor->nvars<=atensor->ncomps);
  assert (atensor->vars);
  for (i=0; i<atensor->ncomps; ++i) {
    assert (atensor->vars[i]>=0 && atensor->vars[i]<atensor->nvars);
    assert (abs(atensor->parity[i]) <= 1);
  }
  assert (atensor->comps);
  for (n=0; n<atensor->nvars; ++n) {
    assert (atensor->comps[n]>=0 && atensor->comps[n]<atensor->ncomps);
    assert (atensor->vars[atensor->comps[n]] == n);
  }
}
