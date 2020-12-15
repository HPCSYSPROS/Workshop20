/* $Header$ */

#ifndef TENSORTYPES_H
#define TENSORTYPES_H

#ifdef __cplusplus
#  undef restrict
#  ifdef CCTK_CXX_RESTRICT
#    define restrict CCTK_CXX_RESTRICT
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif



#include "cctk.h"



/*
 * A tensor type has a name, which is an arbitrary descriptive string.
 * It has a dimension and a rank.
 * It has dim^rank components, numbered in Fortran array order.
 * A tensor has symmetries, which means that not all of its components
 *    are stored.  The stored components are called "vars".
 * Each component is either stored as var,
 *    or is equivalent to a var times a parity (+1 or -1),
 *    or is zero, indicated by a parity of 0.
 * There is a mapping from var to the component,
 *    and another mapping back from component to var and parity.
 */



/* A tensor description */
struct tensor
{
  const char * name;            /* description */
  int dim;                      /* dimension */
  int rank;                     /* rank */
  int ncomps;                   /* dim^rank */
  const int * restrict vars;    /* map component to variable */
  const int * restrict parity;  /* parity for the above */
  int nvars;                    /* depends on symmetries */
  const int * restrict comps;   /* map variable to component */
};



/* Pre-defined tensor types */
extern struct tensor const TT_scalar;       /* []              T          */
extern struct tensor const TT_vector;       /* u               T^i        */
extern struct tensor const TT_tensor;       /* ud              T^i_j      */
extern struct tensor const TT_symmtensor;   /* dd_sym          T_(ij)     */
extern struct tensor const TT_symmtensor3a; /* dd_sym_d        T_(ij)k    */
extern struct tensor const TT_symmtensor3b; /* u_dd_sym        T^i_(jk)   */
extern struct tensor const TT_symmtensor3c; /* d_dd_sym        T_i(jk)    */
extern struct tensor const TT_symmtensor4;  /* dd_sym_dd_sym   T_(ij)(kl) */

extern struct tensor const TT_4scalar;      /* []              T          */
extern struct tensor const TT_4vector;      /* u               T^i        */
extern struct tensor const TT_4symmtensor;  /* dd_sym          T_(ij)     */

extern struct tensor const TT_weylscalars;  /*                 Psi_n      */
extern struct tensor const TT_weylscalars_real; /*             Psi_n      */

/* All pre-defined tensor types */
extern struct tensor const * const TT_alltensors[];
extern int const TT_numtensors;



/* Find the tensor type for a grid function */
struct tensor const *
TT_varindex2tensortype (int const varindex);

/* Map between variable indices and tensor type vars */
int
TT_varindex2var (int const varindex);
int
TT_var2varindex (int const old_varindex, int const var);

/* Map between tensor components and tensor indices */
void
TT_component2indices (int const dim, int const rank,
                      int component,
                      int * restrict const indices);
int
TT_indices2component (int const dim, int const rank,
                      int const * restrict const indices);

/* Convert interpolator derivative codes to the number of derivatives */
int
TT_derivcode2num_derivs (int const derivcode);

/* Convert interpolator derivative codes to the derivative indices */
void
TT_derivcode2derivs (int const derivcode, int const num_derivs,
                     int * restrict const indices);

/* Find derivatives of tensor types */
struct tensor const *
TT_derivative (struct tensor const * const tensortype, int const num_derivs);



#ifdef __cplusplus
}
#endif

#endif /* TENSORTYPES_H */
