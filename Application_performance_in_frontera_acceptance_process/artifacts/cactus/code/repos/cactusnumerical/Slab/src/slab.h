#ifndef SLAB_H
#define SLAB_H

#ifdef __cplusplus 
extern "C" 
{
#endif

#include <stdio.h>

#include <cctk.h>

/*
 * Slab_Transfer copies a slab from one array into a slab of another
 * array.
 * 
 * The src and dst variables describe the source and the destination
 * slab, respectively.
 * 
 * The variables gsh, lbnd, ash and lsh describe the shape and
 * distribution of the array containing the slab.  They are equivalent
 * to the corresponding quantities in the cGH structure.
 * 
 * off, str, and len describe the location and shape of the slab
 * within the array.  off is the offset, i.e., the location of the
 * "first" corner of the slab.  str is the stride, i.e., the distance
 * between to grid points in the slab.  The stride can be negative.
 * len is the length, i.e., the number of grid points making up the
 * slab.  len does not include the grid points that are skipped if the
 * stride is larger than one.
 * 
 * xpose describes a possible permutation of the coordinate axes
 * between the slabs.  It is source-axis = xpose[destination-axis].
 * 
 * flip describes a possible inversion of the coordinate axes (from
 * the point of view of the destination slab).  It is source-axis =
 * xpose[flip[destination-axis]].
 * 
 * The corresponding lengths of the source and the destination slabs
 * must be equal, i.e., for all d: src.len[xpose[d]] = dst.len[d].
 * 
 * The slabs are copied according to
 * 
 *     dst[dst.off + I * dst.str] = src[src.off + J * src.str]
 * 
 * where the multi-indices I and J have the ranges specified by
 * dst.len and src.len, respectively, and I and J are related by the
 * transposition
 * 
 *     J = xpose[flip[I]]
 * 
 * 
 * 
 * Restrictions:
 * 
 *     dim >= 0
 * 
 *     gsh >= 0
 *     lbnd >= 0
 *     lsh >= 0
 *     ash >= lsh
 *     lbnd + lsh <= gsh
 *     lbbox and ubbox must be booleans, i.e., either 0 or 1
 *     nghostzones >= 0
 * 
 *     len >= 0
 *     str != 0
 *     off >= 0
 *     off < gsh
 *     off + (len-1) * str >= 0
 *     off + (len-1) * str < gsh
 * 
 *     xpose must be a permutation of 0 ... dim-1
 *     flip must be a boolean, i.e., either 0 or 1
 * 
 * The source and the destination arrays may be the same.
 */

#define SLAB_MAXDIM 3

struct slabinfo {
  int gsh;
  int lbnd, lsh;
  int ash;
  int lbbox, ubbox, nghostzones;
  int off, str, len;
};

struct xferinfo {
  struct slabinfo src, dst;
  int xpose;
  int flip;
};

void
print_slabinfo (FILE                  *          const out,
                struct slabinfo const * restrict const slabinfo);

void
print_xferinfo (FILE                  *          const out,
                struct xferinfo const * restrict const xferinfo);



struct slabsetup;

struct slabsetup *
Slab_MultiTransfer_Init
(cGH             const * restrict const cctkGH,
 int                              const dim,
 struct xferinfo const * restrict const xferinfo,
 int                              const options);

int
Slab_MultiTransfer_Apply
(cGH              const                  * restrict const cctkGH,
 struct slabsetup const                  * restrict const slabsetup,
 int                                                const nvars,
 int              const                  * restrict const srctypes,
 void             const * restrict const * restrict const srcptrs,
 int              const                  * restrict const dsttypes,
 void                   * restrict const * restrict const dstptrs);

int
Slab_MultiTransfer_Finalize
(cGH              const * restrict const cctkGH,
 struct slabsetup       * restrict const slabsetup);

int
Slab_MultiTransfer
(cGH             const                  * restrict const cctkGH,
 int                                               const dim,
 struct xferinfo const                  * restrict const xferinfo,
 int                                               const options,
 int                                               const nvars,
 int             const                  * restrict const srctypes,
 void            const * restrict const * restrict const srcptrs,
 int             const                  * restrict const dsttypes,
 void                  * restrict const * restrict const dstptrs);

int
Slab_Transfer
(cGH             const * restrict const cctkGH,
 int                              const dim,
 struct xferinfo const * restrict const xferinfo,
 int                              const options,
 int                              const srctype,
 void            const * restrict const srcptr,
 int                              const dsttype,
 void                  * restrict const dstptr);

#ifdef __cplusplus 
}
#endif

#endif /* defined SLAB_H */
