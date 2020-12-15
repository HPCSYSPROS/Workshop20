 /*@@
   @header    cGH.h
   @date      Wed Feb 17 03:17:47 1999
   @author    Tom Goodale
   @desc 
   The cGH structure.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CGH_H_
#define _CGH_H_ 1

#include "cctk_Types.h"

typedef struct
{
  char storage;
  char comm;
} cGHGroupData;

typedef struct _cGH
{
  int cctk_dim;
  int cctk_iteration;

  /* ...[dim]*/
  int *cctk_gsh;
  int *cctk_lsh;
  int *cctk_lbnd;
  int *cctk_ubnd;

  /* allocated shape */
  int *cctk_ash;

  /* unused */
  int *cctk_to;
  int *cctk_from;
  
  /* The grid spacings */
  CCTK_REAL cctk_delta_time;
  CCTK_REAL *cctk_delta_space;

  /* FIXME we want coordinate registration instead of this */
  CCTK_REAL *cctk_origin_space;

  /* The bounding box - 1 => a real boundary, 0 => a local grid boundary. */
  /* bbox[2*dim] */
  int *cctk_bbox;

  /* The refinement factor over the top level (coarsest) grid. */
  int *cctk_levfac;

  /* Offset between this level's and the coarsest level's origin */
  int *cctk_levoff;
  int *cctk_levoffdenom;

  /* The refinement factor in time over the top level (coarsest) grid. */
  int cctk_timefac;

  /* The convergence level (numbered from zero upwards) */
  int cctk_convlevel;

  /* The (per level) convergence factor */
  int cctk_convfac;

  /* The number of ghostzones in each direction */
  int *cctk_nghostzones;

  /* The coordinate time */
  CCTK_REAL cctk_time;

  /* An identifier string for this hierarchy
    (used to construct names for output files) */
  const char *identity;

  /* data[variable_index][timelevel][ijk]
     (ijk is as calculated e.g. by CCTK_GFINDEX3D) */
  void ***data;

  /* The extension array */
  void **extensions;

  /* All the group data for this GH (storage, comm, etc. */
  cGHGroupData *GroupData;

} cGH;

#endif /* _CGH_H_ */
