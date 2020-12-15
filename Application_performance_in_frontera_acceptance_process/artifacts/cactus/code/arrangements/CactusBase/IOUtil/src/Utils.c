/*@@
   @file      Utils.c
   @date      Tue 4th July 2000
   @author    Gabrielle Allen
   @desc
              Utility routines which may be called by other I/O thorns.
   @enddesc
   @version   $Id$
 @@*/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "util_Table.h"
#include "ioGH.h"
#include "ioutil_Utils.h"
#include "ioutil_CheckpointRecovery.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_IOUtil_Utils_c)

/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
/* uncomment this if you want some debugging output */
/* #define DEBUG_IOUTIL 1 */

/********************************************************************
 ********************    Internal Typedefs   ************************
 ********************************************************************/
typedef struct {
  const cGH *GH;
  ioRequest **request_list;
  const char *method_name;
  const char *parameter_name;
  int out_every_default;
  CCTK_REAL out_dt_default;
} info_t;

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void SetOutputVar(int vindex, const char *optstring, void *arg);

/*@@
  @routine    IOUtil_ParseVarsForOutput
  @date       Fri 26 April 2002
  @author     Thomas Radke
  @desc
              Parses the given 'out_vars' string as a list of variable and/or
              groupnames. If an option string is appended to a name this will
              be parsed for I/O request options.
  @enddesc

  @calls      CCTK_TraverseString

  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        method_name
  @vdesc      name of the originating I/O method
  @vtype      const char *
  @vio        in
  @endvar
  @var        parameter_name
  @vdesc      name of the originating I/O parameter to parse
  @vtype      const char *
  @vio        in
  @endvar
  @var        stop_on_parse_errors
  @vdesc      flag indicating whether stop with a level-0 warning message
              on parsing errors
  @vtype      int
  @vio        in
  @endvar
  @var        out_vars
  @vdesc      string with list of variables and/or group names to parse
  @vtype      const char *
  @vio        in
  @endvar
  @var        out_every_default
  @vdesc      default output frequency to set the I/O request to
  @vtype      int
  @vio        in
  @endvar
  @var        out_dt_default
  @vdesc      default output frequency to set the I/O request to
  @vtype      CCTK_REAL
  @vio        in
  @endvar
  @var        request_list
  @vdesc      list of I/O requests to fill out
  @vtype      ioRequest *[]
  @vio        out
  @endvar
@@*/
void IOUtil_ParseVarsForOutput(const cGH *GH, const char *method_name,
                               const char *parameter_name,
                               int stop_on_parse_errors, const char *out_vars,
                               int out_every_default, CCTK_REAL out_dt_default,
                               ioRequest *request_list[]) {
  int i;
  info_t info;

  /* free current list of I/O requests */
  for (i = CCTK_NumVars() - 1; i >= 0; i--) {
    if (request_list[i]) {
      IOUtil_FreeIORequest(&request_list[i]);
    }
  }

  /* generate new list of I/O requests */
  info.GH = GH;
  info.request_list = request_list;
  info.method_name = method_name;
  info.parameter_name = parameter_name;
  info.out_every_default = out_every_default;
  info.out_dt_default = out_dt_default;
  if (CCTK_TraverseString(out_vars, SetOutputVar, &info, CCTK_GROUP_OR_VAR) <
      0) {
    CCTK_VWarn(stop_on_parse_errors ? 0 : 1, __LINE__, __FILE__,
               CCTK_THORNSTRING, "error while parsing parameter '%s'",
               parameter_name);
  }
}

/*@@
  @routine    IOUtil_FreeIORequest
  @date       Fri 26 April 2002
  @author     Thomas Radke
  @desc
              Frees an I/O request description.
  @enddesc

  @var        request
  @vdesc      pointer to the I/O request structure
  @vtype      ioRequest **
  @vio        inout
  @endvar
@@*/
void IOUtil_FreeIORequest(ioRequest **request) {
  free((*request)->reductions);
  free((*request)->vectors);
  free(*request);
  *request = NULL;
}

/*@@
  @routine    IOUtil_ParseOutputFrequency
  @date       Mon 27 May 2002
  @author     Thomas Radke
  @desc
              Parses the option string for a given variable for the
              'out_every' and 'out_dt' options and sets its output
              frequency.
  @enddesc

  @calls      Util_TableCreateFromString

  @var        method_name
  @vdesc      name of the originating I/O method
  @vtype      const char *
  @vio        in
  @endvar
  @var        parameter_name
  @vdesc      name of the originating I/O parameter to parse
  @vtype      const char *
  @vio        in
  @endvar
  @var        stop_on_parse_errors
  @vdesc      flag indicating whether to stop with level-0 warning message
              on a parsing error
  @vtype      int
  @vio        in
  @endvar
  @var        vindex
  @vdesc      index of the variable to set the output frequency
  @vtype      int
  @vio        in
  @endvar
  @var        optstring
  @vdesc      option string to parse
  @vtype      const char *
  @vio        in
  @endvar
  @var        out_every_ptr
  @vdesc      pointer to the variable's output frequency flag
  @vtype      int *
  @vio        out
  @endvar
  @var        out_dt_ptr
  @vdesc      pointer to the variable's output frequency flag
  @vtype      CCTK_REAL *
  @vio        out
  @endvar
@@*/
void IOUtil_ParseOutputFrequency(const char *method_name,
                                 const char *parameter_name,
                                 int stop_on_parse_errors, int vindex,
                                 const char *optstring, CCTK_INT *out_every_ptr,
                                 CCTK_REAL *out_dt_ptr) {
  int table, iterator, warnlevel;
  char key[128];
  char *fullname;
  CCTK_INT type, nelems;
  DECLARE_CCTK_PARAMETERS

  warnlevel = stop_on_parse_errors ? 0 : 1;
  fullname = CCTK_FullName(vindex);
  table = Util_TableCreateFromString(optstring);
  if (table >= 0) {
    if (Util_TableQueryValueInfo(table, &type, &nelems, "out_every") > 0) {
      if (type == CCTK_VARIABLE_INT && nelems == 1) {
        Util_TableGetInt(table, out_every_ptr, "out_every");
      } else {
        CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'out_every' in option string "
                   "'%s' in parameter '%s' (must be an integer)",
                   optstring, parameter_name);
        CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option will be ignored for %s output of variable '%s'",
                   method_name, fullname);
      }
      Util_TableDeleteKey(table, "out_every");
    }

    if (Util_TableQueryValueInfo(table, &type, &nelems, "out_dt") > 0) {
      if (type == CCTK_VARIABLE_REAL && nelems == 1) {
        Util_TableGetReal(table, out_dt_ptr, "out_dt");
      } else {
        CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'out_dt' in option string "
                   "'%s' in parameter '%s' (must be real)",
                   optstring, parameter_name);
        CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option will be ignored for %s output of variable '%s'",
                   method_name, fullname);
      }
      Util_TableDeleteKey(table, "out_dt");
    }

    /* warn about other options */
    iterator = Util_TableItCreate(table);
    for (iterator = Util_TableItCreate(table);
         Util_TableItQueryIsNonNull(iterator) > 0 &&
             Util_TableItQueryKeyValueInfo(iterator, sizeof(key), key, 0, 0) >
                 0;
         Util_TableItAdvance(iterator)) {
      CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Found option with unrecognized key '%s' in option string "
                 "'%s' in parameter '%s'",
                 key, optstring, parameter_name);
      CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Option will be ignored for %s output of variable '%s'",
                 method_name, fullname);
    }
    Util_TableItDestroy(iterator);

    Util_TableDestroy(table);
  } else {
    CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't parse option string '%s' in parameter '%s'", optstring,
               parameter_name);
    CCTK_VWarn(warnlevel, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Option string will be ignored for %s output of variable '%s'",
               method_name, fullname);
  }
  free(fullname);
}

/*@@
   @routine    IOUtil_1DLines
   @date       July 4 2000
   @author     Gabrielle Allen, Gerd Lanfermann, Thomas Radke
   @desc
               Fills out an array determining where to position 1D lines
               for output on a multidimensional regular Cartesian unigrid.
               The first slot of the array specifies the 1D line direction,
               the second slot fixes the indices of the starting point
               of that line on the grid.
   @enddesc

   @calls      Coord_GetDefaultSystem

   @var        GH
   @vdesc      Pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_dims
   @vdesc      number of dimensions of underlying grid
   @vtype      int
   @vio        in
   @endvar
   @var        origin_index
   @vdesc      origin for 1D lines given by grid point indices
   @vtype      const int [num_dims][num_dims]
   @vio        in
   @endvar
   @var        origin_phys
   @vdesc      origin for 1D lines given by physical coordinates
   @vtype      const CCTK_REAL [num_dims][num_dims]
   @vio        in
   @endvar
   @var        slice_center
   @vdesc      resulting 1D line slice center
   @vtype      int [num_dims][num_dims]
   @vio        out
   @endvar

   @returntype int
   @returndesc
               0  for success
               -1 no coordinate system of given dimensions found
   @endreturndesc
 @@*/
int IOUtil_1DLines(const cGH *GH, int num_dims, int *const *const origin_index,
                   CCTK_REAL *const *const origin_phys,
                   int *const *slice_center) {
  int dim, dir, coord_system_handle = -1, have_coords, have_origin_index, len;
  char *coord_system_name;
  CCTK_INT *coord_handles = NULL;
  CCTK_REAL *lower, *upper, *delta;

  /* get the default coordinate system associated with this grid dimension */
  coord_system_name = NULL;
  have_coords = CCTK_IsFunctionAliased("Coord_GetDefaultSystem");
  if (have_coords) {
    coord_system_handle = Coord_GetDefaultSystem(GH, num_dims);
    coord_handles = malloc(num_dims * sizeof(CCTK_INT));
    have_coords = coord_system_handle >= 0 &&
                  Util_TableGetIntArray(coord_system_handle, num_dims,
                                        coord_handles, "COORDINATES") >= 0;
  }
  if (have_coords) {
    len = Util_TableGetString(coord_system_handle, 0, NULL, "NAME");
    have_coords = len > 0;
    if (have_coords) {
      coord_system_name = malloc(len);
      Util_TableGetString(coord_system_handle, len, coord_system_name, "NAME");
      have_coords = Util_StrCmpi(coord_system_name, "cart");
    }
  }

  /* check that origin_index[] is valid */
  have_origin_index = origin_index != NULL;
  for (dir = 0; dir < num_dims && have_origin_index; dir++) {
    for (dim = 0; dim < num_dims && have_origin_index; dim++) {
      have_origin_index = origin_index[dir][dim] >= 0;
    }
  }

  /* if no coordinate info is available then origin_index[] must be given */
  if (!have_coords && !have_origin_index) {
    CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
               "IOUtil_1DLines: Found no default Cartesian coordinate system "
               "associated with grid variables of dimension %d, and no slice "
               "center index coordinates were given either - slice center "
               "will not be set up for output of 1D lines from %dD variables",
               num_dims, num_dims);
    for (dir = 0; dir < num_dims; dir++) {
      memset(slice_center[dir], 0, num_dims * sizeof(int));
    }

    free(coord_handles);
    free(coord_system_name);

    return (-1);
  }

  /* get the ranges in every direction */
  lower = calloc(3 * num_dims, sizeof(CCTK_REAL));
  upper = lower + 1 * num_dims;
  delta = lower + 2 * num_dims;
  if (!have_origin_index) {
    for (dir = 0; dir < num_dims; dir++) {
      if (Util_TableGetReal(coord_handles[dir], &lower[dir], "COMPMIN") < 0 ||
          Util_TableGetReal(coord_handles[dir], &upper[dir], "COMPMAX") < 0 ||
          Util_TableGetReal(coord_handles[dir], &delta[dir], "DELTA") < 0) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "IOUtil_1DLines: Could not get ranges for %c-direction "
                   "of associated default coordinate system '%s'",
                   'x' + dir, coord_system_name);
      }
    }
  }

  /* now set the slice center for each line
     according to origin_index[] or origin_phys[] */
  for (dir = 0; dir < num_dims; dir++) {
    for (dim = 0; dim < num_dims; dim++) {
      if (dim == dir) {
        /* line always starts at the first point */
        slice_center[dir][dim] = 0;
      } else if (origin_index[dir][dim] >= 0) {
        /* FIXME: check upper index bounds also ?? */
        slice_center[dir][dim] = origin_index[dir][dim];
      } else if (lower[dim] > origin_phys[dir][dim] ||
                 upper[dim] < origin_phys[dir][dim]) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "IOUtil_1DLines: %c-coordinate (%f) for slice center of "
                   "1D lines in %c-direction for %dD variables is out of "
                   "grid coordinates range (%f, %f)",
                   'x' + dim, (double)origin_phys[dir][dim], 'x' + dir,
                   num_dims, (double)lower[dim], (double)upper[dim]);

        slice_center[dir][dim] =
            origin_index[dir][dim] == -1
                ? 0
                : ceil((upper[dim] - lower[dim]) / (2 * delta[dim]) - 1e-6);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "IOUtil_1DLines: slice center will default to %c-index %d",
                   'x' + dir, slice_center[dir][dim]);
      } else {
        /* Find index for first point above the chosen coordinate */
        slice_center[dir][dim] =
            ceil((origin_phys[dir][dim] - lower[dim]) / delta[dim] - 1e-6);

#ifdef DEBUG_IOUTIL
        printf("spxyz for %c-coord of lines in %c-direction is %d\n", 'x' + dim,
               'x' + dir, slice_center[dir][dim]);
#endif
      }
    }
  }

  free(lower);
  free(coord_handles);
  free(coord_system_name);

  return (0);
}

/*@@
   @routine    IOUtil_2DPlanes
   @date       July 4 2000
   @author     Gabrielle Allen, Gerd Lanfermann, Thomas Radke
   @desc
               Fills out an array determining where to position 2D planes
               for output on a multidimensional regular Cartesian unigrid.
   @enddesc

   @calls      Coord_GetDefaultSystem

   @var        GH
   @vdesc      Pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_dims
   @vdesc      number of dimensions of underlying grid
   @vtype      int
   @vio        in
   @endvar
   @var        origin_index
   @vdesc      origin for 2D planes given by grid point indices
   @vtype      const int [num_dims]
   @vio        in
   @endvar
   @var        origin_phys
   @vdesc      origin for 2D planes given by physical coordinates
   @vtype      const CCTK_REAL [num_dims]
   @vio        in
   @endvar
   @var        slice_center
   @vdesc      resulting 2D plane slice center
   @vtype      int [num_dims]
   @vio        out
   @endvar

   @returntype int
   @returndesc
               0  for success
               -1 no coordinate system of given dimensions found
   @endreturndesc
 @@*/
int IOUtil_2DPlanes(const cGH *GH, int num_dims, const int *origin_index,
                    const CCTK_REAL *origin_phys, int *slice_center) {
  int have_coords, have_origin_index, len, coord_system_handle = -1, dir;
  char *coord_system_name;
  CCTK_INT *coord_handles = NULL;
  CCTK_REAL *lower, *upper, *delta;

  /* get the default coordinate system associated with this grid dimension */
  coord_system_name = NULL;
  have_coords = CCTK_IsFunctionAliased("Coord_GetDefaultSystem");
  if (have_coords) {
    coord_system_handle = Coord_GetDefaultSystem(GH, num_dims);
    coord_handles = malloc(num_dims * sizeof(CCTK_INT));
    have_coords = coord_system_handle >= 0 &&
                  Util_TableGetIntArray(coord_system_handle, num_dims,
                                        coord_handles, "COORDINATES") >= 0;
  }
  if (have_coords) {
    len = Util_TableGetString(coord_system_handle, 0, NULL, "NAME");
    have_coords = len > 0;
    if (have_coords) {
      coord_system_name = malloc(len);
      Util_TableGetString(coord_system_handle, len, coord_system_name, "NAME");
      have_coords = Util_StrCmpi(coord_system_name, "cart");
    }
  }

  /* check that origin_index[] is valid */
  have_origin_index = origin_index != NULL;
  for (dir = 0; dir < num_dims && have_origin_index; dir++) {
    have_origin_index = origin_index[dir] >= 0;
  }

  /* if no coordinate system is available then origin_index[] must be given */
  if (!have_coords && !have_origin_index) {
    CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
               "IOUtil_2DPlanes: Found no default Cartesian coordinate system "
               "associated with grid variables of dimension %d, and no slice "
               "center index coordinates were given either - slice center "
               "will not be set up for output of 2D planes from %dD variables",
               num_dims, num_dims);

    free(coord_handles);
    free(coord_system_name);

    return (-1);
  }

  /* get the ranges in every direction */
  lower = calloc(3 * num_dims, sizeof(CCTK_REAL));
  upper = lower + 1 * num_dims;
  delta = lower + 2 * num_dims;
  if (have_coords) {
    for (dir = 0; dir < num_dims; dir++) {
      if (Util_TableGetReal(coord_handles[dir], &lower[dir], "COMPMIN") < 0 ||
          Util_TableGetReal(coord_handles[dir], &upper[dir], "COMPMAX") < 0 ||
          Util_TableGetReal(coord_handles[dir], &delta[dir], "DELTA") < 0) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "IOUtil_2DPlanes: Could not get ranges/deltas for "
                   "%c-direction of associated coordinate system '%s'",
                   'x' + (num_dims - dir - 1), coord_system_name);
      }
    }
  }

  /* now set the slice center for each line
     according to origin_index[] or origin_phys[] */
  for (dir = 0; dir < num_dims; dir++) {
    if (origin_index && origin_index[dir] >= 0) {
      slice_center[dir] = origin_index[dir];
    } else if (lower[num_dims - 1 - dir] > origin_phys[dir] ||
               upper[num_dims - 1 - dir] < origin_phys[dir]) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "IOUtil_2DPlanes: %c-coordinate for slice center of 2D "
                 "planes (%f) is out of grid coordinates range (%f, %f)",
                 'x' + (num_dims - dir - 1), (double)origin_phys[dir],
                 (double)lower[num_dims - 1 - dir],
                 (double)upper[num_dims - 1 - dir]);

      slice_center[dir] =
          origin_index[dir] == -1
              ? 0
              : ceil((upper[num_dims - 1 - dir] - lower[num_dims - 1 - dir]) /
                         (2 * delta[num_dims - 1 - dir]) -
                     1e-6);
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "IOUtil_2DPlanes: slice center will default to %c-index %d",
                 'x' + dir, slice_center[dir]);
    } else {
      /* Find index for first point above the chosen coordinate */
      slice_center[dir] = ceil((origin_phys[dir] - lower[num_dims - 1 - dir]) /
                                   delta[num_dims - 1 - dir] -
                               1e-6);
#ifdef DEBUG_IOUTIL
      printf("sp2xyz for planes perpendicular to %d-direction is %d\n", dir,
             slice_center[dir]);
#endif
    }
  }

  free(lower);
  free(coord_handles);
  free(coord_system_name);

  return (0);
}

/*@@
  @routine    IOUtil_PrintTimings
  @date       Wed Jun 28 2000
  @author     Thomas Radke
  @desc
              Gets the timing information for the given timers and prints it
              as INFO messages to screen.
  @enddesc
  @history
  @endhistory
  @var        description
  @vdesc      description of the timers
  @vtype      const char *
  @vio        in
  @endvar
  @var        ntimers
  @vdesc      number of timers passed in
  @vtype      int
  @vio        in
  @endvar
  @var        timers
  @vdesc      array of timers
  @vtype      const int [ntimers]
  @vio        in
  @endvar
  @var        timer_descriptions
  @vdesc      array of timer descriptions
  @vtype      const char *const [ntimers]
  @vio        in
  @endvar
@@*/

void IOUtil_PrintTimings(const char *description, int ntimers,
                         const int *timers,
                         const char *const *const timer_descriptions) {
  int i, j;
  cTimerData *info;

  info = CCTK_TimerCreateData();
  if (info) {
    CCTK_INFO(description);

    for (i = 0; i < info->n_vals; i++) {
      for (j = 0; j < ntimers; j++) {
        CCTK_TimerI(timers[j], info);
        if (j == 0) {
          CCTK_VInfo(CCTK_THORNSTRING, "  %s:", info->vals[i].heading);
        }
        switch (info->vals[i].type) {
        case val_int:
          CCTK_VInfo(CCTK_THORNSTRING, "    %s %5d %s", timer_descriptions[j],
                     info->vals[i].val.i, info->vals[i].units);
          break;

        case val_long:
          CCTK_VInfo(CCTK_THORNSTRING, "    %s %5d %s", timer_descriptions[j],
                     (int)info->vals[i].val.l, info->vals[i].units);
          break;

        case val_double:
          CCTK_VInfo(CCTK_THORNSTRING, "    %s %5.1f %s", timer_descriptions[j],
                     info->vals[i].val.d, info->vals[i].units);
          break;

        default:
          CCTK_WARN(1, "Unknown data type for timer info");
          break;
        }
      }
    }
    CCTK_INFO("-----------------------------------------");
    CCTK_TimerDestroyData(info);
  } else {
    CCTK_WARN(1, "Couldn't create timer info structure ! "
                 "No timing output available.");
  }
}

/*@@
  @routine    IOUtil_CreateDirectory
  @date       Fri 10 Aug 2001
  @author     Thomas Radke
  @desc
              Creates an output directory path and makes sure it is visible
              on all I/O processors.
              It is assumed that processor 0 is always an I/O processor.
              If there are other I/O processors, they will also try and
              create the directory themselfs. This guarantees the directory
              be created on all nodes in case the I/O processors don't share
              a common filesystem.
  @enddesc
  @calls      CCTK_MyProc
              CCTK_CreateDirectory
              CCTK_Barrier

  @var        GH
  @vdesc      pointer to the GH extensions
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        dirname
  @vdesc      the directory to create
  @vtype      const char *
  @vio        in
  @endvar
  @var        multiple_io_procs
  @vdesc      flag indicating that there will be I/O from multiple I/O procs
  @vtype      int
  @vio        in
  @endvar
  @var        ioproc
  @vdesc      I/O processor associated with this processor
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              0 for non-I/O processors, or
              return code of @seeroutine CCTK_CreateDirectory
  @endreturndesc
@@*/
int IOUtil_CreateDirectory(const cGH *GH, const char *dirname,
                           int multiple_io_procs, int ioproc) {
  int myproc, retval;

  /* default return value for non-I/O processors */
  retval = 0;

  /* first, processor 0 creates the directory */
  myproc = CCTK_MyProc(GH);
  if (myproc == 0) {
    retval = CCTK_CreateDirectory(0755, dirname);
  }

  if (multiple_io_procs) {
    /* now the other I/O processors create the directory
       after syncing with processor 0 */
    CCTK_Barrier(GH);
    if (myproc == ioproc || ioproc != 0) {
      retval = CCTK_CreateDirectory(0755, dirname);
    }
  }

  return (retval);
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void SetOutputVar(int vindex, const char *optstring, void *arg) {
  int table, iterator;
  char key[128];
  CCTK_INT type, nelems;
  CCTK_INT *refinement_levels;
  char *fullname;
  const info_t *info = arg;
  ioRequest *request;
  DECLARE_CCTK_PARAMETERS

  /* allocate a new I/O request structure and initialize it with defaults */
  request = IOUtil_DefaultIORequest(info->GH, vindex, info->out_every_default,
                                    info->out_dt_default);

  if (!optstring) {
    info->request_list[vindex] = request;
    return;
  }

  fullname = CCTK_FullName(vindex);

  table = Util_TableCreateFromString(optstring);
  if (table >= 0) {
    /* check for option 'out_every' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "out_every") > 0) {
      if (type == CCTK_VARIABLE_INT && nelems == 1) {
        Util_TableGetInt(table, &request->out_every, "out_every");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'out_every' in option string "
                   "'%s' in parameter '%s' (must be an integer)",
                   optstring, info->parameter_name);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'out_every' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "out_every");
    }

    /* check for option 'out_dt' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "out_dt") > 0) {
      if (type == CCTK_VARIABLE_REAL && nelems == 1) {
        Util_TableGetReal(table, &request->out_dt, "out_dt");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'out_dt' in option string "
                   "'%s' in parameter '%s' (must be real)",
                   optstring, info->parameter_name);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'out_dt' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "out_dt");
    }

    /* check for boolean option 'out_unchunked' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "out_unchunked") > 0) {
      int unchunked = -1;
      if (type == CCTK_VARIABLE_CHAR && nelems >= (int)sizeof("no") - 1 &&
          nelems <= (int)sizeof("false") - 1) {
        char value[sizeof("false")];
        Util_TableGetString(table, sizeof(value), value, "out_unchunked");
        if (CCTK_Equals(value, "yes") || CCTK_Equals(value, "true")) {
          unchunked = 1;
        } else if (CCTK_Equals(value, "no") || CCTK_Equals(value, "false")) {
          unchunked = 0;
        }
      }
      if (unchunked >= 0) {
        request->out_unchunked = unchunked;
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'out_unchunked' in option string "
                   "'%s' in parameter '%s' (must be a boolean)",
                   optstring, info->parameter_name);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'out_unchunked' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "out_unchunked");
    }

    /* check for hyperslab option 'direction' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "direction") > 0) {
      if (type == CCTK_VARIABLE_INT && nelems > 0 &&
          nelems <= request->vdim * request->vdim &&
          nelems % request->vdim == 0) {
        request->hdim = nelems / request->vdim;
        Util_TableGetIntArray(table, nelems, request->direction, "direction");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'direction' in option string "
                   "'%s' in parameter '%s' (must be an integer array with "
                   "%d x hdim elements, 1 <= hdim <= %d)",
                   optstring, info->parameter_name, (int)request->vdim,
                   (int)request->vdim);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'direction' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "direction");
    }

    /* check for hyperslab option 'origin' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "origin") > 0) {
      if (type == CCTK_VARIABLE_INT && nelems == request->vdim) {
        Util_TableGetIntArray(table, nelems, request->origin, "origin");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'origin' in option string "
                   "'%s' in parameter '%s' (must be an integer array with "
                   "%d elements)",
                   optstring, info->parameter_name, (int)request->vdim);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'origin' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "origin");
    }

    /* check for hyperslab option 'extent' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "extent") > 0) {
      if (type == CCTK_VARIABLE_INT && nelems == request->hdim) {
        Util_TableGetIntArray(table, nelems, request->extent, "extent");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'extent' in option string "
                   "'%s' in parameter '%s' (must be an integer array with "
                   "%d elements)",
                   optstring, info->parameter_name, (int)request->hdim);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'extent' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "extent");
    }

    /* check for hyperslab option 'downsample' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "downsample") > 0) {
      if (type == CCTK_VARIABLE_INT && nelems == request->hdim) {
        Util_TableGetIntArray(table, nelems, request->downsample, "downsample");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'downsample' in option string "
                   "'%s' in parameter '%s' (must be an integer array with "
                   "%d elements)",
                   optstring, info->parameter_name, (int)request->hdim);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'downsample' will be ignored for %s output of "
                   "variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "downsample");
    }

    /* check for option 'refinement_levels' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "refinement_levels") >
        0) {
      request->refinement_levels = 0;
      if (type == CCTK_VARIABLE_INT && nelems > 0) {
        refinement_levels = malloc(nelems * sizeof(int));
        Util_TableGetIntArray(table, nelems, refinement_levels,
                              "refinement_levels");
        while (--nelems >= 0) {
          if (refinement_levels[nelems] < 0) {
            request->refinement_levels = 0;
            break;
          }
          request->refinement_levels |= 1 << refinement_levels[nelems];
        }
        free(refinement_levels);
      }
      if (request->refinement_levels == 0) {
        if (CCTK_GroupTypeFromVarI(vindex) == CCTK_GF) {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Invalid value for option 'refinement_levels' in option "
                     "string '%s' in parameter '%s' (must be an array of "
                     "positive integers)",
                     optstring, info->parameter_name);
        } else {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Meaningless option 'refinement_levels' in option "
                     "string '%s' in parameter '%s' (only makes sense for "
                     "grid functions)",
                     optstring, info->parameter_name);
        }
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'refinement_levels' will be ignored for %s output "
                   "of variable '%s'",
                   info->method_name, fullname);
        request->refinement_levels = -1;
      }
      Util_TableDeleteKey(table, "refinement_levels");
    }

    /* check for option 'compression_level' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "compression_level") >
        0) {
      request->compression_level = -1;
      if (type == CCTK_VARIABLE_INT && nelems == 1) {
        Util_TableGetInt(table, &request->compression_level,
                         "compression_level");
        if (request->compression_level < 0 || request->compression_level > 9) {
          type = -1;
          request->compression_level = -1;
        }
      }
      if (!(type == CCTK_VARIABLE_INT && nelems == 1)) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'compression_level' in option "
                   "string '%s' in parameter '%s' (must be an integer in the "
                   "range [0-9])",
                   optstring, info->parameter_name);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'compression_level' will be ignored for %s output "
                   "of variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "compression_level");
    }

    /* check for option 'reductions' */
    if (Util_TableQueryValueInfo(table, &type, &nelems, "reductions") > 0) {
      request->reductions = NULL;
      if (type == CCTK_VARIABLE_CHAR && nelems >= 0) {
        request->reductions = malloc(nelems + 1);
        Util_TableGetString(table, nelems + 1, request->reductions,
                            "reductions");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid value for option 'reductions' in option "
                   "string '%s' in parameter '%s' (must be a string)",
                   optstring, info->parameter_name);
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Option 'reductions' will be ignored for %s output "
                   "of variable '%s'",
                   info->method_name, fullname);
      }
      Util_TableDeleteKey(table, "reductions");
    }

    /* warn about unrecognized options */
    iterator = Util_TableItCreate(table);
    for (iterator = Util_TableItCreate(table);
         Util_TableItQueryIsNonNull(iterator) > 0 &&
             Util_TableItQueryKeyValueInfo(iterator, sizeof(key), key, 0, 0) >
                 0;
         Util_TableItAdvance(iterator)) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Found option with unrecognized key '%s' in option string "
                 "'%s' in parameter '%s'",
                 key, optstring, info->parameter_name);
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Option will be ignored for %s output of variable '%s'",
                 info->method_name, fullname);
    }
    Util_TableItDestroy(iterator);

    Util_TableDestroy(table);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't parse option string '%s' in parameter '%s'", optstring,
               info->parameter_name);
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Option string will be ignored for %s output of variable '%s'",
               info->method_name, fullname);
  }

  /* assign the new I/O request */
  info->request_list[vindex] = request;

  /* clean up */
  free(fullname);
}

/* return the default I/O request description structure for a variable */
ioRequest *IOUtil_DefaultIORequest(const cGH *GH, int vindex,
                                   int out_every_default,
                                   CCTK_REAL out_dt_default) {
  int vdim;
  int *extent_int;
  ioRequest *request;
  const ioGH *myGH;
  DECLARE_CCTK_PARAMETERS

  /* check for valid variable index */
  vdim = CCTK_GroupDimFromVarI(vindex);
  if (vdim < 0) {
    return (NULL);
  }

  myGH = CCTK_GHExtension(GH, "IO");

  /* allocate a new I/O request structure */
  request = malloc(sizeof(ioRequest));

  /* fill out the basics */
  request->vindex = vindex;
  request->timelevel = 0;
  request->check_exist = myGH->recovered;
  request->out_every = out_every_default;
  request->out_dt = out_dt_default;
  request->out_unchunked = out_unchunked;
  request->reductions = NULL;
  request->with_ghostzones = 0;
  request->refinement_levels = -1;
  request->compression_level = -1;

  /* get the I/O request datatype (will be single-precision if requested) */
  request->hdatatype = CCTK_VarTypeI(vindex);
  if (myGH->out_single && CCTK_GroupTypeFromVarI(vindex) != CCTK_SCALAR) {
    if (request->hdatatype == CCTK_VARIABLE_REAL) {
      request->hdatatype = CCTK_VARIABLE_REAL4;
    } else if (request->hdatatype == CCTK_VARIABLE_COMPLEX) {
      request->hdatatype = CCTK_VARIABLE_COMPLEX8;
    }
#ifdef CCTK_INT2
    else if (request->hdatatype == CCTK_VARIABLE_INT) {
      request->hdatatype = CCTK_VARIABLE_INT2;
    }
#endif
  }

  /* get the variable's dimension and extents */
  request->vdim = vdim;
  extent_int = NULL;
  if (request->vdim > 0) {
    extent_int = malloc(request->vdim * sizeof(int));
    CCTK_GroupgshVI(GH, request->vdim, extent_int, vindex);
  }

  /* allocate the arrays all in one go
     only initialize those which are mandatory for the Hyperslab API */
  request->vectors =
      calloc((request->vdim + 6) * request->vdim + 1, sizeof(CCTK_INT));
  request->hoffset = request->vectors + 0 * request->vdim;
  request->hsize = request->vectors + 1 * request->vdim;
  request->hsize_chunk = request->vectors + 2 * request->vdim;
  request->origin = request->vectors + 3 * request->vdim + 1;
  request->extent = request->vectors + 4 * request->vdim + 1;
  request->downsample = request->vectors + 5 * request->vdim + 1;
  request->direction = request->vectors + 6 * request->vdim + 1;

  for (request->hdim = 0; request->hdim < request->vdim; request->hdim++) {
    request->extent[request->hdim] = extent_int[request->hdim];
    request->direction[request->hdim * (request->vdim + 1)] = 1;

    /* take the downsampling parameters from IOUtil */
    if (request->hdim == 0) {
      request->downsample[request->hdim] = out_downsample_x;
    } else if (request->hdim == 1) {
      request->downsample[request->hdim] = out_downsample_y;
    } else if (request->hdim == 2) {
      request->downsample[request->hdim] = out_downsample_z;
    } else {
      request->downsample[request->hdim] = 1;
    }
  }

  /* clean up */
  if (extent_int) {
    free(extent_int);
  }

  return (request);
}
