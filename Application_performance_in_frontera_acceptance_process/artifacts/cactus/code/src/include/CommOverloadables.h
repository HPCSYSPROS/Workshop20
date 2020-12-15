 /*@@
   @header    CommOverloadables.h
   @date      Thu Feb  4 08:11:41 1999
   @author    Tom Goodale
   @desc
              The overloadable functions for the comm layer.
              See OverloadMacros.h to see how to use these.
   @enddesc
   @version   $Header$
 @@*/

#include "cctk_Flesh.h"
#include "cctk_GroupsOnGH.h"

#ifdef OVERLOADABLE_CALL
#undef OVERLOADABLE_CALL
#endif

#ifdef OVERLOABLE_PREFIX
#undef OVERLOADABLE_PREFIX
#endif

#ifdef OVERLOABLE_DUMMY_PREFIX
#undef OVERLOADABLE_DUMMY_PREFIX
#endif

#define OVERLOADABLE_CALL CCTK_
#define OVERLOADABLE_PREFIX CCTK_
#define OVERLOADABLE_DUMMY_PREFIX CCTKi_Dummy

#ifdef ARGUMENTS
#undef ARGUMENTS
#endif

#ifdef USE_ARGUMENTS
#undef USE_ARGUMENTS
#endif

#ifdef RETURN_TYPE
#undef RETURN_TYPE
#endif

#ifdef ATTRIBUTES
#undef ATTRIBUTES
#endif

#define RETURN_TYPE int
#define ARGUMENTS const cGH *GH, const char *group
#define USE_ARGUMENTS GH = GH; group = group;
#define ATTRIBUTES
OVERLOADABLE(SyncGroup)

OVERLOADABLE(EnableGroupStorage)
OVERLOADABLE(DisableGroupStorage)

OVERLOADABLE(EnableGroupComm)
OVERLOADABLE(DisableGroupComm)

#undef ARGUMENTS
#define ARGUMENTS const cGH *GH, int num_groups, const int *groups, const int *directions
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH; num_groups = num_groups; groups = groups; directions = directions;
OVERLOADABLE(SyncGroupsByDirI)

#undef ARGUMENTS
#define ARGUMENTS const cGH *GH
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH;
OVERLOADABLE(Barrier)
OVERLOADABLE(MyProc)
OVERLOADABLE(nProcs)

#undef ARGUMENTS
#define ARGUMENTS cGH *GH
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH;
OVERLOADABLE(ParallelInit)

#undef ARGUMENTS
#define ARGUMENTS cGH *GH, int retval
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH; retval = retval;
#undef RETURN_TYPE
#define RETURN_TYPE int
#undef ATTRIBUTES
#define ATTRIBUTES CCTK_ATTRIBUTE_NORETURN
OVERLOADABLE(Exit)
OVERLOADABLE(Abort)

#undef ARGUMENTS
#define ARGUMENTS tFleshConfig *config, int convergence_level
#undef USE_ARGUMENTS
#define USE_ARGUMENTS config = config; convergence_level = convergence_level;
#undef RETURN_TYPE
#define RETURN_TYPE cGH *
#undef ATTRIBUTES
#define ATTRIBUTES
OVERLOADABLE(SetupGH)

#undef ARGUMENTS
#define ARGUMENTS const cGH *GH, int dir, int group, const char *groupname
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH; dir = dir; group = group; groupname = groupname;
#undef RETURN_TYPE
#define RETURN_TYPE const int *
OVERLOADABLE(ArrayGroupSizeB)

#undef ARGUMENTS
#define ARGUMENTS const cGH *GH, int group, const char *groupname
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH; group = group; groupname = groupname;
#undef RETURN_TYPE
#define RETURN_TYPE int
OVERLOADABLE(QueryGroupStorageB)

#undef ARGUMENTS
#define ARGUMENTS const cGH *GH, int group, cGroupDynamicData *data
#undef USE_ARGUMENTS
#define USE_ARGUMENTS GH = GH; group = group; data = data;
#undef RETURN_TYPE
#define RETURN_TYPE int
OVERLOADABLE(GroupDynamicData)

#undef ARGUMENTS
#undef USE_ARGUMENTS
#undef RETURN_TYPE

#define RETURN_TYPE int
#define ARGUMENTS const cGH *GH, int n_groups,const int *groups,const int *timelevels, int *status
#define USE_ARGUMENTS GH = GH; n_groups=n_groups; groups = groups; timelevels = timelevels; status = status;
OVERLOADABLE(GroupStorageIncrease)
OVERLOADABLE(GroupStorageDecrease)

#undef ARGUMENTS
#undef USE_ARGUMENTS
#undef RETURN_TYPE

#define RETURN_TYPE int
#define ARGUMENTS const cGH *GH, int n_groups,const int *groups, int *status
#define USE_ARGUMENTS GH = GH; n_groups=n_groups; groups = groups; status = status;
OVERLOADABLE(QueryMaxTimeLevels)
#undef ARGUMENTS
#undef USE_ARGUMENTS
#undef RETURN_TYPE

/* overloadable routine CCTK_InterpGridArrays() */
#define RETURN_TYPE int
#define ARGUMENTS const cGH *GH,                                              \
                  int N_dims,                                                 \
                  int local_interp_handle,                                    \
                  int param_table_handle,                                     \
                  int coord_system_handle,                                    \
                  int N_interp_points,                                        \
                    int interp_coords_type,                                   \
                    const void *const interp_coords[],                        \
                  int N_input_arrays,                                         \
                    const CCTK_INT input_array_indices[],                     \
                  int N_output_arrays,                                        \
                    const CCTK_INT output_array_types[],                      \
                    void *const output_arrays[]
#define USE_ARGUMENTS (void) (GH + 0);                                        \
                      (void) (N_dims + 0);                                    \
                      (void) (local_interp_handle + 0);                       \
                      (void) (param_table_handle + 0);                        \
                      (void) (coord_system_handle + 0);                       \
                      (void) (N_interp_points + 0);                           \
                      (void) (interp_coords_type + 0);                        \
                      (void) (interp_coords + 0);                             \
                      (void) (N_input_arrays + 0);                            \
                      (void) (input_array_indices + 0);                       \
                      (void) (N_output_arrays + 0);                           \
                      (void) (output_array_types + 0);                        \
                      (void) (output_arrays + 0);
OVERLOADABLE(InterpGridArrays)

#undef ARGUMENTS
#undef USE_ARGUMENTS
#undef RETURN_TYPE
#undef ATTRIBUTES

#undef OVERLOADABLE_CALL
#undef OVERLOADABLE_PREFIX
#undef OVERLOADABLE_DUMMY_PREFIX
