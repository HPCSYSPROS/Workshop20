/*@@
  @file     CoordBase.c
  @date     23 July 2003
  @author   Gabrielle Allen and David Rideout
  @desc
            Functions which implement the coordinate API.
  @enddesc
  @version  $Id$
@@*/

#include <string.h>
#include <stdlib.h>

#include "cctk.h"
#include "util_Table.h"
#include "util_Hash.h"
#include "util_String.h"

#define COORD_IN_COORDBASE 1
#include "CoordBase.h"
#undef COORD_IN_COORDBASE

#include "coordbaseGH.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_CoordBase_CoordBase_c);

/********************************************************************
 *********************  Aliased Routine Prototypes  *****************
 ********************************************************************/

/********************************************************************
 *********************  Local Data  *********************************
 ********************************************************************/

static unsigned int longest_coordname = 0;
static int longest_systemname = 0;

/********************************************************************
 *********************  (to be Aliased?) Routines  ******************
 ********************************************************************/

/*@@
  @routine    CoordBase_SystemRegister
  @date       24 July 2003
  @author     Gabrielle Allen and David Rideout
  @desc
              Registers a coordinate system.
  @enddesc
  @calls      CCTK_GHExtension, Util_HashData, CCTK_VWarn, Util_TableCreate,
              Util_TableSetString, Util_TableSetInt, Util_HashAdd

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        dim
  @vdesc      coordinate system dimension
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        systemname
  @vdesc      coordinate system name
  @vtype      CCTK_STRING
  @vio        in
  @endvar

  @returntype int
  @returndesc
  coordinate system handle, or:
  COORDERROR_INVALIDDIM   invalid dimension passed in
  COORDERROR_INVALIDNAME  invalid name passed in
  COORDERROR_TABLEERROR   error from key-value or hash tables in flesh
  COORDERROR_SYSTEMEXISTS coordinate system of this name already exists
  @endreturndesc
@@*/

CCTK_INT CoordBase_SystemRegister(CCTK_POINTER_TO_CONST GH, CCTK_INT dim,
                                  CCTK_STRING systemname) {
  int len, retval;
  int *ptr;
  const coordbaseGH *GHex;
  uHash *hash;

  /* Initialize variables */
  len = retval = 0;
  hash = NULL;

  /* Check input arguments */
  if (dim <= 0) {
    retval = COORDERROR_INVALIDDIM;
  } else if (!systemname) {
    retval = COORDERROR_INVALIDNAME;
  } else {
    /* Get coordinate system name / handle hash table from GH Extension */
    GHex = CCTK_GHExtension(GH, "CoordBase");
    hash = GHex->coordsystems;

    /* Has a coordinate system of this name already been registered? */
    len = strlen(systemname);
    if (len > longest_systemname) {
      longest_systemname = len + 1;
    }
    if (Util_HashData(hash, len, systemname, 0)) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Coordinate system %s already registered", systemname);
      retval = COORDERROR_SYSTEMEXISTS;
    }
  }

  if (retval == 0) {
    /* Now make a table for this system */
    ptr = malloc(sizeof(int));
    retval = *ptr = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (retval < 0) {
      CCTK_WARN(0, "Could not create table");
      retval = COORDERROR_TABLEERROR;
    }

    /* Add coordinate system information */
    if (Util_TableSetString(*ptr, systemname, "NAME") < 0 ||
        Util_TableSetInt(*ptr, dim, "DIMENSION") < 0) {
      CCTK_WARN(0, "Internal error when adding key/value table entries");
      retval = COORDERROR_TABLEERROR;
    }

    /* Register the system with the Coordinate GH Extension */
    if (Util_HashAdd(hash, len, systemname, 0, ptr) < 0) {
      CCTK_WARN(0, "Internal error when storing coordinate system handle");
      retval = COORDERROR_TABLEERROR;
    }
  }

  return retval;
}

/*@@
  @routine    CoordBase_SystemHandle
  @date       24 July 2003
  @author     Gabrielle Allen and David Rideout
  @desc
              Returns the coordinate system handle for a given name.
  @enddesc
  @calls      CCTK_GHExtension, Util_HashData

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        systemname
  @vdesc      coordinate system name
  @vtype      CCTK_STRING
  @vio        in
  @endvar

  @returntype int
  @returndesc
  coordinate system handle, or:
  COORDERROR_TABLEERROR   error from hash table
  COORDERROR_NOSYSTEM     no coordinate system of this name is registered
  @endreturndesc
@@*/

CCTK_INT CoordBase_SystemHandle(CCTK_POINTER_TO_CONST GH,
                                CCTK_STRING systemname) {
  int *handle_ptr, retval;
  const coordbaseGH *GHex;
  uHash *hash;

  /* Get coordinate system name / handle hash table from GH Extension */
  GHex = CCTK_GHExtension(GH, "CoordBase");
  hash = GHex->coordsystems;

  /* Get coordinate handle from hash table */
  handle_ptr = Util_HashData(hash, strlen(systemname), systemname, 0);
  if (!handle_ptr) {
    retval = COORDERROR_NOSYSTEM;
  } else {
    retval = *handle_ptr >= 0 ? *handle_ptr : COORDERROR_TABLEERROR;
  }

  return retval;
}

/*@@
  @routine    CoordBase_CoordRegister
  @date       24 July 2003
  @author     Gabrielle Allen and David Rideout
  @desc
              Registers a coordinate within a coordinate system.
  @enddesc
  @calls      Util_TableGetInt, CCTK_VWarn, Util_TableGetIntArray,
              Util_TableSetIntArray, Util_TableGetString, Util_StrCmpi,
              Util_TableCreate, Util_TableSetInt, Util_TableSetString

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        systemhandle
  @vdesc      coordinate system handle
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        direction
  @vdesc      index within coodinate basis
  @vtype      CCTK_INT
  @vio        in
  @vcomment   ranges from 1 to the dimension of the coordinate system
  @endvar
  @var        coordname
  @vdesc      coordinate system name
  @vtype      CCTK_STRING
  @vio        in
  @endvar

  @returntype int
  @returndesc
  coordinate handle, or:
  COORDERROR_INVALIDDIM        invalid 'direction'
  COORDERROR_INVALIDHANDLE     invalid handle passed in / coordinate system
                                 does not exist
  COORDERROR_TABLEERROR        error from hash or key-value tables in flesh
  COORDERROR_COORDINATEEXISTS  coordinate already exists for this 'direction'
  COORDERROR_DUPLICATENAME     coordinate of this name already exists in
                                 this system
  @endreturndesc
@@*/

CCTK_INT CoordBase_CoordRegister(CCTK_POINTER_TO_CONST GH,
                                 CCTK_INT systemhandle, CCTK_INT direction,
                                 CCTK_STRING coordname) {
  int coordinate_handle;
  CCTK_INT system_dimension;
  CCTK_INT *coordinate_tables;
  CCTK_INT i;
  char *coordname_buffer;

  (void)(&GH + 0);

  coordinate_handle = 0;

  /* Check input arguments */
  /* Get coordinate system dimension, and check for valid systemhandle */
  if (Util_TableGetInt(systemhandle, &system_dimension, "DIMENSION") < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Invalid system handle "
                                                        "%d",
               (int)systemhandle);
    coordinate_handle = COORDERROR_INVALIDHANDLE;
  }

  /* Check for valid 'direction' */
  if (direction < 1 || direction > system_dimension) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid 'direction'.  "
               "It must be an integer from 1 to the coordinate system "
               "dimension %d.",
               (int)system_dimension);
    coordinate_handle = COORDERROR_INVALIDDIM;
  }

  /* Get coordinate tables from system table */
  coordinate_tables = malloc(system_dimension * sizeof(CCTK_INT));
  if (Util_TableGetIntArray(systemhandle, system_dimension, coordinate_tables,
                            "COORDINATES") < 0) {
    /* The only possibility is that this is a bad key, so add this entry */
    for (i = 0; i < system_dimension; ++i) {
      coordinate_tables[i] = -1;
    }
    if (Util_TableSetIntArray(systemhandle, system_dimension, coordinate_tables,
                              "COORDINATES") < 0) {
      CCTK_WARN(0, "Error adding coordinate tables to system table");
      coordinate_handle = COORDERROR_TABLEERROR;
    }
  }

  /* Has a coordinate been registered for this 'direction'? */
  if (coordinate_tables[direction - 1] >= 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Coordinate already "
               "registered for direction %d",
               (int)direction);
    coordinate_handle = COORDERROR_COORDINATEEXISTS;
  }

  /* Has this coordinate name been used already in this system? */
  coordname_buffer = malloc(longest_coordname * sizeof(char));
  for (i = 0; i < system_dimension; ++i) {
    if (coordinate_tables[i] >= 0) {
      if (Util_TableGetString(coordinate_tables[i], longest_coordname,
                              coordname_buffer, "NAME") < 0) {
        CCTK_WARN(0, "Error reading coordinate name from its table");
        coordinate_handle = COORDERROR_TABLEERROR;
      } else if (Util_StrCmpi(coordname_buffer, coordname) == 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Coordinate '%s' "
                   "is already registered for this system",
                   coordname_buffer);
        coordinate_handle = COORDERROR_DUPLICATENAME;
        break;
      }
    }
  }

  /* Everything checks out */
  if (coordinate_handle == 0) {
    /* Is this coordinate name the longest yet? */
    if (strlen(coordname) + 1 > longest_coordname) {
      longest_coordname = strlen(coordname) + 1;
    }

    /* Create coordinate table */
    coordinate_handle = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (coordinate_handle < 0) {
      CCTK_WARN(0, "Could not create coordinate table");
      coordinate_handle = COORDERROR_TABLEERROR;
    }
    if (Util_TableSetInt(coordinate_handle, systemhandle, "SYSTEM") < 0 ||
        Util_TableSetString(coordinate_handle, coordname, "NAME") < 0 ||
        Util_TableSetInt(coordinate_handle, direction, "DIRECTION") < 0) {
      CCTK_WARN(0, "Error writing to coordinate table");
      coordinate_handle = COORDERROR_TABLEERROR;
    }

    /* Register coordinate table with system table */
    coordinate_tables[direction - 1] = coordinate_handle;
    if (Util_TableSetIntArray(systemhandle, system_dimension, coordinate_tables,
                              "COORDINATES") < 0) {
      CCTK_WARN(0, "Error writing coordinate table to system table");
      coordinate_handle = COORDERROR_TABLEERROR;
    }
  }

  free(coordinate_tables);
  free(coordname_buffer);
  return coordinate_handle;
}

/*@@
  @routine    CoordBase_CoordHandle
  @date       25 July 2003
  @author     Gabrielle Allen and David Rideout
  @desc
              Returns the coordinate handle for a given coordinate name and
              system name
  @enddesc
  @calls      CCTK_GHExtension, Util_HashData, CCTK_WARN, Util_TableGetInt,
              Util_TableGetIntArray, Util_TableGetString, Util_StrCmpi

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        coordname
  @vdesc      coordinate name
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @var        systemname
  @vdesc      coordinate system name
  @vtype      CCTK_STRING
  @vio        in
  @endvar

  @returntype int
  @returndesc
  coordinate table handle, or:
  COORDERROR_NOSYSTEM     no coordinate system of this name is registered
  COORDERROR_TABLEERROR   error from hash table
  COORDERROR_NOSUCHCOORD  no coordinate of the name is registered for this
                            system
  @endreturndesc
@@*/

CCTK_INT CoordBase_CoordHandle(CCTK_POINTER_TO_CONST GH, CCTK_STRING coordname,
                               CCTK_STRING systemname) {
  int coord_handle;
  int *system_handle_ptr;
  CCTK_INT system_dimension;
  CCTK_INT *coordinate_table_handles;
  int coordname_length;
  int i;
  const coordbaseGH *GHex;
  uHash *hash;
  char *coordname_buffer;

  /* Get coordinate system name / handle hash table from GH Extension */
  GHex = CCTK_GHExtension(GH, "CoordBase");
  hash = GHex->coordsystems;

  coord_handle = COORDERROR_NOSYSTEM;

  /* Get system table handle */
  system_handle_ptr = Util_HashData(hash, strlen(systemname), systemname, 0);
  if (!system_handle_ptr) {
    CCTK_WARN(1, "No such coordinate system registered");
  } else {
    /* Get system dimension */
    if (Util_TableGetInt(*system_handle_ptr, &system_dimension, "DIMENSION") <
        0) {
      CCTK_WARN(0, "Error reading dimension from system table");
      coord_handle = COORDERROR_TABLEERROR;
    }

    /* Get coordinate table handle */
    coordinate_table_handles = malloc(system_dimension * sizeof(CCTK_INT));
    if (Util_TableGetIntArray(*system_handle_ptr, system_dimension,
                              coordinate_table_handles, "COORDINATES") < 0) {
      CCTK_WARN(0, "Error reading coordinate tables from system table");
      coord_handle = COORDERROR_TABLEERROR;
    }

    /* Search for coordname */
    coordname_length = strlen(coordname) + 1;
    coordname_buffer = malloc(coordname_length * sizeof(char));
    for (i = 0; i < system_dimension; ++i) {
      if (coordinate_table_handles[i] >= 0) {
        if (Util_TableGetString(coordinate_table_handles[i], coordname_length,
                                coordname_buffer, "NAME") < 0) {
          CCTK_WARN(0, "Error reading coordinate name");
          coord_handle = COORDERROR_TABLEERROR;
        } else if (Util_StrCmpi(coordname_buffer, coordname) == 0) {
          coord_handle = coordinate_table_handles[i];
          break;
        }
      }
      if (i == system_dimension - 1) {
        coord_handle = COORDERROR_NOSUCHCOORD;
      }
    }

    free(coordinate_table_handles);
    free(coordname_buffer);
  }
  return coord_handle;
}

/*@@
  @routine    CoordBase_GroupSystem
  @date       28 July 2003
  @author     David Rideout
  @desc
              Returns the coordinate system handle associated with a group of
              grid variables
  @enddesc
  @calls      CCTK_GroupTagsTable, CCTK_VWarn, Util_TableGetString,
              CoordBase_SystemHandle, CCTK_GHExtension, CCTK_GroupIndex,
              CCTK_GroupDimI

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        groupname
  @vdesc      variable group name
  @vtype      CCTK_STRING
  @vio        in
  @endvar

  @returntype int
  @returndesc
  coordinate system handle, or:
  COORDERROR_INVALIDGROUPNAME no such group exists
  COORDERROR_NOCOORDSYS       no coordinate system is associated with the
                                group
  @endreturndesc
@@*/

CCTK_INT CoordBase_GroupSystem(CCTK_POINTER_TO_CONST GH,
                               CCTK_STRING groupname) {
  int system_handle;
  int group_tags_table;
  int gindex; /* group index */
  int gdim;   /* group dimension */
  const coordbaseGH *GHex;
  char *system_name;

  system_handle = -9999;

  /* Get tags table associated with this group */
  group_tags_table = CCTK_GroupTagsTable(groupname);
  if (group_tags_table < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group name '%s'", groupname);
    system_handle = COORDERROR_INVALIDGROUPNAME;
  } else {
    /* Get associated coordinate system */
    system_name = malloc(longest_systemname * sizeof(char));
    if (Util_TableGetString(group_tags_table, longest_systemname, system_name,
                            "COORDSYSTEM") >= 0) {
      system_handle = CoordBase_SystemHandle(GH, system_name);
      if (system_handle < 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Coordinate system "
                   "'%s' associated with group '%s' does not exist",
                   system_name, groupname);
        system_handle = COORDERROR_INVALIDNAME;
      }
    }

    /* Check for a default for all GF variables of this dimension */
    gindex = CCTK_GroupIndex(groupname);
    if (gindex < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, "CoordBase",
                 "group index is out of bounds");
    if (system_handle == -9999 && CCTK_GroupTypeI(gindex) == CCTK_GF) {
      GHex = CCTK_GHExtension(GH, "CoordBase");
      gdim = CCTK_GroupDimI(gindex);
      system_handle = GHex->default_coord_systems[gdim - 1];
      if (system_handle < 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "No coordinate "
                   "system is associated with group '%s'",
                   groupname);
        system_handle = COORDERROR_NOCOORDSYS;
      }
    }
    free(system_name);
  }

  return system_handle;
}

/*@@
  @routine    CoordBase_SetDefaultSystem
  @date       28 July 2003
  @author     David Rideout
  @desc
              Set this coordinate system to be the default for grid functions
              of the same dimension
  @enddesc
  @calls      CCTK_GHExtension, CoordBase_SystemHandle, CCTK_VWarn,
              Util_TableGetInt

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        systemname
  @vdesc      coordinate system name
  @vtype      CCTK_STRING
  @vio        in
  @endvar

  @returntype int
  @returndesc
  coordinate system handle, or
  COORDERROR_INVALIDNAME   no coordinate system of this name has been
                             registered
  COORDERROR_NODIMENSION   coordinate system does not have a valid dimension
  COORDERROR_DEFAULTEXISTS grid variables of this dimension already have a
                             default coordinate system registered
  @endreturndesc
@@*/

CCTK_INT CoordBase_SetDefaultSystem(CCTK_POINTER_TO_CONST GH,
                                    CCTK_STRING systemname) {
  int system_handle;
  CCTK_INT system_dimension;
  int retval;
  const coordbaseGH *GHex;
#ifdef DEBUG
  int i;
#endif

  /* Retrieve array from GH extension */
  GHex = CCTK_GHExtension(GH, "CoordBase");

  /* Get coordinate system handle */
  system_handle = CoordBase_SystemHandle(GH, systemname);
  if (system_handle < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Invalid system name "
                                                        "'%s'",
               systemname);
    retval = COORDERROR_INVALIDNAME;
  } else {
    /* Get system dimension */
    if (Util_TableGetInt(system_handle, &system_dimension, "DIMENSION") < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "System '%s' does not "
                 "have a registered dimension",
                 systemname);
      retval = COORDERROR_NODIMENSION;
    } else {
      /* Set the default for this dimension */
      if (GHex->default_coord_systems[system_dimension - 1] >= 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Default coordinate "
                   "system for %d dimensional grid variables is already set.  "
                   "Overwriting to '%s'.",
                   (int)system_dimension, systemname);
        retval = COORDERROR_DEFAULTEXISTS;
      } else {
        retval = system_handle;
      }
      GHex->default_coord_systems[system_dimension - 1] = system_handle;
    }
  }

#ifdef DEBUG
  for (i = 0; i < system_dimension; ++i) {
    printf("def_coord_sys[%d] = %d\n", i, GHex->default_coord_systems[i]);
  }
#endif

  return retval;
}

/*@@
  @routine    CoordBase_GetDefaultSystem
  @date       9 May 2004
  @author     Thomas Radke
  @desc
              Get the default coordinate system for grid variables
              of a given dimension
  @enddesc

  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        systemdim
  @vdesc      dimension of coordinate system to query
  @vtype      CCTK_INT
  @vio        in
  @vcomment   ranges from 1 to the maximum number of dimensions
  @endvar

  @returntype int
  @returndesc
  coordinate system handle, or
  COORDERROR_INVALIDDIM   given dimension is invalid
  COORDERROR_NOSYSTEM     given dimension does not have a default coordinate
                          system associated
  @endreturndesc
@@*/

CCTK_INT CoordBase_GetDefaultSystem(CCTK_POINTER_TO_CONST GH,
                                    CCTK_INT systemdim) {
  CCTK_INT retval;
  const coordbaseGH *myGH;

  myGH = CCTK_GHExtension(GH, "CoordBase");

  retval = COORDERROR_INVALIDDIM;
  if (systemdim >= 1 && systemdim <= CCTK_MaxDim()) {
    retval = myGH->default_coord_systems[systemdim - 1];
  }

  return (retval);
}
