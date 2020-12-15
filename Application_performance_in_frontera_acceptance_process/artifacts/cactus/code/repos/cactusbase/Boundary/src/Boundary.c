/*@@
  @file      Boundary.c
  @date      Sat Oct 26 22:39:40 CEST 2002
  @author    David Rideout
  @desc
             Implements the new boundary specification.
  @enddesc
  @version   $Header$
@@*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include "util_String.h"
#include "Boundary.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_Boundary_Boundary_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/
/* #define DEBUG 1 */

/* Linked list, called a var list, for holding variables selected for a bc:
 * Entries are sorted in the order they appear in struct BCVAR,
 * i.e. the var index varies more rapidly that table handle.
 * (Currently no sorting is done on faces specification.)
 *
 * There will be one such linked list for each type of boundary
 * condition selected (i.e. one for each bc_name).
 */
struct BCVAR {
  struct BCVAR *next; /* pointer to next entry in list */
  int faces;          /* set of faces for this application of bc */
  int width;          /* width of the boundary, if it is equal for all faces */
  int table;          /* table handle holding extra arguments */
  int var;            /* index of grid variable to which to apply the bc */
};

/*
 * Linked list of var lists, one for each type of requested bc
 * (i.e. one for each bc_name).
 *
 * Here is also recorded how many of each bc type have
 * been selected so far, so that the GetSelectedBCs doesn't have to be
 * run through twice; once simply to get the number of selected vars.
 *
 * This list is sorted by bc_name.  Alternatively one could sort it by
 * associated function pointer, but this seems a bit obtuse.
 */
struct BCDATA {
  struct BCDATA *next;    /* pointer to next element of this list */
  struct BCVAR *var_list; /* pointer to first element of a var list */
  const char *bc_name;    /* name of bc */
  int num;                /* number of entries for this bc in var list */
};

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int entry_greater_than(struct BCVAR *new, struct BCVAR *current);
#ifdef DEBUG
static void print_selections_database(void);
static void print_selected_faces(void);
#endif

/********************************************************************
 ***************** Aliased Routine Prototypes ***********************
 ********************************************************************/

CCTK_INT Bdry_Boundary_RegisterPhysicalBC(CCTK_POINTER_TO_CONST _GH,
                                          phys_bc_fn_ptr fn_pointer,
                                          CCTK_STRING bc_name);
CCTK_INT Bdry_Boundary_SelectVarForBC(CCTK_POINTER_TO_CONST _GH, CCTK_INT faces,
                                      CCTK_INT boundary_width,
                                      CCTK_INT table_handle,
                                      CCTK_STRING var_name,
                                      CCTK_STRING bc_name);
CCTK_INT Bdry_Boundary_SelectVarForBCI(CCTK_POINTER_TO_CONST _GH,
                                       CCTK_INT faces, CCTK_INT boundary_width,
                                       CCTK_INT table_handle,
                                       CCTK_INT var_index, CCTK_STRING bc_name);
CCTK_INT Bdry_Boundary_SelectGroupForBC(CCTK_POINTER_TO_CONST _GH,
                                        CCTK_INT faces, CCTK_INT boundary_width,
                                        CCTK_INT table_handle,
                                        CCTK_STRING group_name,
                                        CCTK_STRING bc_name);
CCTK_INT
Bdry_Boundary_SelectGroupForBCI(CCTK_POINTER_TO_CONST _GH, CCTK_INT faces,
                                CCTK_INT boundary_width, CCTK_INT table_handle,
                                CCTK_INT group_index, CCTK_STRING bc_name);
CCTK_INT Bdry_Boundary_SelectedGVs(CCTK_POINTER_TO_CONST _GH,
                                   CCTK_INT array_size, CCTK_INT *var_indices,
                                   CCTK_INT *faces, CCTK_INT *boundary_widths,
                                   CCTK_INT *table_handles,
                                   CCTK_STRING bc_name);

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void Boundary_ApplyPhysicalBCs(CCTK_ARGUMENTS);
void Boundary_ClearSelection(CCTK_ARGUMENTS);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/* Table for holding function pointers associated with each boundary condition:
 * This table has
 *        key   = boundary condition name (eg "Radiation")
 *        value = a CCTK_FPOINTER pointing to a function to implement that BC
 */
static int physbc_table_handle = -1;

/* Linked list for storing data associated with selections list itself */
static struct BCDATA *bcdata_list = NULL;

/* Array of (number of variables) faces specifications, for checking
   for duplicate bc selection */
static CCTK_INT *selected_faces = NULL;
static int num_cctk_vars = 0;

/* 'The' GH, i.e. to check that there is not more than one... */
static CCTK_POINTER_TO_CONST theGH = NULL;

/********************************************************************
 *********************     Aliased Routines   **********************
 ********************************************************************/

/*@@
  @routine    Bdry_Boundary_RegisterPhysicalBC
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Used to register physical boundary conditions with the boundary
              thorn.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        fn_pointer
  @vdesc      pointer to function which implements boundary condition bc_name
  @vtype      phys_bc_fn_ptr
  @vio        in
  @endvar
  @var        bc_name
  @vdesc      name of boundary condition
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
              0 success
             -1 error creating table to hold bcs
             -2 duplicate registration of bc_name
             -3 error adding bc to table
  @endreturndesc
@@*/
CCTK_INT Bdry_Boundary_RegisterPhysicalBC(CCTK_POINTER_TO_CONST _GH,
                                          phys_bc_fn_ptr fn_pointer,
                                          CCTK_STRING bc_name) {
  int retval;
  const cGH *GH = _GH;

#ifdef DEBUG
  printf("Boundary_RegisterPhysicalBC: called with GH=%p\n", GH);
#endif

  /* Check to see if this is a new GH */
  if (!theGH) /* This is the first valid GH passed to a Boundary routine */
  {
    theGH = GH;
  } else if (GH != theGH) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "New GH passed to Boundary_RegisterPhysicalBC.  "
               "Thorn CactusBase/Boundary does not yet handle multiple GHs "
               "properly.");
  }

  /* Check input arguments */
  if (!fn_pointer) {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Null pointer passed to Boundary_RegisterPhysicalBC.  "
               "Is this intentional?");
  }

  /* Check if NULL has been passed for fn_pointer */
  if (!fn_pointer) {
    /* Use dummy function if NULL function registered (e.g. for
       non-local physical bcs) */
    fn_pointer = (phys_bc_fn_ptr)&BndNone;
  }

  /* Create the registered routines table if necessary */
  if (physbc_table_handle == -1) {
    physbc_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (physbc_table_handle < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error creating table to hold registered physical boundary "
                 "conditions");
      retval = -1;
    }
  }

  /* Check if boundary condition has already been registered under this name */
  if (Util_TableGetFnPointer(physbc_table_handle, NULL, bc_name) >= 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There already exists a physical boundary condition "
               "registered under the name \"%s\"",
               bc_name);
    retval = -2;
  } else {
    /* Add boundary condition to table */
    if (Util_TableSetFnPointer(physbc_table_handle, (CCTK_FPOINTER)fn_pointer,
                               bc_name) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error adding boundary condition to table");
      retval = -3;
    }
    retval = 0;
  }

  return retval;
}

/*@@
  @routine    Bdry_Boundary_SelectVarForBC
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Used to select a Cactus variable to have boundary
              conditions applied, using its variable name.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        faces
  @vdesc      set of faces to which to apply the boundary condition
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        width
  @vdesc      if >=0, width of boundary in all directions
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        table_handle
  @vdesc      handle of table which holds arguments to be passed to bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        var_name
  @vdesc      name of variable to which to apply bc
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @var        bc_name
  @vdesc      name of bc to apply
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
   0 success
  -11 invalid variable name
  or the returncode of @seeroutine Bdry_Boundary_SelectVarForBCI
  @endreturndesc
@@*/
CCTK_INT Bdry_Boundary_SelectVarForBC(CCTK_POINTER_TO_CONST _GH, CCTK_INT faces,
                                      CCTK_INT width, CCTK_INT table_handle,
                                      CCTK_STRING var_name,
                                      CCTK_STRING bc_name) {
  int retval, var_index;
  const cGH *GH = _GH;

#ifdef DEBUG
  printf("Boundary_SelectVarForBC:\n");
  printf("  called with faces=%d, width=%d, table_handle=%d, var_name=%s, "
         "bc_name=%s\n",
         faces, width, table_handle, var_name, bc_name);
#endif

  retval = 0;

  var_index = CCTK_VarIndex(var_name);
  if (var_index < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable name");
    retval = -11;
  } else {
    retval = Bdry_Boundary_SelectVarForBCI(GH, faces, width, table_handle,
                                           var_index, bc_name);
  }

  return retval;
}

/*@@
  @routine    Bdry_Boundary_SelectVarForBCI
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Used to select a Cactus variable to have boundary
              conditions applied, using its var index.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        faces
  @vdesc      set of faces to which to apply the boundary condition
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        width
  @vdesc      if >=0, width of boundary in all directions
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        table_handle
  @vdesc      handle of table which holds arguments to be passed to bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        var_index
  @vdesc      index of variable to which to apply bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        bc_name
  @vdesc      name of bc to apply
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
   0 success
  -2 no such physical boundary condition registered
  -3 faces already selected for bc
  -5 new value for GH passed in
  -7 invalid variable index
  @endreturndesc
@@*/
CCTK_INT Bdry_Boundary_SelectVarForBCI(CCTK_POINTER_TO_CONST _GH,
                                       CCTK_INT faces, CCTK_INT width,
                                       CCTK_INT table_handle,
                                       CCTK_INT var_index,
                                       CCTK_STRING bc_name) {
  int retval;
  struct BCVAR *new_entry;
  struct BCVAR *current;
  struct BCVAR *previous;
  struct BCDATA *current_bcdata;
  struct BCDATA *previous_bcdata;
  struct BCDATA *new_bcdata;
  const cGH *GH = _GH;

  retval = 0;
  current = NULL;
  previous = NULL;
  previous_bcdata = NULL;

  if (var_index < 0 || var_index >= CCTK_NumVars()) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable index");
    retval = -7;
    return retval;
  }

#ifdef DEBUG
  printf("Boundary_SelectVarForBCI:\n");
  printf("  called with GH=%p, faces=%d, width=%d, table_handle=%d, "
         "var_index=%d, bc_name=%s\n",
         GH, faces, width, table_handle, var_index, bc_name);
  printf("  vi %d corresponds to %s\n", var_index, CCTK_VarName(var_index));
#endif

  /* Check to see if this is a new GH */
  if (!theGH) /* This is the first valid GH passed to a Boundary routine */
  {
    theGH = GH;
  } else if (GH != theGH) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "New GH passed to Boundary_SelectVarForBCI.  "
               "Thorn CactusBase/Boundary does not yet handle multiple GHs "
               "properly.");
    retval = -5;
  }

  /* Check that this request is allowed
     ---------------------------------- */

  /* Has some function implementing bc_name been registered? */
  if (!Util_TableQueryValueInfo(physbc_table_handle, NULL, NULL, bc_name)) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There is no function implementing the physical boundary "
               "condition %s",
               bc_name);
    retval = -2;
  }

/* Have any of these faces already been selected for a bc? */
#ifdef DEBUG
  print_selected_faces();
#endif
  if (selected_faces && selected_faces[var_index] & faces) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "%s has already been selected for a bc",
               CCTK_VarName(var_index));
    retval = -3;
  }

  /* Honor request
     ------------- */
  if (!retval) {

    /* allocate memory for new entry in database */
    new_entry = malloc(sizeof(struct BCVAR));
    if (!new_entry) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Unable to allocate memory for entry into "
                 "'selected for bcs' database");
      /*retval = -1;*/
    }

    /* populate new entry with data */
    new_entry->faces = faces;
    new_entry->width = width;
    new_entry->table = table_handle;
    new_entry->var = var_index;
    /* new_entry -> next will be filled in later */

    /* Into which of bcdata's lists should this variable be placed? */
    for (current_bcdata = bcdata_list; current_bcdata;
         previous_bcdata = current_bcdata,
        current_bcdata = current_bcdata->next) {
#ifdef DEBUG
      printf("Boundary_SelectVarForBCI: looping through bcdata list, at "
             "current_bcdata for %s\n",
             current_bcdata->bc_name);
#endif

      if (CCTK_Equals(current_bcdata->bc_name, bc_name)) {
        current = current_bcdata->var_list;
        current_bcdata->num++;
#ifdef DEBUG
        printf("Boundary_SelectVarForBCI: var %s brings bc %s to %d vars\n",
               CCTK_VarName(var_index), bc_name, current_bcdata->num);
#endif
        break; /* now that current is set we don't need to look at any
                  more bcdata entries */
      }
    }

    /* If current_bcdata is NULL, we got to the end of the above loop, this is
     * a new bc_name that does not appear in the bcdata list.
     */
    if (!current_bcdata) /* bc_name was not found in bcdata_list */
    {
#ifdef DEBUG
      printf("Boundary_SelectVarForBCI: adding new entry to bcdata list\n");
#endif

      /* new bc_name.  Create new entry for bcdata list. */
      new_bcdata = malloc(sizeof(struct BCDATA));
      if (!new_bcdata) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unable to allocate memory for internal 'bcdata' list");
        /*retval = -4;*/
      }
      /* new_bcdata->next is set below, when inserting into bcdata_list */
      new_bcdata->var_list = NULL;
      new_bcdata->bc_name = Util_Strdup(bc_name);
      /* new_bcdata->num is set further below, when adding first entry to
         var_list */

      /* Place new entry into bcdata list, maintaining case independent sort */
      for (current_bcdata = bcdata_list, previous_bcdata = NULL; current_bcdata;
           previous_bcdata = current_bcdata,
          current_bcdata = current_bcdata->next) {
#ifdef DEBUG
        printf("  looping through bcdata list, to insert new entry; at %s\n",
               current_bcdata->bc_name);
#endif
        if (Util_StrCmpi(bc_name, current_bcdata->bc_name) < 0) {
          /* bc_name precedes current->bc_name; place new entry here */
          if (!previous_bcdata) /* goes at start of bcdata list */
          {
#ifdef DEBUG
            printf("  new entry goes at beginning of bcdata list\n");
#endif
            bcdata_list = new_bcdata;
            new_bcdata->next = current_bcdata;
          } else {
            new_bcdata->next = current_bcdata;
            previous_bcdata->next = new_bcdata;
          }
          break;
        }
      }

      /* If current_bcdata still NULL, this is the last entry in the list */
      if (!current_bcdata) {
        if (!bcdata_list) /* list is empty */
        {
          bcdata_list = new_bcdata;
        } else {
          previous_bcdata->next = new_bcdata;
        }
        new_bcdata->next = NULL;
      }

      /* Set current_bcdata to new_bcdata, so the new bcdata entry will be
         filled in below */
      current_bcdata = new_bcdata;
    }

#ifdef DEBUG
    printf("  Finished sorting out which bcdata to use.  Now add entry to var "
           "list.\n");
    printf("  Preparing to loop through elements of var list.  current is %p"
           " previous is %p\n",
           current, previous);
#endif

    if (!current) /* This is the first element in the var_list */
    {
#ifdef DEBUG
      printf("  New element of var list.\n");
#endif
      current_bcdata->var_list = new_entry;
      new_entry->next = NULL;
      current_bcdata->num = 1;
    }

    /* Enter new_entry into correct location in linked list.
     * Note that this loop is skipped if new_entry was already inserted as
     * first element of a new var list above (since in that case current will
     * be NULL) */
    for (; /* starting value for current is selected using the bcdata
              list, above */
         current && entry_greater_than(new_entry, current) > 0;
         /* continue if not at end of list, and new_entry is greater than
            current entry */
         previous = current,
         current = current->next)
    /* store previous value for later use */
    {
    }

    /* The possibilities:
     * 1 nothing NULL: new_entry goes between previous and current
     * 2 previous NULL, but current non-null : new_entry goes at the start of
     *   the list
     * 3 current NULL, previous non-NULL : new_entry goes at the end of list
     * 4 both NULL: selections list is empty, case already caught above
     */
    if (previous) /* case 1 or 3 */
    {
      if (current) /* case 1 : goes in middle of list */
      {
        previous->next = new_entry;
        new_entry->next = current;
      } else /* case 3 : goes at end of list */
      {
        previous->next = new_entry;
        new_entry->next = NULL;
      }
    } else /* case 2 or 4 */
    {
      if (current) /* case 2 : goes at start of list */
      {
        current_bcdata->var_list = new_entry;
        new_entry->next = current;
      } /* case 4 : starts list, this case has already been handled above */
    }

    /* Record that this variable has been selected for a bc on these
       faces, for duplicate physical bc checking */
    if (!selected_faces) {
      num_cctk_vars = CCTK_NumVars();
      selected_faces = calloc(num_cctk_vars, sizeof(CCTK_INT));
      if (!selected_faces) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unable to allocate memory for internal 'selected_faces' "
                   "array");
      }
    }
    selected_faces[var_index] |= faces;
  }

#ifdef DEBUG
  print_selections_database();
  printf("\n");
#endif
  return retval;
}

/*@@
  @routine    Bdry_Boundary_SelectGroupForBC
  @date       Thu Dec 26 21:45:34 CET 2002
  @author     David Rideout
  @desc
              Used to select a Cactus variable group to have boundary
              conditions applied, using the group name.
              Table handle and faces must be the same for each member of the
              group.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        faces
  @vdesc      set of faces to which to apply bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        width
  @vdesc      if >=0, width of boundary in all directions
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        table_handle
  @vdesc      handle of table which holds arguments to be passed to bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        group_name
  @vdesc      name of group to which to apply bc
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @var        bc_name
  @vdesc      name of bc to apply
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
   0 success
  -6 invalid group name
   or the returncode of @seeroutine Bdry_Boundary_SelectGroupForBCI
  @endreturndesc
@@*/
CCTK_INT Bdry_Boundary_SelectGroupForBC(CCTK_POINTER_TO_CONST _GH,
                                        CCTK_INT faces, CCTK_INT width,
                                        CCTK_INT table_handle,
                                        CCTK_STRING group_name,
                                        CCTK_STRING bc_name) {
  int retval, gi;
  const cGH *GH = _GH;

#ifdef DEBUG
  printf("Boundary_SelectGroupForBC: called for group %s\n", group_name);
#endif

  /* get group index */
  gi = CCTK_GroupIndex(group_name);
  if (gi < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group name %s\n", group_name);
    retval = -6;
  } else {
    /* call Bdry_Boundary_SelectGroupForBCI() */
    retval = Bdry_Boundary_SelectGroupForBCI(GH, faces, width, table_handle, gi,
                                             bc_name);
  }

  return retval;
}

/*@@
  @routine    Bdry_Boundary_SelectGroupForBCI
  @date       Thu Dec 26 21:45:34 CET 2002
  @author     David Rideout
  @desc
              Used to select a Cactus variable group to have boundary
              conditions applied, using the group index.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        faces
  @vdesc      set of faces to which to apply bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        width
  @vdesc      if >=0, width of boundary in all directions
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        table_handle
  @vdesc      handle of table which holds arguments to be passed to bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        group_index
  @vdesc      index of group to which to apply bc
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        bc_name
  @vdesc      name of bc to apply
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
   0 success
  -7 invalid group index
  -8 group has zero vars
   or the returncode of @seeroutine Bdry_Boundary_SelectVarForBCI
  @endreturndesc
@@*/
CCTK_INT Bdry_Boundary_SelectGroupForBCI(CCTK_POINTER_TO_CONST _GH,
                                         CCTK_INT faces, CCTK_INT width,
                                         CCTK_INT table_handle,
                                         CCTK_INT group_index,
                                         CCTK_STRING bc_name) {
  int num_vars, vi, max_vi, retval;
  const cGH *GH = _GH;

  retval = -8;

  /* Get var indices from group name */
  num_vars = CCTK_NumVarsInGroupI(group_index);
  if (num_vars < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Invalid group index");
    retval = -7;
    return retval;
  }

#ifdef DEBUG
  printf("Boundary_SelectGroupForBCI: group %s has %d vars\n",
         CCTK_GroupName(group_index), num_vars);
#endif

  /* loop over variables in group */
  vi = CCTK_FirstVarIndexI(group_index);
  max_vi = vi + num_vars;
  for (; vi < max_vi; ++vi) {
    retval = Bdry_Boundary_SelectVarForBCI(GH, faces, width, table_handle, vi,
                                           bc_name);
  }

  return retval;
}

/*@@
  @routine    Bdry_Boundary_SelectedGVs
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Returns list of variable indices and table handles of
              variables selected for boundary conditions.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        array_size
  @vdesc      size of arrays to which var_indices and table_handles point
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        var_indices
  @vdesc      array into which selected variable indices will be placed
  @vtype      CCTK_INT
  @vio        out
  @endvar
  @var        faces
  @vdesc      array into which a set of selected faces for variables selected
              for bc will be placed
  @vtype      CCTK_INT
  @vio        out
  @endvar
  @var        widths
  @vdesc      array into which boundary widths of selected variables will be
              placed
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        table_handles
  @vdesc      array into which table_handles for variables selected for bc
              will be placed
  @vtype      CCTK_INT
  @vio        out
  @endvar
  @var        bc_name
  @vdesc      name of bc for which to get the selected vars,
              NULL returns all selected vars for all bcs
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
  -1 no boundary condition registered under bc_name
  -5 new value passed for GH
   number of variables selected for bc_name
  @endreturndesc
@@*/
CCTK_INT Bdry_Boundary_SelectedGVs(CCTK_POINTER_TO_CONST _GH,
                                   CCTK_INT array_size, CCTK_INT *var_indices,
                                   CCTK_INT *faces, CCTK_INT *widths,
                                   CCTK_INT *table_handles,
                                   CCTK_STRING bc_name) {
  int retval, i, j;
  struct BCVAR *current;
  struct BCDATA *current_bcdata;
  const cGH *GH = _GH;

  current = NULL;
  retval = 0;

  /* Check to see if this is a new GH */
  if (!theGH) /* This is the first valid GH passed to a Boundary routine */
  {
    theGH = GH;
  } else if (GH != theGH) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "New GH passed to Boundary_SelectedGVs.  "
               "Thorn CactusBase/Boundary does not yet handle multiple GHs "
               "properly.");
    retval = -5;
  }

#ifdef DEBUG
  printf("Boundary_SelectedGVs: called with bc_name=\"%s\" array_size=%d\n",
         bc_name, array_size);
  fflush(stdout);
#endif

  i = 0; /* i indexes the location in the returned arrays */

  /* Step through bcdata list */
  for (current_bcdata = bcdata_list; current_bcdata;
       current_bcdata = current_bcdata->next) {
#ifdef DEBUG
    printf("  looping through bcdata list, at bcdata entry for %s bc\n",
           current_bcdata->bc_name);
#endif

    if (!bc_name || CCTK_Equals(current_bcdata->bc_name, bc_name)) {

      /* Add these selected vars to return value */
      retval += current_bcdata->num;
      current = current_bcdata->var_list;

      /* Loop through var list */
      for (j = 0; /* j counts the bcs to check internal consistency */
           i < array_size && current; current = current->next, ++i, ++j) {

#ifdef DEBUG
        printf("    looping through selected vars, at current->var_index = "
               "%d\n",
               current->var);
        printf("      current->next is %p\n", current->next);
#endif

        if (faces) {
          faces[i] = current->faces;
        }
        if (widths) {
          widths[i] = current->width;
        }
        if (table_handles) {
          table_handles[i] = current->table;
        }
        if (var_indices) {
          var_indices[i] = current->var;
        }
      }
      if (j > current_bcdata->num)
        CCTK_WARN(0, "internal error");
      if (i != array_size && j != current_bcdata->num)
        CCTK_WARN(0, "internal error");
    }
  }

  /* Warn if there is no bc registered under this name */
  if (bc_name &&
      !Util_TableQueryValueInfo(physbc_table_handle, NULL, NULL, bc_name)) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There is no boundary condition registered under the name %s",
               bc_name);
    retval = -1;
  }

  return retval;
}

/********************************************************************
 *********************     Scheduled Routines   **********************
 ********************************************************************/

/*@@
  @routine    Boundary_ApplyPhysicalBCs
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              This will apply all requested physical boundary conditions.
  @enddesc
  @calls
  @history
  @endhistory
  @var        CCTK_ARGUMENTS
  @vdesc      Cactus argument list
  @vtype      CCTK_*
  @vio        in
  @endvar
  @returntype void
  @returndesc
  @endreturndesc
@@*/

void Boundary_ApplyPhysicalBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  phys_bc_fn_ptr bc_fn;
  int num_vars, max_num_vars, err;
  CCTK_INT *vars, *faces, *widths, *tables;
  struct BCDATA *current_bcdata;

  max_num_vars = 0;
  vars = NULL;   /* avoids a compiler warning */
  faces = NULL;  /* avoids a compiler warning */
  widths = NULL; /* avoids a compiler warning */
  tables = NULL; /* avoids a compiler warning */

  /* Warning: This function does not consider which GH it is called on */

  /* Step through each requested physical boundary condition */
  for (current_bcdata = bcdata_list; current_bcdata;
       current_bcdata = current_bcdata->next) {
#ifdef DEBUG
    printf("Boundary_ApplyPhysicalBCs: looping through bcdata list, at bcdata "
           "entry for %s bc\n",
           current_bcdata->bc_name);
    printf("  max_num_vars is %d\n", max_num_vars);
#endif

    /* Allocate memory to hold selected vars */
    num_vars = current_bcdata->num;
#ifdef DEBUG
    printf("  num_vars is %d\n", num_vars);
#endif
    if (num_vars > max_num_vars) {
      max_num_vars = num_vars; /* store new maximum */
      vars = realloc(vars, num_vars * sizeof(CCTK_INT));
      faces = realloc(faces, num_vars * sizeof(CCTK_INT));
      widths = realloc(widths, num_vars * sizeof(CCTK_INT));
      tables = realloc(tables, num_vars * sizeof(CCTK_INT));
    }

    /* get selected vars for this bc_name */
    err = Bdry_Boundary_SelectedGVs(cctkGH, num_vars, vars, faces, widths,
                                    tables, current_bcdata->bc_name);
    if (err < 0) /* This is a redundant test for now,
                    Bdry_Boundary_SelectedGVs never returns <0 */
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error in Boundary_SelectedGVs for %s boundary condition",
                 current_bcdata->bc_name);
    } else if (err != num_vars) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Boundary_SelectedGVs returned %d selected variables for "
                 "\"%s\" boundary condition, but %d expected\n",
                 err, current_bcdata->bc_name, num_vars);
    }

/* Get the fn ptr for the bc */
#ifdef DEBUG
    printf("Boundary_ApplyPhysicalBCs: current_bcdata->bc_name=\"%s\"\n",
           current_bcdata->bc_name);
#endif
    err = Util_TableGetFnPointer(physbc_table_handle, (CCTK_FPOINTER *)&bc_fn,
                                 current_bcdata->bc_name);
    if (err < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Boundary_ApplyPhysicalBCs: Util_TableGetFnPointer "
                 "returned %d",
                 err);
    }

    /* check that all variables have storage at least on the current level
     * individual boundary conditions (eg static) may check for more */
    err = 0;
    for (int i = 0; i < num_vars; i++) {
      if (CCTK_VarDataPtrI(cctkGH, 0, vars[i]) == NULL) {
        char *fullname = CCTK_FullName(vars[i]);
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Boundary_ApplyPhysicalBCs: variable \"%s\" has no storage "
                   "when attempting to apply boundary condition \"%s\".",
                   fullname, current_bcdata->bc_name);
        err += 1;
        free(fullname);
      }
    }
    if (err) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Boundary_ApplyPhysicalBCs: boundary conditions were "
                  "requested for %d variables that do not have storage",
                  err);
    }

/* Apply bc to vi */
#ifdef DEBUG
    printf(
        "Boundary_ApplyPhysicalBCs: Attempting to call boundary condition\n"
        "                           Using function pointer %p with arguments\n"
        "                           cctkGH %p, num_vars %d, vars, tables\n",
        (void *)bc_fn, (const void *)cctkGH, num_vars);
#endif
    err = (*bc_fn)(cctkGH, num_vars, vars, faces, widths, tables);
    if (err < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Function associated with boundary condition %s returned %d",
                 current_bcdata->bc_name, err);
    }
  }

  /* Free data */
  free(vars);
  free(faces);
  free(widths);
  free(tables);
}

/*@@
  @routine    Boundary_ClearSelection
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Clears all boundary condition selections.
  @enddesc
  @calls
  @history
  @endhistory
@@*/

void Boundary_ClearSelection(CCTK_ARGUMENTS) {
  struct BCVAR *current, *next;
  struct BCDATA *current_bcdata;

/* Warning: This function does not consider which GH it is called on, which
   will be a bug in Cactus 4.1. */

/* Step through bcdata list */
#ifdef DEBUG
  printf("Boundary_ClearSelection: looping through bcdata list for freeing\n");
#endif
  for (current_bcdata = bcdata_list; current_bcdata;
       current_bcdata = current_bcdata->next) {
#ifdef DEBUG
    printf("Boundary_ClearSelection: freeing var list rooted at %p\n",
           current_bcdata);
#endif

    /* Free selections list */
    current = next = current_bcdata->var_list;
    for (; current; current = next) {
      next = current->next;
      free(current);
    }

    current_bcdata->var_list = NULL;
    current_bcdata->num = 0;
  }

  /* Clear list of selected faces */
  memset(selected_faces, 0, num_cctk_vars * sizeof(CCTK_INT));
}

/*@@
  @routine    Boundary_MakeSureThatTheSelectionIsEmpty
  @date       Mon Jul  7 21:51:37 CET 2003
  @author     Erik Schnetter
  @desc
              Abort if the selections is not empty.
              This routine is currently unused, but is very
              helpful for debugging.
  @enddesc
  @calls
  @history
  @endhistory
@@*/

void Boundary_MakeSureThatTheSelectionIsEmpty(void);
void Boundary_MakeSureThatTheSelectionIsEmpty(void) {
  struct BCDATA *current_bcdata;
  struct BCVAR *current;
  int is_empty;

  is_empty = 1;
  for (current_bcdata = bcdata_list; current_bcdata;
       current_bcdata = current_bcdata->next) {
    for (current = current_bcdata->var_list; current; current = current->next) {
      char *fullname = CCTK_FullName(current->var);
      if (is_empty) {
        is_empty = 0;
        printf("The following boundary conditions are currently selected for "
               "the following variables:\n");
      }
      printf("Boundary condition %s for variable %s\n", current_bcdata->bc_name,
             fullname);
      free(fullname);
    }
  }

  if (!is_empty) {
    CCTK_WARN(0, "Someone thinks that the boundary selection list should be "
                 "empty at this point.  Alas, it is not; I better abort.");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/*@@
  @routine    entry_greater_than
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Sorts entries in selections list.
              Returns non-zero value if new > current.
  @enddesc
  @calls
  @history
  @endhistory
  @var        new
  @vdesc      new entry to be inserted into selections list
  @vtype      struct BCVAR *
  @vio        in
  @endvar
  @var        current
  @vdesc      entry in selections list to which to compare new
  @vtype      struct BCVAR *
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
             -1 new equals current
              0 new 'less than' current
              1 new->table_handle > current->table_handle
              2 new and current have same table_handle, but new has greater
                var_index
  @endreturndesc
@@*/
static int entry_greater_than(struct BCVAR *new, struct BCVAR *current) {
  int retval;

  /* can assume both arguments are valid (non-null) */
  if (new->table > current->table) {
    retval = 1;
  } else if (new->table < current->table) {
    retval = 0;
  } else if (new->var > current->var) {
    retval = 2;
  } else if (new->var < current->var) {
    retval = 0;
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "%s has already been selected for this boundary condition!",
               CCTK_VarName(new->var));
    retval = -1;
  }

  return retval;
}

/*@@
  @routine    print_selections_database
  @date       13 March 2003
  @author     David Rideout
  @desc
              Prints selected variables database, for debugging only
  @enddesc
  @calls
  @history
  @endhistory
@@*/

#ifdef DEBUG

static void print_selections_database(void) {
  struct BCDATA *current_bcdata;
  struct BCVAR *current;

  printf("Current list of selected vars:\n");
  for (current_bcdata = bcdata_list; current_bcdata;
       current_bcdata = current_bcdata->next) {
    printf("%d entries for %s:\n", current_bcdata->num,
           current_bcdata->bc_name);
    printf(" vi  gi      var name  table handle\n");
    for (current = current_bcdata->var_list; current; current = current->next) {
      printf("%3d  %2d  %12s  %2d\n", current->var,
             CCTK_GroupIndexFromVarI(current->var), CCTK_VarName(current->var),
             current->table);
    }
  }
}

#endif

/*@@
  @routine    print_selected_faces
  @date       13 May 2003
  @author     David Rideout
  @desc
              Prints selected faces database, for debugging only
  @enddesc
  @calls
  @history
  @endhistory
@@*/

#ifdef DEBUG

static void print_selected_faces(void) {
  int i;

  if (selected_faces) {
    for (i = 0; i < num_cctk_vars; ++i) {
      printf("selected_faces[%d] = %d\n", i, selected_faces[i]);
    }
  } else {
    printf("no selected faces database yet\n");
  }
}

#endif
