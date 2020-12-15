/*@@
  @file      Register.c
  @date      Sat Oct 26 22:39:40 CEST 2002
  @author    David Rideout
  @desc
             Register implemented boundary conditions.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "Boundary.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_Boundary_Register_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

void Boundary_RegisterBCs(CCTK_ARGUMENTS);

/********************************************************************
 ***************** Aliased Routine Prototypes ***********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     Aliased Routines   ***********************
 ********************************************************************/

/********************************************************************
 *********************     Scheduled Routines   *********************
 ********************************************************************/

/*@@
  @routine    Boundary_RegisterBCs
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
              Register all boundary conditions implemented by this thorn.
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
@@*/

void Boundary_RegisterBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  int err;

  if (register_scalar) {
    err = Boundary_RegisterPhysicalBC((CCTK_POINTER)cctkGH,
                                      (phys_bc_fn_ptr)&BndScalar, "Scalar");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Scalar\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_flat) {
    err = Boundary_RegisterPhysicalBC((CCTK_POINTER)cctkGH,
                                      (phys_bc_fn_ptr)&BndFlat, "Flat");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Flat\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_radiation) {
    err = Boundary_RegisterPhysicalBC(
        (CCTK_POINTER)cctkGH, (phys_bc_fn_ptr)&BndRadiative, "Radiation");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Radiation\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_copy) {
    err = Boundary_RegisterPhysicalBC((CCTK_POINTER)cctkGH,
                                      (phys_bc_fn_ptr)&BndCopy, "Copy");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Copy\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_robin) {
    err = Boundary_RegisterPhysicalBC((CCTK_POINTER)cctkGH,
                                      (phys_bc_fn_ptr)&BndRobin, "Robin");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Robin\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_static) {
    err = Boundary_RegisterPhysicalBC((CCTK_POINTER)cctkGH,
                                      (phys_bc_fn_ptr)&BndStatic, "Static");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Static\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_none) {
    err = Boundary_RegisterPhysicalBC((CCTK_POINTER)cctkGH,
                                      (phys_bc_fn_ptr)&BndNone, "None");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"None\" "
                 "boundary condition",
                 err);
    }
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
