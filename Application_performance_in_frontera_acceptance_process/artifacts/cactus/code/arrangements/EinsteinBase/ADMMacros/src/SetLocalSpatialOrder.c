 /*@@
   @file      SetLocalSpatialOrder.c
   @date      Fri Jul  4 12:41:00 2003
   @author    Ian Hawke
   @desc 
   Sets the local_spatial_order to the spatial_order parameter.
   @enddesc 
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMMacros_SetLocalSpatialOrder_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void ADMMacros_SetLocalSpatialOrder(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    ADMMacros_SetLocalSpatialOrder
   @date       Fri Jul  4 12:42:15 2003
   @author     Ian Hawke
   @desc 
   Sets the local_spatial_order to the spatial_order parameter.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

 @@*/

void ADMMacros_SetLocalSpatialOrder(CCTK_ARGUMENTS)
{
 
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *local_spatial_order = spatial_order;

  if (CCTK_MyProc(cctkGH)==0)
    CCTK_VInfo(CCTK_THORNSTRING, "Spatial finite differencing order: %d",
               (int)spatial_order);

  return;
}
