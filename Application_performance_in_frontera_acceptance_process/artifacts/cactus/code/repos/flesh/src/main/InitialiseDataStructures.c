 /*@@
   @file      InitialiseDataStructures.c
   @date      Wed Jan 13 20:28:08 1999
   @author    Tom Goodale
   @desc 
   Initialise various datastructures.
   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk_Flesh.h"
#include "cctk_ActiveThorns.h"
#include "cctki_ActiveThorns.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_InitialiseDataStructures_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

int CCTKi_RegisterDefaultTimerFunctions(void);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_InitialiseDataStructures
   @date       Wed Jan 20 09:27:56 1999
   @author     Tom Goodale
   @desc 
   
   @enddesc 
   @calls     CCTKi_ActivateThorn
   @calledby   
   @history 
 
   @endhistory 
   @var     ConfigData
   @vdesc   Flesh configuration data
   @vtype   tFleshConfig
   @vio     inout
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc
   0  - success
   @endreturndesc
@@*/

int CCTKi_InitialiseDataStructures(tFleshConfig *ConfigData)
{

  CCTKi_RegisterDefaultTimerFunctions();

  ConfigData->nGHs = 0;
  ConfigData->GH = NULL;

#if 0
  ConfigData->timer[INITIALISATION] = CactusNewTimer();
  ConfigData->timer[EVOLUTION] = CactusNewTimer();
  ConfigData->timer[ELLIPTIC] = CactusNewTimer();
#endif

  CCTKi_ActivateThorn("Cactus");

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
