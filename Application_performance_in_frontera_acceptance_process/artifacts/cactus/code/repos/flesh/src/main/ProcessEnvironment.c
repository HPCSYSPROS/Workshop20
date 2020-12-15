 /*@@
   @file      ProcessEnvironment.c
   @date      Fri Feb 26 11:20:15 1999
   @author    Tom Goodale
   @desc 
   Checks the environment for various settings, and acts on them.
   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk_Capabilities.h"
#include "cctk_Flesh.h"

#ifdef HAVE_CAPABILITY_MPI
#  include <mpi.h>
#endif

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_ProcessEnvironment_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/
int CCTKi_ProcessEnvironment (int *argc, char ***argv,tFleshConfig *ConfigData);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

#ifdef HAVE_CAPABILITY_MPI
#define CACTUS_MPI_ERROR(xf)                                                  \
          do                                                                  \
          {                                                                   \
            int errcode;                                                      \
                                                                              \
                                                                              \
            if((errcode = xf) != MPI_SUCCESS)                                 \
            {                                                                 \
              char mpi_error_string[MPI_MAX_ERROR_STRING+1];                  \
              int resultlen;                                                  \
                                                                              \
                                                                              \
              MPI_Error_string(errcode, mpi_error_string, &resultlen);        \
              fprintf(stderr, "MPI Call %s returned error code %d (%s)\n",    \
                              #xf, errcode, mpi_error_string);                \
              fprintf(stderr, "At line %d of file %s\n",                      \
                              __LINE__, __FILE__);                            \
            }                                                                 \
          } while (0)
#endif

#ifdef HAVE_CAPABILITY_MPI
char cctki_MPI_Active = 0;      /* We're using MPI */
char cctki_MPI_Managing = 0;    /* We manage MPI */
#endif


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_ProcessEnvironment
   @date       Fri Feb 26 11:20:15 1999
   @author     Tom Goodale
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     argc
   @vdesc   Number of arguments
   @vtype   int *
   @vio     inout
   @vcomment 
 
   @endvar 
   @var     argv
   @vdesc   Argument list
   @vtype   char ***
   @vio     inout
   @vcomment 
 
   @endvar 
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
int CCTKi_ProcessEnvironment(int *argcp, char ***argvp,
                             tFleshConfig *ConfigData)
{
  /* avoid compiler warnings about unused arguments */
  argcp = argcp;
  argvp = argvp;
  ConfigData = ConfigData;

  /* Check if MPI compiled in but choosing not to use MPI. */  
#ifdef HAVE_CAPABILITY_MPI
  if (!getenv("CACTUS_NOMPI"))
  {
    cctki_MPI_Active = 1;
    int initialized;
    CACTUS_MPI_ERROR(MPI_Initialized(&initialized));
    if (!initialized)
    {
      cctki_MPI_Managing = 1;
      /* CACTUS_MPI_ERROR(MPI_Init(argcp, argvp)); */
      int provided;
      CACTUS_MPI_ERROR(MPI_Init_thread(argcp, argvp,
                                       MPI_THREAD_SERIALIZED, &provided));
    }
  }
#endif

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
