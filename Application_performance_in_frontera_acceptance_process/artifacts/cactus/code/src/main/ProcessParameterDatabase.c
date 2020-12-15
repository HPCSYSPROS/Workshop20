 /*@@
   @file      ProcessParameterDatabase.c
   @date      Thu Sep 24 10:34:46 1998
   @author    Tom Goodale
   @desc 
              Routines to determine the parameters and store them.
   @enddesc 
   @version   $Id$
 @@*/

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "cctk_Flesh.h"
#include "cctk_Parameter.h"
#include "cctk_WarnLevel.h"

#include "cctki_Parameter.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_ProcessParameterDatabase_c);

/********************************************************************
 *********************  Macro Definitions  **************************
 ********************************************************************/

/* some systems (eg. Windows NT) don't define this macro */
#ifndef S_ISDIR
#define S_ISDIR(mode)   (((mode) & S_IFMT) == S_IFDIR)
#endif

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/
int ParseFile (FILE *ifp, 
               int (*set_function) (const char *, const char *, int),
               tFleshConfig *ConfigData);
void CCTKi_SetParameterSetMask (int mask);


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_ProcessParameterDatabase
   @date       Thu Sep 24 10:37:07 1998
   @author     Tom Goodale
   @desc 
   
   @enddesc 
   @calls      CCTKi_SetParameterSetMask
               CCTKi_NumParameterFileErrors
               ParseFile
 
   @var        ConfigData
   @vdesc      Flesh configuration data
   @vtype      tFleshConfig *
   @vio        inout
   @endvar 

   @returntype int
   @returndesc
                0 - success, or<BR>
               -1 - unable to open parameter file
   @endreturndesc
@@*/
int CCTKi_ProcessParameterDatabase (tFleshConfig *ConfigData)
{
  int parse_errors;
  int major, minor;
  FILE *parameter_file;
  struct stat statbuf;


  CCTKi_SetParameterSetMask (PARAMETER_RECOVERY_PRE);

  if (! strcmp (ConfigData->parameter_file_name, "-"))
  {
    parameter_file = stdin;
  }
  else if (!stat(ConfigData->parameter_file_name, &statbuf))
  {
    if(S_ISDIR(statbuf.st_mode))
    {
      parameter_file = NULL;
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "Cannot open parameter file '%s': it is a directory", 
                  ConfigData->parameter_file_name);
    }
    else
    {
      parameter_file = fopen (ConfigData->parameter_file_name, "r");
    }
  }
  else
  {
    /* Stat failed */
    parameter_file = NULL;

    switch(errno)
    {
      case ENOENT:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s': file doesn't exist", 
                    ConfigData->parameter_file_name);
        break;
      case ENOTDIR:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s': the path is invalid", 
                    ConfigData->parameter_file_name);
        break;
#ifdef ELOOP
      case ELOOP:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s': too many symbolic links", 
                    ConfigData->parameter_file_name);
        break;
#endif /* ELOOP */
      case EACCES:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s': permission denied",
                    ConfigData->parameter_file_name);
        break;
      case ENAMETOOLONG:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s': filename too long",
                    ConfigData->parameter_file_name);
        break;
      case ENOMEM:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s': out of system memory",
                    ConfigData->parameter_file_name);
        break;
      default:
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "Cannot open parameter file '%s'",
                    ConfigData->parameter_file_name);
    }
  }
  
  if (parameter_file == NULL)
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "Cannot open parameter file '%s'.  "
                "(This can also be an MPI problem; "
                "e.g. you may be using the wrong version of 'mpirun', "
                "or may have forgotten to call 'lamboot'.)",
                ConfigData->parameter_file_name);
  }

  if (parameter_file)
  {
    parse_errors = ParseFile (parameter_file, CCTKi_SetParameter, ConfigData);

    if (strcmp (ConfigData->parameter_file_name, "-"))
    {
      fclose (parameter_file);
    }

    if (parse_errors)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "CCTKi_SetParameterSetMask: %d parsing errors in "
                  "parameter file", parse_errors);
    }      

    minor = CCTKi_NumParameterFileErrors(1);
    if (minor)
    {
      if (minor > 1) 
      {
        CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                    "CCTKi_SetParameterSetMask: %d minor errors in "
                    "parameter file", CCTKi_NumParameterFileErrors (1));
      }
      else
      {
        CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                    "CCTKi_SetParameterSetMask: %d minor error in "
                    "parameter file", CCTKi_NumParameterFileErrors (1));
      }

    }

    major = CCTKi_NumParameterFileErrors (0);
    if (major)
    {
      if (major > 1) 
      {
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "CCTKi_SetParameterSetMask: %d major errors in "
                    "parameter file", CCTKi_NumParameterFileErrors (0));
      }
      else
      {
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "CCTKi_SetParameterSetMask: %d major error in "
                    "parameter file", CCTKi_NumParameterFileErrors (0));
      }

    }
  }
  else
  {
    CCTK_VWarn(0, __LINE__, __FILE__, "Cactus", 
               "Unable to open parameter file '%s'\n", 
               ConfigData->parameter_file_name);
  }
      
  return (parameter_file ? 0 : -1);
}
