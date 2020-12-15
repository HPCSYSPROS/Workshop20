#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Termination.h"
#include "cctk_Timers.h"
#include "util_String.h"


void TriggerTerminationFile (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char filename[10000];
  FILE *file;
  int terminate;

  /* only one processor needs to check the file */
  if (CCTK_MyProc (cctkGH) != 0)
  {
    return;
  }

  if (cctk_iteration % check_file_every != 0)
  {
    return;
  }

  if (termination_file[0] == '/')
  {
    /* file name begins with a slash, do not use IO::out_dir */
    Util_snprintf (filename, sizeof filename, "%s", termination_file);
  }
  else
  {
    /* add IO::out_dir to filename */
    Util_snprintf (filename, sizeof filename,
                   "%s/%s", out_dir, termination_file);
  }

  file = fopen (filename, "r");

  if (file != NULL) {
    terminate = 0;
    fscanf (file, "%d", &terminate);
    fclose (file);
    
    if (terminate == 1) {
      CCTK_INFO ("OH MY GOD! Found termination signal in termination file. Termination NOW!!! ");

      CCTK_TerminateNext (cctkGH);
    }
   
  }
}
