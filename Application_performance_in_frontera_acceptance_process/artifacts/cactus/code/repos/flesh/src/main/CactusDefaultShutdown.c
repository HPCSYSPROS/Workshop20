 /*@@
   @file      CactusDefaultShutdown.c
   @date      Fri Feb 26 16:53:58 1999
   @author    Tom Goodale
   @desc 
   The default shutdown routines.
   @enddesc 
 @@*/

#include <stdio.h>

#include "cctk_Flesh.h"
#include "cctk_Comm.h"
#include "CactusMainDefaults.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_CactusDefaultShutdown_c);

 /*@@
   @routine    CactusDefaultShutdown
   @date       Tue Sep 29 12:45:04 1998
   @author     Tom Goodale 
   @desc 
   Default shutdown routine.
   @enddesc 
   @calls     
   @calledby   
   @history introducing CCTK_SHUTDOWN scheduling [03/00  Gerd Lanfermann]
 
   @endhistory 

@@*/
int CactusDefaultShutdown(tFleshConfig *config)
{
  unsigned int conv_level;

  /* Execute termination for all convergence levels */
  for(conv_level = 0 ; conv_level < config->nGHs;  conv_level++) 
  {    
    CCTK_Traverse(config->GH[conv_level], "CCTK_TERMINATE"); 
  }
 
  /* Execute shutdown for all convergence levels */
  for(conv_level = 0 ; conv_level < config->nGHs;  conv_level++) 
  {    
    CCTK_Traverse(config->GH[conv_level], "CCTK_SHUTDOWN"); 
  }

  return 0;
}
