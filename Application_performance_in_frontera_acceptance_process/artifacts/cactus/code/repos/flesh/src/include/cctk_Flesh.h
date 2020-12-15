 /*@@
   @header    cctk_Flesh.h
   @date      Thu Sep 24 10:18:52 1998
   @author    Tom Goodale
   @desc 
              Header file for flesh functions.
   @enddesc 
   @version   $Header$
 @@*/

#ifndef _CCTK_FLESH_H_
#define _CCTK_FLESH_H_

#include "cGH.h"


/*  Typedefs */

typedef struct
{
  char *parameter_file_name;

  /* Array of pointers to cactus grid hierarchies. */
  cGH **GH;
  unsigned int nGHs;

  /* flag telling whether we restart from a checkpoint or not */
  int recovered;

  /*  cTimer *timer[3];*/
} tFleshConfig;


/* Function prototypes */

#ifdef __cplusplus
extern "C" 
{
#endif

#define CCTK_FILEVERSION(file) const char *CCTKi_version_##file (void);       \
                               const char *CCTKi_version_##file (void)        \
                               { return (rcsid); }

int CCTK_Traverse (cGH *GH, const char *where);
int CCTKi_ProcessCommandLine (int *argc, char ***argv, tFleshConfig *config);
int CCTKi_ProcessEnvironment (int *argc, char ***argv, tFleshConfig *config);
int CCTKi_InitialiseDataStructures (tFleshConfig *config);
int CCTKi_ProcessParameterDatabase (tFleshConfig *config);
int CCTKi_CallStartupFunctions (tFleshConfig *config);
int CCTKi_AddGH (tFleshConfig *config, unsigned int convergence_level, cGH *GH);
int CCTKi_InitialiseCactus (int *argc, char ***argv, tFleshConfig *config);
int CCTKi_ShutdownCactus (tFleshConfig *config);
int CCTKi_DummyExit (cGH *GH, int retval);

#ifdef __cplusplus
}
#endif

#endif  /* _CCTK_FLESH_H_ */
