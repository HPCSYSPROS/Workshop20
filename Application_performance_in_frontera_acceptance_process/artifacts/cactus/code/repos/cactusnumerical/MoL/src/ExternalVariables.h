 /*@@
   @file      ExternalVariables.h
   @date      Wed May 22 02:32:10 2002
   @author    Ian Hawke
   @desc 
   The header file containing the local variables used across routines.
   These are the arrays containing GF indexes for all types of variables,
   and the number of each type of variable currently in use (the 
   parameters only give the maximum possible number).
   No function prototypes are defined in this file, so we do not protect
   it with an ifdef so that we can do inclusion within multiple routines
   in the same file.
   @enddesc 
   @version   $Header$
 @@*/


#include <cctk.h>


extern CCTK_INT *EvolvedVariableIndex;
extern CCTK_INT *EvolvedVariableIndexSlow;
extern CCTK_INT *RHSVariableIndex;
extern CCTK_INT *RHSVariableIndexSlow;
extern CCTK_INT *ConstrainedVariableIndex;
extern CCTK_INT *SandRVariableIndex;


extern CCTK_INT MoLNumEvolvedVariables;
extern CCTK_INT MoLNumEvolvedVariablesSlow;
extern CCTK_INT MoLNumConstrainedVariables;
extern CCTK_INT MoLNumSandRVariables;



extern CCTK_INT *EvolvedArrayVariableIndex;
extern CCTK_INT *RHSArrayVariableIndex;
extern CCTK_INT *ConstrainedArrayVariableIndex;
extern CCTK_INT *SandRArrayVariableIndex;


extern CCTK_INT MoLNumEvolvedArrayVariables;
extern CCTK_INT MoLNumConstrainedArrayVariables;
extern CCTK_INT MoLNumSandRArrayVariables;



extern CCTK_INT ScheduleStatus;


extern CCTK_REAL *ArrayScratchSpace;
extern CCTK_INT *ArrayScratchSizes;
extern CCTK_INT CurrentArrayScratchSize;

extern CCTK_INT MoLMaxNumRegisteredVariables;
