 /*@@
   @file      MoLFunctions.h
   @date      Thu Jun 13 22:17:32 2002
   @author    Ian Hawke
   @desc 
   A header file with the definitions of the functions to be 
   called by other thorns. This will exist until the function
   aliasing is sorted out.
   @enddesc 
 @@*/

#ifndef MOL_FUNCTIONS_H
#define MOL_FUNCTIONS_H

CCTK_INT MoL_RegisterEvolved(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);
CCTK_INT MoL_RegisterEvolvedSlow(CCTK_INT EvolvedIndex, CCTK_INT RHSIndexSlow);
CCTK_INT MoL_RegisterConstrained(CCTK_INT ConstrainedIndex);
CCTK_INT MoL_RegisterSaveAndRestore(CCTK_INT SandRIndex);
CCTK_INT MoL_ChangeToEvolved(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);
CCTK_INT MoL_ChangeToConstrained(CCTK_INT ConstrainedIndex);
CCTK_INT MoL_ChangeToSaveAndRestore(CCTK_INT SandRIndex);
CCTK_INT MoL_ChangeToNone(CCTK_INT RemoveIndex);
CCTK_INT MoL_RegisterEvolvedGroup(CCTK_INT EvolvedGroupIndex, 
                                  CCTK_INT RHSGroupIndex);
CCTK_INT MoL_RegisterEvolvedGroupSlow(CCTK_INT EvolvedGroupIndexSlow,
                                      CCTK_INT RHSGroupIndexSlow);
CCTK_INT MoL_RegisterConstrainedGroup(CCTK_INT ConstrainedGroupIndex);
CCTK_INT MoL_RegisterSaveAndRestoreGroup(CCTK_INT SandRGroupIndex);

#endif
