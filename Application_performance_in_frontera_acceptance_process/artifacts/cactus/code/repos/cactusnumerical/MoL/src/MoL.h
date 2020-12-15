 /*@@
   @file      MoL.h
   @date      Thu Jun 13 22:17:32 2002
   @author    Ian Hawke
   @desc 
   A header file with the definitions of the functions to be 
   called by other thorns. This will exist until the function
   aliasing is sorted out.
   @enddesc 
 @@*/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MOL_ACTIVE
#define MOL_ACTIVE
#endif

#ifdef CCODE

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

/* Old functions from MoL1 for compatibility */

CCTK_INT MoL_RegisterVar(CCTK_INT molvarindex,CCTK_INT molrhsvarindex);
CCTK_INT MoL_RegisterPrimitive(CCTK_INT primitiveindex);
CCTK_INT MoL_RegisterDepends(CCTK_INT dependsindex);
CCTK_INT MoL_RegisterVarGroup(CCTK_INT groupindex,CCTK_INT rhsgroupindex);
CCTK_INT MoL_RegisterPrimitiveGroup(CCTK_INT groupindex);
CCTK_INT MoL_RegisterDependsGroup(CCTK_INT groupindex);
CCTK_INT MoL_ChangeVarToEvolved(CCTK_INT varindex, CCTK_INT rhsindex);
CCTK_INT MoL_ChangeVarToDependent(CCTK_INT dependsindex);
CCTK_INT MoL_ChangeVarToPrimitive(CCTK_INT primitiveindex);
CCTK_INT MoL_ChangeVarToNone(CCTK_INT removeindex);

#endif

#ifdef __cplusplus
}
#endif
