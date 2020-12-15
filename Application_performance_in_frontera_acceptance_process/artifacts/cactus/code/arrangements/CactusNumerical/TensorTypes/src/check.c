/* $Header$ */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "tensortypes_imp.h"



void
CheckTensorTypes (CCTK_ARGUMENTS)
{
  int groupindex;
  int numvars;
  int firstvar;
  
  struct tensor const * tensortype;
  
  int allisgood;
  
  int n;
  
  /* Check all internal tensor type definitions */
  for (n=0; n<TT_numtensors; ++n) {
    CheckTensorType (TT_alltensors[n]);
  }
  
  /* Check tensor types of all groups */
  allisgood = 1;
  for (groupindex=0; groupindex<CCTK_NumGroups(); ++groupindex) {
    
    numvars = CCTK_NumVarsInGroupI (groupindex);
    if (numvars != 0) {
      assert (numvars>0);
      firstvar = CCTK_FirstVarIndexI (groupindex);
      assert (firstvar>=0);
      
      tensortype = TT_varindex2tensortype (firstvar);
      allisgood = allisgood && tensortype;
    }
  }
  
  if (! allisgood) {
    CCTK_PARAMWARN ("Some grid variable groups had incorrect tensor type declarations");
  }
}
