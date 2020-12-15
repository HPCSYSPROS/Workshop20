/* $Header$ */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "tensortypes_imp.h"



struct tensor const *
TT_varindex2tensortype (int const varindex)
{
  int groupindex;
  cGroup group;
  char tensortypealias[100];
  int ierr;
  
  assert (varindex>=0 && varindex<CCTK_NumVars());
  groupindex = CCTK_GroupIndexFromVarI (varindex);
  assert (groupindex>=0);
  ierr = CCTK_GroupData (groupindex, &group);
  assert (! ierr);
  
  ierr = Util_TableGetString
    (group.tagstable,
     sizeof tensortypealias, tensortypealias, "tensortypealias");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    char * const groupname = CCTK_GroupName (groupindex);
    assert (groupname);
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Tensor type alias not declared for group \"%s\" -- assuming a scalar",
                groupname);
    free (groupname);
    strcpy (tensortypealias, "");
  } else if (ierr<0) {
    char * const groupname = CCTK_GroupName (groupindex);
    assert (groupname);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Error in tensor type alias declaration for group \"%s\"",
                groupname);
    free (groupname);
    return NULL;
  }
  
  if (CCTK_EQUALS (tensortypealias, "")
      || CCTK_EQUALS (tensortypealias, "scalar"))
  {
    /* scalar */
    if (group.numvars != 1) {
      char * const groupname = CCTK_GroupName (groupindex);
      assert (groupname);
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                  groupname);
      free (groupname);
    }
    return &TT_scalar;
  } else if (CCTK_EQUALS (tensortypealias, "u")) {
    /* vector */
    if (group.numvars != 3) {
      char * const groupname = CCTK_GroupName (groupindex);
      assert (groupname);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the tensor type alias \"u\", but contains %d instead of 3 elements",
                  groupname, group.numvars);
      free (groupname);
      return NULL;
    }
    return &TT_vector;
  } else if (CCTK_EQUALS (tensortypealias, "dd_sym")) {
    /* symmetric tensor */
    if (group.numvars != 6) {
      char * const groupname = CCTK_GroupName (groupindex);
      assert (groupname);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the tensor type alias \"dd_sym\", but contains %d instead of 6 elements",
                  groupname, group.numvars);
      free (groupname);
      return NULL;
    }
    return &TT_symmtensor;
  }
  
  char * const groupname = CCTK_GroupName (groupindex);
  assert (groupname);
  CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
              "Illegal tensor type alias for group \"%s\"",
              groupname);
  free (groupname);
  return NULL;
}



int
TT_varindex2var (int const varindex)
{
  struct tensor const * tensortype;
  int groupindex;
  int numvars;
  int firstvarindex;
  int var;
  
  assert (varindex>=0 && varindex<CCTK_NumVars());
  
  tensortype = TT_varindex2tensortype (varindex);
  if (! tensortype) {
    char * const fullname = CCTK_FullName (varindex);
    assert (fullname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot determine tensor type for variable \"%s\"", fullname);
    free (fullname);
  }
  
  groupindex = CCTK_GroupIndexFromVarI (varindex);
  assert (groupindex>=0);
  
  numvars = CCTK_NumVarsInGroupI (groupindex);
  assert (numvars>=0);
  assert (numvars > 0);
  
  firstvarindex = CCTK_FirstVarIndexI (groupindex);
  assert (firstvarindex>=0);
  assert (firstvarindex <= varindex);
  
  var = varindex - firstvarindex;
  
  /* scalars are a special case */
  if (tensortype == &TT_scalar) {
    var = 0;
  }
  
  assert (var>=0 && var<tensortype->ncomps);
  
  return var;
}



int
TT_var2varindex (int const old_varindex, int const var)
{
  struct tensor const * tensortype;
  int groupindex;
  int numvars;
  int firstvarindex;
  int varindex;
  
  assert (old_varindex>=0 && old_varindex<CCTK_NumVars());
  
  tensortype = TT_varindex2tensortype (old_varindex);
  if (! tensortype) {
    char * const fullname = CCTK_FullName (old_varindex);
    assert (fullname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot determine tensor type for variable \"%s\"", fullname);
    free (fullname);
  }
  
  assert (var>=0 && var<tensortype->ncomps);
  
  groupindex = CCTK_GroupIndexFromVarI (old_varindex);
  assert (groupindex>=0);
  
  numvars = CCTK_NumVarsInGroupI (groupindex);
  assert (numvars>=0);
  assert (numvars > 0);
  
  firstvarindex = CCTK_FirstVarIndexI (groupindex);
  assert (firstvarindex>=0);
  assert (firstvarindex <= old_varindex);
  
  varindex = firstvarindex + var;
  
  /* scalars are a special case */
  if (tensortype == &TT_scalar) {
    varindex = old_varindex;
  }
  
  return varindex;
}



void
TT_component2indices (int const dim, int const rank,
                      int component,
                      int * restrict const indices)
{
  int ncomps;
  int r;
  
  assert (dim>=0);
  assert (rank>=0);
  ncomps = floor(pow(dim, rank) + 0.5);
  assert (component>=0 && component<ncomps);
  assert (rank==0 || indices);
  
#if 1
  for (r=0; r<rank; ++r) {
    indices[r] = component % dim;
    component /= dim;
  }
#else
  for (r=rank-1; r>=0; --r) {
    indices[r] = component % dim;
    component /= dim;
  }
#endif
  assert (component == 0);
}



int
TT_indices2component (int const dim, int const rank,
                      int const * restrict const indices)
{
  int component;
  int r;
  
  assert (dim>=0);
  assert (rank>=0);
  assert (rank==0 || indices);
  
  component = 0;
#if 1
  for (r=rank-1; r>=0; --r) {
    component = component * dim + indices[r];
  }
#else
  for (r=0; r<rank; ++r) {
    component = component * dim + indices[r];
  }
#endif
  
  return component;
}



int
TT_derivcode2num_derivs (int const derivcode)
{
  if (derivcode == 0) return 0;
  if (derivcode < 10) return 1;
  if (derivcode < 100) return 2;
  assert (0);
  return -1;
}



void
TT_derivcode2derivs (int derivcode, int const num_derivs,
                     int * restrict const indices)
{
  int r;
  
  assert (derivcode >= 0);
  assert (num_derivs == TT_derivcode2num_derivs (derivcode));
  if (num_derivs > 0) assert (indices);
  
  for (r=0; r<num_derivs; ++r) {
    indices[r] = derivcode % 10 - 1;
    derivcode /= 10;
  }
}



struct tensor const *
TT_derivative (struct tensor const * const tensortype, int const num_derivs)
{
  assert (tensortype);
  assert (num_derivs >= 0);
  
  switch (num_derivs) {
    
  case 0:
    return tensortype;
    
  case 1:
    if (tensortype == &TT_scalar    ) return &TT_vector;
    if (tensortype == &TT_vector    ) return &TT_tensor;
    if (tensortype == &TT_symmtensor) return &TT_symmtensor3a;
    break;
    
  case 2:
    if (tensortype == &TT_scalar    ) return &TT_symmtensor;
    if (tensortype == &TT_vector    ) return &TT_symmtensor3b;
    if (tensortype == &TT_symmtensor) return &TT_symmtensor4;
    break;
  }
  
  return NULL;
}
