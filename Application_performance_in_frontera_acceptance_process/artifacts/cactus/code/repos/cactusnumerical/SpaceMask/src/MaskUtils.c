/*@@
  @file      MaskUtils.c
  @date      October 2002
  @author    Denis Pollney
  @desc
             Utilities for registering/setting/checking the mask.
  @desc
  @version   $Header$
@@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_FortranString.h"

#include "SpaceMask.h"

static const char* rcsid = "$Header$";

CCTK_FILEVERSION(CACTUSEINSTEIN_SPACEMASK_MaskUtils_c);

SpaceMask_Registry* spacemask_registry = NULL;

/* prototypes for functions local to this file */
static
  int SpaceMask_get_bit_nbr(int nstates);
static
  CCTK_INT SpaceMask_get_free_bits(int nbits);
static
  CCTK_INT* SpaceMask_determine_state_mask(CCTK_INT allocated_bits,
                                           int nstates);
static
  SpaceMask_State* SpaceMask_setup_new_state(CCTK_INT state_mask,
                                             const char* name);
static
  SpaceMask_Type*
  SpaceMask_setup_new_type(CCTK_INT new_bits, const char* type_name,
                           int nstates, const char* const state_list[], 
                           const CCTK_INT state_mask[]);
static
  void SpaceMask_append_type_to_registry(SpaceMask_Type* new_type);

/* prototypes for Fortran wrapper functions defined in this file */
void CCTK_FCALL CCTK_FNAME(SpaceMask_GetTypeBits)(CCTK_INT* type_bits,
                                                  ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(SpaceMask_GetStateBits)(CCTK_INT* state_bits,
                                                   TWO_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(SpaceMask_SetState)(CCTK_INT* mask,
                                               CCTK_INT* point, TWO_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(SpaceMask_CheckState)(int* retval,
                                                 CCTK_INT* mask,
                                                 CCTK_INT* point,
                                                 TWO_FORTSTRING_ARG);

/*@@
  @routine    SpaceMask_get_bit_nbr
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Return the number of bits required to represent a given number.
              Cheesy divide-by-two method, maybe something else is quicker.
  @enddesc 
@@*/
static 
int
SpaceMask_get_bit_nbr(int nstates)
{
  int bstates;
  int bit_nbr;

  bstates = nstates-1;

  for (bit_nbr=0; bstates>0; ++bit_nbr)
    bstates /= 2;

  return bit_nbr;
}

/*@@
  @routine    SpaceMask_get_free_bits
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Determine the mask bits which have not yet been allocated.
              The return value is a bitmask with the requested number of
              free bits set to 1. If there are not enough free bits left
              to satisfy the request, stop with an error.
  @enddesc 
@@*/
static
CCTK_INT
SpaceMask_get_free_bits(int nbits)
{
  CCTK_INT used_bits;
  CCTK_INT new_bits;
  int i;
  int j;
  int n;

  used_bits = 0;
  if (spacemask_registry != NULL)
    {
      for (i=0; i<spacemask_registry->ntypes; ++i)
        used_bits |= spacemask_registry->type_list[i]->bitmask;
    }

  n=1;
  new_bits = 0;
  j = 0;
  /* Do not use the 2 top bits.  Signed integer datatypes in C do not
     have the overflow semantics that many people expect.  With care,
     we could use the second topmost bit as well, but this doesn't
     seem important at the moment.  */
  for (i=0; i < (int) sizeof(CCTK_INT)*8-2 && j<nbits; ++i)
    {
      if (!(n & used_bits))
        {
          ++j;
          new_bits |= n;
        }
      n *= 2;
    }

  if (j<nbits)
    {
      CCTK_WARN (1, "Cannot allocate mask: Not enough free bits.");
      return 0;
    }

  return new_bits;
}

/*@@
  @routine    SpaceMask_determine_state_mask
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Determine appropriate bitmasks to represent a number of
              states, using the allocated bit mask.
              That is, if allocated_bits is set to 00111000, then
              the returned list of states are permutations of the three
              non-zero (active) bits.
  @enddesc 
@@*/
static
CCTK_INT*
SpaceMask_determine_state_mask(CCTK_INT allocated_bits, int nstates)
{
  CCTK_INT* state_mask;
  int n;
  int bit;
  int i;
  int j;

  state_mask =
    (CCTK_INT*) malloc (nstates*sizeof(CCTK_INT));

  for (j=0; j<nstates; ++j)
    state_mask[j] = 0;

  n=1;
  bit=1;

  /* Do not use the 2 top bits.  Signed integer datatypes in C do not
     have the overflow semantics that many people expect.  With care,
     we could use the second topmost bit as well, but this doesn't
     seem important at the moment.  */
  for (i=0; i < (int) sizeof(CCTK_INT)*8-2; ++i)
    {
      if (n & allocated_bits)
        {
          for (j=0; j<nstates; ++j)
            {
              if (bit & j)
                state_mask[j] += n;
            }
          bit *= 2;
        }
      n *= 2;
    }

  return state_mask;
}

/*@@
  @routine    SpaceMask_setup_new_state
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Allocate a new SpaceMask_State, giving it a name and
              a bit mask.
  @enddesc 
@@*/
static
SpaceMask_State*
SpaceMask_setup_new_state(CCTK_INT state_mask, const char* name)
{
  SpaceMask_State* new_state;

  new_state = (SpaceMask_State*) malloc (sizeof(SpaceMask_State));
  new_state->name = name;
  new_state->bitmask = state_mask;

  return new_state;
}

/*@@
  @routine    SpaceMask_setup_new_type
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Allocate a new SpaceMask_Type. This involves allocating
              a set of bits within the mask to the new type, and creating
              appropriate bitmasks (using the allocated bits) for each of
              the requested states.
  @enddesc 
@@*/
static
SpaceMask_Type*
SpaceMask_setup_new_type(CCTK_INT new_bits, const char* type_name,
                         int nstates, const char* const state_list[], 
                         const CCTK_INT state_mask[])
{
  SpaceMask_Type* new_type;
  int j;

  new_type = (SpaceMask_Type*) malloc (sizeof(SpaceMask_Type));

  if (new_type == NULL)
    return NULL;

  new_type->bitmask = new_bits;
  new_type->nstates = nstates;
  new_type->name = type_name;
  new_type->state_list =
    (SpaceMask_State**) malloc (nstates*sizeof(SpaceMask_State*));

  if (state_list==NULL)
    {
      for (j=0; j<nstates; ++j)
        new_type->state_list[j] = SpaceMask_setup_new_state(state_mask[j],
                                                            NULL);
    }
  else
    {
      for (j=0; j<nstates; ++j)
        new_type->state_list[j] = SpaceMask_setup_new_state(state_mask[j],
                                                            state_list[j]);
    }

  return new_type;
}

/*@@
  @routine    SpaceMask_append_type_to_registry
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Adds a new type to the spacemask_registry.
  @enddesc 
@@*/
static
void
SpaceMask_append_type_to_registry(SpaceMask_Type* new_type)
{
  SpaceMask_Type** new_type_list;
  int ntypes;
  int i;

  if (spacemask_registry == NULL)
    {
      spacemask_registry =
        (SpaceMask_Registry*) malloc (sizeof(SpaceMask_Registry));
      ntypes = 1;
      new_type_list = (SpaceMask_Type**) malloc (sizeof(SpaceMask_Type*));
      new_type_list[0] = new_type;
    }
  else
    {
      ntypes = spacemask_registry->ntypes + 1;
      new_type_list = 
        (SpaceMask_Type**) malloc (ntypes*sizeof(SpaceMask_Type*));
      for (i=0; i<ntypes-1; ++i)
        new_type_list[i] = spacemask_registry->type_list[i];
      new_type_list[ntypes-1] = new_type;
      free(spacemask_registry->type_list);
    }

  spacemask_registry->ntypes = ntypes;
  spacemask_registry->type_list = new_type_list;
}

/*@@
  @routine    SpaceMask_RegisterType
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Allocates a set of bits of the SpaceMask to represent
              a set of states of the named type.
  @enddesc 
@@*/

int
SpaceMask_RegisterType(const char* type_name, int nstates,
                       const char* const state_list[])
{
  SpaceMask_Type* new_type;
  CCTK_INT new_bits;
  CCTK_INT* state_mask;
  int bit_nbr;
 
  bit_nbr = SpaceMask_get_bit_nbr(nstates);
  new_bits = SpaceMask_get_free_bits(bit_nbr);

  if (new_bits == 0)
    return -1;

  state_mask = SpaceMask_determine_state_mask(new_bits, nstates);

  new_type = SpaceMask_setup_new_type(new_bits, type_name, nstates, 
                                      state_list, state_mask);
  if (new_type == NULL)
    return -1;

  SpaceMask_append_type_to_registry (new_type);

  return 0;
}

/*@@
  @routine    SpaceMask_AppendStatesToType
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Adds a new set of possible states to an already existing
              type. If required, new bits of the mask are allocated to the
              given state. If bits are not available for this, then an
              error is returned.
  @enddesc 
@@*/

int
SpaceMask_AppendStatesToType(const char* type_name, int nstates,
                             const char* const state_list[])
{
  CCTK_INT new_bits;
  CCTK_INT allocated_bits;
  CCTK_INT* state_mask;
  SpaceMask_Type* old_type;
  SpaceMask_Type* new_type;
  const char** new_state_list;
  int i;
  int j;
  int old_type_idx;
  int new_bit_nbr;
  int total_nstates;
  
  old_type = NULL;
  for (i=0; i<spacemask_registry->ntypes; ++i)
    {
      if (!strcmp(spacemask_registry->type_list[i]->name, type_name))
        {
          old_type = spacemask_registry->type_list[i];
          old_type_idx = i;
        }
    }

  if (old_type == NULL)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Mask type \"%s\" has not been registered", type_name);
      return -1;
    }

  total_nstates = old_type->nstates + nstates;
  new_bit_nbr = SpaceMask_get_bit_nbr(total_nstates) - 
    SpaceMask_get_bit_nbr(old_type->nstates);

  new_bits = SpaceMask_get_free_bits(new_bit_nbr);

  if (new_bits == 0)
    return -1;

  allocated_bits = old_type->bitmask | new_bits;

  state_mask = SpaceMask_determine_state_mask(allocated_bits, total_nstates);

  new_state_list = (const char**) malloc (total_nstates*sizeof(char*));

  i = 0;
  for (j=0; j<old_type->nstates; ++j, ++i)
    new_state_list[i] = old_type->state_list[j]->name;
  for (j=0; j<nstates; ++j, ++i)
    new_state_list[i] = state_list[j];

  new_type = SpaceMask_setup_new_type(allocated_bits, type_name, total_nstates,
                                      new_state_list, state_mask);
  if (new_type == NULL)
    return -1;

  free(old_type);

  spacemask_registry->type_list[old_type_idx] = new_type;

  return 0;
}

/*@@
  @routine    SpaceMask_GetTypeBits
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Returns the bitmask corresponding to the given type. Bits
              which are allocated to the given type are set to 1, all other
              bits are 0.
  @enddesc 
@@*/

CCTK_INT
SpaceMask_GetTypeBits(const char* type_name)
{
  int i;

  if (spacemask_registry != NULL)
    {
      for (i=0; i<spacemask_registry->ntypes; ++i)
        {
          if (!strcmp(spacemask_registry->type_list[i]->name, type_name))
            return spacemask_registry->type_list[i]->bitmask;
        }
    }

  CCTK_VInfo (CCTK_THORNSTRING, "Type \"%s\" has not been registered.\n",
              type_name);

  return 0;
}

/*@@
  @routine    SpaceMask_GetStateBits
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Returns the bitmask corresponding to the given state, or
              -1 if the state does not exist.
  @enddesc 
@@*/

CCTK_INT
SpaceMask_GetStateBits(const char* type_name, const char* state_name)
{
  SpaceMask_Type* type;
  int i, j;

  if (spacemask_registry != NULL)
    {
      for (i=0; i<spacemask_registry->ntypes; ++i)
        {
          if (!strcmp(spacemask_registry->type_list[i]->name, type_name))
            {
              type = spacemask_registry->type_list[i];
              for (j=0; j<type->nstates; ++j)
                {
                  if (!strcmp(type->state_list[j]->name, state_name))
                    return type->state_list[j]->bitmask;
                }
            }
	}
    }

  CCTK_VInfo (CCTK_THORNSTRING, 
              "Requested state \"%s\" could not be found in type \"%s\"\n",
              state_name, type_name);

  return -1;
}

/*@@
  @routine    SpaceMask_GetStateBitsList
  @author     Denis Pollney
  @date       25 May 2003
  @desc 
              Returns a list of all of the state masks allocated
              for a given type. The number of elements in the list
              is equal to the number of states of the type.
              If the type name is not recognised, NULL is returned.
  @enddesc 
@@*/

CCTK_INT*
SpaceMask_GetStateBitsList(const char* type_name)
{
  SpaceMask_Type* type;
  CCTK_INT* state_bits_list;
  int i, j;

  for (i=0; i<spacemask_registry->ntypes; ++i)
    {
      if (!strcmp(spacemask_registry->type_list[i]->name, type_name))
        {
          type = spacemask_registry->type_list[i];
          state_bits_list = (CCTK_INT*) 
            malloc (sizeof(CCTK_INT)*(type->nstates));

          for (j=0; j<type->nstates; ++j)
            state_bits_list[j] = type->state_list[j]->bitmask;

          return state_bits_list;
        }
    }

  CCTK_VInfo (CCTK_THORNSTRING, 
              "Requested type \"%s\" could not be found\n", type_name);

  return NULL;
}

/*@@
  @routine    SpaceMask_SetState
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Sets the mask at a point to the given state.
  @enddesc 
@@*/

void
SpaceMask_SetState(CCTK_INT* mask, int point, const char* type_name, const char* state)
{
  CCTK_INT type_bits;
  CCTK_INT state_bits;

  type_bits = SpaceMask_GetTypeBits(type_name);
  state_bits = SpaceMask_GetStateBits(type_name, state);

  SpaceMask_SetStateBits(mask, point, type_bits, state_bits);
}

/*@@
  @routine    SpaceMask_CheckState
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Checks that the mask at a point has the given state, in which
              case return 1, otherwise return 0.
  @enddesc 
@@*/

CCTK_INT
SpaceMask_CheckState(const CCTK_INT* mask, int point, const char* type_name, const char* state)
{
  CCTK_INT type_bits;
  CCTK_INT state_bits;

  type_bits = SpaceMask_GetTypeBits(type_name);
  state_bits = SpaceMask_GetStateBits(type_name, state);

  return SpaceMask_CheckStateBits(mask, point, type_bits, state_bits);
}


/********************************************************************
 *********************  Fortran Wrappers  ***************************
 ********************************************************************/

/* How do you pass the char** state list from fortran???
void
CCTK_FCALL CCTK_FNAME(SpaceMask_RegisterType)(int* ierr, char* type_name,
                                              int* nstates, char** state_list)
{
  *ierr = SpaceMask_RegisterType(type_name, *nstates, state_list);
}
*/

void
CCTK_FCALL CCTK_FNAME(SpaceMask_GetTypeBits)(CCTK_INT* type_bits,
                                             ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(type_name)

  *type_bits = SpaceMask_GetTypeBits(type_name);

  free(type_name);
}

void
CCTK_FCALL CCTK_FNAME(SpaceMask_GetStateBits)(CCTK_INT* state_bits,
                                              TWO_FORTSTRING_ARG)
{

  TWO_FORTSTRING_CREATE(type_name, state_name)

  *state_bits = SpaceMask_GetStateBits(type_name, state_name);

  free(type_name);
  free(state_name);
}

void
CCTK_FCALL CCTK_FNAME(SpaceMask_SetState)(CCTK_INT* mask,
                                          CCTK_INT* point, TWO_FORTSTRING_ARG)
{
  TWO_FORTSTRING_CREATE(type_name, state)

  SpaceMask_SetState(mask, *point, type_name, state);

  free(type_name);
  free(state);
}

void
CCTK_FCALL CCTK_FNAME(SpaceMask_CheckState)(int* retval,
                                            CCTK_INT* mask,
                                            CCTK_INT* point,
                                            TWO_FORTSTRING_ARG)
{
  TWO_FORTSTRING_CREATE(type_name, state)

  *retval = SpaceMask_CheckState(mask, *point, type_name, state);

  free(type_name);
  free(state);
}
