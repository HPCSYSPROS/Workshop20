/*@@
  @file      MaskUtils.c
  @date      October 2002
  @author    Denis Pollney
  @desc
             Function prototypes and macros mask utilities.
  @desc
  @version   $Header$
@@*/

#ifndef __CACTUSEINSTEIN_SPACEMASK_H__
#define __CACTUSEINSTEIN_SPACEMASK_H__

#define SPACEMASK_DEBUG 0

/********************************************************************
 *********************  Routine Prototypes  *************************
 ********************************************************************/
#ifdef CCODE

#ifdef __cplusplus
extern "C"
{
#endif 

/* publicly visible routines */
int SpaceMask_RegisterType(const char* type_name, int nstates,
                           const char* const state_list[]);
int SpaceMask_AppendStatesToType(const char* type_name, int nstates,
                                 const char* const state_list[]);
CCTK_INT SpaceMask_GetTypeBits(const char* type_name);
CCTK_INT SpaceMask_GetStateBits(const char* type_name, const char* state_name);
void SpaceMask_SetState(CCTK_INT* mask, int point,
                        const char* type_name, const char* state);
CCTK_INT SpaceMask_CheckState(const CCTK_INT* mask, int point,
                         const char* type_name, const char* state);

#ifdef __cplusplus
}
#endif

#endif

/********************************************************************
 *********************     Local Data Types    **********************
 ********************************************************************/

#ifdef CCODE

typedef struct
{
  const char* name;
  CCTK_INT bitmask;
} SpaceMask_State;

typedef struct
{
  const char* name;
  int nstates;
  CCTK_INT bitmask;
  SpaceMask_State** state_list;
} SpaceMask_Type;

typedef struct
{
  int ntypes;
  SpaceMask_Type** type_list;
} SpaceMask_Registry;

extern SpaceMask_Registry* spacemask_registry;

#endif

/********************************************************************
 *********************     Macros      ******************************
 ********************************************************************/

/*@@
  @routine    SpaceMask_SetStateBits
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Sets the mask at a point to the given state, as specified
              using a bitmask.
  @enddesc 
@@*/

#ifdef FCODE
#define SpaceMask_SetStateBitsF90(mask, i,j,k, type_bits, state_bits) \
  mask(i,j,k) = ior(iand(mask(i,j,k), not(type_bits)), state_bits)
#endif

#ifdef CCODE
#define SpaceMask_SetStateBits(mask, ijk, type_bits, state_bits) \
  ((mask)[ijk] = (((mask)[ijk] & ~(type_bits)) | (state_bits)))
#endif

/*@@
  @routine    SpaceMask_CheckStateBits
  @author     Denis Pollney
  @date       15 October 2002
  @desc 
              Checks that the mask at a point has the specified state,
              in which case return 1, otherwise return 0.
  @enddesc 
@@*/

#ifdef FCODE
#define SpaceMask_CheckStateBitsF90(mask,i,j,k,type_bits,state_bits) \
  (iand(mask((i),(j),(k)),(type_bits)).eq.(state_bits))
#endif

#ifdef CCODE
#define SpaceMask_CheckStateBits(mask, ijk, type_bits, state_bits) \
  (((mask)[ijk] & (type_bits)) == (state_bits))
#endif

#endif /* __CACTUSEINSTEIN_SPACEMASK_H__ */
