#ifndef CARPET_PUBLIC_H
#define CARPET_PUBLIC_H

#include <cctk.h>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

/* Tell thorns that the Carpet routines exist */
#define HAVE_CARPET

#ifdef __cplusplus
namespace Carpet {
extern "C" {
#endif

/* Carpet's GH extension */
struct CarpetGH {

#if 0
      /* Maximum number of refinement levels */
      int maxreflevels;
      
      /* Refinement levels */
      int reflevels;
      
      /* Refinement factor */
      int reffact;
      
      /* Refinement factor on finest possible grid */
      int maxreflevelfact;
      
      /* Base multigrid level */
      int basemglevel;
      
      /* Multigrid levels */
      int mglevels;
      
      /* Multigrid factor */
      int mgfact;
      
      /* Multigrid factor on coarsest grid */
      int maxmglevelfact;
#endif

  /* Maps */
  int maps;

#if 0
      /* Current position on the grid hierarchy */
      int reflevel;
      int mglevel;
#endif
  int map;
#if 0
      int component;
      
      /* Current refinement factor */
      int reflevelfact;
      
      /* Current multigrid factor */
      int mglevelfact;
#endif

#if 0
      /* Number of buffer zones */
      int const * nbufferzones;
#endif
};

struct CarpetGH const *GetCarpetGH(const cGH *const cgh) CCTK_ATTRIBUTE_PURE;

/* Prolongation management */
CCTK_INT CarpetEnableProlongating(const CCTK_INT flag);
CCTK_INT CarpetQueryProlongating(void);

/* Grid function access */
CCTK_POINTER Carpet_VarDataPtrI(CCTK_POINTER_TO_CONST const cctkGH,
                                CCTK_INT const m, CCTK_INT const rl,
                                CCTK_INT const c, CCTK_INT const tl,
                                CCTK_INT const varindex);

/* Call a schedule group */
int CallScheduleGroup(cGH *const cgh, const char *const group);

/* Call a local function */
int CallLocalFunction(cGH *const cgh, void (*const function)(cGH *const cgh));
int CallSinglemapFunction(cGH *const cgh,
                          void (*const function)(cGH *const cgh));
int CallLevelFunction(cGH *const cgh, void (*const function)(cGH *const cgh));
int CallGlobalFunction(cGH *const cgh, void (*const function)(cGH *const cgh));
int CallMetaFunction(cGH *const cgh, void (*const function)(cGH *const cgh));

/* Helper functions */
MPI_Comm CarpetMPIComm(void) CCTK_ATTRIBUTE_PURE;
MPI_Datatype CarpetMPIDatatype(int vartype) CCTK_ATTRIBUTE_PURE;
MPI_Datatype CarpetSimpleMPIDatatype(int vartype) CCTK_ATTRIBUTE_PURE;
int CarpetSimpleMPIDatatypeLength(int vartype) CCTK_ATTRIBUTE_PURE;

#ifdef __cplusplus
} /* extern "C" */
} /* namespace Carpet */
#endif

#endif /* !defined(CARPET_PUBLIC_H) */
