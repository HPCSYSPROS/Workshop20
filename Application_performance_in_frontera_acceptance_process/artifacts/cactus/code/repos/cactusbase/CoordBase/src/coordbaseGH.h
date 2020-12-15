/*@@
  @header    coordbaseGH.h
  @date      Sun Sept 22 2002
  @author    Gabrielle Allen
  @desc
             The extensions to the GH structure from CoordBase
  @version   $Header$
@@*/

#ifndef _COORDBASE_COORDBASEGH_H_
#define _COORDBASE_COORDBASEGH_H_ 1

#include "util_Hash.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct COORDBASEGH {
  uHash *coordsystems;        /* Hash table mapping coord systems to tables */
  int *default_coord_systems; /* Default coord system for each dimension */
  /* int max_dimension;          Allocated size of default_coord_system
                                 array (not used currently) */
} coordbaseGH;

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* _COORDBASE_COORDBASEGH_H_ */
