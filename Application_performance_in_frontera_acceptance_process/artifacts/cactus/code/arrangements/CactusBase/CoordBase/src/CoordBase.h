/*@@
   @header  CoordBase.h
   @date    Wed June 26 2002
   @author  Gabrielle Allen
   @desc
            Defines and prototypes for the coordinate base thorn
   @enddesc
   @version $Header$
 @@*/

#ifndef _COORDBASE_H_
#define _COORDBASE_H_

#ifndef COORD_IN_COORDBASE

/* Since some people may still be including this header file
 * rather than using the function aliasing, we need to
 * allow them access to the internal functions.
 */

#define Coord_SystemRegister CoordBase_SystemRegister
#define Coord_SystemHandle CoordBase_SystemHandle
#define Coord_CoordRegister CoordBase_CoordRegister
#define Coord_CoordHandle CoordBase_CoordHandle
#define Coord_GroupSystem CoordBase_GroupSystem
#define Coord_SetDefaultSystem CoordBase_SetDefaultSystem
#define Coord_GetDefaultSystem CoordBase_GetDefaultSystem

#endif /* COORD_IN_COORDBASE */

#ifdef __cplusplus
extern "C" {
#endif

CCTK_INT CoordBase_SystemRegister(CCTK_POINTER_TO_CONST GH, CCTK_INT dim,
                                  CCTK_STRING systemname);

CCTK_INT CoordBase_SystemHandle(CCTK_POINTER_TO_CONST GH,
                                CCTK_STRING systemname);

CCTK_INT CoordBase_CoordRegister(CCTK_POINTER_TO_CONST GH,
                                 CCTK_INT systemhandle, CCTK_INT direction,
                                 CCTK_STRING coordname);

CCTK_INT CoordBase_CoordHandle(CCTK_POINTER_TO_CONST GH, CCTK_STRING coordname,
                               CCTK_STRING systemname);

CCTK_INT CoordBase_GroupSystem(CCTK_POINTER_TO_CONST GH, CCTK_STRING groupname);

CCTK_INT CoordBase_SetDefaultSystem(CCTK_POINTER_TO_CONST GH,
                                    CCTK_STRING systemname);

CCTK_INT CoordBase_GetDefaultSystem(CCTK_POINTER_TO_CONST GH,
                                    CCTK_INT systemdim);

#ifdef __cplusplus
}
#endif

#define COORDERROR_SYSTEMEXISTS -1
#define COORDERROR_INVALIDDIM -2
#define COORDERROR_INVALIDNAME -3
#define COORDERROR_TABLEERROR -4
#define COORDERROR_NOSYSTEM -5
#define COORDERROR_INVALIDHANDLE -6
#define COORDERROR_COORDINATEEXISTS -7
#define COORDERROR_DUPLICATENAME -8
#define COORDERROR_NOSUCHCOORD -9
#define COORDERROR_INVALIDGROUPNAME -10
#define COORDERROR_NOCOORDSYS -11
#define COORDERROR_NODIMENSION -12
#define COORDERROR_DEFAULTEXISTS -13

#endif /* _COORDBASE_H_ */
