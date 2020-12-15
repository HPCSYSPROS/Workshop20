 /*@@
   @header    cctk_Constants.h
   @date      Fri Oct 15 21:29:23 CEST 1999
   @author    Gabrielle Allen
   @desc 
   Constants used by Cactus
   @enddesc
   @version $Header$
 @@*/

#ifndef _CCTK_CONSTANTS_H_
#define _CCTK_CONSTANTS_H_

#define CCTK_VARIABLE_VOID             100
#define CCTK_VARIABLE_BYTE             110
#define CCTK_VARIABLE_INT              120
#define CCTK_VARIABLE_INT1             121
#define CCTK_VARIABLE_INT2             122
#define CCTK_VARIABLE_INT4             123
#define CCTK_VARIABLE_INT8             124
#define CCTK_VARIABLE_INT16            125
#define CCTK_VARIABLE_REAL             130
#define CCTK_VARIABLE_REAL4            131
#define CCTK_VARIABLE_REAL8            132
#define CCTK_VARIABLE_REAL16           133
#define CCTK_VARIABLE_COMPLEX          140
#define CCTK_VARIABLE_COMPLEX8         141
#define CCTK_VARIABLE_COMPLEX16        142
#define CCTK_VARIABLE_COMPLEX32        143
#define CCTK_VARIABLE_CHAR             150
#define CCTK_VARIABLE_STRING           151
#define CCTK_VARIABLE_POINTER          160
#define CCTK_VARIABLE_POINTER_TO_CONST 161
#define CCTK_VARIABLE_FPOINTER         162

/* DEPRECATED IN BETA 12 */
#define CCTK_VARIABLE_FN_POINTER CCTK_VARIABLE_FPOINTER

/* steerable status of parameters */
#define CCTK_STEERABLE_NEVER   200
#define CCTK_STEERABLE_ALWAYS  201
#define CCTK_STEERABLE_RECOVER 202

/* group distributions */
#define CCTK_DISTRIB_CONSTANT 301
#define CCTK_DISTRIB_DEFAULT  302

/* group types */
#define CCTK_SCALAR 401
#define CCTK_GF     402
#define CCTK_ARRAY  403

/* group scopes */
#define CCTK_PRIVATE   501
#define CCTK_PROTECTED 502
#define CCTK_PUBLIC    503

/* constants for CCTK_TraverseString() */
#define CCTK_VAR          601
#define CCTK_GROUP        602
#define CCTK_GROUP_OR_VAR 603


#endif /* _CCTK_CONSTANTS_ */

