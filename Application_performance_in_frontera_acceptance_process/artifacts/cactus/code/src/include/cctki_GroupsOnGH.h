 /*@@
   @header    cctki_GroupsOnGH.h
   @date      2012-10-26
   @author    Erik Schnetter
   @desc 
   Prototypes and constants for internal group functions which use GH structure.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTKI_GROUPSONGH_H_
#define _CCTKI_GROUPSONGH_H_

/* Prototypes */

#ifdef __cplusplus 
extern "C" {
#endif

void *CCTKi_VarDataPtrI(const cGH *GH, int timelevel, int varindex);

#ifdef __cplusplus 
}
#endif

#endif
