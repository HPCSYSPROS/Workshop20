 /*@@
   @header    cctk_File.h
   @date      Mon Oct 18 00:42:35 CEST 1999
   @author    Gabrielle Allen
   @desc 
   Prototypes and constants for File routines
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTK_FILE_H_
#define _CCTK_FILE_H_

/* Prototypes */

#ifdef __cplusplus 
extern "C" {
#endif

int CCTK_CreateDirectory(int mode, const char *pathname);

#ifdef __cplusplus 
}
#endif

#endif



