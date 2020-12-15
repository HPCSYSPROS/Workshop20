 /*@@
   @header    cctk_Version.h
   @date      Fri Dec 15 18:44:56 2000
   @author    Tom Goodale
   @desc 
   Version info about this Cactus
   @enddesc
   @version $Header$
 @@*/

#ifndef _CCTK_VERSION_H_
#define _CCTK_VERSION_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

/* Version Stuff */
const char *CCTK_FullVersion(void);

const char *CCTK_MajorVersion(void);

const char *CCTK_MinorVersion(void);

const char *CCTK_MicroVersion(void);

/* Compile date and time stuff */

void CCTKi_DateStamp(void); 

const char *CCTK_CompileTime(void);

const char *CCTK_CompileDate(void);

const char *CCTK_CompileDateTime(void);

#ifdef __cplusplus
}
#endif

#endif /* _CCTK_VERSION_H_ */
