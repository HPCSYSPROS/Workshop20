 /*@@
   @header    cctk_WarnLevel.h
   @date      Wed Feb 17 00:53:55 1999
   @author    Tom Goodale
   @desc 
   Header for the warning functions.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTK_WARNLEVEL_H_
#define _CCTK_WARNLEVEL_H_

#include <cctk_Config.h>

#ifdef CCODE

#ifdef __cplusplus 
extern "C" 
{
#endif

int CCTK_Warn(int level, 
               int line, 
               const char *file, 
               const char *thorn, 
               const char *message);
int CCTK_VWarn(int level, 
                int line, 
                const char *file, 
                const char *thorn, 
                const char *format, ...)
  CCTK_ATTRIBUTE_FORMAT(printf, 5, 6);

void CCTK_Error(int line, 
                const char *file, 
                const char *thorn, 
                const char *message)
  CCTK_ATTRIBUTE_NORETURN;

void CCTK_VError(int line, 
                 const char *file, 
                 const char *thorn, 
                 const char *format, ...)
  CCTK_ATTRIBUTE_FORMAT(printf, 4, 5)
  CCTK_ATTRIBUTE_NORETURN;
int CCTK_VParamWarn (const char *thorn,
                     const char *format,
                     ...)
  CCTK_ATTRIBUTE_FORMAT(printf, 2, 3);
int CCTK_ParamWarn(const char *thorn, const char *message);
int CCTK_Info(const char *thorn, const char *message);
int CCTK_VInfo(const char *thorn, const char *format, ...)
  CCTK_ATTRIBUTE_FORMAT(printf, 2, 3);

/* prototypes for warn/info callback routines */

typedef void (*cctk_warnfunc)(int level,
                              int line,
                              const char *file, 
                              const char *thorn,
                              const char *message,
                              void *data);
                            
typedef void (*cctk_infofunc)(const char *thorn,
                              const char *message,
                              void *data);                            

/* prototypes for warn/info registration routines */
int CCTK_WarnCallbackRegister(int minlevel,
                              int maxlevel,
                              void *data,
                              cctk_warnfunc callback);
                              

int CCTK_InfoCallbackRegister(void *data, 
                              cctk_infofunc callback);
                              

#ifdef __cplusplus 
}
#endif

#endif  /* CCODE */

/* suggested values for warning levels (courtesy of Steve, PR#1742) */
#define CCTK_WARN_ABORT    0    /* abort the Cactus run */
#define CCTK_WARN_ALERT    1    /* the results of this run will be wrong, */
                                /* and this will surprise the user, */
                                /* but we can still continue the run */
#define CCTK_WARN_COMPLAIN 2    /* the user should know about this, */
                                /* but the problem is not terribly */
                                /* surprising */
#define CCTK_WARN_PICKY    3    /* this is for small problems that can */
                                /* probably be ignored, but that careful */
                                /* people may want to know about */
#define CCTK_WARN_DEBUG    4    /* these messages are probably useful */
                                /* only for debugging purposes */

#endif  /* ! _CCTK_WARNLEVEL_H_ */
