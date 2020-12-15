 /*@@
   @header    cctk_MemAlloc.h
   @date      Thu Jan 20 2000
   @author    Gerd Lanfermann
   @desc
   Prototypes for Cactus MemAlloc functions.
   @enddesc
   @version $Header$
 @@*/

#ifndef _CCTK_MEMALLOC_H_
#define _CCTK_MEMALLOC_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

void *CCTK_Malloc(size_t size, int line, const char *file);
void  CCTK_Free(void *pointer);
unsigned long int CCTK_TotalMemory(void);

#ifdef __cplusplus
}
#endif

#endif /* _CCTK_MEMALLOC_H_ */


