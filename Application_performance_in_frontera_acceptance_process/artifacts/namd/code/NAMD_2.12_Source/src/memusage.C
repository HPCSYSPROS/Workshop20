/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "converse.h"
#ifndef WIN32
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#else
int sbrk(int) { return 0; }
#endif

int memusageinit::initialized;
unsigned long memusageinit::sbrkval;

memusageinit::memusageinit() {
  if ( initialized == 0 ) {
    sbrkval = (unsigned long) sbrk(0);
    initialized = 1;
  }
}

unsigned long memusageinit::memusage_sbrk() {
  unsigned long newval = (unsigned long) sbrk(0);
  return ( newval - memusageinit::sbrkval );
}


#ifdef WIN32
#define MEMUSAGE_USE_SBRK
#endif

#ifdef _NO_MALLOC_H
#define NO_MALLINFO
#endif

#ifdef MEMUSAGE_USE_SBRK
#define NO_MALLINFO
#define NO_PS
#endif


#ifdef MEMUSAGE_USE_MSTATS

#include <malloc/malloc.h>

unsigned long memusage_mstats() {
  struct mstats ms = mstats();
  unsigned long memtotal = ms.bytes_used;
  return memtotal;
}

#else
inline unsigned long memusage_mstats() { return 0; }
#endif


#ifndef NO_MALLINFO

#include <malloc.h>

unsigned long memusage_mallinfo() {
  struct mallinfo mi = mallinfo();
  // unsigned long memtotal = mi.usmblks + mi.uordblks + mi.hblkhd;
  unsigned long memtotal = (unsigned int) mi.uordblks;
  unsigned long memtotal2 = (unsigned int) mi.usmblks;
  memtotal2 += (unsigned int) mi.hblkhd;
  if ( memtotal2 > memtotal ) memtotal = memtotal2;

  // printf("mallinfo %d %d %d\n", mi.usmblks, mi.uordblks, mi.hblkhd);

  return memtotal;
}

#else
inline unsigned long memusage_mallinfo() { return 0; }
#endif


inline unsigned long memusage_ps() {
#ifdef NO_PS
  return 0;
#else
  char pscmd[100];
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  sprintf(pscmd, "/bin/ps -o rss= -p %d", getpid());
#else
  sprintf(pscmd, "/bin/ps -o vsz= -p %d", getpid());
#endif
  unsigned long vsz = 0;
  FILE *p = popen(pscmd, "r");
  if ( p ) {
    fscanf(p, "%ld", &vsz);
    pclose(p);
  }
  return ( vsz * (unsigned long) 1024 );
#endif
}

#if CMK_BLUEGENEQ
/* Report the memusage according to
 * https://wiki.alcf.anl.gov/parts/index.php/Blue_Gene/Q#CNK_Memory_Information
 */
#include <spi/include/kernel/memory.h>
inline long long memusage_bgq(){
    uint64_t heap;
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
    return heap;
}
#endif

#if CMK_BLUEGENEP
/* Report the memusage according to
 * http://www.alcf.anl.gov/resource-guides/determining-memory-use
 */
#include <malloc.h>
inline long long memusage_bgp(){
    struct mallinfo m = mallinfo();
    long long hblkhd = (unsigned int) m.hblkhd;
    long long uordblks = (unsigned int) m.uordblks;
    return hblkhd + uordblks;
}
#endif

inline unsigned long memusage_proc_self_stat() {
#ifdef NO_PS
  return 0;
#else
  static int failed_once = 0;
  if ( failed_once ) return 0;  // no point in retrying

  FILE *f = fopen("/proc/self/stat","r");
  if ( ! f ) { failed_once = 1; return 0; }
  for ( int i=0; i<22; ++i ) fscanf(f,"%*s");
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  // skip vss, next value is rss in pages
  fscanf(f,"%*s");
#endif
  unsigned long vsz = 0;  // should remain 0 on failure
  fscanf(f,"%lu",&vsz);
  fclose(f);
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  vsz *= sysconf(_SC_PAGESIZE);
#endif
  if ( ! vsz ) failed_once = 1;
  // printf("/proc/self/stat reports %d MB\n", vsz/(1024*1024));
  return vsz;
#endif
}


unsigned long memusage(const char **source) {

  unsigned long memtotal = 0;
  const char* s = "ERROR";

  if ( ! CmiMemoryIs(CMI_MEMORY_IS_OS) ) {
    memtotal = CmiMemoryUsage();  s = "CmiMemoryUsage";
  }

#if CMK_BLUEGENEQ
  if( ! memtotal) { memtotal = memusage_bgq(); s="Kernel_GetMemorySize on BG/Q"; }
#endif

#if CMK_BLUEGENEP
  if( ! memtotal) { memtotal = memusage_bgp(); s="mallinfo on BG/P"; }
#endif

#if defined(WIN32) && !defined(__CYGWIN__)
  if ( ! memtotal ) {
    memtotal = CmiMemoryUsage();  s = "GetProcessMemoryInfo";
  }
#endif

  if ( ! memtotal ) {
    memtotal = memusage_proc_self_stat();  s = "/proc/self/stat";
  }

  if ( ! memtotal ) { memtotal = memusage_mstats(); s = "mstats"; }

  if ( ! memtotal ) { memtotal = memusage_mallinfo(); s = "mallinfo"; }

  if ( ! memtotal ) { memtotal = memusageinit::memusage_sbrk(); s = "sbrk"; }

  if ( ! memtotal ) { memtotal = memusage_ps(); s = "ps"; }

  if ( ! memtotal ) { memtotal = CmiMemoryUsage();  s = "CmiMemoryUsage"; }

  if ( ! memtotal ) s = "nothing";

  if ( source ) *source = s;

  return memtotal;

}

