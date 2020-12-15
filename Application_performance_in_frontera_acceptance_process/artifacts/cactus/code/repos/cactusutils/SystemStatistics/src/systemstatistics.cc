
#include <stdio.h> 
#include <string.h> 
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include "assert.h"

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

// There doesn't seem to be a Cactus macro defined for the Mach header
// files we actually need, so use the mach_time.h macro instead
#ifdef HAVE_MACH_MACH_TIME_H
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

#ifndef HAVE_MALLINFO

// Provide a dummy mallinfo function if none is available
struct mallinfo_t {
  int arena;
  int ordblks;
  int smblks;
  int hblks;
  int hblkhd;
  int usmblks;
  int fsmblks;
  int uordblks;
  int fordblks;
  int keepcost;
};

struct mallinfo_t mallinfo()
{
  struct mallinfo_t m;
  m.arena = 0;
  m.ordblks = 0;
  m.smblks = 0;
  m.hblks = 0;
  m.hblkhd = 0;
  m.usmblks = 0;
  m.fsmblks = 0;
  m.uordblks = 0;
  m.fordblks = 0;
  m.keepcost = 0;
  return m;
}

#endif

#ifndef _MACH_INIT_

static unsigned long int get_rss()
{
  unsigned int size=0; //       total program size
  unsigned int resident=0;//   resident set size
  unsigned int share=0;//      shared pages
  unsigned int text=0;//       text (code)
  unsigned int lib=0;//        library
  unsigned int data=0;//       data/stack

  int page_size = sysconf(_SC_PAGESIZE);

  char buf[30];
  snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  // If the /proc filesystem does not exist, this file will not be
  // found and the function will return 0
  if (pf) 
  {
    if (fscanf(pf, "%u %u %u %u %u %u",
               &size, &resident, &share, &text, &lib, &data) != 6)
    {
      CCTK_WARN(1, "Error while reading memory statistics (rss); results will be invalid");
      fclose(pf);
      return 0;
    }
    
    fclose(pf);
  }
  return (unsigned long int ) resident * (unsigned long int) page_size;
}

#else

// The code to get the RSS from Mac OS has been modified from
//   http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html

static unsigned long int get_rss()
{
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    task_t task = MACH_PORT_NULL;

    if (task_for_pid(current_task(), getpid(), &task) != KERN_SUCCESS)
        abort();

    task_info(task, TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
    return t_info.resident_size;
}

#endif

static unsigned int get_majflt()
{
  int pid;
  char exe[256];
  char  state;
  int dummyi;
  unsigned long int dummyu;
  unsigned long int majflt;

  unsigned int page_size = sysconf(_SC_PAGESIZE);

  char buf[30];
  snprintf(buf, 30, "/proc/%u/stat", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  // If the /proc filesystem does not exist, this file will not be
  // found and the function will return 0
  if (pf) 
  {
    if (fscanf(pf, "%d %s %c %d %d %d %d %d %lu %lu %lu %lu", 
           &pid, exe, &state, &dummyi,  &dummyi, &dummyi, &dummyi, 
               &dummyi, &dummyu, &dummyu, &dummyu, &majflt) != 12)
    {
      CCTK_WARN(1, "Error while reading memory statistics (majflt); results will be invalid");
      fclose(pf);
      return 0;
    }

    fclose(pf);
  }
  return majflt * page_size;
}

static long long int get_swap_kB()
{
  FILE *f = fopen("/proc/meminfo", "r");
  if (f == 0)
    return -1;
  const int buf_len = 100;
  char buffer[buf_len];
  long long int swap_total = 0;
  long long int swap_free = 0;
  bool read_swap_total = false;
  bool read_swap_free = false;

  while(!feof(f) && fgets(buffer, buf_len, f) != NULL)
  {
    char key[100];
    char unit[100];
    long long int val = -1;
    sscanf(buffer, "%s %lld %s", key, &val, unit);
    if (strcmp(key, "SwapTotal:") == 0)
    {
      assert(strcmp(unit, "kB") == 0);
      swap_total = val;
      read_swap_total = true;
    }
    else if (strcmp(key, "SwapFree:") == 0)
    {
      assert(strcmp(unit, "kB") == 0);
      swap_free = val;
      read_swap_free = true;
    }
  }

  if (!read_swap_total || !read_swap_free)
  {
    CCTK_WARN(1, "Unable to read swap usage from /proc/meminfo");
    swap_total = 0; swap_free = 0;
  }

  fclose(f);
  return swap_total - swap_free;
}

extern "C" void SystemStatistics_Collect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const int mb = 1024*1024;
  const int kb = 1024;

  *maxrss = get_rss();
  *majflt = get_majflt();
  *arena = mallinfo().arena;
  *ordblks = mallinfo().ordblks;
  *hblks = mallinfo().hblks;
  *hblkhd = mallinfo().hblkhd;
  *uordblks = mallinfo().uordblks;
  *fordblks = mallinfo().fordblks;
  *keepcost = mallinfo().keepcost;
  *swap_used = get_swap_kB() / 1024.0;

  *maxrss_mb = get_rss() / mb;
  *majflt_mb = *majflt / mb;
  *arena_mb = *arena / mb;
  *ordblks_mb = *ordblks / mb;
  *hblks_mb = *hblks / mb;
  *hblkhd_mb = *hblkhd / mb;
  *uordblks_mb = *uordblks / mb;
  *fordblks_mb = *fordblks / mb;
  *keepcost_mb = *keepcost / mb;
  *swap_used_mb = *swap_used / mb;

  *maxrss_kb = get_rss() / kb;
  *majflt_kb = *majflt / kb;
  *arena_kb = *arena / kb;
  *ordblks_kb = *ordblks / kb;
  *hblks_kb = *hblks / kb;
  *hblkhd_kb = *hblkhd / kb;
  *uordblks_kb = *uordblks / kb;
  *fordblks_kb = *fordblks / kb;
  *keepcost_kb = *keepcost / kb;
  *swap_used_kb = *swap_used / kb;

}
