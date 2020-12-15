/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdio.h>
#include <string.h>
#include "converse.h"
#include <errno.h>

#if CMK_WHEN_PROCESSOR_IDLE_USLEEP
#include <sys/types.h>
#include <sys/time.h>
#endif

#if CMK_TIMER_USE_TIMES
#include <sys/times.h>
#include <limits.h>
#ifndef WIN32
#include <unistd.h>
#endif
#endif

#if CMK_TIMER_USE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#endif

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)

#ifdef __cplusplus
extern "C" {
#endif
void CApplicationInit(void);
void CApplicationDepositData(char *data);
void CApplicationDepositNode0Data(char *data);
#ifdef __cplusplus
}
#endif

#endif
