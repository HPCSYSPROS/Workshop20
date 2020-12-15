/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef DEBUG_H
#define DEBUG_H

#ifndef MIN_DEBUG_LEVEL
  #define MIN_DEBUG_LEVEL 0
#endif
#ifndef MAX_DEBUG_LEVEL
  #define MAX_DEBUG_LEVEL 10
#endif
#ifndef STDERR_LEVEL
  /* anything >= this error level goes to stderr */
  #define STDERR_LEVEL 5
#endif


/*****************************************************************
 *  DebugM(): function to display a debug message.
 *  Messages have different levels.  The low numbers are low severity
 *  while the high numbers are really important.  Very high numbers
 *  are sent to stderr rather than stdout.
 *  The default severity scale is from 0 to 10.
 *     0 = plain message
 *     4 = important message
 *     5 = warning (stderr)
 *     10 = CRASH BANG BOOM error (stderr)
 *  The remaining args are like printf: a format string and some args.
 *  This function can be turned off by compiling without the DEBUGM flag
 *  No parameters to this function should have a side effect!
 *  No functions should be passed as parameters!  (including inline)
 *****************************************************************/

 #ifdef DEBUGM

#include "InfoStream.h"

  #define Debug(x) (x)
  #define DebugM(level,format) \
	{ \
	  if ((level >= MIN_DEBUG_LEVEL) && (level <= MAX_DEBUG_LEVEL)) \
	  { \
	    infostream Dout; \
	    if (level >= STDERR_LEVEL)	Dout << "[ERROR " << level << "] "; \
	    else if (level > 0) Dout << "[Debug " << level << "] "; \
	    Dout << iPE << ' ' << iFILE; \
	    Dout << format << endi; \
	  } \
	}

 #else
  /* make a void function. */
  /* parameters with side effects will be removed! */
  #define Debug(x) ;
  #define DebugM(x,y)	;

 #endif /* DEBUGM */

#ifdef PROCTRACE_DEBUG
#include "charm++.h"
#include <stdarg.h>
#include "ProcessorPrivate.h"
class DebugFileTrace{
private:
    char *fname;
    FILE *fp;
public:
    inline static DebugFileTrace *Instance(char *fn){
        if(CkpvAccess(DebugFileTrace_instance)==0){
            CkpvAccess(DebugFileTrace_instance) = new DebugFileTrace(fn);
        }
        return CkpvAccess(DebugFileTrace_instance);
    }
    inline static DebugFileTrace *Object(){
        return CkpvAccess(DebugFileTrace_instance);
    }
    DebugFileTrace(char *fn){
        if(fn==NULL) {
            fname = NULL;
            fp = stdout;
            return;
        }else{
            char tmp[128];
            memset(tmp, 0, 128*sizeof(char));
            sprintf(tmp, "%s.%d", fn, CkMyPe());
            fname = new char[strlen(tmp)+1];
            memcpy(fname, tmp, strlen(tmp)+1);
            fp = fopen(fname, "w");
            fclose(fp);
        }
    }
    ~DebugFileTrace(){
        delete [] fname;
    }
    inline void writeTrace(const char *msg, ...){
        va_list argList;
        va_start(argList, msg);
        vfprintf(fp, msg, argList);
        va_end(argList);
    }   
    inline int openTrace(){ 
        if(fname==NULL)  return 0;
        fp = fopen(fname, "a"); 
        if(fp==NULL)
            return 1;
        else
            return 0;
    }
    inline int closeTrace(){ 
        if(fname==NULL)  return 0;
        return fclose(fp); 
    }
    inline int flushTrace(){
        return fflush(fp);
    }
};
#endif

#endif /* DEBUG_H */

