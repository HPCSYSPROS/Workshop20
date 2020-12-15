 /*@@
   @header    CommandLine.h
   @date      Wed Feb 17 00:53:55 1999
   @author    Tom Goodale
   @desc 
   Header for the warning functions.
   @enddesc 
 @@*/

#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_

#ifdef __cplusplus 
extern "C" {
#endif

typedef enum {REDIRECT_NONE, REDIRECT_NONROOT, REDIRECT_ALL} redirect_t;

void CCTKi_CommandLineTestThornCompiled(const char *argument);
void CCTKi_CommandLineDescribeAllParameters(const char *argument);
void CCTKi_CommandLineDescribeParameter(const char *argument);
void CCTKi_CommandLineTestParameters(const char *argument);
void CCTKi_CommandLineLoggingLevel(const char *argument);
void CCTKi_CommandLineWarningLevel(const char *argument);
void CCTKi_CommandLineErrorLevel(const char *argument);
void CCTKi_CommandLineParameterLevel(const char *argument);
void CCTKi_CommandLineRedirect(const char *argument, const redirect_t type);
void CCTKi_CommandLineLogDir(const char *argument);
void CCTKi_CommandLineSetBuffering(const char *argument);
void CCTKi_CommandLinePrintSchedule(void);
void CCTKi_CommandLineListThorns(void);
void CCTKi_CommandLineVersion(void);
void CCTKi_CommandLineHelp(void);
void CCTKi_CommandLineUsage(void);
void CCTKi_CommandLineFinished(void);

int CCTKi_ExitAfterParamCheck(void);

#ifdef __cplusplus 
}
#endif

#endif
