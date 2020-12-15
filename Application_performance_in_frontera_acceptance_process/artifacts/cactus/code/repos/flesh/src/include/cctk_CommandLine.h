 /*@@
   @header    cctk_CommandLine.h
   @date      Thu Jan 20 2000
   @author    Gabrielle Allen
   @desc 
   Prototypes and constants for command line functions
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTK_COMMANDLINE_H_
#define _CCTK_COMMANDLINE_H_

#ifdef __cplusplus 
extern "C" 
{
#endif

int CCTK_CommandLine(char ***outargv);
int CCTK_ParameterFilename(int len, char *filename);

#ifdef __cplusplus 
}
#endif

#endif /* _CCTK_COMMANDLINE_H_ */
