/*@@
  @header    ioutil_CheckpointRecovery.h
  @date      Tue 19 Sep 2000
  @author    Thomas Radke
  @desc
             Typedefs and function prototypes for checkpointing and recovery.
  @history
  @endhistory
  @version $Header$
@@*/

#ifndef _IOUTIL_CHECKPOINTRECOVERY_H_
#define _IOUTIL_CHECKPOINTRECOVERY_H_

#ifdef __cplusplus
extern "C" {
#endif

/* enums for checkpointing/recovery and filereader functions */
enum {
  CP_INITIAL_DATA,
  CP_EVOLUTION_DATA,
  CP_RECOVER_PARAMETERS,
  CP_RECOVER_DATA,
  FILEREADER_DATA
};

/************************************************************************
 *
 *  Function prototypes
 *
 ************************************************************************/

/* create a checkpoint filename */
char *IOUtil_AssembleFilename(const cGH *GH, const char *basefilename,
                              const char *postfix, const char *extension,
                              int called_from, int file_ioproc,
                              int file_unchunked);

/* register a new recovery method */
int IOUtil_RegisterRecover(const char *name,
                           int (*func)(cGH *, const char *, int));

/* recover variables from a list of data files */
int IOUtil_RecoverVarsFromDatafiles(cGH *GH, const char *in_files,
                                    const char *in_vars);

/* generic recovery routine called by other IO thorns */
int IOUtil_RecoverParameters(int (*recoverFn)(cGH *GH, const char *basefilename,
                                              int called_from),
                             const char *fileExtension, const char *fileType);

/* return the parameter database as one single string */
char *IOUtil_GetAllParameters(const cGH *GH, int all);

/* set all parameters contained in the given string */
void IOUtil_SetAllParameters(const char *parameters);

/* print the checkpoint/recovery timings to stdout */
void IOUtil_PrintTimings(const char *description, int ntimers,
                         const int *timers,
                         const char *const *const timer_descriptions);

#ifdef __cplusplus
}
#endif

#endif /* _IOUTIL_CHECKPOINTRECOVERY_H_ */
