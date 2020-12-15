/*@@
   @file      CheckpointRecovery.c
   @date      Jun 04 1999
   @author    Thomas Radke
   @desc
              Utility routines for checkpointing/recovery and the filereader
              The actual work is done by the IO thorns.
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"
#include "util_Table.h"

#include "ioGH.h"
#include "ioutil_CheckpointRecovery.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#include "util_String.h"

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_IOUtil_CheckpointRecovery_c)

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
void IOUtil_RecoverGH(cGH *GH);
void IOUtil_RecoverIDFromDatafiles(cGH *GH);
void CCTK_FCALL CCTK_FNAME(IOUtil_RecoverVarsFromDatafiles)(int *result,
                                                            cGH **GH,
                                                            TWO_FORTSTRING_ARG);

/********************************************************************
 ********************    Static Variables   *************************
 ********************************************************************/
/* table to store registered recovery functions */
static int recover_fn_table = -1;
static int checkpoint_file_exists = 0;

/********************************************************************
 ********************    Internal Typedefs   ************************
 ********************************************************************/
typedef struct {
  char *basename;
  int iteration;
} filelist_t;

typedef int (*recover_fn_t)(cGH *GH, const char *basefilename, int called_from);

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static char *EncodeString(const char *string);
static char *DecodeString(const char *string);
static int OctalDigitToInt(const char digit);
static void SetInputFlag(int vindex, const char *optstring, void *arg);
#ifdef HAVE_DIRENT_H
static int CompareFiles(const void *a, const void *b);
#endif

/************************************************************************
 *
 *  Registration functions for Restoring from IO Files
 *
 ************************************************************************/

/*@@
  @routine    IOUtil_RegisterRecover
  @date       Monday 21 June 1999
  @author     Gabrielle Allen
  @desc
              Registers a new recovery method
  @enddesc
  @calls      Util_GetHandle
              Util_NewHandle

  @var        name
  @vdesc      The name of the function for recovery
  @vtype      const char *
  @vio        in
  @endvar
  @var        recover_fn
  @vdesc      the recovery function to register
  @vtype      recover_fn_t
  @vio        in
  @endvar

  @returntype int
  @returndesc
              0  for success, or<BR>
              -1 method with this name already registered
  @endreturndesc
@@*/
int IOUtil_RegisterRecover(const char *name, recover_fn_t recover_fn) {
  int retval;

  /* create the function pointer table on the first time through */
  if (recover_fn_table < 0) {
    recover_fn_table = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  }

  /* check that the method hasn't already been registered */
  retval =
      Util_TableQueryValueInfo(recover_fn_table, NULL, NULL, name) ? -1 : 0;
  if (!retval) {
    Util_TableSetFnPointer(recover_fn_table, (CCTK_FPOINTER)recover_fn, name);
  }

  return (retval);
}

/*@@
  @routine    IOUtil_AssembleFilename
  @date       30 June 2004
  @author     Thomas Radke
  @desc
              This routine assembles the full basename for checkpoint/recovery
              and filereader files, paying attention to the different types:

                * it returns the full filename (directory+filename)
                * for cp files it prepends the iteration number as "it_%d"
                * for chunked files it prepends the file number as "file_%d"

              This routine replaces the old routine IOUtil_PrepareFilename()
              which should not be used anymore.
  @enddesc

  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        basefilename
  @vdesc      basename of the file
  @vtype      const char *
  @vio        in
  @endvar
  @var        postfix
  @vdesc      optional postfix to append to the basename
              (eg. ".tmp" for temporary files, or "" for no postfix)
  @vtype      const char *
  @vio        in
  @var        extension
  @vdesc      file extension to append to the basefilename
              (should be given with the dot, eg. as ".h5")
  @vtype      const char *
  @vio        in
  @endvar
  @var        called_from
  @vdesc      indicates the caller function:
                * either Filereader (FILEREADER_DATA)
                * or IOUtil_RecoverGH() (CP_RECOVER_DATA)
  @vtype      int
  @vio        in
  @endvar
  @var        file_ioproc
  @vdesc      the IO processor number (for chunked files)
  @vtype      int
  @vio        in
  @endvar
  @var        file_unchunked
  @vdesc      flag to indicate whether file mode is unchunked or not
  @vtype      int
  @vio        in
  @endvar

  @returntype char *
  @returndesc
              the filename (must be freed after use)
  @endreturndesc
@@*/
char *IOUtil_AssembleFilename(const cGH *GH, const char *basefilename,
                              const char *postfix, const char *extension,
                              int called_from, int file_ioproc,
                              int file_unchunked) {
  DECLARE_CCTK_PARAMETERS
  size_t filenamelen;
  char *filename;
  const char *dir, *basename;
  char iteration_postfix[32], filenumber_postfix[32];

  /* get the right parameters */
  dir = basename = NULL;
  switch (called_from) {
  case CP_INITIAL_DATA:
    dir = checkpoint_dir;
    basename = checkpoint_ID_file;
    if (CCTK_CreateDirectory(0755, checkpoint_dir) < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Checkpoint directory %s could not be created",
                 checkpoint_dir);
    break;

  case CP_EVOLUTION_DATA:
    dir = checkpoint_dir;
    basename = checkpoint_file;
    if (CCTK_CreateDirectory(0755, checkpoint_dir) < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Checkpoint directory %s could not be created",
                 checkpoint_dir);
    break;

  case CP_RECOVER_DATA:
  case CP_RECOVER_PARAMETERS:
    dir = recover_dir;
    basename = basefilename ? basefilename : recover_file;
    break;
  case FILEREADER_DATA:
    dir = filereader_ID_dir;
    basename = basefilename ? basefilename : recover_file;
    break;

  default:
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "IOUtil_AssembleFilename: unknown calling mode %d", called_from);
    break;
  }

  /* if checkpoint filename, merge in the iteration number
     and for chunked files also the file number */
  iteration_postfix[0] = 0;
  if (called_from == CP_INITIAL_DATA || called_from == CP_EVOLUTION_DATA) {
    Util_snprintf(iteration_postfix, sizeof(iteration_postfix), ".it_%d",
                  (int)GH->cctk_iteration);
  }

  /* If not one unchunked file give a file number */
  filenumber_postfix[0] = 0;
  if (!file_unchunked) {
    Util_snprintf(filenumber_postfix, sizeof(filenumber_postfix), ".file_%d",
                  file_ioproc);
  }

  if (!(dir && basename && postfix && extension))
    CCTK_VWarn(0, __LINE__, __FILE__, "IOUtil",
               "Directory, basename, postfix or extension error");
  filenamelen = strlen(dir) + strlen(basename) + strlen(postfix) +
                strlen(iteration_postfix) + strlen(filenumber_postfix) +
                strlen(extension) + 2;
  filename = malloc(filenamelen);
  if (!(filename))
    CCTK_VWarn(0, __LINE__, __FILE__, "IOUtil", "File name error");
  Util_snprintf(filename, filenamelen, "%s/%s%s%s%s%s", dir, basename, postfix,
                iteration_postfix, filenumber_postfix, extension);

  return (filename);
}

/*@@
  @routine    IOUtil_RecoverFromFile
  @date       Jun 14 1999
  @author     Thomas Radke
  @desc
              Recover from a given file.
              This routine loops through all XXX_RecoverGH routines
              registered by IO thorns.
  @enddesc
  @calls      Util_GetHandledData
              <registered RecoverGH routines>

  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      cGH *
  @vio        in
  @endvar
  @var        basefilename
  @vdesc      basefilename of the file(s) to recover from
  @vtype      const char *
  @vio        in
  @endvar
  @var        called_from
  @vdesc      indicates the caller function:
                * either Filereader (FILEREADER_DATA)
                * or IOUtil_RecoverGH() (CP_RECOVER_DATA)
  @vtype      int
  @vio        in
  @endvar
@@*/
static int IOUtil_RecoverFromFile(cGH *GH, const char *basefilename,
                                  int called_from) {
  int keylen, ihandle, retval;
  char *key;
  CCTK_FPOINTER func;
  recover_fn_t recover_fn;

  if (Util_TableQueryNKeys(recover_fn_table) <= 0) {
    CCTK_WARN(1, "IOUtil_RecoverFromFile: No recovery routines "
                 "were registered");
    return (-1);
  }

  keylen = Util_TableQueryMaxKeyLength(recover_fn_table) + 1;
  key = malloc(keylen);

  for (ihandle = Util_TableItCreate(recover_fn_table), retval = -1;
       Util_TableItQueryIsNonNull(ihandle) > 0 && retval < 0;
       Util_TableItAdvance(ihandle)) {
    Util_TableItQueryKeyValueInfo(ihandle, keylen, key, NULL, NULL);
    Util_TableGetFnPointer(recover_fn_table, &func, key);
    if (func) {
      recover_fn = (recover_fn_t)func;
      retval = recover_fn(GH, basefilename, called_from);
    }
  }
  free(key);
  Util_TableItDestroy(ihandle);

  return (retval);
}

/*@@
  @routine    IOUtil_RecoverGH
  @date       Jun 14 1999
  @author     Thomas Radke
  @desc
              The rfr-registered recovery routine.
              Just calls IOUtil_RecoverFromFile()
              with called_from == CP_RECOVER_DATA.
  @enddesc
  @calls      IOUtil_RecoverFromFile

  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      cGH *
  @vio        in
  @endvar
@@*/
void IOUtil_RecoverGH(cGH *GH) {
  if (checkpoint_file_exists) {
    /* stop if recovery failed */
    if (IOUtil_RecoverFromFile(GH, NULL, CP_RECOVER_DATA) < 0) {
      CCTK_WARN(0, "Failed to restart from recovery !");
    }
  }
}

/*@@
  @routine    IOUtil_RecoverVarsFromDatafiles
  @date       Wed Apr 19 2000
  @author     Thomas Radke
  @desc
              IOUtil's function interface to recover variables from data files.
              Any thorn can call this routine with a list of datafile names
              and a list of variables to read in from these files. <BR>

              It just calls IOUtil_RecoverFromFile() with
              called_from == FILEREADER_DATA for each data file
              from the given file list.
  @enddesc

  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      cGH *
  @vio        in
  @endvar
  @var        in_files
  @vdesc      list of filenames to process (separated by spaces)
  @vtype      const char *
  @vio        in
  @endvar
  @var        in_vars
  @vdesc      list of variable fullnames to read in from the files
  @vtype      const char *
  @vio        in
  @endvar

  @returntype int
  @returndesc
              total number of recovered variables
  @endreturndesc
@@*/
int IOUtil_RecoverVarsFromDatafiles(cGH *GH, const char *in_files,
                                    const char *in_vars) {
  int retval, num_recovered_vars;
  ioGH *myGH;
  char *basefilename, *delim, delim_char;
  DECLARE_CCTK_PARAMETERS

  myGH = CCTK_GHExtension(GH, "IO");
  if (myGH == NULL) {
    CCTK_WARN(0, "invalid cctkGH pointer passed into "
                 "IOUtil_RecoverVarsFromDatafiles()");
  }

  if (CCTK_NumVars() > 0) {
    void *calldata[2];
    myGH->do_inVars = calloc(CCTK_NumVars(), sizeof(CCTK_INT));
    myGH->alias = calloc(CCTK_NumVars(), sizeof(char *));
    calldata[0] = myGH->do_inVars;
    calldata[1] = myGH->alias;
    if (CCTK_TraverseString(in_vars, SetInputFlag, calldata,
                            CCTK_GROUP_OR_VAR) < 0) {
      CCTK_WARN(myGH->stop_on_parse_errors ? 0 : 1,
                "error while parsing parameter 'IO::filereader_ID_vars'");
    }
  } else {
    myGH->do_inVars = NULL;
    myGH->alias = NULL;
  }

  num_recovered_vars = 0;

  /* duplicate the filename list and parse it */
  basefilename = strdup(in_files);
  while (basefilename && *basefilename) {
    /* skip leading spaces */
    while (isspace((int)*basefilename)) {
      basefilename++;
    }
    if (!*basefilename) {
      break;
    }

    /* find delimiter for current filename and cut there */
    for (delim = basefilename + 1; !isspace((int)*delim) && *delim; delim++)
      ;
    delim_char = *delim;
    *delim = 0;

    if (!CCTK_Equals(verbose, "none")) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Reading variables from file with base name '%s'",
                 basefilename);
    }

    retval = IOUtil_RecoverFromFile(GH, basefilename, FILEREADER_DATA);
    if (retval >= 0) {
      num_recovered_vars += retval;
    } else {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Failed to read variables from data file with base name '%s'",
                 basefilename);
    }

    *delim = delim_char;
    basefilename = delim;
  }

  if (basefilename) {
    free(basefilename - strlen(in_files));
  }

  if (myGH->do_inVars) {
    free(myGH->do_inVars);
    myGH->do_inVars = NULL;
  }

  if (myGH->alias) {
    const int numVars = CCTK_NumVars();
    for (int i = 0; i < numVars; i++) {
      if (myGH->alias[i])
        free((void *)myGH->alias[i]);
    }
    free(myGH->alias);
    myGH->alias = NULL;
  }

  return (num_recovered_vars);
}

void CCTK_FCALL
    CCTK_FNAME(IOUtil_RecoverVarsFromDatafiles)(int *result, cGH **GH,
                                                TWO_FORTSTRING_ARG) {
  TWO_FORTSTRING_CREATE(in_files, in_vars)
  *result = IOUtil_RecoverVarsFromDatafiles(*GH, in_files, in_vars);
  free(in_files);
  free(in_vars);
}

/*@@
  @routine    IOUtil_RecoverIDFromDatafiles
  @date       Wed Apr 19 2000
  @author     Thomas Radke
  @desc
              The rfr-registered initial data recovery routine.
              Just calls IOUtil's generic routine
              IOUtil_RecoverVarsFromDatafiles() with the filereader_ID_XXX
              parameters.
  @enddesc
  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      cGH *
  @vio        in
  @endvar
@@*/
void IOUtil_RecoverIDFromDatafiles(cGH *GH) {
  ioGH *myGH;
  DECLARE_CCTK_PARAMETERS

  myGH = CCTK_GHExtension(GH, "IO");
  myGH->stop_on_parse_errors = strict_io_parameter_check;
  IOUtil_RecoverVarsFromDatafiles(GH, filereader_ID_files, filereader_ID_vars);
  myGH->stop_on_parse_errors = 0;
}

/*@@
  @routine    IOUtil_RecoverParameters
  @date       Apr 22 2000
  @author     Thomas Radke
  @desc
              The generic parameter recovery routine.
              It is called by the IO thorns' parameter recovery routines
              scheduled at CCTK_RECOVER_PARAMETERS, and simply calls
              the given callback routine with its arguments
              plus a checkpoint filename.
  @enddesc

  @var        recover_fn
  @vdesc      callback function for recovery of parameters
              from a given checkpoint file
  @vtype      int (*) (cGH *, const char *, int)
  @vio        in
  @endvar
  @var        fileExtension
  @vdesc      extension of valid checkpoint files for given callback
  @vtype      const char *
  @vio        in
  @endvar
  @var        fileType
  @vdesc      string to describe the type of checkpoint file
              (used for warning/info messages)
  @vtype      const char *
  @vio        in
  @endvar

  @returntype int
  @returndesc
               0 if in "autoprobe" mode and no cp files were found, or<BR>
              +1 if parameter recovery was successful for some cp file,<BR>
              -1 if in "auto" mode and no checkpoint files were found,
                 or if parameter recovery failed for some cp file,<BR>
              -2 if in "auto*" mode and recovery dir doesn't exist
  @endreturndesc
@@*/
int IOUtil_RecoverParameters(int (*recover_fn)(cGH *GH,
                                               const char *basefilename,
                                               int called_from),
                             const char *fileExtension, const char *fileType) {
  int retval;
#ifdef HAVE_DIRENT_H
  int len, recover_file_len;
  unsigned int num_files;
  const char *p, *q;
  DIR *dir;
  struct dirent *file;
  filelist_t *filelist, *tmp;
#endif
  DECLARE_CCTK_PARAMETERS
#ifdef HAVE_DIRENT_H
  bool verbose_full = CCTK_Equals(verbose, "full");
#endif

  if (CCTK_Equals(recover, "auto") || CCTK_Equals(recover, "autoprobe")) {
#ifdef HAVE_DIRENT_H
    if (verbose_full) {
      CCTK_VInfo(CCTK_THORNSTRING, "Searching for %s checkpoint files "
                                   "with basefilename '%s' in directory '%s'",
                 fileType, recover_file, recover_dir);
    }

    dir = opendir(recover_dir);
    if (!dir) {
      /* The recovery directory does not exist */
      if (CCTK_Equals(recover, "autoprobe")) {
        /* This is harmless when "autoprobe" is used */
        CCTK_VInfo(CCTK_THORNSTRING, "Recovery directory '%s' doesn't exist",
                   recover_dir);
        return 0;
      } else {
        /* This is an error when "auto" is used */
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Recovery directory '%s' doesn't exist", recover_dir);
        return -2;
      }
    }

    /* get the list of potential recovery files */
    recover_file_len = strlen(recover_file);
    num_files = 0;
    filelist = NULL;

    while ((file = readdir(dir)) != NULL) {
      /* first check the file prefix */
      if (strncmp(file->d_name, recover_file, recover_file_len) ||
          strncmp(file->d_name + recover_file_len, ".it_", 4)) {
        continue;
      }

      /* now check if there is an iteration number following the file prefix */
      for (p = file->d_name + recover_file_len + 4; *p && *p != '.'; p++) {
        if (!isdigit((int)*p)) {
          break;
        }
      }

      /* check for a '.file_<processor>' suffix for chunked output */
      q = p;
      if (!strncmp(q, ".file_", 6)) {
        /* read past the <processor> number */
        for (q = q + 6; *q && *q != '.'; q++) {
          if (!isdigit((int)*q)) {
            break;
          }
        }
      }

      /* finally check the file extension suffix */
      if (*q != '.' || strcmp(q, fileExtension)) {
        continue;
      }

      /* found a recovery file by that basename */
      if (num_files == 0) {
        tmp = malloc(sizeof(filelist_t));
      } else {
        tmp = realloc(filelist, (num_files + 1) * sizeof(filelist_t));
      }
      if (tmp == NULL) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Failed to allocate memory for file list");
        continue;
      }
      filelist = tmp;
      filelist[num_files].basename = strdup(file->d_name);
      /* cut the filename after the iteration number field */
      len = p - file->d_name;
      filelist[num_files].basename[len] = 0;
      filelist[num_files].iteration = atoi(file->d_name + recover_file_len + 4);

      num_files++;
    }
    closedir(dir);

    retval = CCTK_Equals(recover, "auto") ? -1 : 0;
    if (num_files) {
      /* sort the list according to their iteration numbers */
      qsort(filelist, num_files, sizeof(filelist_t), CompareFiles);

      /* loop over all recovery files found and call the callback routine;
         skip all following files after the first successful recovery (when
         recover_fn() returned a positive value) */
      while (num_files--) {
        if (retval <= 0) {
          retval = recover_fn(NULL, filelist[num_files].basename,
                              CP_RECOVER_PARAMETERS);
        }
        free(filelist[num_files].basename);
      }
      free(filelist);
    } else {
      CCTK_VWarn(retval ? 1 : 3, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "No %s checkpoint files with basefilename '%s' found in "
                 "recovery directory '%s'",
                 fileType, recover_file, recover_dir);
    }
#else
    fileExtension = fileExtension;

    /* No opendir(3) ??? It's probably a Windows box, so just give up ! */
    CCTK_WARN(0, "You cannot use 'IO::recover = \"auto*\"' on "
                 "this architecture because it doesn't provide opendir(3) to "
                 "automatically look for checkpoint files.\n"
                 "Please use 'IO::recover = \"manual\"' instead !");
    retval = -1;
#endif
  } else {
    /* just call the recovery routine */
    retval = (*recover_fn)(NULL, recover_file, CP_RECOVER_PARAMETERS);
  }

  if (retval < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not recover parameters from %s checkpoint file(s) "
               "with basefilename '%s' in recovery directory '%s'",
               fileType, recover_file, recover_dir);
  }

  /* remember parameter recovery status for later evaluation in
     IOUtil_RecoverGH() */
  checkpoint_file_exists = retval > 0;

  return (retval);
}

/*@@
  @routine    IOUtil_GetAllParameters
  @date       Mon Apr 10 2000
  @author     Thomas Radke
  @desc
              Collect all parameters of active implementations
              into a single string which can then be dumped as an attribute.

              Each "<parameter> = <value>" pair in the resulting string
              is separated by a newline character from the next pair.
              The "<value>" substring is encoded such that non-printable
              characters are contained as escape sequences, and a single
              '\' character is escaped by another '\' character.
  @enddesc

  @var        GH
  @vdesc      pointer to CCTK grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        all
  @vdesc      flag indicating whether to save all parameters or just the ones
              which have been set before
  @vtype      int
  @vio        in
  @endvar

  @returntype char *
  @returndesc the allocated string or NULL
  @endreturndesc
@@*/
char *IOUtil_GetAllParameters(const cGH *GH, int all) {
  int i, first, add_len, current_len, max_len;
  char *param, *value, *encoded_value, *tmp, *retval;
  const char *thorn;
  const cParamData *pdata;
  DECLARE_CCTK_PARAMETERS

  /* avoid compiler warning about unused parameter */
  GH = GH;

  /* preallocate buffer for parameters, to prevent too many subsequent reallocs
     (causing unnecessary memory fragmentation) */
  max_len = 1024 * 1024;
  retval = malloc(max_len);
  if (!retval) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to allocate buffer of %d bytes", max_len);
  }
  current_len = 0;

  /* loop over all thorns */
  for (i = CCTK_NumCompiledThorns() - 1; i >= 0; i--) {
    thorn = CCTK_CompiledThorn(i);

    /* skip all inactive thorns */
    if (!CCTK_IsThornActive(thorn)) {
      continue;
    }

    /* now walk through all parameters of given thorn */
    first = 1;
    while (CCTK_ParameterWalk(first, thorn, &param, &pdata) == 0) {
      first = 0;

      /* if 'all' is false get only parameters which have been set */
      if (all || pdata->n_set > 0) {
        value = CCTK_ParameterValString(pdata->name, pdata->thorn);
        if (value == NULL) {
          CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Couldn't get value for parameter '%s'", param);
        } else {
          encoded_value = EncodeString(value);
          add_len = strlen(param) + strlen(encoded_value) + 5;
          if (current_len + add_len >= max_len) {
            /* double new buffer length until buffer is large enough */
            while (current_len + add_len >= max_len) {
              max_len *= 2;
            }
            tmp = realloc(retval, max_len);
            if (!tmp) {
              CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Failed to allocate buffer of %d bytes", max_len);
            }
            retval = tmp;
          }

          sprintf(retval + current_len, "%s = %s\n", param, encoded_value);
          current_len += add_len - 1;

          free(encoded_value);
          free(value);
        }
      }

      free(param);

    } /* end of loop walking over all parameters of given thorn */
  }   /* end of looping over all thorns */

  return (retval);
}

/*@@
  @routine    IOUtil_SetAllParameters
  @date       Mon Apr 10 2000
  @author     Thomas Radke
  @desc
              Parse the given string for parameters
              and call CCTK_SetParameter() for each.

              The input string is expected to follow the format of
              IOUtil_GetAllParameters() (see above) where each
              "<parameter> = <value>" pair is separated by a newline character
              from the next pair.
              The "<value>" substring is encoded such that non-printable
              characters are contained as escape sequences, and a single
              '\' character is escaped by another '\' character. This
              encoding is undone before setting the a value of a parameter.
  @enddesc

  @var        parameters
  @vdesc      the parameter string
  @vtype      const char *
  @vio        in
  @endvar
@@*/
void IOUtil_SetAllParameters(const char *parameters) {
  char *tmp, *nextparam, *avalue, *decoded_value, *param;
  char oldchar;
  char *name, *thorn_impl, *parameter_string, *free_me;
  const char *thorn;
  const cParamData *paramdata;
  int ierr, steered_non_steerable;
  DECLARE_CCTK_PARAMETERS
  bool verbose_full = CCTK_Equals(verbose, "full");

  parameter_string = free_me = strdup(parameters);
  steered_non_steerable = 0;
  while (*parameter_string) {
    /* find the end of the current "<parameter> = <value>" pair substring */
    nextparam = parameter_string;
    while (*nextparam != '\n' && *nextparam) {
      nextparam++;
    }
    oldchar = *nextparam;
    *nextparam = 0;

    /* find the end of "<parameter>" */
    tmp = parameter_string;
    while (*tmp != ' ') {
      tmp++;
    }
    *tmp = 0;

    param = parameter_string;
    /* skip the " = " delimiter to get to "<value>" */
    avalue = tmp + 3;

    name = thorn_impl = NULL;
    if (Util_SplitString(&thorn_impl, &name, param, "::") == 0) {
      /* find out the implementing thorn of the parameter given */
      thorn = CCTK_ActivatingThorn(thorn_impl);
      if (!thorn) {
        thorn = thorn_impl;
      }

      /* set parameter only if it belongs to an active implementation
         and is not an accumulator parameter */
      if (CCTK_IsThornActive(thorn)) {
        paramdata = CCTK_ParameterData(name, thorn);
        if (paramdata && !paramdata->accumulator_expression) {
          decoded_value = DecodeString(avalue);
          ierr = CCTK_ParameterSet(name, thorn, decoded_value);
          if (ierr == -10) {
            steered_non_steerable = 1;
          } else if (ierr < 0) {
            CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Couldn't set parameter '%s' to '%s'", param,
                       decoded_value);
          }
          free(decoded_value);
        }
      } else if (verbose_full) {
        CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Ignoring inactive parameter '%s' for recovery", param);
      }

      if (name) {
        free(name);
      }
      if (thorn_impl) {
        free(thorn_impl);
      }
    }

    *nextparam = oldchar;
    parameter_string = nextparam;
    if (*parameter_string == '\n') {
      parameter_string++;
    }
  }

  free(free_me);

  if (steered_non_steerable) {
    CCTK_ERROR("Attempted to steer non steerable parameters.");
  }
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/* encodes a C string with all the non-printable characters converted
   into a 4-char-long escape sequence representation "\xxx" and
   the single '\' character converted into "\\" */
static char *EncodeString(const char *string) {
  size_t i, len = strlen(string);
  char *p, *retval;

  p = retval = malloc(4 * len + 1);
  if (!retval) {
    CCTK_WARN(0, "Out of memory !");
  }

  for (i = 0; i < len; i++) {
    if (!isprint(string[i])) {
      sprintf(p, "\\%.03o", (unsigned char)string[i]);
      p += 4;
    } else {
      *p++ = string[i];
      if (string[i] == '\\') {
        *p++ = string[i];
      }
    }
  }
  *p = 0;

  return (retval);
}

/* decodes a string encoded by EncodeString() back into a C string */
static char *DecodeString(const char *string) {
  size_t i, len = strlen(string);
  char *p, *retval;
  int o1, o2, o3;

  p = retval = malloc(len + 1);
  if (!retval) {
    CCTK_WARN(0, "Out of memory !");
  }

  for (i = 0; i < len; i++) {
    if (string[i] != '\\') {
      *p++ = string[i];
    } else {
      if (i + 1 >= len)
        goto bail_out;
      if (string[i + 1] == '\\') {
        *p++ = '\\';
        i++;
      } else {
        if (i + 3 >= len)
          goto bail_out;
        o1 = OctalDigitToInt(string[i + 1]);
        o2 = OctalDigitToInt(string[i + 2]);
        o3 = OctalDigitToInt(string[i + 3]);
        if (o1 < 0 || o2 < 0 || o3 < 0)
          goto bail_out;
        *p++ = (char)(8 * 8 * o1 + 8 * o2 + o3);
        i += 3;
      }
    }
  }
  *p = 0;

  return (retval);

bail_out:
  free(retval);
  CCTK_WARN(1, "Malformed string");
  return NULL;
}

/* callback for CCTK_TraverseString() to set the input flag
   for the given variable */
static void SetInputFlag(int vindex, const char *optstring, void *calldata) {
  int table, iterator;
  char key[128];
  CCTK_INT type, nelems;
  void **flags = (void **)calldata;
  CCTK_INT *do_inVars = (CCTK_INT *)flags[0];
  const char **alias = (const char **)flags[1];

  /* default -1 is to read the last iteration number from the file */
  do_inVars[vindex] = -1;

  /* default is NULL to use CCTK_FullName() */
  alias[vindex] = NULL;

  if (optstring) {
    table = Util_TableCreateFromString(optstring);
    if (table >= 0) {
      if (Util_TableQueryValueInfo(table, &type, &nelems, "cctk_iteration") >
          0) {
        if (type == CCTK_VARIABLE_INT && nelems == 1) {
          Util_TableGetInt(table, &do_inVars[vindex], "cctk_iteration");

          /* if a specific iteration number was given then increment it
             to keep 0 as disabling value */
          if (do_inVars[vindex] >= 0) {
            do_inVars[vindex]++;
          }
        } else {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Invalid value for option 'cctk_iteration' in option "
                     "string '%s' (must be an integer)",
                     optstring);
          CCTK_WARN(1, "Option will be ignored by file reader routines");
        }
        Util_TableDeleteKey(table, "cctk_iteration");
      }

      if (Util_TableQueryValueInfo(table, &type, &nelems, "alias") > 0) {
        if (type == CCTK_VARIABLE_CHAR) {
          char *buf = malloc(nelems + 1);
          if (!(buf)) {
            CCTK_VError(__LINE__, __FILE__, "IOUtil",
                        "Could not allocate %d bytes", (int)nelems + 1);
          }
          int len = Util_TableGetString(table, nelems + 1, buf, "alias");
          if (len != nelems) {
            CCTK_VError(__LINE__, __FILE__, "IOUtil",
                        "Unexpected length %d of 'alias' option value in '%s'. "
                        "Expected length %d.",
                        len, optstring, (int)nelems);
          }

          alias[vindex] = buf;
        } else {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Invalid value for option 'alias' in option "
                     "string '%s' (must be a string)",
                     optstring);
          CCTK_WARN(1, "Option will be ignored by file reader routines");
        }
        Util_TableDeleteKey(table, "alias");
      }

      /* warn about other options */
      for (iterator = Util_TableItCreate(table);
           Util_TableItQueryIsNonNull(iterator) > 0 &&
               Util_TableItQueryKeyValueInfo(iterator, sizeof(key), key, 0, 0) >
                   0;
           Util_TableItAdvance(iterator)) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Found option with unrecognized key '%s' in option string "
                   "'%s'",
                   key, optstring);
        CCTK_WARN(1, "Option will be ignored by file reader routines");
      }
      Util_TableItDestroy(iterator);

      Util_TableDestroy(table);
    } else {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Couldn't parse option string '%s'", optstring);
      CCTK_WARN(1, "Option will be ignored by file reader routines");
    }
  }
}

/* Convert an octal digit to an integer.  Return -1 if this is not an
   octal digit.  Note that this routine does not assume ASCII encoding
   (and thus cannot calculate digit-'0'), and does not use
   non-standard library function such as isoctal.  */
static int OctalDigitToInt(const char digit) {
  switch (digit) {
  case '0':
    return 0;
  case '1':
    return 1;
  case '2':
    return 2;
  case '3':
    return 3;
  case '4':
    return 4;
  case '5':
    return 5;
  case '6':
    return 6;
  case '7':
    return 7;
  }
  return -1;
}

#ifdef HAVE_DIRENT_H
/* callback for qsort() to sort the list of recovery files found */
static int CompareFiles(const void *a, const void *b) {
  return (((const filelist_t *)a)->iteration -
          ((const filelist_t *)b)->iteration);
}
#endif
