/*@@
  @file      Startup.c
  @date      Sat Feb 6 1999
  @author    Gabrielle Allen
  @desc
             Startup routines for IOUtil.
  @enddesc
  @version   $Id$
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Misc.h"
#include "cctk_Version.h"
#include "cctk_Parameters.h"
#include "util_Network.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#include "ioGH.h"
#include "ioutil_Utils.h"
#include "ioutil_AdvertisedFiles.h"

/* the rcs ID and its dummy funtion to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_IOUtil_Startup_c)

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
int IOUtil_Startup(void);
void IOUtil_UpdateParFile(CCTK_ARGUMENTS);

/********************************************************************
 *********************     Internal Typedefs   **********************
 ********************************************************************/
/* structure defining a single-linked list of parameters
   (stores the name of the parameter, its stringified value, and its type) */
typedef struct t_param_list {
  char *name;
  char *value;
  const char *format_string;
  struct t_param_list *next;
} t_param_list;

/********************************************************************
 **********************     Static Routines     *********************
 ********************************************************************/
static void *SetupGH(tFleshConfig *config, int convergence_level, cGH *GH);
static int CopyParFile(int recovered);
static int GenerateParFile(int recovered);
static int DumpParameters(FILE *outfile);

/*@@
  @routine    IOUtil_Startup
  @date       Sat Feb 6 1999
  @author     Gabrielle Allen
  @desc
              The startup registration routine for IOUtil.<BR>
              It registers the GH extension "IO" for IOUtil, along with its
              setup routine.
              It also advertises the original parameter file.
  @enddesc
  @calls      CCTK_RegisterGHExtension
              CCTK_RegisterGHExtensionSetupGH
              CCTK_ParameterFilename
              IOUtil_AdvertiseFile
@@*/
int IOUtil_Startup(void) {
  char parfile[1024];
  ioAdvertisedFileDesc advertised_file;
  DECLARE_CCTK_PARAMETERS

  CCTK_RegisterGHExtensionSetupGH(CCTK_RegisterGHExtension("IO"), SetupGH);

  /* advertise the parameter file */
  parfile[0] = 0;

  /* TODO: check the returned size; the function silently truncates
   * the filename if it doesn't fit in the buffer; I can't think of a
   * situation where this is useful. Even better, modify
   * CCTK_ParameterFilename to allocate memory for the string, like
   * many other Cactus functions do. */

  CCTK_ParameterFilename(sizeof(parfile), parfile);
  advertised_file.slice = "";
  advertised_file.thorn = CCTK_THORNSTRING;
  advertised_file.varname = "";
  advertised_file.description = "Parameter File";
  advertised_file.mimetype = "text/plain";
  IOUtil_AdvertiseFile(NULL, parfile, &advertised_file);

  return 0;
}

/*@@
  @routine    IOUtil_UpdateParFile
  @date       Tue 18 Dec 2001
  @author     Thomas Radke
  @desc
              Updates the parameter file for all parameters which got steered
              during its last invocation.
  @enddesc

  @var        GH
  @vdesc      pointer to grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar

  @returntype int
  @returndesc 0 for success
  @endreturndesc
@@*/
void IOUtil_UpdateParFile(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  const cGH *GH = cctkGH;

  /* check if it's time to update */
  if (CCTK_MyProc(GH) > 0 || GH->cctk_iteration % parfile_update_every) {
    return;
  }

  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
             "Not yet implemented: updating parameter file for steered "
             "parameters up to iteration %d",
             GH->cctk_iteration);
}

/*@@
  @routine    IOUtil_TruncateOutputFiles
  @date       Wed 17 November 2004
  @author     Thomas Radke
  @desc
              Aliased function 'IO_TruncateOutputFiles' which can be called
              by I/O methods to check whether to truncate existing output
              files or not.
  @enddesc

  @var        GH
  @vdesc      pointer to grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar

  @returntype CCTK_INT
  @returndesc 1 or 0 for truncate or not
  @endreturndesc
@@*/
CCTK_INT IOUtil_TruncateOutputFiles(const CCTK_POINTER_TO_CONST GH) {
  const ioGH *myGH;
  DECLARE_CCTK_PARAMETERS

  myGH = CCTK_GHExtension(GH, "IO");
  if (!myGH) {
    CCTK_WARN(CCTK_WARN_ABORT, "IOUtil_TruncateOutputFiles called before being "
                               "fully initialized. Please make sure to only "
                               "call me after the STARTUP bin has finished and "
                               "the driver called CCTKi_InitGHExtensions.");
  }

  return (myGH->recovered ? truncate_files_after_recovering : truncate_files);
}

/********************************************************************
 ***********************    Local Functions   ***********************
 ********************************************************************/
/*@@
  @routine    SetupGH
  @date       Tue May 09 2000
  @author     Thomas Radke
  @desc
              The GH allocation and initialization routine for IOUtil.<BR>
              Necessary output dirs are created, checkpoint/recovery timers
              are created if timing information is wanted, and a parameter
              file is written to 'IO::out_dir' if requested.
  @enddesc
  @calls      IOUtil_CreateDirectory
              CopyParFile
              GenerateParFile

  @var        config
  @vdesc      flesh configuration structure (unused)
  @vtype      tFleshConfig *
  @vio        in
  @endvar
  @var        convergence_level
  @vdesc      convergence level (unused)
  @vtype      int
  @vio        in
  @endvar
  @var        GH
  @vdesc      pointer to grid hierarchy
  @vtype      cGH *
  @vio        in
  @endvar

  @returntype void *
  @returndesc the pointer to IOUtil's GH extension structure, or<BR>
              NULL if memory allocation failed
  @endreturndesc
@@*/
static void *SetupGH(tFleshConfig *config, int convergence_level, cGH *GH) {
  int i, maxdim, myproc;
  ioGH *myGH;
#ifdef HAVE_DIRENT_H
  DIR *dir;
  struct dirent *file;
#endif
  DECLARE_CCTK_PARAMETERS

  /* avoid compiler warnings about unused parameters */
  convergence_level = convergence_level;

  myproc = CCTK_MyProc(GH);
  myGH = calloc(1, sizeof(ioGH));
  if (!myGH) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't allocate GH extension structure for IOUtil");
    return (NULL);
  }

  if (CCTK_Equals(out_mode, "proc")) {
    myGH->ioproc = myproc;
    myGH->nioprocs = CCTK_nProcs(GH);
    myGH->ioproc_every = 1;
  } else if (CCTK_Equals(out_mode, "np")) {
    myGH->ioproc_every = out_proc_every;
    if (out_proc_every > CCTK_nProcs(GH)) {
      myGH->ioproc_every = CCTK_nProcs(GH);
      CCTK_VInfo(CCTK_THORNSTRING, "Reducing 'IO::out_proc_every' to %d",
                 myGH->ioproc_every);
    }

    myGH->nioprocs = CCTK_nProcs(GH) / myGH->ioproc_every +
                     (CCTK_nProcs(GH) % myGH->ioproc_every ? 1 : 0);
    myGH->ioproc = myproc - (myproc % myGH->ioproc_every);
  } else /* IO::out_mode = "onefile" */
  {
    myGH->ioproc = 0;
    myGH->nioprocs = 1;
    myGH->ioproc_every = CCTK_nProcs(GH);
  }

  /* For now we can only have unchunked for a single output file */
  myGH->unchunked = 0;
  if (out_unchunked || CCTK_nProcs(GH) == 1) {
    if (myGH->ioproc_every >= CCTK_nProcs(GH)) {
      myGH->unchunked = 1;
    } else {
      CCTK_INFO("Unchunked output not supported for multiple "
                "output files. Output will be chunked.");
    }
  }

  myGH->stop_on_parse_errors = 0;

  /* create the default output and checkpoint directories */
  i = IOUtil_CreateDirectory(GH, out_dir, !CCTK_Equals(out_mode, "onefile"),
                             myGH->ioproc);
  if (i < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Problem creating default output directory '%s'", out_dir);
  } else if (i > 0) {
    if (CCTK_Equals(verbose, "full")) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "default output directory '%s' already exists", out_dir);
    }
    if (require_empty_output_directory) {
#ifdef HAVE_DIRENT_H
      i = 0;
      dir = opendir(out_dir);
      while ((file = readdir(dir)) != NULL) {
        if (strcmp(file->d_name, ".") == 0 || strcmp(file->d_name, "..") == 0) {
          continue;
        }
        CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Found file '%s' in already existing output directory '%s'",
                   file->d_name, out_dir);
        i++;
      }
      closedir(dir);
      if (i) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "You have set IO::out_dir = '%s' and "
                   "IO::require_empty_output_directory = 'yes'. "
                   "This output directory already exists and is non-empty.",
                   out_dir);
      }
#else
      /* No opendir(3) ??? It's probably a Windows box, so just give up ! */
      CCTK_WARN(0, "You cannot use 'IO::require_empty_output_directory = "
                   "\"yes\"' on this architecture because it doesn't provide "
                   "opendir(3) to browse the IO::out_dir  directory.\n"
                   "Please use 'IO::require_empty_output_directories = "
                   "\"no\"' instead !");
#endif
    }
  }
  i = IOUtil_CreateDirectory(GH, checkpoint_dir,
                             !CCTK_Equals(out_mode, "onefile"), myGH->ioproc);
  if (i < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Problem creating checkpoint directory '%s'", checkpoint_dir);
  } else if (i > 0 && CCTK_Equals(verbose, "full")) {
    CCTK_VInfo(CCTK_THORNSTRING, "checkpoint directory '%s' already exists",
               checkpoint_dir);
  }

  /* save downsampling parameters in ioUtilGH because they are temporarily
     reset during checkpointing */
  /* for now we have only parameters for the first 3 dimensions
     the rest is constantly initialized to 1 */
  maxdim = CCTK_MaxDim();
  myGH->downsample = malloc(maxdim * sizeof(int));

  switch (maxdim > 3 ? 3 : maxdim) {
  case 3:
    myGH->downsample[2] = out_downsample_z;
  case 2:
    myGH->downsample[1] = out_downsample_y;
  case 1:
    myGH->downsample[0] = out_downsample_x;
  }

  for (i = 3; i < maxdim; i++) {
    myGH->downsample[i] = 1;
  }

/* evaluate the 'IO::out_single_precision' parameter only
   if Cactus was compiled with double precision */
#ifdef SINGLE_PRECISION
  myGH->out_single = 0;
#else
  myGH->out_single = out_single_precision;
#endif

  /* copy the 'recovered' flag to this GH extension */
  myGH->recovered = config->recovered;

  /* reset the flags array for the file reader */
  myGH->do_inVars = NULL;

  /* write the parameter file if requested */
  if (myproc == 0) {
    if (CCTK_Equals(parfile_write, "copy")) {
      CopyParFile(myGH->recovered);
    } else if (CCTK_Equals(parfile_write, "generate")) {
      GenerateParFile(myGH->recovered);
    }
  }

  return (myGH);
}

/*@@
  @routine    CopyParFile
  @date       Tue 18 Dec 2001
  @author     Thomas Radke
  @desc
              Copies the original parameter file to a new one in 'IO::out_dir'.
              Note that the new parameter file will usually overwrite an
              existing file unless
                - the new parameter file is identical with the original one
                - this is a recovery run
  @enddesc
  @calls      CCTK_ParameterFilename

  @returntype int
  @returndesc
               0 for success, or<BR>
              -1 if original parameter file couldn't be opened for reading,<BR>
              -2 if new parameter file couldn't be opened for writing<BR>
  @endreturndesc
@@*/
static int CopyParFile(int recovered) {
  int in_parfile, out_parfile, bytes, bytes_written, flags;
  char *out_parfilename, buffer[1024];
  const char *my_parfile_name;
  struct stat in_stat_buf, out_stat_buf;
  DECLARE_CCTK_PARAMETERS

  /* get the name of the original parfile and open it for reading */
  CCTK_ParameterFilename(sizeof(buffer), buffer);
  in_parfile = open(buffer, O_RDONLY);
  if (in_parfile < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't open original parameter file '%s' (%s)", buffer,
               strerror(errno));
    return (-1);
  }
  /* use the file mode of the original parfile to create the new one */
  if (fstat(in_parfile, &in_stat_buf)) {
    in_stat_buf.st_mode = 0644;
  }

  /* build the name of the output parfile */
  my_parfile_name = parfile_name;
  if (!*my_parfile_name) {
    /* cut off any path names */
    my_parfile_name = strrchr(buffer, '/');
    if (my_parfile_name) {
      my_parfile_name++;
    } else {
      my_parfile_name = buffer;
    }
  }
  out_parfilename = malloc(strlen(out_dir) + strlen(my_parfile_name) + 2);
  sprintf(out_parfilename, "%s/%s", out_dir, my_parfile_name);

  /* check whether input and output files are identical */
  if (!stat(out_parfilename, &out_stat_buf) &&
      in_stat_buf.st_ino == out_stat_buf.st_ino) {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Parameter file '%s' to be written into directory '%s' is "
               "identical with original parameter file. Parameter file will "
               "not be copied.",
               my_parfile_name, out_dir);
    out_parfile = 0;
  } else {
    /* binary-copy the input parfile to the output parfile */
    flags = O_CREAT | O_TRUNC | O_WRONLY;
    if (recovered) {
      flags |= O_EXCL;
    }
    out_parfile = open(out_parfilename, flags, in_stat_buf.st_mode);
    if (out_parfile >= 0) {
      while ((bytes = read(in_parfile, buffer, sizeof(buffer))) > 0) {
        bytes_written = write(out_parfile, buffer, bytes);
        if (bytes_written != bytes) {
          if (bytes_written < 0) {
            CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Error while writing parameter file '%s' (%s)",
                       out_parfilename, strerror(errno));
          } else {
            CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Error while writing parameter file '%s'",
                       out_parfilename);
          }
        }
      }
      close(out_parfile);
    } else {
      CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Couldn't write parameter file '%s' (%s)", out_parfilename,
                 strerror(errno));
    }
  }

  /* clean up */
  close(in_parfile);
  free(out_parfilename);

  return (out_parfile >= 0 ? 0 : -2);
}

/*@@
  @routine    DumpParameters
  @date       Tue 18 Dec 2001
  @author     Thomas Radke
  @desc
              Generates a new parameter file in 'IO::out_dir' with an
              ActiveThorns list and a sorted list of all active parameters
              which have been set in the original parameter file.<BR>
              Note that the new parameter file will usually overwrite an
              existing file unless
                - the new parameter file is identical with the original one
                - this is a recovery run
  @enddesc
  @calls      CCTK_ParameterFilename

  @var        recovered
  @vdesc      flag indicating whether this is a recovery run
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
               0 for success, or<BR>
              -1 if original parameter file couldn't be stat(2),<BR>
              -2 if new parameter file couldn't be opened for writing<BR>
  @endreturndesc
@@*/
static int GenerateParFile(int recovered) {
  FILE *outfile;
  int out_parfile, flags;
  char *out_parfilename, buffer[1024];
  const char *my_parfile_name;
  struct stat in_stat_buf, out_stat_buf;
  DECLARE_CCTK_PARAMETERS

  /* get the name of the original parfile and stat(2) it */
  CCTK_ParameterFilename(sizeof(buffer), buffer);
  if (stat(buffer, &in_stat_buf)) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't stat(2) original parameter file '%s' (%s)", buffer,
               strerror(errno));
    return (-1);
  }

  /* build the name of the output parfile */
  my_parfile_name = parfile_name;
  if (!*my_parfile_name) {
    /* cut off any path names */
    my_parfile_name = strrchr(buffer, '/');
    if (my_parfile_name) {
      my_parfile_name++;
    } else {
      my_parfile_name = buffer;
    }
  }
  out_parfilename = malloc(strlen(out_dir) + strlen(my_parfile_name) + 2);
  sprintf(out_parfilename, "%s/%s", out_dir, my_parfile_name);

  /* check whether input and output files are identical */
  if (!stat(out_parfilename, &out_stat_buf) &&
      in_stat_buf.st_ino == out_stat_buf.st_ino) {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Parameter file '%s' to be written into directory '%s' is "
               "identical with original parameter file. Parameter file will "
               "not be generated.",
               my_parfile_name, out_dir);
    out_parfile = 0;
  } else {
    flags = O_CREAT | O_TRUNC | O_WRONLY;
    if (recovered) {
      flags |= O_EXCL;
    }

    out_parfile = open(out_parfilename, flags, in_stat_buf.st_mode);
    if (out_parfile >= 0) {
      outfile = fdopen(out_parfile, "w");
    }
    if (out_parfile >= 0 && outfile) {
      fprintf(outfile, "# '%s' automatically generated by Cactus version %s\n",
              out_parfilename, CCTK_FullVersion());
      fprintf(outfile, "# Original parameter file was '%s'\n", buffer);
      Util_CurrentTime(sizeof(buffer), buffer);
      fprintf(outfile, "# Run time/date was %s ", buffer);
      Util_CurrentDate(sizeof(buffer), buffer);
      fprintf(outfile, "%s ", buffer);
      Util_GetHostName(buffer, sizeof(buffer));
      fprintf(outfile, "on host '%s' with %d processor(s)\n\n", buffer,
              CCTK_nProcs(NULL));
      DumpParameters(outfile);
      fclose(outfile);
    } else {
      CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Couldn't write parameter file '%s' (%s)", out_parfilename,
                 strerror(errno));
    }
    if (out_parfile >= 0) {
      close(out_parfile);
    }
  }

  /* clean up */
  free(out_parfilename);

  return (out_parfile >= 0 ? 0 : -2);
}

/*@@
  @routine    DumpParameters
  @date       Thu Oct 25 2001
  @author     Thomas Radke
  @desc
              Dumps parameter settings of all active thorns into a file.
  @enddesc
  @calls      CCTK_NumCompiledThorns
              CCTK_CompiledThorn
              CCTK_IsThornActive
              CCTK_ParameterWalk
              CCTK_ParameterValString

  @var        parfile
  @vdesc      open file descriptor for the output file
  @vtype      FILE *
  @vio        out
  @endvar

  @returntype int
  @returndesc
               0 for success, or<BR>
              -1 if file descriptor is invalid
  @endreturndesc
@@*/
static int DumpParameters(FILE *outfile) {
  int thorn, first, num_thorns;
  int len, maxname_len;
  const char *thornname;
  char *param;
  const cParamData *pdata;
  t_param_list *list, *last, *new;
  const char *quoted_format = "%s::%-*s = \"%s\"\n",
             *non_quoted_format = "%s::%-*s = %s\n";

  /* check passed file descriptor */
  if (outfile == NULL) {
    return (-1);
  }

  num_thorns = CCTK_NumCompiledThorns();

  /* loop over all thorns to print the 'ActiveThorns' list */
  first = 1;
  fputs("ActiveThorns = \"", outfile);
  for (thorn = 0; thorn < num_thorns; thorn++) {
    thornname = CCTK_CompiledThorn(thorn);

    /* skip all inactive thorns and "Cactus" */
    if (strcmp("Cactus", thornname) && CCTK_IsThornActive(thornname)) {
      if (!first) {
        putc(' ', outfile);
      }
      fputs(thornname, outfile);
      first = 0;
    }
  }
  fputs("\"\n", outfile);

  /* loop over all thorns */
  for (thorn = 0; thorn < num_thorns; thorn++) {
    thornname = CCTK_CompiledThorn(thorn);

    /* skip all inactive thorns */
    if (!CCTK_IsThornActive(thornname)) {
      continue;
    }

    /* now walk through all parameters of given thorn */
    maxname_len = 0;
    list = last = NULL;
    first = 1;
    while (CCTK_ParameterWalk(first, thornname, &param, &pdata) == 0) {
      first = 0;

      /* skip the parameters which weren't explicitely set */
      if (pdata->n_set > 0) {
        new = malloc(sizeof(t_param_list));
        if (new) {
          new->value = CCTK_ParameterValString(pdata->name, pdata->thorn);
          if (new->value) {
            new->format_string =
                pdata->type == PARAMETER_INT || pdata->type == PARAMETER_REAL
                    ? non_quoted_format
                    : quoted_format;
            new->name = pdata->name;
            len = strlen(new->name);
            if (len > maxname_len) {
              maxname_len = len;
            }
            if (last) {
              last->next = new;
              last = new;
            } else {
              list = last = new;
            }
            last->next = NULL;
          } else {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Couldn't get value for parameter '%s'", param);
            free(new);
          }
        } else {
          CCTK_WARN(1, "Couldn't allocate list element");
        }
      }
      free(param);
    } /* end of loop walking over all parameters of given thorn */

    /* finally dump out the list for this thorn */
    if (list) {
      fprintf(outfile, "\n# Parameters of thorn %s (implementing %s)\n",
              thornname, CCTK_ThornImplementation(thornname));
    }
    while (list) {
      fprintf(outfile, list->format_string, thornname, maxname_len, list->name,
              list->value);
      free(list->value);
      last = list->next;
      free(list);
      list = last;
    }
  } /* end of looping over all thorns */

  return (0);
}
