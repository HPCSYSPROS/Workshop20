 /*@@
   @file      CommandLine.c
   @date      Wed Feb 17 00:11:26 1999
   @author    Tom Goodale
   @desc
              Routines to deal with command line arguments.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CommandLine.h"
#include "cctk_Flesh.h"
#include "cGH.h"

#include "util_String.h"
#include "cctk_Version.h"
#include "cctk_ActiveThorns.h"
#include "cctk_CommandLine.h"
#include "cctk_Comm.h"
#include "cctk_File.h"
#include "cctk_Misc.h"
#include "cctk_ParamCheck.h"
#include "cctk_WarnLevel.h"

#include "cctki_ActiveThorns.h"
#include "cctki_WarnLevel.h"


#define NEED_PARAMETER_SCOPE_STRINGS
#define NEED_PARAMETER_TYPE_STRINGS

#include "cctk_Parameter.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_CommandLine_c);

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
static void CommandLinePrintParameter (const cParamData *properties);


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/
static char* logdir = NULL;
static redirect_t requested_stdout_redirection = REDIRECT_NONE;
static redirect_t requested_stderr_redirection = REDIRECT_NONE;
static int buffering_type = 0;
/* buffering: 0=default, 1=unbuffered, 2=line, 3=fully */
static int paramchecking = 0;
/* used to keep track of whether to inform the user that we changed the logging
 * or warning levels */
static int requested_change_in_logginglevel = 0;
static int requested_change_in_warninglevel = 0;
static int requested_change_in_errorlevel = 0;
static int requested_logginglevel = 0;
static int requested_warninglevel = 0;
static int requested_errorlevel = 0;


/********************************************************************
 *********************     Global Data   *****************************
 ********************************************************************/
int cctki_paramchecking;
int cctki_paramcheck_nprocs;

int cctki_onlyprintschedule = 0;

/********************************************************************
 *********************        Defines          **********************
 ********************************************************************/

/*
 * See comments in  CCTKi_ProcessCommandLine()  for a list of all the
 * places you have to update if you add a new command-line option.
 */
#define CACTUS_COMMANDLINE_OPTIONS                                      \
        "[-h] [-O] [-o paramname] [-L n] [-W n] [-E n] "                \
        "[-r[o|e|oe|eo]] [--logdir <dir>] "                             \
        "[-b <no|line|full>] "                                          \
        "[-S] [-T] [-t name] [--parameter-level <level>] [-v] "         \
        "<parameter_file_name>"


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_CommandLineTestThorncompiled
   @date       Wed Feb 17 10:25:30 1999
   @author     Gabrielle Allen
   @desc
               Tests if a given thorn has been compiled.
   @enddesc
   @calls      CCTK_IsThornCompiled

   @var        argument
   @vdesc      thorn name
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineTestThornCompiled (const char *argument)
{
  int retval;


  retval = CCTK_IsThornCompiled (argument);
  if (CCTK_MyProc (NULL) == 0)
  {
    printf ("Thorn '%s' %savailable.\n", argument, retval ? "" : "un");
  }

  CCTK_Exit (NULL, retval);
}


 /*@@
   @routine    CCTKi_CommandLineDescribeAllParameters
   @date       Tue Apr 18 15:00:12 2000
   @author     Tom Goodale
   @desc
               Describe all the parameters
   @enddesc
   @calls      CCTK_NumCompiledThorns
               CCTK_CompiledThorn
               CCTK_ParameterWalk
               CommandLinePrintParameter

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineDescribeAllParameters (const char *argument)
{
  int first, n_thorns, thorn;
  char *param;
  const char *thornname;
  const cParamData *properties;


  if (CCTK_MyProc (NULL) == 0)
  {
    n_thorns = CCTK_NumCompiledThorns ();

    for (thorn = 0; thorn < n_thorns; thorn++)
    {
      thornname = CCTK_CompiledThorn (thorn);
      printf ("\nParameters of thorn '%s' providing implementation '%s':\n",
              thornname, CCTK_ThornImplementation (thornname));

      first = 1;
      while (CCTK_ParameterWalk (first, thornname, &param, &properties) == 0)
      {
        if (argument)
        {
          switch (*argument)
          {
            case 'v':
              CommandLinePrintParameter (properties);
              break;
            default :
              fprintf (stderr, "Unknown verbosity option %s\n", argument);
              CCTK_Exit (NULL, 2);
          }
        }
        else
        {
          printf ("%s\n", param);
        }

        free (param);
        first = 0;
      }
    }
  }

  CCTK_Exit (NULL, 0);
}


 /*@@
   @routine    CCTKi_CommandLineDescribeParameter
   @date       Tue Apr 18 15:00:33 2000
   @author     Tom Goodale
   @desc
               Describe a particular parameter.
   @enddesc
   @calls      Util_SplitString
               CCTK_ParameterData
               CCTK_ImplementationThorn
               CommandLinePrintParameter

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineDescribeParameter (const char *argument)
{
  int retcode;
  char *thorn, *param;
  const char *cthorn;
  const cParamData *properties;


  Util_SplitString (&thorn, &param, argument, "::");

  if (! param)
  {
    properties = CCTK_ParameterData (argument, NULL);
  }
  else
  {
    properties = CCTK_ParameterData (param, thorn);
    if (! properties)
    {
      cthorn = CCTK_ImplementationThorn (thorn);
      properties = CCTK_ParameterData (param, cthorn);
    }

    free (thorn);
    free (param);
  }

  if(properties)
  {
    if (CCTK_MyProc (NULL) == 0)
    {
      CommandLinePrintParameter (properties);
    }
    retcode = 0;
  }
  else
  {
    if (CCTK_MyProc (NULL) == 0)
    {
      fprintf(stderr, "No such parameter\n");
    }
    retcode = 1;
  }

  CCTK_Exit (NULL, retcode);
}


 /*@@
   @routine    CCTKi_CommandLineTestParameters
   @date       Tue Apr 18 15:00:45 2000
   @author     Tom Goodale
   @desc

   @enddesc
   @calls      CCTKi_CommandLineUsage

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineTestParameters (const char *argument)
{
  int nprocs;
  char *endptr;


  if (argument == NULL)
  {
    nprocs = 1;
  }
  else
  {
    nprocs = strtol (argument, &endptr, 10);
    if (! (endptr && *endptr == 0))
    {
      CCTKi_CommandLineUsage ();
    }
  }

  paramchecking = 1;
  cctki_paramchecking = 1;
  cctki_paramcheck_nprocs = nprocs;
}


/*@@
   @routine    CCTKi_CommandLineLoggingLevel
   @date       Thu Dec 30 2004
   @author     Erik Schnetter
   @desc
               Sets the CCTK logging level from a command line argument.
   @enddesc
   @calls      CCTKi_SetLoggingLevel
               CCTKi_CommandLineUsage

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineLoggingLevel (const char *argument)
{
  char *endptr;
  long int logginglevel;


  logginglevel = strtol (argument, &endptr, 10);
  if (endptr && *endptr == 0)
  {
    requested_logginglevel = (int)logginglevel;
    requested_change_in_logginglevel = CCTKi_SetLogLevel (logginglevel);
  }
  else
  {
    CCTKi_CommandLineUsage ();
  }
}


/*@@
   @routine    CCTKi_CommandLineWarningLevel
   @date       Wed Feb 17 00:58:56 1999
   @author     Tom Goodale
   @desc
               Sets the CCTK warning level from a command line argument.
   @enddesc
   @calls      CCTKi_SetWarningLevel
               CCTKi_CommandLineUsage

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineWarningLevel (const char *argument)
{
  char *endptr;
  long int warninglevel;


  warninglevel = strtol (argument, &endptr, 10);
  if (endptr && *endptr == 0)
  {
    requested_warninglevel = (int)warninglevel;
    requested_change_in_warninglevel = CCTKi_SetWarnLevel (warninglevel);
  }
  else
  {
    CCTKi_CommandLineUsage ();
  }
}


 /*@@
   @routine    CCTKi_CommandLineErrorLevel
   @date       Wed Feb 17 00:58:56 1999
   @author     Tom Goodale
   @desc
               Sets the CCTK error level from a command line argument.
   @enddesc
   @calls      CCTKi_SetErrorLevel
               CCTKi_CommandLineUsage

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineErrorLevel (const char *argument)
{
  char *endptr;
  long int errorlevel;


  errorlevel = strtol (argument, &endptr, 10);
  if (endptr && *endptr == 0)
  {
    if (errorlevel < 0)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "Error level cannot be negative, but %d was requested.",
                  (int)errorlevel);
      CCTK_Exit (NULL, 1);
    }
    else
    {
      requested_errorlevel = (int)errorlevel;
      requested_change_in_errorlevel = CCTKi_SetErrorLevel (errorlevel);
    }
  }
  else
  {
    CCTKi_CommandLineUsage ();
  }
}


 /*@@
   @routine    CCTKi_CommandLineParameterLevel
   @date       Wed Feb 21 2001
   @author     Gabrielle Allen
   @desc
               Sets the parameter checking level from a command line argument.
   @enddesc
   @calls      CCTKi_SetParameterLevel

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineParameterLevel (const char *argument)
{
  int parameterlevel;


  if (CCTK_Equals (argument, "strict"))
  {
    parameterlevel = CCTK_PARAMETER_STRICT;
  }
  else if (CCTK_Equals (argument, "normal"))
  {
    parameterlevel = CCTK_PARAMETER_NORMAL;
  }
  else if (CCTK_Equals (argument, "relaxed"))
  {
    parameterlevel = CCTK_PARAMETER_RELAXED;
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "Parameter checking level '%s' not recognised, "
                "defaulting to normal", argument);
    parameterlevel = CCTK_PARAMETER_NORMAL;
  }

  CCTKi_SetParameterLevel (parameterlevel);
}


 /*@@
   @routine    CCTKi_CommandLineRedirect
   @date       Fri Jul 23 11:32:46 1999
   @author     Tom Goodale
   @desc
               Redirect stdout/stderr on non-root processors into a file

               Redirection is defered until all command line options
               (including a possible '-logdir <dir>') have been processed.
   @enddesc

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
   @var        type
   @vdesc      redirect type (1: all MPI ranks > 0, 2: all MPI ranks)
   @vtype      const int
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineRedirect (const char *argument, const redirect_t type)
{
  if (!argument || strchr(argument,'o')) /* redirect stdout */
  {
    requested_stdout_redirection = type;
  }
  if (argument && strchr(argument,'e')) /* redirect stderr */
  {
    requested_stderr_redirection = type;
  }
}


 /*@@
   @routine    CCTKi_CommandLineLogDir
   @date       Thu 22 July 2006
   @author     Thomas Radke
   @desc
               Set the output directory for redirected stdout/stderr logfiles.
   @enddesc

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineLogDir (const char *argument)
{
  logdir = Util_Strdup (argument);
}
 /*@@
   @routine    CCTKi_CommandLineSetBuffering
   @date       2006-05-27
   @author     Erik Schnetter
   @desc
               Set stdout buffering.  (stderr is always unbuffered.)
   @enddesc

   @var        argument
   @vdesc      option argument
   @vtype      const char *
   @vio        in
   @endvar
@@*/
void CCTKi_CommandLineSetBuffering (const char *argument)
{
  if (! strcmp (argument, "no"))
  {
    /* Switch to unbuffered mode (best for debugging) */
    buffering_type = 1;
  }
  else if (! strcmp (argument, "line"))
  {
    /* Switch to line buffered mode (good for screen output) */
    buffering_type = 2;
  }
  else if (! strcmp (argument, "full"))
  {
    /* Switch to fully buffered mode (fastest) */
    buffering_type = 3;
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "Stdout buffering mode '%s' not recognised, "
                "not changing the default setting", argument);
  }
}


 /*@@
   @routine    CCTKi_CommandLinePrintSchedule
   @date       2005-06-10
   @author     Erik Schnetter
   @desc
               Set a flag that makes the run abort after printing the
               schedule tree.
   @enddesc
@@*/
void CCTKi_CommandLinePrintSchedule (void)
{
  cctki_onlyprintschedule = 1;
}


 /*@@
   @routine    CCTKi_CommandLineListThorns
   @date       Tue Apr 18 15:05:00 2000
   @author     Tom Goodale
   @desc
               List the thorns which are compiled in.
   @enddesc
   @calls      CCTKi_PrintThorns
@@*/
void CCTKi_CommandLineListThorns (void)
{
  if (CCTK_MyProc (NULL) == 0)
  {
    printf ("\n---------------Compiled Thorns-------------\n");
    CCTKi_PrintThorns (stdout, "  %s\n", 0);
    printf ("-------------------------------------------\n\n");
  }

  CCTK_Exit (NULL, 1);
}


 /*@@
   @routine    CCTKi_CommandLineVersion
   @date       Fri Jul 23 12:57:45 1999
   @author     Tom Goodale
   @desc
               Prints version info.
   @enddesc
   @calls      CCTK_FullVersion
               CCTK_CompileDate
               CCTK_CompileTime
@@*/
void CCTKi_CommandLineVersion(void)
{
  char **argv;
  const char *version;


  if (CCTK_MyProc (NULL) == 0)
  {
    CCTK_CommandLine (&argv);
    version = CCTK_FullVersion ();

    printf ("%s: Version %s.  Compiled on %s at %s\n", argv[0], version,
            CCTK_CompileDate (), CCTK_CompileTime ());
  }

  CCTK_Exit (NULL, 1);
}


 /*@@
   @routine    CCTKi_CommandLineHelp
   @date       Fri Jul 23 12:57:23 1999
   @author     Tom Goodale
   @desc
               Prints a help message.
   @enddesc
   @calls      CCTK_CommandLine
               CCTK_CompileDate
               CCTK_CompileTime
@@*/
void CCTKi_CommandLineHelp (void)
{
  char **argv;

  /*
   * See comments in  CCTKi_ProcessCommandLine()  for a list of all the
   * places you have to update if you add a new command-line option.
   */
  const char *commandline_options_description =
    "-h, --help                           : gets this help.\n"
    "-O[v], --describe-all-parameters     : describes all the parameters.\n"
    "                                       v makes this verbose, i.e., it gives\n"
    "                                       a verbose description of all parameters.\n"
    "-o, --describe-parameter <paramname> : describe the given parameter.\n"
    "-L, --logging-level <n>              : Sets the logging level to n.\n"
    "-W, --warning-level <n>              : Sets the warning level to n.\n"
    "-E, --error-level <n>                : Sets the error level to n.\n"
    "-r, --redirect[o|e|oe|eo]            : Redirects standard output and/or\n"
    "                                       standard error to files;\n"
    "                                       only MPI-non-root processes.\n"
    "-R, --Redirect[o|e|oe|eo]            : Redirects standard output and/or\n"
    "                                       standard error to files;\n"
    "                                       all MPI processes (no screen output).\n"
    "    --logdir <dir>                   : Sets the output directory for logfiles\n"
    "                                       created by the '-r' option\n"
    "-b, --buffering <no|line|full>       : Set stdout buffering mode.\n"
    "-S, --print-schedule                 : Print the schedule tree, then exit.\n"
    "-P, --exit-after-param-check         : Check the parameters, then exit.\n"
    "-T, --list-thorns                    : Lists the compiled-in thorns.\n"
    "-t, --test-thorn-compiled <name>     : Tests for the presence of thorn <name>.\n"
    "    --parameter-level <level>        : Sets the amount of parameter checking, \n"
    "                                       level can be strict, normal, relaxed.\n"
    "-v, --version                        : Prints the version.\n"
    "-i, --ignore-next                    : Ignores the next argument.\n";

  /* test-parameter option to be added back when implemented 
    "-x, --test-parameters [nprocs]       : does a quick test of the parameter file\n"
    "                                       pretending to be on nprocs processors,\n"
    "                                       or 1 if not given.\n"
   */

  if (CCTK_MyProc (NULL) == 0)
  {
    CCTK_CommandLine (&argv);

    printf ("%s, compiled on %s at %s\n",
            argv[0], CCTK_CompileDate(), CCTK_CompileTime());
    printf ("Usage: %s %s\n", argv[0], CACTUS_COMMANDLINE_OPTIONS);
    printf ("\nValid options:\n%s", commandline_options_description);
  }

  CCTK_Exit (NULL, 1);
}


 /*@@
   @routine    CCTKi_CommandLineUsage
   @date       Fri Jul 23 12:57:04 1999
   @author     Tom Goodale
   @desc
               Prints a usage message.
   @enddesc
   @calls      CCTK_CommandLine
@@*/
void CCTKi_CommandLineUsage (void)
{
  char **argv;


  if (CCTK_MyProc (NULL) == 0)
  {
    CCTK_CommandLine (&argv);

    printf ("Usage: %s %s\n", argv[0], CACTUS_COMMANDLINE_OPTIONS);
  }

  CCTK_Exit (NULL, 1);
}


 /*@@
   @routine    CCTKi_CommandLineFinished
   @date       Fri Jul 23 12:55:39 1999
   @author     Tom Goodale
   @desc
               Subroutine to do anything which has to be done based upon the
               commandline, but needs to be have a default.
   @enddesc
@@*/
void CCTKi_CommandLineFinished (void)
{
  int myproc;
  char *logfilename;
  FILE *newfile;


  /* Are we in a paramcheck run ? */
  if (! paramchecking)
  {
    cctki_paramchecking = 0;
  }

  /* redirect stdout/stderr on non-root processors */
  if (logdir && requested_stdout_redirection == REDIRECT_NONE &&
                requested_stderr_redirection == REDIRECT_NONE)
  {
    CCTK_VWarn (CCTK_WARN_PICKY, __LINE__, __FILE__, "Cactus",
                "Specifying the '-logdir' option without the '-r' or '-R' "
                "option is a no-op and will be ignored.");
  }
  myproc = CCTK_MyProc (NULL);
  /* if specified on the command line, create the output directory
     for redirected stdout/stderr logfiles */
  if (logdir)
  {
    if (requested_stdout_redirection != REDIRECT_NONE ||
        requested_stderr_redirection != REDIRECT_NONE)
    {
      if (CCTK_CreateDirectory (0755, logdir) < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                    "Could not create output directory '%s' for "
                    "stdout/stderr logfiles ! Falling back to using the "
                    "current working directory...", logdir);
        free (logdir);
        logdir = Util_Strdup (".");
      }
    }
  }
  else
  {
    /* make cwd the default logdir */
    logdir = Util_Strdup (".");
  }

  /* if redirection was requested on the command line
     send stdout/stderr messages to <logdir>/CCTK_Proc<id>.{out,err}
     otherwise redirect stdout to the NULL device */
  logfilename = malloc (strlen (logdir) + 32);
  if ( (myproc && requested_stdout_redirection != REDIRECT_NONE) ||
       requested_stdout_redirection == REDIRECT_ALL )
  {
    if (myproc == 0 && requested_stdout_redirection == REDIRECT_ALL)
      printf("Redirection of all stdout to file(s) was requested. This means "
             "that there will be no output to the screen. In order to see the "
             "redirected output you will need to look at these files, e.g., "
             "by using \"tail -f %s/CCTK_Proc0.out\".\n", logdir);
    sprintf (logfilename, "%s/CCTK_Proc%u.out", logdir, myproc);
    newfile = freopen (logfilename, "w", stdout);
    if (! newfile)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Could not redirect stdout to logfile '%s'", logfilename);
    }
  }
  else if (myproc && requested_stdout_redirection == REDIRECT_NONE)
  {
    newfile = freopen (NULL_DEVICE, "w", stdout);
    if (! newfile)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Could not disable stdout "
                  "(was trying to redirect it to '%s')", NULL_DEVICE);
    }
  }

  if ( (myproc && requested_stderr_redirection != REDIRECT_NONE) ||
       requested_stderr_redirection == REDIRECT_ALL )
  {
    sprintf (logfilename, "%s/CCTK_Proc%u.err", logdir, myproc);
    newfile = freopen (logfilename, "w", stderr);
    if (! newfile)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Could not redirect stderr to logfile '%s'", logfilename);
    }
  }
  free (logfilename);
  free (logdir);
  
  /* if requested, change buffering mode of stdout */
  switch (buffering_type)
  {
  case 0:
    /* Keep default */
    break;
  case 1:
    /* Switch to unbuffered mode (best for debugging) */
    setvbuf (stdout, NULL, _IONBF, 0);
    break;
  case 2:
    /* Switch to line buffered mode (good for screen output) */
    setvbuf (stdout, NULL, _IOLBF, 0);
    break;
  case 3:
    /* Switch to fully buffered mode (fastest) */
    setvbuf (stdout, NULL, _IOFBF, 0);
    break;
  }
  
  /* ensure that stderr is unbuffered (best for debugging) */
  setvbuf (stderr, NULL, _IONBF, 0);

  /* output information about changed log levels now that log output is set up
   * correctly */
  if (requested_change_in_warninglevel)
  {
    CCTK_VInfo ("Cactus", "%s warning level to %d",
                requested_change_in_warninglevel > 0 ? "Increased" :
                "Decreased", requested_warninglevel);
  }
  if (requested_change_in_errorlevel)
  {
    CCTK_VInfo ("Cactus", "%s error level to %d",
                requested_change_in_errorlevel > 0 ? "Increased" :
                "Decreased", requested_errorlevel);
  }
  if (requested_change_in_logginglevel)
  {
    CCTK_VInfo ("Cactus", "%s logging level to %d",
                requested_change_in_logginglevel > 0 ? "Increased" :
                "Decreased", requested_logginglevel);
  }
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    CommandLinePrintParameter
   @date       Sun Oct 17 22:11:31 1999
   @author     Tom Goodale
   @desc
               Prints info for a parameter.
   @enddesc
   @calls      CCTK_ThornImplementation

   @var        properties
   @vdesc      Parameter properties
   @vtype      const cParamData *
   @vio        in
   @endvar
@@*/
static void CommandLinePrintParameter (const cParamData *properties)
{
  t_range *range;


  if (properties)
  {
    printf ("Parameter:   %s::%s", properties->thorn, properties->name);
    if (properties->scope != SCOPE_PRIVATE)
    {
      printf (", %s::%s", CCTK_ThornImplementation (properties->thorn),
                          properties->name);
    }
    printf ("\n");
    printf ("Description: \"%s\"\n", properties->description);
    printf ("Type:        %s\n", cctk_parameter_type_names[properties->type-PARAMETER_FIRST]);
    printf ("Default:     %s\n", properties->defval);
    printf ("Scope:       %s\n", cctk_parameter_scopes[properties->scope-SCOPE_FIRST]);

    for (range = properties->range; range; range = range->next)
    {
      printf ("  Range:     %s\n", range->range);
      printf ("    Origin:      %s\n", range->origin);
      printf ("    Description: %s\n", range->description);
    }
  }
}
