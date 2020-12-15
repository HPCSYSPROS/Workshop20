 /*@@
   @file      ProcessCommandLine.c
   @date      Thu Sep 24 10:32:28 1998
   @author    Tom Goodale
   @desc 
   Routines to deal with the command line arguments.
   @enddesc 
   @version $Header$
 @@*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk_Flesh.h"
#include "cctk_GNU.h"
#include "cctk_Misc.h"
#include "cctk_CommandLine.h"

#include "CommandLine.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_ProcessCommandLine_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

void CCTK_FCALL CCTK_FNAME(CCTK_ParameterFilename)
     (int *retval, int *len, char *name);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static int exit_after_param_check = 0;

static char *parameter_file_name=NULL;

static int argc = 0;

static char **argv = NULL;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_ProcessCommandLine
   @date       Thu Sep 24 10:33:31 1998
   @author     Tom Goodale
   @desc 
   Processes the command line arguments.
   @enddesc 
   @calls    CCTKi_CommandLineTestThornCompiled CCTKi_CommandLineDescribeAllParameters CCTKi_CommandLineDescribeParameter CCTKi_CommandLineTestParameters CCTKi_CommandLineLoggingLevel CCTKi_CommandLineWarningLevel CCTKi_CommandLineErrorLevel CCTKi_CommandLineParameterLevel CCTKi_CommandLineRedirectStdout CCTKi_CommandLineListThorns() CCTKi_CommandLineVersion() CCTKi_CommandLineHelp
   @calledby   
   @history 
 
   @endhistory 
   @var     inargc
   @vdesc   Number of runtime arguments
   @vtype   int *
   @vio     inout
   @vcomment 
 
   @endvar 
   @var     inargv
   @vdesc   Command line arguments
   @vtype   char ***
   @vio     inout
   @vcomment 
 
   @endvar 
   @var     ConfigData
   @vdesc   Flesh configuration data
   @vtype   tFleshConfig
   @vio     inout
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc
   0 -- success
   @endreturndesc

@@*/
int CCTKi_ProcessCommandLine(int *inargc, char ***inargv, tFleshConfig *ConfigData)
{
  /*
   * If you add a new command-line option, you must update (at least)
   * the following different places in the code:
   * - the definition of  long_options[]  in this function
   * - the 3rd argument in the call to  getopt_long_only()  in this function
   * - the  switch (c)  statement in this function
   * - the #define of CACTUS_COMMANDLINE_OPTIONS near the top of
   *   src/main/CommandLine.c
   * - the help message printed by CCTKi_CommandLineHelp()
   *   (also in  src/main/CommandLine.c )
   * You should also update the description of command-line options in the
   * Cactus Users' Guide, in  doc/UsersGuide/RunningCactus.tex .
   */

  int option_index = 0;
  int c;
  int ignore;
  /* constants for identifying each options by the return value
     from getopt_long_only() */
  enum
  {
    help_option                   = 'h',
    describe_all_paramters_option = 'O',
    describe_parameter_option     = 'o',
    test_parameters_option        = 'x',
    logging_level_option          = 'L',
    warning_level_option          = 'W',
    error_level_option            = 'E',
    parameter_level_option        = 256,  /* no short option */
    redirect_option               = 'r',
    Redirect_option               = 'R',
    logdir_option                 = 257,  /* no short option */
    buffering_option              = 'b',
    print_schedule_option         = 'S',
    list_thorns_option            = 'T',
    test_thorns_compiled_option   = 't',
    version_option                = 'v',
    exit_after_param_check_option = 'P',
    ignore_next_option            = 'i'
  };
  /* the longopts argument passed into getopt_long_only() */
  const struct option long_options[] =
  {
    {"help",                    no_argument,       NULL, help_option},
    {"describe-all-parameters", optional_argument, NULL, describe_all_paramters_option},
    {"describe-parameter",      required_argument, NULL, describe_parameter_option},
    /*{"test-parameters",         optional_argument, NULL, test_parameters_option},*/
    {"logging-level",           required_argument, NULL, logging_level_option},
    {"warning-level",           required_argument, NULL, warning_level_option},
    {"error-level",             required_argument, NULL, error_level_option},
    {"parameter-level",         required_argument, NULL, parameter_level_option},
    {"redirect",                optional_argument, NULL, redirect_option},
    {"Redirect",                optional_argument, NULL, Redirect_option},
    {"logdir",                  required_argument, NULL, logdir_option},
    {"buffering",               required_argument, NULL, buffering_option},
    {"print-schedule",          no_argument,       NULL, print_schedule_option},
    {"list-thorns",             no_argument,       NULL, list_thorns_option},
    {"test-thorn-compiled",     required_argument, NULL, test_thorns_compiled_option},
    {"exit-after-param-check",  no_argument,       NULL, exit_after_param_check_option},
    {"version",                 no_argument,       NULL, version_option},
    {"ignore-next",             no_argument,       NULL, ignore_next_option},
    {0, 0, 0, 0}
  };


  /* Store the command line */
  argc = *inargc;
  argv = *inargv;

  ignore = 0;

  if(argc>1)
  {
    while (1)
    {
      c = getopt_long_only (argc, argv, "hO::o:x::L:W:E:r::R::b:STt:Pvi",
                            long_options, &option_index);
      if (c == -1)
        break;
  
      if(!ignore)
      {
        switch (c)
        {
          case describe_all_paramters_option:
            CCTKi_CommandLineDescribeAllParameters(optarg); break;
          case describe_parameter_option:
            CCTKi_CommandLineDescribeParameter(optarg); break;
          case test_parameters_option:
            CCTKi_CommandLineTestParameters(optarg); break;
          case logging_level_option:
            CCTKi_CommandLineLoggingLevel(optarg); break;
          case warning_level_option:
            CCTKi_CommandLineWarningLevel(optarg); break;
          case error_level_option:
            CCTKi_CommandLineErrorLevel(optarg); break;
          case parameter_level_option:
            CCTKi_CommandLineParameterLevel(optarg); break;
          case redirect_option:
            CCTKi_CommandLineRedirect(optarg, REDIRECT_NONROOT); break;
          case Redirect_option:
            CCTKi_CommandLineRedirect(optarg, REDIRECT_ALL); break;
          case logdir_option:
            CCTKi_CommandLineLogDir(optarg); break;
          case buffering_option:
            CCTKi_CommandLineSetBuffering(optarg); break;
          case print_schedule_option:
            CCTKi_CommandLinePrintSchedule(); break;
          case list_thorns_option:
            CCTKi_CommandLineListThorns(); break;
          case test_thorns_compiled_option:
            CCTKi_CommandLineTestThornCompiled(optarg); break;
          case version_option:
            CCTKi_CommandLineVersion(); break;
          case ignore_next_option:
            ignore = 1; break;
          case exit_after_param_check_option:
            exit_after_param_check = 1; break;
          case help_option: 
          case '?':
            CCTKi_CommandLineHelp(); break;
          default:
            printf ("?? getopt returned character code 0%o ??\n", c);
        }
      }
      else
      {
        printf("Ignoring option\n");
        ignore = 0;
      }
    }

    if(argc > optind)
    {
      ConfigData->parameter_file_name = argv[optind];
      parameter_file_name = ConfigData->parameter_file_name;
    }
    else
    {
      CCTKi_CommandLineUsage();
    }
  }
  else
  {
    CCTKi_CommandLineUsage();
  }

  CCTKi_CommandLineFinished();

  return 0;
}


 /*@@
   @routine    CCTK_CommandLine
   @date       Wed Feb 17 00:19:30 1999
   @author     Tom Goodale
   @desc 
   Gets the command line arguments.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     outargv
   @vdesc   Place to dump the command line arguments
   @vtype   char ***
   @vio     out
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc
   The number of command line arguments.
   @endreturndesc
@@*/
int CCTK_CommandLine(char ***outargv)
{
  *outargv = argv;

  return argc;
}


 /*@@
   @routine    CCTK_ParameterFilename
   @date       Tue Oct 3 2000
   @author     Gabrielle Allen
   @desc 
   Returns the parameter filename
   @enddesc 
   @calls    CCTK_Equals 
   @calledby   
   @history 
 
   @endhistory 
   @var     len
   @vdesc   The length of the incoming string
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     filename
   @vdesc   String to contain the filename
   @vtype   char *
   @vio     out
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc
   The length of the returned string.
   @endreturndesc

@@*/
int CCTK_ParameterFilename(int len, char *filename)
{
  int retval;
  const char *copy_string;


  if (CCTK_Equals(parameter_file_name,"-"))
  {
    copy_string = "STDIN";
  }
  else
  {
    copy_string = parameter_file_name;
  }
  retval = strlen (copy_string);
  if (retval > len - 1)
  {
    retval = len - 1;
  }
  strncpy (filename, copy_string, retval);
  filename[retval] = 0;
  return retval;
}

void CCTK_FCALL CCTK_FNAME(CCTK_ParameterFilename)
     (int *retval, int *len, char *name)
{
  *retval = CCTK_ParameterFilename(*len,name);
}

int CCTKi_ExitAfterParamCheck(void)
{
    return exit_after_param_check;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
