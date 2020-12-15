/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "converse.h"
#include "common.h"
#include "BackEnd.h"
#include "InfoStream.h"
#include "Broadcasts.h"

#include "NamdState.h"
#include "Node.h"
#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR _chdir
#define GETCWD _getcwd
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#else
#include <unistd.h>
#define CHDIR chdir
#define GETCWD getcwd
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif
#include <sys/stat.h>
#include "ConfigList.h"
#include "ScriptTcl.h"


void after_backend_init(int argc, char **argv);

#ifdef MEM_OPT_VERSION
//record the working directory when reading the configuration file
//for Parallel IO Input --Chao Mei
char *gWorkDir = NULL;
#endif

int main(int argc, char **argv) {
#ifdef NAMD_TCL
  if ( argc > 2 && ! strcmp(argv[1],"+tclsh") ) {
    // pass all remaining arguments to script
    return ScriptTcl::tclsh(argc-2,argv+2);
  }
#endif
  BackEnd::init(argc,argv);
  after_backend_init(argc, argv);
  return 0;
}

void after_backend_init(int argc, char **argv){
#define CWDSIZE 1024
  char origcwd_buf[CWDSIZE];
  char currentdir_buf[CWDSIZE];

  ScriptTcl *script = new ScriptTcl;
  Node::Object()->setScript(script);

  for(argc = 0; argv[argc]; ++argc);
  if ( argc < 2 ) {
#if defined(WIN32) && !defined(__CYGWIN__)
    CkPrintf("\nFATAL ERROR: No simulation config file specified on command line.\n");
    CkPrintf("\nNOTE: NAMD has no graphical interface and must be run from a command line.\n");
    int nsleep = 10;
    CkPrintf("\nSleeping %d seconds before exiting...\n", nsleep);
    fflush(stdout);
    sleep(nsleep);
    CkPrintf("\n");
#endif
    NAMD_die("No simulation config file specified on command line.");
  }
  char *origcwd = GETCWD(origcwd_buf,CWDSIZE);
  if ( ! origcwd ) NAMD_err("getcwd");
#ifdef NAMD_TCL
  for(int i = 1; i < argc; ++i) {
  if ( strstr(argv[i],"--") == argv[i] ) {
    char buf[1024];
    if ( i + 1 == argc ) {
      sprintf(buf, "missing argument for command line option %s", argv[i]);
      NAMD_die(buf);
    }
    if ( ! strcmp(argv[i],"--tclmain") ) {
      // pass all remaining arguments to script
      iout << iINFO << "Command-line argument is";
      for ( int j=i; j<argc; ++j ) { iout << " " << argv[j]; }
      iout << "\n" << endi;
      script->tclmain(argc-i-1,argv+i+1);
      BackEnd::exit();
      return;
    }
    sprintf(buf, "%s %s", argv[i]+2, argv[i+1]);
    iout << iINFO << "Command-line argument is --" << buf << "\n" << endi;
    script->eval(buf);
    ++i;
    continue;
  }
  char *confFile = argv[i];
#else
  char *confFile = argv[argc-1];
#endif
  iout << iINFO << "Configuration file is " << confFile << "\n" << endi;

  char *currentdir=confFile;
  char *tmp;
  for(tmp=confFile;*tmp;++tmp); // find final null
  for( ; tmp != confFile && *tmp != PATHSEP; --tmp); // find last '/'
#if defined(WIN32) && !defined(__CYGWIN__)
  if (tmp == confFile) {
    // in case this is under cygwin, search for '/' as well
    for(tmp=confFile;*tmp;++tmp); // find final null
    for( ; tmp != confFile && *tmp != '/'; --tmp); // find last '/'
  }
#endif
  if ( CHDIR(origcwd) ) NAMD_err(origcwd);
  if ( tmp != confFile )
  {
    *tmp = 0; confFile = tmp + 1;
    if ( CHDIR(currentdir) ) NAMD_err(currentdir);
    struct stat statBuf;
    if (stat(confFile, &statBuf)) {
      char buf[1024];
      sprintf(buf,"Unable to access config file %s%c%s",currentdir,PATHSEP,confFile);
      NAMD_die(buf);
    }
    iout << iINFO << "Changed directory to " << currentdir << "\n" << endi;
    currentdir = GETCWD(currentdir_buf,CWDSIZE);
    if ( ! currentdir ) NAMD_err("getcwd after chdir");
  }
  else{
      if ( *tmp == PATHSEP ){ // config file in / is odd, but it might happen
          if ( CHDIR(PATHSEPSTR) ) NAMD_err(PATHSEPSTR);
          struct stat statBuf;
          if (stat(confFile, &statBuf)) {
            char buf[1024];
            sprintf(buf,"Unable to access config file %s",confFile);
            NAMD_die(buf);
          }
      }else{ // just a config file name, so the path is the current working path
          struct stat statBuf;
          if (stat(confFile, &statBuf)) {
            char buf[1024];
            if ( confFile[0] == '-' || confFile[0] == '+' ) {
              sprintf(buf,"Unknown command-line option %s",confFile);
            } else {
              sprintf(buf,"Unable to access config file %s",confFile);
            }
            NAMD_die(buf);
          }
          char tmpcurdir[3];
          tmpcurdir[0] = '.';
          tmpcurdir[1] = PATHSEP;
          tmpcurdir[2] = 0;
          currentdir = tmpcurdir;
          iout << iINFO << "Working in the current directory " << origcwd << "\n" << endi;
      }
  }

#ifdef MEM_OPT_VERSION
    int dirlen = strlen(currentdir);
    gWorkDir = new char[dirlen+1];
    gWorkDir[dirlen]=0;
    memcpy(gWorkDir, currentdir, dirlen);
#endif

  currentdir = NULL;

#ifdef NAMD_TCL
  script->load(confFile);
#else
  script->run(confFile);
#endif

#ifdef NAMD_TCL
}
  script->run();
#endif

  BackEnd::exit();
}

