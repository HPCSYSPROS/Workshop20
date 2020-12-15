/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/


#ifdef WIN32
#ifndef NO_SOCKET
#define NO_SOCKET
#endif
#endif

#ifndef NO_SOCKET

// to track usage, define target host and port in compile line, e.g.:
//    -DTBSOFT_TRACK_HOST=\"127.0.0.1\" -DTBSOFT_TRACK_PORT=3141
// #define TBSOFT_TRACK_HOST   "127.0.0.1"
// #define TBSOFT_TRACK_PORT   3141

#ifndef TBSOFT_TRACK_MAXLEN
#define TBSOFT_TRACK_MAXLEN 1024              /* maximum message length */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <pwd.h>

#endif

#include "InfoStream.h"
#include "memusage.h"

#include "Lattice.h"
#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "main.decl.h"
#include "main.h"

#ifndef NO_SOCKET

int send_dgram(const char *host_addr, int port, const char *buf, int buflen) {
  struct sockaddr_in addr;
  int sockfd;

#ifndef NOHOSTNAME
  if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
    return -1;
  } 

  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  addr.sin_addr.s_addr = inet_addr(host_addr);

  sendto(sockfd, buf, buflen, 0, (struct sockaddr *)&addr, sizeof(addr));

  close(sockfd);
#endif

  return 0;
}                     


int tbsoft_sendusage(const char *program, 
                     const char *versionnum,
                     const char *platform,
                     const char *numcpus,
                     const char *miscinfo) {

#ifndef NOHOSTNAME
#ifdef TBSOFT_TRACK_HOST
  iout << iINFO
    << "Sending usage information to " << TBSOFT_TRACK_HOST
    << ":" << TBSOFT_TRACK_PORT << " via UDP.  Sent data is:\n";
#endif

  char sendbuf[TBSOFT_TRACK_MAXLEN];
  char host[128];
  struct passwd *pw;
  char user[128];

  memset(sendbuf, 0, sizeof(sendbuf));

  gethostname(host, 128);  host[127] = 0;
  pw = getpwuid(getuid());
  if ( pw && pw->pw_name ) {
    strncpy(user, pw->pw_name, 127);  user[127] = 0;
  } else {
    sprintf(user,"%d",getuid());
  }

  sprintf(sendbuf, "1 %s  %s  %s  %s  %s  %s  %s", 
    program, versionnum, platform, numcpus, miscinfo, host, user);
  iout << iINFO << sendbuf << "\n" << endi;
#ifdef TBSOFT_TRACK_HOST
  send_dgram(TBSOFT_TRACK_HOST, TBSOFT_TRACK_PORT, sendbuf, strlen(sendbuf));
#endif

#endif
  return 0;
}

#endif

extern const char *namd_build_date;
extern const char *namd_build_user;
extern const char *namd_build_machine;

class main : public Chare
{
public:
  main(CkArgMsg *)
  {

    // print banner
    iout << iINFO << "NAMD " << NAMD_VERSION << " for " << NAMD_PLATFORM
         << "\n"
#ifdef MEM_OPT_VERSION
         << iWARN << "\n"
         << iWARN << "       ***  EXPERIMENTAL MEMORY OPTIMIZED VERSION  ***\n"
         << iWARN << "\n"
#endif
#if 0
         << iWARN << "\n"
         << iWARN << "          ***  UNRELEASED EXPERIMENTAL VERSION  ***\n"
         << iWARN << "\n"
#endif
#ifdef SPEC_DISABLED_VERSION

         << iINFO << "\n"
         << iINFO << "NAMD is a parallel, object-oriented molecular dynamics\n"
         << iINFO << "code designed for high-performance simulation of large\n"
         << iINFO << "biomolecular systems.  NAMD is distributed free of\n"
         << iINFO << "charge and includes source code.  For more information\n" 
         << iINFO << "please visit http://www.ks.uiuc.edu/Research/namd/\n"
         << iINFO << "\n"
         << iINFO << "*********************************************************\n"
         << iINFO << "This version of NAMD may be distributed only as a part of\n"
         << iINFO << "the SPEC Workstation Benchmark and all other distribution\n"
         << iINFO << "is prohibited.  Any use of this software is bound by\n"
         << iINFO << "the terms of the NAMD License, which is available at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/license.html\n"
         << iINFO << "The NAMD development team will not provide support for\n"
         << iINFO << "any version of NAMD unless you have first registered\n"
         << iINFO << "and downloaded the latest version of NAMD available at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/\n"
         << iINFO << "*********************************************************\n"
#else
         << iINFO << "\n"
         << iINFO << "Please visit http://www.ks.uiuc.edu/Research/namd/\n"
         << iINFO << "for updates, documentation, and support information.\n"
#endif
<< iINFO << "\n"
<< iINFO << "Please cite Phillips et al., J. Comp. Chem. 26:1781-1802 (2005)\n"
<< iINFO << "in all publications reporting results obtained with NAMD.\n"
<< iINFO << "\n"
         << endi;

    char charm_version[64];
    sprintf(charm_version,"%d",CHARM_VERSION);

#if CHARM_VERSION < 60500
#error "Charm++ 6.5.1 or later is required to build NAMD"
#endif

    iout << iINFO << "Based on Charm++/Converse " << charm_version
         << " for " << CMK_MACHINE_NAME << "\n" << endi;

    iout << iINFO << "Built " << namd_build_date << " by "
         << namd_build_user << " on " << namd_build_machine << "\n"
         << endi;
#ifndef NO_SOCKET
    char numcpus[512];
    sprintf(numcpus,"%d",CkNumPes());
    tbsoft_sendusage("NAMD",NAMD_VERSION,NAMD_PLATFORM,numcpus,"");
#endif

#if CMK_BLUEGENE_CHARM
    iout << iINFO << "Running on BigSim using " << CmiNumPes() << " real processors.\n" << endi;
#endif
    iout << iINFO << "Running on " << CkNumPes() << " processors, "
         << CmiNumNodes() << " nodes, "
         << CmiNumPhysicalNodes() << " physical nodes.\n" << endi;
    iout << iINFO << "CPU topology information " << (CmiCpuTopologyEnabled()?"available":"unavailable") << ".\n" << endi;
    iout << iINFO << "Charm++/Converse parallel runtime startup completed at "
	 << CmiWallTimer() << " s\n"<< endi;
    const char* memsource;
    memusage(&memsource);
    iout << iINFO << memusage_MB() << " MB of memory in use"
	 << " based on " << memsource << "\n";

#if CMK_SMP
    if ( CmiNumNodes() > 1 && CkNumPes() == CmiNumNodes() ) {
      NAMD_die("SMP build launched as multiple single-thread processes.  Use ++ppn to set number of worker threads per process to match available cores, reserving one core per process for communication thread.");
    }
#endif
  }
};

#include "main.def.h"

