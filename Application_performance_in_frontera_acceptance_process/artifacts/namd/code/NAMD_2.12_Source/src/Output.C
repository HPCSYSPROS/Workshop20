/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   This object outputs the data collected on the master node
*/

#include "largefiles.h"  // must be first!

#include <string.h>
#include <stdlib.h>

#include "InfoStream.h"
#include "IMDOutput.h"
#include "Output.h"
#include "dcdlib.h"
#include "strlib.h"
#include "Molecule.h"
#include "Node.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Vector.h"
#include "structures.h"
#include "MStream.h"
#include "Communicate.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "ScriptTcl.h"
#include "Lattice.h"
#include "DataExchanger.h"
#include <fcntl.h>
#include <sys/stat.h>
#ifdef WIN32
#include <io.h>
#define access(PATH,MODE) _access(PATH,00)
#endif

#if defined(WIN32) && !defined(__CYGWIN__)
#define PATHSEPSTR "\\"
#define MKDIR(X) mkdir(X)
#else
#define PATHSEPSTR "/"
#define MKDIR(X) mkdir(X,0777)
#endif

#define NAMD_open NAMD_open64
#define NAMD_write NAMD_write64
#define NAMD_close NAMD_close64

#ifndef O_LARGEFILE
#define O_LARGEFILE 0x0
#endif

// same as open, only does error checking internally
int NAMD_open(const char *fname) {
  int fd;

  //  open the file and die if the open fails
#ifdef WIN32
  while ( (fd = _open(fname, O_WRONLY|O_CREAT|O_EXCL|O_BINARY|O_LARGEFILE,_S_IREAD|_S_IWRITE)) < 0) {
#else
#ifdef NAMD_NO_O_EXCL
  while ( (fd = open(fname, O_WRONLY|O_CREAT|O_TRUNC|O_LARGEFILE,
#else
  while ( (fd = open(fname, O_WRONLY|O_CREAT|O_EXCL|O_LARGEFILE,
#endif
                           S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0) {
#endif
    if ( errno != EINTR ) {
      char errmsg[1024];
      sprintf(errmsg, "Unable to open binary file %s", fname);
      NAMD_err(errmsg);
    }
  }

  return fd;
}

// same as write, only does error checking internally
void NAMD_write(int fd, const char *buf, size_t count, const char *errmsg="NAMD_write64") {
  double firsttime = 0.;
  while ( count ) {
#if defined(WIN32) && !defined(__CYGWIN__)
    long retval = _write(fd,buf,count);
#else
    ssize_t retval = write(fd,buf,count);
#endif
    if ( retval < 0 && errno == EINTR ) retval = 0;
    if ( retval < 0 && errno == ENOMEM ) {
      if ( firsttime == 0. ) firsttime = CmiWallTimer();
      if ( (CmiWallTimer() - firsttime) < 300. ) retval = 0;
    }
    if ( retval < 0 ) NAMD_err(errmsg);
    if ( retval > count ) NAMD_bug("extra bytes written in NAMD_write64()");
    buf += retval;
    count -= retval;
  }
  if ( firsttime != 0. ) {
    iout << iWARN << errmsg << ": NAMD_write64() retried for " << (CmiWallTimer() - firsttime) << " seconds.\n" << endi;
  }
}

// same as close, only does error checking internally
void NAMD_close(int fd, const char *fname) {
#ifdef WIN32
  while ( _close(fd) ) {
#else
  while ( close(fd) ) {
#endif
    if ( errno != EINTR ) {
      char errmsg[1024];
      sprintf(errmsg, "Error on closing file %s", fname);
      NAMD_err(errmsg);
    }
  }
}


#define seek_dcdfile NAMD_seek

// These make the NAMD 1 names work in NAMD 2
#define namdMyNode Node::Object()
#define simParams simParameters
#define pdbData pdb


static void lattice_to_unitcell(const Lattice *lattice, double *unitcell) {
   if (lattice && lattice->a_p() && lattice->b_p() && lattice->c_p()) {
      const Vector &a=lattice->a();
      const Vector &b=lattice->b();
      const Vector &c=lattice->c();
      unitcell[0] = a.length();
      unitcell[2] = b.length();
      unitcell[5] = c.length();
      double cosAB = (a*b)/(unitcell[0]*unitcell[2]);
      double cosAC = (a*c)/(unitcell[0]*unitcell[5]);
      double cosBC = (b*c)/(unitcell[2]*unitcell[5]);
      if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
      if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
      if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;
      unitcell[1] = cosAB;
      unitcell[3] = cosAC;
      unitcell[4] = cosBC;
   } else {
      unitcell[0] = unitcell[2] = unitcell[5] = 1.0;
      unitcell[1] = unitcell[3] = unitcell[4] = 0.0;
   }
}


/************************************************************************/
/*                  */
/*      FUNCTION Output          */
/*                  */
/*  This is the constructor for the Ouput class.  It just sets   */
/*  up some values for the VMD connection.        */
/*                  */
/************************************************************************/

Output::Output() : replicaDcdActive(0) { }

/*      END OF FUNCTION Output        */

/************************************************************************/
/*                  */
/*      FUNCTION ~Output        */
/*                  */
/************************************************************************/

Output::~Output() { }

/*      END OF FUNCTION ~Output        */

/************************************************************************/
/*                  */
/*      FUNCTION coordinate        */
/*                  */
/*   INPUTS:                */
/*  timestep - Timestep coordinates were accumulated for    */
/*  n - This is the number of coordinates accumulated.    */
/*  vel - Array of Vectors containing the velocities    */
/*                  */
/*  This function receives the coordinates accumulated for a given  */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output coordinates information    */
/*   should be called from here.          */
/*                  */
/************************************************************************/

int Output::coordinateNeeded(int timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if(simParams->benchTimestep) return 0;

  int positionsNeeded = 0;

  if ( timestep >= 0 ) {

    //  Output a DCD trajectory 
    if ( simParams->dcdFrequency &&
       ((timestep % simParams->dcdFrequency) == 0) )
    { positionsNeeded |= 1; }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    { positionsNeeded |= 2; }

    //  Iteractive MD
    if ( simParams->IMDon &&
       ( ((timestep % simParams->IMDfreq) == 0) ||
         (timestep == simParams->firstTimestep) ) )
      { positionsNeeded |= 1; }

  }

  //  Output final coordinates
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN ||
					timestep == EVAL_MEASURE)
  {
    positionsNeeded |= 2;
  }

  return positionsNeeded;
}

template <class xVector, class xDone>
void wrap_coor_int(xVector *coor, Lattice &lattice, xDone *done) {
  SimParameters *simParams = Node::Object()->simParameters;
  if ( *done ) return;
  *done = 1;
  if ( ! ( simParams->wrapAll || simParams->wrapWater ) ) return;
  const int wrapNearest = simParams->wrapNearest;
  const int wrapAll = simParams->wrapAll;
  Molecule *molecule = Node::Object()->molecule;
  int n = molecule->numAtoms;
  int i;
#ifndef MEM_OPT_VERSION
  Position *con = new Position[n];
  for ( i = 0; i < n; ++i ) {
    con[i] = 0;
    int ci = molecule->get_cluster(i);
    con[ci] += coor[i];
  }
  for ( i = 0; i < n; ++i ) {
    if ( ! wrapAll && ! molecule->is_water(i) ) continue;
    int ci = molecule->get_cluster(i);
    if ( ci == i ) {
      Vector coni = con[i] / molecule->get_clusterSize(i);
      Vector trans = ( wrapNearest ?
	lattice.wrap_nearest_delta(coni) : lattice.wrap_delta(coni) );
      con[i] = trans;
    }
    coor[i] = coor[i] + con[ci];
  }
  delete [] con;
#endif
}

void wrap_coor(Vector *coor, Lattice &lattice, double *done) {
  wrap_coor_int(coor,lattice,done);
};

void wrap_coor(FloatVector *coor, Lattice &lattice, float *done) {
  wrap_coor_int(coor,lattice,done);
};

void Output::coordinate(int timestep, int n, Vector *coor, FloatVector *fcoor,
							Lattice &lattice)
{
  SimParameters *simParams = Node::Object()->simParameters;
  double coor_wrapped = 0;
  float fcoor_wrapped = 0;

  if ( timestep >= 0 ) {

    //  Output a DCD trajectory 
    if ( simParams->dcdFrequency &&
       ((timestep % simParams->dcdFrequency) == 0) )
    {
      wrap_coor(fcoor,lattice,&fcoor_wrapped);
      output_dcdfile(timestep, n, fcoor, 
          simParams->dcdUnitCell ? &lattice : NULL);
    }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    {
      iout << "WRITING COORDINATES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
      wrap_coor(coor,lattice,&coor_wrapped);
      output_restart_coordinates(coor, n, timestep);
      iout << "FINISHED WRITING RESTART COORDINATES\n" <<endi;
      fflush(stdout);
    }

    //  Interactive MD
    if ( simParams->IMDon &&
       ( ((timestep % simParams->IMDfreq) == 0) ||
         (timestep == simParams->firstTimestep) ) )
    {
      IMDOutput *imd = Node::Object()->imd;
      wrap_coor(fcoor,lattice,&fcoor_wrapped);
      if (imd != NULL) imd->gather_coordinates(timestep, n, fcoor);
    }

  }

  if (timestep == EVAL_MEASURE)
  {
#ifdef NAMD_TCL
    wrap_coor(coor,lattice,&coor_wrapped);
    Node::Object()->getScript()->measure(coor);
#endif
  }

  //  Output final coordinates
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    int realstep = ( timestep == FILE_OUTPUT ?
          simParams->firstTimestep : simParams->N );
    iout << "WRITING COORDINATES TO OUTPUT FILE AT STEP "
				<< realstep << "\n" << endi;
    fflush(stdout);
    wrap_coor(coor,lattice,&coor_wrapped);
    output_final_coordinates(coor, n, realstep);
  }

  //  Close trajectory files
  if (timestep == END_OF_RUN)
  {
    if (simParams->dcdFrequency) output_dcdfile(END_OF_RUN,0,0, 
        simParams->dcdUnitCell ? &lattice : NULL);
  }

}
/*    END OF FUNCTION coordinate        */

/************************************************************************/
/*                  */
/*      FUNCTION velocity        */
/*                  */
/*   INPUTS:                */
/*  timestep - Timestep velocities were accumulated for    */
/*  n - This is the number of velocities accumulated.    */
/*  vel - Array of Vectors containing the velocities    */
/*                  */
/*  This function receives the velocities accumulated for a given   */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output velocity information should*/
/*   be called from here.            */
/*                  */
/************************************************************************/

int Output::velocityNeeded(int timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if(simParams->benchTimestep) return 0;

  int velocitiesNeeded = 0;

  if ( timestep >= 0 ) {

    //  Output a velocity DCD trajectory
    if ( simParams->velDcdFrequency &&
       ((timestep % simParams->velDcdFrequency) == 0) )
      { velocitiesNeeded |= 1; }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
      { velocitiesNeeded |= 2; }

  }

  //  Output final velocities
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    velocitiesNeeded |= 2;
  }

  return velocitiesNeeded;
}

void Output::velocity(int timestep, int n, Vector *vel)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if ( timestep >= 0 ) {

    //  Output velocity DCD trajectory
    if ( simParams->velDcdFrequency &&
       ((timestep % simParams->velDcdFrequency) == 0) )
    {
      output_veldcdfile(timestep, n, vel);
    }

  //  Output restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    {
      iout << "WRITING VELOCITIES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
      output_restart_velocities(timestep, n, vel);
      iout << "FINISHED WRITING RESTART VELOCITIES\n" <<endi;
      fflush(stdout);
    }

  }

  //  Output final velocities
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    int realstep = ( timestep == FILE_OUTPUT ?
          simParams->firstTimestep : simParams->N );
    iout << "WRITING VELOCITIES TO OUTPUT FILE AT STEP "
				<< realstep << "\n" << endi;
    fflush(stdout);
    output_final_velocities(realstep, n, vel);
  }

  //  Close trajectory files
  if (timestep == END_OF_RUN)
  {
    if (simParams->velDcdFrequency) output_veldcdfile(END_OF_RUN,0,0);
    // close force dcd file here since no final force output below
    if (simParams->forceDcdFrequency) output_forcedcdfile(END_OF_RUN,0,0);
  }

}
/*      END OF FUNCTION velocity      */

/************************************************************************/
/*                  */
/*      FUNCTION force */
/*                  */
/*   INPUTS:                */
/*  timestep - Timestep forces were accumulated for    */
/*  n - This is the number of forces accumulated.    */
/*  frc - Array of Vectors containing the forces */
/*                  */
/*  This function receives the forces accumulated for a given   */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output force information should*/
/*   be called from here.            */
/*                  */
/************************************************************************/

int Output::forceNeeded(int timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if(simParams->benchTimestep) return 0;

  int forcesNeeded = 0;

  if ( timestep >= 0 ) {

    //  Output a force DCD trajectory
    if ( simParams->forceDcdFrequency &&
       ((timestep % simParams->forceDcdFrequency) == 0) )
      { forcesNeeded |= 1; }

  }

  //  Output forces
  if (timestep == FORCE_OUTPUT)
  {
    forcesNeeded |= 2;
  }

  return forcesNeeded;
}

void Output::force(int timestep, int n, Vector *frc)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if ( timestep >= 0 ) {

    //  Output force DCD trajectory
    if ( simParams->forceDcdFrequency &&
       ((timestep % simParams->forceDcdFrequency) == 0) )
    {
      output_forcedcdfile(timestep, n, frc);
    }

  }

  //  Output forces
  if (timestep == FORCE_OUTPUT)
  {
    int realstep = simParams->firstTimestep;
    iout << "WRITING FORCES TO OUTPUT FILE AT STEP "
				<< realstep << "\n" << endi;
    fflush(stdout);
    output_forces(realstep, n, frc);
  }

  //  Trajectory file closed by velocity() above

}
/*      END OF FUNCTION force */

/************************************************************************/
/*                  */
/*      FUNCTION output_restart_coordinates    */
/*                  */
/*   INPUTS:                */
/*  coor - Array of vectors containing current positions    */
/*  n - Number of coordinates to output        */
/*  timestep - Timestep for which the coordinates are being written */
/*                  */
/*  This function writes out the current positions of all the atoms */
/*   in PDB format to the restart file specified by the user in the   */
/*   configuration file.            */
/*                  */
/************************************************************************/

void Output::output_restart_coordinates(Vector *coor, int n, int timestep)

{
  char comment[128];    //  Comment for header of PDB file
  char timestepstr[20];

  int baselen = strlen(namdMyNode->simParams->restartFilename);
  char *restart_name = new char[baselen+26];
  const char *bsuffix = ".old";

  strcpy(restart_name, namdMyNode->simParams->restartFilename);
  if ( namdMyNode->simParams->restartSave ) {
    sprintf(timestepstr,".%d",timestep);
    strcat(restart_name, timestepstr);
    bsuffix = ".BAK";
  }
  strcat(restart_name, ".coor");

  NAMD_backup_file(restart_name,bsuffix);

  //  Check to see if we should generate a binary or PDB file
  if (!namdMyNode->simParams->binaryRestart)
  {
    //  Generate a PDB restart file
    sprintf(comment, "RESTART COORDINATES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    namdMyNode->pdbData->set_all_positions(coor);
    namdMyNode->pdbData->write(restart_name, comment);
  }
  else
  {
    //  Generate a binary restart file
    write_binary_file(restart_name, n, coor);
  }

  delete [] restart_name;

  if ( namdMyNode->simParams->restartSaveDcd ) {
    if ( ! output_dcdfile(END_OF_RUN, 0, 0, 0) ) { // close old file
      const char *old_name = namdMyNode->simParams->dcdFilename;
      int old_len = strlen(old_name);
      char *new_name = new char[old_len+26];
      strcpy(new_name, old_name);
      if ( old_len >= 4 && ! strcmp(new_name+old_len-4,".dcd") ) {
        old_len -= 4;
        new_name[old_len] = 0;
      }
      sprintf(timestepstr,".%d",timestep);
      strcat(new_name, timestepstr);
      strcat(new_name, ".dcd");
      iout << "RENAMING COORDINATE DCD FILE " << old_name << " TO " << new_name << "\n" << endi;
      NAMD_backup_file(new_name,".BAK");
      while ( rename(old_name, new_name) ) {
        if ( errno == EINTR || errno == EXDEV ) continue;
        char err_msg[257];
        sprintf(err_msg, "Unable to rename DCD file %s to %s", old_name, new_name);
        NAMD_err(err_msg);
      }
      delete [] new_name;
    }
  }

}
/*      END OF FUNCTION output_restart_coordinates  */

/************************************************************************/
/*                  */
/*      FUNCTION output_restart_velocities    */
/*                  */
/*   INPUTS:                */
/*  vel - Array of vectors containing current velocities    */
/*  timestep - Timestep for which the velocities are being written  */
/*                  */
/*  This function writes out the current velocites of all the atoms */
/*   in PDB format to the restart file specified by the user in the   */
/*   configuration file.            */
/*                  */
/************************************************************************/

void Output::output_restart_velocities(int timestep, int n, Vector *vel)

{
  char comment[128];    //  comment for the header of PDB file
  char timestepstr[20];

  int baselen = strlen(namdMyNode->simParams->restartFilename);
  char *restart_name = new char[baselen+26];
  const char *bsuffix = ".old";

  strcpy(restart_name, namdMyNode->simParams->restartFilename);
  if ( namdMyNode->simParams->restartSave ) {
    sprintf(timestepstr,".%d",timestep);
    strcat(restart_name, timestepstr);
    bsuffix = ".BAK";
  }
  strcat(restart_name, ".vel");

  NAMD_backup_file(restart_name,bsuffix);

  //  Check to see if we should write out a PDB or a binary file
  if (!namdMyNode->simParams->binaryRestart)
  {
    //  Write the coordinates to a PDB file.  Multiple them by 20
    //  first to make the numbers bigger
    sprintf(comment, "RESTART VELOCITIES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    scale_vels(vel, n, PDBVELFACTOR);
    namdMyNode->pdbData->set_all_positions(vel);
    namdMyNode->pdbData->write(restart_name, comment);
    scale_vels(vel, n, PDBVELINVFACTOR);
  }
  else
  {
    //  Write the velocities to a binary file
    write_binary_file(restart_name, n, vel);
  }

  delete [] restart_name;
}
/*      END OF FUNCTION output_restart_velocities  */



// This is here so it can access Output.h and Node.h
void SimParameters::close_dcdfile() {

  Output *output = Node::Object()->output;
  if ( ! output ) return;

  output->output_dcdfile(END_OF_RUN, 0, 0, 0);

}

void Output::setReplicaDcdIndex(int index) {
  replicaDcdActive = 1;
  replicaDcdIndex = index;
}

void Output::replicaDcdInit(int index, const char *filename) {
  replicaDcdActive = 1;
  replicaDcdIndex = index;
  int msgsize = sizeof(ReplicaDcdInitMsg) + strlen(filename);
  ReplicaDcdInitMsg *msg = (ReplicaDcdInitMsg *) CmiAlloc(msgsize);
  msg->srcPart = CmiMyPartition();
  msg->dcdIndex = replicaDcdIndex;
  strcpy(msg->data, filename);
  sendReplicaDcdInit(abs(replicaDcdIndex) % CmiNumPartitions(), msg, msgsize);
}

void Output::recvReplicaDcdInit(ReplicaDcdInitMsg *msg) {
  replicaDcdFile &f = replicaDcdFiles[msg->dcdIndex];
  if ( f.fileid ) {
    iout << "CLOSING REPLICA DCD FILE " << msg->dcdIndex << " " << f.filename.c_str() << "\n" << endi;
    close_dcd_write(f.fileid);
    f.fileid = 0;
  }
  f.filename = (const char*) msg->data;
  sendReplicaDcdAck(msg->srcPart, (ReplicaDcdAckMsg*) CmiAlloc(sizeof(ReplicaDcdAckMsg)));
}

void Output::recvReplicaDcdData(ReplicaDcdDataMsg *msg) {
  if ( ! replicaDcdFiles.count(msg->dcdIndex) ) {
    char err_msg[257];
    sprintf(err_msg, "Unknown replicaDcdFile identifier %d\n", msg->dcdIndex);
    NAMD_die(err_msg);
  }
  replicaDcdFile &f = replicaDcdFiles[msg->dcdIndex];

  if ( ! f.fileid ) {
    //  Open the DCD file
    iout << "OPENING REPLICA DCD FILE " << msg->dcdIndex << " " << f.filename.c_str() << "\n" << endi;

    f.fileid=open_dcd_write(f.filename.c_str());

    if (f.fileid == DCD_FILEEXISTS) {
      char err_msg[257];
      sprintf(err_msg, "DCD file %s already exists!!", f.filename.c_str());
      NAMD_err(err_msg);
    } else if (f.fileid < 0) {
      char err_msg[257];
      sprintf(err_msg, "Couldn't open DCD file %s", f.filename.c_str());
      NAMD_err(err_msg);
    } else if (! f.fileid) {
      NAMD_bug("Output::recvReplicaDcdData open_dcd_write returned fileid of zero");
    }

    //  Write out the header
    int ret_code = write_dcdheader(f.fileid, f.filename.c_str(),
        msg->numAtoms, msg->NFILE, msg->NPRIV, msg->NSAVC, msg->NSTEP,
        msg->DELTA, msg->with_unitcell);

    if (ret_code<0) {
      NAMD_err("Writing of DCD header failed!!");
    }
  }

  //  Write out the values for this timestep
  iout << "WRITING TO REPLICA DCD FILE " << msg->dcdIndex << " " << f.filename.c_str() << "\n" << endi;
  float *msgx = (float*) msg->data;
  float *msgy = msgx + msg->numAtoms;
  float *msgz = msgy + msg->numAtoms;
  int ret_code = write_dcdstep(f.fileid, msg->numAtoms, msgx, msgy, msgz,
                                   msg->with_unitcell ? msg->unitcell : 0);
  if (ret_code < 0) NAMD_err("Writing of DCD step failed!!");

  sendReplicaDcdAck(msg->srcPart, (ReplicaDcdAckMsg*) CmiAlloc(sizeof(ReplicaDcdAckMsg)));
}


/************************************************************************/
/*                  */
/*      FUNCTION output_dcdfile        */
/*                  */
/*   INPUTS:                */
/*  timestep - Current timestep          */
/*  n - Number of atoms in simulation        */
/*  coor - Coordinate vectors for all atoms        */
/*  lattice - periodic cell data; NULL if not to be written */
/*                  */
/*  This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.      */
/*                  */
/************************************************************************/

#define RAD2DEG 180.0/3.14159265359

int Output::output_dcdfile(int timestep, int n, FloatVector *coor,
    const Lattice *lattice)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file

  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  static int n_alloc;  // allocated size
  
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  SimParameters *simParams = namdMyNode->simParams;

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( timestep == END_OF_RUN ) {
    for ( std::map<int,replicaDcdFile>::iterator it = replicaDcdFiles.begin();
          it != replicaDcdFiles.end(); ++it ) {
      replicaDcdFile &f = it->second;
      if ( f.fileid ) {
        iout << "CLOSING REPLICA DCD FILE " << it->first << " " << f.filename.c_str() << "\n" << endi;
        close_dcd_write(f.fileid);
        f.fileid = 0;
      }
    }
    int rval = 0;
    if ( ! first ) {
      iout << "CLOSING COORDINATE DCD FILE " << simParams->dcdFilename << "\n" << endi;
      close_dcd_write(fileid);
    } else {
      iout << "COORDINATE DCD FILE " << simParams->dcdFilename << " WAS NOT CREATED\n" << endi;
      rval = -1;
    }
    first = 1;
    fileid = 0;
    return rval;
  }

  if ( replicaDcdActive ) {
    int msgsize = sizeof(ReplicaDcdDataMsg) + 3*n*sizeof(float);
    ReplicaDcdDataMsg *msg = (ReplicaDcdDataMsg *) CmiAlloc(msgsize);
    float *msgx = (float*) msg->data;
    float *msgy = msgx + n;
    float *msgz = msgy + n;
    for (i=0; i<n; i++) { msgx[i] = coor[i].x; }
    for (i=0; i<n; i++) { msgy[i] = coor[i].y; }
    for (i=0; i<n; i++) { msgz[i] = coor[i].z; }
    msg->numAtoms = n;
    lattice_to_unitcell(lattice,msg->unitcell);
    msg->with_unitcell = lattice ? 1 : 0;
    msg->NSAVC = simParams->dcdFrequency;
    msg->NPRIV = timestep;
    msg->NSTEP = msg->NPRIV - msg->NSAVC;
    msg->NFILE = 0;
    msg->DELTA = simParams->dt/TIMEFACTOR;
    msg->srcPart = CmiMyPartition();
    msg->dcdIndex = replicaDcdIndex;
    sendReplicaDcdData(abs(replicaDcdIndex) % CmiNumPartitions(), msg, msgsize);
    return 0;
  }

  if (first)
  {
    //  Allocate x, y, and z arrays since the DCD file routines
    //  need them passed as three independant arrays to be
    //  efficient
    if ( n > n_alloc ) {
      delete [] x;  x = new float[3*n];
      y = x + n;
      z = x + 2*n;
      n_alloc = n;
    }

    //  Open the DCD file
    iout << "OPENING COORDINATE DCD FILE\n" << endi;

    fileid=open_dcd_write(simParams->dcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "DCD file %s already exists!!",
        simParams->dcdFilename);

      NAMD_err(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open DCD file %s",
        simParams->dcdFilename);

      NAMD_err(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->dcdFrequency;
    NPRIV = timestep;
    NSTEP = NPRIV - NSAVC;
    NFILE = 0;

    //  Write out the header
    ret_code = write_dcdheader(fileid, 
        simParams->dcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR, lattice != NULL);


    if (ret_code<0)
    {
      NAMD_err("Writing of DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = coor[i].x;
    y[i] = coor[i].y;
    z[i] = coor[i].z;
  }

  //  Write out the values for this timestep
  iout << "WRITING COORDINATES TO DCD FILE " << simParams->dcdFilename << " AT STEP "
	<< timestep << "\n" << endi;
  fflush(stdout);
  if (lattice) {
    double unitcell[6];
    lattice_to_unitcell(lattice,unitcell);
    ret_code = write_dcdstep(fileid, n, x, y, z, unitcell);
  } else {
    ret_code = write_dcdstep(fileid, n, x, y, z, NULL);
  }
  if (ret_code < 0)
  {
    NAMD_err("Writing of DCD step failed!!");
  }

  return 0;
}
/*      END OF FUNCTION output_dcdfile      */

/************************************************************************/
/*                  */
/*      FUNCTION output_final_coordinates    */
/*                  */
/*   INPUTS:                */
/*  coor - Array of vectors containing final coordinates    */
/*  n - Number of coordinates to output        */
/*  timestep - Timestep that coordinates are being written in  */
/*                  */
/*  This function writes out the final coordinates for the    */
/*   simulation in PDB format to the file specified in the config  */
/*   file.                */
/*                  */
/************************************************************************/

void Output::output_final_coordinates(Vector *coor, int n, int timestep)

{
  char output_name[140];  //  Output filename
  char comment[128];    //  comment for PDB header

  //  Built the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".coor");

  NAMD_backup_file(output_name);

  //  Check to see if we should write out a binary file or a
  //  PDB file
  if (!namdMyNode->simParams->binaryOutput)
  {
    sprintf(comment, "FINAL COORDINATES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    namdMyNode->pdbData->set_all_positions(coor);
    namdMyNode->pdbData->write(output_name, comment);
  }
  else
  {
    //  Write the velocities to a binary file
    write_binary_file(output_name, n, coor);
  }
}
/*    END OF FUNCTION output_final_coordinates    */

/************************************************************************/
/*                  */
/*      FUNCTION output_final_velocities    */
/*                  */
/*   INPUTS:                */
/*  vel - Array of vectors containing final velocities    */
/*  timestep - Timestep that vleocities are being written in  */
/*                  */
/*  This function writes out the final vleocities for the    */
/*   simulation in PDB format to the file specified in the config  */
/*   file.                */
/*                  */
/************************************************************************/

void Output::output_final_velocities(int timestep, int n, Vector *vel)

{
  char output_name[140];  //  Output filename
  char comment[128];    //  Comment for PDB header

  //  Build the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".vel");

  NAMD_backup_file(output_name);

  //  Check to see if we should write a PDB or binary file
  if (!(namdMyNode->simParams->binaryOutput))
  {
    //  Write the final velocities to a PDB file
    sprintf(comment, "FINAL VELOCITIES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    scale_vels(vel, n, PDBVELFACTOR);
    namdMyNode->pdbData->set_all_positions(vel);
    namdMyNode->pdbData->write(output_name, comment);
    scale_vels(vel, n, PDBVELINVFACTOR);
  }
  else
  {
    //  Write the coordinates to a binary file
    write_binary_file(output_name, n, vel);
  }

}
/*      END OF FUNCTION output_final_velocities    */

/************************************************************************/
/*                  */
/*      FUNCTION output_veldcdfile      */
/*                  */
/*   INPUTS:                */
/*  timestep - Current timestep          */
/*  n - Number of atoms in simulation        */
/*  coor - velocity vectors for all atoms        */
/*                  */
/*  This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.  This fucntion   */
/*   writes out the velocity vectors in DCD format.      */
/*                  */
/************************************************************************/

void Output::output_veldcdfile(int timestep, int n, Vector *vel)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  static int n_alloc;  // allocated size
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  SimParameters *simParams = Node::Object()->simParameters;

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( timestep == END_OF_RUN ) {
    if ( ! first ) {
      iout << "CLOSING VELOCITY DCD FILE\n" << endi;
      close_dcd_write(fileid);
    } else {
      iout << "VELOCITY DCD FILE WAS NOT CREATED\n" << endi;
    }
    return;
  }

  if (first)
  {
    //  Allocate x, y, and z arrays since the DCD file routines
    //  need them passed as three independant arrays to be
    //  efficient
    if ( n > n_alloc ) {
      delete [] x;  x = new float[3*n];
      y = x + n;
      z = x + 2*n;
      n_alloc = n;
    }

    //  Open the DCD file
    iout << "OPENING VELOCITY DCD FILE\n" << endi;

    fileid=open_dcd_write(namdMyNode->simParams->velDcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Velocity DCD file %s already exists!!",
        namdMyNode->simParams->velDcdFilename);

      NAMD_err(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open velocity DCD file %s",
        namdMyNode->simParams->velDcdFilename);

      NAMD_err(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->velDcdFrequency;
    NPRIV = timestep;
    NSTEP = NPRIV - NSAVC;
    NFILE = 0;

    //  Write out the header
    const int with_unitcell = 0;
    ret_code = write_dcdheader(fileid, 
        simParams->velDcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR, with_unitcell);


    if (ret_code<0)
    {
      NAMD_err("Writing of velocity DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = vel[i].x;
    y[i] = vel[i].y;
    z[i] = vel[i].z;
  }

  //  Write out the values for this timestep
  iout << "WRITING VELOCITIES TO DCD FILE AT STEP "
	<< timestep << "\n" << endi;
  fflush(stdout);
  ret_code = write_dcdstep(fileid, n, x, y, z, NULL);

  if (ret_code < 0)
  {
    NAMD_err("Writing of velocity DCD step failed!!");
  }

}
/*      END OF FUNCTION output_veldcdfile    */

/************************************************************************/
/*                  */
/*      FUNCTION output_forces    */
/*                  */
/*   INPUTS:                */
/*  frc - Array of vectors containing final forces */
/*  timestep - Timestep that vleocities are being written in  */
/*                  */
/*  This function writes out the final vleocities for the    */
/*   simulation in PDB format to the file specified in the config  */
/*   file.                */
/*                  */
/************************************************************************/

void Output::output_forces(int timestep, int n, Vector *frc)

{
  char output_name[140];  //  Output filename
  char comment[128];    //  Comment for PDB header

  //  Build the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".force");

  NAMD_backup_file(output_name);

  //  Check to see if we should write a PDB or binary file
  if (!(namdMyNode->simParams->binaryOutput))
  {
    //  Write the forces to a PDB file
    sprintf(comment, "FORCES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    namdMyNode->pdbData->set_all_positions(frc);
    namdMyNode->pdbData->write(output_name, comment);
  }
  else
  {
    //  Write the coordinates to a binary file
    write_binary_file(output_name, n, frc);
  }

}
/*      END OF FUNCTION output_forces */

/************************************************************************/
/*                  */
/*      FUNCTION output_forcedcdfile      */
/*                  */
/*   INPUTS:                */
/*  timestep - Current timestep          */
/*  n - Number of atoms in simulation        */
/*  frc - force vectors for all atoms        */
/*                  */
/*  This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.  This fucntion   */
/*   writes out the force vectors in DCD format.      */
/*                  */
/************************************************************************/

void Output::output_forcedcdfile(int timestep, int n, Vector *frc)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  static int n_alloc;  // allocated size
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  SimParameters *simParams = Node::Object()->simParameters;

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( timestep == END_OF_RUN ) {
    if ( ! first ) {
      iout << "CLOSING FORCE DCD FILE\n" << endi;
      close_dcd_write(fileid);
    } else {
      iout << "FORCE DCD FILE WAS NOT CREATED\n" << endi;
    }
    return;
  }

  if (first)
  {
    //  Allocate x, y, and z arrays since the DCD file routines
    //  need them passed as three independant arrays to be
    //  efficient
    if ( n > n_alloc ) {
      delete [] x;  x = new float[3*n];
      y = x + n;
      z = x + 2*n;
      n_alloc = n;
    }

    //  Open the DCD file
    iout << "OPENING FORCE DCD FILE\n" << endi;

    fileid=open_dcd_write(namdMyNode->simParams->forceDcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Force DCD file %s already exists!!",
        namdMyNode->simParams->forceDcdFilename);

      NAMD_err(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open force DCD file %s",
        namdMyNode->simParams->forceDcdFilename);

      NAMD_err(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->forceDcdFrequency;
    NPRIV = timestep;
    NSTEP = NPRIV - NSAVC;
    NFILE = 0;

    //  Write out the header
    const int with_unitcell = 0;
    ret_code = write_dcdheader(fileid, 
        simParams->forceDcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR, with_unitcell);


    if (ret_code<0)
    {
      NAMD_err("Writing of force DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = frc[i].x;
    y[i] = frc[i].y;
    z[i] = frc[i].z;
  }

  //  Write out the values for this timestep
  iout << "WRITING FORCES TO DCD FILE AT STEP "
	<< timestep << "\n" << endi;
  fflush(stdout);
  ret_code = write_dcdstep(fileid, n, x, y, z, NULL);

  if (ret_code < 0)
  {
    NAMD_err("Writing of force DCD step failed!!");
  }

}
/*      END OF FUNCTION output_forcedcdfile    */

/************************************************************************/
/*                  */
/*      FUNCTION write_binary_file      */
/*                  */
/*   INPUTS:                */
/*  fname - file name to write velocities to      */
/*  n - Number of atoms in system          */
/*  vels - Array of vectors            */
/*                  */
/*  This function writes out vectors in binary format to    */
/*   the specified file.            */
/*                  */
/************************************************************************/

void Output::write_binary_file(char *fname, int n, Vector *vecs)

{
  char errmsg[256];
  int fd;    //  File descriptor
  int32 n32 = n;

  fd = NAMD_open(fname);

  sprintf(errmsg, "Error on write to binary file %s", fname);

  //  Write out the number of atoms and the vectors
  NAMD_write(fd, (char *) &n32, sizeof(int32), errmsg);
  NAMD_write(fd, (char *) vecs, sizeof(Vector)*n, errmsg);

  NAMD_close(fd, fname);
}
/*      END OF FUNCTION write_binary_file    */

/************************************************************************/
/*                  */
/*      FUNCTION scale_vels        */
/*                  */
/*   INPUTS:                */
/*  v - Array of velocity vectors          */
/*  n - Number of atoms in system          */
/*  fact - Scaling factor            */
/*                  */
/*  This function scales all the vectors passed in by a constant  */
/*   factor.  This is used before writing out velocity vectors that  */
/*   need to be resized.            */
/*                  */
/************************************************************************/

void Output::scale_vels(Vector *v, int n, Real fact)

{
  int i;

  for (i=0; i<n; i++)
  {
    v[i].x *= fact;
    v[i].y *= fact;
    v[i].z *= fact;
  }
}
/*      END OF FUNCTION scale_vels      */


#ifdef MEM_OPT_VERSION
//////Beginning of Functions related with velocity output//////
void ParOutput::velocityMaster(int timestep, int n){
    SimParameters *simParams = Node::Object()->simParameters;

    if ( timestep >= 0 ) {

      //  Output velocity DCD trajectory
      if ( simParams->velDcdFrequency &&
         ((timestep % simParams->velDcdFrequency) == 0) )
      {         
        output_veldcdfile_master(timestep, n);        
      }

    //  Output restart file
      if ( simParams->restartFrequency &&
         ((timestep % simParams->restartFrequency) == 0) )
      {
        iout << "WRITING VELOCITIES TO RESTART FILE AT STEP "
                  << timestep << "\n" << endi;
        output_restart_velocities_master(timestep, n);
        iout << "FINISHED WRITING RESTART VELOCITIES\n" <<endi;
        fflush(stdout);
      }

    }

    //  Output final velocities
    if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
    {
      int realstep = ( timestep == FILE_OUTPUT ?
          simParams->firstTimestep : simParams->N );
      iout << "WRITING VELOCITIES TO OUTPUT FILE AT STEP "
                  << realstep << "\n" << endi;
      fflush(stdout);
      output_final_velocities_master(n);
    }

    //  Close trajectory files
    if (timestep == END_OF_RUN)
    {
      if (simParams->velDcdFrequency) output_veldcdfile_master(END_OF_RUN,0);
      if (simParams->forceDcdFrequency) output_forcedcdfile_master(END_OF_RUN,0);
    }

}

//output atoms' velocities from id fID to tID.
void ParOutput::velocitySlave(int timestep, int fID, int tID, Vector *vecs){
    SimParameters *simParams = Node::Object()->simParameters;

    if ( timestep >= 0 ) {

      //  Output velocity DCD trajectory
      if ( simParams->velDcdFrequency &&
         ((timestep % simParams->velDcdFrequency) == 0) )
      {         
        output_veldcdfile_slave(timestep, fID, tID, vecs);
      }

    //  Output restart file
      if ( simParams->restartFrequency &&
         ((timestep % simParams->restartFrequency) == 0) )
      {          
          int64 offset = sizeof(int)+sizeof(Vector)*((int64)fID);
          output_restart_velocities_slave(timestep, fID, tID, vecs, offset);     
      }

    }

    //  Output final velocities
    if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
    {
        int64 offset = sizeof(int)+sizeof(Vector)*((int64)fID);
        output_final_velocities_slave(fID, tID, vecs, offset);
    }

    //  Close trajectory files
    if (timestep == END_OF_RUN)
    {
      if (simParams->velDcdFrequency) output_veldcdfile_slave(END_OF_RUN, 0, 0, NULL);
      if (simParams->forceDcdFrequency) output_forcedcdfile_slave(END_OF_RUN, 0, 0, NULL);
    }
}

void ParOutput::output_veldcdfile_master(int timestep, int n){
    int ret_code;    //  Return code from DCD calls
    SimParameters *simParams = Node::Object()->simParameters;

    //  If this is the last time we will be writing coordinates,
    //  close the file before exiting
    if ( timestep == END_OF_RUN ) {
      if ( ! veldcdFirst ) {
        iout << "CLOSING VELOCITY DCD FILE\n" << endi;
        close_dcd_write(veldcdFileID);
      } else {
        iout << "VELOCITY DCD FILE WAS NOT CREATED\n" << endi;
      }
      return;      
    }

    if (veldcdFirst)
    {
      //  Open the DCD file
      iout << "OPENING VELOCITY DCD FILE\n" << endi;

#ifndef OUTPUT_SINGLE_FILE
#error OUTPUT_SINGLE_FILE not defined!
#endif
	  
	#if OUTPUT_SINGLE_FILE
	  char *veldcdFilename = simParams->velDcdFilename;
	#else	  
	  char *veldcdFilename = buildFileName(veldcdType);
	#endif

      veldcdFileID=open_dcd_write(veldcdFilename);

      if (veldcdFileID == DCD_FILEEXISTS)
      {
        char err_msg[257];
        sprintf(err_msg, "Velocity DCD file %s already exists!",veldcdFilename);
        NAMD_err(err_msg);
      }
      else if (veldcdFileID < 0)
      {
        char err_msg[257];
        sprintf(err_msg, "Couldn't open velocity DCD file %s",veldcdFilename);
        NAMD_err(err_msg);
      }

	#if !OUTPUT_SINGLE_FILE
	  // Write out extra fields as MAGIC number etc.
	  int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	  float tmpFlt = OUTPUT_FILE_VERSION;
	  NAMD_write(veldcdFileID, (char *) &tmpInt, sizeof(int32));
	  NAMD_write(veldcdFileID, (char *) &tmpFlt, sizeof(float));
	  tmpInt = simParams->numoutputprocs;
	  NAMD_write(veldcdFileID, (char *) &tmpInt, sizeof(int32));
	#endif

      int NSAVC, NFILE, NPRIV, NSTEP;
      NSAVC = simParams->velDcdFrequency;
      NPRIV = timestep;
      NSTEP = NPRIV - NSAVC;
      NFILE = 0;

      //  Write out the header
      const int with_unitcell = 0;
      ret_code = write_dcdheader(veldcdFileID, 
          veldcdFilename,
          n, NFILE, NPRIV, NSAVC, NSTEP,
          simParams->dt/TIMEFACTOR, with_unitcell);


      if (ret_code<0)
      {
        NAMD_err("Writing of velocity DCD header failed!!");
      }

	#if !OUTPUT_SINGLE_FILE
	  //In this case, the file name is dynamically allocated
	  delete [] veldcdFilename;
	#endif

      veldcdFirst = FALSE;
    }

    //  Write out the values for this timestep
    iout << "WRITING VELOCITIES TO DCD FILE AT STEP "
      << timestep << "\n" << endi;
    fflush(stdout);

	//In the case of writing to multiple files, only the header
	//of the dcd file needs to be updated. Note that the format of
	//the new dcd file has also changed! -Chao Mei 

#if OUTPUT_SINGLE_FILE
    //write X,Y,Z headers
    int totalAtoms = namdMyNode->molecule->numAtoms;
    write_dcdstep_par_XYZUnits(veldcdFileID, totalAtoms);
#endif

    //update the header
    update_dcdstep_par_header(veldcdFileID);    

}

void ParOutput::output_veldcdfile_slave(int timestep, int fID, int tID, Vector *vecs){
    int ret_code;    //  Return code from DCD calls
    SimParameters *simParams = Node::Object()->simParameters;

    //  If this is the last time we will be writing coordinates,
    //  close the file before exiting
    if ( timestep == END_OF_RUN ) {
      if ( ! veldcdFirst ) {        
        close_dcd_write(veldcdFileID);
      }
#if OUTPUT_SINGLE_FILE
      delete [] veldcdX;
      delete [] veldcdY;
      delete [] veldcdZ;
#endif
      return;
    }

    int parN = tID-fID+1;

    if (veldcdFirst)
    {

	#if OUTPUT_SINGLE_FILE
	  char *veldcdFilename = namdMyNode->simParams->velDcdFilename;
	#else	  
	  char *veldcdFilename = buildFileName(veldcdType);
	#endif
      veldcdFileID=open_dcd_write_par_slave(veldcdFilename);
      if(veldcdFileID < 0)
      {
        char err_msg[257];
        sprintf(err_msg, "Couldn't open velocity DCD file %s",veldcdFilename);
        NAMD_err(err_msg);
      }
	#if OUTPUT_SINGLE_FILE
	  //If outputting to a single file, dcd files conforms to the old format
	  //as data are organized as 3 seperate arrays of X,Y,Z, while in the new
	  //format used in outputing multiple files, the data are organized as an
	  //array of XYZs.
      veldcdX = new float[parN];
      veldcdY = new float[parN];
      veldcdZ = new float[parN];
	  //seek to beginning of X,Y,Z sections which means skipping header. 
	  //Cell data is not needed because this is velocity trajectory
	  int skipbytes = get_dcdheader_size();
	  seek_dcdfile(veldcdFileID, skipbytes, SEEK_SET);
	#endif

	#if !OUTPUT_SINGLE_FILE
	  delete [] veldcdFilename;
	  // Write out extra fields as MAGIC number etc.
	  int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	  float tmpFlt = OUTPUT_FILE_VERSION;
	  NAMD_write(veldcdFileID, (char *) &tmpInt, sizeof(int32));
	  NAMD_write(veldcdFileID, (char *) &tmpFlt, sizeof(float));	  
	  NAMD_write(veldcdFileID, (char *) &outputID, sizeof(int));
	  NAMD_write(veldcdFileID, (char *) &fID, sizeof(int));
	  NAMD_write(veldcdFileID, (char *) &tID, sizeof(int));
	#endif
	  
      veldcdFirst = FALSE;
    }

#if OUTPUT_SINGLE_FILE
    //The following seek will set the stream position to the
    //beginning of the place where a new timestep output should
    //be performed.
    CmiAssert(sizeof(off_t)==8);
    int totalAtoms = namdMyNode->molecule->numAtoms;

    for(int i=0; i<parN; i++){
        veldcdX[i] = vecs[i].x;
        veldcdY[i] = vecs[i].y;
        veldcdZ[i] = vecs[i].z;
    }

    write_dcdstep_par_slave(veldcdFileID, fID, tID, totalAtoms, veldcdX, veldcdY, veldcdZ);

	//same with the slave output for coordiantes trajectory file
	//but cell data is not needed because this is velocity trajectory
	int atomsRemains = (totalAtoms-1)-(tID+1)+1;
	off_t offset = ((off_t)atomsRemains)*sizeof(float)+1*sizeof(int);
	seek_dcdfile(veldcdFileID, offset, SEEK_CUR);

#else
	//write the timestep
	NAMD_write(veldcdFileID, (char *)&timestep, sizeof(int));
	//write the values for this timestep
	NAMD_write(veldcdFileID, (char *)vecs, sizeof(Vector)*parN);
#endif
}

void ParOutput::output_restart_velocities_master(int timestep, int n){
#if OUTPUT_SINGLE_FILE
	char timestepstr[20];

    int baselen = strlen(namdMyNode->simParams->restartFilename);
    char *restart_name = new char[baselen+26];

    strcpy(restart_name, namdMyNode->simParams->restartFilename);
    if ( namdMyNode->simParams->restartSave ) {
      sprintf(timestepstr,".%d",timestep);
      strcat(restart_name, timestepstr);
    }
    strcat(restart_name, ".vel");
#else
	char *restart_name = NULL;
	if ( namdMyNode->simParams->restartSave )
		restart_name = buildFileName(velType);
	else
		restart_name = buildFileName(velType,timestep);	
#endif

    NAMD_backup_file(restart_name,".old");

    //Always output a binary file
    write_binary_file_master(restart_name, n);

    delete [] restart_name;
}

void ParOutput::output_restart_velocities_slave(int timestep, int fID, int tID, Vector *vecs, int64 offset){    
#if OUTPUT_SINGLE_FILE
    char timestepstr[20];

    int baselen = strlen(namdMyNode->simParams->restartFilename);
    char *restart_name = new char[baselen+26];

    strcpy(restart_name, namdMyNode->simParams->restartFilename);
    if ( namdMyNode->simParams->restartSave ) {
      sprintf(timestepstr,".%d",timestep);
      strcat(restart_name, timestepstr);
    }
    strcat(restart_name, ".vel");    
#else
	char *restart_name = NULL;
	if ( namdMyNode->simParams->restartSave )
		restart_name = buildFileName(velType);
	else
		restart_name = buildFileName(velType,timestep);	

	NAMD_backup_file(restart_name,".old");
#endif

    //Always output a binary file	
    write_binary_file_slave(restart_name, fID, tID, vecs, offset);

    delete [] restart_name;
}

void ParOutput::output_final_velocities_master(int n){
#if OUTPUT_SINGLE_FILE
    char *output_name = new char[strlen(namdMyNode->simParams->outputFilename)+8];
    //  Build the output filename
    strcpy(output_name, namdMyNode->simParams->outputFilename);
    strcat(output_name, ".vel");
#else	
	char *output_name = buildFileName(velType);
#endif

    NAMD_backup_file(output_name);

    //Write the velocities to a binary file
    write_binary_file_master(output_name, n);
}

void ParOutput::output_final_velocities_slave(int fID, int tID, Vector *vecs, int64 offset){
#if OUTPUT_SINGLE_FILE
    char *output_name = new char[strlen(namdMyNode->simParams->outputFilename)+8];
    //  Build the output filename
    strcpy(output_name, namdMyNode->simParams->outputFilename);
    strcat(output_name, ".vel");
#else
	char *output_name = buildFileName(velType);
	NAMD_backup_file(output_name);
#endif
    
    //Write the velocities to a binary file
    write_binary_file_slave(output_name, fID, tID, vecs, offset);

	delete [] output_name;
}

void ParOutput::write_binary_file_master(char *fname, int n){
    char errmsg[256];
    int fd;    //  File descriptor
    int32 n32 = n;

    fd = NAMD_open(fname);

    sprintf(errmsg, "Error on write to binary file %s", fname);

  #if !OUTPUT_SINGLE_FILE
	// Write out extra fields as MAGIC number etc.
	int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	float tmpFlt = OUTPUT_FILE_VERSION;
	NAMD_write(fd, (char *) &tmpInt, sizeof(int32), errmsg);
	NAMD_write(fd, (char *) &tmpFlt, sizeof(float), errmsg);
	tmpInt = namdMyNode->simParams->numoutputprocs;
	NAMD_write(fd, (char *) &tmpInt, sizeof(int32), errmsg);
  #endif

    //  Write out the number of atoms and the vectors
    NAMD_write(fd, (char *) &n32, sizeof(int32), errmsg);

    NAMD_close(fd, fname);
}

void ParOutput::write_binary_file_slave(char *fname, int fID, int tID, Vector *vecs, int64 offset){
    char errmsg[256];

#if OUTPUT_SINGLE_FILE
	//the mode has to be "r+" because the file already exists
	FILE *ofp = fopen(fname, "rb+");
        if ( ! ofp ) {
          sprintf(errmsg, "Error on opening binary file %s", fname);
          NAMD_err(errmsg);
        }

	//if the output is a single file, then the file position needs to be set correctly
#ifdef WIN32
        if ( _fseeki64(ofp, offset, SEEK_SET) )
#else
        if ( fseeko(ofp, offset, SEEK_SET) )
#endif
        {
          sprintf(errmsg, "Error on seeking binary file %s", fname);
          NAMD_err(errmsg);
        }
#else
	//the mode has to be "w+" because the file doesn't exist yet
	FILE *ofp = fopen(fname, "wb+"); 
        if ( ! ofp ) {
          sprintf(errmsg, "Error on opening binary file %s", fname);
          NAMD_err(errmsg);
        }

	// Write out extra fields as MAGIC number etc.
	int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	float tmpFlt = OUTPUT_FILE_VERSION;	
	fwrite(&tmpInt, sizeof(int32), 1, ofp);	
	fwrite(&tmpFlt, sizeof(float), 1, ofp);
	fwrite(&outputID, sizeof(int), 1, ofp);
	fwrite(&fID, sizeof(int), 1, ofp);
	fwrite(&tID, sizeof(int), 1, ofp);
#endif

	int parN = tID-fID+1;
        if ( fwrite(vecs, sizeof(Vector), parN, ofp) != parN ) {
          sprintf(errmsg, "Error on writing to binary file %s", fname);
          NAMD_err(errmsg);
        }

	if ( fclose(ofp) ) {
          sprintf(errmsg, "Error on closing binary file %s", fname);
          NAMD_err(errmsg);
        }
}
//////End of Functions related with velocity output//////


//////Beginning of Functions related with force output//////
void ParOutput::forceMaster(int timestep, int n){
    SimParameters *simParams = Node::Object()->simParameters;

    if ( timestep >= 0 ) {

      //  Output force DCD trajectory
      if ( simParams->forceDcdFrequency &&
         ((timestep % simParams->forceDcdFrequency) == 0) )
      {         
        output_forcedcdfile_master(timestep, n);        
      }

    }

    //  Output forces
    if (timestep == FORCE_OUTPUT)
    {
      int realstep = simParams->firstTimestep;
      iout << "WRITING FORCES TO OUTPUT FILE AT STEP "
                  << realstep << "\n" << endi;
      fflush(stdout);
      output_forces_master(n);
    }

    //  Close trajectory files in velocityMaster above
}

//output atoms' forces from id fID to tID.
void ParOutput::forceSlave(int timestep, int fID, int tID, Vector *vecs){
    SimParameters *simParams = Node::Object()->simParameters;

    if ( timestep >= 0 ) {

      //  Output force DCD trajectory
      if ( simParams->forceDcdFrequency &&
         ((timestep % simParams->forceDcdFrequency) == 0) )
      {         
        output_forcedcdfile_slave(timestep, fID, tID, vecs);
      }

    }

    //  Output forces
    if (timestep == FORCE_OUTPUT)
    {
        int64 offset = sizeof(int)+sizeof(Vector)*((int64)fID);
        output_forces_slave(fID, tID, vecs, offset);
    }

    //  Close trajectory files in velocitySlave above
}

void ParOutput::output_forcedcdfile_master(int timestep, int n){
    int ret_code;    //  Return code from DCD calls
    SimParameters *simParams = Node::Object()->simParameters;

    //  If this is the last time we will be writing coordinates,
    //  close the file before exiting
    if ( timestep == END_OF_RUN ) {
      if ( ! forcedcdFirst ) {
        iout << "CLOSING FORCE DCD FILE\n" << endi;
        close_dcd_write(forcedcdFileID);
      } else {
        iout << "FORCE DCD FILE WAS NOT CREATED\n" << endi;
      }
      return;      
    }

    if (forcedcdFirst)
    {
      //  Open the DCD file
      iout << "OPENING FORCE DCD FILE\n" << endi;
	  
	#if OUTPUT_SINGLE_FILE
	  char *forcedcdFilename = simParams->forceDcdFilename;
	#else	  
	  char *forcedcdFilename = buildFileName(forcedcdType);
	#endif

      forcedcdFileID=open_dcd_write(forcedcdFilename);

      if (forcedcdFileID == DCD_FILEEXISTS)
      {
        char err_msg[257];
        sprintf(err_msg, "Force DCD file %s already exists!",forcedcdFilename);
        NAMD_err(err_msg);
      }
      else if (forcedcdFileID < 0)
      {
        char err_msg[257];
        sprintf(err_msg, "Couldn't open force DCD file %s",forcedcdFilename);
        NAMD_err(err_msg);
      }

	#if !OUTPUT_SINGLE_FILE
	  // Write out extra fields as MAGIC number etc.
	  int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	  float tmpFlt = OUTPUT_FILE_VERSION;
	  NAMD_write(forcedcdFileID, (char *) &tmpInt, sizeof(int32));
	  NAMD_write(forcedcdFileID, (char *) &tmpFlt, sizeof(float));
	  tmpInt = simParams->numoutputprocs;
	  NAMD_write(forcedcdFileID, (char *) &tmpInt, sizeof(int32));
	#endif

      int NSAVC, NFILE, NPRIV, NSTEP;
      NSAVC = simParams->forceDcdFrequency;
      NPRIV = timestep;
      NSTEP = NPRIV - NSAVC;
      NFILE = 0;

      //  Write out the header
      const int with_unitcell = 0;
      ret_code = write_dcdheader(forcedcdFileID, 
          forcedcdFilename,
          n, NFILE, NPRIV, NSAVC, NSTEP,
          simParams->dt/TIMEFACTOR, with_unitcell);


      if (ret_code<0)
      {
        NAMD_err("Writing of force DCD header failed!!");
      }

	#if !OUTPUT_SINGLE_FILE
	  //In this case, the file name is dynamically allocated
	  delete [] forcedcdFilename;
	#endif

      forcedcdFirst = FALSE;
    }

    //  Write out the values for this timestep
    iout << "WRITING FORCES TO DCD FILE AT STEP "
      << timestep << "\n" << endi;
    fflush(stdout);

	//In the case of writing to multiple files, only the header
	//of the dcd file needs to be updated. Note that the format of
	//the new dcd file has also changed! -Chao Mei 

#if OUTPUT_SINGLE_FILE
    //write X,Y,Z headers
    int totalAtoms = namdMyNode->molecule->numAtoms;
    write_dcdstep_par_XYZUnits(forcedcdFileID, totalAtoms);
#endif

    //update the header
    update_dcdstep_par_header(forcedcdFileID);    

}

void ParOutput::output_forcedcdfile_slave(int timestep, int fID, int tID, Vector *vecs){
    int ret_code;    //  Return code from DCD calls
    SimParameters *simParams = Node::Object()->simParameters;

    //  If this is the last time we will be writing coordinates,
    //  close the file before exiting
    if ( timestep == END_OF_RUN ) {
      if ( ! forcedcdFirst ) {        
        close_dcd_write(forcedcdFileID);
      }
#if OUTPUT_SINGLE_FILE
      delete [] forcedcdX;
      delete [] forcedcdY;
      delete [] forcedcdZ;
#endif
      return;
    }

    int parN = tID-fID+1;

    if (forcedcdFirst)
    {

	#if OUTPUT_SINGLE_FILE
	  char *forcedcdFilename = namdMyNode->simParams->forceDcdFilename;
	#else	  
	  char *forcedcdFilename = buildFileName(forcedcdType);
	#endif
      forcedcdFileID=open_dcd_write_par_slave(forcedcdFilename);
      if(forcedcdFileID < 0)
      {
        char err_msg[257];
        sprintf(err_msg, "Couldn't open force DCD file %s",forcedcdFilename);
        NAMD_err(err_msg);
      }
	#if OUTPUT_SINGLE_FILE
	  //If outputting to a single file, dcd files conforms to the old format
	  //as data are organized as 3 seperate arrays of X,Y,Z, while in the new
	  //format used in outputing multiple files, the data are organized as an
	  //array of XYZs.
      forcedcdX = new float[parN];
      forcedcdY = new float[parN];
      forcedcdZ = new float[parN];
	  //seek to beginning of X,Y,Z sections which means skipping header. 
	  //Cell data is not needed because this is force trajectory
	  int skipbytes = get_dcdheader_size();
	  seek_dcdfile(forcedcdFileID, skipbytes, SEEK_SET);
	#endif

	#if !OUTPUT_SINGLE_FILE
	  delete [] forcedcdFilename;
	  // Write out extra fields as MAGIC number etc.
	  int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	  float tmpFlt = OUTPUT_FILE_VERSION;
	  NAMD_write(forcedcdFileID, (char *) &tmpInt, sizeof(int32));
	  NAMD_write(forcedcdFileID, (char *) &tmpFlt, sizeof(float));	  
	  NAMD_write(forcedcdFileID, (char *) &outputID, sizeof(int));
	  NAMD_write(forcedcdFileID, (char *) &fID, sizeof(int));
	  NAMD_write(forcedcdFileID, (char *) &tID, sizeof(int));
	#endif
	  
      forcedcdFirst = FALSE;
    }

#if OUTPUT_SINGLE_FILE
    //The following seek will set the stream position to the
    //beginning of the place where a new timestep output should
    //be performed.
    CmiAssert(sizeof(off_t)==8);
    int totalAtoms = namdMyNode->molecule->numAtoms;

    for(int i=0; i<parN; i++){
        forcedcdX[i] = vecs[i].x;
        forcedcdY[i] = vecs[i].y;
        forcedcdZ[i] = vecs[i].z;
    }

    write_dcdstep_par_slave(forcedcdFileID, fID, tID, totalAtoms, forcedcdX, forcedcdY, forcedcdZ);
	//same with the slave output for coordiantes trajectory file
	//but cell data is not needed because this is force trajectory
	int atomsRemains = (totalAtoms-1)-(tID+1)+1;
	off_t offset = ((off_t)atomsRemains)*sizeof(float)+1*sizeof(int);
	seek_dcdfile(forcedcdFileID, offset, SEEK_CUR);
#else
	//write the timestep
	NAMD_write(forcedcdFileID, (char *)&timestep, sizeof(int));
	//write the values for this timestep
	NAMD_write(forcedcdFileID, (char *)vecs, sizeof(Vector)*parN);
#endif
}

void ParOutput::output_forces_master(int n){
#if OUTPUT_SINGLE_FILE
    char *output_name = new char[strlen(namdMyNode->simParams->outputFilename)+8];
    //  Build the output filename
    strcpy(output_name, namdMyNode->simParams->outputFilename);
    strcat(output_name, ".force");
#else	
	char *output_name = buildFileName(forceType);
#endif

    NAMD_backup_file(output_name);

    //Write the force to a binary file
    write_binary_file_master(output_name, n);
}

void ParOutput::output_forces_slave(int fID, int tID, Vector *vecs, int64 offset){
#if OUTPUT_SINGLE_FILE
    char *output_name = new char[strlen(namdMyNode->simParams->outputFilename)+8];
    //  Build the output filename
    strcpy(output_name, namdMyNode->simParams->outputFilename);
    strcat(output_name, ".force");
#else
	char *output_name = buildFileName(forceType);
	NAMD_backup_file(output_name);
#endif
    
    //Write the forces to a binary file
    write_binary_file_slave(output_name, fID, tID, vecs, offset);

	delete [] output_name;
}
//////End of Functions related with force output//////


//////Beginning of Functions related with coordinate output//////
void ParOutput::coordinateMaster(int timestep, int n, Lattice &lat){
    SimParameters *simParams = Node::Object()->simParameters;

    if ( timestep >= 0 ) {
      //  Output a DCD trajectory 
      if ( simParams->dcdFrequency &&
         ((timestep % simParams->dcdFrequency) == 0) )
      {        
        output_dcdfile_master(timestep, n, 
            simParams->dcdUnitCell ? &lat : NULL);
      }

      //  Output a restart file
      if ( simParams->restartFrequency &&
         ((timestep % simParams->restartFrequency) == 0) )
      {
        iout << "WRITING COORDINATES TO RESTART FILE AT STEP "
                  << timestep << "\n" << endi;
        output_restart_coordinates_master(timestep, n);
        iout << "FINISHED WRITING RESTART COORDINATES\n" <<endi;
        fflush(stdout);
      }

/*  Interactive MD is not supported in Parallel IO
      //  Interactive MD
      if ( simParams->IMDon &&
         ( ((timestep % simParams->IMDfreq) == 0) ||
           (timestep == simParams->firstTimestep) ) )
      {
        IMDOutput *imd = Node::Object()->imd;
        wrap_coor(fcoor,lattice,&fcoor_wrapped);
        if (imd != NULL) imd->gather_coordinates(timestep, n, fcoor);
      }
*/
    }

/*  EVAL_MEASURE of a timestep is not supported in Parallel IO 
    if (timestep == EVAL_MEASURE)
    {
  #ifdef NAMD_TCL
      wrap_coor(coor,lattice,&coor_wrapped);
      Node::Object()->getScript()->measure(coor);
  #endif
    }
*/
    //  Output final coordinates
    if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
    {
      int realstep = ( timestep == FILE_OUTPUT ?
          simParams->firstTimestep : simParams->N );
      iout << "WRITING COORDINATES TO OUTPUT FILE AT STEP "
                  << realstep << "\n" << endi;
      fflush(stdout);      
      output_final_coordinates_master(n);
    }

    //  Close trajectory files
    if (timestep == END_OF_RUN)
    {
      if (simParams->dcdFrequency) output_dcdfile_master(END_OF_RUN,0,NULL);
    }
}

void ParOutput::coordinateSlave(int timestep, int fID, int tID, Vector *vecs, FloatVector *fvecs){
    SimParameters *simParams = Node::Object()->simParameters;

    if ( timestep >= 0 ) {
      //  Output a DCD trajectory 
      if ( simParams->dcdFrequency &&
         ((timestep % simParams->dcdFrequency) == 0) )
      {        
        output_dcdfile_slave(timestep, fID, tID, fvecs);
      }

      //  Output a restart file
      if ( simParams->restartFrequency &&
         ((timestep % simParams->restartFrequency) == 0) )
      {
        int64 offset = sizeof(int)+sizeof(Vector)*((int64)fID);
        output_restart_coordinates_slave(timestep, fID, tID, vecs, offset);
      }

/*  Interactive MD is not supported in Parallel IO
      //  Interactive MD
      if ( simParams->IMDon &&
         ( ((timestep % simParams->IMDfreq) == 0) ||
           (timestep == simParams->firstTimestep) ) )
      {
        IMDOutput *imd = Node::Object()->imd;
        wrap_coor(fcoor,lattice,&fcoor_wrapped);
        if (imd != NULL) imd->gather_coordinates(timestep, n, fcoor);
      }
*/
    }

/*  EVAL_MEASURE of a timestep is not supported in Parallel IO 
    if (timestep == EVAL_MEASURE)
    {
  #ifdef NAMD_TCL
      wrap_coor(coor,lattice,&coor_wrapped);
      Node::Object()->getScript()->measure(coor);
  #endif
    }
*/
    //  Output final coordinates
    if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
    {
      int64 offset = sizeof(int)+sizeof(Vector)*((int64)fID);
      output_final_coordinates_slave(fID, tID, vecs, offset);
    }

    //  Close trajectory files
    if (timestep == END_OF_RUN)
    {
      if (simParams->dcdFrequency) output_dcdfile_slave(END_OF_RUN,0,0,NULL);
    }
}

void ParOutput::output_dcdfile_master(int timestep, int n, const Lattice *lattice){
    int ret_code;    //  Return code from DCD calls
    SimParameters *simParams = namdMyNode->simParams;

    //  If this is the last time we will be writing coordinates,
    //  close the file before exiting
    if ( timestep == END_OF_RUN ) {
      if ( ! dcdFirst ) {
        iout << "CLOSING COORDINATE DCD FILE\n" << endi;
        close_dcd_write(dcdFileID);
      } else {
        iout << "COORDINATE DCD FILE WAS NOT CREATED\n" << endi;
      }
      return;
    }

    if (dcdFirst)
    {
      //  Open the DCD file
      iout << "OPENING COORDINATE DCD FILE\n" << endi;

	#if OUTPUT_SINGLE_FILE
	  char *dcdFilename = simParams->dcdFilename;
	#else	  
	  char *dcdFilename = buildFileName(dcdType);
	#endif


      dcdFileID=open_dcd_write(dcdFilename);

      if (dcdFileID == DCD_FILEEXISTS)
      {
        char err_msg[257];
        sprintf(err_msg, "DCD file %s already exists!!",dcdFilename);
        NAMD_err(err_msg);
      }
      else if (dcdFileID < 0)
      {
        char err_msg[257];
        sprintf(err_msg, "Couldn't open DCD file %s",dcdFilename);
        NAMD_err(err_msg);
      }

	#if !OUTPUT_SINGLE_FILE
	  // Write out extra fields as MAGIC number etc.
	  int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	  float tmpFlt = OUTPUT_FILE_VERSION;
	  NAMD_write(dcdFileID, (char *) &tmpInt, sizeof(int32));
	  NAMD_write(dcdFileID, (char *) &tmpFlt, sizeof(float));
	  tmpInt = simParams->numoutputprocs;
	  NAMD_write(dcdFileID, (char *) &tmpInt, sizeof(int32));
	#endif

      int NSAVC, NFILE, NPRIV, NSTEP;
      NSAVC = simParams->dcdFrequency;
      NPRIV = timestep;
      NSTEP = NPRIV - NSAVC;
      NFILE = 0;

      //  Write out the header
      ret_code = write_dcdheader(dcdFileID, 
          dcdFilename,
          n, NFILE, NPRIV, NSAVC, NSTEP,
          simParams->dt/TIMEFACTOR, lattice != NULL);


      if (ret_code<0)
      {
        NAMD_err("Writing of DCD header failed!!");
      }

	  #if !OUTPUT_SINGLE_FILE
	  //dcdFilename needs to be freed as it is dynamically allocated
	  delete [] dcdFilename;
	  #endif

      dcdFirst = FALSE;
    }

    //  Write out the values for this timestep
    iout << "WRITING COORDINATES TO DCD FILE AT STEP "
      << timestep << "\n" << endi;
    fflush(stdout);

	//In the case of writing to multiple files, the header of the
	//dcd file needs to be updated. In addition, the lattice data
	//needs to be written if necessary. Note that the format of	
	//the new dcd file has also changed! -Chao Mei

    // Write out the Cell data
    if (lattice) {
      double unitcell[6];
      lattice_to_unitcell(lattice,unitcell);
      write_dcdstep_par_cell(dcdFileID, unitcell);
    }
            
#if OUTPUT_SINGLE_FILE
    //write X,Y,Z headers
    int totalAtoms = namdMyNode->molecule->numAtoms;
    write_dcdstep_par_XYZUnits(dcdFileID, totalAtoms);
#endif

    //update the header
    update_dcdstep_par_header(dcdFileID);
}
void ParOutput::output_dcdfile_slave(int timestep, int fID, int tID, FloatVector *fvecs){
    int ret_code;    //  Return code from DCD calls
    SimParameters *simParams = Node::Object()->simParameters;

    //  If this is the last time we will be writing coordinates,
    //  close the file before exiting
    if ( timestep == END_OF_RUN ) {
      if ( ! dcdFirst ) {        
        close_dcd_write(dcdFileID);
      }
#if OUTPUT_SINGLE_FILE
      delete [] dcdX;
      delete [] dcdY;
      delete [] dcdZ; 
#endif	       
      return;
    }

    int parN = tID-fID+1;

    if (dcdFirst)
    {

	#if OUTPUT_SINGLE_FILE
	  char *dcdFilename = simParams->dcdFilename;
	#else
	  char *dcdFilename = buildFileName(dcdType);
	#endif
      dcdFileID=open_dcd_write_par_slave(dcdFilename);
      if(dcdFileID < 0)
      {
        char err_msg[257];
        sprintf(err_msg, "Couldn't open DCD file %s", dcdFilename);
        NAMD_err(err_msg);
      }
	
	#if OUTPUT_SINGLE_FILE
      dcdX = new float[parN];
      dcdY = new float[parN];
      dcdZ = new float[parN];
	  //seek to beginning of X,Y,Z sections which means skipping header 
	  //skip the cell data if necessary
	  int skipbytes = get_dcdheader_size();
	  if(simParams->dcdUnitCell) {
		  skipbytes += sizeof(int)*2 + 6*sizeof(double);
	  }
	  seek_dcdfile(dcdFileID, skipbytes, SEEK_SET);
	#endif

	#if !OUTPUT_SINGLE_FILE
	  delete [] dcdFilename;

	  // Write out extra fields as MAGIC number etc.
	  int32 tmpInt = OUTPUT_MAGIC_NUMBER;
	  float tmpFlt = OUTPUT_FILE_VERSION;
	  NAMD_write(dcdFileID, (char *) &tmpInt, sizeof(int32));
	  NAMD_write(dcdFileID, (char *) &tmpFlt, sizeof(float));
	  NAMD_write(dcdFileID, (char *) &outputID, sizeof(int));
	  NAMD_write(dcdFileID, (char *) &fID, sizeof(int));
	  NAMD_write(dcdFileID, (char *) &tID, sizeof(int));
	#endif
      dcdFirst = FALSE;
    }

#if OUTPUT_SINGLE_FILE
    //The following seek will set the stream position to the
    //beginning of the place where a new timestep output should
    //be performed.
    CmiAssert(sizeof(off_t)==8);
    int totalAtoms = namdMyNode->molecule->numAtoms;

    for(int i=0; i<parN; i++){
        dcdX[i] = fvecs[i].x;
        dcdY[i] = fvecs[i].y;
        dcdZ[i] = fvecs[i].z;
    }

    write_dcdstep_par_slave(dcdFileID, fID, tID, totalAtoms, dcdX, dcdY, dcdZ);

	//1. Always foward to the beginning position of X,Y,Z sections in the
	//next timeframe.	
	//2. SHOULD AVOID USING SEEK_END: although the file size is updated on the
	//master proc, slave procs should rely on this update by seeking 
	//from SEEK_END because such info may not be updated in a timely manner
	//when this slave proc needs it.

	//We know the current position of file after slave writing is at the
	//end of Z output for atom "tID" (tID is the atom id indexed from 0)
	//(totalAtoms-1) is the last atom id, (tID+1) is the next atom id
	int atomsRemains = (totalAtoms-1)-(tID+1)+1;
	off_t offset = ((off_t)atomsRemains)*sizeof(float)+1*sizeof(int);
	//then skip the cell data if necessary
	if(simParams->dcdUnitCell) {
		offset += sizeof(int)*2 + 6*sizeof(double);		
	}
	seek_dcdfile(dcdFileID, offset, SEEK_CUR);
#else
	//write the timestep
	NAMD_write(dcdFileID, (char *)&timestep, sizeof(int));
	//write the values for this timestep
	NAMD_write(dcdFileID, (char *)fvecs, sizeof(FloatVector)*parN);
#endif
}

void ParOutput::output_restart_coordinates_master(int timestep, int n){
#if OUTPUT_SINGLE_FILE
	char timestepstr[20];

    int baselen = strlen(namdMyNode->simParams->restartFilename);
    char *restart_name = new char[baselen+26];

    strcpy(restart_name, namdMyNode->simParams->restartFilename);
    if ( namdMyNode->simParams->restartSave ) {
      sprintf(timestepstr,".%d",timestep);
      strcat(restart_name, timestepstr);
    }
    strcat(restart_name, ".coor");
#else
	char *restart_name = NULL;
	if ( namdMyNode->simParams->restartSave )
		restart_name = buildFileName(coorType);
	else
		restart_name = buildFileName(coorType,timestep);
#endif

    NAMD_backup_file(restart_name,".old");

    //  Generate a binary restart file
    write_binary_file_master(restart_name, n);

    delete [] restart_name;
}
void ParOutput::output_restart_coordinates_slave(int timestep, int fID, int tID, Vector *vecs, int64 offset){
#if OUTPUT_SINGLE_FILE
	char timestepstr[20];

    int baselen = strlen(namdMyNode->simParams->restartFilename);
    char *restart_name = new char[baselen+26];

    strcpy(restart_name, namdMyNode->simParams->restartFilename);
    if ( namdMyNode->simParams->restartSave ) {
      sprintf(timestepstr,".%d",timestep);
      strcat(restart_name, timestepstr);
    }
    strcat(restart_name, ".coor");
#else
	char *restart_name = NULL;
	if ( namdMyNode->simParams->restartSave )
		restart_name = buildFileName(coorType);
	else
		restart_name = buildFileName(coorType,timestep);

	NAMD_backup_file(restart_name,".old");
#endif

    //  Generate a binary restart file
    write_binary_file_slave(restart_name, fID, tID, vecs, offset);

    delete [] restart_name;
}

void ParOutput::output_final_coordinates_master(int n){
#if OUTPUT_SINGLE_FILE
	char *output_name = new char[strlen(namdMyNode->simParams->outputFilename)+8];

    //  Built the output filename
    strcpy(output_name, namdMyNode->simParams->outputFilename);
    strcat(output_name, ".coor");
#else	
	char *output_name = buildFileName(coorType);
#endif

    NAMD_backup_file(output_name);

    //  Write the coordinates to a binary file
    write_binary_file_master(output_name, n);

	delete [] output_name;
}
void ParOutput::output_final_coordinates_slave(int fID, int tID, Vector *vecs, int64 offset){
#if OUTPUT_SINGLE_FILE
    char *output_name = new char[strlen(namdMyNode->simParams->outputFilename)+8];

    //  Built the output filename
    strcpy(output_name, namdMyNode->simParams->outputFilename);
    strcat(output_name, ".coor");
#else	
	char *output_name = buildFileName(coorType);	
	NAMD_backup_file(output_name);
#endif

    //  Write the coordinates to a binary file
    write_binary_file_slave(output_name, fID, tID, vecs, offset);

	delete [] output_name;
}
//////End of Functions related with coordinate output//////

//////Beginning of Utility Functions for ParOutput//////
#if !OUTPUT_SINGLE_FILE
char *ParOutput::buildFileName(OUTPUTFILETYPE type, int timestep){
	char *filename = NULL;
	const char *typeName = NULL;
	switch(type) {
	case dcdType:
		typeName = "dcd";
		break;
	case forcedcdType:
		typeName = "forcedcd";
		break;
	case veldcdType:
		typeName = "veldcd";
		break;
	case coorType:
		typeName = "coor";
		break;
	case forceType:
		typeName = "force";
		break;
	case velType:
		typeName = "vel";
		break;
	default:
		typeName = "invalid";
		break;
	}
	int baselen = strlen(namdMyNode->simParams->outputFilename);	
	filename = new char[baselen+32];
	memset(filename, 0, baselen+32);

#if 0
	strcpy(filename, namdMyNode->simParams->outputFilename);

	//check if the directory exists or not
	if(access(filename, F_OK)!=0) {
		int ret = MKDIR(filename);
		if(ret!=0) {
			char errmsg[512];
			sprintf(errmsg, "Error in creating top-level directory %s!", filename);
			NAMD_err(errmsg);
		}
	}

	strcat(filename, PATHSEPSTR);
	strcat(filename, typeName);

	//check if the directory exists or not	
	if(access(filename, F_OK)!=0) {
		int ret = MKDIR(filename);
		if(ret!=0) {
			char errmsg[512];
			sprintf(errmsg, "Error in creating middle-level directory %s!", filename);
			NAMD_err(errmsg);
		}
	}
#else
	sprintf(filename, "%s%s%s", namdMyNode->simParams->outputFilename, PATHSEPSTR, typeName);
#endif

	char tmpstr[20];
	if(outputID == -1) {
		//indicating the output from master			
		if(timestep!=-9999) {
			//not the default value
			sprintf(tmpstr, "%smeta.%d", PATHSEPSTR, timestep);
		}else{
			sprintf(tmpstr, "%smeta", PATHSEPSTR);
		}
	}else{
		//indicating the output from slave		
		sprintf(tmpstr, "%s%d", PATHSEPSTR, outputID);
		strcat(filename, tmpstr);
		#if 0
		if(access(filename, F_OK)!=0) {
			int ret = MKDIR(filename);
			if(ret!=0) {
				char errmsg[512];
				sprintf(errmsg, "Error in creating last-level directory %s!", filename);
				NAMD_err(errmsg);
			}
		}
		#endif
		if(timestep!=-9999) {
			//not the default value
			sprintf(tmpstr, "%s%s_%d.%d", PATHSEPSTR,typeName,outputID,timestep);
		}else{
			sprintf(tmpstr,"%s%s_%d", PATHSEPSTR,typeName,outputID);
		}	
	}

	strcat(filename, tmpstr);
	return filename;
}
#endif
//////End of Utility Functions for ParOutput//////
#endif
