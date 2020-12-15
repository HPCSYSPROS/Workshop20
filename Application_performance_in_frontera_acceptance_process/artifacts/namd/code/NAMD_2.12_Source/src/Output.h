/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   This object outputs the data collected on the master  node
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include "common.h"
#include <string>
#include <map>

class Vector;
class FloatVector;
class Lattice;
class ReplicaDcdInitMsg;
class ReplicaDcdDataMsg;

// semaphore "steps", must be negative
#define FILE_OUTPUT -1
#define END_OF_RUN -2
#define EVAL_MEASURE -3
#define FORCE_OUTPUT -4

#define OUTPUT_SINGLE_FILE 1
#define OUTPUT_MAGIC_NUMBER 123456
#define OUTPUT_FILE_VERSION 1.00

enum OUTPUTFILETYPE {
	dcdType,
	forcedcdType,
	veldcdType,
	coorType,
	forceType,
	velType
};

class Output 
{

private:

friend class SimParameters;

   //  output coords to dcd file
   //  Pass non-NULL Lattice to include unit cell in the timesteps.
   int output_dcdfile(int, int, FloatVector *, const Lattice *); 
   void output_veldcdfile(int, int, Vector *); 	//  output velocities to
						//  dcd file
   void output_forcedcdfile(int, int, Vector *); //  output forces to
						//  dcd file

   void output_restart_coordinates(Vector *, int, int);
						//  output coords to 
						//  restart file
   void output_restart_velocities(int, int, Vector *);
						//  output velocities to 
						//  restart file
   void output_final_coordinates(Vector *, int, int);//  output final coordinates
   void output_final_velocities(int, int, Vector *);	//  output final coordinates
   void output_forces(int, int, Vector *);	//  output forces

   void scale_vels(Vector *, int, Real);	//  scale velocity vectors before output
   void write_binary_file(char *, int, Vector *); // Write a binary restart file with
						//  coordinates or velocities

   struct replicaDcdFile {
     std::string filename;
     int fileid;
     replicaDcdFile() : fileid(0) { ; }
   };
   std::map<int,replicaDcdFile> replicaDcdFiles;
   int replicaDcdActive;
   int replicaDcdIndex;

public :
   Output();					//  Constructor
   ~Output();					//  Destructor
   void energy(int, BigReal *);			//  Output energies

   static int coordinateNeeded(int);
   void coordinate(int, int, Vector *, FloatVector *, Lattice &);
						//  Produce appropriate 
						//  coordinate output for 
						//  the current timestep
   static int velocityNeeded(int);
   void velocity(int, int, Vector *);		//  Produce appropriate velocity
						//  output for the current 
						//  timestep
   static int forceNeeded(int);
   void force(int, int, Vector *);		//  Produce appropriate force
						//  output for the current 
						//  timestep

  void replicaDcdOff() { replicaDcdActive = 0; }
  void setReplicaDcdIndex(int index);
  void replicaDcdInit(int index, const char *filename);
  void recvReplicaDcdInit(ReplicaDcdInitMsg *msg);
  void recvReplicaDcdData(ReplicaDcdDataMsg *msg);
};

#ifdef MEM_OPT_VERSION
class ParOutput{
private:
    void output_veldcdfile_master(int timestep, int n);
    void output_veldcdfile_slave(int timestep, int fID, int tID, Vector *vecs);
    void output_restart_velocities_master(int timestep, int n);
    void output_restart_velocities_slave(int timestep, int fID, int tID, Vector *vecs, int64 offset);
    void output_final_velocities_master(int n);
    void output_final_velocities_slave(int fID, int tID, Vector *vecs, int64 offset);

    void output_forcedcdfile_master(int timestep, int n);
    void output_forcedcdfile_slave(int timestep, int fID, int tID, Vector *vecs);
    void output_forces_master(int n);
    void output_forces_slave(int fID, int tID, Vector *vecs, int64 offset);

    void output_dcdfile_master(int timestep, int n, const Lattice *lat);
    void output_dcdfile_slave(int timestep, int fID, int tID, FloatVector *fvecs);
    void output_restart_coordinates_master(int timestep, int n);
    void output_restart_coordinates_slave(int timestep, int fID, int tID, Vector *vecs, int64 offset);
    void output_final_coordinates_master(int n);
    void output_final_coordinates_slave(int fID, int tID, Vector *vecs, int64 offset);

    void write_binary_file_master(char *fname, int n);
    void write_binary_file_slave(char *fname, int fID, int tID, Vector *vecs, int64 offset);

	#if !OUTPUT_SINGLE_FILE
	char *buildFileName(OUTPUTFILETYPE type, int timestep=-9999);
	#endif

    int dcdFileID;
    Bool dcdFirst;
    float *dcdX, *dcdY, *dcdZ;

    int veldcdFileID;
    Bool veldcdFirst;    
    float *veldcdX, *veldcdY, *veldcdZ;

    int forcedcdFileID;
    Bool forcedcdFirst;    
    float *forcedcdX, *forcedcdY, *forcedcdZ;

	int outputID; //the sequence of this output

public:
    ParOutput(int oid=-1){
        dcdFileID=veldcdFileID=-99999;
        forcedcdFileID=-99999;
        dcdFirst=veldcdFirst=TRUE;
        forcedcdFirst=TRUE;
        dcdX=dcdY=dcdZ=veldcdX=veldcdY=veldcdZ=NULL;
        forcedcdX=forcedcdY=forcedcdZ=NULL;
		outputID=oid;
    }
    ~ParOutput() {}

    void velocityMaster(int timestep, int n);
    void velocitySlave(int timestep, int fID, int tID, Vector *vecs);

    void forceMaster(int timestep, int n);
    void forceSlave(int timestep, int fID, int tID, Vector *vecs);

    void coordinateMaster(int timestep, int n, Lattice &lat);
    void coordinateSlave(int timestep, int fID, int tID, Vector *vecs, FloatVector *fvecs);
};
#endif

#endif

