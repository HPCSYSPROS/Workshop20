/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "largefiles.h"  // must be first!

#include "InfoStream.h"
#include "CollectionMaster.h"
#include "ProcessorPrivate.h"
#include "SimParameters.h"
#include "packmsg.h"
#include "CollectionMaster.decl.h"
#include "Molecule.h"

#include "memusage.h"

// #define DEBUGM
#include "Debug.h"

CollectionMaster::CollectionMaster()
{
  if (CkpvAccess(CollectionMaster_instance) == 0) {
    CkpvAccess(CollectionMaster_instance) = this;
  } else {
    DebugM(1, "CollectionMaster::CollectionMaster() - another instance of CollectionMaster exists!\n");
  }
  dataStreamFile = 0;

#ifdef MEM_OPT_VERSION
  wrapCoorDoneCnt = 0;
  posDoneCnt = 0;
  velDoneCnt = 0;
  parOut = new ParOutput();
#endif

  posTimings = 10;  velTimings = forceTimings = 5;
}


CollectionMaster::~CollectionMaster(void)
{
}

void CollectionMaster::receivePositions(CollectVectorMsg *msg)
{
#ifndef MEM_OPT_VERSION
  positions.submitData(msg,Node::Object()->molecule->numAtoms);
  delete msg;
  
  CollectVectorInstance *c;
  while ( ( c = positions.removeReady() ) ) { disposePositions(c); }
#endif
}

void CollectionMaster::enqueuePositions(int seq, Lattice &lattice)
{
  positions.enqueue(seq,lattice);

#ifndef MEM_OPT_VERSION
  CollectVectorInstance *c;
  while ( ( c = positions.removeReady() ) ) { disposePositions(c); }
#else
  checkPosReady();
#endif
}

void CollectionMaster::disposePositions(CollectVectorInstance *c)
{
#ifndef MEM_OPT_VERSION
    DebugM(3,"Collected positions at " << c->seq << std::endl);
    int seq = c->seq;
    int size = c->data.size();
    if ( ! size ) size = c->fdata.size();
    Vector *data = c->data.begin();
    FloatVector *fdata = c->fdata.begin();
    double exectime = CmiWallTimer();
    double mem = memusage_MB();
    Node::Object()->output->coordinate(seq,size,data,fdata,c->lattice);
    c->free();
    exectime = CmiWallTimer()-exectime;
    if ( posTimings ) {
      CkPrintf("The last position output (seq=%d) takes %.3f seconds, %.3f MB of memory in use\n", seq, exectime, mem);
      --posTimings;
    }
#endif
}


void CollectionMaster::receiveVelocities(CollectVectorMsg *msg)
{
#ifndef MEM_OPT_VERSION
  velocities.submitData(msg,Node::Object()->molecule->numAtoms);
  delete msg;

  CollectVectorInstance *c;
  while ( ( c = velocities.removeReady() ) ) { disposeVelocities(c); }
#endif
}

void CollectionMaster::enqueueVelocities(int seq)
{
  Lattice dummy;
  velocities.enqueue(seq,dummy);
#ifndef MEM_OPT_VERSION
  CollectVectorInstance *c;
  while ( ( c = velocities.removeReady() ) ) { disposeVelocities(c); }
#else
  checkVelReady();
#endif
}

void CollectionMaster::disposeVelocities(CollectVectorInstance *c)
{
#ifndef MEM_OPT_VERSION
    DebugM(3,"Collected velocities at " << c->seq << std::endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = c->data.begin();
    double exectime = CmiWallTimer();
    double mem = memusage_MB();
    Node::Object()->output->velocity(seq,size,data);
    c->free();
    exectime = CmiWallTimer()-exectime;
    if ( velTimings ) {
      CkPrintf("The last velocity output (seq=%d) takes %.3f seconds, %.3f MB of memory in use\n", seq, exectime, mem);
      --velTimings;
    }
#endif
}


void CollectionMaster::receiveForces(CollectVectorMsg *msg)
{
#ifndef MEM_OPT_VERSION
  forces.submitData(msg,Node::Object()->molecule->numAtoms);
  delete msg;

  CollectVectorInstance *c;
  while ( ( c = forces.removeReady() ) ) { disposeForces(c); }
#endif
}

void CollectionMaster::enqueueForces(int seq)
{
  Lattice dummy;
  forces.enqueue(seq,dummy);
#ifndef MEM_OPT_VERSION
  CollectVectorInstance *c;
  while ( ( c = forces.removeReady() ) ) { disposeForces(c); }
#else
  checkForceReady();
#endif
}

void CollectionMaster::disposeForces(CollectVectorInstance *c)
{
#ifndef MEM_OPT_VERSION
    DebugM(3,"Collected forces at " << c->seq << std::endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = c->data.begin();
    double exectime = CmiWallTimer();
    double mem = memusage_MB();
    Node::Object()->output->force(seq,size,data);
    c->free();
    exectime = CmiWallTimer()-exectime;
    if ( forceTimings ) {
      CkPrintf("The last force output (seq=%d) takes %.3f seconds, %.3f MB of memory in use\n", seq, exectime, mem);
      --forceTimings;
    }
#endif
}


void CollectionMaster::receiveDataStream(DataStreamMsg *msg) {
    if ( ! dataStreamFile ) {
      char *fname = Node::Object()->simParameters->auxFilename;
      // iout has large file linking issues on AIX
      // iout << iINFO << "OPENING AUXILIARY DATA STREAM FILE "
      // 				<< fname << "\n" << endi;
      CkPrintf("Info: OPENING AUXILIARY DATA STREAM FILE %s\n", fname);
      NAMD_backup_file(fname);
      dataStreamFile = fopen(fname,"w");
      if ( ! dataStreamFile )
		NAMD_die("Can't open auxiliary data stream file!");
    }
    fprintf(dataStreamFile,"%s",msg->data.begin());
    fflush(dataStreamFile);
    delete msg;
}

PACK_MSG(DataStreamMsg,
  PACK_RESIZE(data);
)

//The timesteps on CollectionMaster will be always in increasing order
//because they are enqueued on PE0 by controller in order. -Chao Mei
//
//The computation of wrap_coor is also serialized for the sake of easy
//implementation (easier management of output)
void CollectionMaster::receiveOutputPosReady(int seq){
#ifdef MEM_OPT_VERSION
    positions.submitData(seq);
    checkPosReady();
#endif
}

void CollectionMaster::receiveOutputVelReady(int seq){
#ifdef MEM_OPT_VERSION
    velocities.submitData(seq);
    checkVelReady();
#endif
}

void CollectionMaster::receiveOutputForceReady(int seq){
#ifdef MEM_OPT_VERSION
    forces.submitData(seq);
    checkForceReady();
#endif
}


void CollectionMaster::startNextRoundOutputPos(double totalT){
#ifdef MEM_OPT_VERSION

	if(totalT > posIOTime) posIOTime = totalT;

#ifndef OUTPUT_SINGLE_FILE
#error OUTPUT_SINGLE_FILE not defined!
#endif

#if OUTPUT_SINGLE_FILE
    if(++posDoneCnt < Node::Object()->simParameters->numoutputwrts)  return;
#else
	if(++posDoneCnt < Node::Object()->simParameters->numoutputprocs)  return;
#endif

	posDoneCnt = 0;

    //retrieve the last ready instance
    CollectVectorInstance *c = positions.getReady();
    int seq = c->seq;
    CmiAssert(c->status == IN_PROCESS);
    double mem = memusage_MB();
    positions.removeFirstReady();
    c->free();
    posOutTime = CmiWallTimer()-posOutTime;
    if ( posTimings ) {
      CkPrintf("The last position output (seq=%d) takes %.3f seconds(file I/O: %.3f secs), %.3f MB of memory in use\n", seq, posOutTime, posIOTime, mem);
      --posTimings;
    }

    //Actually the c->status doesn't need to be checked because it is
    //certain that the new ready one will not be in  IN_PROCESS status 
    checkPosReady();
#endif
}

void CollectionMaster::startNextRoundOutputVel(double totalT){
#ifdef MEM_OPT_VERSION
	
	if(totalT > velIOTime) velIOTime = totalT;

#if OUTPUT_SINGLE_FILE
    if(++velDoneCnt < Node::Object()->simParameters->numoutputwrts)  return;
#else
	if(++velDoneCnt < Node::Object()->simParameters->numoutputprocs)  return;
#endif

    velDoneCnt = 0;

    //retrieve the last ready instance
    CollectVectorInstance *c = velocities.getReady();
    int seq = c->seq;
    CmiAssert(c->status == IN_PROCESS);
    double mem = memusage_MB();
    velocities.removeFirstReady();
    c->free();
    velOutTime = CmiWallTimer()-velOutTime;
    if ( velTimings ) {
      CkPrintf("The last velocity output (seq=%d) takes %.3f seconds(file I/O: %.3f secs), %.3f MB of memory in use\n", seq, velOutTime, velIOTime, mem);
      --velTimings;
    }

    //Actually the c->status doesn't need to be checked because it is
    //certain that the new ready one will not be in  IN_PROCESS status 
    checkVelReady();
#endif
}

void CollectionMaster::startNextRoundOutputForce(double totalT){
#ifdef MEM_OPT_VERSION
	
	if(totalT > forceIOTime) forceIOTime = totalT;

#if OUTPUT_SINGLE_FILE
    if(++forceDoneCnt < Node::Object()->simParameters->numoutputwrts)  return;
#else
	if(++forceDoneCnt < Node::Object()->simParameters->numoutputprocs)  return;
#endif

    forceDoneCnt = 0;

    //retrieve the last ready instance
    CollectVectorInstance *c = forces.getReady();
    int seq = c->seq;
    CmiAssert(c->status == IN_PROCESS);
    double mem = memusage_MB();
    forces.removeFirstReady();
    c->free();
    forceOutTime = CmiWallTimer()-forceOutTime;
    if ( forceTimings ) {
      CkPrintf("The last force output (seq=%d) takes %.3f seconds(file I/O: %.3f secs), %.3f MB of memory in use\n", seq, forceOutTime, forceIOTime, mem);
      --forceTimings;
    }

    //Actually the c->status doesn't need to be checked because it is
    //certain that the new ready one will not be in  IN_PROCESS status 
    checkForceReady();
#endif
}


void CollectionMaster::wrapCoorFinished(){
#ifdef MEM_OPT_VERSION
    if(++wrapCoorDoneCnt == Node::Object()->simParameters->numoutputprocs){
        wrapCoorDoneCnt = 0;

		//count the wrapping-coor time into master writing time
		posIOTime = CmiWallTimer()-posOutTime; 

        //it's ready to output positions
        CollectVectorInstance *c = positions.getReady();

		CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
		ParallelIOMgr *ioMgr = io.ckLocalBranch();

#if OUTPUT_SINGLE_FILE
        //notify output procs to do Token based output
        int grpsize = ioMgr->numOutputProcs / ioMgr->numOutputWrts;
        int remains = ioMgr->numOutputProcs % ioMgr->numOutputWrts;
        int outrank = 0;
        int i;
        for(i=0; i<remains; i++){
            io[ioMgr->outputProcArray[outrank]].disposePositions(c->seq, posIOTime);
            outrank += (grpsize+1);
        }
        for(; i<ioMgr->numOutputWrts; i++){
            io[ioMgr->outputProcArray[outrank]].disposePositions(c->seq, posIOTime);
            outrank += grpsize;
        }
#else
		//output multiple files
		for(int i=0; i<ioMgr->numOutputProcs; i++) {
			io[ioMgr->outputProcArray[i]].disposePositions(c->seq, posIOTime);
		}
#endif

    }
#endif
}

#ifdef MEM_OPT_VERSION
void CollectionMaster::checkPosReady(){
    CollectVectorInstance *c;
    if((c = positions.getReady())){
        if(c->status == IN_PROCESS){
            //indicating in the process of outputing coordinates
            return;
        }        
        c->status = IN_PROCESS;

        posOutTime = CmiWallTimer();
        SimParameters *simParam = Node::Object()->simParameters;
        CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
        ParallelIOMgr *ioMgr = io.ckLocalBranch();
        if(simParam->wrapAll || simParam->wrapWater){
            for(int i=0; i<ioMgr->numOutputProcs; i++){
                io[ioMgr->outputProcArray[i]].wrapCoor(c->seq, c->lattice);
            }
            //write the header to overlap with the computation of 
            //wrapping coordinates
            parOut->coordinateMaster(c->seq,Node::Object()->molecule->numAtoms,c->lattice);
        }else{
            //write the header 
            parOut->coordinateMaster(c->seq,Node::Object()->molecule->numAtoms,c->lattice);
			posIOTime = CmiWallTimer() - posOutTime;            

		#if OUTPUT_SINGLE_FILE
            int grpsize = ioMgr->numOutputProcs / ioMgr->numOutputWrts;
            int remains = ioMgr->numOutputProcs % ioMgr->numOutputWrts;
            int outrank = 0;
            int i;
            for(i=0; i<remains; i++){
                io[ioMgr->outputProcArray[outrank]].disposePositions(c->seq, posIOTime);
                outrank += (grpsize+1);
            }
            for(; i<ioMgr->numOutputWrts; i++){
                io[ioMgr->outputProcArray[outrank]].disposePositions(c->seq, posIOTime);
                outrank += grpsize;
            }
		#else
			//output multiple files
			for(int i=0; i<ioMgr->numOutputProcs; i++) {
				io[ioMgr->outputProcArray[i]].disposePositions(c->seq, posIOTime);
			}
		#endif
        }
        //this instance c is freed in the next round of output invocation.
    }
}

void CollectionMaster::checkVelReady(){
    CollectVectorInstance *c;
    if((c = velocities.getReady())){
        if(c->status == IN_PROCESS){
            //indicating in the process of outputing velocities
            return;
        }

        c->status = IN_PROCESS;

        velOutTime = CmiWallTimer();
        //write the header
        parOut->velocityMaster(c->seq, Node::Object()->molecule->numAtoms);
		velIOTime = CmiWallTimer() - velOutTime;

        //notify output procs to do Token based output
        CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
        ParallelIOMgr *ioMgr = io.ckLocalBranch();

	#if OUTPUT_SINGLE_FILE
        int grpsize = ioMgr->numOutputProcs / ioMgr->numOutputWrts;
        int remains = ioMgr->numOutputProcs % ioMgr->numOutputWrts;
        int outrank = 0;
        int i;
        for(i=0; i<remains; i++){
            io[ioMgr->outputProcArray[outrank]].disposeVelocities(c->seq, velIOTime);
            outrank += (grpsize+1);
        }
        for(; i<ioMgr->numOutputWrts; i++){
            io[ioMgr->outputProcArray[outrank]].disposeVelocities(c->seq, velIOTime);
            outrank += grpsize;
        }
	#else
		//output multiple files
		for(int i=0; i<ioMgr->numOutputProcs; i++) {
			io[ioMgr->outputProcArray[i]].disposeVelocities(c->seq, velIOTime);
		}
	#endif
        //this instance c is freed in the next round of output invocation.        
    }
}

void CollectionMaster::checkForceReady(){
    CollectVectorInstance *c;
    if((c = forces.getReady())){
        if(c->status == IN_PROCESS){
            //indicating in the process of outputing forces
            return;
        }

        c->status = IN_PROCESS;

        forceOutTime = CmiWallTimer();
        //write the header
        parOut->forceMaster(c->seq, Node::Object()->molecule->numAtoms);
		forceIOTime = CmiWallTimer() - forceOutTime;

        //notify output procs to do Token based output
        CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
        ParallelIOMgr *ioMgr = io.ckLocalBranch();

	#if OUTPUT_SINGLE_FILE
        int grpsize = ioMgr->numOutputProcs / ioMgr->numOutputWrts;
        int remains = ioMgr->numOutputProcs % ioMgr->numOutputWrts;
        int outrank = 0;
        int i;
        for(i=0; i<remains; i++){
            io[ioMgr->outputProcArray[outrank]].disposeForces(c->seq, forceIOTime);
            outrank += (grpsize+1);
        }
        for(; i<ioMgr->numOutputWrts; i++){
            io[ioMgr->outputProcArray[outrank]].disposeForces(c->seq, forceIOTime);
            outrank += grpsize;
        }
	#else
		//output multiple files
		for(int i=0; i<ioMgr->numOutputProcs; i++) {
			io[ioMgr->outputProcArray[i]].disposeForces(c->seq, forceIOTime);
		}
	#endif
        //this instance c is freed in the next round of output invocation.        
    }
}


void CollectionMidMaster::disposePositions(int seq)
{
    CollectMidVectorInstance *c = positions.getReady(seq);
    CmiAssert(c!=NULL);
    parOut->coordinateSlave(seq,c->fromAtomID,c->toAtomID,
                               c->data.begin(),c->fdata.begin());
    c->free();
}

void CollectionMidMaster::disposeVelocities(int seq)
{
    CollectMidVectorInstance *c = velocities.getReady(seq);
    CmiAssert(c!=NULL);    
    parOut->velocitySlave(seq,c->fromAtomID,c->toAtomID,c->data.begin()); 
    c->free();
}

void CollectionMidMaster::disposeForces(int seq)
{
    CollectMidVectorInstance *c = forces.getReady(seq);
    CmiAssert(c!=NULL);    
    parOut->forceSlave(seq,c->fromAtomID,c->toAtomID,c->data.begin()); 
    c->free();
}
#endif

#include "CollectionMaster.def.h"

