/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COLLECTIONMASTER_H
#define COLLECTIONMASTER_H

#include "charm++.h"
#include "main.h"
#include "NamdTypes.h"
#include "Lattice.h"
#include "ProcessorPrivate.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "CollectionMaster.decl.h"
#include <stdio.h>

#include "Node.h"
#include "ParallelIOMgr.decl.h"
#include "ParallelIOMgr.h"
#include "Output.h"

class CollectVectorMsg : public CMessage_CollectVectorMsg
{
public:

  int seq;
  int aid_size;
  int data_size;
  int fdata_size;
  AtomID *aid;
  Vector *data;
  FloatVector *fdata;
};

class DataStreamMsg;

class CollectionMaster : public Chare
{
public:
  static CollectionMaster *Object() { 
    return CkpvAccess(CollectionMaster_instance); 
  }
  CollectionMaster();
  ~CollectionMaster(void);

  void receivePositions(CollectVectorMsg *msg);
  void receiveVelocities(CollectVectorMsg *msg);
  void receiveForces(CollectVectorMsg *msg);

  void receiveDataStream(DataStreamMsg *msg);

  void enqueuePositions(int seq, Lattice &lattice);
  void enqueueVelocities(int seq);
  void enqueueForces(int seq);

  class CollectVectorInstance;
  void disposePositions(CollectVectorInstance *c);
  void disposeVelocities(CollectVectorInstance *c);
  void disposeForces(CollectVectorInstance *c);

  void blockPositions() { positions.block(); }
  void unblockPositions() { positions.unblock(); }

  /////Beginning of declarations for comm with CollectionMidMaster/////
  void receiveOutputPosReady(int seq);
  void receiveOutputVelReady(int seq);
  void receiveOutputForceReady(int seq);
  //totalT is the time taken to do file I/O for each output workflow -Chao Mei  
  void startNextRoundOutputPos(double totalT);
  void startNextRoundOutputVel(double totalT);
  void startNextRoundOutputForce(double totalT);

  void wrapCoorFinished();

  enum OperationStatus {NOT_PROCESSED, IN_PROCESS, HAS_PROCESSED};
  /////End of declarations for comm with CollectionMidMaster/////
#ifdef MEM_OPT_VERSION
  class CollectVectorInstance
  {
  public:
    CollectVectorInstance(void) : seq(-10) { ; }

    CollectVectorInstance(int s) { reset(s); }

    void free() { seq = -10; status = HAS_PROCESSED; }
    int notfree() { return ( seq != -10 ); }

    void reset(int s) {
      if ( s == -10 ) NAMD_bug("seq == free in CollectionMaster");
      seq = s;
      CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);      
      ParallelIOMgr *ioMgr = io.ckLocalBranch();
      remaining = ioMgr->getNumOutputProcs();
      status = NOT_PROCESSED;          
    }
    
    void append(){ --remaining; }
    
    int ready(void) { return ( ! remaining ); }

    int seq;
    Lattice lattice;
    //mainly used for tracking the progress of wrap_coor operation
    //the write to files will not happen until the wrap_coor is finished,
    //and the CollectVectorInstance is freed. -Chao Mei
    OperationStatus status; 
  private:
    int remaining;

  }; //end of declaration for CollectionMaster::CollectVectorInstance

  class CollectVectorSequence
  {
  public:

    void submitData(int seq){
      CollectVectorInstance **c = data.begin();
      CollectVectorInstance **c_e = data.end();
      for( ; c != c_e && (*c)->seq != seq; ++c );
      if ( c == c_e )
      {
        c = data.begin();
        for( ; c != c_e && (*c)->notfree(); ++c );
        if ( c == c_e ) {
          data.add(new CollectVectorInstance(seq));
          c = data.end() - 1;
        }
        (*c)->reset(seq);
      }
      (*c)->append();
    }

    void enqueue(int seq, Lattice &lattice) {
      queue.add(seq);
      latqueue.add(lattice);
    }

    CollectVectorInstance* removeReady(void)
    {
      //it is equal to
      //if(getReady()) removeFirstReady();
      //But it should be not used in memory
      //optimized version as the ready instance
      //is delayed to be freed at the start of
      //next round of output -Chao Mei
      return NULL;
    }
    void block() { ; }  // unimplemented
    void unblock() { ; }  // unimplemented

    //only get the ready instance, not remove their info
    //from timestep queue and lattice queue
    CollectVectorInstance* getReady(void)
    {
      CollectVectorInstance *o = 0;
      if ( queue.size() )
      {
        int seq = queue[0];
        CollectVectorInstance **c = data.begin();
        CollectVectorInstance **c_e = data.end();
        for( ; c != c_e && (*c)->seq != seq; ++c );
        if ( c != c_e && (*c)->ready() )
        {
          o = *c;
          o->lattice = latqueue[0];
        }
      }
      return o;
    }

    //the function is intended to be used after "getReady"
    //to remove the info regarding the timestep and lattice.
    //So, it removes the front ready one. -Chao Mei
    int removeFirstReady(){
      int seq = queue[0];
      queue.del(0,1);
      latqueue.del(0,1);
      return seq;
    }

    ResizeArray<CollectVectorInstance*> data;
    ResizeArray<int> queue;
    ResizeArray<Lattice> latqueue;
  }; //end of declaration for CollectionMaster::CollectVectorSequence
#else
  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-10) { ; }

    CollectVectorInstance(int s) { reset(s); }

    void free() { seq = -10; }
    int notfree() { return ( seq != -10 ); }

    void reset(int s) {
	if ( s == -10 ) NAMD_bug("seq == free in CollectionMaster");
        seq = s;
        remaining = (PatchMap::Object())->numNodesWithPatches(); 
	data.resize(0);
	fdata.resize(0);
    }

    // true -> send it and delete it!
    void append(CollectVectorMsg *msg, int max_index)
    {
      AtomID *a = msg->aid;
      Vector *d = msg->data;
      FloatVector *fd = msg->fdata;
      int size = msg->aid_size;
      if ( msg->data_size ) {
        data.resize(max_index);
        Vector *ptr = data.begin();
	for( int i = 0; i < size; ++i ) { ptr[a[i]] = d[i]; }
      }
      if ( msg->fdata_size ) {
        fdata.resize(max_index);
        FloatVector *ptr = fdata.begin();
	for( int i = 0; i < size; ++i ) { ptr[a[i]] = fd[i]; }
      }
      --remaining;
    }

    int ready(void) { return ( ! remaining ); }

    int seq;
    Lattice lattice;

    ResizeArray<Vector> data;
    ResizeArray<FloatVector> fdata;

  private:
    int remaining;

  };

  class CollectVectorSequence
  {
  public:

    void submitData(CollectVectorMsg *msg, int max_index)
    {
      int seq = msg->seq;
      CollectVectorInstance **c = data.begin();
      CollectVectorInstance **c_e = data.end();
      for( ; c != c_e && (*c)->seq != seq; ++c );
      if ( c == c_e )
      {
        c = data.begin();
        for( ; c != c_e && (*c)->notfree(); ++c );
        if ( c == c_e ) {
          data.add(new CollectVectorInstance(seq));
          c = data.end() - 1;
        }
        (*c)->reset(seq);
      }
      (*c)->append(msg, max_index);
    }

    void enqueue(int seq, Lattice &lattice) {
      queue.add(seq);
      latqueue.add(lattice);
    }

    CollectVectorInstance* removeReady(void)
    {
      CollectVectorInstance *o = 0;
      if ( queue.size() && ! blocked )
      {
        int seq = queue[0];
        CollectVectorInstance **c = data.begin();
        CollectVectorInstance **c_e = data.end();
        for( ; c != c_e && (*c)->seq != seq; ++c );
        if ( c != c_e && (*c)->ready() )
        {
	  o = *c;
	  o->lattice = latqueue[0];
	  queue.del(0,1);
	  latqueue.del(0,1);
        }
      }
      return o;
    }

    void block() { blocked = 1; }
    void unblock() { blocked = 0; }
    CollectVectorSequence() : blocked(0) { ; }
    
    ResizeArray<CollectVectorInstance*> data;
    ResizeArray<int> queue;
    ResizeArray<Lattice> latqueue;
    int blocked;

  };
#endif
private:

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;
  int posTimings, velTimings, forceTimings;
  FILE *dataStreamFile;

#ifdef MEM_OPT_VERSION
  int wrapCoorDoneCnt;
  ParOutput *parOut;
  double posOutTime; //record the output time
  double velOutTime; //record the output time
  double forceOutTime; //record the output time
  double posIOTime; //record the max time spent on real file IO for one output
  double velIOTime; //record the max time spent on real file IO for one output
  double forceIOTime; //record the max time spent on real file IO for one output
  
  //for the sake of simultaneous writing to the same file
  int posDoneCnt;
  int velDoneCnt;
  int forceDoneCnt;

  void checkPosReady();
  void checkVelReady();
  void checkForceReady();
#endif
};

class DataStreamMsg : public CMessage_DataStreamMsg {
public:

  ResizeArray<char> data;

  static void* pack(DataStreamMsg* msg);
  static DataStreamMsg* unpack(void *ptr);

};

//Use varsize message to be more SMP safe 
class CollectVectorVarMsg : public CMessage_CollectVectorVarMsg
{
public:
  enum DataStatus {VectorValid, FloatVectorValid, BothValid};
public:
  int seq;
  int size;
  DataStatus status;
  AtomID *aid;
  Vector *data;
  FloatVector *fdata;
};

#ifdef MEM_OPT_VERSION
class CollectMidVectorInstance{
  public:

    CollectMidVectorInstance(void) : seq(-10) {
      CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
      ParallelIOMgr *ioMgr = io.ckLocalBranch();
      ioMgr->getMyAtomsRangeOnOutput(fromAtomID, toAtomID);    
    }

    CollectMidVectorInstance(int s) { 
      CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
      ParallelIOMgr *ioMgr = io.ckLocalBranch();
      ioMgr->getMyAtomsRangeOnOutput(fromAtomID, toAtomID);
      reset(s); 
    }

    void free() { seq = -10; }
    int notfree() { return ( seq != -10 ); }

    void reset(int s) {
      if ( s == -10 ) NAMD_bug("seq == free in CollectionMidMaster");
      seq = s;
      remaining = toAtomID-fromAtomID+1;
      data.resize(0);
      fdata.resize(0);
    }

    // if 1 is returned, indicates all the expected atoms for a 
    // timestep have been received
    int append(int size, AtomID *a, Vector *d, FloatVector *fd){      
      if (d) {
        if ( size ) data.resize(toAtomID-fromAtomID+1);
        Vector *ptr = data.begin();
        for(int i = 0; i < size; ++i) { ptr[a[i]-fromAtomID] = d[i]; }
      }
      if (fd) {
        if ( size ) fdata.resize(toAtomID-fromAtomID+1);
        FloatVector *ptr = fdata.begin();
        for(int i = 0; i < size; ++i) { ptr[a[i]-fromAtomID] = fd[i]; }
      }
      remaining -= size;

      return ready();
    }

    int ready(void) { return ( ! remaining ); }

    int seq;    
    Lattice lattice;
    ResizeArray<Vector> data;
    ResizeArray<FloatVector> fdata;

    //indicates the range of atoms this object is responsible for collecting
    int fromAtomID; 
    int toAtomID;

  private:
    int remaining;
};//end of declaration for CollectMidVectorInstance


//An object of this class will be a member of Parallel IO Manager
//It is responsible to buffer the position/coordinates data that is
//going to be written to file system.
//In particular, the instance of this class will be on output procs.
//It will communicate via parallel IO manager with the CollectionMaster
//object on PE0 to be notified on when and what to write --Chao Mei
class CollectionMidMaster{
public:

  CollectionMidMaster(ParallelIOMgr *pIO_) : pIO(pIO_) { parOut = new ParOutput(pIO_->myOutputRank); }
  ~CollectionMidMaster(void) { delete parOut; }

  int receivePositions(CollectVectorVarMsg *msg) {return positions.submitData(msg);}
  int receiveVelocities(CollectVectorVarMsg *msg) {return velocities.submitData(msg);}
  int receiveForces(CollectVectorVarMsg *msg) {return forces.submitData(msg);}

  void disposePositions(int seq);
  void disposeVelocities(int seq);
  void disposeForces(int seq);
  
  CollectMidVectorInstance *getReadyPositions(int seq) { return positions.getReady(seq); }

  //containing an array of CollectVectorInstance and their corresponding
  //timestep value and lattice value
  class CollectVectorSequence{    
  public:
    int submitData(CollectVectorVarMsg *msg){
      int seq = msg->seq;
      CollectMidVectorInstance **c = data.begin();
      CollectMidVectorInstance **c_e = data.end();
      for( ; c != c_e && (*c)->seq != seq; ++c );
      if ( c == c_e ){
        c = data.begin();
        for( ; c != c_e && (*c)->notfree(); ++c );
        if ( c == c_e ) {
          data.add(new CollectMidVectorInstance(seq));
          c = data.end() - 1;
        }
        (*c)->reset(seq);
      }      
      AtomID *i = msg->aid;
      Vector *d = msg->data;
      FloatVector *fd = msg->fdata;
      if(msg->status==CollectVectorVarMsg::VectorValid) {
        fd = NULL;
      }else if(msg->status==CollectVectorVarMsg::FloatVectorValid){
        d = NULL;
      }
      return (*c)->append(msg->size,i,d,fd);
    }

    CollectMidVectorInstance* getReady(int seq){
      CollectMidVectorInstance **c = data.begin();
      CollectMidVectorInstance **c_e = data.end();
      for( ; c != c_e && (*c)->seq != seq; ++c );
      CmiAssert(c != c_e);
      return *c;      
    }

    ResizeArray<CollectMidVectorInstance*> data;
  };//end of declaration for CollectionMidMaster::CollectVectorSequence

private:
  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;
  ParallelIOMgr *pIO; 
  ParOutput *parOut; 
}; //end of declaration for CollectionMidMaster
#endif

#endif

