/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "CollectionMgr.decl.h"
#include "CollectionMgr.h"
#include "CollectionMaster.decl.h"
#include "CollectionMaster.h"
#include "Node.h"
#include "SimParameters.h"

#include "ParallelIOMgr.decl.h"
#include "ParallelIOMgr.h"

//#define DEBUGM
#include "Debug.h"

CollectionMgr::CollectionMgr(SlaveInitMsg *msg) : master(msg->master)
{
  delete msg;
  if (CkpvAccess(CollectionMgr_instance) == 0) {
    CkpvAccess(CollectionMgr_instance) = this;
  } else {
    DebugM(1, "CollectionMgr::CollectionMgr() - another instance of CollectionMgr exists!\n");
  }
}


CollectionMgr::~CollectionMgr(void)
{
}

#ifdef MEM_OPT_VERSION
//1. record the dest output rank of each atom
//2. distribute the atoms to the corresponding output procs
//1 and 2 are both needed for positions and velocities

void CollectionMgr::submitPositions(int seq, FullAtomList &a,
				Lattice l, int prec)
{
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  ResizeArray<int> oRank(numAtoms);
  PositionList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    oRank[i] = a[i].outputRank;
    d[i] = l.reverse_transform(a[i].position,a[i].transform);
  }
  CollectVectorInstance *c;
  if ( ( c = positions.submitData(seq,aid,oRank,d,prec) ) )
  {    
    CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
    ParallelIOMgr *ioMgr = io.ckLocalBranch();

    //construct per output proc atoms list
    AtomIDList *perOList = new AtomIDList[ioMgr->numOutputProcs];
    for(int i=0; i<c->aid.size(); i++){
        perOList[c->outRank[i]].add(i);
    }    
    CollectVectorVarMsg::DataStatus vstatus;
    if(c->data.size()==0){
        vstatus = CollectVectorVarMsg::FloatVectorValid;
    }else if(c->fdata.size()==0){
        vstatus = CollectVectorVarMsg::VectorValid;
    }else{
        vstatus = CollectVectorVarMsg::BothValid;
    }
    //send msg to output proc if there's one    
    for(int i=0; i<ioMgr->numOutputProcs; i++){
        int numAtoms = perOList[i].size();
        if(!numAtoms) continue;
        CollectVectorVarMsg *msg;
        if( vstatus == CollectVectorVarMsg::VectorValid){
            msg = new(numAtoms, numAtoms, 0, 0)CollectVectorVarMsg;
            for(int j=0; j<numAtoms; j++){
                int lIdx = perOList[i][j];
                msg->aid[j] = c->aid[lIdx];
                msg->data[j] = c->data[lIdx];
            }
        }else if(vstatus == CollectVectorVarMsg::FloatVectorValid){
            msg = new(numAtoms, 0, numAtoms, 0)CollectVectorVarMsg;
            for(int j=0; j<numAtoms; j++){
                int lIdx = perOList[i][j];
                msg->aid[j] = c->aid[lIdx];
                msg->fdata[j] = c->fdata[lIdx];
            }
        }else{
            msg = new(numAtoms, numAtoms, numAtoms, 0)CollectVectorVarMsg;
            for(int j=0; j<numAtoms; j++){
                int lIdx = perOList[i][j];
                msg->aid[j] = c->aid[lIdx];
                msg->data[j] = c->data[lIdx];
                msg->fdata[j] = c->fdata[lIdx];
            }
        }
        msg->seq = c->seq;
        msg->size = numAtoms;
        msg->status = vstatus;
        io[ioMgr->outputProcArray[i]].receivePositions(msg);
    }
    c->free();
    delete [] perOList;
  }
}

void CollectionMgr::submitVelocities(int seq, int zero, FullAtomList &a)
{
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  ResizeArray<int> oRank(numAtoms);
  PositionList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    oRank[i] = a[i].outputRank;
    if ( zero ) d[i] = 0.;
    else d[i] = a[i].velocity;
  }
  CollectVectorInstance *c;
  if ( ( c = velocities.submitData(seq,aid,oRank,d) ) )
  {
      CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
      ParallelIOMgr *ioMgr = io.ckLocalBranch();

      //construct per output proc atoms list
      AtomIDList *perOList = new AtomIDList[ioMgr->numOutputProcs];
      for(int i=0; i<c->aid.size(); i++){
          perOList[c->outRank[i]].add(i);
      }    
      CollectVectorVarMsg::DataStatus vstatus = CollectVectorVarMsg::VectorValid;
      //send msg to output proc if there's one    
        for(int i=0; i<ioMgr->numOutputProcs; i++){
            int numAtoms = perOList[i].size();
            if(!numAtoms) continue;
            CollectVectorVarMsg *msg;            
            msg = new(numAtoms, numAtoms, 0, 0)CollectVectorVarMsg;
            msg->seq = c->seq;
            msg->size = numAtoms;
            msg->status = vstatus;
            for(int j=0; j<numAtoms; j++){
                int lIdx = perOList[i][j];
                msg->aid[j] = c->aid[lIdx];
                msg->data[j] = c->data[lIdx];
            }
            io[ioMgr->outputProcArray[i]].receiveVelocities(msg);            
        }
        c->free();
        delete [] perOList;         
  }
}

void CollectionMgr::submitForces(int seq, FullAtomList &a, int maxForceUsed, ForceList *f)
{
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  ResizeArray<int> oRank(numAtoms);
  ForceList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    oRank[i] = a[i].outputRank;
    d[i] = 0.;
  }
  for ( int j=0; j<=maxForceUsed; ++j ) {
    Force *fptr = f[j].begin();
    for ( int i=0; i<numAtoms; ++i ) {
      d[i] += fptr[i];
    }
  }
  if ( Node::Object()->simParameters->fixedAtomsOn && !Node::Object()->simParameters->fixedAtomsForceOutput) {
    for ( int i=0; i<numAtoms; ++i ) {
      if ( a[i].atomFixed ) d[i] = 0.;
    }
  }
  CollectVectorInstance *c;
  if ( ( c = forces.submitData(seq,aid,oRank,d) ) )
  {
      CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
      ParallelIOMgr *ioMgr = io.ckLocalBranch();

      //construct per output proc atoms list
      AtomIDList *perOList = new AtomIDList[ioMgr->numOutputProcs];
      for(int i=0; i<c->aid.size(); i++){
          perOList[c->outRank[i]].add(i);
      }    
      CollectVectorVarMsg::DataStatus vstatus = CollectVectorVarMsg::VectorValid;
      //send msg to output proc if there's one    
        for(int i=0; i<ioMgr->numOutputProcs; i++){
            int numAtoms = perOList[i].size();
            if(!numAtoms) continue;
            CollectVectorVarMsg *msg;            
            msg = new(numAtoms, numAtoms, 0, 0)CollectVectorVarMsg;
            msg->seq = c->seq;
            msg->size = numAtoms;
            msg->status = vstatus;
            for(int j=0; j<numAtoms; j++){
                int lIdx = perOList[i][j];
                msg->aid[j] = c->aid[lIdx];
                msg->data[j] = c->data[lIdx];
            }
            io[ioMgr->outputProcArray[i]].receiveForces(msg);            
        }
        c->free();
        delete [] perOList;         
  }
}

#else
void CollectionMgr::submitPositions(int seq, FullAtomList &a,
				Lattice l, int prec)
{  
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  PositionList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    d[i] = l.reverse_transform(a[i].position,a[i].transform);
  }
  CollectVectorInstance *c;
  if ( ( c = positions.submitData(seq,aid,d,prec) ) )
  {
    int aid_size = c->aid.size();
    int data_size = c->data.size();
    int fdata_size = c->fdata.size();
    CollectVectorMsg *msg
      = new (aid_size, data_size, fdata_size,0) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid_size = aid_size;
    msg->data_size = data_size;
    msg->fdata_size = fdata_size;
    memcpy(msg->aid,c->aid.begin(),aid_size*sizeof(AtomID));
    memcpy(msg->data,c->data.begin(),data_size*sizeof(Vector));
    memcpy(msg->fdata,c->fdata.begin(),fdata_size*sizeof(FloatVector));
    CProxy_CollectionMaster cm(master);
    cm.receivePositions(msg);
    c->free();
  }
}

void CollectionMgr::submitVelocities(int seq, int zero, FullAtomList &a)
{
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  PositionList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    if ( zero ) d[i] = 0.;
    else d[i] = a[i].velocity;
  }
  CollectVectorInstance *c;
  if ( ( c = velocities.submitData(seq,aid,d) ) )
  {
    int aid_size = c->aid.size();
    int data_size = c->data.size();
    CollectVectorMsg *msg = new (aid_size, data_size, 0, 0) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid_size = aid_size;
    msg->data_size = data_size;
    msg->fdata_size = 0;
    memcpy(msg->aid,c->aid.begin(),aid_size*sizeof(AtomID));
    memcpy(msg->data,c->data.begin(),data_size*sizeof(Vector));
    CProxy_CollectionMaster cm(master);
    cm.receiveVelocities(msg);
    c->free();
  }
}
 
void CollectionMgr::submitForces(int seq, FullAtomList &a, int maxForceUsed, ForceList *f)
{
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  ForceList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    d[i] = 0.;
  }
  for ( int j=0; j<=maxForceUsed; ++j ) {
    Force *fptr = f[j].begin();
    for ( int i=0; i<numAtoms; ++i ) {
      d[i] += fptr[i];
    }
  }
  if ( Node::Object()->simParameters->fixedAtomsOn && !Node::Object()->simParameters->fixedAtomsForceOutput) {
    for ( int i=0; i<numAtoms; ++i ) {
      if ( a[i].atomFixed ) d[i] = 0.;
    }
  }
  CollectVectorInstance *c;
  if ( ( c = forces.submitData(seq,aid,d) ) )
  {
    int aid_size = c->aid.size();
    int data_size = c->data.size();
    CollectVectorMsg *msg = new (aid_size, data_size, 0, 0) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid_size = aid_size;
    msg->data_size = data_size;
    msg->fdata_size = 0;
    memcpy(msg->aid,c->aid.begin(),aid_size*sizeof(AtomID));
    memcpy(msg->data,c->data.begin(),data_size*sizeof(Vector));
    CProxy_CollectionMaster cm(master);
    cm.receiveForces(msg);
    c->free();
  }
}
#endif

void CollectionMgr::sendDataStream(const char *data) {
  DataStreamMsg *msg = new DataStreamMsg;
  msg->data.resize(strlen(data)+1);
  strcpy(msg->data.begin(),data);
  CProxy_CollectionMaster cm(master);
  cm.receiveDataStream(msg);
}

#include "CollectionMgr.def.h"

