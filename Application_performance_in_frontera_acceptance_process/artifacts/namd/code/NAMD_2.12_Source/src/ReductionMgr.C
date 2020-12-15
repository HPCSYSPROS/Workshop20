/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   The order of execution is expected to be:
	    0. instantiate object
   ------------------ (processing barrier)
   (mode 0) 1. register() and subscribe()
   ------------------
   (mode 1) 2. submit() and request()
   ------------------
   (mode 2) 3. unregister() and unsubscribe()
   ------------------ (processing barrier)
            4. destroy object
   Doing this out-of-order will cause errors.

   Assumes that *only* one thread will require() a specific sequence's data.
*/

#include <stdlib.h>
#include <stdio.h>

#include "InfoStream.h"
#include "PatchMap.h"	// for patchMap

#include "Node.h"
#include "SimParameters.h"

#include "ReductionMgr.decl.h"
#include "ReductionMgr.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

// Used to register and unregister reductions to downstream nodes
class ReductionRegisterMsg : public CMessage_ReductionRegisterMsg {
public:
  int reductionSetID;
  int dataSize;
  int sourceNode;
};

// Used to send reduction data to downstream nodes
class ReductionSubmitMsg : public CMessage_ReductionSubmitMsg {
public:
  int reductionSetID;
  int sourceNode;
  int sequenceNumber;
  int dataSize;
  BigReal *data;
};

ReductionSet::ReductionSet(int setID, int size, int numChildren) {
  if ( setID == REDUCTIONS_BASIC || setID == REDUCTIONS_AMD ) {
    if ( size != -1 ) {
      NAMD_bug("ReductionSet size specified for REDUCTIONS_BASIC or REDUCTIONS_AMD.");
    }
    size = REDUCTION_MAX_RESERVED;
  }
  if ( size == -1 ) NAMD_bug("ReductionSet size not specified.");
  dataSize = size;
  reductionSetID = setID;
  nextSequenceNumber = 0;
  submitsRegistered = 0;
  dataQueue = 0;
  requireRegistered = 0;
  threadIsWaiting = 0;
  addToRemoteSequenceNumber = new int[numChildren];
}

ReductionSet::~ReductionSet() {

  ReductionSetData *current = dataQueue;

  while ( current ) {
    ReductionSetData *next = current->next;
    delete current;
    current = next;
  }
  delete [] addToRemoteSequenceNumber;
}

// possibly create and return data for a particular seqNum
ReductionSetData* ReductionSet::getData(int seqNum) {

  ReductionSetData **current = &dataQueue;

  while ( *current ) {
    if ( (*current)->sequenceNumber == seqNum ) return *current;
    current = &((*current)->next);
  }

//iout << "seq " << seqNum << " created on " << CkMyPe() << "\n" << endi;
  *current = new ReductionSetData(seqNum, dataSize);
  return *current;
}

// possibly delete data for a particular seqNum
ReductionSetData* ReductionSet::removeData(int seqNum) {

  ReductionSetData **current = &dataQueue;

  while ( *current ) {
    if ( (*current)->sequenceNumber == seqNum ) break;
    current = &((*current)->next);
  }

  if ( ! *current ) { NAMD_die("ReductionSet::removeData on missing seqNum"); }

  ReductionSetData *toremove = *current;
  *current = (*current)->next;
  return toremove;
}

void ReductionMgr::buildSpanTree(const int pe, 
                                 const int max_intranode_children,
                                 const int max_internode_children,
                                 int* parent, 
                                 int* num_children, 
                                 int** children)
{
  // If pe is a first-node, children are same-node pes and perhaps some
  // other first-nodes, and parents are other first-nodes. If pe is not a 
  // first-node, build the spanning tree among the children, and the parent
  // is the corresponding first-node

  // No matter what, build list of PEs on my node first
  const int num_pes = CkNumPes();
  const int num_node_pes = CmiNumPesOnPhysicalNode(CmiPhysicalNodeID(pe)); 
  int *node_pes = new int[num_node_pes];
  int pe_index = -1;
  const int first_pe = CmiGetFirstPeOnPhysicalNode(CmiPhysicalNodeID(pe));
  int num_nodes = 0;
  int *node_ids = new int[num_pes];
  int first_pe_index = -1;
  int my_parent_index;
  
  // Make sure PE 0 is a first-node
  if (pe == 0 && first_pe != pe) {
    NAMD_die("PE 0 is not the first physical node. This shouldn't happen");
  }
  // Get all the PEs on my node, and also build the list of all first-nodes
  int i;
  int node_pe_count=0;
  for (i = 0; i < num_pes; i++) {
    // Save first-nodes
    if (CmiGetFirstPeOnPhysicalNode(CmiPhysicalNodeID(i)) == i) {
      node_ids[num_nodes] = i;
      if (i == first_pe)
        first_pe_index = num_nodes;
      num_nodes++;
    }

    // Also, find pes on my node
    const int i1 = (i + first_pe) % num_pes;
    if (CmiPeOnSamePhysicalNode(first_pe,i1)) {
      if ( node_pe_count == num_node_pes )
        NAMD_bug("ReductionMgr::buildSpanTree found inconsistent physical node data from Charm++ runtime");
      node_pes[node_pe_count] = i1;
      if (pe == i1)
        pe_index = node_pe_count;
      node_pe_count++;
    }
  }
  if ( pe_index < 0 || first_pe_index < 0 )
    NAMD_bug("ReductionMgr::buildSpanTree found inconsistent physical node data from Charm++ runtime");
  
  // Any PE might have children on the same node, plus, if its a first-node,
  // it may have several children on other nodes

  int first_loc_child_index = pe_index * max_intranode_children + 1;
  int last_loc_child_index 
    = first_loc_child_index + max_intranode_children - 1;
  if (first_loc_child_index > num_node_pes) {
    first_loc_child_index = num_node_pes;
    last_loc_child_index = num_node_pes;
  } else {
    if (last_loc_child_index >= num_node_pes) 
      last_loc_child_index = num_node_pes-1;
  }
//  CkPrintf("Local [%d] firstpe %d max %d num %d firstloc %d lastloc %d\n",
//           pe,pe_index,max_intranode_children,num_node_pes,
//           first_loc_child_index,last_loc_child_index);
  
  int first_rem_child_index = num_nodes;
  int last_rem_child_index = num_nodes;
  int rem_children=0;
  int *rem_child_index = new int[max_internode_children];
  
  if (first_pe != pe) {
    // I'm not a first_pe, so I have no more children, and my parent
    // is someone else on my node
    my_parent_index = (pe_index-1)/max_intranode_children;
    *parent = node_pes[my_parent_index];
  } else {
    // I am a first_pe, so I may have additional children
    // on other nodes, and my parent will be on another node

    int range_begin = 0;
    int range_end = num_nodes;

    if (pe == 0) {
      my_parent_index = -1;
      *parent = -1;
    } else {
      my_parent_index = 0;
      while ( first_pe_index != range_begin ) {
        my_parent_index = range_begin;
        ++range_begin;
        for ( int i = 0; i < max_internode_children; ++i ) {
          int split = range_begin + ( range_end - range_begin ) / ( max_internode_children - i );
          if ( first_pe_index < split ) { range_end = split; break; } 
          else { range_begin = split; }
        }
      }
      *parent = node_ids[my_parent_index];
    }

    // now we know parent and need only repeat calculation of children
    int prev_child_index = range_begin;
    ++range_begin;
    for ( int i = 0; i < max_internode_children; ++i ) {
      if ( range_begin >= range_end ) break;
      if ( range_begin > prev_child_index ) {
        rem_child_index[rem_children++] = prev_child_index = range_begin;
      }
      range_begin += ( range_end - range_begin ) / ( max_internode_children - i );
    }
  }

  *num_children = 0;
  //CkPrintf("TREE pe %d my_parent %d %d\n",pe,my_parent_index,*parent);

  int loc_children=0;
  if (first_loc_child_index != num_node_pes) {
    loc_children = last_loc_child_index - first_loc_child_index + 1;
    *num_children += loc_children;
//    CkPrintf("TREE pe %d %d local children\n",pe,loc_children);
//  } else {
//    CkPrintf("TREE pe %d No local children\n",pe);
  }

  if (rem_children) {
    *num_children += rem_children;
//    CkPrintf("TREE pe %d %d rem children\n",pe,rem_children);
//  } else {
//    CkPrintf("TREE pe %d No rem children\n",pe);
  }
  if (*num_children == 0)
    *children = 0;
  else {
    *children = new int[*num_children];
//    CkPrintf("TREE pe %d children %d\n",pe,*num_children);
    int k;
    int child=0;
    if (loc_children > 0) {
      for(k=first_loc_child_index; k <= last_loc_child_index; k++) {
//        CkPrintf("TREE pe %d loc child[%d,%d] %d\n",pe,child,k,node_pes[k]);
        (*children)[child++]=node_pes[k];
      }
    }
    if (rem_children > 0) {
      for(k=0; k < rem_children; k++)  {
//        CkPrintf("TREE pe %d rem child[%d,%d] %d\n",pe,child,k,node_ids[rem_child_index[k]]);
        (*children)[child++]=node_ids[rem_child_index[k]];
      }
    }
  }
  delete [] rem_child_index;
  delete [] node_ids;
  delete [] node_pes;
}

// constructor
ReductionMgr::ReductionMgr() {
    if (CkpvAccess(ReductionMgr_instance) == 0) {
      CkpvAccess(ReductionMgr_instance) = this;
    } else {
      DebugM(1, "ReductionMgr::ReductionMgr() - another instance exists!\n");
    }
    
    buildSpanTree(CkMyPe(),REDUCTION_MAX_CHILDREN,REDUCTION_MAX_CHILDREN,
                  &myParent,&numChildren,&children);
    
//    CkPrintf("TREE [%d] parent %d %d children\n",
//      CkMyPe(),myParent,numChildren);
//    if (numChildren > 0) {
//      for(int i=0; i < numChildren; i++)  {
//        CkPrintf("TREE [%d] child %d %d\n",CkMyPe(),i,children[i]);
//      }
//    }
    
    // fill in the spanning tree fields
#if 0  // Old spanning tree
    if (CkMyPe() == 0) {
      myParent = -1;
    } else {
      myParent = (CkMyPe()-1)/REDUCTION_MAX_CHILDREN;
    }
    firstChild = CkMyPe()*REDUCTION_MAX_CHILDREN + 1;
    if (firstChild > CkNumPes()) firstChild = CkNumPes();
    lastChild = firstChild + REDUCTION_MAX_CHILDREN;
    if (lastChild > CkNumPes()) lastChild = CkNumPes();
#endif

    // initialize data
    for(int i=0; i<REDUCTION_MAX_SET_ID; i++) {
      reductionSets[i] = 0;
    }

    DebugM(1,"ReductionMgr() instantiated.\n");
}

// destructor
ReductionMgr::~ReductionMgr() {
    if (children != 0)
      delete [] children;
    for(int i=0; i<REDUCTION_MAX_SET_ID; i++) {
      delete reductionSets[i];
    }

}

// possibly create and return reduction set
ReductionSet* ReductionMgr::getSet(int setID, int size) {
  if ( reductionSets[setID] == 0 ) {
    reductionSets[setID] = new ReductionSet(setID,size,numChildren);
    if ( ! isRoot() ) {
      ReductionRegisterMsg *msg = new ReductionRegisterMsg;
      msg->reductionSetID = setID;
      msg->dataSize = size;
      msg->sourceNode = CkMyPe();
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteRegister(msg);
    }
  } else if ( setID == REDUCTIONS_BASIC || setID == REDUCTIONS_AMD ) {
    if ( size != -1 ) NAMD_bug("ReductionMgr::getSet size set");
  } else if ( size < 0 || reductionSets[setID]->dataSize != size ) {
    NAMD_bug("ReductionMgr::getSet size mismatch");
  }
  return reductionSets[setID];
}

// possibly delete reduction set
void ReductionMgr::delSet(int setID) {
  ReductionSet *set = reductionSets[setID];
  if ( set && ! set->submitsRegistered & ! set->requireRegistered ) {
    if ( ! isRoot() ) {
      ReductionRegisterMsg *msg = new ReductionRegisterMsg;
      msg->reductionSetID = setID;
      msg->sourceNode = CkMyPe();
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteUnregister(msg);
    }
    delete set;
    reductionSets[setID] = 0;
  }
}

// register local submit
SubmitReduction* ReductionMgr::willSubmit(int setID, int size) {
  ReductionSet *set = getSet(setID, size);
  ReductionSetData *data = set->getData(set->nextSequenceNumber);
  if ( data->submitsRecorded ) {
    NAMD_die("ReductionMgr::willSubmit called while reductions outstanding!");
  }

  set->submitsRegistered++;

  SubmitReduction *handle = new SubmitReduction;
  handle->reductionSetID = setID;
  handle->sequenceNumber = set->nextSequenceNumber;
  handle->master = this;
  handle->data = data->data;

  return handle;
}

// unregister local submit
void ReductionMgr::remove(SubmitReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("SubmitReduction deleted while reductions outstanding!");
  }

  set->submitsRegistered--;

  delSet(setID);
}

// local submit
void ReductionMgr::submit(SubmitReduction* handle) {
  int setID = handle->reductionSetID;
  int seqNum = handle->sequenceNumber;
  ReductionSet *set = reductionSets[setID];
  ReductionSetData *data = set->getData(seqNum);

  data->submitsRecorded++;
  if ( data->submitsRecorded == set->submitsRegistered ) {
    mergeAndDeliver(set,seqNum);
  }

  handle->sequenceNumber = ++seqNum;
  handle->data = set->getData(seqNum)->data;
}

// register submit from child
void ReductionMgr::remoteRegister(ReductionRegisterMsg *msg) {

  int setID = msg->reductionSetID;
  int size = msg->dataSize;
  ReductionSet *set = getSet(setID,size);
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("ReductionMgr::remoteRegister called while reductions outstanding on parent!");
  }

  set->submitsRegistered++;
  set->addToRemoteSequenceNumber[childIndex(msg->sourceNode)]
					= set->nextSequenceNumber;
//  CkPrintf("[%d] reduction register received from node[%d] %d\n",
//    CkMyPe(),childIndex(msg->sourceNode),msg->sourceNode);
    
  delete msg;
}

// unregister submit from child
void ReductionMgr::remoteUnregister(ReductionRegisterMsg *msg) {

  int setID = msg->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("SubmitReduction deleted while reductions outstanding on parent!");
  }

  set->submitsRegistered--;

  delSet(setID);
  delete msg;
}

// data submitted from child
void ReductionMgr::remoteSubmit(ReductionSubmitMsg *msg) {
  int setID = msg->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = msg->sequenceNumber
	+ set->addToRemoteSequenceNumber[childIndex(msg->sourceNode)];

//iout << "seq " << seqNum << " from " << msg->sourceNode << " received on " << CkMyPe() << "\n" << endi;
  int size = msg->dataSize;
  if ( size != set->dataSize ) {
    NAMD_bug("ReductionMgr::remoteSubmit data sizes do not match.");
  }

  BigReal *newData = msg->data;
  ReductionSetData *data = set->getData(seqNum);
  BigReal *curData = data->data;
#ifdef ARCH_POWERPC
#pragma disjoint (*curData,  *newData)
#pragma unroll(4)
#endif
  if ( setID == REDUCTIONS_MINIMIZER ) {
    for ( int i = 0; i < size; ++i ) {
      if ( newData[i] > curData[i] ) {
        curData[i] = newData[i];
      }
    }
  } else {
    for ( int i = 0; i < size; ++i ) {
      curData[i] += newData[i];
    }
  }
//  CkPrintf("[%d] reduction Submit received from node[%d] %d\n",
//    CkMyPe(),childIndex(msg->sourceNode),msg->sourceNode);
  delete msg;

  data->submitsRecorded++;
  if ( data->submitsRecorded == set->submitsRegistered ) {
    mergeAndDeliver(set,seqNum);
  }
}

// common code for submission and delivery
void ReductionMgr::mergeAndDeliver(ReductionSet *set, int seqNum) {

//iout << "seq " << seqNum << " complete on " << CkMyPe() << "\n" << endi;
 
    set->nextSequenceNumber++; // should match all clients

    ReductionSetData *data = set->getData(seqNum);
    if ( data->submitsRecorded != set->submitsRegistered ) {
      NAMD_bug("ReductionMgr::mergeAndDeliver not ready to deliver.");
    }

    if ( isRoot() ) {
      if ( set->requireRegistered ) {
	if ( set->threadIsWaiting && set->waitingForSequenceNumber == seqNum) {
	  // awaken the thread so it can take the data
	  CthAwaken(set->waitingThread);
	}
      } else {
	NAMD_die("ReductionSet::deliver will never deliver data");
      }
    } else {
      // send data to parent
      ReductionSubmitMsg *msg = new(set->dataSize) ReductionSubmitMsg;
      msg->reductionSetID = set->reductionSetID;
      msg->sourceNode = CkMyPe();
      msg->sequenceNumber = seqNum;
      msg->dataSize = set->dataSize;
      for ( int i = 0; i < msg->dataSize; ++i ) {
        msg->data[i] = data->data[i];
      }
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteSubmit(msg);
      delete set->removeData(seqNum);
    }

}

// register require
RequireReduction* ReductionMgr::willRequire(int setID, int size) {
  ReductionSet *set = getSet(setID,size);
  set->requireRegistered++;
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("ReductionMgr::willRequire called while reductions outstanding!");
  }

  RequireReduction *handle = new RequireReduction;
  handle->reductionSetID = setID;
  handle->sequenceNumber = set->nextSequenceNumber;
  handle->master = this;

  return handle;
}

// unregister require
void ReductionMgr::remove(RequireReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("RequireReduction deleted while reductions outstanding!");
  }

  set->requireRegistered--;

  delSet(setID);
}

// require the data from a thread
void ReductionMgr::require(RequireReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = handle->sequenceNumber;
  ReductionSetData *data = set->getData(seqNum);
  if ( data->submitsRecorded < set->submitsRegistered ) {
    set->threadIsWaiting = 1;
    set->waitingForSequenceNumber = seqNum;
    set->waitingThread = CthSelf();
//iout << "seq " << seqNum << " waiting\n" << endi;
    CthSuspend();
  }
  set->threadIsWaiting = 0;

//iout << "seq " << seqNum << " consumed\n" << endi;
  delete handle->currentData;
  handle->currentData = set->removeData(seqNum);
  handle->data = handle->currentData->data;
  handle->sequenceNumber = ++seqNum;
}


#include "ReductionMgr.def.h"
// nothing should be placed below here

