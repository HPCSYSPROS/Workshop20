/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "ComputeGlobalMsgs.h"
#include "ComputeMgr.h"
#include "NamdTypes.h"
#include "InfoStream.h"
#include "Node.h"
#include "SimParameters.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"
#include <stdlib.h>
#include <vector>
#include <algorithm>

using namespace std;
struct position_index {
  int atomID;
  int index;
};
void GlobalMasterServer::addClient(GlobalMaster *newClient) {
  DebugM(3,"Adding client\n");
  clientList.add(newClient);
  newClient->setLattice(&lattice);
  DebugM(2,"Added.\n");
}

void GlobalMasterServer::recvData(ComputeGlobalDataMsg *msg) {
  DebugM(3,"Storing data (" << msg->aid.size() << " positions) on master\n");

  if ( msg->step != -1 ) step = msg->step;

  /* get the beginning and end of the lists */
  AtomIDList::iterator a_i = msg->aid.begin();
  AtomIDList::iterator a_e = msg->aid.end();
  PositionList::iterator p_i = msg->p.begin();
  PositionList::iterator g_i = msg->gcom.begin();
  PositionList::iterator g_e = msg->gcom.end();
  BigRealList::iterator gm_i = msg->gmass.begin();

  /* iterate over each member of the atom lists */
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    receivedAtomIDs.add(*a_i);
    receivedAtomPositions.add(*p_i);
  }
  
  /* iterate over each member of "total force" lists */
  a_e = msg->fid.end();
  ForceList::iterator f_i=msg->tf.begin();
  for (a_i=msg->fid.begin() ; a_i!=a_e; ++a_i,++f_i) {
    receivedForceIDs.add(*a_i);
    receivedTotalForces.add(*f_i);
  }

  /* iterate over each member of the group position list */
  int i=0;
  for ( ; g_i != g_e; ++g_i, ++gm_i ) {
    DebugM(1,"Received center of mass "<<*g_i<<"\n");
    if(i >= totalGroupsRequested) NAMD_bug("Received too many groups.");
    receivedGroupPositions[i] += (*g_i);
    receivedGroupMasses[i] += (*gm_i);
    i++;
  }
  if(i!=totalGroupsRequested) NAMD_bug("Received too few groups.");

  /* iterate over each member of group total force lists */
  int ntf = msg->gtf.size();
  if ( ntf && ntf != receivedGroupTotalForces.size() ) NAMD_bug("Received wrong number of group forces.");
  ForceList::iterator gf_i=msg->gtf.begin();
  ForceList::iterator gf_e=msg->gtf.end();
  for ( i=0 ; gf_i != gf_e; ++gf_i, ++i ) {
    receivedGroupTotalForces[i] += (*gf_i);
  }

  if ( msg->lat.size() ) {
    if ( latticeCount ) {
      NAMD_bug("GlobalMasterServer::recvData received lattice twice.");
    }
    lattice = msg->lat[0];
    latticeCount = 1;
  }

  recvCount += msg->count;

  /* done with the message, delete it */
  delete msg;

  /* check whether we've gotten all the expected messages */
  if(recvCount > numDataSenders + numForceSenders + 1) {
    NAMD_bug("GlobalMasterServer::recvData recvCount too high.");
  }
  if(recvCount == numDataSenders + numForceSenders + 1) {
    if ( ! latticeCount ) {
      NAMD_bug("GlobalMasterServer::recvData did not receive lattice.");
    }

    int oldTotalGroupsRequested = totalGroupsRequested;

    DebugM(3,"received messages from each of the ComputeGlobals\n");
    int resendCoordinates = callClients();

    /* now restart */
    step = -1;
    receivedAtomIDs.resize(0);
    receivedAtomPositions.resize(0);
    receivedGroupPositions.resize(totalGroupsRequested);
    receivedGroupPositions.setall(Vector(0,0,0));
    receivedGroupMasses.resize(totalGroupsRequested);
    receivedGroupMasses.setall(0);
    latticeCount = 0;
    recvCount = 0;
    if ( resendCoordinates ) {
      recvCount += numForceSenders;
    } else {
      receivedGroupTotalForces.resize(totalGroupsRequested);
      receivedGroupTotalForces.setall(0);
      receivedForceIDs.resize(0);
      receivedTotalForces.resize(0);
    }
  }
}

void GlobalMasterServer::resetAtomList(AtomIDList &atomsRequested) {
  atomsRequested.resize(0);

  /* iterate over all of the masters */
  GlobalMaster **m_i = clientList.begin();
  GlobalMaster **m_e = clientList.end();
  while(m_i != m_e) {
    /* add all of the atoms in this master */
    int i;
    for(i=0;i<(*m_i)->requestedAtoms().size();i++) {
      atomsRequested.add((*m_i)->requestedAtoms()[i]);
    }

    /* go to next master */
    m_i++;
  }

  /* remove duplicates */
  sort(atomsRequested.begin(), atomsRequested.end());
  AtomIDList::iterator it = unique(atomsRequested.begin(), atomsRequested.end());
  atomsRequested.resize( distance(atomsRequested.begin(), it) );
}

void GlobalMasterServer::resetGroupList(AtomIDList &groupsRequested,
					int *numGroups) {
  DebugM(3,"Rebuilding the group list\n");
  groupsRequested.resize(0);
  *numGroups = 0;

  /* iterate over all of the masters */
  GlobalMaster **m_i = clientList.begin();
  GlobalMaster **m_e = clientList.end();
  while(m_i != m_e) {

    /* add all of the groups requested by this master */
    int i;
    GlobalMaster *master = *m_i;
    for(i=0;i<master->requestedGroups().size();i++) {
      /* add all of the atoms in this group, then add a -1 */
      int j;
      const AtomIDList &atoms_in_group = master->requestedGroups()[i];
      (*numGroups) ++;
      DebugM(1,"adding group "<<*numGroups<<"\n");
      for(j=0;j<atoms_in_group.size();j++) {
	groupsRequested.add(atoms_in_group[j]); // add an atom
      }
      DebugM(1,"here\n");
      groupsRequested.add(-1); // add a -1 to separate the groups (yuck)
    }

    /* go to next master */
    m_i++;
  }
}

void GlobalMasterServer::resetForceList(AtomIDList &atomsForced,
					ForceList &forces,
					ForceList &groupForces) {
  DebugM(1,"Restting forces\n");
  atomsForced.resize(0);
  forces.resize(0);
  groupForces.resize(0);
  lastAtomsForced.resize(0);
  lastForces.resize(0);

  /* iterate over all of the masters */
  GlobalMaster **m_i = clientList.begin();
  GlobalMaster **m_e = clientList.end();
  while(m_i != m_e) {
    (*m_i)->check(); // just in case!

    /* add all of the atoms in this master */
    int i;
    DebugM(1,"Adding atom forces\n");
    for(i=0;i<(*m_i)->forcedAtoms().size();i++) {
      atomsForced.add((*m_i)->forcedAtoms()[i]);
      forces.add((*m_i)->appliedForces()[i]);
      lastAtomsForced.add((*m_i)->forcedAtoms()[i]);
      lastForces.add((*m_i)->appliedForces()[i]);
    }

    /* add all of the group forces for this master */
    DebugM(1,"Adding "<<(*m_i)->requestedGroups().size()<<" group forces\n");
    for(i=0;i<(*m_i)->groupForces().size();i++) {
      groupForces.add((*m_i)->groupForces()[i]);
    }

    /* go to next master */
    DebugM(1,"Next master...\n");
    m_i++;
  }
  DebugM(1,"Done restting forces\n");
}

struct atomID_less {
      bool operator ()(position_index const& a, position_index const& b) const {
        if (a.atomID < b.atomID) return true;
        if (a.atomID > b.atomID) return false;

        return false;
}
};

int GlobalMasterServer::callClients() {
  DebugM(3,"Calling clients\n");
  bool forceSendActive = false;
  {
    GlobalMaster **m_i = clientList.begin();
    GlobalMaster **m_e = clientList.end();
    for ( ; m_i != m_e; ++m_i ) {
      if ( (*m_i)->changedAtoms() ) firstTime = 1;
      if ( (*m_i)->changedGroups() ) firstTime = 1;
      if ( (*m_i)->requestedTotalForces() ) forceSendActive = true;
      (*m_i)->clearChanged();
    }
  }

  if(firstTime) {
 
    /* the first time we just get the requested atom ids from the
       clients and send messages to the compute globals requesting
       those atoms, so they can have coordinates for the first time
       step. */
    DebugM(1,"first time.\n");
    ComputeGlobalResultsMsg *msg = new ComputeGlobalResultsMsg;
    resetAtomList(msg->newaid); // add any atom IDs made in constructors
    // resetForceList(msg->aid,msg->f,msg->gforce); // same for forces
    resetGroupList(msg->newgdef,&totalGroupsRequested);
    msg->resendCoordinates = 1;
    msg->reconfig = 1;
    msg->totalforces = forceSendActive;
    totalAtomsRequested = msg->newaid.size(); // record the atom total

    numDataSenders = totalAtomsRequested;
    AtomIDList::iterator g_i = msg->newgdef.begin();
    AtomIDList::iterator g_e = msg->newgdef.end();
    for ( ; g_i != g_e; ++g_i ) {
      if ( *g_i != -1 ) ++numDataSenders;
    }

    DebugM(3,"Sending configure ("<<totalAtomsRequested<<" atoms, "
	   <<totalGroupsRequested<<" groups)\n");
    myComputeManager->sendComputeGlobalResults(msg);

    firstTime = 0;
    return 1;
  }

  /* check to make sure we've got everything */
  if(receivedAtomIDs.size() != totalAtomsRequested) {
    DebugM(3,"Requested " << totalAtomsRequested << " atoms.\n");
    NAMD_bug("Got the wrong number of atoms");
  }
  if(receivedGroupPositions.size() != totalGroupsRequested) {
    DebugM(3,"Requested " << totalGroupsRequested << " groups.\n");
    DebugM(3,"Got " << receivedGroupPositions.size() << " group positions.\n");
    NAMD_bug("Got the wrong number of groups");
  }
  if(receivedGroupMasses.size() != totalGroupsRequested) {
    DebugM(3,"Requested " << totalGroupsRequested << " groups.\n");
    DebugM(3,"Got " << receivedGroupMasses.size() << " group masses.\n");
    NAMD_bug("Got the wrong number of group masses");
  }

  /* get the beginning and end of the lists */
  AtomIDList::iterator a_i = receivedAtomIDs.begin();
  AtomIDList::iterator a_e = receivedAtomIDs.end();
  PositionList::iterator p_i = receivedAtomPositions.begin();
  PositionList::iterator g_i = receivedGroupPositions.begin();
  PositionList::iterator g_e = receivedGroupPositions.end();
  ForceList::iterator gtf_i = receivedGroupTotalForces.begin();
  ForceList::iterator gtf_e = receivedGroupTotalForces.end();
  AtomIDList::iterator forced_atoms_i = lastAtomsForced.begin();
  AtomIDList::iterator forced_atoms_e = lastAtomsForced.end();
  ForceList::iterator forces_i = lastForces.begin();
  GlobalMaster **m_i = clientList.begin();
  GlobalMaster **m_e = clientList.end();
    AtomIDList::iterator f_i = receivedForceIDs.begin();
    AtomIDList::iterator f_e = receivedForceIDs.end();

  /* renormalize group center of mass by group total mass */
  ResizeArray<BigReal>::iterator gm_i = receivedGroupMasses.begin();
  for ( ; g_i < g_e; ++g_i, ++gm_i ) { *g_i /= *gm_i; }
  g_i = receivedGroupPositions.begin();
  gm_i = receivedGroupMasses.begin();

  /* use these to check whether anything has changed for any master */
  bool requested_atoms_changed=false;
  bool requested_forces_changed=false;
  bool requested_groups_changed=false;
  forceSendActive = false;
  
  vector <position_index> positions;
  for (int j = 0; a_i != a_e; ++a_i, ++j) {
    position_index temp;
    temp.atomID = *a_i;
    temp.index = j;
    positions.push_back(temp); 
  }
  sort(positions.begin(), positions.end(), atomID_less());

  /* call each of the masters with the coordinates */
  while(m_i != m_e) {
    int num_atoms_requested, num_groups_requested;
    
    /* get the masters information */
    GlobalMaster *master = *m_i;
    num_atoms_requested = master->requestedAtoms().size();
    num_groups_requested = master->requestedGroups().size();

    AtomIDList   clientAtomIDs;
    PositionList clientAtomPositions;
    AtomIDList   clientReceivedForceIDs;
    ForceList    clientReceivedTotalForces;

    if (num_atoms_requested) {
      vector <int> rqAtoms;
      for (int i = 0; i < master->requestedAtoms().size(); i++){
        rqAtoms.push_back(master->requestedAtoms()[i]);
      }
      sort(rqAtoms.begin(), rqAtoms.end());
      int j = 0;
      for (int i = 0; i < positions.size(); i++){
        if (positions[i].atomID == rqAtoms[j]){
          clientAtomPositions.add(receivedAtomPositions[positions[i].index]);
          clientAtomIDs.add(rqAtoms[j]);
          --i;  // rqAtoms may contain repeats
          if ( ++j == num_atoms_requested ) break;
        }
      }
      if ( j != num_atoms_requested ) NAMD_bug(
        "GlobalMasterServer::callClients() did not find all requested atoms");
    }

    AtomIDList::iterator ma_i = clientAtomIDs.begin();
    AtomIDList::iterator ma_e = clientAtomIDs.end();
    PositionList::iterator mp_i = clientAtomPositions.begin();
    AtomIDList::iterator mf_i = clientReceivedForceIDs.begin();
    AtomIDList::iterator mf_e = clientReceivedForceIDs.end();
    ForceList::iterator mtf_i = clientReceivedTotalForces.begin();

    /* check to make sure we have some atoms left.  This must work for
     zero requested atoms, as well! */
   // if(a_i+num_atoms_requested > a_e)
   //   NAMD_bug("GlobalMasterServer ran out of atom IDs!");

    /* update this master */
    master->step = step;
    master->processData(ma_i,ma_e,
			mp_i,g_i,g_i+num_groups_requested,
			gm_i,gm_i+num_groups_requested,
			gtf_i,gtf_i+(numForceSenders?master->old_num_groups_requested:0),
			forced_atoms_i,forced_atoms_e,forces_i,
      receivedForceIDs.begin(),receivedForceIDs.end(),receivedTotalForces.begin());

    a_i = receivedAtomIDs.begin();
    p_i = receivedAtomPositions.begin();
//}

    /* check to see if anything changed */
    if(master->changedAtoms()) {
      requested_atoms_changed = true;
    }
    if(master->changedForces()) {
      requested_forces_changed = true;
    }
    if(master->changedGroups()) {
      requested_groups_changed = true;
    }
    master->clearChanged();
    if(master->requestedTotalForces()) forceSendActive = true;

    /* go to next master */
    m_i++;

    g_i += num_groups_requested;
    gm_i += num_groups_requested;
    if ( numForceSenders ) gtf_i += master->old_num_groups_requested;
    master->old_num_groups_requested = master->requestedGroups().size();  // include changes
  } 

  if ( numForceSenders && gtf_i != gtf_e ) NAMD_bug("GlobalMasterServer::callClients bad group total force count");

  /* make a new message */
  ComputeGlobalResultsMsg *msg = new ComputeGlobalResultsMsg;

  /* build an atom list, if necessary */
  if(requested_atoms_changed || requested_groups_changed) {
    resetAtomList(msg->newaid); // add all of the atom IDs
    totalAtomsRequested = msg->newaid.size();
    msg->reconfig = 1; // request a reconfig
    resetGroupList(msg->newgdef,&totalGroupsRequested); // add all of the group IDs
    numDataSenders = totalAtomsRequested;
    AtomIDList::iterator g_i = msg->newgdef.begin();
    AtomIDList::iterator g_e = msg->newgdef.end();
    for ( ; g_i != g_e; ++g_i ) {
      if ( *g_i != -1 ) ++numDataSenders;
    }
  }
  msg->totalforces = forceSendActive;
  numForceSenders = (forceSendActive ? numDataSenders : 0);
  resetForceList(msg->aid,msg->f,msg->gforce); // could this be more efficient?

  /* get group acceleration by renormalizing group net force by group total mass */
  ForceList::iterator gf_i = msg->gforce.begin();
  ForceList::iterator gf_e = msg->gforce.end();
  gm_i = receivedGroupMasses.begin();
  for ( ; gf_i < gf_e; ++gf_i, ++gm_i ) { *gf_i /= *gm_i; }
  
  /* now send the results */
  DebugM(3,"Sending results ("<<totalAtomsRequested<<" atoms requested, "
	 <<totalGroupsRequested<<" groups requested, "
	 <<msg->f.size()<<" forces set)\n");
  myComputeManager->sendComputeGlobalResults(msg);
  DebugM(3,"Sent.\n");
  return 0;
}

GlobalMasterServer::GlobalMasterServer(ComputeMgr *m,
				       int theNumDataSenders) {
  DebugM(3,"Constructing GlobalMasterServer\n");
  myComputeManager = m;
  numDataSenders = 0; // theNumDataSenders;
  numForceSenders = 0;
  latticeCount = 0;
  lattice = Node::Object()->simParameters->lattice;
  recvCount = 0; /* we haven't gotten any messages yet */
  firstTime = 1; /* XXX temporary */
  step = -1;
  totalAtomsRequested = 0;
  totalGroupsRequested = 0;
  forceSendEnabled = 0;
  if ( Node::Object()->simParameters->tclForcesOn ) forceSendEnabled = 1;
  if ( Node::Object()->simParameters->colvarsOn ) forceSendEnabled = 1;
  DebugM(3,"done constructing\n");
}

GlobalMasterServer::~GlobalMasterServer() {
  GlobalMaster *m_i = *clientList.begin();
  GlobalMaster *m_e = *clientList.end();
  
  /* delete each of the masters */
  while(m_i != m_e) {
    delete m_i;
    m_i++;
  }
}


 
