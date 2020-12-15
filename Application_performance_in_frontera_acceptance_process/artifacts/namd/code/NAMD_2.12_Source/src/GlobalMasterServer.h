/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* A GlobalMasterServer is responsible for keeping track of all
   of the ComputeGlobalMasters present on the master node.  It must
   relay atom data to them, and pass messages to the nodes on their
   behalf. */

#ifndef GLOBALMASTERSERVER_H
#define GLOBALMASTERSERVER_H

#include "Lattice.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeMgr;

class GlobalMasterServer {
 public:
  /* initializes this to be a GlobalMasterServer that has no
     masters to serve yet.  GlobalMasterServer will wait for
     <theNumDataSenders> data messages before responding */
  GlobalMasterServer(ComputeMgr *m, int theNumDataSenders);

  virtual ~GlobalMasterServer();

  /* passes atom coordinate data to the GlobalMasters once this has
     been called by each of the ComputeGlobals, and potentially
     generates response messages. */
  void recvData(ComputeGlobalDataMsg *);
  
  /* gives this control over <newClient> */
  void addClient(GlobalMaster *newClient);
 private:
  int forceSendEnabled; // are total forces received?
  int numDataSenders; // the number of expected messages each cycle
  int numForceSenders; // the number of expected force messages each cycle
  int latticeCount; // is lattice received so far this cycle
  int recvCount; // the number of messages so far this cycle
  int firstTime; // used to be compatible with the ComputeGlobals
  int totalAtomsRequested; // the total number of atoms requested
                           // (initially zero)
  int totalGroupsRequested; // the total number of groups requested

  /* the receivedAtomIDs and receivedAtomPositions lists give
     correspond to each other: element i of the receivedAtomIDs is the
     ID of an atom that has position given by element i of the
     receivedAtomPositions. The receivedForceIDs and receivedTotalForces
     lists have similar relationship. This data is built up as messages are
     received, and cleared after being passed off to the Masters. */
  AtomIDList receivedAtomIDs;
  PositionList receivedAtomPositions;
  PositionList receivedGroupPositions; // the group positions
  BigRealList receivedGroupMasses; // the group positions
  ForceList receivedGroupTotalForces;
  AtomIDList receivedForceIDs;
  ForceList receivedTotalForces;

  int step;  // current timestep received from patches
  Lattice lattice;  // current lattice received from patches

  /* the compute manager responsible for my message delivery */
  ComputeMgr *myComputeManager;

  /* the list of global compute masters that this server is
     responsible for serving */
  ResizeArray<GlobalMaster *> clientList;

  /* passes atom data to the clients and generates responses.  The
   first time we just get the requested atom ids from the clients and
   send messages to the ComputeGlobals requesting those atoms, so we
   can have coordinates ready for the first time step.  All future
   times this is called we ask the masters for forces, as well.  XXX
   this is a little weird...why should the first timestep be special?
   Only because we are dealing with the unwritten rule made by the
   ComputeGlobals that they must not be sent forces before they are
   "configured" */
  int callClients();

  /* puts all the info requested by the clients into a single list */
  void resetAtomList(AtomIDList &atomsRequested);
  void resetForceList(AtomIDList &atomsForced, ForceList &forces,
		      ForceList &groupforces);

  /* the group list is ugly - it is, as far as I can tell, a list of
     the atoms in each group, separated by -1s.  So if the Masters
     request groups {1,2,3} and {2,3,4}, <groupsRequested> will be set
     to the list (1,2,3,-1,2,3,4,-1).  The number of groups sent is
     stored in the variable <numGroups>*/
  void resetGroupList(AtomIDList &groupsRequested, int *numGroups);

  /* stores the version of forces that was found by resetForceList */
  AtomIDList lastAtomsForced;
  ForceList lastForces;

};

#endif

