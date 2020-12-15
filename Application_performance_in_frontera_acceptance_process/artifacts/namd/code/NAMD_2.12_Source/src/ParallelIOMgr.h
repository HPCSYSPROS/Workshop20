#ifndef PARALLELIOMGR_H
#define PARALLELIOMGR_H

#include "ProcessorPrivate.h"
#include "charm++.h"
#include "BOCgroup.h"
#include "common.h"
#include "CompressPsf.h"
#include "Hydrogen.h"
#include "Vector.h"
#include "NamdState.h"
#include "Node.decl.h"
#include "PatchMgr.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "Molecule.h"

#define COLLECT_PERFORMANCE_DATA 0

#ifdef MEM_OPT_VERSION
class CollectionMgr;
class CollectionMaster;
class CollectionMidMaster;
class CollectMidVectorInstance;
#endif

class CollectVectorVarMsg;
class PatchMap;

#include "ParallelIOMgr.decl.h"

///////Beginning of msg declarations related to parallel input////////
class MolInfoMsg : public CMessage_MolInfoMsg
{
public:
    int numBonds;
    int numCalcBonds;
    int numAngles;
    int numCalcAngles;
    int numDihedrals;
    int numCalcDihedrals;
    int numImpropers;
    int numCalcImpropers;
    int numCrossterms;
    int numCalcCrossterms;
    int numExclusions;
    int numCalcExclusions;
    int numCalcFullExclusions;
    int numLJPairs;
    int numCalcLJPairs;
    
    int numRigidBonds;

    //used for "comMov" is false in the SimParameter object
    //which is usually true -Chao Mei
    BigReal totalMass;
    Vector totalMV;

    BigReal totalCharge;
};

class HydroBasedMsg : public CMessage_HydroBasedMsg
{
//Info inside this message needs to be calculated when the hydrogen
//group is constructed
public:
    int numFixedRigidBonds;
    int numFixedGroups;
};

class MoveInputAtomsMsg : public CMessage_MoveInputAtomsMsg
{
public:
    int length;
    InputAtom *atomList;
};

class AtomsCntPerPatchMsg: public CMessage_AtomsCntPerPatchMsg
{
public:
    //one-to-one mapping between two lists
    int length;
    PatchID *pidList;
    unsigned short *atomsCntList;
    unsigned short *fixedAtomsCntList;
};

class MovePatchAtomsMsg: public CMessage_MovePatchAtomsMsg
{
public:
    //the total size of allAtoms could be calculated from sizeList
    int from;
    int patchCnt;
    PatchID *pidList;
    int *sizeList;
    FullAtom *allAtoms;
};
///////End of msg declarations related to parallel input////////

///////Beginning of data struct declarations related to parallel output////////
struct ClusterElem {
    int clusterId;
    int atomsCnt; 

    ClusterElem() : clusterId(-1), atomsCnt(0) {}
    ClusterElem(int cid) : clusterId(cid), atomsCnt(0) {}
    int hash() const {
        return clusterId;
    }
    int operator==(const ClusterElem &a) const {
        return clusterId == a.clusterId;
    }
};
typedef UniqueSet<ClusterElem> ClusterSet;
typedef UniqueSetIter<ClusterElem> ClusterSetIter;

class ClusterSizeMsg : public CMessage_ClusterSizeMsg
{
public:
    int srcRank;
    int clusterId;
    int atomsCnt;
};
typedef ResizeArray<ClusterSizeMsg *> ClusterSizeMsgBuffer;

struct ClusterCoorElem{
    int clusterId;    
    Vector dsum;

    ClusterCoorElem(): clusterId(-1), dsum(0.0) {}
    ClusterCoorElem(int cid): clusterId(cid), dsum(0.0) {}
    int hash() const {
        return clusterId;
    }
    int operator==(const ClusterCoorElem &a) const {
        return clusterId == a.clusterId;
    }
};
typedef UniqueSet<ClusterCoorElem> ClusterCoorSet;
typedef UniqueSetIter<ClusterCoorElem> ClusterCoorSetIter;

class ClusterCoorMsg : public CMessage_ClusterCoorMsg
{
public:
    int srcRank;
    int clusterId;    
    Vector dsum;    
};
typedef ResizeArray<ClusterCoorMsg *> ClusterCoorMsgBuffer;
///////End of data struct declarations related to parallel output////////

class ParallelIOMgr : public CBase_ParallelIOMgr
{
#ifdef MEM_OPT_VERSION
    friend class CollectionMgr;
    friend class CollectionMaster;
    friend class CollectionMidMaster;
    friend class CollectMidVectorInstance;
#endif

private:
    SimParameters *simParameters;
    Molecule *molecule;

///////Beginning of fields related to parallel input////////
    int numInputProcs;
    int *inputProcArray;
    //the index to the inputProcArray i.e.
    //inputProcArray[myInputRank] == CkMyPe();
    //if it is not a input proc, the rank is -1;
    int myInputRank;

    //Initially this atom list contains the initially assigned
    //atoms, later it will contain atoms from other input processors
    //based on the migration group
    InputAtomList initAtoms;

    //This atom list contains the migrated atoms from other input
    //processors based on the migration group. Once the migration
    //finishes, atoms in this list is added to initAtoms, and its
    //space is freed.
    InputAtomList tmpRecvAtoms;


    //This variable indicates whether this processor is ready
    //to receive atoms of HomePatches to be created on this
    //processor
    bool isOKToRecvHPAtoms;
    FullAtomList *hpAtomsList;
    ResizeArray<int> hpIDList; //in increasing order

    //tmp variables
    int procsReceived; //used at updateMolInfo and recvAtomsCntPerPatch
    int hydroMsgRecved; //used at recvHydroBasedCounter
    Vector totalMV; //used to remove center of mass motion
    BigReal totalMass; //used to remove center of mass motion
    BigReal totalCharge;
///////End of fields related to parallel input////////

///////Beginning of fields related to parallel output////////
    int numOutputProcs;
    int *outputProcArray;
    char *outputProcFlags;
    //the index to the outputProcArray i.e.
    //outputProcArray[myOutputRank] == CkMyPe();
    //if it is not a output proc, the rank is -1;
    int myOutputRank;
    //the number of simutaneous writers 
    //output procs with rank distance of numOutputProcs/numOutputWrts do the
    //output at a time
    int numOutputWrts;

    //both arrays are of size #local atoms on this output proc
    int *clusterID;
    int *clusterSize;
    //record the number of atoms that a remote cluster has on this
    //output processor
    ClusterSet remoteClusters;
    //record the number of clusters that have atoms on other output procs
    //on this output proc. Should be remoteClusters.size();
    int numRemoteClusters;
    //TEMP var to indicate how many msgs from remote proc have been recved
    //for updating cluster sizes in my local repository (linked with 
    //numRemoteClusters.
    int numCSMAck;

    ClusterSizeMsgBuffer csmBuf; //used to buffer cluster size msgs
    //record the number of remote cluster info queries for this output proc.
    //i.e. SUM(for a particular cluster on this local output proc, 
    //the number of remote output procs that has some atoms belonging 
    //to this cluster). Should be csmBuf.size();  
    int numRemoteReqs;
    //TEMP var to record the number of remote cluster info queries that
    //has received (linked with numRemoteReqs)
    int numReqRecved;

    //used to store the caculated centralized coordinates for each cluster
    //on this local output proc.
    ClusterCoorSet remoteCoors; //similar to remoteClusters
    ClusterCoorMsgBuffer ccmBuf; //similar to csmBuf but for cluster coor msgs
    Position *tmpCoorCon;     
    //the array is of size #local atoms on this output proc
    char *isWater;

#ifdef MEM_OPT_VERSION
    CollectMidVectorInstance *coorInstance;
    CollectionMidMaster *midCM;
#endif

    CkChareID mainMaster;

///////End of fields related to parallel output////////

#if COLLECT_PERFORMANCE_DATA
    int numFixedAtomLookup;
#endif    

private:
    void readCoordinatesAndVelocity();
    //create atom lists that are used for creating home patch
    void prepareHomePatchAtomList();
    //returns the index in hpIDList which points to pid
    int binaryFindHPID(int pid);

    void readInfoForParOutput();

    void integrateClusterCoor();

    int numMyAtoms(int rank, int numProcs);
    //returns the number of atoms INITIALLY assigned on this input processor
    inline int numInitMyAtomsOnInput() {
        return numMyAtoms(myInputRank, numInputProcs);
    }
    inline int numMyAtomsOnOutput() {
        return numMyAtoms(myOutputRank, numOutputProcs);
    }

    int atomRank(int atomID, int numProcs);
    //returns the rank of the input proc that the atom resides on INITIALLY
    inline int atomInitRankOnInput(int atomID) {
        return atomRank(atomID, numInputProcs);
    }
    inline int atomRankOnOutput(int atomID) {
        return atomRank(atomID, numOutputProcs);
    }

    void getMyAtomsRange(int &lowerIdx, int &upperIdx, int rank, int numProcs);
    //get the range of atoms to be read based on the initial distribution
    //i.e. atoms from [lowerIdx ... upperIdx] are going to be loaded
    inline void getMyAtomsInitRangeOnInput(int &lowerIdx, int &upperIdx) {
        return getMyAtomsRange(lowerIdx,upperIdx,myInputRank,numInputProcs);
    }
    inline void getMyAtomsRangeOnOutput(int &lowerIdx, int &upperIdx) {
        return getMyAtomsRange(lowerIdx,upperIdx,myOutputRank,numOutputProcs);
    }
    inline void getAtomsRangeOnOutput(int &lowerIdx, int &upperIdx, int rank) {
        return getMyAtomsRange(lowerIdx,upperIdx,rank,numOutputProcs);
    }

    //This function is only valid for the inital distribution of input atoms
    //mySAId: the atom id this input proc starts with
    //regAId: the atom id to look up
    inline int isAtomFixed(int mySAId, int reqAId){
        int localIdx = reqAId-mySAId;
        if(localIdx>=0 && localIdx<initAtoms.size()){
            //atom "thisAId" is on this input proc now!
            return initAtoms[localIdx].atomFixed;            
        } else{
        #if COLLECT_PERFORMANCE_DATA
            numFixedAtomLookup++;
        #endif
            //atom "thisAId" is NOT on this input proc now!
            return molecule->is_atom_fixed(reqAId);
        }
    }   
    
    //This function returns the highest rank of this output group procs
    inline int getMyOutputGroupHighestRank(){
        int step = numOutputProcs/numOutputWrts;
        int remains = numOutputProcs%numOutputWrts;
        //so "remains" output groups contain "step+1" output procs;
        //"numOutputWrts-remains" output groups contain "step" output procs;
        int limit = (step+1)*remains;
        if(myOutputRank<limit){
            int idx = myOutputRank/(step+1);
            return (idx+1)*(step+1)-1;
        }else{
            int idx = (myOutputRank-limit)/step;
            return limit+(idx+1)*step-1;
        }
    }
public:
    ParallelIOMgr();
    ~ParallelIOMgr();

    void initialize(Node *node);

    //read per-atom files including the binary compressed psf
    //file, coordinate and velocity files
    void readPerAtomInfo();

    //migrate the initally assigned atoms to appropriate input processors
    //based on the migration group
    void migrateAtomsMGrp();
    void recvAtomsMGrp(MoveInputAtomsMsg *msg);
    void integrateMigratedAtoms();

    //Reduce counters for Tuples and Exclusions in Molecule globally
    void updateMolInfo();
    void recvMolInfo(MolInfoMsg *msg);
    void bcastMolInfo(MolInfoMsg *msg);
    void recvHydroBasedCounter(HydroBasedMsg *msg);
    void bcastHydroBasedCounter(HydroBasedMsg *msg);

    //calculate #atoms in each patch and reduce to proc 0
    void calcAtomsInEachPatch();
    void recvAtomsCntPerPatch(AtomsCntPerPatchMsg *msg);

    //distribute atoms to their homepatch processors
    CthThread sendAtomsThread;
    void sendAtomsToHomePatchProcs();
    int numAcksOutstanding;
    void ackAtomsToHomePatchProcs();
    void recvAtomsToHomePatchProcs(MovePatchAtomsMsg *msg);

    //create home patches on this processor
    void createHomePatches();

    //free the space occupied by atoms' names etc.
    void freeMolSpace();

    //used in parallel IO output
    int getNumOutputProcs() { return numOutputProcs; }
    bool isOutputProcessor(int pe);

    void recvClusterSize(ClusterSizeMsg *msg);
    void integrateClusterSize();
    void recvFinalClusterSize(ClusterSizeMsg *msg);

    void receivePositions(CollectVectorVarMsg *msg);
    void receiveVelocities(CollectVectorVarMsg *msg);
    void receiveForces(CollectVectorVarMsg *msg);
    void disposePositions(int seq, double prevT);
    void disposeVelocities(int seq, double prevT);
    void disposeForces(int seq, double prevT);

    void wrapCoor(int seq, Lattice lat);
    void recvClusterCoor(ClusterCoorMsg *msg);
    void recvFinalClusterCoor(ClusterCoorMsg *msg);
};

#endif
