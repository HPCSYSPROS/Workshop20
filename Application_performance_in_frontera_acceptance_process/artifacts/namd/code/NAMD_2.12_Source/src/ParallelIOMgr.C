
#include "largefiles.h"  // must be first!

#include <stdio.h>
#include "BOCgroup.h"
#include "Molecule.h"
#include "Node.h"
#include "Node.decl.h"
#include "NamdState.h"
#include "WorkDistrib.h"
#include "PDB.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "packmsg.h"
#include "HomePatch.h"
#include "InfoStream.h"
#include "CollectionMaster.h"
#include "CollectionMgr.h"

#include "ParallelIOMgr.decl.h"
#include "ParallelIOMgr.h"

#include "Output.h"
#include "Random.h"

#include <algorithm>
using namespace std;

ParallelIOMgr::ParallelIOMgr()
{
    CkpvAccess(BOCclass_group).ioMgr = thisgroup;

    numInputProcs=-1;
    inputProcArray = NULL;
    numOutputProcs=-1;
    outputProcArray = NULL;

    procsReceived=0;
    hydroMsgRecved=0;

    totalMV.x = totalMV.y = totalMV.z = 0.0;
    totalMass = 0.0;
    totalCharge = 0.0;

    isOKToRecvHPAtoms = false;
    hpAtomsList = NULL;

    clusterID = NULL;
    clusterSize = NULL;

#ifdef MEM_OPT_VERSION
    midCM = NULL;
#endif

    isWater = NULL;

    numCSMAck = 0;
    numReqRecved = 0;

    sendAtomsThread = 0;

#if COLLECT_PERFORMANCE_DATA
    numFixedAtomLookup = 0;
#endif
}

ParallelIOMgr::~ParallelIOMgr()
{
    delete [] inputProcArray;
    delete [] outputProcArray;
    delete [] clusterID;
    delete [] clusterSize;

#ifdef MEM_OPT_VERSION
    delete midCM;
#endif

    delete [] isWater;
}

#ifndef OUTPUT_SINGLE_FILE
#error OUTPUT_SINGLE_FILE not defined!
#endif

// initialization needed for the parallel IO manager. Entry point through Scripttcl
void ParallelIOMgr::initialize(Node *node)
{
    simParameters = node->simParameters;
    molecule = node->molecule;

    numInputProcs = simParameters->numinputprocs;
    numOutputProcs = simParameters->numoutputprocs;
    numOutputWrts = simParameters->numoutputwrts;


    if(!CkMyPe()) {
        iout << iINFO << "Running with " <<numInputProcs<<" input processors.\n"<<endi;
        #if OUTPUT_SINGLE_FILE
        iout << iINFO << "Running with " <<numOutputProcs<<" output processors ("<<numOutputWrts<<" of them will output simultaneously).\n"<<endi;
        #else
        iout << iINFO << "Running with " <<numOutputProcs<<" output processors, and each of them will output to its own separate file.\n"<<endi;
        #endif
    }

    //build inputProcArray
   {
    inputProcArray = new int[numInputProcs];
    myInputRank = -1;
    for(int i=0; i<numInputProcs; ++i) {
      inputProcArray[i] = WorkDistrib::peDiffuseOrdering[(1+numOutputProcs+i)%CkNumPes()];
    }
    std::sort(inputProcArray, inputProcArray+numInputProcs);
    for(int i=0; i<numInputProcs; ++i) {
      if ( CkMyPe() == inputProcArray[i] ) {
        if ( myInputRank != -1 ) NAMD_bug("Duplicate input proc");
        myInputRank = i;
      }
    }

    if(!CkMyPe()) {
      iout << iINFO << "INPUT PROC LOCATIONS:";
      int i;
      for ( i=0; i<numInputProcs && i < 10; ++i ) {
        iout << " " << inputProcArray[i];
      }
      if ( i<numInputProcs ) iout << " ... " << inputProcArray[numInputProcs-1];
      iout << "\n" << endi;
    }
   }

    if(myInputRank!=-1) {
        //NOTE: this could further be optimized by pre-allocate the memory
        //for incoming atoms --Chao Mei
        int numMyAtoms = numInitMyAtomsOnInput();
        initAtoms.resize(numMyAtoms+100);  // extra space for orphan hydrogens
        initAtoms.resize(numMyAtoms);
        tmpRecvAtoms.resize(0);
    } else {
        initAtoms.resize(0);
        tmpRecvAtoms.resize(0);
    }
    hpIDList.resize(0);

    //build outputProcArray
    //spread the output processors across all the processors
   {
    outputProcArray = new int[numOutputProcs];
    outputProcFlags = new char[CkNumPes()];
    myOutputRank = -1;
    for(int i=0; i<numOutputProcs; ++i) {
      outputProcArray[i] = WorkDistrib::peDiffuseOrdering[(1+i)%CkNumPes()];
    }
    std::sort(outputProcArray, outputProcArray+numOutputProcs);
    for(int i=0; i<CkNumPes(); ++i) {
      outputProcFlags[i] = 0;
    }
    for(int i=0; i<numOutputProcs; ++i) {
      outputProcFlags[outputProcArray[i]] = 1;
      if ( CkMyPe() == outputProcArray[i] ) {
        if ( myOutputRank != -1 ) NAMD_bug("Duplicate output proc");
        myOutputRank = i;
      }
    }

    if(!CkMyPe()) {
      iout << iINFO << "OUTPUT PROC LOCATIONS:";
      int i;
      for ( i=0; i<numOutputProcs && i < 10; ++i ) {
        iout << " " << outputProcArray[i];
      }
      if ( i<numOutputProcs ) iout << " ... " << outputProcArray[numOutputProcs-1];
      iout << "\n" << endi;
    }
   }

#ifdef MEM_OPT_VERSION
    if(myOutputRank!=-1) {
        midCM = new CollectionMidMaster(this);
    }
    remoteClusters.clear();
    csmBuf.resize(0);
    remoteCoors.clear();
    ccmBuf.resize(0);

    mainMaster = CollectionMgr::Object()->getMasterChareID();
#endif
}

bool ParallelIOMgr::isOutputProcessor(int pe) {
   return outputProcFlags[pe];
}

int isOutputProcessor(int pe){ 
  return CProxy_ParallelIOMgr::ckLocalBranch(CkpvAccess(BOCclass_group).ioMgr)->isOutputProcessor(pe);
}



void ParallelIOMgr::readPerAtomInfo()
{
#ifdef MEM_OPT_VERSION
    if(myInputRank!=-1) {
        int myAtomLIdx, myAtomUIdx;
        getMyAtomsInitRangeOnInput(myAtomLIdx, myAtomUIdx);

        //1. read the file that contains per-atom info such as signature index
        molecule->read_binary_atom_info(myAtomLIdx, myAtomUIdx, initAtoms);

        //2. read coordinates and velocities of each atom if the velocity file
        //exists, otherwise, the velocity of each atom is randomly generated.
        //This has to be DONE AFTER THE FIRST STEP as the atom mass is required
        //if the velocity is generated randomly.
        readCoordinatesAndVelocity();

        //3. set every atom's output processor rank, i.e. the dest pe this
        //atom will be sent for writing positions and velocities etc.
        int oRank=atomRankOnOutput(myAtomLIdx);
        for(int i=oRank; i<numOutputProcs; i++) {
            int lIdx, uIdx; //indicates the range of atom ids outputProcArray[i] has
            getAtomsRangeOnOutput(lIdx, uIdx, i);
            if(lIdx > myAtomUIdx) break;
            int fid = lIdx>myAtomLIdx?lIdx:myAtomLIdx;
            int tid = uIdx>myAtomUIdx?myAtomUIdx:uIdx;
            for(int j=fid; j<=tid; j++) initAtoms[j-myAtomLIdx].outputRank = i;
        }
    }

    //read clusters
    if(myOutputRank!=-1) {
        //only when wrapAll or wrapWater is set, cluster info is required
        if(!(simParameters->wrapAll || simParameters->wrapWater)) return;
        readInfoForParOutput();
    }
#endif
}

void ParallelIOMgr::readCoordinatesAndVelocity()
{
#ifdef MEM_OPT_VERSION
    int needFlip = 0;
    int myAtomLIdx, myAtomUIdx;
    getMyAtomsInitRangeOnInput(myAtomLIdx, myAtomUIdx);
    int myNumAtoms = myAtomUIdx-myAtomLIdx+1;

    //contains the data for Position and Velocity
    Vector *tmpData = new Vector[myNumAtoms];

    //begin to read coordinates
    //step1: open the file
    FILE *ifp = fopen(simParameters->binCoorFile, "rb");
    if(!ifp) {
        char s[256];
        sprintf(s, "The binary coordinate file %s cannot be opened on proc %d\n", simParameters->binCoorFile, CkMyPe());
        NAMD_err(s);
    }
    //step2: check whether flip is needed
    int filelen;
    fread(&filelen, sizeof(int32),1,ifp);
    char lenbuf[sizeof(int32)];
    memcpy(lenbuf, (const char *)&filelen, sizeof(int32));
    flipNum(lenbuf, sizeof(int32), 1);
    if(!memcmp(lenbuf, (const char *)&filelen, sizeof(int32))) {
        iout << iWARN << "Number of atoms in binary file " << simParameters->binCoorFile
             <<" is palindromic, assuming same endian.\n" << endi;
    }
    if(filelen!=molecule->numAtoms) {
        needFlip = 1;
        memcpy((void *)&filelen, lenbuf,sizeof(int32));
    }
    if(filelen!=molecule->numAtoms) {
        char s[256];
        sprintf(s, "Incorrect atom count in binary file %s", simParameters->binCoorFile);
        NAMD_die(s);
    }
    //step3: read the file specified by the range
    int64 offsetPos = ((int64)myAtomLIdx)*sizeof(Position);
#ifdef WIN32
    if ( _fseeki64(ifp, offsetPos, SEEK_CUR) )
#else
    if ( fseeko(ifp, offsetPos, SEEK_CUR) )
#endif
    {
        char s[256];
        sprintf(s, "Error in seeking binary file %s on proc %d",  simParameters->binCoorFile, CkMyPe());
        NAMD_err(s);
    }
    size_t totalRead = fread(tmpData, sizeof(Vector), myNumAtoms, ifp);
    if(totalRead!=myNumAtoms) {
        char s[256];
        sprintf(s, "Error in reading binary file %s on proc %d",  simParameters->binCoorFile, CkMyPe());
        NAMD_err(s);
    }
    if(needFlip) flipNum((char *)tmpData, sizeof(BigReal), myNumAtoms*3);
    fclose(ifp);
    for(int i=0; i<myNumAtoms; i++) initAtoms[i].position = tmpData[i];

    //begin to read velocity
    //step1: generate velocity randomly or open the file
    if(!simParameters->binVelFile) {
        //generate velocity randomly
        Node::Object()->workDistrib->random_velocities_parallel(simParameters->initialTemp, initAtoms);
    } else {
        ifp = fopen(simParameters->binVelFile, "rb");
        if(!ifp) {
            char s[256];
            sprintf(s, "The binary velocity file %s cannot be opened on proc %d\n", simParameters->binVelFile, CkMyPe());
            NAMD_err(s);
        }
        //step2: check whether flip is needed
        fread(&filelen, sizeof(int32),1,ifp);
        memcpy(lenbuf, (const char *)&filelen, sizeof(int32));
        flipNum(lenbuf, sizeof(int32), 1);
        if(!memcmp(lenbuf, (const char *)&filelen, sizeof(int32))) {
            iout << iWARN << "Number of atoms in binary file " << simParameters->binVelFile
                 <<" is palindromic, assuming same endian.\n" << endi;
        }
        if(filelen!=molecule->numAtoms) {
            needFlip = 1;
            memcpy((void *)&filelen, lenbuf,sizeof(int32));
        }
        if(filelen!=molecule->numAtoms) {
            char s[256];
            sprintf(s, "Incorrect atom count in binary file %s", simParameters->binVelFile);
            NAMD_die(s);
        }

        //step3: read the file specified by the range
        int64 offsetPos = ((int64)myAtomLIdx)*sizeof(Velocity);
#ifdef WIN32
        if ( _fseeki64(ifp, offsetPos, SEEK_CUR) )
#else
        if ( fseeko(ifp, offsetPos, SEEK_CUR) )
#endif
        {
            char s[256];
            sprintf(s, "Error in seeking binary file %s on proc %d",  simParameters->binVelFile, CkMyPe());
            NAMD_err(s);
        }
        totalRead = fread(tmpData, sizeof(Vector), myNumAtoms, ifp);
        if(totalRead!=myNumAtoms) {
            char s[256];
            sprintf(s, "Error in reading binary file %s on proc %d",  simParameters->binVelFile, CkMyPe());
            NAMD_err(s);
        }
        if(needFlip) flipNum((char *)tmpData, sizeof(BigReal), myNumAtoms*3);
        fclose(ifp);
        for(int i=0; i<myNumAtoms; i++) initAtoms[i].velocity = tmpData[i];
    }

    //begin to read reference coordinates
    //step1: use initial positions or open the file
    if(!simParameters->binRefFile) {
        for(int i=0; i<myNumAtoms; i++) initAtoms[i].fixedPosition = initAtoms[i].position;
    } else {
        ifp = fopen(simParameters->binRefFile, "rb");
        if(!ifp) {
            char s[256];
            sprintf(s, "The binary reference coordinate file %s cannot be opened on proc %d\n", simParameters->binRefFile, CkMyPe());
            NAMD_err(s);
        }
        //step2: check whether flip is needed
        fread(&filelen, sizeof(int32),1,ifp);
        memcpy(lenbuf, (const char *)&filelen, sizeof(int32));
        flipNum(lenbuf, sizeof(int32), 1);
        if(!memcmp(lenbuf, (const char *)&filelen, sizeof(int32))) {
            iout << iWARN << "Number of atoms in binary file " << simParameters->binRefFile
                 <<" is palindromic, assuming same endian.\n" << endi;
        }
        if(filelen!=molecule->numAtoms) {
            needFlip = 1;
            memcpy((void *)&filelen, lenbuf,sizeof(int32));
        }
        if(filelen!=molecule->numAtoms) {
            char s[256];
            sprintf(s, "Incorrect atom count in binary file %s", simParameters->binRefFile);
            NAMD_die(s);
        }

        //step3: read the file specified by the range
        int64 offsetPos = ((int64)myAtomLIdx)*sizeof(Position);
#ifdef WIN32
        if ( _fseeki64(ifp, offsetPos, SEEK_CUR) )
#else
        if ( fseeko(ifp, offsetPos, SEEK_CUR) )
#endif
        {
            char s[256];
            sprintf(s, "Error in seeking binary file %s on proc %d",  simParameters->binRefFile, CkMyPe());
            NAMD_err(s);
        }
        totalRead = fread(tmpData, sizeof(Vector), myNumAtoms, ifp);
        if(totalRead!=myNumAtoms) {
            char s[256];
            sprintf(s, "Error in reading binary file %s on proc %d",  simParameters->binRefFile, CkMyPe());
            NAMD_err(s);
        }
        if(needFlip) flipNum((char *)tmpData, sizeof(BigReal), myNumAtoms*3);
        fclose(ifp);
        for(int i=0; i<myNumAtoms; i++) initAtoms[i].fixedPosition = tmpData[i];
    }

    delete [] tmpData;
#endif
}

void ParallelIOMgr::readInfoForParOutput()
{
    int fromIdx, toIdx; //atoms' range
    getMyAtomsRangeOnOutput(fromIdx,toIdx);
    int numMyAtoms = toIdx-fromIdx+1;

    clusterID = new int[numMyAtoms];
    clusterSize = new int[numMyAtoms];

    //Since the input proc also reads this file, all the checks
    //(file version check, atom record size check) are omitted here
    FILE *ifp = fopen(simParameters->binAtomFile, "rb");
    //read magic number to set needFlip
    int needFlip = 0;
    int magicNum;
    fread(&magicNum, sizeof(int), 1, ifp);
    if (magicNum!=COMPRESSED_PSF_MAGICNUM) {
        needFlip = 1;
    }

    //need to load isWater info
    isWater = new char[numMyAtoms];
    //seek from the end of the file (note offset is negative!)
    int64 offset = sizeof(char)*((int64)(fromIdx-molecule->numAtoms));
#ifdef WIN32
    if ( _fseeki64(ifp, offset, SEEK_END) )
#else
    if ( fseeko(ifp, offset, SEEK_END) )
#endif
    {
        char s[256];
        sprintf(s, "Error in seeking binary file %s on proc %d",  simParameters->binAtomFile, CkMyPe());
        NAMD_err(s);
    }
    fread(isWater, sizeof(char), numMyAtoms, ifp);
    //there's no need for flipping as it's a char array

    //seek from the end of the file (note offset is negative!)
    offset = sizeof(int)*((int64)(fromIdx-molecule->numAtoms))
                   - sizeof(char)*((int64)(molecule->numAtoms));
#ifdef WIN32
    if ( _fseeki64(ifp, offset, SEEK_END) )
#else
    if ( fseeko(ifp, offset, SEEK_END) )
#endif
    {
        char s[256];
        sprintf(s, "Error in seeking binary file %s on proc %d",  simParameters->binAtomFile, CkMyPe());
        NAMD_err(s);
    }
    fread(clusterID, sizeof(int), numMyAtoms, ifp);
    if(needFlip) flipNum((char *)clusterID, sizeof(int), numMyAtoms);
    fclose(ifp);

    //calculate cluster size (communication may be neccessary)
    ClusterElem one;
    for(int i=0; i<numMyAtoms; i++) {
        clusterSize[i] = 0;        
        int cid = clusterID[i];
        //check if the cluster id (i.e. the header of the atom of
        //the cluster) is in my local repository. If the cluster id
        //is not in my local repository, then it must appear in output
        //processors that are in front of this one because of the way
        //of determing cluster info for the molecular system.
        CmiAssert(cid<=toIdx);
        if(cid<fromIdx) {
            //on output procs ahead of me
            one.clusterId = cid;
            ClusterElem *ret = remoteClusters.find(one);
            if(ret==NULL) {
                one.atomsCnt = 1;
                remoteClusters.add(one);
            } else {
                ret->atomsCnt++;
            }
        } else {
            int lidx = cid-fromIdx;
            CmiAssert(lidx<=i);
            clusterSize[lidx]++;
        }
    }

    //Prepare to send msgs to remote output procs to reduce the cluster size
    //Since the expected number of msgs to be very small, msgs to the same proc
    //are not aggregated. --Chao Mei
#if 0
    printf("output[%d]=%d: prepare to send %d remote msgs for cluster size\n",
           myOutputRank, CkMyPe(), remoteClusters.size());
#endif

    numRemoteClusters = remoteClusters.size();
    numCSMAck = 0; //set to 0 to prepare recving the final cluster size update
    CProxy_ParallelIOMgr pIO(thisgroup);
    ClusterSetIter iter(remoteClusters);
    for(iter=iter.begin(); iter!=iter.end(); iter++) {
        ClusterSizeMsg *msg = new ClusterSizeMsg;
        msg->srcRank = myOutputRank;
        msg->clusterId = iter->clusterId;
        msg->atomsCnt = iter->atomsCnt;
        int dstRank = atomRankOnOutput(iter->clusterId);
        pIO[outputProcArray[dstRank]].recvClusterSize(msg);
    }
}

void ParallelIOMgr::recvClusterSize(ClusterSizeMsg *msg)
{
    csmBuf.add(msg); //added to buffer for reuse to send back to src

    //update cluster size has to be delayed to integration to prevent
    //data racing where the clusterSize has not been created!
}

void ParallelIOMgr::integrateClusterSize()
{
    if(myOutputRank==-1) return;
    if(!(simParameters->wrapAll || simParameters->wrapWater)) return;

    int fromIdx, toIdx; //atoms' range
    getMyAtomsRangeOnOutput(fromIdx,toIdx);

    //calculated the final cluster size
    for(int i=0; i<csmBuf.size(); i++) {
        ClusterSizeMsg *msg = csmBuf[i];
        int lidx = msg->clusterId - fromIdx;
        clusterSize[lidx] += msg->atomsCnt;
    }

    CProxy_ParallelIOMgr pIO(thisgroup);
    for(int i=0; i<csmBuf.size(); i++) {
        ClusterSizeMsg *msg = csmBuf[i];
        int lidx = msg->clusterId - fromIdx;
        msg->atomsCnt = clusterSize[lidx];
        pIO[outputProcArray[msg->srcRank]].recvFinalClusterSize(msg);
    }
    numRemoteReqs = csmBuf.size();
    csmBuf.resize(0);
    
    //There's a possible msg race problem here that recvFinalClusterSize 
    //executes before integrateClusterSize because other proc finishes faster
    //in calculating the cluster size. The recvFinalClusterSize should be
    //executed after integrateClusterSize. To avoid this, a self message is
    //sent to participate the reduction.
    if(numRemoteClusters!=0){
        recvFinalClusterSize(NULL);
    }else{
        //this output proc already has the final cluster size for each atom
        int numMyAtoms = toIdx-fromIdx+1;
        for(int i=0; i<numMyAtoms; i++) {
            int lidx = clusterID[i]-fromIdx;
            clusterSize[i] = clusterSize[lidx];
        }
        
        #if 0 //write out cluster debug info
        char fname[128];
        sprintf(fname, "cluster.par.%d", CkMyPe());
        FILE *ofp = fopen(fname, "w");
        for(int i=0; i<numMyAtoms; i++) {
            fprintf(ofp, "%d: %d: %d\n", i+fromIdx, clusterID[i], clusterSize[i]);
        }
        fclose(ofp);
        #endif
    }
}

void ParallelIOMgr::recvFinalClusterSize(ClusterSizeMsg *msg)
{
    //only process the message sent by other procs
    if(msg!=NULL) {
        //indicating a message from other procs
        ClusterElem one(msg->clusterId);
        ClusterElem *ret = remoteClusters.find(one);
        CmiAssert(ret!=NULL);
        ret->atomsCnt = msg->atomsCnt;
    }
    delete msg;

    //include a msg sent by itself for reduction
    if(++numCSMAck == (numRemoteClusters+1)) {
        //recved all the msgs needed to update the cluster size for each atom finally
        int fromIdx, toIdx; //atoms' range
        getMyAtomsRangeOnOutput(fromIdx,toIdx);
        int numMyAtoms = toIdx-fromIdx+1;
        ClusterElem tmp;
        for(int i=0; i<numMyAtoms; i++) {
            int cid = clusterID[i];
            int lidx = cid-fromIdx;
            if(lidx<0) {
                //this cid should be inside remoteClusters
                tmp.clusterId = cid;
                ClusterElem *fone = remoteClusters.find(tmp);
                clusterSize[i] = fone->atomsCnt;
            } else {
                clusterSize[i] = clusterSize[lidx];
            }
        }
        numCSMAck = 0;
        remoteClusters.clear();

#if 0 //write out cluster debug info
        char fname[128];
        sprintf(fname, "cluster.par.%d", CkMyPe());
        FILE *ofp = fopen(fname, "w");
        for(int i=0; i<numMyAtoms; i++) {
            fprintf(ofp, "%d: %d: %d\n", i+fromIdx, clusterID[i], clusterSize[i]);
        }
        fclose(ofp);
#endif

    }
}

void ParallelIOMgr::migrateAtomsMGrp()
{
    if(myInputRank==-1) return;

    //1. first get the list of atoms to be migrated
    //which should be few compared with the number of atoms
    //initially assigned to this input proc.
    AtomIDList toMigrateList; //the list of atoms to be migrated
    //the max distance from this processor of atoms to be sent
    int maxOffset = 0;
    for(int i=0; i<initAtoms.size(); i++) {
        //returns the proc id on which atom MPID resides on
        int parentRank = atomInitRankOnInput(initAtoms[i].MPID);
        if(parentRank != myInputRank) {
            toMigrateList.add(i);
            initAtoms[i].isValid = false;
            int tmp = parentRank - myInputRank;
            tmp = tmp>0 ? tmp : -tmp;
            if(tmp > maxOffset) maxOffset = tmp;
        }
    }

    //2. prepare atom migration messages
    //the messages are indexed as [-maxOffset,..., -1,0,1,..., maxOffset]
    //where the "0" is not used at all. It is added for the sake of
    //computing the index easily.
    InputAtomList *migLists = new InputAtomList[2*maxOffset+1];
    for(int i=0; i<toMigrateList.size(); i++) {
        int idx = toMigrateList[i];
        int parentRank = atomInitRankOnInput(initAtoms[idx].MPID);
        //decide which migList to put this atom
        int offset = parentRank - myInputRank + maxOffset;
        migLists[offset].add(initAtoms[idx]);
    }

    CProxy_ParallelIOMgr pIO(thisgroup);
    for(int i=0; i<2*maxOffset+1; i++) {
        int migLen = migLists[i].size();
        if(migLen>0) {
            MoveInputAtomsMsg *msg = new (migLen, 0)MoveInputAtomsMsg;
            msg->length = migLen;
            memcpy(msg->atomList, migLists[i].begin(), sizeof(InputAtom)*migLen);
            int destRank = i-maxOffset+myInputRank;
            pIO[inputProcArray[destRank]].recvAtomsMGrp(msg);
            migLists[i].clear();
        }
    }
    
    toMigrateList.clear();
    delete [] migLists;
}

void ParallelIOMgr::recvAtomsMGrp(MoveInputAtomsMsg *msg)
{
    for(int i=0; i<msg->length; i++) {
        tmpRecvAtoms.add((msg->atomList)[i]);
    }
    delete msg;
}

void ParallelIOMgr::integrateMigratedAtoms()
{
    if(myInputRank==-1) return;

    for(int i=0; i<tmpRecvAtoms.size(); i++) {
        tmpRecvAtoms[i].isValid = true;
        initAtoms.add(tmpRecvAtoms[i]);
    }
    tmpRecvAtoms.clear();

    //sort atom list based on hydrogenList value
    std::sort(initAtoms.begin(), initAtoms.end());

    //now compute the counters inside Molecule such as numFixedRigidBonds
    //which is based on the hydrogen group info

    int numFixedRigidBonds = 0;
    if(molecule->numRigidBonds){
        int parentIsFixed = 0;
        for(int i=0; i<initAtoms.size(); i++) {
            InputAtom *one = &(initAtoms[i]);
            if(!one->isValid) continue;
            if(one->isGP) {
                parentIsFixed = one->atomFixed;
                InputAtom *a1 = &(initAtoms[i+1]);
                InputAtom *a2 = &(initAtoms[i+2]);
                if((one->rigidBondLength>0.0) &&
                   a1->atomFixed && a2->atomFixed) {
                    numFixedRigidBonds++;
                }
            }else{
                if((one->rigidBondLength>0.0) &&
                   one->atomFixed && parentIsFixed) {
                    numFixedRigidBonds++;
                }
            }
        }
    }

    int numFixedGroups = 0;
    if(molecule->numFixedAtoms){        
        for(int i=0; i<initAtoms.size();) {
            InputAtom *one = &(initAtoms[i]);
            if(!one->isValid){
                i++;
                continue;
            }
            if(one->isGP) {
                int allFixed = 1;                
                for(int j=0; j<one->hydrogenGroupSize; j++){
                    InputAtom *a1 = &(initAtoms[i+j]);
                    allFixed = allFixed & a1->atomFixed;
                    if(!allFixed) break;
                }
                if(allFixed) numFixedGroups++;                
                i += one->hydrogenGroupSize;
            }
        }
    }
    
    CProxy_ParallelIOMgr pIO(thisgroup);
    HydroBasedMsg *msg = new HydroBasedMsg;
    msg->numFixedGroups = numFixedGroups;
    msg->numFixedRigidBonds = numFixedRigidBonds;
    pIO[0].recvHydroBasedCounter(msg);
}

void ParallelIOMgr::updateMolInfo()
{
#ifdef MEM_OPT_VERSION
    if(myInputRank==-1) return;

    CProxy_ParallelIOMgr pIO(thisgroup);

    MolInfoMsg *msg = new MolInfoMsg;
    msg->numBonds = msg->numCalcBonds = 0;
    msg->numAngles = msg->numCalcAngles = 0;
    msg->numDihedrals = msg->numCalcDihedrals = 0;
    msg->numImpropers = msg->numCalcImpropers = 0;
    msg->numCrossterms = msg->numCalcCrossterms = 0;
    msg->numExclusions = msg->numCalcExclusions = 0;
    int numFullExclusions = msg->numCalcFullExclusions = 0;
    // JLai
    msg->numLJPairs = msg->numCalcLJPairs = 0;
    // End of JLai
    msg->numRigidBonds = 0;
    msg->totalMass = 0.0;
    msg->totalCharge = 0.0;

    //calculate the tuples this input processor have
    AtomSignature *atomSigPool = molecule->atomSigPool;
    ExclusionSignature *exclSigPool = molecule->exclSigPool;
    for(int i=0; i<initAtoms.size(); i++) {
        AtomSignature *thisSig = &atomSigPool[initAtoms[i].sigId];
        msg->numBonds += thisSig->bondCnt;
        msg->numAngles += thisSig->angleCnt;
        msg->numDihedrals += thisSig->dihedralCnt;
        msg->numImpropers += thisSig->improperCnt;
        msg->numCrossterms += thisSig->crosstermCnt;
	// JLai
	msg->numLJPairs += thisSig->gromacsPairCnt;
	// End of JLai

        ExclusionSignature *exclSig = &exclSigPool[initAtoms[i].exclId];
        msg->numExclusions += (exclSig->fullExclCnt + exclSig->modExclCnt);
        numFullExclusions += exclSig->fullExclCnt;

        if(initAtoms[i].rigidBondLength > 0.0) msg->numRigidBonds++;

        msg->totalMass += initAtoms[i].mass;
        msg->totalCharge += initAtoms[i].charge;
    }

    //deal with numCalc* which is related with fixed atoms!
    if(molecule->numFixedAtoms>0 && ! simParameters->fixedAtomsForces) {
        //if there's fixed atoms, calcExclusions needs to be calculated
        //Since it's possible the atom inside the this exclusion set is on
        //another input processor, we have to resort to the global fixed atoms
        //info inside the Molecule object. The number of such accesses should
        //be very small! --Chao Mei
        int sAId = initAtoms[0].id;
        int remoteCnt=0; //stats info
        for(int i=0; i<initAtoms.size(); i++) {
            //When all the atoms in the set are fixed, the elem (Bond etc.)
            //is not counted as a calc*.
            int myAId = initAtoms[i].id;
            AtomSignature *thisSig = &atomSigPool[initAtoms[i].sigId];
            ExclusionSignature *exclSig = &exclSigPool[initAtoms[i].exclId];
            if(!initAtoms[i].atomFixed) {
                msg->numCalcBonds += thisSig->bondCnt;                
                msg->numCalcAngles += thisSig->angleCnt;
                msg->numCalcDihedrals += thisSig->dihedralCnt;
                msg->numCalcImpropers += thisSig->improperCnt;
                msg->numCalcCrossterms += thisSig->crosstermCnt;
                msg->numCalcExclusions+=(exclSig->fullExclCnt+exclSig->modExclCnt);
                msg->numCalcFullExclusions+=(exclSig->fullExclCnt);
                continue;
            }
                       
            //1. Bonds
            for(int j=0; j<thisSig->bondCnt; j++) {            
                TupleSignature *bsig = &(thisSig->bondSigs[j]);
                int a1 = myAId + bsig->offset[0];
                if(!isAtomFixed(sAId, a1)) msg->numCalcBonds++;
            }
            
            //2. Angles
            for(int j=0; j<thisSig->angleCnt; j++) {            
                TupleSignature *bsig = &(thisSig->angleSigs[j]);
                int a1 = myAId + bsig->offset[0];
                int a2 = myAId + bsig->offset[1];
                if(!isAtomFixed(sAId, a1) || !isAtomFixed(sAId, a2)) 
                    msg->numCalcAngles++;
            }

            //3. Dihedrals
            for(int j=0; j<thisSig->dihedralCnt; j++) {            
                TupleSignature *bsig = &(thisSig->dihedralSigs[j]);
                int a1 = myAId + bsig->offset[0];
                int a2 = myAId + bsig->offset[1];
                int a3 = myAId + bsig->offset[2];
                if(!isAtomFixed(sAId, a1) || 
                   !isAtomFixed(sAId, a2) ||
                   !isAtomFixed(sAId, a3)) 
                    msg->numCalcDihedrals++;
            }

            //4. Impropers
            for(int j=0; j<thisSig->improperCnt; j++) {            
                TupleSignature *bsig = &(thisSig->improperSigs[j]);
                int a1 = myAId + bsig->offset[0];
                int a2 = myAId + bsig->offset[1];
                int a3 = myAId + bsig->offset[2];
                if(!isAtomFixed(sAId, a1) || 
                   !isAtomFixed(sAId, a2) ||
                   !isAtomFixed(sAId, a3)) 
                    msg->numCalcImpropers++;
            }

            //5. Crossterms
            for(int j=0; j<thisSig->crosstermCnt; j++) {            
                TupleSignature *bsig = &(thisSig->crosstermSigs[j]);
                int a1 = myAId + bsig->offset[0];
                int a2 = myAId + bsig->offset[1];
                int a3 = myAId + bsig->offset[2];
                int a4 = myAId + bsig->offset[3];
                int a5 = myAId + bsig->offset[4];
                int a6 = myAId + bsig->offset[5];
                int a7 = myAId + bsig->offset[6];

                if(!isAtomFixed(sAId, a1) || 
                   !isAtomFixed(sAId, a2) ||
                   !isAtomFixed(sAId, a3) ||
                   !isAtomFixed(sAId, a4) ||
                   !isAtomFixed(sAId, a5) ||
                   !isAtomFixed(sAId, a6) ||
                   !isAtomFixed(sAId, a7)) 
                    msg->numCalcDihedrals++;
            }
            
            //6: Exclusions            
            //this atom is fixed, check atoms in the exclusion set
            for(int j=0; j<exclSig->fullExclCnt; j++) {
                int thisAId = exclSig->fullOffset[j]+myAId;
                if(!isAtomFixed(sAId, thisAId)) { msg->numCalcExclusions++; msg->numCalcFullExclusions++; }
            }
            for(int j=0; j<exclSig->modExclCnt; j++) {
                int thisAId = exclSig->modOffset[j]+myAId;
                if(!isAtomFixed(sAId, thisAId)) msg->numCalcExclusions++;
            }

	    //7: GromacsPair
	    for(int j=0; j<thisSig->gromacsPairCnt; j++) {
		TupleSignature *bsig = &(thisSig->gromacsPairSigs[j]);
		int a1 = myAId + bsig->offset[0];
		int a2 = myAId + bsig->offset[1];
                if(!isAtomFixed(sAId, a1) || 
                   !isAtomFixed(sAId, a2))
                    msg->numCalcLJPairs++;
	    }
        }
#if COLLECT_PERFORMANCE_DATA
        printf("Num fixedAtom lookup on proc %d is %d\n", CkMyPe(), numFixedAtomLookup);
#endif
    } else {
        //no fixed atoms, numCalc* is same with numExclusions
        msg->numCalcBonds = msg->numBonds;
        msg->numCalcAngles = msg->numAngles;
        msg->numCalcDihedrals = msg->numDihedrals;
        msg->numCalcImpropers = msg->numImpropers;
        msg->numCalcCrossterms = msg->numCrossterms;
        msg->numCalcExclusions = msg->numExclusions;
        msg->numCalcFullExclusions = numFullExclusions;
    }


    if(!simParameters->comMove) {
        //to remove the center of mass motion from a molecule.
        //first calculate the values on every input proc, then reduce.
        //For more info, refer to WorkDistrib::remove_com_motion
        //-Chao Mei
        (msg->totalMV).x = 0.0;
        (msg->totalMV).y = 0.0;
        (msg->totalMV).z = 0.0;
        for (int i=0; i<initAtoms.size(); i++) {            
            msg->totalMV += initAtoms[i].mass * initAtoms[i].velocity;
        }
    }

    //always send to the master processor (proc 0)
    pIO[0].recvMolInfo(msg);
#endif
}

//only executed on proc 0
void ParallelIOMgr::recvMolInfo(MolInfoMsg *msg)
{
    molecule->numBonds += msg->numBonds;
    molecule->numCalcBonds += msg->numCalcBonds;
    molecule->numAngles += msg->numAngles;
    molecule->numCalcAngles += msg->numCalcAngles;
    molecule->numDihedrals += msg->numDihedrals;
    molecule->numCalcDihedrals += msg->numCalcDihedrals;
    molecule->numImpropers += msg->numImpropers;
    molecule->numCalcImpropers += msg->numCalcImpropers;
    molecule->numCrossterms += msg->numCrossterms;
    molecule->numCalcCrossterms += msg->numCalcCrossterms;
    molecule->numTotalExclusions += msg->numExclusions;
    molecule->numCalcExclusions += msg->numCalcExclusions;
    molecule->numCalcFullExclusions += msg->numCalcFullExclusions;
    molecule->numRigidBonds += msg->numRigidBonds;

    totalMass += msg->totalMass;
    totalCharge += msg->totalCharge;

    if(!simParameters->comMove) {
        totalMV += msg->totalMV;        
    }

    if(++procsReceived == numInputProcs) {
        //received all the counters
        msg->numBonds = molecule->numBonds;
        msg->numCalcBonds = molecule->numCalcBonds;
        msg->numAngles = molecule->numAngles;
        msg->numCalcAngles = molecule->numCalcAngles;
        msg->numDihedrals = molecule->numDihedrals;
        msg->numCalcDihedrals = molecule->numCalcDihedrals;
        msg->numImpropers = molecule->numImpropers;
        msg->numCalcImpropers = molecule->numCalcImpropers;
        msg->numCrossterms = molecule->numCrossterms;
        msg->numCalcCrossterms = molecule->numCalcCrossterms;
        msg->numExclusions = molecule->numTotalExclusions/2;
        msg->numCalcExclusions = molecule->numCalcExclusions/2;
        msg->numCalcFullExclusions = molecule->numCalcFullExclusions/2;
        msg->numRigidBonds = molecule->numRigidBonds;

        msg->totalMass = totalMass;
        msg->totalCharge = totalCharge;

        if(!simParameters->comMove) {
            msg->totalMV = totalMV;            
        }

        CProxy_ParallelIOMgr pIO(thisgroup);
        pIO.bcastMolInfo(msg);

        //reset to 0 for the next p2p-based reduction on input procs
        procsReceived = 0;
    } else delete msg;
}

void ParallelIOMgr::bcastMolInfo(MolInfoMsg *msg)
{
#ifdef MEM_OPT_VERSION
    if(myInputRank!=-1) {
        if(!simParameters->comMove) {
            //needs to remove the center of mass motion from a molecule
            Vector val = msg->totalMV / msg->totalMass;
            for (int i=0; i<initAtoms.size(); i++) initAtoms[i].velocity -= val;
        }
    }

    //only the rank 0 in the SMP node update the Molecule object
    if(CmiMyRank()) {
        delete msg;
        return;
    }

    molecule->numBonds = msg->numBonds;
    molecule->numCalcBonds = msg->numCalcBonds;
    molecule->numAngles = msg->numAngles;
    molecule->numCalcAngles = msg->numCalcAngles;
    molecule->numDihedrals = msg->numDihedrals;
    molecule->numCalcDihedrals = msg->numCalcDihedrals;
    molecule->numImpropers = msg->numImpropers;
    molecule->numCalcImpropers = msg->numCalcImpropers;
    molecule->numCrossterms = msg->numCrossterms;
    molecule->numCalcCrossterms = msg->numCalcCrossterms;

    molecule->numTotalExclusions = msg->numExclusions;
    molecule->numCalcExclusions = msg->numCalcExclusions;
    molecule->numCalcFullExclusions = msg->numCalcFullExclusions;

    molecule->numRigidBonds = msg->numRigidBonds;
    delete msg;

    if(!CkMyPe()) {
        iout << iINFO << "LOADED " << molecule->numTotalExclusions << " TOTAL EXCLUSIONS\n" << endi;
        if(!simParameters->comMove) {
            iout << iINFO << "REMOVING COM VELOCITY "
                 << (PDBVELFACTOR * (msg->totalMV / msg->totalMass))<< "\n" <<endi;
        }
    }
#endif
}

//only called on PE0
void ParallelIOMgr::recvHydroBasedCounter(HydroBasedMsg *msg){
    molecule->numFixedRigidBonds += msg->numFixedRigidBonds;
    molecule->numFixedGroups += msg->numFixedGroups;

    if(++hydroMsgRecved == numInputProcs){
        msg->numFixedRigidBonds = molecule->numFixedRigidBonds;
        msg->numFixedGroups = molecule->numFixedGroups;
        CProxy_ParallelIOMgr pIO(thisgroup);
        pIO.bcastHydroBasedCounter(msg);
        hydroMsgRecved = 0;
    }else delete msg;
}

void ParallelIOMgr::bcastHydroBasedCounter(HydroBasedMsg *msg){
#ifdef MEM_OPT_VERSION
    //only the rank 0 in the SMP node update the Molecule object
    if(CmiMyRank()) {
        delete msg;
        return;
    }
    molecule->numFixedRigidBonds = msg->numFixedRigidBonds;
    molecule->numFixedGroups = msg->numFixedGroups;
    delete msg;

    if(!CkMyPe()) {
        iout << iINFO << "****************************\n";
        iout << iINFO << "STRUCTURE SUMMARY:\n";
        iout << iINFO << molecule->numAtoms << " ATOMS\n";
        iout << iINFO << molecule->numBonds << " BONDS\n";
        iout << iINFO << molecule->numAngles << " ANGLES\n";
        iout << iINFO << molecule->numDihedrals << " DIHEDRALS\n";
        iout << iINFO << molecule->numImpropers << " IMPROPERS\n";
        iout << iINFO << molecule->numCrossterms << " CROSSTERMS\n";
        iout << iINFO << molecule->numExclusions << " EXCLUSIONS\n";

            //****** BEGIN CHARMM/XPLOR type changes
        if ((molecule->numMultipleDihedrals) && (simParameters->paraTypeXplorOn)){
            iout << iINFO << molecule->numMultipleDihedrals 
             << " DIHEDRALS WITH MULTIPLE PERIODICITY (BASED ON PSF FILE)\n";
        }
        if ((molecule->numMultipleDihedrals) && (simParameters->paraTypeCharmmOn)){
            iout << iINFO << molecule->numMultipleDihedrals 
         << " DIHEDRALS WITH MULTIPLE PERIODICITY IGNORED (BASED ON PSF FILE) \n";
            iout << iINFO  
         << " CHARMM MULTIPLICITIES BASED ON PARAMETER FILE INFO! \n";
        }
            //****** END CHARMM/XPLOR type changes

        if (molecule->numMultipleImpropers){
            iout << iINFO << molecule->numMultipleImpropers 
                 << " IMPROPERS WITH MULTIPLE PERIODICITY\n";
        }

        if (simParameters->fixedAtomsOn)
           iout << iINFO << molecule->numFixedAtoms << " FIXED ATOMS\n";
        

        if (simParameters->rigidBonds)        
           iout << iINFO << molecule->numRigidBonds << " RIGID BONDS\n";        

        if (simParameters->fixedAtomsOn && simParameters->rigidBonds)        
           iout << iINFO << molecule->numFixedRigidBonds <<
                " RIGID BONDS BETWEEN FIXED ATOMS\n";

        iout << iINFO << molecule->num_deg_freedom(1)
             << " DEGREES OF FREEDOM\n";

        iout << iINFO << molecule->numHydrogenGroups << " HYDROGEN GROUPS\n";
        iout << iINFO << molecule->maxHydrogenGroupSize
            << " ATOMS IN LARGEST HYDROGEN GROUP\n";
        iout << iINFO << molecule->numMigrationGroups << " MIGRATION GROUPS\n";
        iout << iINFO << molecule->maxMigrationGroupSize
            << " ATOMS IN LARGEST MIGRATION GROUP\n";
        if (simParameters->fixedAtomsOn)
        {
           iout << iINFO << molecule->numFixedGroups <<
                " HYDROGEN GROUPS WITH ALL ATOMS FIXED\n";
        }
        
        iout << iINFO << "TOTAL MASS = " << totalMass << " amu\n"; 
        iout << iINFO << "TOTAL CHARGE = " << totalCharge << " e\n"; 

        BigReal volume = simParameters->lattice.volume();
        if ( volume ) {
            iout << iINFO << "MASS DENSITY = "
                << ((totalMass/volume) / 0.6022) << " g/cm^3\n";
            iout << iINFO << "ATOM DENSITY = "
                << (molecule->numAtoms/volume) << " atoms/A^3\n";
        }
    
        iout << iINFO << "*****************************\n";
        iout << endi;
        fflush(stdout);               
    }
#endif
}

void ParallelIOMgr::calcAtomsInEachPatch()
{
    if(myInputRank==-1) return;

    PatchMap *patchMap = PatchMap::Object();
    int numPatches = patchMap->numPatches();

    patchMap->initTmpPatchAtomsList();

    //each list contains the atom index to the initAtoms
    vector<int> *eachPatchAtomList = patchMap->getTmpPatchAtomsList();

    CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
    PatchMgr *patchMgr = pm.ckLocalBranch();

    int pid=0;
    const Lattice lattice = simParameters->lattice;
    for(int i=0; i<initAtoms.size(); i++) {
        InputAtom *atom = &(initAtoms[i]);
        if(!atom->isValid) continue;
        if(atom->isMP) {
            pid = patchMap->assignToPatch(atom->position, lattice);
        }
        eachPatchAtomList[pid].push_back(i);
    }

    CProxy_ParallelIOMgr pIO(thisgroup);

    int patchCnt = 0;
    for(int i=0; i<numPatches; i++) {
        int cursize = eachPatchAtomList[i].size();
        if(cursize>0) patchCnt++;
    }

    AtomsCntPerPatchMsg *msg = NULL;
    if(simParameters->fixedAtomsOn) {
        msg = new (patchCnt, patchCnt, patchCnt, 0)AtomsCntPerPatchMsg;
    } else {
        msg = new (patchCnt, patchCnt, 0, 0)AtomsCntPerPatchMsg;
    }

    msg->length = patchCnt;
    patchCnt = 0;
    for(int i=0; i<numPatches; i++) {
        int cursize = eachPatchAtomList[i].size();
        if(cursize>0) {
            if ( cursize > USHRT_MAX ) {
              char errstr[512];
              sprintf(errstr, "Patch %d exceeds %d atoms.", i, USHRT_MAX);
              NAMD_die(errstr);
            }
            msg->pidList[patchCnt] = i;
            msg->atomsCntList[patchCnt] = cursize;
            patchCnt++;
        }
    }

    if(simParameters->fixedAtomsOn) {
        patchCnt = 0;
        for(int i=0; i<numPatches; i++) {
            int cursize = eachPatchAtomList[i].size();
            if(cursize>0) {
                int fixedCnt = 0;
                for(int j=0; j<cursize; j++) {
                    int aid = eachPatchAtomList[i][j];
                    //atomFixed is either 0 or 1
                    fixedCnt += initAtoms[aid].atomFixed;
                }
                msg->fixedAtomsCntList[patchCnt] = fixedCnt;
                patchCnt++;
            }
        }
    }

    pIO[0].recvAtomsCntPerPatch(msg);

}

void ParallelIOMgr::recvAtomsCntPerPatch(AtomsCntPerPatchMsg *msg)
{
#ifdef MEM_OPT_VERSION
    PatchMap *patchMap = PatchMap::Object();
    for(int i=0; i<msg->length; i++) {
        int pid = msg->pidList[i];
        int oldNum = patchMap->numAtoms(pid);
        if ( oldNum + msg->atomsCntList[i] > USHRT_MAX ) {
          char errstr[512];
          sprintf(errstr, "Patch %d exceeds %d atoms.", pid, USHRT_MAX);
          NAMD_die(errstr);
        }
        patchMap->setNumAtoms(pid, oldNum+msg->atomsCntList[i]);
        if(simParameters->fixedAtomsOn) {
            oldNum = patchMap->numFixedAtoms(pid);
            patchMap->setNumFixedAtoms(pid, oldNum+msg->fixedAtomsCntList[i]);
        }
    }
    delete msg;

    if(++procsReceived == numInputProcs) {
        //print max PATCH info
        int maxAtoms = -1;
        int maxPatch = -1;
        int totalAtoms = 0;
        for(int i=0; i<patchMap->numPatches(); i++) {
            int cnt = patchMap->numAtoms(i);
            totalAtoms += cnt;
            if(cnt>maxAtoms) {
                maxAtoms = cnt;
                maxPatch = i;
            }
        }
        procsReceived = 0;
        iout << iINFO << "LARGEST PATCH (" << maxPatch <<
             ") HAS " << maxAtoms << " ATOMS\n" << endi;
        if ( totalAtoms !=  Node::Object()->molecule->numAtoms ) {
          char errstr[512];
          sprintf(errstr, "Incorrect atom count in void ParallelIOMgr::recvAtomsCntPerPatch: %d vs %d", totalAtoms, Node::Object()->molecule->numAtoms);
          NAMD_die(errstr);
        }
    }
#endif
}

void call_sendAtomsToHomePatchProcs(void *arg)
{
  ((ParallelIOMgr*)arg)->sendAtomsToHomePatchProcs();
}

void ParallelIOMgr::sendAtomsToHomePatchProcs()
{
#ifdef MEM_OPT_VERSION
    if(myInputRank==-1) return;

    if ( sendAtomsThread == 0 ) {
      sendAtomsThread = CthCreate((CthVoidFn)call_sendAtomsToHomePatchProcs,this,0);
      CthAwaken(sendAtomsThread);
      return;
    }
    sendAtomsThread = 0;
    numAcksOutstanding = 0;

    PatchMap *patchMap = PatchMap::Object();
    int numPatches = patchMap->numPatches();
    vector<int> *eachPatchAtomList = patchMap->getTmpPatchAtomsList();

    //each element (proc) contains the list of ids of patches which will stay
    //on that processor
    ResizeArray<int> *procList = new ResizeArray<int>[CkNumPes()];
    ResizeArray<int> pesToSend;
    for(int i=0; i<numPatches; i++) {
        if(eachPatchAtomList[i].size()==0) continue;
        int onPE = patchMap->node(i);
        if ( procList[onPE].size() == 0 ) pesToSend.add(onPE);
        procList[onPE].add(i);
    }

    Random(CkMyPe()).reorder(pesToSend.begin(),pesToSend.size());
    //CkPrintf("Pe %d ParallelIOMgr::sendAtomsToHomePatchProcs sending to %d pes\n",CkMyPe(),pesToSend.size());

    //go over every processor to send a message if necessary
    //TODO: Optimization for local home patches to save temp memory usage??? -CHAOMEI
    CProxy_ParallelIOMgr pIO(thisgroup);
    for(int k=0; k<pesToSend.size(); k++) {
        const int i = pesToSend[k];
        int len = procList[i].size();
        if(len==0) continue;

        if ( numAcksOutstanding >= 10 ) {
          //CkPrintf("Pe %d ParallelIOMgr::sendAtomsToHomePatchProcs suspending at %d of %d pes\n",CkMyPe(),k,pesToSend.size());
          //fflush(stdout);
          sendAtomsThread = CthSelf();
          CthSuspend();
        }
        ++numAcksOutstanding;

        //prepare a message to send
        int patchCnt = len;
        int totalAtomCnt = 0;
        for(int j=0; j<len; j++) {
            int pid = procList[i][j];
            int atomCnt = eachPatchAtomList[pid].size();
            totalAtomCnt += atomCnt;
        }

        MovePatchAtomsMsg *msg = new (patchCnt, patchCnt, totalAtomCnt, 0)MovePatchAtomsMsg;
        msg->from = CkMyPe();
        msg->patchCnt = patchCnt;
        int atomIdx = 0;
        for(int j=0; j<len; j++) {
            int pid = procList[i][j];
            int atomCnt = eachPatchAtomList[pid].size();
            msg->pidList[j] = pid;
            msg->sizeList[j] = atomCnt;
            for(int k=0; k<atomCnt; k++, atomIdx++) {
                int aid = eachPatchAtomList[pid][k];
                FullAtom one = initAtoms[aid];
                //HACK to re-sort the atom list after receiving the atom list on
                //home patch processor -Chao Mei
                one.hydVal = initAtoms[aid].hydList;
                msg->allAtoms[atomIdx] = one;
            }
        }
        pIO[i].recvAtomsToHomePatchProcs(msg);

        procList[i].clear();
    }

    //clean up to free space
    delete [] procList;
    patchMap->delTmpPatchAtomsList();

    //free the space occupied by the list that contains the input atoms
    initAtoms.clear();
#endif
}

void ParallelIOMgr::ackAtomsToHomePatchProcs()
{
  --numAcksOutstanding;
  if ( sendAtomsThread ) {
    CthAwaken(sendAtomsThread);
    sendAtomsThread = 0;
  }
}

void ParallelIOMgr::recvAtomsToHomePatchProcs(MovePatchAtomsMsg *msg)
{
    CProxy_ParallelIOMgr pIO(thisgroup);
    pIO[msg->from].ackAtomsToHomePatchProcs();

    if(!isOKToRecvHPAtoms) {
        prepareHomePatchAtomList();
        isOKToRecvHPAtoms = true;
    }

    int numRecvPatches = msg->patchCnt;
    int aid = 0;
    for(int i=0; i<numRecvPatches; i++) {
        int pid = msg->pidList[i];
        int size = msg->sizeList[i];
        int idx = binaryFindHPID(pid);
        for(int j=0; j<size; j++, aid++) {
            hpAtomsList[idx].add(msg->allAtoms[aid]);
        }
    }
    //CkPrintf("Pe %d recvAtomsToHomePatchProcs for %d patches %d atoms\n",CkMyPe(),numRecvPatches,aid);
    delete msg;
}

void ParallelIOMgr::prepareHomePatchAtomList()
{
    PatchMap *patchMap = PatchMap::Object();
    for(int i=0; i<patchMap->numPatches(); i++) {
        if(patchMap->node(i)==CkMyPe()) {
            hpIDList.add(i);
        }
    }
    if(hpIDList.size()>0)
        hpAtomsList = new FullAtomList[hpIDList.size()];
}

int ParallelIOMgr::binaryFindHPID(int pid)
{
    //hpIDList should be in increasing order!!
    int lIdx, rIdx;
    int retIdx = -1;

    rIdx=0;
    lIdx=hpIDList.size()-1;

    while(rIdx<=lIdx ) {
        int idx = (rIdx+lIdx)/2;
        int curPid = hpIDList[idx];
        if(pid>curPid) {
            //in the left
            rIdx = idx+1;
        } else if(pid<curPid) {
            //in the right
            lIdx = idx-1;
        } else {
            //found!
            retIdx = idx;
            break;
        }
    }
    CmiAssert(retIdx!=-1);
    return retIdx;
}

void ParallelIOMgr::createHomePatches()
{
#ifdef MEM_OPT_VERSION

    int assignedPids = PatchMap::Object()->numPatchesOnNode(CkMyPe());
    int numPids = hpIDList.size();
    if(numPids==0){
        //this node actually contains no homepatches
        if(assignedPids == 0) return; 

        //Entering the rare condition that all the homepatches this node has
        //are empty so that "recvAtomsToHomePatchProcs" is never called!
        //But we still need to create those empty homepatches!
        CmiAssert(isOKToRecvHPAtoms == false);        
        PatchMap *patchMap = PatchMap::Object();
        CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
        PatchMgr *patchMgr = pm.ckLocalBranch();
        for(int i=0; i<patchMap->numPatches(); i++) {
            if(patchMap->node(i)==CkMyPe()) {
                FullAtomList emptyone;
                patchMgr->createHomePatch(i, emptyone);
            }
        }        
        return;
    }

    CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
    PatchMgr *patchMgr = pm.ckLocalBranch();

    //go through the home patch list
    for(int i=0; i<numPids; i++) {
        int pid = hpIDList[i];

        //re-sort the atom list of this patch
        std::sort(hpAtomsList[i].begin(), hpAtomsList[i].end());
        Node::Object()->workDistrib->fillAtomListForOnePatch(pid, hpAtomsList[i]);
        patchMgr->createHomePatch(pid, hpAtomsList[i]);
    }

    hpIDList.clear();   
    delete [] hpAtomsList;

    hpAtomsList = NULL;
#endif
}

void ParallelIOMgr::freeMolSpace()
{
#ifdef MEM_OPT_VERSION
    molecule->delAtomNames();
    molecule->delChargeSpace();

    //???TODO NOT SURE WHETHER freeEnergyOn is support in MEM_OPT_VERSION
    //-CHAOMEI
    if(!CkMyPe() && !simParameters->freeEnergyOn)
        molecule->delMassSpace();

    molecule->delFixedAtoms();
#endif
}

//The initial distribution of atoms:
//(avgNum+1),...,(avgNum+1),avgNum,...,avgNum where
//avgNum equals to the total number of atoms in the system
//divided by number of input processors --Chao Mei
int ParallelIOMgr::numMyAtoms(int rank, int numProcs)
{
    if(rank==-1) return -1;
    int avgNum = molecule->numAtoms/numProcs;
    int remainder = molecule->numAtoms%numProcs;
    if(rank<remainder) return avgNum+1;
    else return avgNum;
}

int ParallelIOMgr::atomRank(int atomID, int numProcs)
{
    int avgNum = molecule->numAtoms/numProcs;
    int remainder = molecule->numAtoms%numProcs;
    int midLimit = remainder*(avgNum+1);
    int idx;
    if(atomID<midLimit) {
        idx = atomID/(avgNum+1);
    } else {
        idx = remainder+(atomID-midLimit)/avgNum;
    }
    return idx;
}

void ParallelIOMgr::getMyAtomsRange(int &lowerIdx, int &upperIdx, int rank, int numProcs)
{
    if(rank==-1) {
        //just to make sure upperIdx-lowerIdx+1 == -1
        lowerIdx=-1;
        upperIdx=-3;
        return;
    }
    int avgNum = molecule->numAtoms/numProcs;
    int remainder = molecule->numAtoms%numProcs;
    if(rank<remainder) {
        lowerIdx = rank*(avgNum+1);
        upperIdx = lowerIdx+avgNum;
    } else {
        int midLimit = remainder*(avgNum+1);
        lowerIdx = midLimit+(rank-remainder)*avgNum;
        upperIdx = lowerIdx+avgNum-1;
    }
}

void ParallelIOMgr::receivePositions(CollectVectorVarMsg *msg)
{
#ifdef MEM_OPT_VERSION
    int ready = midCM->receivePositions(msg);
    if(ready) {
        CProxy_CollectionMaster cm(mainMaster);
        cm.receiveOutputPosReady(msg->seq);
    }
    delete msg;
#endif
}

void ParallelIOMgr::receiveVelocities(CollectVectorVarMsg *msg)
{
#ifdef MEM_OPT_VERSION
    int ready = midCM->receiveVelocities(msg);
    if(ready) {
        CProxy_CollectionMaster cm(mainMaster);
        cm.receiveOutputVelReady(msg->seq);        
    }
    delete msg;
#endif
}

void ParallelIOMgr::receiveForces(CollectVectorVarMsg *msg)
{
#ifdef MEM_OPT_VERSION
    int ready = midCM->receiveForces(msg);
    if(ready) {
        CProxy_CollectionMaster cm(mainMaster);
        cm.receiveOutputForceReady(msg->seq);        
    }
    delete msg;
#endif
}


void ParallelIOMgr::disposePositions(int seq, double prevT)
{
#ifdef MEM_OPT_VERSION
	double iotime = CmiWallTimer();
    midCM->disposePositions(seq);
	iotime = CmiWallTimer()-iotime+prevT;

#if OUTPUT_SINGLE_FILE    
	//Token-based file output
    if(myOutputRank == getMyOutputGroupHighestRank()) {
        //notify the CollectionMaster to start the next round
        CProxy_CollectionMaster cm(mainMaster);
        cm.startNextRoundOutputPos(iotime);
    } else {
        CProxy_ParallelIOMgr io(thisgroup);
        io[outputProcArray[myOutputRank+1]].disposePositions(seq, iotime);
    }
#else
	//notify the CollectionMaster to start the next round
	CProxy_CollectionMaster cm(mainMaster);
	cm.startNextRoundOutputPos(iotime);
#endif

#endif
}

void ParallelIOMgr::disposeVelocities(int seq, double prevT)
{
#ifdef MEM_OPT_VERSION
	double iotime = CmiWallTimer();
    midCM->disposeVelocities(seq);
	iotime = CmiWallTimer()-iotime+prevT;
    
#if OUTPUT_SINGLE_FILE
	//Token-based file output
    if(myOutputRank==getMyOutputGroupHighestRank()) {
        //notify the CollectionMaster to start the next round
        CProxy_CollectionMaster cm(mainMaster);
        cm.startNextRoundOutputVel(iotime);
    } else {
        CProxy_ParallelIOMgr io(thisgroup);
        io[outputProcArray[myOutputRank+1]].disposeVelocities(seq, iotime);
    }
#else
	//notify the CollectionMaster to start the next round
	CProxy_CollectionMaster cm(mainMaster);
	cm.startNextRoundOutputVel(iotime);	
#endif

#endif
}

void ParallelIOMgr::disposeForces(int seq, double prevT)
{
#ifdef MEM_OPT_VERSION
	double iotime = CmiWallTimer();
    midCM->disposeForces(seq);
	iotime = CmiWallTimer()-iotime+prevT;
    
#if OUTPUT_SINGLE_FILE
	//Token-based file output
    if(myOutputRank==getMyOutputGroupHighestRank()) {
        //notify the CollectionMaster to start the next round
        CProxy_CollectionMaster cm(mainMaster);
        cm.startNextRoundOutputForce(iotime);
    } else {
        CProxy_ParallelIOMgr io(thisgroup);
        io[outputProcArray[myOutputRank+1]].disposeForces(seq, iotime);
    }
#else
	//notify the CollectionMaster to start the next round
	CProxy_CollectionMaster cm(mainMaster);
	cm.startNextRoundOutputForce(iotime);	
#endif

#endif
}


void ParallelIOMgr::wrapCoor(int seq, Lattice lat)
{
#ifdef MEM_OPT_VERSION
    coorInstance = midCM->getReadyPositions(seq);

    coorInstance->lattice = lat; //record the lattice to use for wrapAll/Water!
    int fromAtomID = coorInstance->fromAtomID;
    int toAtomID = coorInstance->toAtomID;

    //only reference copies
    ResizeArray<Vector> &data = coorInstance->data;
    ResizeArray<FloatVector> &fdata = coorInstance->fdata;
    //if both data and fdata are not empty, they contain exact values, the only
    //difference lies in their precisions. Therefore, we only need to compute 
    //the higher precision coordinate array. -Chao Mei
    int dsize = data.size();    
    int numMyAtoms = toAtomID-fromAtomID+1;
    tmpCoorCon = new Vector[numMyAtoms];    
    ClusterCoorElem one;
    //1. compute wrapped coordinates locally 
    for(int i=0; i<numMyAtoms; i++){
        tmpCoorCon[i] = 0.0;
        int cid = clusterID[i];
        if(cid<fromAtomID){
            //on output procs ahead of me
            one.clusterId = cid;
            ClusterCoorElem *ret = remoteCoors.find(one);
            if(ret==NULL){
                if(dsize==0) 
                    one.dsum = fdata[i];
                else 
                    one.dsum = data[i];
                                
                remoteCoors.add(one);                 
            }else{
                if(dsize==0) 
                    ret->dsum += fdata[i];
                else 
                    ret->dsum += data[i];               
            }
        }else{
            if(dsize==0) 
                tmpCoorCon[cid-fromAtomID] += fdata[i];
            else 
                tmpCoorCon[cid-fromAtomID] += data[i];
        }
    }

    //2. Prepare to send msgs to remote output procs to reduce coordinates 
    //values of a cluster
    CmiAssert(numRemoteClusters == remoteCoors.size());
    numCSMAck = 0; //set to 0 to prepare recving the final coor update
    CProxy_ParallelIOMgr pIO(thisgroup);
    ClusterCoorSetIter iter(remoteCoors);
    for(iter=iter.begin(); iter!=iter.end(); iter++){
        ClusterCoorMsg *msg = new ClusterCoorMsg;
        msg->srcRank = myOutputRank;
        msg->clusterId = iter->clusterId;
        msg->dsum = iter->dsum;
        int dstRank = atomRankOnOutput(iter->clusterId);
        pIO[outputProcArray[dstRank]].recvClusterCoor(msg);
    }
    
    //Just send a local NULL msg to indicate the local wrapping
    //coordinates has finished.
    recvClusterCoor(NULL);
#endif
}

//On the output proc, it's possible (could be a rare case) that recvClusterCoor
//is executed before wrapCoor, so we need to make sure the local tmpCoorCon has
//been calculated. This is why (numRemoteReqs+1) msgs are expected as the 
//additional one is sent by itself when it finishes wrapping coordinates.
// --Chao Mei
void ParallelIOMgr::recvClusterCoor(ClusterCoorMsg *msg){
    //only add the msg from remote procs
    if(msg!=NULL) ccmBuf.add(msg);

    //include a msg sent by itself
    if(++numReqRecved == (numRemoteReqs+1)){
        numReqRecved = 0;
        integrateClusterCoor();
    }
}

void ParallelIOMgr::integrateClusterCoor(){
#ifdef MEM_OPT_VERSION
    int fromIdx = coorInstance->fromAtomID;
    int toIdx = coorInstance->toAtomID;
    for(int i=0; i<ccmBuf.size(); i++){
        ClusterCoorMsg *msg = ccmBuf[i];
        int lidx = msg->clusterId - fromIdx;        
        tmpCoorCon[lidx] += msg->dsum;
    }

    //send back those msgs
    CProxy_ParallelIOMgr pIO(thisgroup);
    for(int i=0; i<ccmBuf.size(); i++){
        ClusterCoorMsg *msg = ccmBuf[i];
        int lidx = msg->clusterId - fromIdx;        
        if(simParameters->wrapAll || isWater[lidx]) {
            Lattice *lat = &(coorInstance->lattice);
            Vector coni = tmpCoorCon[lidx]/clusterSize[lidx];
            msg->dsum = (simParameters->wrapNearest ?
                  lat->wrap_nearest_delta(coni) : lat->wrap_delta(coni));
        }else{
            msg->dsum = 0.0;
        }
        pIO[outputProcArray[msg->srcRank]].recvFinalClusterCoor(msg);
    }
    ccmBuf.resize(0);    

    //It's possible that recvFinalClusterCoor is executed before integrateClusterCoor
    //on this processor (the other proc executes faster in doing wrapCoor). So avoid
    //this msg race, do send a message to itself to participate the reduction.
    if(numRemoteClusters!=0){
        recvFinalClusterCoor(NULL);
    } else {
        //this output proc is ready has the sum of coordinates of each cluster
        //on it, so it is ready to do the final wrap coor computation        
        int numMyAtoms = toIdx-fromIdx+1;
        ResizeArray<Vector> &data = coorInstance->data;
        ResizeArray<FloatVector> &fdata = coorInstance->fdata;
        for(int i=0; i<numMyAtoms; i++){
            if(!simParameters->wrapAll && !isWater[i]) continue;
            int lidx = clusterID[i]-fromIdx;
            if(lidx==i){
                //the head atom of the cluster
                Lattice *lat = &(coorInstance->lattice);
                Vector coni = tmpCoorCon[lidx]/clusterSize[lidx];
                tmpCoorCon[lidx] = (simParameters->wrapNearest ?
                  lat->wrap_nearest_delta(coni) : lat->wrap_delta(coni));
            }
            if(data.size()) data[i] += tmpCoorCon[lidx]; 
            //no operator += (FloatVector, Vector)
            if(fdata.size()) fdata[i] = fdata[i] + tmpCoorCon[lidx]; 
        }
        
        delete [] tmpCoorCon;
        tmpCoorCon = NULL;
        CProxy_CollectionMaster cm(mainMaster);
        cm.wrapCoorFinished();
    }
#endif    
}

void ParallelIOMgr::recvFinalClusterCoor(ClusterCoorMsg *msg){
#ifdef MEM_OPT_VERSION
    if(msg!=NULL){
        //only process the message sent from other procs!
        ClusterCoorElem one(msg->clusterId);
        ClusterCoorElem *ret = remoteCoors.find(one);
        ret->dsum = msg->dsum;
        delete msg;
    }
    
    if(++numCSMAck == (numRemoteClusters+1)){        
        //final wrap coor computation
        int fromIdx = coorInstance->fromAtomID;
        int toIdx = coorInstance->toAtomID;
        int numMyAtoms = toIdx-fromIdx+1;
        ResizeArray<Vector> &data = coorInstance->data;
        ResizeArray<FloatVector> &fdata = coorInstance->fdata;
        ClusterCoorElem tmp;
        for(int i=0; i<numMyAtoms; i++){
            if(!simParameters->wrapAll && !isWater[i]) continue;
            int cid = clusterID[i];
            int lidx = cid-fromIdx;
            if(lidx<0){
                //this cid should be inside remoteCoors
                tmp.clusterId = cid;
                ClusterCoorElem *fone = remoteCoors.find(tmp);
                if(data.size()) data[i] += fone->dsum; 
                if(fdata.size()) fdata[i] = fdata[i] + fone->dsum; 
            }else{
                if(lidx==i){
                    Lattice *lat = &(coorInstance->lattice);
                    Vector coni = tmpCoorCon[lidx]/clusterSize[lidx];
                    tmpCoorCon[lidx] = (simParameters->wrapNearest ?
                    lat->wrap_nearest_delta(coni) : lat->wrap_delta(coni));
                }
                if(data.size()) data[i] += tmpCoorCon[lidx]; 
                if(fdata.size()) fdata[i] = fdata[i] + tmpCoorCon[lidx];
            }
        }

        delete [] tmpCoorCon;
        tmpCoorCon = NULL;
        CProxy_CollectionMaster cm(mainMaster);
        cm.wrapCoorFinished();
        numCSMAck = 0;
        remoteCoors.clear();
    }
#endif
}
#include "ParallelIOMgr.def.h"
