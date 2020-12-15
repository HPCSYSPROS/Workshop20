/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PROXYMGR_H
#define PROXYMGR_H


#include "charm++.h"
#include "envelope.h"

#include "main.h"
#include "NamdTypes.h"
#include "PatchTypes.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "ProcessorPrivate.h"
#include "ProxyMgr.decl.h"

extern int proxySendSpanning, proxyRecvSpanning;
extern int proxySpanDim;
extern int inNodeProxySpanDim;

#if CMK_PERSISTENT_COMM
#define USE_PERSISTENT_TREE                  0
#endif

class ProxyGBISP1ResultMsg: public CMessage_ProxyGBISP1ResultMsg {
  public:
    int destPe;
    int origPe;
    PatchID patch;
    GBReal *psiSum;
    int psiSumLen;// = numAtoms
};
class ProxyGBISP2DataMsg: public CMessage_ProxyGBISP2DataMsg {
public:
    int origPe;
    int destPe;
    PatchID patch;
    Real *bornRad;//numAtoms
    int bornRadLen;
  };
class ProxyGBISP2ResultMsg: public CMessage_ProxyGBISP2ResultMsg {
public:
    int destPe;
    int origPe;
    PatchID patch;
    GBReal *dEdaSum;
    int dEdaSumLen;//numAtoms
  };
class ProxyGBISP3DataMsg: public CMessage_ProxyGBISP3DataMsg {
public:
    int origPe;
    int destPe;
    PatchID patch;
    Real *dHdrPrefix;//numAtoms
    int dHdrPrefixLen;
  };

class RegisterProxyMsg : public CMessage_RegisterProxyMsg {
public:
  NodeID node;
  PatchID patch;
};

class UnregisterProxyMsg : public CMessage_UnregisterProxyMsg {
public:
  NodeID node;
  PatchID patch;
};

//1. This class represents for both msg types: one that
//is originally known as ProxyAllMsg which is sent
//at the step where atoms migrate; and the other is
//sent during the steps between two migrations.
//2. In the case of memory optimized version, the scenario
//becomes tricky as load balancer will move compute objects
//around so that new ProxyPatches will be created where
//the CompAtomExt list information is not available. If
//the step immediately after the load balancing is a normal
//step, then the CompAtomExt list info has to be resent by
//the HomePatch. Because of the current Proxy msg communication
//scheme where msg is sent to ProxyMgr first, and then retransmitted
//to ProxyPatches, there's overhead when we want to resend CompAtomExt
//list as not all the ProxyPatches that are managed by ProxyMgr are
//newly created ProxyPatches. 
//--Chao Mei
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
class ProxyDataMsg_not_32_byte : public CMessage_ProxyDataMsg {
#else
class ProxyDataMsg : public CMessage_ProxyDataMsg {
#endif
public:
  PatchID patch;
  Flags flags;

  int plLen;

  CompAtom *positionList;
  int avgPlLen;
  CompAtom *avgPositionList;
  // BEGIN LA
  int vlLen;
  CompAtom *velocityList;
  // END LA

  Real *intRadList;// gbis atom intrinsic radii

  int *lcpoTypeList;// LCPO atom type

  //1. The following field will be only
  //useful for memory optimized version.
  //2. In normal case, adding this field only
  //increases the msg length by 4 bytes which
  //can be ignored considering the current fast
  //communication network
  //--Chao Mei
  int plExtLen;
  CompAtomExt *positionExtList;
  CudaAtom *cudaAtomList;

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
  //In smp layer, the couter for msg creation and process of communication
  //thread is not included in the quiescence detection process. In addition,
  //the immediate messages from other nodes are executed on the communication
  //thread. If inside the process of immediate messages, some normal Charm++
  //messages sent out which will be processed on worker threads. Then QD will
  //be a problem that the process of the normal messages sent from communication
  //thread is recorded, but the creation of such messages (although recorded
  //in the comm thread) is virtually not recorded, i.e., not visible the 
  //QD process. So we need to artificially increase the QD counter to 
  //compensate for aforementioned msg creation loss.
  //The idea is to use the following variable to indicate the normal message
  //is sent from the communication thread inside a processing of immediate
  //message. If the variable is set, then we should increase the QD counter.
  //Chao Mei
  char isFromImmMsgCall; //hack for imm msg with QD in SMP 
#endif

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int numWaterAtoms;  // Number of atoms in positionList (from start)
	                //   that are part of water hydrogen groups.
  #endif

};

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
class ProxyDataMsg : public ProxyDataMsg_not_32_byte {
  // Adding padding bytes to make sure that positionList is
  // 32-byte aligned which usually gives better cache performance.
  char padding[(32-(sizeof(envelope)+sizeof(ProxyDataMsg_not_32_byte))%32)%32];
};
class assert_ProxyDataMsg {
  char assert_sizeof_envelope_is_multiple_of_ptr[(sizeof(envelope)%sizeof(void*))?-1:1];
  char assert_sizeof_ProxyDataMsg_is_multiple_of_32[((sizeof(envelope)+sizeof(ProxyDataMsg))%32)?-1:1];
};
#endif


class ProxyResultMsg : public CMessage_ProxyResultMsg {
public:
  NodeID node;
  PatchID patch;
  ForceList *forceList[Results::maxNumForces];
  static void* pack(ProxyResultMsg *msg);
  static ProxyResultMsg* unpack(void *ptr);
private:
  ForceList forceListInternal[Results::maxNumForces];
};

class ProxyResultVarsizeMsg: public CMessage_ProxyResultVarsizeMsg{
public:
    NodeID node;
    PatchID patch;
    int flLen[Results::maxNumForces];   

    Force *forceArr;
    //Indicate the position of the force list that has zero value
    //which is not recorded in the above force array.
    char *isZero;

    //add padding bytes to make sure the beginning 
    //of force arrays is 8-byte aligned as it is originally.
    //Therefore, we have to put the forceArr field as
    //the first variable of varsize array type
    char padding[(8-(sizeof(envelope)+sizeof(NodeID)+sizeof(PatchID)+sizeof(int)*Results::maxNumForces+2*sizeof(void *))%8)%8];   

    //The length of "fls" is Results::maxNumForces
    static ProxyResultVarsizeMsg *getANewMsg(NodeID nid, PatchID pid, int prioSize, ForceList *fls); 
};

class assert_ProxyResultVarsizeMsg {
  char assert_sizeof_envelope_is_multiple_of_ptr[(sizeof(envelope)%sizeof(void*))?-1:1];
  char assert_sizeof_ProxyResultVarsizeMsg_is_multiple_of_8[((sizeof(envelope)+sizeof(ProxyResultVarsizeMsg))%8)?-1:1];
};

class ProxyNodeAwareSpanningTreeMsg: public CMessage_ProxyNodeAwareSpanningTreeMsg{
public:
    PatchID patch;
    NodeID procID;
    int numNodesWithProxies;
    int *numPesOfNode;
    int *allPes;

    static ProxyNodeAwareSpanningTreeMsg *getANewMsg(PatchID pid, NodeID nid, proxyTreeNode *tree, int size);

    //For debug
    void printOut(char *tag);
};

class ProxyCombinedResultRawMsg : public CMessage_ProxyCombinedResultRawMsg {
public:
        int nodeSize;
        NodeID *nodes;

        #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
        //since this msg may be processed by comm thread in the smp mode,
        //this variable helps comm thread to find which proc will actually process it.
        NodeID destPe;
        #if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
        //Mainly for QD in the presence of the optimization of using immediate
        //message. Refer to the explanation from ProxyDataMsg for the same 
        //variable. --Chao Mei
        char isFromImmMsgCall;
        #endif
        #endif
        PatchID patch;

        int flLen[Results::maxNumForces];
        char *isForceNonZero;
        //The beginning address of this variable should be 8-byte aligned!!! -Chao Mei
        Force *forceArr;
};

class ProxyCombinedResultMsg : public CMessage_ProxyCombinedResultMsg {
public:
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  //since this msg may be processed by comm thread in the smp mode,
  //this variable helps comm thread to find which proc will actually process it.
  NodeID destPe;
  #if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
  //Mainly for QD in the presence of the optimization of using immediate
  //message. Refer to the explanation from ProxyDataMsg for the same 
  //variable. --Chao Mei
  char isFromImmMsgCall;
  #endif
  #endif
  PatchID patch;
  NodeIDList nodes;
  ForceList *forceList[Results::maxNumForces];
  static ProxyCombinedResultRawMsg* toRaw(ProxyCombinedResultMsg *msg);
  static ProxyCombinedResultMsg* fromRaw(ProxyCombinedResultRawMsg *msg);
private:
  ForceList forceListInternal[Results::maxNumForces];
};

class ProxySpanningTreeMsg : public CMessage_ProxySpanningTreeMsg {
public:
  PatchID patch;
  NodeID  node;
  NodeIDList tree;
  static void* pack(ProxySpanningTreeMsg *msg);
  static ProxySpanningTreeMsg* unpack(void *ptr);
};

class ProxyPatch;
class PatchMap;

struct ProxyElem {
  ProxyElem() : proxyPatch(0) { };
  ProxyElem(PatchID pid) : patchID(pid), proxyPatch(0) { };
  ProxyElem(PatchID pid, ProxyPatch *p) : patchID(pid), proxyPatch(p) { };

  int hash() const { return patchID; }
  int operator==(const ProxyElem & pe) const { return patchID == pe.patchID; }

  PatchID patchID;
  ProxyPatch *proxyPatch;
};

typedef UniqueSet<ProxyElem> ProxySet;
typedef UniqueSetIter<ProxyElem> ProxySetIter;

class ProxyTree {       // keep track of the spanning trees
  public:
    int proxyMsgCount;
    NodeIDList *proxylist;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    //a node-aware spanning tree array, each element of which
    //is a spanning tree for all proxies of a patch
    proxyTreeNodeList *naTrees;
#else
    NodeIDList *trees;
    int *sizes;
#endif
    
  public:
    ProxyTree() {
      proxyMsgCount = 0;
      proxylist = NULL;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
      naTrees = NULL;
#else
      trees = NULL;
      sizes = NULL;
#endif      
    }
    ~ProxyTree() {
    }
};

class ProxyMgr : public CBase_ProxyMgr
{
public:
  ProxyMgr();
  ~ProxyMgr();

  void removeProxies(void);
  void removeUnusedProxies(void);
  void createProxies(void);

  void createProxy(PatchID pid);
  void removeProxy(PatchID pid);

  void registerProxy(PatchID pid);
  void recvRegisterProxy(RegisterProxyMsg *);

  void unregisterProxy(PatchID pid);
  void recvUnregisterProxy(UnregisterProxyMsg *);

  void setSendSpanning();
  int  getSendSpanning();

  void setRecvSpanning();
  int  getRecvSpanning();

  void setProxyTreeBranchFactor(int dim);

  void buildProxySpanningTree();
  void sendSpanningTrees();
  void sendSpanningTreeToHomePatch(int pid, int *tree, int n);
  void recvSpanningTreeOnHomePatch(int pid, int *tree, int n);
  void sendSpanningTree(ProxySpanningTreeMsg *);
  void recvSpanningTree(ProxySpanningTreeMsg *);

  void sendNodeAwareSpanningTreeToHomePatch(int pid, proxyTreeNode *tree, int n);
  void recvNodeAwareSpanningTreeOnHomePatch(ProxyNodeAwareSpanningTreeMsg *msg);
  void sendNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *);
  void recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *);
  //set the proxy patch's parent field
  void recvNodeAwareSTParent(int patch, int parent);

  void buildProxySpanningTree2();               // centralized version
  void sendProxies(int pid, int *list, int n);
  void recvProxies(int pid, int *list, int n);
  void recvPatchProxyInfo(PatchProxyListMsg *msg);

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  void buildNodeAwareSpanningTree0();
  static void buildSinglePatchNodeAwareSpanningTree(PatchID pid, NodeIDList &proxyList, 
                                                    proxyTreeNodeList &ptnTree);
#else
  void buildSpanningTree0();
#endif

  void sendResults(ProxyResultVarsizeMsg *);
  void recvResults(ProxyResultVarsizeMsg *);
  void sendResults(ProxyResultMsg *);
  void recvResults(ProxyResultMsg *);
  void sendResults(ProxyCombinedResultMsg *);

  void sendResult(ProxyGBISP1ResultMsg *);//psiSum
  void recvResult(ProxyGBISP1ResultMsg *);
  void recvData(   ProxyGBISP2DataMsg *);
  void sendResult(ProxyGBISP2ResultMsg *);//dEdaSum
  void recvResult(ProxyGBISP2ResultMsg *);
  void recvData(   ProxyGBISP3DataMsg *);

  void recvResults(ProxyCombinedResultRawMsg *);
  void recvImmediateResults(ProxyCombinedResultRawMsg *);

  void sendProxyData(ProxyDataMsg *, int, int*);
  void recvImmediateProxyData(ProxyDataMsg *);
  void recvProxyData(ProxyDataMsg *);

  void sendProxyAll(ProxyDataMsg *, int, int*);
  void recvImmediateProxyAll(ProxyDataMsg *);
  void recvProxyAll(ProxyDataMsg *);

  static ProxyMgr *Object() { return CkpvAccess(ProxyMgr_instance); }
  
  int numProxies() { return proxySet.size(); }

  static int nodecount;
  ProxyTree &getPtree();
 
private:
  ProxySet proxySet;
  ProxyTree ptree;

  void printProxySpanningTree();
};

struct ProxyListInfo{
	int patchID;
	int numProxies;
	int *proxyList; //record which PE the proxy is on 
};

class PatchProxyListMsg: public CMessage_PatchProxyListMsg {
public:
	int numPatches;
	int *patchIDs;
	int *proxyListLen;
	int *proxyPEs;

public:
	PatchProxyListMsg(int num) { numPatches = num; }
	static PatchProxyListMsg *createPatchProxyListMsg(PatchProxyListMsg **bufs, int bufSize, ProxyListInfo *info, int size);
};

class NodeProxyMgr : public CBase_NodeProxyMgr
{
private:
/*The following vars are for node-aware spanning tree if NodeProxyMgr is used*/
	//Potential TODO: change it to a hashtable if necessary
    proxyTreeNode **proxyInfo;
    int numPatches;

    CkGroupID localProxyMgr; //a charm Group variable
    PatchMap **localPatchMaps;

/* The following vars are for managing sending proxy list of a patch to PE 0*/
	int parentNode; //-1 for root node
	int numKidNodes;
	int kidRecved;
	//about local home patches
	int numHomePatches; //the number of homepatches on this node
	int homepatchRecved;
	ProxyListInfo *localProxyLists;
	PatchProxyListMsg **remoteProxyLists;
	CmiNodeLock localDepositLock;
	CmiNodeLock remoteDepositLock;

public:
    NodeProxyMgr(){
        proxyInfo = NULL;
        numPatches = 0;
        localPatchMaps = new PatchMap *[CkMyNodeSize()];

		parentNode = -1;
		numKidNodes = 0;
		kidRecved = 0;
		numHomePatches = 0;
		homepatchRecved = 0;
		localProxyLists = NULL;
		remoteProxyLists = NULL;
		localDepositLock = CmiCreateLock();
		remoteDepositLock = CmiCreateLock();
    }
    ~NodeProxyMgr(){
        for(int i=0; i<numPatches; i++) {
            delete proxyInfo[i];
        }
        delete [] proxyInfo;
        delete [] localPatchMaps;

		CmiDestroyLock(localDepositLock);
		CmiDestroyLock(remoteDepositLock);
    }

    void createProxyInfo(int numPs){
        numPatches = numPs;
        proxyInfo = new proxyTreeNode *[numPs];
        memset(proxyInfo, 0, sizeof(proxyTreeNode *)*numPs);
    }
    void registerPatch(int patchID, int numPes, int *pes);
    proxyTreeNode *getPatchProxyInfo(int patchID){
        return proxyInfo[patchID];
    }

    void registerLocalProxyMgr(CkGroupID one){
        localProxyMgr = one;
    }
    const CkGroupID &getLocalProxyMgr(){
        return localProxyMgr;
    }
    void registerLocalPatchMap(int rank, PatchMap *one){
        localPatchMaps[rank] = one;
    }
    PatchMap *getLocalPatchMap(int rank){
        return localPatchMaps[rank];
    }   

    void recvImmediateProxyData(ProxyDataMsg *msg);
    void recvImmediateProxyAll(ProxyDataMsg *msg);
    void recvImmediateResults(ProxyCombinedResultRawMsg *);

	//initialize the spanning tree of home patches
	void createSTForHomePatches(PatchMap *pmap);
	//direct call from local home patches
	void sendProxyList(int pid, int *plist, int size);
	//remote call to send this node's proxy list info
	void sendProxyListInfo(PatchProxyListMsg *msg);
	void contributeToParent();
};

#endif /* PATCHMGR_H */

