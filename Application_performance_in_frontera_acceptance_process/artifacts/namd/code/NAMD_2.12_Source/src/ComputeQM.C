/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "dcdlib.h"

// #define DEBUGM
// #define DEBUG_QM

#ifdef DEBUG_QM
  #define DEBUGM
#endif

#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeQM.h"
#include "ComputeQMMgr.decl.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "ComputePme.h"
#include "ComputePmeMgr.decl.h"

#include <fstream>
#include <iomanip>

#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define mkdir(X,Y) _mkdir(X)
#define S_ISDIR(X) ((X) & S_IFDIR)
#endif

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif

#define QMLENTYPE 1
#define QMRATIOTYPE 2

#define QMPCSCHEMENONE 1
#define QMPCSCHEMEROUND 2
#define QMPCSCHEMEZERO 3

#define QMPCSCALESHIFT 1
#define QMPCSCALESWITCH 2

struct ComputeQMAtom {
    Position position;
    float charge;
    int id;
    Real qmGrpID;
    int homeIndx;
    int vdwType;
    ComputeQMAtom() : position(0), charge(0), id(-1), qmGrpID(-1), 
                        homeIndx(-1), vdwType(0) {}
    ComputeQMAtom(Position posInit, float chrgInit, int idInit, 
                  Real qmInit, int hiInit, int newvdwType) {
        position = posInit;
        charge = chrgInit;
        id = idInit;
        qmGrpID = qmInit;
        homeIndx = hiInit;
        vdwType = newvdwType;
    }
    ComputeQMAtom(const ComputeQMAtom &ref) {
        position = ref.position;
        charge = ref.charge;
        id = ref.id;
        qmGrpID = ref.qmGrpID;
        homeIndx = ref.homeIndx;
        vdwType = ref.vdwType;
    }
};

#define PCMODEUPDATESEL 1
#define PCMODEUPDATEPOS 2
#define PCMODECUSTOMSEL 3

class QMCoordMsg : public CMessage_QMCoordMsg {
public:
  int sourceNode;
  int numAtoms;
  int timestep;
  int numPCIndxs;
  int pcSelMode;
  ComputeQMAtom *coord;
  int *pcIndxs;
};

struct pntChrgDist {
    int index;
    Real dist;
    pntChrgDist() : index(-1), dist(0) {};
    pntChrgDist(int newIndx, Real newDist) {
        index = newIndx;
        dist = newDist;
    }
    int operator <(const pntChrgDist& v) {return dist < v.dist;}
};

struct ComputeQMPntChrg {
    Position position;
    float charge;
    int id;
    Real qmGrpID;
    int homeIndx;
    Real dist;
    Mass mass;
    int vdwType;
    ComputeQMPntChrg() : position(0), charge(0), id(-1), qmGrpID(-1), 
        homeIndx(-1), dist(0), mass(0), vdwType(0) {}
    ComputeQMPntChrg(Position posInit, float chrgInit, int idInit, 
                     Real qmInit, int hiInit, Real newDist, 
                     Mass newM, int newvdwType) {
        position = posInit;
        charge = chrgInit;
        id = idInit;
        qmGrpID = qmInit;
        homeIndx = hiInit;
        dist = newDist;
        mass = newM;
        vdwType = newvdwType;
    }
    ComputeQMPntChrg(const ComputeQMPntChrg &ref) {
        position = ref.position;
        charge = ref.charge;
        id = ref.id;
        qmGrpID = ref.qmGrpID;
        homeIndx = ref.homeIndx;
        dist = ref.dist;
        mass = ref.mass;
        vdwType = ref.vdwType;
    }
    
    bool operator <(const ComputeQMPntChrg& ref) {
        return (id < ref.id);
    }
    bool operator ==(const ComputeQMPntChrg& ref) {
        return (id == ref.id) ;
    }
};

class QMPntChrgMsg : public CMessage_QMPntChrgMsg {
public:
  int sourceNode;
  int numAtoms;
  ComputeQMPntChrg *coord;
};

struct QMForce {
  int replace;
  Force force;
  int homeIndx;
  float charge;
  int id;
  QMForce() : replace(0), force(0), homeIndx(-1), charge(0), id(-1) {;}
};

class QMGrpResMsg : public CMessage_QMGrpResMsg {
public:
  int grpIndx;
  BigReal energyOrig; // Original QM Energy
  BigReal energyCorr; // Corrected energy due to PME
  BigReal virial[3][3];
  int numForces;
  QMForce *force;
};

class QMForceMsg : public CMessage_QMForceMsg {
public:
  bool PMEOn;
  BigReal energy;
  BigReal virial[3][3];
  int numForces;
  int numLssDat;
  QMForce *force;
  LSSSubsDat *lssDat;
};

# define QMATOMTYPE_DUMMY 0
# define QMATOMTYPE_QM 1
# define QMPCTYPE_IGNORE 0
# define QMPCTYPE_CLASSICAL 1
# define QMPCTYPE_EXTRA 2

#define QMLSSQMRES 1
#define QMLSSCLASSICALRES 2

struct idIndxStr {
    int ID; // Global atom ID
    int indx; // Atom index in the vector received from other ranks.
    
    idIndxStr() {
        ID = -1;
        indx = -1;
    };
    idIndxStr(int newID, int newIndx) {
        ID = newID;
        indx = newIndx;
    }
    idIndxStr(const idIndxStr& ref) {
        ID = ref.ID ;
        indx = ref.indx;
    }
    
    idIndxStr &operator=(const idIndxStr& ref) {
        ID = ref.ID ;
        indx = ref.indx;
        return *this;
    }
    
    bool operator==(const idIndxStr& ref) const {
        return (ID == ref.ID && indx == ref.indx);
    }
    bool operator!=(const idIndxStr& ref) const {
        return !(*this == ref);
    }
    
    // We sort the atoms based on their global ID so that the same
    // order is kept when replacing classical and QM residues.
    // We are assuming here that every time a residue is added to a PDB
    // file, its atoms are placed in the same order.
    bool operator<(const idIndxStr& ref) {return ID < ref.ID;}
};

struct lssDistSort {
    int type;
    Real dist;
//     std::vector<int> idVec;
//     std::vector<int> indxVec;
    SortedArray<idIndxStr> idIndx;
    
    lssDistSort() {
        idIndx.resize(0);
    };
    lssDistSort(int newType, Real newDist) {
        type = newType;
        dist = newDist;
        idIndx.resize(0);
    }
    lssDistSort(const lssDistSort& ref) {
        type = ref.type;
        dist = ref.dist;
        idIndx.resize(idIndx.size());
        for (int i=0; i<ref.idIndx.size(); i++) idIndx.insert(ref.idIndx[i]) ;
        idIndx.sort();
    }
    
    lssDistSort& operator=(const lssDistSort& ref) {
        type = ref.type;
        dist = ref.dist;
        idIndx.resize(idIndx.size());
        for (int i=0; i<ref.idIndx.size(); i++) idIndx.insert(ref.idIndx[i]) ;
        idIndx.sort();
        
        return *this;
    }
    
    bool operator==(const lssDistSort& ref) {
        bool returnVal = true;
        
        if (! (type == ref.type && dist == ref.dist))
            return false;
        
        if (idIndx.size() != ref.idIndx.size())
            return false;
        
        for (int i=0; i<ref.idIndx.size(); i++) {
            if (idIndx[i] != ref.idIndx[i])
                return false;
        }
        
        return returnVal;
    }
    
    bool operator<(const lssDistSort& ref) {return dist < ref.dist;}
    
} ;

struct QMAtomData {
    Position position;
    float charge;
    int id;
    int bountToIndx; // So as to implement an "exclude 1-2" functionality.
    int type;
    char element[3];
    Real dist;
    QMAtomData() : position(0), charge(0), id(-1), bountToIndx(-1), 
                   type(-1), dist(0) {}
    QMAtomData(Position posInit, float chrgInit, int idInit, 
               int bountToIndxInit, int newType, 
               char *elementInit, Real newDist) {
        position = posInit;
        charge = chrgInit;
        id = idInit;
        bountToIndx = bountToIndxInit;
        type = newType;
        strncpy(element,elementInit,3);
        dist = newDist;
    }
};

class QMGrpCalcMsg: public CMessage_QMGrpCalcMsg {
public:
    int grpIndx;
    int peIter;
    int numQMAtoms ;
    int numAllAtoms ; // Including dummy atoms.
    int numRealPntChrgs;
    int numAllPntChrgs; // Inlcuding point charges created to handle QM-MM bonds.
    Real charge, multiplicity;
    BigReal constants;
    bool secProcOn ;
    bool prepProcOn ;
    bool PMEOn;
    bool switching ;
    int switchType;
    Real cutoff, swdist;
    int pcScheme ;
    BigReal PMEEwaldCoefficient;
    int qmAtmChrgMode;
    char baseDir[256], execPath[256], secProc[256], prepProc[256];
    QMAtomData *data;
    char *configLines;
};

struct dummyData {
    // Position of dummy atom.
    Position pos ;
    // Indicates to which QM atom is this dummy atom bound.
    // This is indexed in the context of the number of QM atoms in a QM group.
    int qmInd ;
    // Indicates the index of this bond, used to get the element of the dummy atom.
    int bondIndx;
    dummyData(Position newPos, int newQmInd, int newIndx) { 
        pos = newPos;
        qmInd = newQmInd;
        bondIndx = newIndx;
    }
} ;

struct LSSDataStr {
    int resIndx ;
    Mass mass ;
    LSSDataStr(int newResIndx, Mass newMass) {
        resIndx = newResIndx;
        mass = newMass;
    }
} ;

typedef std::pair<int, LSSDataStr> atmLSSData;
typedef std::map<int, LSSDataStr> LSSDataMap;

typedef std::pair<int, Mass> refLSSData;
typedef std::map<int, Mass> LSSRefMap;

typedef std::vector<ComputeQMPntChrg> QMPCVec ;
typedef std::vector<ComputeQMAtom> QMAtmVec ;

class ComputeQMMgr : public CBase_ComputeQMMgr {
public:
    ComputeQMMgr();
    ~ComputeQMMgr();

    void setCompute(ComputeQM *c) { QMCompute = c; }
    
    // Function called on pe ZERO by all pe's.
    // Receives partial QM coordinates form the QM atoms in each pe and aggregates
    // them all in a single structure.
    // This function initializes some variables and allcoates memory.
    void recvPartQM(QMCoordMsg*) ;
    
    void recvFullQM(QMCoordMsg*) ;
    
    void recvPntChrg(QMPntChrgMsg*) ;
    
    void calcMOPAC(QMGrpCalcMsg *) ;
    void calcORCA(QMGrpCalcMsg *) ;
    void calcUSR(QMGrpCalcMsg *) ;
    
    void storeQMRes(QMGrpResMsg *) ;
    
    void procQMRes();
    
    void recvForce(QMForceMsg*) ;
    
    SortedArray<LSSSubsDat> &get_subsArray() { return subsArray; } ;
private:
    ComputeQMMgr_SDAG_CODE
    
    CProxy_ComputeQMMgr QMProxy;
    ComputeQM *QMCompute;

    int numSources;
    int numArrivedQMMsg, numArrivedPntChrgMsg, numRecQMRes;
    QMCoordMsg **qmCoordMsgs;
    QMPntChrgMsg **pntChrgCoordMsgs;
    
    std::vector<int> qmPEs ;
    
    ComputeQMAtom *qmCoord;
    QMPCVec pntChrgMsgVec;
    QMForce *force;

    int numQMGrps;

    int numAtoms;

    int qmAtmIter;

    int numQMAtms;
    
    Bool noPC ;
    int meNumMMIndx;
    
    int qmPCFreq;
    // In case we are skiping steps between point charge re-assignment (i.e. qmPCFreq > 0),
    // we keep a list in rank zero whihch is replaced every X steps, and 
    // then re-sent to all PEs every "X+1" steps. This is done because atoms
    // may move between pathes and PEs between PC reassignments, if they are 
    // sufficiently long or unlucky.
    int *pcIDList ;
    int pcIDListSize;
    Bool resendPCList;
    SortedArray<ComputeQMPntChrg> pcDataSort;
    
    int replaceForces;
    int bondValType;
    SimParameters * simParams;
    Molecule *molPtr;
    // ID for each QM group.
    Real *grpID;
    Real *grpChrg, *grpMult;
    
    Bool qmLSSOn ;
    int qmLSSFreq;
    int qmLSSResSize;
    int *qmLSSSize;
    int *qmLSSIdxs;
    Mass *qmLSSMass;
    int lssTotRes;
    // This will hold one map per QM group. Each will map an atom ID with
    // the local residue index for the solvent molecule it belongs to, 
    // and its atomic mass (in case COM calculation is being performed).
    ResizeArray<LSSDataMap> grpIDResNum ;
    int *qmLSSRefIDs;
    int *qmLSSRefSize;
    // This will hold one map per QM group. Each will map an atom ID with
    // its atomic mass. This is used to calculate the reference COM to which
    // solvent molecules will be compared for distance calculation.
    ResizeArray<LSSRefMap> lssGrpRefMass ;
    // Current positions of each LSS RESIDUE
    Position *lssPos;
    Mass lssResMass;
    SortedArray<LSSSubsDat> subsArray;
    
    BigReal ewaldcof, pi_ewaldcof;
    
    // Accumulates the enery of different QM groups. This is the energy passed to NAMD.
    BigReal totalEnergy;
    BigReal totVirial[3][3];

    int dcdOutFile, dcdPosOutFile;
    Real *outputData ;
    int timeStep ;
    
    void procBonds(int numBonds,
                     const int *const qmGrpBondsGrp,
                     const int *const qmMMBondedIndxGrp,
                     const int *const *const chargeTarget,
                     const int *const numTargs,
                     const QMPCVec grpPntChrgVec,
                     QMPCVec &grpAppldChrgVec) ;
    
    void pntChrgSwitching(QMGrpCalcMsg* msg) ;
    
    void lssPrepare() ;
    void lssUpdate(int grpIter, QMAtmVec &grpQMAtmVec, QMPCVec &grpPntChrgVec);
    
    #ifdef DEBUG_QM
    void Write_PDB(std::string Filename, const QMGrpCalcMsg *dataMsg);
    void Write_PDB(std::string Filename, const QMCoordMsg *dataMsg);
    #endif
};

ComputeQMMgr::ComputeQMMgr() :
  QMProxy(thisgroup), QMCompute(0), numSources(0), numArrivedQMMsg(0), 
  numArrivedPntChrgMsg(0), numRecQMRes(0), qmCoordMsgs(0), pntChrgCoordMsgs(0),
  qmCoord(0), force(0), numAtoms(0), dcdOutFile(0), outputData(0), pcIDListSize(0) {
      
  CkpvAccess(BOCclass_group).computeQMMgr = thisgroup;
  
}

ComputeQMMgr::~ComputeQMMgr() {
    
    if (qmCoordMsgs != 0) {
        for ( int i=0; i<numSources; ++i ) {
            if (qmCoordMsgs[i] != 0)
                delete qmCoordMsgs[i];
        }
        delete [] qmCoordMsgs;
    }
    
    if (pntChrgCoordMsgs != 0) {
        for ( int i=0; i<numSources; ++i ) {
            if (pntChrgCoordMsgs[i] != 0)
                delete pntChrgCoordMsgs[i];
        }
        delete [] pntChrgCoordMsgs;
    }
    
    if (qmCoord != 0)
        delete [] qmCoord;
    
    if (force != 0)
        delete [] force;
    
    
    if (dcdOutFile != 0)
        close_dcd_write(dcdOutFile);
    
    if (dcdPosOutFile != 0)
        close_dcd_write(dcdPosOutFile);
    
    if (outputData != 0)
        delete [] outputData;
    
    if (lssPos != NULL)
        delete [] lssPos;
}

SortedArray<LSSSubsDat> &lssSubs(ComputeQMMgr *mgr) { 
    return mgr->get_subsArray(); 
} ;

ComputeQM::ComputeQM(ComputeID c) :
  ComputeHomePatches(c), oldForces(0)
{
  CProxy_ComputeQMMgr::ckLocalBranch(
	CkpvAccess(BOCclass_group).computeQMMgr)->setCompute(this);

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  
}

ComputeQM::~ComputeQM()
{
    if (oldForces != 0)
        delete [] oldForces;
}

void ComputeQM::initialize()
{
    ComputeHomePatches::initialize();
    
    // Get some pointers that will be used in the future,
    simParams = Node::Object()->simParameters ;
    molPtr = Node::Object()->molecule;
    // Future comes fast...
    qmAtomGroup = molPtr->get_qmAtomGroup() ;
    
    noPC = molPtr->get_noPC();
    if (noPC) {
        meNumMMIndx = molPtr->get_qmMeNumBonds();
        meMMindx = molPtr->get_qmMeMMindx();
        meQMGrp =molPtr->get_qmMeQMGrp();
        
        for (int i=0; i<meNumMMIndx; i++)
            meQMBonds.add(meMMQMGrp(meMMindx[i], meQMGrp[i]));
    }
    
    numQMAtms = molPtr->get_numQMAtoms();
    qmAtmChrg = molPtr->get_qmAtmChrg() ;
    qmAtmIndx = molPtr->get_qmAtmIndx() ;
    
    numQMGrps = molPtr->get_qmNumGrps();
    qmGrpIDArray = molPtr->get_qmGrpID() ;
    
    cutoff = simParams->cutoff;
    
    customPC = simParams->qmCustomPCSel;
    if (customPC) {
        
        customPCLists.resize(numQMGrps);
        
        int totCustPCs = molPtr->get_qmTotCustPCs();
        const int *custPCSizes = molPtr->get_qmCustPCSizes();
        const int *customPCIdxs = molPtr->get_qmCustomPCIdxs();
        
        int minI = 0, maxI = 0, grpIter = 0;
        
        while (grpIter < numQMGrps) {
            
            maxI += custPCSizes[grpIter];
            
            for (size_t i=minI; i < maxI; i++) {
                
                customPCLists[grpIter].add( customPCIdxs[i] ) ;
            }
            
            // Minimum index gets the size of the current group
            minI += custPCSizes[grpIter];
            // Group iterator gets incremented
            grpIter++;
            
        }
    }
}


void ComputeQM::doWork()
{
    CProxy_ComputeQMMgr QMProxy(CkpvAccess(BOCclass_group).computeQMMgr);
    
    ResizeArrayIter<PatchElem> ap(patchList);
    
    int timeStep ;
    
    #ifdef DEBUG_QM
    DebugM(4,"----> Initiating QM work on rank " << CkMyPe() <<
    " with " << patchList.size() << " patches." << std::endl );
    #endif
    
    patchData.clear();
    
    // Determines hou many QM atoms are in this node.
    int numLocalQMAtoms = 0, localMM1atoms = 0;
    for (ap = ap.begin(); ap != ap.end(); ap++) {
        CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
        int localNumAtoms = (*ap).p->getNumAtoms();
        
        patchData.push_back(patchDataStrc((*ap).positionBox,
                    (*ap).positionBox->open(),
                    (*ap).p) );
        
        for(int i=0; i<localNumAtoms; ++i) {
            if ( qmAtomGroup[xExt[i].id] > 0 ) {
                numLocalQMAtoms++;
            }
            else if (meNumMMIndx > 0) {
                // If we have a mechanical embedding simulation with no 
                // point charges, but with QM-MM bonds, we look for the MM1 atoms 
                // here and pass them directly. No need for an extra round of 
                // comms just to get atoms we already know we need.
                
                auto retIt = meQMBonds.find(meMMQMGrp(xExt[i].id));
                
                if (retIt != NULL) {
                    localMM1atoms++;
                }
            }
        }
        
        timeStep = (*ap).p->flags.step ;
    }
    
    DebugM(4, "Node " << CkMyPe() << " has " << numLocalQMAtoms
    << " QM atoms." << std::endl) ;
    #ifdef DEBUG_QM
    if (localMM1atoms > 0)
        DebugM(4, "Node " << CkMyPe() << " has " << localMM1atoms
            << " MM1 atoms." << std::endl) ;
    #endif
    // Send QM atoms
    
    // This pointer will be freed in rank zero.
    QMCoordMsg *msg = new (numLocalQMAtoms + localMM1atoms,0, 0) QMCoordMsg;
    msg->sourceNode = CkMyPe();
    msg->numAtoms = numLocalQMAtoms + localMM1atoms;
    msg->timestep = timeStep;
    ComputeQMAtom *data_ptr = msg->coord;
    
    // Build a message with coordinates, charge, atom ID and QM group
    // for all QM atoms in the node, then send it to node zero.
    int homeCounter = 0;
    for (auto pdIt = patchData.begin(); pdIt != patchData.end(); pdIt++ ) {
        
        CompAtom *x = (*pdIt).compAtomP;
        const CompAtomExt *xExt =  (*pdIt).homePatchP->getCompAtomExtInfo();
        int localNumAtoms =  (*pdIt).homePatchP->getNumAtoms();
        const FullAtom* fullAtms = (*pdIt).homePatchP->getAtomList().const_begin() ;
        
        for(int i=0; i<localNumAtoms; ++i) {
            
            if ( qmAtomGroup[xExt[i].id] > 0 ) {
                
                Real charge = 0;
                
                for (int qi=0; qi<numQMAtms; qi++) {
                    if (qmAtmIndx[qi] == xExt[i].id) {
                        charge = qmAtmChrg[qi];
                        break;
                    }
                }
                
                data_ptr->position = patchList[0].p->lattice.reverse_transform(x[i].position, 
                                                          fullAtms[i].transform) ;
//                 data_ptr->charge = x[i].charge;
                data_ptr->charge = charge;
                data_ptr->id = xExt[i].id;
                data_ptr->qmGrpID = qmAtomGroup[xExt[i].id] ;
                data_ptr->homeIndx = homeCounter;
                data_ptr->vdwType = fullAtms[i].vdwType;
                
                ++data_ptr;
            }
            else if (meNumMMIndx > 0) {
                // If we have a mechanical embedding simulation with no 
                // point charges, but with QM-MM bonds, we look for the MM1 atoms 
                // connected here and pass them directly.
                
                auto retIt = meQMBonds.find(meMMQMGrp(xExt[i].id));
                
                if (retIt != NULL) {
                    
                    DebugM(3,"Found atom " << retIt->mmIndx << "," << retIt->qmGrp << "\n" );
                    
                    data_ptr->position = patchList[0].p->lattice.reverse_transform(x[i].position, 
                                                          fullAtms[i].transform) ;
                    // We force the charge to be zero since we would only get
                    // here if mechanical embeding was selected.
                    data_ptr->charge = 0;
                    data_ptr->id = xExt[i].id;
                    data_ptr->qmGrpID = retIt->qmGrp ;
                    data_ptr->homeIndx = homeCounter;
                    
                    // We re-use this varaible to indicate this is an MM1 atom.
                    data_ptr->vdwType = -1;
                    
                    ++data_ptr;
                    
                }
                
            }
                
            homeCounter++;
            
        }
        
    }
    
    if (noPC) {
        for (auto pdIt = patchData.begin(); pdIt != patchData.end(); pdIt++ ) {
            CompAtom *x = (*pdIt).compAtomP;
            (*pdIt).posBoxP->close(&x);
        }
    }
    
    QMProxy[0].recvPartQM(msg);
    
}


void ComputeQMMgr::recvPartQM(QMCoordMsg*msg)
{
    // In the first (ever) step of the simulation, we allocate arrays that
    // will be used throughout the simulation, and get some important numbers.
    if ( ! numSources ) {
        DebugM(4,"Initializing ComputeQMMgr variables." << std::endl);
        numSources = (PatchMap::Object())->numNodesWithPatches();
        
        DebugM(4,"There are " << numSources << " nodes with patches." << std::endl);
        qmCoordMsgs = new QMCoordMsg*[numSources];
        for ( int i=0; i<numSources; ++i ) { 
            qmCoordMsgs[i] = 0;
        }
        
        // Prepares the allocation for the recvPntChrg function.
        
        DebugM(4,"Getting data from molecule and simParameters." << std::endl);
        
        molPtr = Node::Object()->molecule;
        simParams = Node::Object()->simParameters;
        
        numAtoms = molPtr->numAtoms;
        force = new QMForce[numAtoms];
        
        numQMAtms = molPtr->get_numQMAtoms();
        qmAtmIter = 0;
        
        noPC = simParams->qmNoPC ;
        meNumMMIndx = molPtr->get_qmMeNumBonds();
        if (noPC && meNumMMIndx == 0) {
            pntChrgCoordMsgs = NULL;
        }
        else {
            pntChrgCoordMsgs = new QMPntChrgMsg*[numSources];
            for ( int i=0; i<numSources; ++i ) { 
                pntChrgCoordMsgs[i] = 0;
            }
        }
        
        qmPCFreq = molPtr->get_qmPCFreq();
        resendPCList = false;
        
        grpID = molPtr->get_qmGrpID() ;
        bondValType = simParams->qmBondValType;
        
        numQMGrps = molPtr->get_qmNumGrps();
        
        grpChrg = molPtr->get_qmGrpChrg() ;
        
        grpMult = molPtr->get_qmGrpMult() ;
        
        qmLSSOn = simParams->qmLSSOn ;
        if (qmLSSOn) {
            qmLSSFreq = molPtr->get_qmLSSFreq() ;
            qmLSSSize = molPtr->get_qmLSSSize() ;
            qmLSSIdxs = molPtr->get_qmLSSIdxs() ;
            qmLSSMass = molPtr->get_qmLSSMass() ;
            qmLSSResSize = molPtr->get_qmLSSResSize() ;
            qmLSSRefIDs = molPtr->get_qmLSSRefIDs() ;
            qmLSSRefSize = molPtr->get_qmLSSRefSize() ;
            
            lssPrepare();
        }
        
        numArrivedQMMsg = 0 ;
        numArrivedPntChrgMsg = 0 ;
        
        qmCoord = new ComputeQMAtom[numQMAtms];
        
        replaceForces = 0;
        if (molPtr->get_qmReplaceAll()) {
            replaceForces = 1;
        }
        
        DebugM(4,"Initializing DCD file for charge information." << std::endl);
        
        // Initializes output DCD file for charge information.
        if (simParams->qmOutFreq > 0) {
            std::string dcdOutFileStr;
            dcdOutFileStr.clear();
            dcdOutFileStr.append(simParams->outputFilename) ;
            dcdOutFileStr.append(".qdcd") ;
            dcdOutFile = open_dcd_write(dcdOutFileStr.c_str()) ;
            
            if (dcdOutFile == DCD_FILEEXISTS) {
                iout << iERROR << "DCD file " << dcdOutFile << " already exists!!\n" << endi;
                NAMD_err("Could not write QM charge DCD file.");
            } else if (dcdOutFile < 0) {
                iout << iERROR << "Couldn't open DCD file " << dcdOutFile << ".\n" << endi;
                NAMD_err("Could not write QM charge DCD file.");
            } else if (! dcdOutFile) {
                NAMD_bug("ComputeQMMgr::recvPartQM open_dcd_write returned fileid of zero");
            }
            
            timeStep = simParams->firstTimestep;
            
            int NSAVC, NFILE, NPRIV, NSTEP;
            NSAVC = simParams->qmOutFreq;
            NPRIV = timeStep;
            NSTEP = NPRIV - NSAVC;
            NFILE = 0;
            
            //  Write out the header
            int ret_code = write_dcdheader(dcdOutFile, dcdOutFileStr.c_str(),
            numQMAtms, NFILE, NPRIV, NSAVC, NSTEP,
            simParams->dt/TIMEFACTOR, 0);

            if (ret_code<0) {
                NAMD_err("Writing of DCD header failed!!");
            }
            
            // The DCD write function takes 3 independent arrays for X, Y and Z
            // coordinates, but we allocate one and send in the pieces.
            outputData = new Real[3*numQMAtms];
        }
        
        DebugM(4,"Initializing DCD file for position information." << std::endl);
        // Initializes output DCD file for position information.
        if (simParams->qmPosOutFreq > 0) {
            std::string dcdPosOutFileStr;
            dcdPosOutFileStr.clear();
            dcdPosOutFileStr.append(simParams->outputFilename) ;
            dcdPosOutFileStr.append(".QMonly.dcd") ;
            dcdPosOutFile = open_dcd_write(dcdPosOutFileStr.c_str()) ;
            
            if (dcdPosOutFile == DCD_FILEEXISTS) {
                iout << iERROR << "DCD file " << dcdPosOutFile << " already exists!!\n" << endi;
                NAMD_err("Could not write QM charge DCD file.");
            } else if (dcdPosOutFile < 0) {
                iout << iERROR << "Couldn't open DCD file " << dcdPosOutFile << ".\n" << endi;
                NAMD_err("Could not write QM charge DCD file.");
            } else if (! dcdPosOutFile) {
                NAMD_bug("ComputeQMMgr::recvPartQM open_dcd_write returned fileid of zero");
            }
            
            timeStep = simParams->firstTimestep;
            
            int NSAVC, NFILE, NPRIV, NSTEP;
            NSAVC = simParams->qmOutFreq;
            NPRIV = timeStep;
            NSTEP = NPRIV - NSAVC;
            NFILE = 0;
            
            //  Write out the header
            int ret_code = write_dcdheader(dcdPosOutFile, dcdPosOutFileStr.c_str(),
            numQMAtms, NFILE, NPRIV, NSAVC, NSTEP,
            simParams->dt/TIMEFACTOR, 0);

            if (ret_code<0) {
                NAMD_err("Writing of DCD header failed!!");
            }
            
            // The DCD write function takes 3 independent arrays for X, Y and Z
            // coordinates, but we allocate one and send in the pieces.
            outputData = new Real[3*numQMAtms];
        }
        
        // Prepares list of PEs which will run the QM software
        int simsPerNode = simParams->qmSimsPerNode ;
        int zeroNodeSize = CmiNumPesOnPhysicalNode(0);
        
        // Check if the node has enought PEs to run the requested number of simulations.
        if ( zeroNodeSize < simsPerNode ) {
            iout << iERROR << "There are " << zeroNodeSize << " cores pernode, but "
            << simsPerNode << " QM simulations per node were requested.\n" << endi ;
            NAMD_die("Error preparing QM simulations.");
        }
        
        DebugM(4,"Preparing PE list for QM execution.\n");
        qmPEs.clear(); // Making sure its empty.
        
        int numNodes = CmiNumPhysicalNodes();
        int numPlacedQMGrps = 0;
        int placedOnNode = 0;
        int nodeIt = 0 ;
        
        // The default is to only run on rank zero.
        if ( simsPerNode <= 0 ) {
            qmPEs.push_back(0);
            numPlacedQMGrps = 1;
        }
        
        while ( (numPlacedQMGrps < numQMGrps) && (simsPerNode > 0) ) {
            
            // If we searched all nodes, break the loop.
            if (nodeIt == numNodes) {
                break;
            }
            
            int *pelist;
            int nodeSize;
            CmiGetPesOnPhysicalNode(CmiPhysicalNodeID(nodeIt), &pelist, &nodeSize);
            
            DebugM(4,"Checking node " << nodeIt +1 << " out of " << numNodes
            << " (" << nodeSize << " PEs: " << pelist[0] << " to " 
            << pelist[nodeSize-1] << ")." << std::endl );
            
            for ( int i=0; i<nodeSize; ++i ) {
                
                // Check if the PE has patches. We only run on PEs with patches!
                if ( (PatchMap::Object())->numPatchesOnNode(pelist[i]) == 0 ) {
                    DebugM(1,"PE " << pelist[i] << " has no patches! We'll skip it." 
                    << std::endl);
                    continue ;
                }
                
                // Add the target PE on the target node to the list
                // of PEs which will carry QM simulations.
                qmPEs.push_back(pelist[i]);
                
                DebugM(1,"Added PE to QM execution list: " << pelist[i]  << "\n");
                
                numPlacedQMGrps++;
                placedOnNode++;
                
                if (placedOnNode == simsPerNode) {
                    DebugM(1,"Placed enought simulations on this node.\n");
                    break;
                }
                
                
            }
            
            nodeIt++ ;
            placedOnNode = 0;
        }
        
        if ( numPlacedQMGrps < numQMGrps ) {
            iout << iWARN << "Could not compute all QM groups in parallel.\n" << endi ;
        }
        
        iout << iINFO << "List of ranks running QM simulations: " << qmPEs[0] ;
        for (int i=1; i < qmPEs.size(); i++) {
            iout << ", " << qmPEs[i] ;
        }
        iout << ".\n" << endi;
        
    }
    
    DebugM(1,"Receiving from rank " << msg->sourceNode
    << " a total of " << msg->numAtoms << " QM atoms." << std::endl);
    
    // In case we are NOT using point charges but there are QM-MM bonds,
    // test each QM message for MM1 atoms.
    if (meNumMMIndx > 0) {
        
        ResizeArray< ComputeQMAtom > tmpAtm;
        ComputeQMAtom *data_ptr = msg->coord;
        
        for (int i=0; i<msg->numAtoms; i++) {
            if (data_ptr[i].vdwType < 0) {
                tmpAtm.add(data_ptr[i]) ;
            }
        }
        
        QMPntChrgMsg *pntChrgMsg = new (tmpAtm.size(), 0) QMPntChrgMsg;
        pntChrgMsg->sourceNode = msg->sourceNode ;
        pntChrgMsg->numAtoms = tmpAtm.size() ;
        ComputeQMPntChrg* newPCData = pntChrgMsg->coord ;
        
        QMCoordMsg *newMsg = msg;
        
        if (tmpAtm.size() > 0) {
            
            newMsg = new (msg->numAtoms - tmpAtm.size(),0, 0) QMCoordMsg;
            newMsg->sourceNode = msg->sourceNode ;
            newMsg->numAtoms = msg->numAtoms - tmpAtm.size() ;
            newMsg->timestep = msg->timestep ;
            ComputeQMAtom *newMsgData = newMsg->coord;
            
            for (int i=0; i<msg->numAtoms; i++) {
                if (data_ptr[i].vdwType < 0) {
                    newPCData->position = data_ptr[i].position ;
                    newPCData->charge = data_ptr[i].charge ;
                    newPCData->id = data_ptr[i].id ;
                    newPCData->qmGrpID = data_ptr[i].qmGrpID ;
                    newPCData->homeIndx = data_ptr[i].homeIndx ;
                    newPCData->dist = 0 ;
                    newPCData->mass = 0 ;
                    newPCData->vdwType = 0 ;
                    newPCData++;
                }
                else {
                    *newMsgData = data_ptr[i] ;
                    newMsgData++;
                }
            }
            
            delete msg;
            
        }
        
        qmCoordMsgs[numArrivedQMMsg] = newMsg;
        ++numArrivedQMMsg;
        
        pntChrgCoordMsgs[numArrivedPntChrgMsg] = pntChrgMsg;
        ++numArrivedPntChrgMsg;
    }
    else {
        qmCoordMsgs[numArrivedQMMsg] = msg;
        ++numArrivedQMMsg;
    }
    
    if ( numArrivedQMMsg < numSources ) 
        return;
    
    // Now that all messages arrived, get all QM positions.
    for (int msgIt=0; msgIt < numArrivedQMMsg; msgIt++){
        
        DebugM(1, "Getting positions for " << qmCoordMsgs[msgIt]->numAtoms 
        << " QM atoms in this message." << std::endl);
        
        for ( int i=0; i < qmCoordMsgs[msgIt]->numAtoms; ++i ) {
            qmCoord[qmAtmIter] = qmCoordMsgs[msgIt]->coord[i];
            qmAtmIter++;
        }
    }
    
    if (qmAtmIter != numQMAtms) {
        iout << iERROR << "The number of QM atoms received (" << qmAtmIter
        << ") is different than expected: " << numQMAtms << "\n" << endi;
        NAMD_err("Problems broadcasting QM atoms.");
    }
    
    // Resets the counter for the next step.
    numArrivedQMMsg = 0;
    qmAtmIter = 0;
    
    timeStep = qmCoordMsgs[0]->timestep;
    
    // Makes sure there is no trash or old info in the force array.
    // This is where we accumulate forces from the QM software and our own
    // Coulomb forces. It will have info on QM atoms and point charges only.
    for (int i=0; i<numAtoms; ++i ) {
        force[i].force = 0;  // this assigns 0 (ZERO) to all componenets of the vector.
        force[i].replace = 0;  // By default, no atom has all forces replaced.
        force[i].homeIndx = -1;  // To prevent errors from sliping.
        force[i].charge = 0;
        force[i].id = i;        // Initializes the variable since not all forces are sent to all ranks.
    }
    
    
    for (int i=0; i<numQMAtms; i++) {
        // Each force receives the home index of its atom with respect to the 
        // local set of atoms in each node.
        if (force[qmCoord[i].id].homeIndx != -1
            && force[qmCoord[i].id].homeIndx != qmCoord[i].homeIndx
        ) {
            iout << iERROR << "Overloading QM atom " 
            << qmCoord[i].id << "; home index: " 
            << force[qmCoord[i].id].homeIndx << " with " << qmCoord[i].homeIndx
            << "\n" << endi ;
            NAMD_die("Error preparing QM simulations.");
        }
        
        force[qmCoord[i].id].homeIndx = qmCoord[i].homeIndx;
        // Each force on QM atoms has an indicator so NAMD knows if all 
        // NAMD forces should be erased and completely overwritten by the
        // external QM forces.
        force[qmCoord[i].id].replace = replaceForces;
    }
    
    if (noPC) {
        // this pointer should be freed in rank zero, after receiving it.
        QMPntChrgMsg *pntChrgMsg = new (0, 0) QMPntChrgMsg;
        pntChrgMsg->sourceNode = CkMyPe();
        pntChrgMsg->numAtoms = 0;
        
        QMProxy[0].recvPntChrg(pntChrgMsg);
    }
    else {
        // The default mode is to update the poitn charge selection at every step.
        int pcSelMode = PCMODEUPDATESEL;
        
        int msgPCListSize = 0;
        // We check wether point charges are to be selected with a stride.
        if ( qmPCFreq > 0 ) {
            if (timeStep % qmPCFreq == 0) {
                // Marks that the PC list determined in this step will
                // be broadcast on the *next* step, and kept for the following
                // qmPCFreq-1 steps.
                resendPCList = true;
                
                // Clears the list since we don't know how many charges 
                // will be selected this time.
                delete [] pcIDList;
            }
            else {
                // If the PC selection is not to be updated in this step,
                // inform the nodes that the previous list (or the updated
                // list being sent in this message) is to be used and only 
                // updated positions will be returned.
                pcSelMode = PCMODEUPDATEPOS;
                
                // If this is the first step after a PC re-selection, all 
                // ranks receive the total list, since atoms may move between
                // PEs in between PC re-assignments (if they are far enought apart
                // or if you are unlucky)
                if (resendPCList) {
                    msgPCListSize = pcIDListSize;
                    resendPCList = false;
                }
            }
        }
        
        // In case we are using custom selection of point charges, indicate the mode.
        if (simParams->qmCustomPCSel)
            pcSelMode = PCMODECUSTOMSEL;
        
        DebugM(1,"Broadcasting current positions of QM atoms.\n");
        for ( int j=0; j < numSources; ++j ) {
            // This pointer will be freed in the receiving rank.
            QMCoordMsg *qmFullMsg = new (numQMAtms, msgPCListSize, 0) QMCoordMsg;
            qmFullMsg->sourceNode = CkMyPe();
            qmFullMsg->numAtoms = numQMAtms;
            qmFullMsg->pcSelMode = pcSelMode;
            qmFullMsg->numPCIndxs =  msgPCListSize;
            ComputeQMAtom *data_ptr = qmFullMsg->coord;
            int *msgPCListP = qmFullMsg->pcIndxs;
            
            for (int i=0; i<numQMAtms; i++) {
                data_ptr->position = qmCoord[i].position;
                data_ptr->charge = qmCoord[i].charge;
                data_ptr->id = qmCoord[i].id;
                data_ptr->qmGrpID = qmCoord[i].qmGrpID;
                data_ptr->homeIndx = -1; // So there is no mistake that there is no info here.
                data_ptr++;
            }
            
            for (int i=0; i<msgPCListSize; i++) {
                msgPCListP[i] = pcIDList[i] ;
            }
            
            #ifdef DEBUG_QM
            if (j == 0)
                Write_PDB("/home/melomcr/Research/NAMD_QMMM/TestSystem/qmMsg.pdb", qmFullMsg) ;
            #endif
            
            // The messages are deleted later, we will need them.
            QMProxy[qmCoordMsgs[j]->sourceNode].recvFullQM(qmFullMsg);
        }
    }
}

void ComputeQMMgr::recvFullQM(QMCoordMsg* qmFullMsg) {
    
    if (subsArray.size() > 0)
        subsArray.clear();
    
    QMCompute->processFullQM(qmFullMsg);
}

typedef std::map<Real,BigReal> GrpDistMap ;

void ComputeQM::processFullQM(QMCoordMsg* qmFullMsg) {
    
    ResizeArrayIter<PatchElem> ap(patchList); 
    
    // Dynamically accumulates point charges as they are found.
    // The same MM atom may be added to this vector more than once if it is 
    // within the cutoff region of two or more QM regions.
    ResizeArray<ComputeQMPntChrg> compPCVec ;
    
    // This will keep the set of QM groups with which each 
    // point charge should be used. It is re-used for each atom in this
    // patch so we can controll the creation of the compPCVec vector.
    // The first item in the pair is the QM Group ID, the second is the shortest
    // distance between the point charge and an atom of this QM group.
    GrpDistMap chrgGrpMap ;
    
    DebugM(4,"Rank " << CkMyPe() << " receiving from rank " << qmFullMsg->sourceNode
    << " a total of " << qmFullMsg->numAtoms << " QM atoms and " 
    << qmFullMsg->numPCIndxs << " PC IDs." << std::endl);
    
    // If this is the firts step after a step of PC re-selection,
    // we store all PC IDs that all PEs found, in case some atoms end up
    // here untill the next re-selection step.
    if (qmFullMsg->numPCIndxs) {
        
        pcIDSortList.clear();
        
        int *pcIndx = qmFullMsg->pcIndxs;
        for (int i=0; i< qmFullMsg->numPCIndxs;i++) {
            pcIDSortList.load(pcIndx[i]);
        }
        
        pcIDSortList.sort();
    }
    
    int totalNodeAtoms = 0;
    int atomIter = 0;
    int uniqueCounter = 0;
    
    switch ( qmFullMsg->pcSelMode ) {
    
    case PCMODEUPDATESEL:
    {
        
    DebugM(4,"Updating PC selection.\n")
    
    // Loops over all atoms in this node and checks if any MM atom is within 
    // the cutof radius from a QM atom
    for (auto pdIt = patchData.begin(); pdIt != patchData.end(); pdIt++ ) {
        
        CompAtom *x = (*pdIt).compAtomP;
        const CompAtomExt *xExt =  (*pdIt).homePatchP->getCompAtomExtInfo();
        int localNumAtoms =  (*pdIt).homePatchP->getNumAtoms();
        const FullAtom* fullAtms = (*pdIt).homePatchP->getAtomList().const_begin() ;
        
        totalNodeAtoms += localNumAtoms;
        
        // Iterates over the local atoms in a PE, all patches.
        for(int i=0; i<localNumAtoms; ++i) {
            
            const ComputeQMAtom *qmDataPtn = qmFullMsg->coord;
            
            // A new subset for each atom in this node.
            chrgGrpMap.clear();
            
            Position posMM = patchList[0].p->lattice.reverse_transform(x[i].position, 
                                                          fullAtms[i].transform) ;
            
            Real pcGrpID = qmAtomGroup[xExt[i].id];
            
            Real charge = x[i].charge;
            
            // If the atom is a QM atom, there will be no charge info in the 
            // atom box. Therefore, we take the charge in the previous step
            // in case this atom is a point charge for a neighboring QM region.
            if (pcGrpID > 0) {
                for (int qi=0; qi<numQMAtms; qi++) {
                    if (qmAtmIndx[qi] == xExt[i].id) {
                        charge = qmAtmChrg[qi];
                        break;
                    }
                    
                }
            }
            
            for(int j=0; j<qmFullMsg->numAtoms; j++, ++qmDataPtn) {
                
                Real qmGrpID = qmDataPtn->qmGrpID;
                
                // We check If a QM atom and the node atom "i"
                // belong to the same QM group. If not, atom "i" is either an MM
                // atom or a QM atom from another group, which will be seen
                // as an MM point charge.
                // If they belong to the same group, just skip the distance 
                // calculation and move on to the next QM atom.
                // This loop needs to go over all QM atoms since a point charge
                // may me be close to two different QM groups, in which case it 
                // will be sent to both.
                if (qmGrpID == pcGrpID) {
                    continue;
                }
                
                Position posQM = qmDataPtn->position;
                
                Vector distV = posMM - posQM;
                BigReal dist = distV.length();
                
                if ( dist < cutoff ){
                    
                    auto ret = chrgGrpMap.find(qmGrpID) ;
                    
                    // If 'ret' points to end, it means the item was not found
                    // and will be added, which means a new QM region has this 
                    // atom within its cutoff region,
                    // which means a new point charge will be 
                    // created in the QM system in question.
                    // 'ret' means a lot to me!
                    if ( ret ==  chrgGrpMap.end()) {
                        chrgGrpMap.insert(std::pair<Real,BigReal>(qmGrpID, dist));
                    }
                    else {
                        // If the current distance is smaller than the one in 
                        // the map, we update it so that we have the smallest
                        // distance between the point charge and any QM atom in 
                        // this group.
                        if (dist < ret->second) {
                            ret->second = dist;
                        }
                    }
                    
                }
            }
            
            for(auto mapIt = chrgGrpMap.begin(); 
                mapIt != chrgGrpMap.end(); mapIt++) {
                
                // We now add the final info about this point charge 
                // to the vector, repeating it for each QM group that has it
                // within its cuttoff radius.
                compPCVec.add(
                    ComputeQMPntChrg(posMM, charge, xExt[i].id,
                                  mapIt->first, atomIter, mapIt->second, 
                                  fullAtms[i].mass, fullAtms[i].vdwType)
                               );
            }
            
            // Counts how many atoms are seens as point charges, by one or more
            // QM groups.
            if (chrgGrpMap.size() > 0)
                uniqueCounter++;
            
            atomIter++;
        }
        
        (*pdIt).posBoxP->close(&x);
    }
    
    }break; // End case PCMODEUPDATESEL
    
    case PCMODEUPDATEPOS: 
    {
    
    DebugM(4,"Updating PC positions ONLY.\n")
    
    for (auto pdIt = patchData.begin(); pdIt != patchData.end(); pdIt++ ) {
        
        CompAtom *x = (*pdIt).compAtomP;
        const CompAtomExt *xExt =  (*pdIt).homePatchP->getCompAtomExtInfo();
        int localNumAtoms =  (*pdIt).homePatchP->getNumAtoms();
        const FullAtom* fullAtms = (*pdIt).homePatchP->getAtomList().const_begin() ;
        
        totalNodeAtoms += localNumAtoms;
        
        // Iterates over the local atoms in a PE, all patches.
        for(int i=0; i<localNumAtoms; ++i) {
            
            if (pcIDSortList.find(xExt[i].id) != NULL ) {
                Position posMM = patchList[0].p->lattice.reverse_transform(
                                                        x[i].position, 
                                                        fullAtms[i].transform) ;
                
                Real pcGrpID = qmAtomGroup[xExt[i].id];
                Real charge = x[i].charge;
                
                // If the atom is a QM atom, there will be no charge info in the 
                // atom box. Therefore, we take the charge in the previous step
                // in case this atom is a point charge for a neighboring QM region.
                if (pcGrpID > 0) {
                    for (int qi=0; qi<numQMAtms; qi++) {
                        if (qmAtmIndx[qi] == xExt[i].id) {
                            charge = qmAtmChrg[qi];
                            break;
                        }
                        
                    }
                }
                
                compPCVec.add(
                    ComputeQMPntChrg(posMM, charge, xExt[i].id,
                                  0, atomIter, 0, fullAtms[i].mass, 
                                  fullAtms[i].vdwType));
            }
            
            atomIter++;
        }
        
        (*pdIt).posBoxP->close(&x);
    }
    }break ; // End case PCMODEUPDATEPOS
    
    case PCMODECUSTOMSEL:
    {
    
    DebugM(4,"Updating PC positions for custom PC selection.\n")
    
    for (auto pdIt = patchData.begin(); pdIt != patchData.end(); pdIt++ ) {
        
        CompAtom *x = (*pdIt).compAtomP;
        const CompAtomExt *xExt =  (*pdIt).homePatchP->getCompAtomExtInfo();
        int localNumAtoms =  (*pdIt).homePatchP->getNumAtoms();
        const FullAtom* fullAtms = (*pdIt).homePatchP->getAtomList().const_begin() ;
        
        totalNodeAtoms += localNumAtoms;
        
        // Iterates over the local atoms in a PE, all patches.
        for(int i=0; i<localNumAtoms; ++i) {
            
            // Checks if the atom ID belongs to a point charge list,
            // which indicates it is to be sent to rank zero for QM calculations.
            // With one list per QM group, the same point charge can be added to 
            // different QM groups, as if it was being dynamically selected.
            for (int grpIndx=0; grpIndx<numQMGrps; grpIndx++) {
                
                if (customPCLists[grpIndx].find(xExt[i].id) != NULL){
                    
                    Position posMM = patchList[0].p->lattice.reverse_transform(
                                                        x[i].position, 
                                                        fullAtms[i].transform) ;
                    
                    Real pcGrpID = qmAtomGroup[xExt[i].id];
                    Real charge = x[i].charge;
                    
                    // If the atom is a QM atom, there will be no charge info in the 
                    // atom box. Therefore, we take the charge in the previous step
                    // in case this atom is a point charge for a neighboring QM region.
                    if (pcGrpID > 0) {
                        for (int qi=0; qi<numQMAtms; qi++) {
                            if (qmAtmIndx[qi] == xExt[i].id) {
                                charge = qmAtmChrg[qi];
                                break;
                            }
                            
                        }
                    }
                    
                    compPCVec.add(
                        ComputeQMPntChrg(posMM, charge, xExt[i].id,
                                      qmGrpIDArray[grpIndx], atomIter, 
                                      0, fullAtms[i].mass, fullAtms[i].vdwType));
                    
                }
                
            }
            
            atomIter++;
        }
        
        (*pdIt).posBoxP->close(&x);
    }
    
    } break;
    }
    
    DebugM(4,"Rank " << CkMyPe() << " found a total of " << compPCVec.size()
    << " point charges, out of " << totalNodeAtoms
    << " atoms in this node. " << uniqueCounter << " are unique." << std::endl);
    
    // Send only the MM atoms within radius
    
    // this pointer should be freed in rank zero, after receiving it.
    QMPntChrgMsg *pntChrgMsg = new (compPCVec.size(), 0) QMPntChrgMsg;
    pntChrgMsg->sourceNode = CkMyPe();
    pntChrgMsg->numAtoms = compPCVec.size();
    
    for (int i=0; i<compPCVec.size(); i++ ) {
        
        // Only sends the positions and charges of atoms within
        // cutoff of a (any) QM atom. Repeats for the number of QM groups
        // this charge is near to, and indicates the QM group it should 
        // be used with.
        pntChrgMsg->coord[i] = compPCVec[i] ;
        
    }
    
    DebugM(4,"Rank " << pntChrgMsg->sourceNode << " sending a total of " 
    << compPCVec.size() << " elements to rank zero." << std::endl);
    
    CProxy_ComputeQMMgr QMProxy(CkpvAccess(BOCclass_group).computeQMMgr);
    QMProxy[0].recvPntChrg(pntChrgMsg);
    
    delete qmFullMsg;
}



void ComputeQMMgr::procBonds(int numBonds,
                             const int *const qmGrpBondsGrp,
                             const int *const qmMMBondedIndxGrp,
                             const int *const *const chargeTarget,
                             const int *const numTargs,
                             const QMPCVec grpPntChrgVec,
                             QMPCVec &grpAppldChrgVec) {
    
    DebugM(1,"Processing QM-MM bonds in rank zero.\n");
    
    // Indices and IDs of MM atoms which will be passed as point charges.
    std::map<int, int> mmPntChrgMap ;
    
    // Build the map of atom ID and charge for this QM group's point charges.
    // They may be changed by different charge distribution schemes and/or
    // dummy aotm placements.
    for (size_t i=0; i<grpPntChrgVec.size(); i++) {
        
        mmPntChrgMap.insert(std::pair<int,int>(grpPntChrgVec[i].id, (int) i) );
        
        grpAppldChrgVec.push_back( grpPntChrgVec[i] ) ;
        
    }
    
    // If we are treating QM-MM bonds and if there are any
    // in this QM group, we have to treat the MM partial charge
    // using some scheme.
    for (int bondIt = 0; bondIt < numBonds; bondIt++) {
        
        // Gets the global index of this particular QM-MM bond.
        int bondIndx = qmGrpBondsGrp[bondIt] ;
        
        auto retIt = mmPntChrgMap.find(qmMMBondedIndxGrp[bondIt]) ;
        
        // Checks if this MM atom was included as a 
        // point charge due proximity to a QM atom.
        if (retIt == mmPntChrgMap.end()) {
            // If it wasn't, there is an error somwhere.
            
            iout << iERROR << "The MM atom " << qmMMBondedIndxGrp[bondIt]
            << " is bound to a QM atom, but it was not selected as a poitn charge."
            << " Check your cutoff radius!\n" << endi ;
            
            NAMD_die("Charge placement error in QM-MM bond.");
        }
        
        // Gets the (local) index of the MM atom in the point charge vector.
        int mmIndex = (*retIt).second;
        // Gets the position of the MM atom and its partial charge.
        Position mmPos = grpAppldChrgVec[mmIndex].position ;
        BigReal mmCharge = grpAppldChrgVec[mmIndex].charge/numTargs[bondIndx] ;
        
        
        // gives part of the MM charge to neighboring atoms
        for (int i=0; i<numTargs[bondIndx]; i++){
            
            int targetIndxGLobal = chargeTarget[bondIndx][i] ;
            
            retIt = mmPntChrgMap.find(targetIndxGLobal);
            
            // If one of the neighboring atoms which should receive charge
            // is not among partial charges, stop and run away.
            if (retIt == mmPntChrgMap.end()) {
                
                iout << iERROR << "The MM atom " << targetIndxGLobal
                << " is bound to the MM atom of a QM-MM bond and is needed for"
                << " the required bond scheme"
                << " but it was not selected as a poitn charge."
                << " Check your cutoff radius!\n" << endi ;
                
                NAMD_die("Charge placement error in QM-MM bond.");
            }
            
            int trgIndxLocal = (*retIt).second;
            
            // Choose charge treatment method and apply it.
            switch (simParams->qmBondScheme) {
            
                // Charge Schiftig scheme.
                case QMSCHEMECS:
                {
                
//                 grpAppldChrgVec[trgIndxLocal].charge += mmCharge ;
                    
                // Charge Shifting Scheme (DOI: 10.1021/ct100530r)
                // Here we create a dipole to counter act the charge movement
                // we created by moving parts of the MM charge to target MM atoms.
                // 
                // The MM1 charge is equally divided and added to all MM2 atoms.
                // Two point charges are created in the direction of the MM1-MM2 bond,
                // one before and one after the MM2 atom.
                // 
                // Below we see the diagram of the positions of new charges along
                // the direction of the bond between the MM1 atom and a target 
                // MM2 atom, wiht respect to the new Dummy atom (a Hydrogen).
                // 
                //  QM --------  MM1(p0) ------------ MM2
                //  QM ------ H ------------ (+p0) - (MM2 +p0) - (-p0)
                
                Position trgPos = grpAppldChrgVec[trgIndxLocal].position ;
                
                // We create a new point charge at the same position so that
                // the fraction of charge from the MM1 atom is "added" to the
                // MM2 position, but without actually changing the charge of 
                // the original MM2 atom. This trick helps keeping the original 
                // charge for electrostatic calculations after the new charges 
                // of QM atoms is received from the QM software.
                grpAppldChrgVec.push_back(
                    ComputeQMPntChrg(trgPos, mmCharge, -1, 0, -1, 0, 0, 0)
                );
                
                Vector bondVec = trgPos - mmPos ;
                
                Vector bondVec1 = bondVec*0.94 ;
                Vector bondVec2 = bondVec*1.06 ;
                
                Position chrgPos1 = mmPos + bondVec1;
                Position chrgPos2 = mmPos + bondVec2;
                
                BigReal trgChrg1 = mmCharge;
                BigReal trgChrg2 = -1*mmCharge;
                
                
                grpAppldChrgVec.push_back(
                    ComputeQMPntChrg(chrgPos1, trgChrg1, -1, 0, -1, 0, 0, 0)
                );
                
                grpAppldChrgVec.push_back(
                    ComputeQMPntChrg(chrgPos2, trgChrg2, -1, 0, -1, 0, 0, 0)
                );
                
                } break;
                
                // Redistributed Charge and Dipole scheme.
                case QMSCHEMERCD:
                {
                
//                 grpAppldChrgVec[trgIndxLocal].charge -= mmCharge ;
                
                // Redistributed Charge and Dipole method (DOI: 10.1021/ct100530r)
                // Here we create a dipole to counter act the charge movement
                // we created by moving parts of the MM charge to target MM atoms.
                // 
                // The MM1 charge is equally divided and subtracted from all MM2 atoms.
                // One point charge is created in the midpoint of the MM1-MM2 bond.
                // 
                // Below we see the diagram of the positions of new charges along
                // the direction of the bond between the MM1 atom and a target 
                // MM2 atom, wiht respect to the new Dummy atom (a Hydrogen).
                // 
                //  QM --------  MM1(p0) -------------- MM2
                //  QM ------ H ------------ (+2*p0) - (MM2 -p0)
                
                Position trgPos = grpAppldChrgVec[trgIndxLocal].position ;
                
                // We create a new point charge at the same position so that
                // the fraction of charge from the MM1 atom is "added" to the
                // MM2 position, but without actually changing the charge of 
                // the original MM2 atom. This trick helps keeping the original 
                // charge for electrostatic calculations after the new charges 
                // of QM atoms is received from the QM software.
                grpAppldChrgVec.push_back(
                    ComputeQMPntChrg(trgPos, -1*mmCharge, -1, 0, -1, 0, 0, 0)
                );
                
                Vector bondVec = trgPos - mmPos ;
                
                Vector bondVec1 = bondVec*0.5 ;
                
                Position chrgPos1 = mmPos + bondVec1;
                
                BigReal trgChrg1 = 2*mmCharge;
                
                grpAppldChrgVec.push_back(
                    ComputeQMPntChrg(chrgPos1, trgChrg1, -1, 0, -1, 0, 0, 0)
                );
                
                
                
                } break;
                
                // The following schemes simply make the charges of the 
                // MM atom(s) equal to zero. No redistribution or 
                // artificial dipole is attmepted.
                //
                // For Z1, only the MM1 atom has its charge set to zero.
                // For Z2, the MM1 and all MM2 atom charges are changed.
                // For Z3, the MM1, all MM2 and all MM3 atom charges are changed.
                //
                // All modification are local ONLY. They do not change the atom's
                // charge for the classical side of the simulation. Only the QM
                // side ceases to see the electrostatic influence of these atoms,
                // and they cease to see the QM region electrostatics as well.
                
                // Z1, Z2 and Z3 schemes do the same thing. The only difference 
                // is the target list, which is created in MoleculeQM.C.
                case QMSCHEMEZ1:
                
                case QMSCHEMEZ2:
                
                case QMSCHEMEZ3:
                    grpAppldChrgVec[trgIndxLocal].charge = 0.0;
                    break;
            }
            
        }
        
        // We keep this "point charge" so we can calculate forces on it later
        // but it will not influence the QM system.
        // We use the qmGrpID variable to send the message that this poitn charge 
        // should be ignored since this variable will not be relevant anymore,
        // all point charges gathered here are for a specific qmGroup.
        grpAppldChrgVec[mmIndex].qmGrpID = -1 ;
    }
    
    return ;
}


void ComputeQMMgr::recvPntChrg(QMPntChrgMsg *msg) {
    
    // All the preparation that used to be in this function was moved to 
    // recvPartQM, which is called first in rank zero.
    
    if (noPC) {
        // Even if we have QM-MM bonds, the point charge messages 
        // are handled in recvPartQM
        delete msg;
    }
    else {
        pntChrgCoordMsgs[numArrivedPntChrgMsg] = msg;
        ++numArrivedPntChrgMsg;

        if ( numArrivedPntChrgMsg < numSources ) 
            return;
    }
    
    // Resets counters for next step.
    numRecQMRes = 0;
    
    totalEnergy = 0;
    
    for ( int k=0; k<3; ++k )
        for ( int l=0; l<3; ++l )
            totVirial[k][l] = 0;
    
    // ALL DATA ARRIVED --- CALCULATE FORCES
    
    const char *const *const elementArray = molPtr->get_qmElements() ;
    const char *const *const dummyElmArray = molPtr->get_qmDummyElement();
    const int *const qmGrpSize = molPtr->get_qmGrpSizes();
    
    const BigReal *const dummyBondVal = molPtr->get_qmDummyBondVal();
    const int *const grpNumBonds = molPtr->get_qmGrpNumBonds() ;
    const int *const *const qmMMBond = molPtr->get_qmMMBond() ;
    const int *const *const qmGrpBonds = molPtr->get_qmGrpBonds() ;
    const int *const *const qmMMBondedIndx = molPtr->get_qmMMBondedIndx() ;
    const int *const *const chargeTarget = molPtr->get_qmMMChargeTarget() ;
    const int *const numTargs = molPtr->get_qmMMNumTargs() ;
    
    BigReal constants = COULOMB*simParams->nonbondedScaling/(simParams->dielectric*4.0*PI) ;
//     BigReal constants = COULOMB*simParams->nonbondedScaling/(simParams->dielectric) ;
    
    if ( qmPCFreq > 0 ) {
        DebugM(4,"Using point charge stride of " << qmPCFreq << "\n")
        // In case we are skiping steps between PC selection, store a main list
        // in rank zero for future steps. Then we only update positions untill
        // we reach a new "qmPCFreq" step, when a new PC selection is made.
        
        if (timeStep % qmPCFreq == 0) {
            DebugM(4,"Loading a new selection of PCs.\n")
            
            // We only re-set this arrya in a step where a new PC selection is made.
            pntChrgMsgVec.clear();
            for (int pcMsgIt=0; pcMsgIt < numArrivedPntChrgMsg; pcMsgIt++) {
                // Accumulates the message point charges into a local vector.
                for ( int i=0; i < pntChrgCoordMsgs[pcMsgIt]->numAtoms; ++i ) {
                    pntChrgMsgVec.push_back(pntChrgCoordMsgs[pcMsgIt]->coord[i]);
                }
            }
            
            // This fast array is created to send the entire PC IDs list to all
            // PEs in the next step.
            pcIDListSize = pntChrgMsgVec.size();
            pcIDList = new int[pcIDListSize] ;
            for (size_t i=0; i<pcIDListSize; i++) {
                pcIDList[i] = pntChrgMsgVec[i].id;
                
                // Loads home indexes of MM atoms received as point charges.
                // Each force receives the home index of its atom with respect to the 
                // local set of atoms in each node.
                force[pntChrgMsgVec[i].id].homeIndx = pntChrgMsgVec[i].homeIndx;
            }
        }
        else {
            DebugM(4,"Updating position/homeIdex of old PC selection.\n")
            
            // We build a sorted array so that we can quickly access the unordered
            // data we just received, and update positions of the same selection
            // of point charges.
            pcDataSort.clear();
            for (int pcMsgIt=0; pcMsgIt < numArrivedPntChrgMsg; pcMsgIt++) {
                // Accumulates the message point charges into a local sorted array.
                for ( int i=0; i < pntChrgCoordMsgs[pcMsgIt]->numAtoms; ++i ) {
                    pcDataSort.load(pntChrgCoordMsgs[pcMsgIt]->coord[i]);
                }
            }
            pcDataSort.sort();
            
            iout << "Loaded new positions in sorted array: " << pcDataSort.size() << "\n" << endi;
            
            // If we are using a set of point charges that was selected in a
            // previous step, we update the PC POSITION ONLY.
            for (size_t i=0; i<pntChrgMsgVec.size(); i++) {
                
                auto pntr = pcDataSort.find(pntChrgMsgVec[i]);
                
                pntChrgMsgVec[i].position = pntr->position ;
                pntChrgMsgVec[i].homeIndx = pntr->homeIndx ;
                
                // Loads home indexes of MM atoms received as point charges.
                // Each force receives the home index of its atom with respect to the 
                // local set of atoms in each node.
                force[pntChrgMsgVec[i].id].homeIndx = pntChrgMsgVec[i].homeIndx;
            }
        }
    }
    else {
        DebugM(4,"Updating PCs at every step.\n")
        
        pntChrgMsgVec.clear();
        for (int pcMsgIt=0; pcMsgIt < numArrivedPntChrgMsg; pcMsgIt++) {
            // Accumulates the message point charges into a local vector.
            for ( int i=0; i < pntChrgCoordMsgs[pcMsgIt]->numAtoms; ++i ) {
                pntChrgMsgVec.push_back(pntChrgCoordMsgs[pcMsgIt]->coord[i]);
            }
        }
        
        // Loads home indexes of MM atoms received as point charges.
        for (size_t i=0; i<pntChrgMsgVec.size(); i++) {
            // Each force receives the home index of its atom with respect to the 
            // local set of atoms in each node.
            
            #ifdef DEBUG_QM
            if (force[pntChrgMsgVec[i].id].homeIndx != -1 
                and force[pntChrgMsgVec[i].id].homeIndx != pntChrgMsgVec[i].homeIndx
            ) {
                iout << iERROR << "Overloading point charge " 
                << pntChrgMsgVec[i].id << "; home index: " 
                << force[pntChrgMsgVec[i].id].homeIndx << " with " << pntChrgMsgVec[i].homeIndx
                << "\n" << endi ;
                NAMD_die("Error preparing QM simulations.");
            }
            #endif
            
            force[pntChrgMsgVec[i].id].homeIndx = pntChrgMsgVec[i].homeIndx;
        }
    }
    
    // Reset counter for next timestep
    numArrivedPntChrgMsg = 0;
    
    DebugM(4,"A total of " << pntChrgMsgVec.size() 
    << " point charges arrived." << std::endl);
    
    DebugM(4,"Starting QM groups processing.\n");
    
    QMAtmVec grpQMAtmVec;
    QMPCVec grpPntChrgVec;
    
    // Final set of charges, created or not, that is sent to the QM software.
    // This set will depend on how QM-MM bonds are processed and presented to the
    // QM region.
    QMPCVec grpAppldChrgVec;
    
    // Vector of dummy atoms created to treat QM-MM bonds.
    std::vector<dummyData> dummyAtoms ;
    
    // This will hold a big sting with all configuration lines the user supplied.
    std::string configLines ;
    StringList *current = Node::Object()->configList->find("QMConfigLine");
    for ( ; current; current = current->next ) {
        std::string confLineStr(current->data);
        configLines.append(confLineStr);
        configLines.append("\n");
    }
    
    // Initializes the loop for receiving the QM results.
    thisProxy[0].recvQMResLoop() ;
    
    // Iterator for target PE where QM simulations will run.
    int peIter = 0;
    
    for (int grpIter = 0; grpIter < numQMGrps; grpIter++) {
        
        grpQMAtmVec.clear();
        grpPntChrgVec.clear();
        grpAppldChrgVec.clear();
        dummyAtoms.clear();
        
        DebugM(4,"Calculating QM group " << grpIter +1 
        << " (ID: " << grpID[grpIter] << ")." << std::endl);
        
        DebugM(4,"Determining point charges...\n");
        
        Real qmTotalCharge = 0;
        // Loads the QM atoms for this group.
        for (int qmIt=0; qmIt<numQMAtms; qmIt++){
            if (qmCoord[qmIt].qmGrpID == grpID[grpIter]) {
                grpQMAtmVec.push_back(qmCoord[qmIt]);
                qmTotalCharge += qmCoord[qmIt].charge;
            }
        }
        if ((fabsf(roundf(qmTotalCharge) - qmTotalCharge) <= 0.001f) ) {
            qmTotalCharge = roundf(qmTotalCharge) ;
        }
        
        Real pcTotalCharge = 0;
        // Loads the point charges to a local vector for this QM group.
        for (auto pntChrgIt = pntChrgMsgVec.begin(); 
             pntChrgIt != pntChrgMsgVec.end(); pntChrgIt++) {
            if ((*pntChrgIt).qmGrpID == grpID[grpIter] ) {
                grpPntChrgVec.push_back( (*pntChrgIt) );
                pcTotalCharge += (*pntChrgIt).charge;
            }
        }
        if ((fabsf(roundf(pcTotalCharge) - pcTotalCharge) <= 0.001f) ) {
            pcTotalCharge = roundf(pcTotalCharge) ;
        }
        
        #ifdef DEBUG_QM
        if (grpQMAtmVec.size() != qmGrpSize[grpIter]) {
            iout << iERROR << "The number of QM atoms received for group " << grpID[grpIter]
            << " does not match the expected: " << grpQMAtmVec.size()
            << " vs " << qmGrpSize[grpIter] << "\n" << endi ;
            
            NAMD_die("Error processing QM group.");
        }
        #endif
        
        DebugM(1,"Found " << grpPntChrgVec.size() << " point charges for QM group " 
        << grpIter << " (ID: " << grpID[grpIter] << "; Num QM atoms: " 
        << grpQMAtmVec.size() <<  "; Num QM-MM bonds: " 
        << grpNumBonds[grpIter] << ")" << std::endl);
        
        DebugM(1,"Total QM charge: " << qmTotalCharge 
        << "; Total point-charge charge: " << pcTotalCharge << std::endl);
        
        // If we have a frequency for LSS update, check if we shoudl do it in 
        // the current time step.
        if ( qmLSSFreq > 0 && ((timeStep + 1) % qmLSSFreq == 0 )) {
            lssUpdate(grpIter, grpQMAtmVec, grpPntChrgVec);
        }
        
        // This function checks data and treats the charge (and existence) of
        // the MM atoms in and around QM-MM bonds. It is only executed in 
        // electrostatic embeding QM/MM simulations.
        if (! noPC ) {
            procBonds(grpNumBonds[grpIter], qmGrpBonds[grpIter], 
                     qmMMBondedIndx[grpIter], 
                     chargeTarget, numTargs, 
                     grpPntChrgVec, grpAppldChrgVec) ;
        }
        else {
            grpAppldChrgVec = grpPntChrgVec;
        }
        
        // For future use, we get the pairs of indexes of QM atoms and MM atoms which are
        // bound in QM-MM bonds.
        std::vector< std::pair<int,int> > qmPCLocalIndxPairs ;
        
        // Create and position dummy atoms.
        Position mmPos(0), qmPos(0);
        for (int dummyIt = 0; dummyIt < grpNumBonds[grpIter]; dummyIt++) {
            
            int qmMMBondIndx = qmGrpBonds[grpIter][dummyIt] ;
            
            BigReal bondVal = dummyBondVal[qmMMBondIndx] ;
            
            int mmAtmIndx = qmMMBond[qmMMBondIndx][0] ;
            int qmAtmIndx = qmMMBond[qmMMBondIndx][1] ;
            
            // Sicne we don't know in which patch/node the QM atom is, or the 
            // order in which they will arrive in rank zero, we have
            // no direct index to it.
            #ifdef DEBUG_QM
            bool missingQM = true, missingMM = true;
            #endif
            size_t qmIt ;
            for (qmIt=0; qmIt<grpQMAtmVec.size(); qmIt++){
                if (grpQMAtmVec[qmIt].id == qmAtmIndx) {
                    qmPos = grpQMAtmVec[qmIt].position;
                    
                    #ifdef DEBUG_QM
                    missingQM = false;
                    #endif
                    
                    break;
                }
            }
            // The same is true about the MM atom to which the QM atom is bound,
            // we must look
            size_t pcIt;
            for (pcIt=0; pcIt < grpPntChrgVec.size(); pcIt++) {
                if (grpPntChrgVec[pcIt].id == mmAtmIndx) {
                    mmPos = grpPntChrgVec[pcIt].position ;
                    
                    #ifdef DEBUG_QM
                    missingMM = false;
                    #endif
                    
                    break;
                }
            }
            
            qmPCLocalIndxPairs.push_back(std::pair<int,int>(qmIt, pcIt) );
            
            #ifdef DEBUG_QM
            // Checks if the MM atom was included as a 
            // point charge due proximity to a QM atom, and if the QM atom arrived.
            if ( missingMM or missingQM ) {
                // If it wasn't, there is an error somwhere.
                
                if (missingMM) {
                    iout << iERROR << "The MM atom " << mmAtmIndx
                    << " is bound to a QM atom, but it was not selected as a poitn charge."
                    << " Check your cutoff radius!\n" << endi ;
                    
                    NAMD_die("Error in QM-MM bond processing.");
                }
                if (missingQM) {
                    iout << iERROR << "The QM atom " << qmAtmIndx
                    << " is bound to an MM atom, but it was not sent to rank zero for processing."
                    << " Check your configuration!\n" << endi ;
                    
                    NAMD_die("Error in QM-MM bond processing.");
                }
            }
            #endif
            
            Vector bondVec = mmPos - qmPos ;
            
            if (bondValType == QMLENTYPE) {
                // If a length is defined by the user, or a default len
                // is used, we calculate the unit vector for the displacement
                // and multiply by the desired length in order 
                // to get the final dummy atom position relative to the
                // QM atom.
                bondVec = bondVec.unit() ;
                bondVec *= bondVal ;
            }
            else if (bondValType == QMRATIOTYPE) {
                // If distance a ratio was defined by the user, then
                // the displacement vector is multiplied by that ratio
                // to get the final dummy atom position relative to the
                // QM atom.
                bondVec *= bondVal ;
            }
            
            Position dummyPos = qmPos + bondVec;
            
            DebugM(1,"Creating dummy atom " << dummyPos << " ; QM ind: " 
            << qmIt << " ; PC ind: " << pcIt << std::endl);
            
            dummyAtoms.push_back(dummyData(dummyPos, qmIt, qmMMBondIndx)) ;
            
        }
        
        DebugM(3, "Creating data for " << grpQMAtmVec.size() << " QM atoms " 
        << dummyAtoms.size() << " dummy atoms " << grpPntChrgVec.size()
        << " real point charges and " << grpAppldChrgVec.size() - grpPntChrgVec.size()
        << " virtual point charges\n") ;
        
        int dataSize = grpQMAtmVec.size() + dummyAtoms.size() + grpAppldChrgVec.size();
        QMGrpCalcMsg *msg = new (dataSize, configLines.size(), 0) QMGrpCalcMsg;
        msg->grpIndx = grpIter;
        msg->peIter = peIter;
        msg->charge = grpChrg[grpIter];
        msg->multiplicity = grpMult[grpIter];
        msg->numQMAtoms = grpQMAtmVec.size();
        msg->numAllAtoms = grpQMAtmVec.size() + dummyAtoms.size();
        msg->numRealPntChrgs = grpPntChrgVec.size(); // The original set of classical atoms.
        msg->numAllPntChrgs = grpAppldChrgVec.size(); // The extended set with virtual point charges.
        msg->secProcOn = simParams->qmSecProcOn ;
        msg->constants = constants;
        msg->PMEOn = simParams->PMEOn ;
        if (msg->PMEOn)
            msg->PMEEwaldCoefficient = simParams->PMEEwaldCoefficient ;
        msg->switching = simParams->qmPCSwitchOn;
        msg->switchType = simParams->qmPCSwitchType;
        msg->cutoff = simParams->cutoff;
        msg->swdist = simParams->switchingDist;
        msg->pcScheme = simParams->qmPCScheme;
        msg->qmAtmChrgMode = simParams->qmChrgMode;
        
        strncpy(msg->baseDir, simParams->qmBaseDir, 256);
        strncpy(msg->execPath, simParams->qmExecPath, 256);
        if (msg->secProcOn)
            strncpy(msg->secProc, simParams->qmSecProc, 256);
        
        if (simParams->qmPrepProcOn && (timeStep == simParams->firstTimestep)) {
            msg->prepProcOn = true;
            strncpy(msg->prepProc, simParams->qmPrepProc, 256);
        } else
            msg->prepProcOn = false;
        
        QMAtomData *dataP = msg->data;
        
        for (int i=0; i<grpQMAtmVec.size(); i++) {
            dataP->position = grpQMAtmVec[i].position ;
            dataP->charge = grpQMAtmVec[i].charge ;
            dataP->id = grpQMAtmVec[i].id ;
            dataP->bountToIndx = -1;
            dataP->type = QMATOMTYPE_QM ;
            strncpy(dataP->element,elementArray[grpQMAtmVec[i].id],3);
            dataP++;
        }
        
        for (int i=0; i<dummyAtoms.size(); i++) {
            dataP->position = dummyAtoms[i].pos ;
            dataP->charge = 0 ;
            dataP->id = -1 ;
            dataP->bountToIndx = dummyAtoms[i].qmInd ;
            dataP->type = QMATOMTYPE_DUMMY ;
            strncpy(dataP->element,dummyElmArray[dummyAtoms[i].bondIndx],3);
            dataP++;
        }
        
        for (int i=0; i<grpAppldChrgVec.size(); i++) {
            dataP->position = grpAppldChrgVec[i].position ;
            dataP->charge = grpAppldChrgVec[i].charge ;
            // Point charges created to handle QM-MM bonds will have an id of -1.
            dataP->id = grpAppldChrgVec[i].id ;
            dataP->bountToIndx = -1;
            dataP->dist = grpAppldChrgVec[i].dist ;
            // If we are loading the classical atoms' charges
            // the point charge type is 1, unless it is from an 
            // atom which is bound to a QM atom.
            if (i < grpPntChrgVec.size()) {
                if (grpAppldChrgVec[i].qmGrpID < 0) {
                    dataP->type = QMPCTYPE_IGNORE ;
                }
                else {
                    dataP->type = QMPCTYPE_CLASSICAL ;
                }
            }
            else {
                // Extra charges are created to handle QM-MM bonds (if they exist).
                dataP->type = QMPCTYPE_EXTRA ;
            }
            dataP++;
        }
        
        QMAtomData *qmP = msg->data ;
        QMAtomData *pcP = msg->data + msg->numAllAtoms ;
        
        // With this, every QM atom knows to which MM atom is is bound, 
        // and vice versa. This will be usefull later on to prevent them from
        // feeling eachother's electrostatic charges AND to place the dummy
        // atom forces on the "real" atoms that form the bond.
        for( auto vecPtr  = qmPCLocalIndxPairs.begin(); 
                  vecPtr != qmPCLocalIndxPairs.end(); 
                  vecPtr++ ) {
            
            int qmLocInd = (*vecPtr).first;
            int pcLocInd = (*vecPtr).second;
            
            qmP[qmLocInd].bountToIndx = pcLocInd ;
            pcP[pcLocInd].bountToIndx = qmLocInd ;
        }
        
        
        strcpy(msg->configLines, configLines.c_str());
        
        int targetPE = qmPEs[peIter] ;
        
        DebugM(4,"Sending QM group " << grpIter << " (ID " << grpID[grpIter] 
        << ") to PE " << targetPE << std::endl);
        
        switch (simParams->qmFormat) {
            // Creates the command line that will be executed according to the 
            // chosen QM software, as well as the input file with coordinates.
            case QMFormatORCA:
                QMProxy[targetPE].calcORCA(msg) ;
                break;
            
            case QMFormatMOPAC:
                QMProxy[targetPE].calcMOPAC(msg) ;
                break;
            
            case QMFormatUSR:
                QMProxy[targetPE].calcUSR(msg) ;
                break;
        }
        
        peIter++;
        
        if (peIter == qmPEs.size())
            peIter = 0;
    }
    
}

void ComputeQMMgr::storeQMRes(QMGrpResMsg *resMsg) {
    
    iout << iINFO << "Storing QM results for region " << resMsg->grpIndx  
    << " (ID: "  << grpID[resMsg->grpIndx] 
    << ") with original energy: " << endi;
    std::cout << std::fixed << std::setprecision(6) << resMsg->energyOrig << endi;
    iout << " Kcal/mol\n" << endi;
    
    if (resMsg->energyCorr != resMsg->energyOrig) {
        iout << iINFO << "PME corrected energy: " << endi;
        std::cout << std::fixed << std::setprecision(6) << resMsg->energyCorr << endi;
        iout << " Kcal/mol\n" << endi;
    }
    
    totalEnergy += resMsg->energyCorr ;
    
    for ( int k=0; k<3; ++k )
        for ( int l=0; l<3; ++l )
            totVirial[k][l] += resMsg->virial[k][l];
    
    QMForce *fres = resMsg->force ;
    Real qmTotalCharge = 0;
    
    for (int i=0; i<resMsg->numForces; i++) {
        
        force[ fres[i].id ].force += fres[i].force;
        
        // Indicates the result is a QM atom, and we should get its updated charge.
        if (fres[i].replace == 1) {
            force[ fres[i].id ].charge =  fres[i].charge;
            qmTotalCharge += fres[i].charge;
        }
    }
    
    if ((fabsf(roundf(qmTotalCharge) - qmTotalCharge) <= 0.001f) ) {
        qmTotalCharge = roundf(qmTotalCharge) ;
    }
    
    DebugM(4,"QM total charge received is " << qmTotalCharge << std::endl);
    
    DebugM(4,"Current accumulated energy is " << totalEnergy << std::endl);
    
    numRecQMRes++;
    
    delete resMsg;
}

void ComputeQMMgr::procQMRes() {
    
    // Writes a DCD file with the charges of all QM atoms at a frequency 
    // defined by the user in qmOutFreq.
    if ( simParams->qmOutFreq > 0 && 
         timeStep % simParams->qmOutFreq == 0 ) {
        
        iout << iINFO << "Writing QM charge output at step " 
        << timeStep <<  "\n" << endi;
    
        Real *x = outputData, 
             *y = outputData + molPtr->get_numQMAtoms(), 
             *z = outputData + 2*molPtr->get_numQMAtoms();
        
        for (int qmIt=0; qmIt<numQMAtms; qmIt++){
            x[qmIt] = qmCoord[qmIt].id;
            y[qmIt] = force[qmCoord[qmIt].id].charge ;
            z[qmIt] = 0;
        }
        
        write_dcdstep(dcdOutFile, numQMAtms, x, y, z, 0) ;
    }
    
    // Writes a DCD file with the charges of all QM atoms at a frequency 
    // defined by the user in qmPosOutFreq.
    if ( simParams->qmPosOutFreq > 0 && 
         timeStep % simParams->qmPosOutFreq == 0 ) {
        
        iout << iINFO << "Writing QM position output at step " 
        << timeStep <<  "\n" << endi;
        
        SortedArray<idIndxStr> idIndx;
        
        for(int i=0; i<numQMAtms;i++) {
            idIndx.insert( idIndxStr(qmCoord[i].id, i) );
        }
        idIndx.sort();
        
        Real *x = outputData, 
             *y = outputData + molPtr->get_numQMAtoms(), 
             *z = outputData + 2*molPtr->get_numQMAtoms();
        
        for (int qmIt=0; qmIt<numQMAtms; qmIt++){
            x[qmIt] = qmCoord[idIndx[qmIt].indx].position.x;
            y[qmIt] = qmCoord[idIndx[qmIt].indx].position.y;
            z[qmIt] = qmCoord[idIndx[qmIt].indx].position.z;
        }
        
        write_dcdstep(dcdPosOutFile, numQMAtms, x, y, z, 0) ;
    }
    
    // distribute forces
    DebugM(4,"Distributing QM forces for all ranks.\n");
    for ( int j=0; j < numSources; ++j ) {
        
        DebugM(1,"Sending forces and charges to source " << j << std::endl);
        
        QMCoordMsg *qmmsg = 0;
        QMPntChrgMsg *pcmsg = 0 ;
        
        int totForces = 0;
        int sourceNode = -1;
        
        if (pntChrgCoordMsgs == NULL) {
            
            qmmsg = qmCoordMsgs[j];
            qmCoordMsgs[j] = 0;
            
            totForces = qmmsg->numAtoms ;
            
            sourceNode = qmmsg->sourceNode;
        }
        else {
            pcmsg = pntChrgCoordMsgs[j];
            pntChrgCoordMsgs[j] = 0;
            
            sourceNode = pcmsg->sourceNode;
            
            // Since we receive two different messages from nodes, there is no 
            // guarantee the two sets of messages will come in the same order.
            // Therefore, we match the messages by comaring their sourceNodes.
            for (int aux=0; aux<numSources; aux++) {
                
                if (qmCoordMsgs[aux] == 0)
                    continue;
                
                qmmsg = qmCoordMsgs[aux];
                
                if (qmmsg->sourceNode == sourceNode) {
                    qmCoordMsgs[aux] = 0;
                    break;
                }
            }
            
            DebugM(1,"Building force mesage for rank " 
            << pcmsg->sourceNode << std::endl);
            
            totForces = qmmsg->numAtoms + pcmsg->numAtoms;
        }
        
        QMForceMsg *fmsg = new (totForces, subsArray.size(), 0) QMForceMsg;
        fmsg->PMEOn = simParams->PMEOn;
        fmsg->numForces = totForces;
        fmsg->numLssDat = subsArray.size();
        
        DebugM(1,"Loading QM forces.\n");
        
        // This iterator is used in BOTH loops!
        int forceIter = 0;
        
        for ( int i=0; i < qmmsg->numAtoms; ++i ) {
            fmsg->force[forceIter] = force[qmmsg->coord[i].id];
            forceIter++;
        }
        
        delete qmmsg;
        
        if (pntChrgCoordMsgs != NULL) {
            DebugM(1,"Loading PntChrg forces.\n");
            
            for ( int i=0; i < pcmsg->numAtoms; ++i ) {
                fmsg->force[forceIter] = force[pcmsg->coord[i].id];
                forceIter++;
            }
            
            delete pcmsg;
        }
        
        DebugM(1,"A total of " << forceIter << " forces were loaded." << std::endl);
        
        for ( int i=0; i < subsArray.size(); ++i ) {
            fmsg->lssDat[i] = subsArray[i];
        }
        
        #ifdef DEBUG_QM
        if (subsArray.size() > 0)
            DebugM(3,"A total of " << subsArray.size() << " LSS data structures were loaded." << std::endl);
        #endif
        
        if ( ! j ) {
            fmsg->energy = totalEnergy;
            for ( int k=0; k<3; ++k )
                for ( int l=0; l<3; ++l )
                    fmsg->virial[k][l] = totVirial[k][l];
        } else {
            fmsg->energy = 0;
            for ( int k=0; k<3; ++k )
                for ( int l=0; l<3; ++l )
                    fmsg->virial[k][l] = 0;
        }
        
        DebugM(4,"Sending forces...\n");
        
        QMProxy[sourceNode].recvForce(fmsg);
        
    }
    
    DebugM(4,"All forces sent from node zero.\n");
}

void ComputeQMMgr::recvForce(QMForceMsg *fmsg) {
    
    if (CkMyPe()) {
        for (int i=0; i<fmsg->numLssDat; i++) {
            subsArray.add(fmsg->lssDat[i]) ;
        }
    }
    
    QMCompute->saveResults(fmsg);
}

void ComputeQM::saveResults(QMForceMsg *fmsg) {
    
    ResizeArrayIter<PatchElem> ap(patchList);
    
    bool callReplaceForces = false;
    
    int numQMAtms = Node::Object()->molecule->get_numQMAtoms();
    const Real * const qmAtomGroup = Node::Object()->molecule->get_qmAtomGroup() ;
    const int *qmAtmIndx = Node::Object()->molecule->get_qmAtmIndx() ;
    Real *qmAtmChrg = Node::Object()->molecule->get_qmAtmChrg() ;
    
    int totalNumbAtoms = 0;
    for (ap = ap.begin(); ap != ap.end(); ap++) {
        totalNumbAtoms += (*ap).p->getNumAtoms();
    }
    
    // This is kept to be deleted in the next step so that the message can be 
    // used in "replaceForces" routine later on, in case there are forces to 
    // be replaced. The "replaceForces" function uses but DOES NOT free the pointer,
    // so we free the data from the previous iteration and allocate a new one for
    // the current iteration (remember that the number of atoms can change in a 
    // home patch between iterations).
    if (oldForces != 0)
        delete [] oldForces;
    oldForces = new ExtForce[totalNumbAtoms] ;
    
    for (int i=0; i < totalNumbAtoms; ++i) {
        oldForces[i].force = Force(0) ;
    }
    
    DebugM(1,"Force array has been created and zeroed in rank " 
    << CkMyPe() << std::endl);
    
    DebugM(1,"Preparing " << fmsg->numForces << " forces in rank " 
    << CkMyPe() << std::endl);
    
    QMForce *results_ptr = fmsg->force;
    // Iterates over the received forces and place them on the full array.
    for (int i=0; i < fmsg->numForces; ++i, ++results_ptr) {
        // For now we may have more than one item in the message acting on the same 
        // atom, such as an MM atom which is a point charge for two or more QM regions.
        
        oldForces[results_ptr->homeIndx].force += results_ptr->force;
        oldForces[results_ptr->homeIndx].replace = results_ptr->replace;
        
        if (results_ptr->replace == 1)
            callReplaceForces = true;
        
        // If the atom is in a QM group, update its charge to the local (this homePatch)
        // copy of the qmAtmChrg array.
        if (qmAtomGroup[results_ptr->id] > 0 && fmsg->PMEOn) {
            
            // Loops over all QM atoms (in all QM groups) comparing their global indices
            for (int qmIter=0; qmIter<numQMAtms; qmIter++) {
                
                if (qmAtmIndx[qmIter] == results_ptr->id) {
                    qmAtmChrg[qmIter] = results_ptr->charge;
                    break;
                }
                
            }
            
        }
        
    }
    
    DebugM(1,"Placing forces on force boxes in rank " 
    << CkMyPe() << std::endl);
    
    // Places the received forces on the force array for each patch.
    int homeIndxIter = 0;
    for (ap = ap.begin(); ap != ap.end(); ap++) {
        Results *r = (*ap).forceBox->open();
        Force *f = r->f[Results::normal];
        const CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
        int localNumAtoms = (*ap).p->getNumAtoms();
        
        for(int i=0; i<localNumAtoms; ++i) {
            
            f[i] += oldForces[homeIndxIter].force;
            
            ++homeIndxIter;
        }
        
        if ( callReplaceForces )
            (*ap).p->replaceForces(oldForces);
        
        (*ap).forceBox->close(&r);
        
    }
    
    DebugM(1,"Forces placed on force boxes in rank " 
    << CkMyPe() << std::endl);
    
    if (fmsg->PMEOn) {
        
        DebugM(1,"PME ON! Accessing PMEmgr in rank " << CkMyPe() << std::endl);
        
        ComputePmeMgr *mgrP = CProxy_ComputePmeMgr::ckLocalBranch(
            CkpvAccess(BOCclass_group).computePmeMgr) ;
        
        DebugM(4, "Initiating ComputePme::doQMWork on rank " << CkMyPe() << " over " 
            << getComputes(mgrP).size() << " pmeComputes." << std::endl) ;
        
        for ( int i=0; i< getComputes(mgrP).size(); ++i ) {
    //         WorkDistrib::messageEnqueueWork(pmeComputes[i]);
            getComputes(mgrP)[i]->doQMWork();
        }
    }
    
    DebugM(1,"Submitting reduction in rank " << CkMyPe() << std::endl);
    
    reduction->item(REDUCTION_MISC_ENERGY) += fmsg->energy;
    reduction->item(REDUCTION_VIRIAL_NORMAL_XX) += fmsg->virial[0][0];
    reduction->item(REDUCTION_VIRIAL_NORMAL_XY) += fmsg->virial[0][1];
    reduction->item(REDUCTION_VIRIAL_NORMAL_XZ) += fmsg->virial[0][2];
    reduction->item(REDUCTION_VIRIAL_NORMAL_YX) += fmsg->virial[1][0];
    reduction->item(REDUCTION_VIRIAL_NORMAL_YY) += fmsg->virial[1][1];
    reduction->item(REDUCTION_VIRIAL_NORMAL_YZ) += fmsg->virial[1][2];
    reduction->item(REDUCTION_VIRIAL_NORMAL_ZX) += fmsg->virial[2][0];
    reduction->item(REDUCTION_VIRIAL_NORMAL_ZY) += fmsg->virial[2][1];
    reduction->item(REDUCTION_VIRIAL_NORMAL_ZZ) += fmsg->virial[2][2];
    reduction->submit();
    
    delete fmsg ;
}

void ComputeQMMgr::calcMOPAC(QMGrpCalcMsg *msg)
{
    
    FILE *inputFile,*outputFile,*chrgFile;
    int iret;
    
    const size_t lineLen = 256;
    char *line = new char[lineLen];
    
    std::string qmCommand, inputFileName, outputFileName, pntChrgFileName;
    
    // For coulomb calculation
    BigReal constants = msg->constants;
    
    double gradient[3];
    
    DebugM(4,"Running MOPAC on PE " << CkMyPe() << std::endl);
    
    if (msg->switching)
        pntChrgSwitching(msg) ;
    
    // For each QM group, create a subdirectory where files will be palced.
    std::string baseDir(msg->baseDir);
    baseDir.append("/") ;
    std::ostringstream itosConv ;
    itosConv << msg->peIter ;
    baseDir += itosConv.str() ;
    
    struct stat info;
    
    if (stat(msg->baseDir, &info) != 0 ) {
        CkPrintf( "Node %d cannot access directory %s\n",
                  CkMyPe(), baseDir.c_str() );
        NAMD_die("QM calculation could not be ran. Check your qmBaseDir!");
    }
    else if (! (stat(baseDir.c_str(), &info) == 0 && S_ISDIR(info.st_mode)) ) {
        DebugM(4,"Creating directory " << baseDir.c_str() << std::endl);
        int retVal = mkdir(baseDir.c_str(), S_IRWXU);
    }
    
    #ifdef DEBUG_QM
    Write_PDB(std::string(baseDir)+"/input.pdb", msg ) ;
    #endif
    
    inputFileName.clear();
    inputFileName.append(baseDir.c_str()) ;
    inputFileName.append("/qmmm_") ;
    inputFileName += itosConv.str() ;
    inputFileName.append(".input") ;
    
    // Opens file for coordinate and parameter input
    inputFile = fopen(inputFileName.c_str(),"w");
    if ( ! inputFile ) {
        iout << iERROR << "Could not open input file for writing: " 
        << inputFileName << "\n" << endi ;
        NAMD_die(strerror(errno));
    }
    
    // Builds the command that will be executed
    qmCommand.clear();
    qmCommand.append("cd ");
    qmCommand.append(baseDir);
    qmCommand.append(" ; ");
    qmCommand.append(msg->execPath) ;
    qmCommand.append(" ") ;
    qmCommand.append(inputFileName) ;
    
    // Builds the file name where MOPAC will place the gradient
    // This is also relative to the input file name.
    outputFileName = inputFileName ;
    outputFileName.append(".aux") ;
    
    if (msg->numAllPntChrgs) {
        // Builds the file name where we will write the point charges.
        pntChrgFileName.clear();
        pntChrgFileName.append(baseDir.c_str()) ;
        pntChrgFileName.append("/mol.in") ;
        
        chrgFile = fopen(pntChrgFileName.c_str(),"w");
        if ( ! chrgFile ) {
            iout << iERROR << "Could not open charge file for writing: " 
            << pntChrgFileName << "\n" << endi ;
            NAMD_die(strerror(errno));
        }
        
        iret = fprintf(chrgFile,"\n%d %d\n", 
                       msg->numQMAtoms, msg->numAllAtoms - msg->numQMAtoms);
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    }
    
    // writes configuration lines to the input file.
    std::stringstream ss(msg->configLines) ;
    std::string confLineStr;
    while (std::getline(ss, confLineStr) ) {
        confLineStr.append("\n");
        iret = fprintf(inputFile,confLineStr.c_str());
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    }
    
    
    iret = fprintf(inputFile,"\n");
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    DebugM(4, "Writing " << msg->numAllAtoms << " QM atom coords in file " 
    << inputFileName.c_str() << " and " << msg->numAllPntChrgs 
    << " point charges in file " << pntChrgFileName.c_str() << "\n");
    
    // write QM and dummy atom coordinates to input file and
    // MM electric field from MM point charges.
    QMAtomData *atmP = msg->data ;
    QMAtomData *pcP = msg->data + msg->numAllAtoms ;
    for (size_t i=0; i<msg->numAllAtoms; ++i, ++atmP ) {
        
        double x = atmP->position.x;
        double y = atmP->position.y;
        double z = atmP->position.z;
        
        iret = fprintf(inputFile,"%s %f 1 %f 1 %f 1\n",
                       atmP->element,x,y,z);
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
        
        if (msg->numAllPntChrgs) {
            BigReal phi = 0;
            
            pcP = msg->data + msg->numAllAtoms ;
            for ( size_t j=0; j < msg->numAllPntChrgs; ++j, ++pcP ) {
                
                if (pcP->type == QMPCTYPE_IGNORE)
                    continue;
                
                double charge = pcP->charge;
                
                double xMM = pcP->position.x;
                double yMM = pcP->position.y;
                double zMM = pcP->position.z;
                
                BigReal rij = sqrt((x-xMM)*(x-xMM) +
                     (y-yMM)*(y-yMM) +
                     (z-zMM)*(z-zMM) ) ;
                     
                phi += charge/rij ;
            }
            
            phi = phi*constants ;
            
            iret = fprintf(chrgFile,"%s %f %f %f %f\n",
                           atmP->element,x,y,z, phi);
            if ( iret < 0 ) { NAMD_die(strerror(errno)); }
        }
    }
    
    DebugM(4,"Closing input and charge file\n");
    
    if (msg->numAllPntChrgs)
        fclose(chrgFile);
    
    fclose(inputFile);
    
    if (msg->prepProcOn) {
        
        std::string prepProc(msg->prepProc) ;
        prepProc.append(" ") ;
        prepProc.append(inputFileName) ;
        iret = system(prepProc.c_str());
        if ( iret == -1 ) { NAMD_die(strerror(errno)); }
        if ( iret ) { NAMD_die("Error running preparation command for QM calculation."); }
    }
    
    // runs QM command
    DebugM(4,"Running command ->" << qmCommand.c_str() << "<-" << std::endl);
    iret = system(qmCommand.c_str());
    
    if ( iret == -1 ) { NAMD_die(strerror(errno)); }
    if ( iret ) { NAMD_die("Error running command for QM forces calculation."); }

    if (msg->secProcOn) {
        
        std::string secProc(msg->secProc) ;
        secProc.append(" ") ;
        secProc.append(inputFileName) ;
        iret = system(secProc.c_str());
        if ( iret == -1 ) { NAMD_die(strerror(errno)); }
        if ( iret ) { NAMD_die("Error running second command for QM calculation."); }
    }
    
    // remove coordinate file
//     iret = remove(inputFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }

    // remove coordinate file
//     iret = remove(pntChrgFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }
    
    // opens output file
    DebugM(4,"Reading QM data from file " << outputFileName.c_str() << std::endl);
    outputFile = fopen(outputFileName.c_str(),"r");
    if ( ! outputFile ) {
        iout << iERROR << "Could not find QM output file!\n" << endi;
        NAMD_die(strerror(errno)); 
    }
    
    // Resets the pointers.
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // Allocates the return message.
    QMGrpResMsg *resMsg = new (msg->numQMAtoms + msg->numRealPntChrgs, 0) QMGrpResMsg;
    resMsg->grpIndx = msg->grpIndx;
    resMsg->numForces = msg->numQMAtoms + msg->numRealPntChrgs;
    resMsg->energyOrig = 0;
    resMsg->energyCorr = 0;
    for ( int k=0; k<3; ++k )
        for ( int l=0; l<3; ++l )
            resMsg->virial[k][l] = 0;
    QMForce *resForce = resMsg->force ;
    
    // We split the data into two pointers so we can "skip" the dummy atoms
    // which have no representation in NAMD.
    for (int i=0; i<resMsg->numForces; i++) {
        resForce[i].force = 0;
        resForce[i].charge = 0 ;
        if (i < msg->numQMAtoms) {
            // We use the replace field to indicate QM atoms,
            // which will have charge information.
            resForce[i].replace = 1 ;
            resForce[i].id = atmP->id;
            atmP++;
        }
        else {
            // We use the replace field to indicate QM atoms,
            // which will have charge information.
            resForce[i].replace = 0 ;
            resForce[i].id = pcP->id;
            pcP++;
        }
    }
    
    // Resets the pointers.
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // Reads the data form the output file created by the QM software.
    // Gradients over the QM atoms, and Charges for QM atoms will be read.
    
    size_t atmIndx = 0, gradCount = 0;
    Bool gradFields = false, chargeFields = false;
    Bool chargesRead = false, gradsRead = false;
    while ( fgets(line, lineLen, outputFile) != NULL){
        // We loop over the lines of the output file untill we find
        // the line that initiates the "atom charges" lines. Then
        // we read all lines and expect to get one or more charges 
        // per line, separated by space(s), untill we find a line that
        // initiates the "gradients", then once more, we expect several
        // numbers separated by space(s). When the "overlap matrix" 
        // string is found, we break the loop and stop reading the file.
        
        
        if ( strstr(line,"TOTAL_ENERGY") != NULL ) {
            
            char strEnergy[14], *endPtr ;
            
            strncpy(strEnergy, line + 17, 13) ;
            strEnergy[13] = '\0';
            
            // We have to convert the notation from FORTRAN double precision to
            // the natural, normal, reasonable and not terribly out dated format.
            resMsg->energyOrig = strtod(strEnergy, &endPtr);
            if ( *endPtr == 'D' ) {
                *endPtr = 'e';
                resMsg->energyOrig = strtod(strEnergy, &endPtr);
            }
            
            // In MOPAC, the total energy is given in EV, so we convert to Kcal/mol
            resMsg->energyOrig *= 23.060 ;
            
//             DebugM(4,"Reading QM energy from file: " << resMsg->energyOrig << "\n");
            
            resMsg->energyCorr = resMsg->energyOrig;
            
            continue;
        }
        
        if ( strstr(line,"ATOM_CHARGES") != NULL ) {
            gradFields = false;
            chargeFields = true;
            atmIndx = 0;
            
            // If the used ask for charges NOT to be read from MOPAC,
            // we skip the charge reading step.
            if (msg->qmAtmChrgMode == QMCHRGNONE) {
                chargeFields = false;
                atmIndx = msg->numAllAtoms;
            }
            else {
                chargeFields = true;
                atmIndx = 0;
            }
            
            // Now we expect the following line(s) to have atom charges
            continue;
        }
        
        if ( strstr(line,"GRADIENTS") != NULL ) {
            
            // Now that all charges have been read, checks if the 
            // numbers match
            if (atmIndx != msg->numAllAtoms) {
                NAMD_die("Error reading QM forces file. Wrong number of atom charges");
            }
            
            chargesRead = true;
            
            // Now we expect the following line(s) to have gradients
            chargeFields = false ;
            gradFields = true;
            gradCount = 0;
            atmIndx = 0;
            
            continue;
        }
        
        if ( strstr(line,"OVERLAP_MATRIX") != NULL ) {
            
            // Now that all gradients have been read, checks if the 
            // numbers match
            if (atmIndx != msg->numAllAtoms) {
                NAMD_die("Error reading QM forces file. Wrong number of gradients");
            }
            
            gradsRead = true;
            
            // Nothing more to read from the ".aux" file
            break;
        }
        
        char result[10] ;
        size_t strIndx = 0;
        
        if (chargeFields) {
            while ((strIndx < (strlen(line)-9)) && (strlen(line)-1 >=9 ) ) {
                
                strncpy(result, line+strIndx,9) ;
                result[9] = '\0';
                
                Real localCharge = atof(result);
                
                // If we are reading charges from QM atoms, store them.
                if (atmIndx < msg->numQMAtoms ) {
                    atmP[atmIndx].charge = localCharge;
                    resForce[atmIndx].charge = localCharge;
                }
                
                // If we are reading charges from Dummy atoms,
                // place them on the appropriate QM atom.
                if ( msg->numQMAtoms <= atmIndx &&
                    atmIndx < msg->numAllAtoms ) {
                    // The dummy atom points to the QM atom to which it is bound.
                    int qmInd = atmP[atmIndx].bountToIndx ;
                    resForce[qmInd].charge += localCharge;
                }
                
                strIndx += 9;
                atmIndx++;
                
                // If we found all charges for QM and dummy atoms, break the loop
                // and stop reading the line.
                if (atmIndx == msg->numAllAtoms) {
                    chargeFields = false;
                    break;
                }
                
            }
            
        }
        
        if (gradFields) {
            while ((strIndx < (strlen(line)-9)) && (strlen(line)-1 >=9 ) ) {
                
                strncpy(result, line+strIndx,9) ;
                result[9] = '\0';
                
                gradient[gradCount] = atof(result);
                if (gradCount == 2) {
                    
                    if (atmIndx < msg->numQMAtoms ) {
                        // Gradients in MOPAC are written in kcal/mol/A.
                        resForce[atmIndx].force.x = -1*gradient[0];
                        resForce[atmIndx].force.y = -1*gradient[1];
                        resForce[atmIndx].force.z = -1*gradient[2];
                    }
                    
                    // If we are reading forces applied on Dummy atoms,
                    // place them on the appropriate QM and MM atoms to conserve energy.
                    
                    // This implementation was based on the description in 
                    // DOI: 10.1002/jcc.20857  :  Equations 30 to 32
                    if ( msg->numQMAtoms <= atmIndx &&
                    atmIndx < msg->numAllAtoms ) {
                        // The dummy atom points to the QM atom to which it is bound.
                        // The QM atom and the MM atom (in a QM-MM bond) point to each other.
                        int qmInd = atmP[atmIndx].bountToIndx ;
                        int mmInd = atmP[qmInd].bountToIndx ;
                        
                        Vector dir = pcP[mmInd].position - atmP[qmInd].position ;
                        Real mmqmDist = dir.length() ;
                        
                        Real linkDist = Vector(atmP[atmIndx].position - 
                                        atmP[qmInd].position).length() ;
                        
                        Force mmForce(0), qmForce(0), 
                            linkForce(gradient[0], gradient[1], gradient[2]);
                        linkForce *= -1;
                        
                        Vector base = (linkDist/(mmqmDist*mmqmDist*mmqmDist))*dir;
                        // Unit vectors
                        Vector xuv(1,0,0), yuv(0,1,0), zuv(0,0,1);
                        Real xDelta = pcP[mmInd].position.x - atmP[qmInd].position.x;
                        Real yDelta = pcP[mmInd].position.y - atmP[qmInd].position.y;
                        Real zDelta = pcP[mmInd].position.z - atmP[qmInd].position.z;
                        
                        qmForce += (linkForce*((1 - linkDist/mmqmDist)*xuv + 
                                    (xDelta)*base) )*xuv;
                        
                        qmForce += (linkForce*((1 - linkDist/mmqmDist)*yuv + 
                                    (yDelta)*base) )*yuv;
                        
                        qmForce += (linkForce*((1 - linkDist/mmqmDist)*zuv + 
                                    (zDelta)*base) )*zuv;
                        
                        
                        mmForce += (linkForce*((linkDist/mmqmDist)*xuv -
                                    (xDelta)*base) )*xuv;
                        
                        mmForce += (linkForce*((linkDist/mmqmDist)*yuv -
                                    (yDelta)*base) )*yuv;
                        
                        mmForce += (linkForce*((linkDist/mmqmDist)*zuv -
                                    (zDelta)*base) )*zuv;
                        
                        resForce[qmInd].force += qmForce;
                        resForce[msg->numQMAtoms + mmInd].force += mmForce;
                    }
                    
                    gradCount = 0;
                    atmIndx++;
                } else {
                    gradCount++;
                }
                
                strIndx += 9;
                
                // If we found all gradients for QM atoms, break the loop
                // and keep the next gradient line from being read, if there
                // is one. Following gradients, if present, will refer to link
                // atoms, and who cares about those?.
                if (atmIndx == msg->numAllAtoms) {
                    gradFields = false;
                    break;
                }
            }
        }
        
    }
    
    delete [] line;
    
    fclose(outputFile);
    
    // In case charges are not to be read form the QM software output,
    // we load the origianl atom charges.
    if (msg->qmAtmChrgMode == QMCHRGNONE) {
        int globalNumQMAtoms = Node::Object()->molecule->get_numQMAtoms();
        const int *qmAtmIndx = Node::Object()->molecule->get_qmAtmIndx() ;
        Real *qmAtmChrg = Node::Object()->molecule->get_qmAtmChrg() ;
        
        atmIndx = 0 ;
        for (; atmIndx < msg->numQMAtoms; atmIndx++) {
            
            // Loops over all QM atoms (in all QM groups) comparing their global indices
            for (int qmIter=0; qmIter<globalNumQMAtoms; qmIter++) {
                
                if (qmAtmIndx[qmIter] == atmP[atmIndx].id) {
                    
                    atmP[atmIndx].charge = qmAtmChrg[qmIter];
                    resForce[atmIndx].charge = qmAtmChrg[qmIter];
                    
                    break;
                }
                
            }
            
        }
    }
    
    // remove force file
//     DebugM(4, "Removing output file: " << outputFileName << std::endl) ;
//     iret = remove(outputFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }
    
    if (! (chargesRead && gradsRead) ) {
        NAMD_die("Error reading QM forces file. Not all data could be read!");
    }
    
    DebugM(4, "Applying forces on " << msg->numRealPntChrgs << " point charges" << std::endl) ;
    
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // The initial point charge index for the force message is the number of QM
    // atoms, since the dummy atoms have no representation in NAMD
    int pcIndx = msg->numQMAtoms;
    
    // We only loop over point charges from real atoms, ignoring the ones 
    // created to handle QM-MM bonds.
    for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
        
        BigReal Force = 0;
        
        BigReal pntCharge = pcP[i].charge;
        
        Position posMM = pcP[i].position ;
        
        for (size_t j=0; j<msg->numQMAtoms; ++j ) {
            
            // Not perfect
            // This prevents the MM point charge of a MM-QM bond from feeling 
            // the influence from the QM atom it is bount to. 
            if ( pcP[i].bountToIndx == j ) {
                DebugM(4,"Skiping charged interaction between " << atmP[j].id 
                << " and " << pcP[i].id << std::endl);
                continue ;
            }
            
            BigReal qmCharge = atmP[j].charge ;
            
            Force = pntCharge*qmCharge*constants ;
            
            Position rVec = posMM - atmP[i].position ;
            
            Force /= rVec.length2();
            
            resForce[pcIndx].force += Force*(rVec*rVec.length());
        }
    }
    
    // Adjusts forces from PME, canceling contributions from the QM and 
    // direct Coulomb forces calculated here.
    if (msg->PMEOn) {
        
        DebugM(1,"Correcting forces and energy for PME.\n");
        
        ewaldcof = msg->PMEEwaldCoefficient;
        BigReal TwoBySqrtPi = 1.12837916709551;
        pi_ewaldcof = TwoBySqrtPi * ewaldcof;
        
        for (size_t i=0; i < msg->numQMAtoms; i++) {
            
            BigReal p_i_charge = atmP[i].charge ;
            Position pos_i = atmP[i].position ;
            
            const BigReal kq_i = p_i_charge * constants;
            
            for (size_t j=i+1; j < msg->numQMAtoms; j++) {
                
                BigReal p_j_charge = atmP[j].charge ;
                
                Position pos_j = atmP[j].position ;
                
                BigReal r = Vector(pos_i - pos_j).length() ;
                
                BigReal tmp_a = r * ewaldcof;
                BigReal tmp_b = erfc(tmp_a);
                BigReal corr_energy = tmp_b;
                BigReal corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
                
//                 BigReal recip_energy = (1-tmp_b)/r = erf(tmp_a)/r;
                BigReal recip_energy = (1-tmp_b)/r;
                
                BigReal recip_gradient = -(1-corr_gradient)/(r*2);
                
                // Final force and energy correction for this pair of atoms.
                BigReal energy = kq_i * p_j_charge * recip_energy ;
                
                Force fixForce = -1*kq_i*p_j_charge*(recip_gradient/r)*(pos_i - pos_j) ;
                
                // The force is *subtracted* from the total force acting on
                // both atoms. The sign on fixForce corrects the orientation
                // of the subtracted force.
//                 DebugM(4,"Old forces for QM " << i << ": " << resForce[i].force
//                     << std::endl);
//                 DebugM(4,"Old forces for QM " << j << ": " << resForce[j].force
//                     << std::endl);
//                 DebugM(4,"Force correction: " << fixForce << std::endl);
//                 DebugM(4,"Energy correction: " << energy << "\n");
                resForce[i].force -= fixForce ;
                resForce[j].force -= -1*fixForce ;
                
                // The energy is *subtracted* from the total energy calculated here.
                resMsg->energyCorr -= energy;
            }
        }
        
//         DebugM(4,"Corrected energy QM-QM interactions: " << resMsg->energyCorr << "\n");
        
        pcIndx = msg->numQMAtoms;
        // We only loop over point charges from real atoms, ignoring the ones 
        // created to handle QM-MM bonds.
        for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
            
            BigReal p_i_charge = pcP[i].charge ;
            Position pos_i = pcP[i].position ;
            
            const BigReal kq_i = p_i_charge * constants;
            
            Force fixForce = 0;
            
            for (size_t j=0; j<msg->numQMAtoms; ++j ) {
                
                // Not perfect
                // This prevents the MM point charge of a MM-QM bond from feeling 
                // the influence from the QM atom it is bount to. 
                if ( pcP[i].bountToIndx == j ) continue ;
                
                BigReal p_j_charge = atmP[j].charge ;
                
                Position pos_j = atmP[j].position ;
                
                BigReal r = Vector(pos_i - pos_j).length() ;
                
                BigReal tmp_a = r * ewaldcof;
                BigReal tmp_b = erfc(tmp_a);
                BigReal corr_energy = tmp_b;
                BigReal corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
                
//                 BigReal recip_energy = (1-tmp_b)/r = erf(tmp_a)/r;
                BigReal recip_energy = (1-tmp_b)/r;
                
                BigReal recip_gradient = -(1-corr_gradient)/(r*2);
                
                // Final force and energy correction for this pair of atoms.
                BigReal energy = kq_i * p_j_charge * recip_energy ;
                
                fixForce += -1*p_j_charge*(recip_gradient/r)*(pos_i - pos_j) ;
                
                // The energy is *subtracted* from the total energy calculated here.
                resMsg->energyCorr -= energy;
            }
            
            // The force is *subtracted* from the total force acting on
            // the point charge..
//             DebugM(4,"Old forces for PC " << pcIndx << ": " << resForce[pcIndx].force
//             << std::endl);
//             DebugM(4,"Force correction: " << fixForce << std::endl);
            resForce[pcIndx].force -= kq_i*fixForce ;
        }
        
    }
    
    
    // Calculates the virial contribution form the forces on QM atoms and 
    // point charges calculated here.
    for (size_t i=0; i < msg->numQMAtoms; i++) {
        // virial used by NAMD is -'ve of normal convention, so reverse it!
        // virial[i][j] in file should be sum of -1 * f_i * r_j
        for ( int k=0; k<3; ++k )
            for ( int l=0; l<3; ++l )
                resMsg->virial[k][l] = -1*resForce[i].force[k]*atmP[i].position[l];
    }
    
    pcIndx = msg->numQMAtoms; // Index in the real PC force array.
    for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
        // virial used by NAMD is -'ve of normal convention, so reverse it!
        // virial[i][j] in file should be sum of -1 * f_i * r_j
        for ( int k=0; k<3; ++k )
            for ( int l=0; l<3; ++l )
                resMsg->virial[k][l] = -1*resForce[pcIndx].force[k]*pcP[i].position[l];
    }
    
    
    
    // Send message to rank zero with results.
    QMProxy[0].recvQMRes(resMsg);
    
    delete msg;
    return ;
}

void ComputeQMMgr::calcORCA(QMGrpCalcMsg *msg)
{
    
    FILE *inputFile,*outputFile,*chrgFile;
    int iret;
    
    const size_t lineLen = 256;
    char *line = new char[lineLen];
    
    std::string qmCommand, inputFileName, outputFileName, pntChrgFileName;
    std::string tmpRedirectFileName, pcGradFileName;
    
    // For coulomb calculation
    BigReal constants = msg->constants;
    
    double gradient[3];
    
    DebugM(4,"Running ORCA on PE " << CkMyPe() << std::endl);
    
    if (msg->switching)
        pntChrgSwitching(msg) ;
    
    // For each QM group, create a subdirectory where files will be palced.
    std::string baseDir(msg->baseDir);
    baseDir.append("/") ;
    std::ostringstream itosConv ;
    itosConv << msg->peIter ;
    baseDir += itosConv.str() ;
    
    struct stat info;
    
    if (stat(msg->baseDir, &info) != 0 ) {
        CkPrintf( "Node %d cannot access directory %s\n",
                  CkMyPe(), baseDir.c_str() );
        NAMD_die("QM calculation could not be ran. Check your qmBaseDir!");
    }
    else if (! (stat(baseDir.c_str(), &info) == 0 && S_ISDIR(info.st_mode)) ) {
        DebugM(4,"Creating directory " << baseDir.c_str() << std::endl);
        int retVal = mkdir(baseDir.c_str(), S_IRWXU);
    }
    
    #ifdef DEBUG_QM
    Write_PDB(std::string(baseDir)+"/input.pdb", msg ) ;
    #endif
    
    inputFileName.clear();
    inputFileName.append(baseDir.c_str()) ;
    inputFileName.append("/qmmm_") ;
    inputFileName += itosConv.str() ;
    inputFileName.append(".input") ;
    
    // Opens file for coordinate and parameter input
    inputFile = fopen(inputFileName.c_str(),"w");
    if ( ! inputFile ) {
        iout << iERROR << "Could not open input file for writing: " 
        << inputFileName << "\n" << endi ;
        NAMD_die(strerror(errno));
    }
    
    // Builds the command that will be executed
    qmCommand.clear();
    qmCommand.append("cd ");
    qmCommand.append(baseDir);
    qmCommand.append(" ; ");
    qmCommand.append(msg->execPath) ;
    qmCommand.append(" ") ;
    qmCommand.append(inputFileName) ;
    
    // Build the redirect file name we need for the screen output.
    // That's the only place where we can get partial charges for QM atoms.
    tmpRedirectFileName = inputFileName ;
    tmpRedirectFileName.append(".TmpOut") ;
    
    qmCommand.append(" > ") ;
    qmCommand.append(tmpRedirectFileName) ;
    
    // Builds the file name where orca will place the gradient
    // This will be relative to the input file
    outputFileName = inputFileName ;
    outputFileName.append(".engrad") ;
    
    QMAtomData *pcP = msg->data + msg->numAllAtoms ;
    if (msg->numAllPntChrgs) {
        // Builds the file name where we will write the point charges.
        pntChrgFileName = inputFileName ;
        pntChrgFileName.append(".pntchrg") ;
        
        pcGradFileName = inputFileName;
        pcGradFileName.append(".pcgrad");
        
        chrgFile = fopen(pntChrgFileName.c_str(),"w");
        if ( ! chrgFile ) {
            iout << iERROR << "Could not open charge file for writing: " 
            << pntChrgFileName << "\n" << endi ;
            NAMD_die(strerror(errno));
        }
        
        int numPntChrgs = 0;
        for (int i=0; i<msg->numAllPntChrgs; i++ ) {
            if (pcP[i].type != QMPCTYPE_IGNORE)
                numPntChrgs++;
        }
        
        iret = fprintf(chrgFile,"%d\n", numPntChrgs);
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    }
    
    // writes configuration lines to the input file.
    std::stringstream ss(msg->configLines) ;
    std::string confLineStr;
    while (std::getline(ss, confLineStr) ) {
        confLineStr.append("\n");
        iret = fprintf(inputFile,confLineStr.c_str());
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    }
    
    if (msg->numAllPntChrgs) {
        iret = fprintf(inputFile,"%%pointcharges \"%s\"\n", pntChrgFileName.c_str());
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    }
    
    iret = fprintf(inputFile,"\n\n%%coords\n  CTyp xyz\n");
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    iret = fprintf(inputFile,"  Charge %f\n",msg->charge);
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    iret = fprintf(inputFile,"  Mult %f\n",msg->multiplicity);
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    iret = fprintf(inputFile,"  Units Angs\n  coords\n\n");
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    DebugM(4, "Writing " << msg->numAllAtoms << " QM atom coords in file " << 
    inputFileName.c_str() << " and " << msg->numAllPntChrgs << 
    " point charges in file " << pntChrgFileName.c_str() << "\n");
    
    // write QM and dummy atom coordinates to input file.
    QMAtomData *atmP = msg->data ;
    for (size_t i=0; i<msg->numAllAtoms; ++i, ++atmP ) {
        
        double x = atmP->position.x;
        double y = atmP->position.y;
        double z = atmP->position.z;
        
        iret = fprintf(inputFile,"  %s %f %f %f\n",
                       atmP->element,x,y,z);
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
        
    }
    
    iret = fprintf(inputFile,"  end\nend\n");
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    if (msg->numAllPntChrgs) {
        // Write point charges to file.
        pcP = msg->data + msg->numAllAtoms ;
        for ( size_t j=0; j < msg->numAllPntChrgs; j++, ++pcP) {
            
            if (pcP->type == QMPCTYPE_IGNORE)
                    continue;
            
            double charge = pcP->charge;
            
            double x = pcP->position.x;
            double y = pcP->position.y;
            double z = pcP->position.z;
            
            iret = fprintf(chrgFile,"%f %f %f %f\n",
                           charge,x,y,z);
            if ( iret < 0 ) { NAMD_die(strerror(errno)); }
        }
        
        fclose(chrgFile);
    }
    
    DebugM(4,"Closing input and charge file\n");
    fclose(inputFile);
    
    if (msg->prepProcOn) {
        
        std::string prepProc(msg->prepProc) ;
        prepProc.append(" ") ;
        prepProc.append(inputFileName) ;
        iret = system(prepProc.c_str());
        if ( iret == -1 ) { NAMD_die(strerror(errno)); }
        if ( iret ) { NAMD_die("Error running preparation command for QM calculation."); }
    }
    
        // runs QM command
    DebugM(4,"Running command ->" << qmCommand.c_str() << "<-" << std::endl);
    iret = system(qmCommand.c_str());
    
    if ( iret == -1 ) { NAMD_die(strerror(errno)); }
    if ( iret ) { NAMD_die("Error running command for QM forces calculation."); }

    if (msg->secProcOn) {
        
        std::string secProc(msg->secProc) ;
        secProc.append(" ") ;
        secProc.append(inputFileName) ;
        iret = system(secProc.c_str());
        if ( iret == -1 ) { NAMD_die(strerror(errno)); }
        if ( iret ) { NAMD_die("Error running second command for QM calculation."); }
    }

    // remove coordinate file
//     iret = remove(inputFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }

    // remove coordinate file
//     iret = remove(pntChrgFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }
    
    // opens output file
    DebugM(4,"Reading QM data from file " << outputFileName.c_str() << std::endl);
    outputFile = fopen(outputFileName.c_str(),"r");
    if ( ! outputFile ) {
        iout << iERROR << "Could not find QM output file!\n" << endi;
        NAMD_die(strerror(errno)); 
    }

    // Resets the pointers.
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // Allocates the return message.
    QMGrpResMsg *resMsg = new (msg->numQMAtoms + msg->numRealPntChrgs, 0) QMGrpResMsg;
    resMsg->grpIndx = msg->grpIndx;
    resMsg->numForces = msg->numQMAtoms + msg->numRealPntChrgs;
    resMsg->energyOrig = 0;
    resMsg->energyCorr = 0;
    for ( int k=0; k<3; ++k )
        for ( int l=0; l<3; ++l )
            resMsg->virial[k][l] = 0;
    QMForce *resForce = resMsg->force ;
    
    // We split the data into two pointers so we can "skip" the dummy atoms
    // which have no representation in NAMD.
    for (int i=0; i<resMsg->numForces; i++) {
        resForce[i].force = 0;
        resForce[i].charge = 0 ;
        if (i < msg->numQMAtoms) {
            // We use the replace field to indicate QM atoms,
            // which will have charge information.
            resForce[i].replace = 1 ;
            resForce[i].id = atmP->id;
            atmP++;
        }
        else {
            // We use the replace field to indicate QM atoms,
            // which will have charge information.
            resForce[i].replace = 0 ;
            resForce[i].id = pcP->id;
            pcP++;
        }
    }
    
    // Resets the pointers.
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    size_t atmIndx = 0, gradCount = 0;
    int numAtomsInFile = 0;
    
    // Reads the data form the output file created by the QM software.
    // Gradients over the QM atoms, and Charges for QM atoms will be read.
    
    // skip lines before number of atoms
    for (int i = 0; i < 3; i++) {
        fgets(line, lineLen, outputFile); 
    }
    
    iret = fscanf(outputFile,"%d\n", &numAtomsInFile);
    if ( iret != 1 ) {
        NAMD_die("Error reading QM forces file.");
    }
    
    #ifdef DEBUG_QM
    if(numAtomsInFile != msg->numAllAtoms) {
        NAMD_die("Error reading QM forces file. Number of atoms in QM output\
        does not match the expected.");
    }
    #endif
    
    // skip lines before energy
    for (int i = 0; i < 3; i++) {
        fgets(line, lineLen, outputFile); 
    }
    
    iret = fscanf(outputFile,"%lf\n", &resMsg->energyOrig);
    if ( iret != 1 ) {
        NAMD_die("Error reading QM forces file.");
    }
//     iout << "Energy in step (Hartree): " << resMsg->energyOrig << "\n" <<  endi;
    // All energies are given in Eh (Hartree)
    // NAMD needs energies in kcal/mol
    // The conversion factor is 627.509469
    resMsg->energyOrig *= 627.509469;
//     iout << "Energy in step (Kcal/mol): " << resMsg->energyOrig << "\n" <<  endi;
    
    resMsg->energyCorr = resMsg->energyOrig;
    
    // skip lines before gradient
    for (int i = 0; i < 3; i++) {
        fgets(line, lineLen, outputFile) ;
    }
    
    // Break the loop when we find all gradients for QM atoms, 
    // and keep the next gradient lines from being read, if there
    // are more. Following gradients, if present, will refer to link
    // atoms.
    atmIndx = 0;
    gradCount = 0;
    for ( size_t i=0; i<msg->numAllAtoms*3; ++i ) {
        
        iret = fscanf(outputFile,"%lf\n", &gradient[gradCount]);
        if ( iret != 1 ) { NAMD_die("Error reading QM forces file."); }
        
        if (gradCount == 2){
            
            // All gradients are given in Eh/a0 (Hartree over Bohr radius)
            // NAMD needs forces in kcal/mol/angstrons
            // The conversion factor is 627.509469/0.529177 = 1185.82151
            if (atmIndx < msg->numQMAtoms ) {
                resForce[atmIndx].force.x = -1*gradient[0]*1185.82151;
                resForce[atmIndx].force.y = -1*gradient[1]*1185.82151;
                resForce[atmIndx].force.z = -1*gradient[2]*1185.82151;
            }
            
            // If we are reading forces applied on Dummy atoms,
            // place them on the appropriate QM and MM atoms to conserve energy.
            
            // This implementation was based on the description in 
            // DOI: 10.1002/jcc.20857  :  Equations 30 to 32
            if ( msg->numQMAtoms <= atmIndx &&
            atmIndx < msg->numAllAtoms ) {
                // The dummy atom points to the QM atom to which it is bound.
                // The QM atom and the MM atom (in a QM-MM bond) point to each other.
                int qmInd = atmP[atmIndx].bountToIndx ;
                int mmInd = atmP[qmInd].bountToIndx ;
                
                Vector dir = pcP[mmInd].position - atmP[qmInd].position ;
                Real mmqmDist = dir.length() ;
                
                Real linkDist = Vector(atmP[atmIndx].position - atmP[qmInd].position).length() ;
                
                Force mmForce(0), qmForce(0), linkForce(gradient[0], gradient[1], gradient[2]);
                linkForce *= -1*1185.82151;
                
                Vector base = (linkDist/(mmqmDist*mmqmDist*mmqmDist))*dir;
                // Unit vectors
                Vector xuv(1,0,0), yuv(0,1,0), zuv(0,0,1);
                Real xDelta = pcP[mmInd].position.x - atmP[qmInd].position.x;
                Real yDelta = pcP[mmInd].position.y - atmP[qmInd].position.y;
                Real zDelta = pcP[mmInd].position.z - atmP[qmInd].position.z;
                
                qmForce += (linkForce*((1 - linkDist/mmqmDist)*xuv + 
                            (xDelta)*base) )*xuv;
                
                qmForce += (linkForce*((1 - linkDist/mmqmDist)*yuv + 
                            (yDelta)*base) )*yuv;
                
                qmForce += (linkForce*((1 - linkDist/mmqmDist)*zuv + 
                            (zDelta)*base) )*zuv;
                
                
                mmForce += (linkForce*((linkDist/mmqmDist)*xuv -
                            (xDelta)*base) )*xuv;
                
                mmForce += (linkForce*((linkDist/mmqmDist)*yuv -
                            (yDelta)*base) )*yuv;
                
                mmForce += (linkForce*((linkDist/mmqmDist)*zuv -
                            (zDelta)*base) )*zuv;
                
                resForce[qmInd].force += qmForce;
                resForce[msg->numQMAtoms + mmInd].force += mmForce;
            }
            
            gradCount = 0;
            atmIndx++ ;
        } else
            gradCount++ ;
        
    }
    
    fclose(outputFile);
    
    // In case charges are not to be read form the QM software output,
    // we load the origianl atom charges.
    if (msg->qmAtmChrgMode == QMCHRGNONE) {
        int globalNumQMAtoms = Node::Object()->molecule->get_numQMAtoms();
        const int *qmAtmIndx = Node::Object()->molecule->get_qmAtmIndx() ;
        Real *qmAtmChrg = Node::Object()->molecule->get_qmAtmChrg() ;
        
        atmIndx = 0 ;
        for (; atmIndx < msg->numQMAtoms; atmIndx++) {
            
            // Loops over all QM atoms (in all QM groups) comparing their global indices
            for (int qmIter=0; qmIter<globalNumQMAtoms; qmIter++) {
                
                if (qmAtmIndx[qmIter] == atmP[atmIndx].id) {
                    
                    atmP[atmIndx].charge = qmAtmChrg[qmIter];
                    resForce[atmIndx].charge = qmAtmChrg[qmIter];
                    
                    break;
                }
                
            }
            
        }
    }
    else {
        // opens redirected output file
        outputFile = fopen(tmpRedirectFileName.c_str(),"r");
        if ( ! outputFile ) {
            iout << iERROR << "Could not find Redirect output file:"
            << tmpRedirectFileName << std::endl << endi;
            NAMD_die(strerror(errno)); 
        }
        
        DebugM(4,"Opened tmeporary output for charge reading: " << tmpRedirectFileName.c_str() << "\n");
        
        int lineState = 0;
        char result[20] ;
        while ( fgets(line, lineLen, outputFile) != NULL){
            
            // We loop over the lines of the output file untill we find
            // the line that initiates the charges lines. Then
            // we read all lines and expect to get one charge
            // per line, untill we find a line that sums all charges.
            
            switch (msg->qmAtmChrgMode) {
                
                case QMCHRGMULLIKEN:
                {
                
                if ( strstr(line,"MULLIKEN ATOMIC CHARGES") != NULL ) {
                    lineState = 1;
                    atmIndx = 0;
                    
                    // Now we expect the following line to have a series of dashes
                    // and the folowing lines to have atom charges (among other info)
                    continue;
                }
                
                if ( strstr(line,"Sum of atomic charges") != NULL ) {
                    
                    // Now that all charges have been read, grabs the total charge
                    // just for fun... and checking purposes.
                    #ifdef DEBUG_QM
                    strncpy(result, line + 31,12) ;
                    result[12] = '\0';
                    
                    DebugM(4,"Total charge of QM region calculated by ORCA is: " 
                    << atof(result) << std::endl )
                    #endif
                    
                    // Nothing more to read from the output file
                    break;
                }
                
                // Line state 1 means we have to skip ONLY one line.
                if (lineState == 1) {
                    lineState = 2;
                    continue;
                }
                
                // Line state 2 means we have to read the line and grab the info.
                if (lineState == 2) {
                    
                    strncpy(result, line+8,12) ;
                    result[12] = '\0';
                    
                    Real localCharge = atof(result);
                    
                    // If we are reading charges from QM atoms, store them.
                    if (atmIndx < msg->numQMAtoms ) {
                        atmP[atmIndx].charge = localCharge;
                        resForce[atmIndx].charge = localCharge;
                    }
                    
                    // If we are reading charges from Dummy atoms,
                    // place the on the appropriate QM atom.
                     if ( msg->numQMAtoms <= atmIndx &&
                        atmIndx < msg->numAllAtoms ) {
                        int qmInd = atmP[atmIndx].bountToIndx ;
                        atmP[qmInd].charge += localCharge;
                        resForce[qmInd].charge += localCharge;
                    }
                    
                    atmIndx++ ;
                    
                    // If we found all charges for QM atoms, change the lineState
                    // untill we reach the "total charge" line.
                    if (atmIndx == msg->numAllAtoms ) {
                        lineState = 0;
                    }
                    
                    continue;
                }
                
                } break ;
                
                case QMCHRGCHELPG :
                {
                
                if ( strstr(line,"CHELPG Charges") != NULL ) {
                    lineState = 1;
                    atmIndx = 0;
                    
                    // Now we expect the following line to have a series of dashes
                    // and the folowing lines to have atom charges (among other info)
                    continue;
                }
                
                if ( strstr(line,"Total charge") != NULL ) {
                    
                    // Now that all charges have been read, grabs the total charge
                    // just for fun... and checking purposes.
                    #ifdef DEBUG_QM
                    strncpy(result, line + 14,13) ;
                    result[13] = '\0';
                    
                    DebugM(4,"Total charge of QM region calculated by ORCA is: " 
                    << atof(result) << std::endl )
                    #endif
                    
                    // Nothing more to read from the output file
                    break;
                }
                
                // Line state 1 means we have to skip ONLY one line.
                if (lineState == 1) {
                    lineState = 2;
                    continue;
                }
                
                // Line state 2 means we have to read the line and grab the info.
                if (lineState == 2) {
                    
                    strncpy(result, line+12,15) ;
                    result[15] = '\0';
                    
                    Real localCharge = atof(result);
                    
                    // If we are reading charges from QM atoms, store them.
                    if (atmIndx < msg->numQMAtoms ) {
                        atmP[atmIndx].charge = localCharge;
                        resForce[atmIndx].charge = localCharge;
                    }
                    
                    // If we are reading charges from Dummy atoms,
                    // place the on the appropriate QM atom.
                     if ( msg->numQMAtoms <= atmIndx &&
                        atmIndx < msg->numAllAtoms ) {
                        int qmInd = atmP[atmIndx].bountToIndx ;
                        atmP[qmInd].charge += localCharge;
                        resForce[qmInd].charge += localCharge;
                    }
                    
                    atmIndx++ ;
                    
                    // If we found all charges for QM atoms, we ignore the following line
                    // untill we reach the "total charge" line.
                    if (atmIndx == msg->numAllAtoms ) {
                        lineState = 1;
                    }
                    
                    continue;
                }
                
                } break;
            }
        }
        
        fclose(outputFile);
    }
    
    delete [] line;
    
    // remove force file
//     DebugM(4, "Removing output file: " << outputFileName << std::endl) ;
//     iret = remove(outputFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }
    
    
    DebugM(4, "Applying forces on " << msg->numRealPntChrgs << " point charges" << std::endl) ;
    
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // The initial point charge index for the force message is the number of QM
    // atoms, since the dummy atoms have no representation in NAMD
    int pcIndx = msg->numQMAtoms;
    
    // We only loop over point charges from real atoms, ignoring the ones 
    // created to handle QM-MM bonds.
    for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
        
        BigReal Force = 0;
        
        BigReal pntCharge = pcP[i].charge;
        
        BigReal xMM = pcP[i].position.x;
        BigReal yMM = pcP[i].position.y;
        BigReal zMM = pcP[i].position.z;
        
        for (size_t j=0; j<msg->numQMAtoms; ++j ) {
            
            // Not perfect
            // This prevents the MM point charge of a MM-QM bond from feeling 
            // the influence from the QM atom it is bount to. 
            if ( pcP[i].bountToIndx == j ) continue ;
            
            BigReal qmCharge = atmP[j].charge ;
            
            Force = pntCharge*qmCharge*constants ;
            
            BigReal xQM = atmP[j].position.x;
            BigReal yQM = atmP[j].position.y;
            BigReal zQM = atmP[j].position.z;
            
            BigReal x_ij = (xMM - xQM);
            BigReal y_ij = (yMM - yQM);
            BigReal z_ij = (zMM - zQM);
            
            BigReal r2 = (x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
            BigReal rNorm = sqrt(r2) ;
            
            Force /= r2;
            
            resForce[pcIndx].force.x += Force*x_ij/rNorm;
            resForce[pcIndx].force.y += Force*y_ij/rNorm;
            resForce[pcIndx].force.z += Force*z_ij/rNorm;
        }
        
    }
    
    // Adjusts forces from PME, canceling contributions from the QM and 
    // direct Coulomb forces calculated here.
    if (msg->PMEOn) {
        
        DebugM(1,"Correcting forces and energy for PME.\n");
        
        ewaldcof = msg->PMEEwaldCoefficient;
        BigReal TwoBySqrtPi = 1.12837916709551;
        pi_ewaldcof = TwoBySqrtPi * ewaldcof;
        
        for (size_t i=0; i < msg->numQMAtoms; i++) {
            
            BigReal p_i_charge = atmP[i].charge ;
            Position pos_i = atmP[i].position ;
            
            for (size_t j=i+1; j < msg->numQMAtoms; j++) {
                
                BigReal p_j_charge = atmP[j].charge ;
                
                Position pos_j = atmP[j].position ;
                
                BigReal r = Vector(pos_i - pos_j).length() ;
                
                BigReal tmp_a = r * ewaldcof;
                BigReal tmp_b = erfc(tmp_a);
                BigReal corr_energy = tmp_b;
                BigReal corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
                
//                 BigReal recip_energy = (1-tmp_b)/r = erf(tmp_a)/r;
                BigReal recip_energy = (1-tmp_b)/r;
                
                BigReal recip_gradient = -(1-corr_gradient)/(r*2);
                
                const BigReal kq_i = p_i_charge * constants;
                
                // Final force and energy correction for this pair of atoms.
                BigReal energy = kq_i * p_j_charge * recip_energy ;
                
                Force fixForce = -1*kq_i*p_j_charge*(recip_gradient/r)*(pos_i - pos_j) ;
                
                // The force is *subtracted* from the total force acting on
                // both atoms. The sign on fixForce corrects the orientation
                // of the subtracted force.
//                 DebugM(4,"Old forces for QM " << i << ": " << resForce[i].force
//                     << std::endl);
//                 DebugM(4,"Old forces for QM " << j << ": " << resForce[j].force
//                     << std::endl);
//                 DebugM(4,"Force correction: " << fixForce << std::endl);
                resForce[i].force -= fixForce ;
                resForce[j].force -= -1*fixForce;
                
                // The energy is *subtracted* from the total energy calculated here.
                resMsg->energyCorr -= energy;
            }
        }
        
        pcIndx = msg->numQMAtoms;
        // We only loop over point charges from real atoms, ignoring the ones 
        // created to handle QM-MM bonds.
        for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
            
            BigReal p_i_charge = pcP[i].charge ;
            Position pos_i = pcP[i].position ;
            
            const BigReal kq_i = p_i_charge * constants;
            
            Force fixForce = 0;
            
            for (size_t j=0; j<msg->numQMAtoms; ++j ) {
                
                // Not perfect
                // This prevents the MM point charge of a MM-QM bond from feeling 
                // the influence from the QM atom it is bount to. 
                if ( pcP[i].bountToIndx == j ) continue ;
                
                BigReal p_j_charge = atmP[j].charge ;
                
                Position pos_j = atmP[j].position ;
                
                BigReal r = Vector(pos_i - pos_j).length() ;
                
                BigReal tmp_a = r * ewaldcof;
                BigReal tmp_b = erfc(tmp_a);
                BigReal corr_energy = tmp_b;
                BigReal corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
                
//                 BigReal recip_energy = (1-tmp_b)/r = erf(tmp_a)/r;
                BigReal recip_energy = (1-tmp_b)/r;
                
                BigReal recip_gradient = -(1-corr_gradient)/(r*2);
                
                // Final force and energy correction for this pair of atoms.
                BigReal energy = kq_i * p_j_charge * recip_energy ;
                
                fixForce += -1*p_j_charge*(recip_gradient/r)*(pos_i - pos_j) ;
                
                // The energy is *subtracted* from the total energy calculated here.
                resMsg->energyCorr -= energy;
                
            }
            
            // The force is *subtracted* from the total force acting on
                // the point charge.
//                 DebugM(4,"Old forces for PC " << pcIndx << ": " << resForce[pcIndx].force
//                     << std::endl);
//                 DebugM(4,"Force correction: " << fixForce << std::endl);
            resForce[pcIndx].force -= kq_i*fixForce ;
        }
        
    }
    
    DebugM(1,"Determining virial...\n");
    
    // Calculates the virial contribution form the forces on QM atoms and 
    // point charges calculated here.
    for (size_t i=0; i < msg->numQMAtoms; i++) {
        // virial used by NAMD is -'ve of normal convention, so reverse it!
        // virial[i][j] in file should be sum of -1 * f_i * r_j
        for ( int k=0; k<3; ++k )
            for ( int l=0; l<3; ++l )
                resMsg->virial[k][l] = -1*resForce[i].force[k]*atmP[i].position[l];
    }
    
    pcIndx = msg->numQMAtoms; // Index in the real PC force array.
    for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
        // virial used by NAMD is -'ve of normal convention, so reverse it!
        // virial[i][j] in file should be sum of -1 * f_i * r_j
        for ( int k=0; k<3; ++k )
            for ( int l=0; l<3; ++l )
                resMsg->virial[k][l] = -1*resForce[pcIndx].force[k]*pcP[i].position[l];
    }
    
    DebugM(1,"End of QM processing. Sending result message.\n");
    
    // Send message to rank zero with results.
    QMProxy[0].recvQMRes(resMsg);
    
    delete msg;
    return ;
}

void ComputeQMMgr::calcUSR(QMGrpCalcMsg *msg) {
    
    FILE *inputFile,*outputFile;
    int iret;
    
    std::string qmCommand, inputFileName, outputFileName;
    
    // For coulomb calculation
    BigReal constants = msg->constants;
    
    DebugM(4,"Running USER DEFINED SOFTWARE on PE " << CkMyPe() << std::endl);
    
    if (msg->switching)
        pntChrgSwitching(msg) ;
    
    // For each QM group, create a subdirectory where files will be palced.
    std::string baseDir(msg->baseDir);
    baseDir.append("/") ;
    std::ostringstream itosConv ;
    itosConv << msg->peIter ;
    baseDir += itosConv.str() ;
    
    struct stat info;
    
    if (stat(msg->baseDir, &info) != 0 ) {
        CkPrintf( "Node %d cannot access directory %s\n",
                  CkMyPe(), baseDir.c_str() );
        NAMD_die("QM calculation could not be ran. Check your qmBaseDir!");
    }
    else if (! (stat(baseDir.c_str(), &info) == 0 && S_ISDIR(info.st_mode)) ) {
        DebugM(4,"Creating directory " << baseDir.c_str() << std::endl);
        int retVal = mkdir(baseDir.c_str(), S_IRWXU);
    }
    
    #ifdef DEBUG_QM
    Write_PDB(std::string(baseDir)+"/input.pdb", msg ) ;
    #endif
    
    inputFileName.clear();
    inputFileName.append(baseDir.c_str()) ;
    inputFileName.append("/qmmm_") ;
    inputFileName += itosConv.str() ;
    inputFileName.append(".input") ;
    
    // Opens file for coordinate and parameter input
    inputFile = fopen(inputFileName.c_str(),"w");
    if ( ! inputFile ) {
        iout << iERROR << "Could not open input file for writing: " 
        << inputFileName << "\n" << endi ;
        NAMD_die(strerror(errno));
    }
    
    // Builds the command that will be executed
    qmCommand.clear();
    qmCommand.append("cd ");
    qmCommand.append(baseDir);
    qmCommand.append(" ; ");
    qmCommand.append(msg->execPath) ;
    qmCommand.append(" ") ;
    qmCommand.append(inputFileName) ;
    
    // Builds the file name where orca will place the gradient
    // This will be relative to the input file
    outputFileName = inputFileName ;
    outputFileName.append(".result") ;
    
    int numPntChrgs = 0;
    QMAtomData *pcP = msg->data + msg->numAllAtoms ;
    for (int i=0; i<msg->numAllPntChrgs; i++ ) {
        if (pcP[i].type != QMPCTYPE_IGNORE)
            numPntChrgs++;
    }
    
    iret = fprintf(inputFile,"%d %d\n",msg->numAllAtoms, numPntChrgs);
    if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    
    DebugM(4, "Writing " << msg->numAllAtoms << " QM atom coords in file " << 
        inputFileName.c_str() << " and " << msg->numAllPntChrgs << 
        " point charges." << std::endl);
    
    // write QM and dummy atom coordinates to input file.
    QMAtomData *atmP = msg->data ;
    for (size_t i=0; i<msg->numAllAtoms; ++i, ++atmP ) {
        
        double x = atmP->position.x;
        double y = atmP->position.y;
        double z = atmP->position.z;
        
        iret = fprintf(inputFile,"%f %f %f %s\n",
                       x,y,z,atmP->element);
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
        
    }
    
    // Write point charges to file.
    pcP = msg->data + msg->numAllAtoms ;
    for ( size_t j=0; j < msg->numAllPntChrgs; j++, ++pcP) {
        
        if (pcP->type == QMPCTYPE_IGNORE)
                continue;
        
        double charge = pcP->charge;
        
        double x = pcP->position.x;
        double y = pcP->position.y;
        double z = pcP->position.z;
        
        iret = fprintf(inputFile,"%f %f %f %f\n",
                       x,y,z,charge);
        if ( iret < 0 ) { NAMD_die(strerror(errno)); }
    }
    
    DebugM(4,"Closing input file\n");
    fclose(inputFile);
    
    if (msg->prepProcOn) {
        
        std::string prepProc(msg->prepProc) ;
        prepProc.append(" ") ;
        prepProc.append(inputFileName) ;
        iret = system(prepProc.c_str());
        if ( iret == -1 ) { NAMD_die(strerror(errno)); }
        if ( iret ) { NAMD_die("Error running preparation command for QM calculation."); }
    }
    
        // runs QM command
    DebugM(4,"Running command ->" << qmCommand.c_str() << "<-" << std::endl);
    iret = system(qmCommand.c_str());
    
    if ( iret == -1 ) { NAMD_die(strerror(errno)); }
    if ( iret ) { NAMD_die("Error running command for QM forces calculation."); }

    if (msg->secProcOn) {
        
        std::string secProc(msg->secProc) ;
        secProc.append(" ") ;
        secProc.append(inputFileName) ;
        iret = system(secProc.c_str());
        if ( iret == -1 ) { NAMD_die(strerror(errno)); }
        if ( iret ) { NAMD_die("Error running second command for QM calculation."); }
    }

    // remove coordinate file
//     iret = remove(inputFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }

    // remove coordinate file
//     iret = remove(pntChrgFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }
    
    // opens output file
    DebugM(4,"Reading QM data from file " << outputFileName.c_str() << std::endl);
    outputFile = fopen(outputFileName.c_str(),"r");
    if ( ! outputFile ) {
        iout << iERROR << "Could not find QM output file!\n" << endi;
        NAMD_die(strerror(errno)); 
    }

    // Resets the pointers.
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // Allocates the return message.
    QMGrpResMsg *resMsg = new (msg->numQMAtoms + msg->numRealPntChrgs, 0) QMGrpResMsg;
    resMsg->grpIndx = msg->grpIndx;
    resMsg->numForces = msg->numQMAtoms + msg->numRealPntChrgs;
    resMsg->energyOrig = 0;
    resMsg->energyCorr = 0;
    for ( int k=0; k<3; ++k )
        for ( int l=0; l<3; ++l )
            resMsg->virial[k][l] = 0;
    QMForce *resForce = resMsg->force ;
    
    // We split the data into two pointers so we can "skip" the dummy atoms
    // which have no representation in NAMD.
    for (int i=0; i<resMsg->numForces; i++) {
        resForce[i].force = 0;
        resForce[i].charge = 0 ;
        if (i < msg->numQMAtoms) {
            // We use the replace field to indicate QM atoms,
            // which will have charge information.
            resForce[i].replace = 1 ;
            resForce[i].id = atmP->id;
            atmP++;
        }
        else {
            // We use the replace field to indicate QM atoms,
            // which will have charge information.
            resForce[i].replace = 0 ;
            resForce[i].id = pcP->id;
            pcP++;
        }
    }
    
    // Resets the pointers.
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // Reads the data form the output file created by the QM software.
    // Gradients over the QM atoms, and Charges for QM atoms will be read.
    
    iret = fscanf(outputFile,"%lf\n", &resMsg->energyOrig);
    if ( iret != 1 ) {
        NAMD_die("Error reading energy from QM results file.");
    }
    
    resMsg->energyCorr = resMsg->energyOrig;
    
    size_t atmIndx;
    double localForce[3];
    double localCharge;
    for (atmIndx = 0; atmIndx < msg->numAllAtoms; atmIndx++) {
        
        iret = fscanf(outputFile,"%lf %lf %lf %lf\n", 
                      localForce+0,
                      localForce+1,
                      localForce+2,
                      &localCharge);
        if ( iret != 4 ) {
            NAMD_die("Error reading gradient and charge from QM results file.");
        }
        
        // If we are reading charges and forces on QM atoms, store
        // them directly.
        if (atmIndx < msg->numQMAtoms ) {
            
            resForce[atmIndx].force.x = localForce[0];
            resForce[atmIndx].force.y = localForce[1];
            resForce[atmIndx].force.z = localForce[2];
            
            atmP[atmIndx].charge = localCharge;
            resForce[atmIndx].charge = localCharge;
        }
        
        // If we are reading forces applied on Dummy atoms,
        // place them on the appropriate QM and MM atoms to conserve energy.
        
        // This implementation was based on the description in 
        // DOI: 10.1002/jcc.20857  :  Equations 30 to 32
        if ( msg->numQMAtoms <= atmIndx &&
        atmIndx < msg->numAllAtoms ) {
            
            // If we are reading charges from Dummy atoms,
            // place them on the appropriate QM atom.
            // The dummy atom points to the QM atom to which it is bound.
            int qmInd = atmP[atmIndx].bountToIndx ;
            resForce[qmInd].charge += localCharge;
            
            // The dummy atom points to the QM atom to which it is bound.
            // The QM atom and the MM atom (in a QM-MM bond) point to each other.
            int mmInd = atmP[qmInd].bountToIndx ;
            
            Vector dir = pcP[mmInd].position - atmP[qmInd].position ;
            Real mmqmDist = dir.length() ;
            
            Real linkDist = Vector(atmP[atmIndx].position - 
                            atmP[qmInd].position).length() ;
            
            Force mmForce(0), qmForce(0), 
                linkForce(localForce[0], localForce[1], localForce[2]);
            
            Vector base = (linkDist/(mmqmDist*mmqmDist*mmqmDist))*dir;
            // Unit vectors
            Vector xuv(1,0,0), yuv(0,1,0), zuv(0,0,1);
            Real xDelta = pcP[mmInd].position.x - atmP[qmInd].position.x;
            Real yDelta = pcP[mmInd].position.y - atmP[qmInd].position.y;
            Real zDelta = pcP[mmInd].position.z - atmP[qmInd].position.z;
            
            qmForce += (linkForce*((1 - linkDist/mmqmDist)*xuv + 
                        (xDelta)*base) )*xuv;
            
            qmForce += (linkForce*((1 - linkDist/mmqmDist)*yuv + 
                        (yDelta)*base) )*yuv;
            
            qmForce += (linkForce*((1 - linkDist/mmqmDist)*zuv + 
                        (zDelta)*base) )*zuv;
            
            
            mmForce += (linkForce*((linkDist/mmqmDist)*xuv -
                        (xDelta)*base) )*xuv;
            
            mmForce += (linkForce*((linkDist/mmqmDist)*yuv -
                        (yDelta)*base) )*yuv;
            
            mmForce += (linkForce*((linkDist/mmqmDist)*zuv -
                        (zDelta)*base) )*zuv;
            
            resForce[qmInd].force += qmForce;
            resForce[msg->numQMAtoms + mmInd].force += mmForce;
        }
    }
    
    fclose(outputFile);
    
    // In case charges are not to be read form the QM software output,
    // we load the origianl atom charges.
    if (msg->qmAtmChrgMode == QMCHRGNONE) {
        int globalNumQMAtoms = Node::Object()->molecule->get_numQMAtoms();
        const int *qmAtmIndx = Node::Object()->molecule->get_qmAtmIndx() ;
        Real *qmAtmChrg = Node::Object()->molecule->get_qmAtmChrg() ;
        
        atmIndx = 0 ;
        for (; atmIndx < msg->numQMAtoms; atmIndx++) {
            
            // Loops over all QM atoms (in all QM groups) comparing their global indices
            for (int qmIter=0; qmIter<globalNumQMAtoms; qmIter++) {
                
                if (qmAtmIndx[qmIter] == atmP[atmIndx].id) {
                    
                    atmP[atmIndx].charge = qmAtmChrg[qmIter];
                    resForce[atmIndx].charge = qmAtmChrg[qmIter];
                    
                    break;
                }
                
            }
            
        }
    }
    
    // remove force file
//     DebugM(4, "Removing output file: " << outputFileName << std::endl) ;
//     iret = remove(outputFileName);
//     if ( iret ) { NAMD_die(strerror(errno)); }
    
    
    DebugM(4, "Applying forces on " << msg->numRealPntChrgs << " point charges" << std::endl) ;
    
    atmP = msg->data ;
    pcP = msg->data + msg->numAllAtoms ;
    
    // The initial point charge index for the force message is the number of QM
    // atoms, since the dummy atoms have no representation in NAMD
    int pcIndx = msg->numQMAtoms;
    
    // We only loop over point charges from real atoms, ignoring the ones 
    // created to handle QM-MM bonds.
    for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
        
        BigReal Force = 0;
        
        BigReal pntCharge = pcP[i].charge;
        
        BigReal xMM = pcP[i].position.x;
        BigReal yMM = pcP[i].position.y;
        BigReal zMM = pcP[i].position.z;
        
        for (size_t j=0; j<msg->numQMAtoms; ++j ) {
            
            // Not perfect
            // This prevents the MM point charge of a MM-QM bond from feeling 
            // the influence from the QM atom it is bount to. 
            if ( pcP[i].bountToIndx == j ) continue ;
            
            BigReal qmCharge = atmP[j].charge ;
            
            Force = pntCharge*qmCharge*constants ;
            
            BigReal xQM = atmP[j].position.x;
            BigReal yQM = atmP[j].position.y;
            BigReal zQM = atmP[j].position.z;
            
            BigReal x_ij = (xMM - xQM);
            BigReal y_ij = (yMM - yQM);
            BigReal z_ij = (zMM - zQM);
            
            BigReal r2 = (x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
            BigReal rNorm = sqrt(r2) ;
            
            Force /= r2;
            
            resForce[pcIndx].force.x += Force*x_ij/rNorm;
            resForce[pcIndx].force.y += Force*y_ij/rNorm;
            resForce[pcIndx].force.z += Force*z_ij/rNorm;
        }
        
    }
    
    // Adjusts forces from PME, canceling contributions from the QM and 
    // direct Coulomb forces calculated here.
    if (msg->PMEOn) {
        
        DebugM(1,"Correcting forces and energy for PME.\n");
        
        ewaldcof = msg->PMEEwaldCoefficient;
        BigReal TwoBySqrtPi = 1.12837916709551;
        pi_ewaldcof = TwoBySqrtPi * ewaldcof;
        
        for (size_t i=0; i < msg->numQMAtoms; i++) {
            
            BigReal p_i_charge = atmP[i].charge ;
            Position pos_i = atmP[i].position ;
            
            for (size_t j=i+1; j < msg->numQMAtoms; j++) {
                
                BigReal p_j_charge = atmP[j].charge ;
                
                Position pos_j = atmP[j].position ;
                
                BigReal r = Vector(pos_i - pos_j).length() ;
                
                BigReal tmp_a = r * ewaldcof;
                BigReal tmp_b = erfc(tmp_a);
                BigReal corr_energy = tmp_b;
                BigReal corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
                
//                 BigReal recip_energy = (1-tmp_b)/r = erf(tmp_a)/r;
                BigReal recip_energy = (1-tmp_b)/r;
                
                BigReal recip_gradient = -(1-corr_gradient)/(r*2);
                
                const BigReal kq_i = p_i_charge * constants;
                
                // Final force and energy correction for this pair of atoms.
                BigReal energy = kq_i * p_j_charge * recip_energy ;
                
                Force fixForce = -1*kq_i*p_j_charge*(recip_gradient/r)*(pos_i - pos_j) ;
                
                // The force is *subtracted* from the total force acting on
                // both atoms. The sign on fixForce corrects the orientation
                // of the subtracted force.
//                 DebugM(4,"Old forces for QM " << i << ": " << resForce[i].force
//                     << std::endl);
//                 DebugM(4,"Old forces for QM " << j << ": " << resForce[j].force
//                     << std::endl);
//                 DebugM(4,"Force correction: " << fixForce << std::endl);
                resForce[i].force -= fixForce ;
                resForce[j].force -= -1*fixForce;
                
                // The energy is *subtracted* from the total energy calculated here.
                resMsg->energyCorr -= energy;
            }
        }
        
        pcIndx = msg->numQMAtoms;
        // We only loop over point charges from real atoms, ignoring the ones 
        // created to handle QM-MM bonds.
        for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
            
            BigReal p_i_charge = pcP[i].charge ;
            Position pos_i = pcP[i].position ;
            
            const BigReal kq_i = p_i_charge * constants;
            
            Force fixForce = 0;
            
            for (size_t j=0; j<msg->numQMAtoms; ++j ) {
                
                // Not perfect
                // This prevents the MM point charge of a MM-QM bond from feeling 
                // the influence from the QM atom it is bount to. 
                if ( pcP[i].bountToIndx == j ) continue ;
                
                BigReal p_j_charge = atmP[j].charge ;
                
                Position pos_j = atmP[j].position ;
                
                BigReal r = Vector(pos_i - pos_j).length() ;
                
                BigReal tmp_a = r * ewaldcof;
                BigReal tmp_b = erfc(tmp_a);
                BigReal corr_energy = tmp_b;
                BigReal corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
                
//                 BigReal recip_energy = (1-tmp_b)/r = erf(tmp_a)/r;
                BigReal recip_energy = (1-tmp_b)/r;
                
                BigReal recip_gradient = -(1-corr_gradient)/(r*2);
                
                // Final force and energy correction for this pair of atoms.
                BigReal energy = kq_i * p_j_charge * recip_energy ;
                
                fixForce += -1*p_j_charge*(recip_gradient/r)*(pos_i - pos_j) ;
                
                // The energy is *subtracted* from the total energy calculated here.
                resMsg->energyCorr -= energy;
                
            }
            
            // The force is *subtracted* from the total force acting on
                // the point charge.
//                 DebugM(4,"Old forces for PC " << pcIndx << ": " << resForce[pcIndx].force
//                     << std::endl);
//                 DebugM(4,"Force correction: " << fixForce << std::endl);
            resForce[pcIndx].force -= kq_i*fixForce ;
        }
        
    }
    
    DebugM(1,"Determining virial...\n");
    
    // Calculates the virial contribution form the forces on QM atoms and 
    // point charges calculated here.
    for (size_t i=0; i < msg->numQMAtoms; i++) {
        // virial used by NAMD is -'ve of normal convention, so reverse it!
        // virial[i][j] in file should be sum of -1 * f_i * r_j
        for ( int k=0; k<3; ++k )
            for ( int l=0; l<3; ++l )
                resMsg->virial[k][l] = -1*resForce[i].force[k]*atmP[i].position[l];
    }
    
    pcIndx = msg->numQMAtoms; // Index in the real PC force array.
    for (size_t i=0; i < msg->numRealPntChrgs; i++, pcIndx++ ) {
        // virial used by NAMD is -'ve of normal convention, so reverse it!
        // virial[i][j] in file should be sum of -1 * f_i * r_j
        for ( int k=0; k<3; ++k )
            for ( int l=0; l<3; ++l )
                resMsg->virial[k][l] = -1*resForce[pcIndx].force[k]*pcP[i].position[l];
    }
    
    DebugM(1,"End of QM processing. Sending result message.\n");
    
    // Send message to rank zero with results.
    QMProxy[0].recvQMRes(resMsg);
    
    delete msg;
    return ;
}


void ComputeQMMgr::pntChrgSwitching(QMGrpCalcMsg *msg) {
    
    // We apply a switching function to the point charges so that there is a 
    // smooth decay of the electrostatic environment seen by the QM system.
    
    BigReal cutoff2 = msg->cutoff*msg->cutoff;
    BigReal swdist = msg->swdist;
    
    SortedArray<pntChrgDist> sortedDists;
    
    DebugM(1,"Initiating point charge switching and processing in rank " 
    << CkMyPe() << "\n" ) ;
    
    QMAtomData *pcP = msg->data + msg->numAllAtoms;
    
#ifdef DEBUG_QM
    Real PCScaleCharge = 0;
    for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
        PCScaleCharge += pcP[i].charge;
    }
    DebugM(4, "The initial total Point-Charge charge is " << PCScaleCharge
            << " before scaling type: " << msg->switchType << "\n" );
#endif
    
    switch (msg->switchType) {
        
        case QMPCSCALESHIFT:
        {
            
            // We store all point charges so that only the furthest away
            // are changed in PC schemes "round" or "zero"
            for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
                sortedDists.add( pntChrgDist(i, pcP[i].dist) ) ;
                
                // https://nmr.cit.nih.gov/xplor-nih/xplorMan/node119.html
                // Applying X-PLOR shifting formula:
                // F = Qi*Qj*(C/(epsilon*r))*(1 - r^2/cutoff^2)^2
                Real dist2 = pcP[i].dist*pcP[i].dist ;
                dist2 /= cutoff2 ;
                Real coef = 1- dist2;
                pcP[i].charge *= coef*coef ;
            }
            
        } break;
        
        case QMPCSCALESWITCH:
        {
            for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
                
                // We store the point charges which are beiond the switching
                // distance.
                if (pcP[i].dist > swdist) {
                    sortedDists.add( pntChrgDist(i, pcP[i].dist) ) ;
                    
                    // https://nmr.cit.nih.gov/xplor-nih/xplorMan/node118.html
                    // Applying the XPLOR switching formula:
                    // (r^2 - cutoff^2)^2*(cutoff^2 + 2*r^2 - 3*swdits^2)/(cutoff^2 - swdits^2)^3
                    Real dist2 = pcP[i].dist*pcP[i].dist ;
                    Real swdist2 = swdist*swdist;
                    Real coef = (dist2 - cutoff2)*(dist2 - cutoff2) ;
                    coef *= (cutoff2 + 2*dist2 - 3*swdist2) ;
                    coef /= (cutoff2 - swdist2)*(cutoff2 - swdist2)*(cutoff2 - swdist2);
                    pcP[i].charge *= coef ;
                }
            }
        } break;
    }
    
#ifdef DEBUG_QM
    PCScaleCharge = 0;
    for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
        PCScaleCharge += pcP[i].charge;
    }
    DebugM(4, "The final total Point-Charge charge is " << PCScaleCharge
            << " after scaling.\n" );
#endif
    
    DebugM(4, sortedDists.size() 
    << " point charges were selected for point charge scheme " << msg->pcScheme << "\n" );
    
    Real totalPCCharge = 0, correction = 0;
    switch (msg->pcScheme) {
        
        case QMPCSCHEMEROUND:
        {
            for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
                totalPCCharge += pcP[i].charge;
            }
            DebugM(4, "The total Point-Charge charge is " << totalPCCharge 
            << "\n" );
            
            if ((fabsf(roundf(totalPCCharge) - totalPCCharge) <= 0.001f) ) {
                DebugM(4, "Charge is already a whole number!\n" );
            } 
            else {
                correction = roundf(totalPCCharge) -totalPCCharge ;
                DebugM(4, "Adding to system the charge: " << correction << "\n" );
            }
        } break;
        
        case QMPCSCHEMEZERO:
        {
            for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
                totalPCCharge += pcP[i].charge;
            }
            DebugM(4, "The total Point-Charge charge is " << totalPCCharge << "\n");
            
            DebugM(4, "Total QM system charge is: " << msg->charge << "\n" );
            
            correction = -1*(totalPCCharge + msg->charge);
            if ((fabsf(correction) <= 0.001f) ) {
                correction = 0;
                DebugM(4, "Total QM + PC charge is already zero!\n" );
            }
            else
                DebugM(4, "Adding a charge of " << correction << " to the system\n");
            
        } break;
    }
    
    if (correction != 0) {
        
        int maxi = sortedDists.size(), mini = sortedDists.size()/2;
        Real fraction = correction/(maxi - mini); 
        
        for (size_t i=mini; i<maxi ; i++) {
            
            pcP[sortedDists[i].index].charge += fraction ;
            
        }
        
        #ifdef DEBUG_QM
        totalPCCharge = 0;
        for ( size_t i=0; i<msg->numRealPntChrgs; i++ ) {
            totalPCCharge += pcP[i].charge;
        }
        DebugM(4, "The total Point-Charge charge is " << totalPCCharge 
        << "\n");
        #endif
    }
    
}

void ComputeQMMgr::lssPrepare() {
    
    lssTotRes = 0;
    lssResMass = 0;
    
    DebugM (4, "Preparing LSS for " << numQMGrps << " QM groups.\n" )
    
    for (int i=0; i<numQMGrps; i++) {
        lssTotRes += qmLSSSize[i];
    }
    
    lssPos = new Position[lssTotRes];
    
    grpIDResNum.resize(numQMGrps);
    
    if (simParams->qmLSSMode == QMLSSMODECOM) {
        
        lssGrpRefMass.resize(numQMGrps);
        
        for (int i=0; i<qmLSSResSize; i++)
            lssResMass += qmLSSMass[i];
        
        DebugM(4, "Total residue mass is " << lssResMass << "\n" )
    }
    
    // Get all atom IDs of solvent QM atoms, per group.
    int solvResBeg = 0, refAtmBeg = 0, locIter = 0, solvIndx = 0;
    int *lssIndxs, *refAtmsIndxs ;
    Mass *lssMasses, *refAtmMasses;
    while (locIter < numQMGrps) {
        lssIndxs = qmLSSIdxs + solvResBeg;
        
        DebugM (4, "Loading atom IDs for QM group " << locIter 
            << " with " << qmLSSSize[locIter]
            << " solvent molecules.\n" )
        
        
        switch (simParams->qmLSSMode) {
        
        case QMLSSMODECOM:
        {
            lssMasses = qmLSSMass + solvResBeg;
            
            // Loads data on QM solvent residues and their atoms.
            for (int i=0; i<qmLSSSize[locIter]; i++) {
                
                lssPos[solvIndx] = 0;
                
                Debug( iout << "Getting atom IDs from QM solvent molecule " 
                << solvIndx << "\n") ;
                for (int j=0; j<qmLSSResSize; j++) {
                    
                    int atmID = lssIndxs[i*qmLSSResSize + j];
                    Mass atmMass = lssMasses[i*qmLSSResSize + j];
                    Debug( iout << atmID << " (" << atmMass << ") ") ;
                    
                    grpIDResNum[locIter].insert(atmLSSData(atmID, LSSDataStr(solvIndx,atmMass)));
                    
                }
                solvIndx++;
                Debug( iout << "\n" << endi );
            }
            
            // Loads data on the mass of QM atoms which will be used to calculate 
            // the COM for solvent selection.
            refAtmsIndxs = qmLSSRefIDs + refAtmBeg;
            refAtmMasses = molPtr->get_qmLSSRefMass() + refAtmBeg;
            for (int i=0; i<qmLSSRefSize[locIter]; i++) {
                lssGrpRefMass[locIter].insert( 
                    refLSSData( refAtmsIndxs[i], refAtmMasses[i] ) 
                ) ;
            }
            refAtmBeg += qmLSSRefSize[locIter] ;
        } break ;
        
        case QMLSSMODEDIST:
        {
            // Loads data on QM solvent residues and their atoms.
            for (int i=0; i<qmLSSSize[locIter]; i++) {
                
                lssPos[solvIndx] = 0;
                
                Debug( iout << "Getting atom IDs from QM solvent molecule " 
                << solvIndx << "\n") ;
                for (int j=0; j<qmLSSResSize; j++) {
                    
                    int atmID = lssIndxs[i*qmLSSResSize + j];
                    Debug( iout << atmID << " ") ;
                    
                    grpIDResNum[locIter].insert(atmLSSData(atmID, LSSDataStr(solvIndx,0)));
                    
                }
                solvIndx++;
                Debug( iout << "\n" << endi );
            }
            
        } break ;
        
        }
        
        solvResBeg += qmLSSSize[locIter]*qmLSSResSize ;
        locIter++;
    }
    
    return ;
}

void ComputeQMMgr::lssUpdate(int grpIter, QMAtmVec& grpQMAtmVec, 
                             QMPCVec& grpPntChrgVec) {
    
    SortedArray<lssDistSort> solvDist;
    
    Position refCOM(0) ;
    Mass totMass = 0;
    
    DebugM(3, "LSS UPDATE...\n")
    
    int solvResBeg = 0 ;
    for (int i=0; i<grpIter; i++)
        solvResBeg += qmLSSSize[i] ;
    
    switch (simParams->qmLSSMode ) {
        
    case QMLSSMODECOM:
    {
        DebugM(3, "Using COM for LSS in group " << grpIter << "\n")
        
        // Determined the reference center of mass for this group.
        for(int i=0; i<grpQMAtmVec.size(); i++) {
            
            auto it = lssGrpRefMass[grpIter].find(grpQMAtmVec[i].id);
            if ( it != lssGrpRefMass[grpIter].end() ) {
                refCOM += grpQMAtmVec[i].position*it->second ;
                totMass += it->second ;
            }
        }
        
        refCOM /= totMass;
        DebugM ( 3, "Reference COM position: " << refCOM << "\n");
        
        // Initialize the positions of all COM of quantum residues 
        for (int solvIter=solvResBeg; solvIter<solvResBeg+qmLSSSize[grpIter]; solvIter++) {
            lssPos[solvIter] = 0 ;
        }
        
        DebugM(3, "Calculating distance of QM solvent COM from reference COM of group.\n")
        
        // Temporary handler of lssDistSort structures, while we accumulate atoms
        // on their respective residues and calculate distances.
        std::map<int,lssDistSort > resQMDist ;
        
        // Initiates COM determination for all QM solvent molecules.
        for(int i=0; i<grpQMAtmVec.size(); i++) {
            auto it = grpIDResNum[grpIter].find(grpQMAtmVec[i].id) ;
            if (it != grpIDResNum[grpIter].end()) {
                lssPos[it->second.resIndx] += grpQMAtmVec[i].position*it->second.mass ;
                
                // tries to find a residue number
                auto itRes = resQMDist.find(it->second.resIndx) ;
                if (itRes == resQMDist.end() ) {
                    resQMDist.insert(std::pair<int,lssDistSort >(
                        it->second.resIndx, lssDistSort(QMLSSQMRES, -1) ));
                }
                
                // For each classical residue ID, we compile a list of atom IDs which 
                // are found among point charges
//                 resQMDist[it->second.resIndx].idVec.push_back(grpQMAtmVec[i].id) ;
//                 resQMDist[it->second.resIndx].indxVec.push_back(i) ;
                
                resQMDist[it->second.resIndx].idIndx.add(idIndxStr(grpQMAtmVec[i].id,i)) ;
            }
        }
        
        DebugM(3, "QM Solvent molecules " << solvResBeg 
        << " to " << solvResBeg+qmLSSSize[grpIter] << "\n")
        
        
        for (int solvIter=solvResBeg; solvIter<solvResBeg+qmLSSSize[grpIter]; solvIter++) {
            Real dist = Vector((lssPos[solvIter]/lssResMass) - refCOM).length() ;
            resQMDist[solvIter].dist = dist ;
            solvDist.add(resQMDist[solvIter]);
        }
        
        #ifdef DEBUG_QM
        DebugM(3, "We loaded the following QM solvent residues and distances:\n")
        for (int i=0; i<solvDist.size(); i++) {
            iout << i << ") type: " << solvDist[i].type 
            << " dist " << solvDist[i].dist
            << " IDs: " ;
            for (int j=0; j<solvDist[i].idIndx.size(); j++) 
                iout << solvDist[i].idIndx[j].ID << " " ;
            iout << "\n" << endi;
        }
        #endif
        
        // Get Center Of Mass distance for all (whole) solvent molecules among 
        // *point charges*, comparing against the non-solvent QM atoms.
        
        // This will count how many PCs are associated with each solvent residue,
        // and associate their atom ID with that residue.
        std::map<int,lssDistSort > resPCSize ;
        
        for(int i=0; i<grpPntChrgVec.size(); i++) {
            
            // This maps aomt IDs with a global classical residue ID.
            auto it = molPtr->get_qmMMSolv().find(grpPntChrgVec[i].id) ;
            if (it != molPtr->get_qmMMSolv().end()) {
                
                // tries to find a residue number
                auto itRes = resPCSize.find(it->second) ;
                if (itRes == resPCSize.end() ) {
                    resPCSize.insert(std::pair<int,lssDistSort >(
                        it->second, lssDistSort(QMLSSCLASSICALRES, -1) ));
                }
                
                // For each classical residue ID, we compile a list of atom IDs which 
                // are found among point charges
//                 resPCSize[it->second].idVec.push_back(grpPntChrgVec[i].id) ;
//                 resPCSize[it->second].indxVec.push_back(i) ;
                
                resPCSize[it->second].idIndx.add(idIndxStr(grpPntChrgVec[i].id,i)) ;
            }
        }
        
        // Now we check which classical solvent molecules are complete,
        // and compute the COM for the complete ones, while ignoring the 
        // incomplete ones.
        for (auto it=resPCSize.begin(); it!=resPCSize.end(); it++) {
            
            if (it->second.idIndx.size() == qmLSSResSize) {
                
                Position currCOM(0);
                Mass totalMass = 0;
                
    //             iout << "Found complete classical residue " << it->first << "\n" << endi;
                
                for (int i=0; i<it->second.idIndx.size(); i++) {
    //                 iout << "AtomID " << it->second.idIndx[i] .ID
    //                 << " indxVec " << it->second.idIndx[i].indx
    //                 << " mass " << grpPntChrgVec[it->second.idIndx[i].indx].mass
    //                 << " position " << grpPntChrgVec[it->second.idIndx[i].indx].position << "\n" << endi;
                    currCOM += grpPntChrgVec[it->second.idIndx[i].indx].position*grpPntChrgVec[it->second.idIndx[i].indx].mass;
                    totalMass += grpPntChrgVec[it->second.idIndx[i].indx].mass;
                }
                currCOM /= totalMass;
                
                Real currDist = Vector(currCOM - refCOM).length() ;
                
                it->second.dist = currDist ;
                
                solvDist.add(it->second) ;
            }
            
        }
    } break ;
    
    case QMLSSMODEDIST:
    {
        DebugM(3, "Using minimal distances for LSS in group " << grpIter << "\n")
        
        DebugM(3, "QM Solvent molecules " << solvResBeg 
        << " to " << solvResBeg+qmLSSSize[grpIter] << "\n")
        
        // List of QM indices which will be used as reference for distance calculation.
        ResizeArray<int> qmRefIndx ;
        
        // Temporary handler of lssDistSort structures, while we accumulate atoms
        // on their respective residues and calculate distances.
        std::map<int,lssDistSort > resQMDist ;
        
        // Initiates COM determination for all QM solvent molecules.
        for(int i=0; i<grpQMAtmVec.size(); i++) {
            
            auto it = grpIDResNum[grpIter].find(grpQMAtmVec[i].id) ;
            if (it != grpIDResNum[grpIter].end()) {
                
                // tries to find a residue number
                auto itRes = resQMDist.find(it->second.resIndx) ;
                if (itRes == resQMDist.end() ) {
                    resQMDist.insert(std::pair<int,lssDistSort >(
                        it->second.resIndx, lssDistSort(QMLSSQMRES, -1) ));
                }
                
                // For each classical residue ID, we compile a list of atom IDs which 
                // are found among point charges
//                 resQMDist[it->second.resIndx].idVec.push_back(grpQMAtmVec[i].id) ;
//                 resQMDist[it->second.resIndx].indxVec.push_back(i) ;
                
                resQMDist[it->second.resIndx].idIndx.add(idIndxStr(grpQMAtmVec[i].id,i)) ;
                
            }
            else {
                qmRefIndx.add(i) ;
            }
        }
        
        // We now calculate the shortest distance between a QM solvent residue 
        // and any non-solvent QM atom.
        for (auto it=resQMDist.begin(); it != resQMDist.end(); it++) {
            
            // We prime the residue distance with the first non-solvent
            // QM atom.
            it->second.dist = Vector(
                    grpQMAtmVec[it->second.idIndx[0].indx].position - 
                    grpQMAtmVec[qmRefIndx[0]].position
                    ).length() ;
            
            for (int i=0; i<it->second.idIndx.size(); i++) {
                
                for(int j=0; j<qmRefIndx.size(); j++) {
                    Real currDist = Vector(
                        grpQMAtmVec[it->second.idIndx[i].indx].position - 
                        grpQMAtmVec[qmRefIndx[j]].position
                        ).length() ;
                    
                    if (currDist < it->second.dist)
                        it->second.dist = currDist;
                }
            }
            
            // After determining the distance of this QM solvent residue,
            // we add it to the sorted list.
            solvDist.add(it->second) ;
        }
        
                
        #ifdef DEBUG_QM
        DebugM(3, "We loaded the following QM solvent residues and distances:\n")
        for (int i=0; i<solvDist.size(); i++) {
            iout << i << ") type: " << solvDist[i].type 
            << " dist " << solvDist[i].dist
            << " IDs: " ;
            for (int j=0; j<solvDist[i].idIndx.size(); j++) 
                iout << solvDist[i].idIndx[j].ID << " " ;
            iout << "\n" << endi;
        }
        #endif
        
        // Get shortest distance for all (whole) solvent molecules among 
        // *point charges*, comparing against the non-solvent QM atoms.
        
        // This will count how many PCs are associated with each solvent residue,
        // and associate their atom ID with that residue.
        std::map<int,lssDistSort > resPCSize ;
        
        for(int i=0; i<grpPntChrgVec.size(); i++) {
            
            // This maps aomt IDs with a global classical residue ID.
            auto it = molPtr->get_qmMMSolv().find(grpPntChrgVec[i].id) ;
            if (it != molPtr->get_qmMMSolv().end()) {
                
                // tries to find a residue number
                auto itRes = resPCSize.find(it->second) ;
                if (itRes == resPCSize.end() ) {
                    resPCSize.insert(std::pair<int,lssDistSort >(
                        it->second, lssDistSort(QMLSSCLASSICALRES, -1) ));
                }
                
                // For each classical residue ID, we compile a list of atom IDs which 
                // are found among point charges
//                 resPCSize[it->second].idVec.push_back(grpPntChrgVec[i].id) ;
//                 resPCSize[it->second].indxVec.push_back(i) ;
                
                resPCSize[it->second].idIndx.add(idIndxStr(grpPntChrgVec[i].id,i)) ;
            }
        }
        
        // Now we check which classical solvent molecules are complete,
        // and compute the COM for the complete ones, while ignoring the 
        // incomplete ones.
        for (auto it=resPCSize.begin(); it!=resPCSize.end(); it++) {
            
            if (it->second.idIndx.size() == qmLSSResSize) {
                
                // We prime the residue distance with the first non-solvent
                // QM atom.
                it->second.dist = Vector(
                        grpPntChrgVec[it->second.idIndx[0].indx].position - 
                        grpQMAtmVec[qmRefIndx[0]].position
                        ).length() ;
                
                for (int i=0; i<it->second.idIndx.size(); i++) {
                    
                    for(int j=0; j<qmRefIndx.size(); j++) {
                        Real currDist = Vector(
                            grpPntChrgVec[it->second.idIndx[i].indx].position - 
                            grpQMAtmVec[qmRefIndx[j]].position
                            ).length() ;
                        
                        if (currDist < it->second.dist)
                            it->second.dist = currDist;
                    }
                }
                
                solvDist.add(it->second) ;
            }
            
        }
        
    } break ;
    
    } // End switch
    
    #ifdef DEBUG_QM
    DebugM(3, "Final selection of solvent residues and distances:\n")
    for (int i=0; i<qmLSSSize[grpIter]; i++) {
        std::string typeS ;
        if (solvDist[i].type != QMLSSCLASSICALRES ) 
            typeS = "QM" ;
        else 
            typeS = "Classical";
        iout << i << ") type: " << typeS
        << " dist " << solvDist[i].dist
        << " IDs: " ;
        for (int j=0; j<solvDist[i].idIndx.size(); j++) 
            iout << solvDist[i].idIndx[j].ID << " " ;
        iout << "\n" << endi;
    }
    #endif
    
    // Compare COM distances of QM and Classical solvent molecules, creating a
    // substitution list, atmID by atmID.
    
    DebugM(3, "Determining residues to be swaped...\n")
    
    ResizeArray<lssDistSort> nearMM, farQM ;
    
    for (int resIter = 0; resIter < qmLSSSize[grpIter] ; resIter++) {
        if (solvDist[resIter].type == QMLSSCLASSICALRES) {
            nearMM.add(solvDist[resIter]);
        }
    }
    
    for (int resIter=qmLSSSize[grpIter]; resIter<solvDist.size(); resIter++) {
        if (solvDist[resIter].type == QMLSSQMRES) {
            farQM.add(solvDist[resIter]);
        }
        
        if (farQM.size() == nearMM.size()) break;
    }
    
    if (farQM.size() != nearMM.size())
        NAMD_die("Could not find complementing residues to be swapped in LSS.") ;
    
    #ifdef DEBUG_QM
    DebugM(3, "Removing the following QM residues:\n")
    for (int i=0; i<farQM.size();i++) {
        std::string typeS ;
        if (farQM[i].type != QMLSSCLASSICALRES ) 
            typeS = "QM" ;
        else 
            typeS = "Classical";
        iout << i << ") type: " << typeS
        << " dist " << farQM[i].dist
        << " IDs: " ;
        for (int j=0; j<farQM[i].idIndx.size(); j++) 
            iout << farQM[i].idIndx[j].ID << " " ;
        iout << "\n" << endi;
    }
    
    DebugM(3, "Replacing with the following Classical residues:\n")
    for (int i=0; i<nearMM.size();i++) {
        std::string typeS ;
        if (nearMM[i].type != QMLSSCLASSICALRES ) 
            typeS = "QM" ;
        else 
            typeS = "Classical";
        iout << i << ") type: " << typeS
        << " dist " << nearMM[i].dist
        << " IDs: " ;
        for (int j=0; j<nearMM[i].idIndx.size(); j++) 
            iout << nearMM[i].idIndx[j].ID << " " ;
        iout << "\n" << endi;
    }
    
    DebugM(3, "Building substitution array...\n")
    #endif
    
    iout << iINFO << "LSS is swapping " << farQM.size() << " solvent residues in QM group "
    << grpIter << "\n" << endi;
    
    // Now we build the array which will be sent to all nodes with force results from
    // this step.
    // Atom reassignment will be done in the next step, and will use this data.
    for (int i=0; i<farQM.size();i++) {
        
        for(int j=0; j<qmLSSResSize; j++) {
            
            int qIndx= farQM[i].idIndx[j].indx;
            int mIndx= nearMM[i].idIndx[j].indx;
            
            subsArray.add( LSSSubsDat( grpQMAtmVec[qIndx].id,
                                       grpPntChrgVec[mIndx].id,
                                       grpPntChrgVec[mIndx].vdwType,
                                       grpPntChrgVec[mIndx].charge ) );
            
            subsArray.add( LSSSubsDat( grpPntChrgVec[mIndx].id,
                                       grpQMAtmVec[qIndx].id,
                                       grpQMAtmVec[qIndx].vdwType,
                                       grpQMAtmVec[qIndx].charge ) );
        }
        
    }
    
    #ifdef DEBUG_QM
    for(int i=0; i<subsArray.size() ;i++) {
        DebugM(3, CkMyPe() << ") Storing LSS atom " << subsArray[i].origID
        << " Which will become " << subsArray[i].newID
        << " - " << subsArray[i].newVdWType
        << " - " << subsArray[i].newCharge << "\n" ); 
    }
    #endif
    
    return;
}


#ifdef DEBUG_QM

void ComputeQMMgr::Write_PDB(std::string Filename, const QMGrpCalcMsg *dataMsg)
{
    std::ofstream OutputTmpPDB ;
    
    try
    {
        
        OutputTmpPDB.open( Filename.c_str(), std::fstream::out );
        
        OutputTmpPDB << "REMARK Information used by NAMD to create the QM simulation." << std::endl ;
        OutputTmpPDB << "REMARK Occupancy: bount to index; Beta: System ID; " << std::endl ;
        
        const QMAtomData *dataP = dataMsg->data;
        
        for ( int i=0 ; i < dataMsg->numAllAtoms + dataMsg->numAllPntChrgs ; i++ )
        {
            // PDB format: http://www.wwpdb.org/documentation/format33/sect9.html
            
            // Exception: The field resName was changed from the coluns 18-20 to the columns 18-21
            // This allows the usage of protonated amino acids and other molecules in the PDB file.
            
            // Exception2: The field extraInfo was changed from the coluns 67 - 76 to the columns 67 - 72
            // This allows the usage of segments in PDB files, necessary for the propper usage of NAMD.
            
            std::string name(" x  ");
            std::string resName ("  uk");
            std::string chainName("X");
            std::string element("") ;
            if (i < dataMsg->numAllAtoms ) {
                name  = dataP[i].element;
                chainName = "q" ;
                element = dataP[i].element;
                if (dataP[i].type == QMATOMTYPE_QM)
                    resName = " qm " ;
                else if (dataP[i].type == QMATOMTYPE_DUMMY)
                    resName = " dm " ;
            }
            else {
                chainName = "c" ;
                if (dataP[i].type == QMPCTYPE_CLASSICAL)
                    resName = " pc ";
                else if (dataP[i].type == QMPCTYPE_IGNORE)
                    resName = "ipc ";
                else if (dataP[i].type == QMPCTYPE_EXTRA)
                    resName = "vpc ";
            }
            
            OutputTmpPDB << "ATOM  " ; // ATOM  1 -  6
            
            OutputTmpPDB.width(5) ; // serial  7 - 11
            OutputTmpPDB.right ;
            OutputTmpPDB << i ;
            
            OutputTmpPDB << " " ; // Spacing
            
            OutputTmpPDB.width(4) ; // name  13 - 16
            OutputTmpPDB << name ;
            
            OutputTmpPDB << " " ; // altLoc  17
            
            OutputTmpPDB.width(4) ; // resName  18 - 21
            OutputTmpPDB << resName ;
            
            OutputTmpPDB.width(1) ; // chainID  22
            OutputTmpPDB << chainName ;
            
            OutputTmpPDB.width(4) ; // Residue Index  23 - 26
            OutputTmpPDB << i ;
            
            OutputTmpPDB << " " ; // iCode  27
            
            OutputTmpPDB << "   " ; // Spacing
            
            OutputTmpPDB.width(8) ; // x  31 - 38
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(3) ;
            OutputTmpPDB << dataP[i].position.x ;
            
            OutputTmpPDB.width(8) ; // y  39 - 46
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(3) ;
            OutputTmpPDB << dataP[i].position.y ;
            
            OutputTmpPDB.width(8) ; // z  47 - 54
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(3) ;
            OutputTmpPDB << dataP[i].position.z ;
            
            OutputTmpPDB.width(6) ; // occupancy 55 - 60
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(2) ;
            OutputTmpPDB << dataP[i].bountToIndx ;
            
            OutputTmpPDB.width(6) ; // tempFactor/Beta 61 - 66
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(2) ;
            OutputTmpPDB << dataP[i].id ;
            
            OutputTmpPDB.width(6) ; // extra information not originaly on PDB format  67 - 72
            OutputTmpPDB << "      " ;
            
            OutputTmpPDB.width(4) ; // segment information from NAMD, not originaly on PDB format  72 - 76
            OutputTmpPDB.left ;
            OutputTmpPDB << "QM  ";
            
            OutputTmpPDB.width(2) ; // element 77 - 78
            OutputTmpPDB.right ;
            OutputTmpPDB << element ;
            
            OutputTmpPDB.width(2) ; // charge 77 - 78
            OutputTmpPDB.right ;
            OutputTmpPDB << dataP[i].charge ;
            
            OutputTmpPDB << std::endl;
            
        }
        
        OutputTmpPDB << "END" << std::endl;
        
        OutputTmpPDB.close();
    }
    catch (...)
    {
        iout << iERROR << "Generic exception at QM write PBD." << endi ;
        NAMD_die("PDB write error");
        throw "Generic exception!" ;
    }
    return ;
}

void ComputeQMMgr::Write_PDB(std::string Filename, const QMCoordMsg *dataMsg)
{
    std::ofstream OutputTmpPDB ;
    
    try
    {
        
        OutputTmpPDB.open( Filename.c_str(), std::fstream::out );
        
        OutputTmpPDB << "REMARK Information used by NAMD to create the QM simulation." << std::endl ;
        OutputTmpPDB << "REMARK Occupancy: bount to index; Beta: System ID; " << std::endl ;
        
        const ComputeQMAtom *dataP = dataMsg->coord;
        
        for ( int i=0 ; i < dataMsg->numAtoms; i++ )
        {
            // PDB format: http://www.wwpdb.org/documentation/format33/sect9.html
            
            // Exception: The field resName was changed from the coluns 18-20 to the columns 18-21
            // This allows the usage of protonated amino acids and other molecules in the PDB file.
            
            // Exception2: The field extraInfo was changed from the coluns 67 - 76 to the columns 67 - 72
            // This allows the usage of segments in PDB files, necessary for the propper usage of NAMD.
            
            std::string name(" x  ");
            std::string resName (" atm");
            std::string chainName("X");
            std::string element("") ;
            
            OutputTmpPDB << "ATOM  " ; // ATOM  1 -  6
            
            OutputTmpPDB.width(5) ; // serial  7 - 11
            OutputTmpPDB.right ;
            OutputTmpPDB << i ;
            
            OutputTmpPDB << " " ; // Spacing
            
            OutputTmpPDB.width(4) ; // name  13 - 16
            OutputTmpPDB << name ;
            
            OutputTmpPDB << " " ; // altLoc  17
            
            OutputTmpPDB.width(4) ; // resName  18 - 21
            OutputTmpPDB << resName ;
            
            OutputTmpPDB.width(1) ; // chainID  22
            OutputTmpPDB << chainName ;
            
            OutputTmpPDB.width(4) ; // Residue Index  23 - 26
            OutputTmpPDB << i ;
            
            OutputTmpPDB << " " ; // iCode  27
            
            OutputTmpPDB << "   " ; // Spacing
            
            OutputTmpPDB.width(8) ; // x  31 - 38
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(3) ;
            OutputTmpPDB << dataP[i].position.x ;
            
            OutputTmpPDB.width(8) ; // y  39 - 46
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(3) ;
            OutputTmpPDB << dataP[i].position.y ;
            
            OutputTmpPDB.width(8) ; // z  47 - 54
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(3) ;
            OutputTmpPDB << dataP[i].position.z ;
            
            OutputTmpPDB.width(6) ; // occupancy 55 - 60
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(2) ;
            OutputTmpPDB << dataP[i].qmGrpID ;
            
            OutputTmpPDB.width(6) ; // tempFactor/Beta 61 - 66
            OutputTmpPDB.right ;
            OutputTmpPDB.setf(std::ios::fixed,std::ios::floatfield) ;
            OutputTmpPDB.precision(2) ;
            OutputTmpPDB << dataP[i].id ;
            
            OutputTmpPDB.width(6) ; // extra information not originaly on PDB format  67 - 72
            OutputTmpPDB << "      " ;
            
            OutputTmpPDB.width(4) ; // segment information from NAMD, not originaly on PDB format  72 - 76
            OutputTmpPDB.left ;
            OutputTmpPDB << "QM  ";
            
            OutputTmpPDB.width(2) ; // element 77 - 78
            OutputTmpPDB.right ;
            OutputTmpPDB << element ;
            
            OutputTmpPDB.width(2) ; // charge 77 - 78
            OutputTmpPDB.right ;
            OutputTmpPDB << dataP[i].charge ;
            
            OutputTmpPDB << std::endl;
            
        }
        
        OutputTmpPDB << "END" << std::endl;
        
        OutputTmpPDB.close();
    }
    catch (...)
    {
        iout << iERROR << "Generic exception at QM write PBD." << endi ;
        NAMD_die("PDB write error");
        throw "Generic exception!" ;
    }
    return ;
}

#endif

#include "ComputeQMMgr.def.h"

