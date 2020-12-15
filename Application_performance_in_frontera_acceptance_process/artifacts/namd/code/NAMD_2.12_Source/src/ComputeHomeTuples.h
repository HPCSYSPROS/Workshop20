/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEHOMETUPLES_H
#define COMPUTEHOMETUPLES_H

#include "NamdTypes.h"
#include "common.h"
#include "structures.h"
#include "Compute.h"
#include "HomePatch.h"

#include "Box.h"
#include "OwnerBox.h"
#include "UniqueSet.h"

#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeHomeTuples.h"
#include "PatchMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "Parameters.h"
#include "ReductionMgr.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "Priorities.h"
#include "LdbCoordinator.h"

class TuplePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    Box<Patch,CompAtom> *positionBox;
    Box<Patch,CompAtom> *avgPositionBox;
    Box<Patch,Results> *forceBox;
    CompAtom *x;
    CompAtomExt *xExt;
    CompAtom *x_avg;
    Results *r;
    Force *f;
    Force *af;

    int hash() const { return patchID; }

  TuplePatchElem(PatchID pid = -1) {
    patchID = pid;
    p = NULL;
    positionBox = NULL;
    avgPositionBox = NULL;
    forceBox = NULL;
    x = NULL;
    xExt = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
    af = NULL;
  }

  TuplePatchElem(Patch *p_param, Compute *cid) {
    patchID = p_param->getPatchID();
    p = p_param;
    positionBox = p_param->registerPositionPickup(cid);
    avgPositionBox = p_param->registerAvgPositionPickup(cid);
    forceBox = p_param->registerForceDeposit(cid);
    x = NULL;
    xExt = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
    af = NULL;
  }
    
  ~TuplePatchElem() {};

  int operator==(const TuplePatchElem &elem) const {
    return (elem.patchID == patchID);
  }

  int operator<(const TuplePatchElem &elem) const {
    return (patchID < elem.patchID);
  }
};

typedef UniqueSet<TuplePatchElem> TuplePatchList;
typedef UniqueSetIter<TuplePatchElem> TuplePatchListIter;

class AtomMap;
class ReductionMgr;

#ifdef MEM_OPT_VERSION
template <class T> struct ElemTraits {
  typedef AtomSignature signature;
  static signature* get_sig_pointer(Molecule *mol) { return mol->atomSigPool; }
  static int get_sig_id(const CompAtomExt &a) { return a.sigId; }
};

template <> struct ElemTraits <ExclElem> {
  typedef ExclusionSignature signature;
  static signature* get_sig_pointer(Molecule *mol) { return mol->exclSigPool; }
  static int get_sig_id(const CompAtomExt &a) { return a.exclId; }
};
#endif

template <class T, class S, class P> class ComputeHomeTuples : public Compute {

  protected:
  
    virtual void loadTuples(void) {
      int numTuples;

      #ifdef MEM_OPT_VERSION
      typename ElemTraits<T>::signature *allSigs;      
      #else
      int32 **tuplesByAtom;
      /* const (need to propagate const) */ S *tupleStructs;
      #endif
      
      const P *tupleValues;
      Node *node = Node::Object();

      #ifdef MEM_OPT_VERSION
      allSigs = ElemTraits<T>::get_sig_pointer(node->molecule);
      #else      
      T::getMoleculePointers(node->molecule,
		    &numTuples, &tuplesByAtom, &tupleStructs);      
      #endif
      
      T::getParameterPointers(node->parameters, &tupleValues);

      tupleList.resize(0);

      LocalID aid[T::size];

      const int lesOn = node->simParameters->lesOn;
      Real invLesFactor = lesOn ? 
                          1.0/node->simParameters->lesFactor :
                          1.0;

      // cycle through each patch and gather all tuples
      TuplePatchListIter ai(tuplePatchList);
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
        // CompAtom *atom = (*ai).x;
        Patch *patch = (*ai).p;
        int numAtoms = patch->getNumAtoms();
	CompAtomExt *atomExt = (*ai).xExt; //patch->getCompAtomExtInfo();
    
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {              
           /* cycle through each tuple */
           #ifdef MEM_OPT_VERSION
           typename ElemTraits<T>::signature *thisAtomSig =
                   &allSigs[ElemTraits<T>::get_sig_id(atomExt[j])];
           TupleSignature *allTuples;
           T::getTupleInfo(thisAtomSig, &numTuples, &allTuples);
           for(int k=0; k<numTuples; k++) {
               T t(atomExt[j].id, &allTuples[k], tupleValues);
           #else
           /* get list of all tuples for the atom */
           int32 *curTuple = tuplesByAtom[atomExt[j].id];
           for( ; *curTuple != -1; ++curTuple) {             
             T t(&tupleStructs[*curTuple],tupleValues);
           #endif            
             register int i;
             aid[0] = atomMap->localID(t.atomID[0]);
             int homepatch = aid[0].pid;
             int samepatch = 1;
             int has_les = lesOn && node->molecule->get_fep_type(t.atomID[0]);
             for (i=1; i < T::size; i++) {
	         aid[i] = atomMap->localID(t.atomID[i]);
	         samepatch = samepatch && ( homepatch == aid[i].pid );
                 has_les |= lesOn && node->molecule->get_fep_type(t.atomID[i]);
             }
             if ( samepatch ) continue;
             t.scale = has_les ? invLesFactor : 1;
             for (i=1; i < T::size; i++) {
	         homepatch = patchMap->downstream(homepatch,aid[i].pid);
             }
             if ( homepatch != notUsed && isBasePatch[homepatch] ) {
      	       TuplePatchElem *p;
               for (i=0; i < T::size; i++) {
      	         t.p[i] = p = tuplePatchList.find(TuplePatchElem(aid[i].pid));
      	         if ( ! p ) {
                     #ifdef MEM_OPT_VERSION
                     iout << iWARN << "Tuple with atoms ";
                     #else
      	           iout << iWARN << "Tuple " << *curTuple << " with atoms ";
                     #endif
      	           int erri;
      	           for( erri = 0; erri < T::size; erri++ ) {
      	             iout << t.atomID[erri] << "(" <<  aid[erri].pid << ") ";
      	           }
      	           iout << "missing patch " << aid[i].pid << "\n" << endi;
      	           break;
      	         }
      	         t.localIndex[i] = aid[i].index;
               }
      	       if ( ! p ) continue;
             #ifdef MEM_OPT_VERSION
               //avoid adding Tuples whose atoms are all fixed
               if(node->simParameters->fixedAtomsOn &&
		  !node->simParameters->fixedAtomsForces) {
                 int allfixed = 1;
                 for(i=0; i<T::size; i++){
                   CompAtomExt *one = &(t.p[i]->xExt[aid[i].index]);
                   allfixed = allfixed & one->atomFixed;
                 }
                 if(!allfixed) tupleList.add(t);
               }else{
                 tupleList.add(t);
               }
             #else
               tupleList.add(t);
             #endif               
             }
           }
        }
      }
    }

    int doLoadTuples;
  
  protected:
  
    ResizeArray<T> tupleList;
    TuplePatchList tuplePatchList;
  
    PatchMap *patchMap;
    AtomMap *atomMap;
    SubmitReduction *reduction;
    int accelMDdoDihe;
    SubmitReduction *pressureProfileReduction;
    BigReal *pressureProfileData;
    int pressureProfileSlabs;
    char *isBasePatch;
  
    ComputeHomeTuples(ComputeID c) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      
      SimParameters *params = Node::Object()->simParameters;
      accelMDdoDihe=false;
      if (params->accelMDOn) {
         if (params->accelMDdihe || params->accelMDdual) accelMDdoDihe=true;
      }
      if (params->pressureProfileOn) {
        pressureProfileSlabs = T::pressureProfileSlabs = 
          params->pressureProfileSlabs;
        int n = T::pressureProfileAtomTypes = params->pressureProfileAtomTypes;
        pressureProfileReduction = ReductionMgr::Object()->willSubmit(
          REDUCTIONS_PPROF_BONDED, 3*pressureProfileSlabs*((n*(n+1))/2));
        int numAtomTypePairs = n*n;
        pressureProfileData = new BigReal[3*pressureProfileSlabs*numAtomTypePairs];
      } else {
        pressureProfileReduction = NULL;
        pressureProfileData = NULL;
      }
      doLoadTuples = false;
      isBasePatch = 0;
    }

    ComputeHomeTuples(ComputeID c, PatchIDList &pids) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      SimParameters *params = Node::Object()->simParameters;
      accelMDdoDihe=false;
      if (params->accelMDOn) {
         if (params->accelMDdihe || params->accelMDdual) accelMDdoDihe=true;
      }
      if (params->pressureProfileOn) {
        pressureProfileSlabs = T::pressureProfileSlabs = 
          params->pressureProfileSlabs;
        int n = T::pressureProfileAtomTypes = params->pressureProfileAtomTypes;
        pressureProfileReduction = ReductionMgr::Object()->willSubmit(
          REDUCTIONS_PPROF_BONDED, 3*pressureProfileSlabs*((n*(n+1))/2));
        int numAtomTypePairs = n*n;
        pressureProfileData = new BigReal[3*pressureProfileSlabs*numAtomTypePairs];
      } else {
        pressureProfileReduction = NULL;
        pressureProfileData = NULL;
      }
      doLoadTuples = false;
      int nPatches = patchMap->numPatches();
      isBasePatch = new char[nPatches];
      int i;
      for (i=0; i<nPatches; ++i) { isBasePatch[i] = 0; }
      for (i=0; i<pids.size(); ++i) { isBasePatch[pids[i]] = 1; }
    }

  public:
  
    virtual ~ComputeHomeTuples() {
      delete reduction;
      delete [] isBasePatch;
      delete pressureProfileReduction;
      delete pressureProfileData;
    }

    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    virtual void initialize(void) {
    
      // Start with empty list
      tuplePatchList.clear();
    
      int nPatches = patchMap->numPatches();
      int pid;
      for (pid=0; pid<nPatches; ++pid) {
        if ( isBasePatch[pid] ) {
          Patch *patch = patchMap->patch(pid);
	  tuplePatchList.add(TuplePatchElem(patch, this));
        }
      }
    
      // Gather all proxy patches (neighbors, that is)
      PatchID neighbors[PatchMap::MaxOneOrTwoAway];
    
      for (pid=0; pid<nPatches; ++pid) if ( isBasePatch[pid] ) {
        int numNeighbors = patchMap->upstreamNeighbors(pid,neighbors);
        for ( int i = 0; i < numNeighbors; ++i ) {
          if ( ! tuplePatchList.find(TuplePatchElem(neighbors[i])) ) {
            Patch *patch = patchMap->patch(neighbors[i]);
	    tuplePatchList.add(TuplePatchElem(patch, this));
          }
        }
      }
      setNumPatches(tuplePatchList.size());
      doLoadTuples = true;

      basePriority = COMPUTE_PROXY_PRIORITY;  // no patch dependence
    }

    //======================================================================
    // atomUpdate() - Method is invoked after anytime that atoms have been
    // changed in patches used by this Compute object.
    //======================================================================
    void atomUpdate(void) {
      doLoadTuples = true;
    }

//-------------------------------------------------------------------
// Routine which is called by enqueued work msg.  It wraps
// actualy Force computation with the apparatus needed
// to get access to atom positions, return forces etc.
//-------------------------------------------------------------------
    virtual void doWork(void) {

      LdbCoordinator::Object()->startWork(ldObjHandle);

      // Open Boxes - register that we are using Positions
      // and will be depositing Forces.
      UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->x = ap->positionBox->open();
	ap->xExt = ap->p->getCompAtomExtInfo();
        if ( ap->p->flags.doMolly ) ap->x_avg = ap->avgPositionBox->open();
        ap->r = ap->forceBox->open();
        ap->f = ap->r->f[Results::normal];
        if (accelMDdoDihe) ap->af = ap->r->f[Results::amdf]; // for dihedral-only or dual-boost accelMD
      } 
    
      BigReal reductionData[T::reductionDataSize];
      int tupleCount = 0;
      int numAtomTypes = T::pressureProfileAtomTypes;
      int numAtomTypePairs = numAtomTypes*numAtomTypes;
    
      for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;
      if (pressureProfileData) {
        memset(pressureProfileData, 0, 3*pressureProfileSlabs*numAtomTypePairs*sizeof(BigReal));
        // Silly variable hiding of the previous iterator
        UniqueSetIter<TuplePatchElem> newap(tuplePatchList);
        newap = newap.begin();
        const Lattice &lattice = newap->p->lattice;
        T::pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
        T::pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;
      }

      if ( ! Node::Object()->simParameters->commOnly ) {
      if ( doLoadTuples ) {
        loadTuples();
        doLoadTuples = false;
      }
      // take triplet and pass with tuple info to force eval
      T *al = tupleList.begin();
      const int ntuple = tupleList.size();
      if ( ntuple ) T::computeForce(al, ntuple, reductionData, pressureProfileData);
      tupleCount += ntuple;
      }
 
    LdbCoordinator::Object()->endWork(ldObjHandle);

      T::submitReductionData(reductionData,reduction);
      reduction->item(T::reductionChecksumLabel) += (BigReal)tupleCount;
      reduction->submit();

      if (pressureProfileReduction) {
        // For ease of calculation we stored interactions between types
        // i and j in (ni+j).  For efficiency now we coalesce the
        // cross interactions so that just i<=j are stored.
        const int arraysize = 3*pressureProfileSlabs;
        const BigReal *data = pressureProfileData;
        for (int i=0; i<numAtomTypes; i++) {
          for (int j=0; j<numAtomTypes; j++) {
            int ii=i;
            int jj=j;
            if (ii > jj) { int tmp=ii; ii=jj; jj=tmp; }
            const int reductionOffset = 
              (ii*numAtomTypes - (ii*(ii+1))/2 + jj)*arraysize;
            for (int k=0; k<arraysize; k++) {
              pressureProfileReduction->item(reductionOffset+k) += data[k];
            }
            data += arraysize;
          }
        }
        pressureProfileReduction->submit();
      }
    
      // Close boxes - i.e. signal we are done with Positions and
      // AtomProperties and that we are depositing Forces
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->positionBox->close(&(ap->x));
        if ( ap->p->flags.doMolly ) ap->avgPositionBox->close(&(ap->x_avg));
        ap->forceBox->close(&(ap->r));
      }
    }
};


#endif

