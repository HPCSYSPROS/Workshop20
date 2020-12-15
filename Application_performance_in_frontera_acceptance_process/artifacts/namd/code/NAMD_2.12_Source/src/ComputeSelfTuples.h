/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTESELFTUPLES_H
#define COMPUTESELFTUPLES_H

#include "ComputeHomeTuples.h"
#include "LdbCoordinator.h"
    
template <class T, class S, class P> class ComputeSelfTuples :
	public ComputeHomeTuples<T,S,P> {

  private:
  
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

      this->tupleList.resize(0);

      LocalID aid[T::size];

      const int lesOn = node->simParameters->lesOn;
      Real invLesFactor = lesOn ?
                          1.0/node->simParameters->lesFactor :
                          1.0;

      // cycle through each patch and gather all tuples
      // There should be only one!
      TuplePatchListIter ai(this->tuplePatchList);

      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
    
        // CompAtomExt *atom = (*ai).x;
        Patch *patch = (*ai).p;
        int numAtoms = patch->getNumAtoms();
	CompAtomExt *atomExt = (*ai).xExt; //patch->getCompAtomExtInfo();
    
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {
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
           /* cycle through each tuple */
           for( ; *curTuple != -1; ++curTuple) {
             T t(&tupleStructs[*curTuple],tupleValues);
           #endif
             register int i;
             aid[0] = this->atomMap->localID(t.atomID[0]);
             int homepatch = aid[0].pid;
             int samepatch = 1;
             int has_les = lesOn && node->molecule->get_fep_type(t.atomID[0]);
             for (i=1; i < T::size; i++) {
	         aid[i] = this->atomMap->localID(t.atomID[i]);
	         samepatch = samepatch && ( homepatch == aid[i].pid );
                 has_les |= lesOn && node->molecule->get_fep_type(t.atomID[i]);
             }
             if ( samepatch ) {
               t.scale = has_les ? invLesFactor : 1;
               TuplePatchElem *p;
               p = this->tuplePatchList.find(TuplePatchElem(homepatch));
               for(i=0; i < T::size; i++) {
                   t.p[i] = p;
                   t.localIndex[i] = aid[i].index;
               }
             #ifdef MEM_OPT_VERSION
               //avoid adding Tuples whose atoms are all fixed
	       if(node->simParameters->fixedAtomsOn &&
                  !node->simParameters->fixedAtomsForces) {
   		 int allfixed = 1;
                 for(i=0; i<T::size; i++){
                   CompAtomExt *one = &(p->xExt[aid[i].index]);
                   allfixed = allfixed & one->atomFixed;
                 }
                 if(!allfixed) this->tupleList.add(t);
               }else{
                 this->tupleList.add(t);
               }
             #else
               this->tupleList.add(t);
             #endif               
             }
           }
        }
      }
    }

    PatchID patchID;

  public:

    ComputeSelfTuples(ComputeID c, PatchID p) : ComputeHomeTuples<T,S,P>(c) {
      patchID = p;
    }

    virtual ~ComputeSelfTuples() {
      UniqueSetIter<TuplePatchElem> ap(this->tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->p->unregisterPositionPickup(this,&(ap->positionBox));
        ap->p->unregisterAvgPositionPickup(this,&(ap->avgPositionBox));
        ap->p->unregisterForceDeposit(this,&(ap->forceBox));
      }
    }


    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    virtual void initialize(void) {
    
      // Start with empty list
      this->tuplePatchList.clear();
    
      this->tuplePatchList.add(TuplePatchElem(ComputeHomeTuples<T,S,P>::patchMap->patch(patchID), this));
    
      this->setNumPatches(this->tuplePatchList.size());
      this->doLoadTuples = true;

      int myNode = CkMyPe();
      if ( PatchMap::Object()->node(patchID) != myNode )
      {
        this->basePriority = COMPUTE_PROXY_PRIORITY + PATCH_PRIORITY(patchID);
      }
      else
      {
        this->basePriority = COMPUTE_HOME_PRIORITY + PATCH_PRIORITY(patchID);
      }
    }

    void doWork(void) {
//      LdbCoordinator::Object()->startWork(this->ldObjHandle);

#ifdef TRACE_COMPUTE_OBJECTS
    double traceObjStartTime = CmiWallTimer();
#endif

      ComputeHomeTuples<T,S,P>::doWork();

#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+this->cid, traceObjStartTime, CmiWallTimer());
#endif

//      LdbCoordinator::Object()->endWork(this->ldObjHandle);
    }

};


#endif

