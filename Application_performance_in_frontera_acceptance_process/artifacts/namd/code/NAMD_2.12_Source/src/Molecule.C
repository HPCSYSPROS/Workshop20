/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Molecule.C is compiled twice!
// MOLECULE2_C undefined only compiles first half of file
// MOLECULE2_C defined only compiles second half of file
// This is shameful but it works.  Molecule needs refactoring badly.

/*
   The class Molecule is used to hold all of the structural information
   for a simulation.  This information is read in from a .psf file and
   cross checked with the Parameters object passed in.  All of the structural
   information is then stored in arrays for use.
*/

#include "largefiles.h"  // must be first!

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "InfoStream.h"
#include "Molecule.h"
#include "strlib.h"
#include "MStream.h"
#include "Communicate.h"
#include "Node.h"
#include "ObjectArena.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Hydrogen.h"
#include "UniqueSetIter.h"
#include "charm++.h"
/* BEGIN gf */
#include "ComputeGridForce.h"
#include "GridForceGrid.h"

#include "MGridforceParams.h"
/* END gf */

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include "CompressPsf.h"
#include "ParallelIOMgr.h"
#include <deque>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef GROMACS_PAIR
#define GROMACS_PAIR 1
#endif

#ifndef GROMACS_EXCLUSIONS
#define GROMACS_EXCLUSIONS 1
#endif

using namespace std;

#ifndef MOLECULE2_C  // first object file

#ifdef MEM_OPT_VERSION
template int lookupCstPool<AtomSignature>(const vector<AtomSignature>&, const AtomSignature&);
template int lookupCstPool<ExclusionSignature>(const vector<ExclusionSignature>&, const ExclusionSignature&);
#endif

int ResidueLookupElem::lookup(
	const char *segid, int resid, int *begin, int *end) const {
    const ResidueLookupElem *elem = this;
    int rval = -1;  // error
    while ( elem && strcasecmp(elem->mySegid,segid) ) elem = elem->next;
    if ( elem && (resid >= elem->firstResid) && (resid <= elem->lastResid) ) {
      *begin = elem->atomIndex[resid - elem->firstResid];
      *end = elem->atomIndex[resid - elem->firstResid + 1];
      rval = 0;  // no error
    }
    return rval;
}

ResidueLookupElem* ResidueLookupElem::append(
	const char *segid, int resid, int aid) {
    ResidueLookupElem *rval = this;
    if ( firstResid == -1 ) {  // nothing added yet
      strcpy(mySegid,segid);
      firstResid = resid;
      lastResid = resid;
      atomIndex.add(aid);
      atomIndex.add(aid+1);
    } else if ( ! strcasecmp(mySegid,segid) ) {  // same segid
      if ( resid == lastResid ) {  // same resid
        atomIndex[lastResid - firstResid + 1] = aid + 1;
      } else if ( resid < lastResid ) {  // error
	// We can work around this by creating a new segment.
	iout << iWARN << "Residue " << resid <<
	  " out of order in segment " << segid <<
	  ", lookup for additional residues in this segment disabled.\n" << endi;
	rval = next = new ResidueLookupElem;
	next->append(segid,resid,aid);
      } else {  // new resid
        for ( ; lastResid < resid; ++lastResid ) atomIndex.add(aid);
        atomIndex[lastResid - firstResid + 1] = aid + 1;
      }
    } else {  // new segid
      rval = next = new ResidueLookupElem;
      next->append(segid,resid,aid);
    }
    return rval;
}


//  Lookup atom id from segment, residue, and name
int Molecule::get_atom_from_name(
	const char *segid, int resid, const char *aname) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }

  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return -1;
  for ( ; i < end; ++i ) {
    #ifdef MEM_OPT_VERSION    
    Index idx = atomNames[i].atomnameIdx;
    if(!strcasecmp(aname, atomNamePool[idx])) return i;
    #else
    if ( ! strcasecmp(aname,atomNames[i].atomname) ) return i;
    #endif
  }
  return -1;
}

//  Lookup number of atoms in residue from segment and residue
int Molecule::get_residue_size(
	const char *segid, int resid) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }
  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return 0;
  return ( end - i );
}

//  Lookup atom id from segment, residue, and index in residue
int Molecule::get_atom_from_index_in_residue(
	const char *segid, int resid, int index) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }
  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return -1;
  if ( index >= 0 && index < ( end - i ) ) return ( index + i );
  return -1;
}

/************************************************************************/
/*                  */
/*      FUNCTION initialize  */
/*                  */
/*  This is the initializer for the Molecule class.  It simply sets */
/*  the counts for all the various parameters to 0 and sets the pointers*/
/*  to the arrays that will store these parameters to NULL, since they  */
/*  have not been allocated yet.          */
/*                  */
/************************************************************************/

void Molecule::initialize(SimParameters *simParams, Parameters *param)
{
  if ( sizeof(int32) != 4 ) { NAMD_bug("sizeof(int32) != 4"); }
  this->simParams = simParams;
  this->params = param;

  /*  Initialize array pointers to NULL  */
  atoms=NULL;
  atomNames=NULL;
  resLookup=NULL;

  // DRUDE
  is_lonepairs_psf = 0;
  is_drude_psf = 0;  // assume not Drude model
  drudeConsts=NULL;
  lphosts=NULL;
  anisos=NULL;
  tholes=NULL;
  lphostIndexes=NULL;
  // DRUDE

  //LCPO
  lcpoParamType = NULL;

  //for compressing molecule info
  atomSegResids=NULL;

  if ( simParams->globalForcesOn ) {
    resLookup = new ResidueLookupElem;
  }

  #ifdef MEM_OPT_VERSION
  eachAtomSig = NULL;
  atomSigPoolSize = 0;
  atomSigPool = NULL;
  massPoolSize = 0;
  atomMassPool = NULL;
  eachAtomMass = NULL;
  chargePoolSize = 0;
  atomChargePool = NULL;
  eachAtomCharge = NULL;
  #else
  bonds=NULL;
  angles=NULL;
  dihedrals=NULL;
  impropers=NULL;
  crossterms=NULL;
  #endif

  donors=NULL;
  acceptors=NULL;
  

  #ifndef MEM_OPT_VERSION      
  tmpArena=NULL;
  exclusions=NULL;
  bondsWithAtom=NULL;
  bondsByAtom=NULL;
  anglesByAtom=NULL;
  dihedralsByAtom=NULL;
  impropersByAtom=NULL;
  crosstermsByAtom=NULL;
  // JLai
  gromacsPairByAtom=NULL;
  // End of JLai
  // DRUDE
  tholesByAtom=NULL;
  anisosByAtom=NULL;
  // DRUDE
  #endif

  #ifdef MEM_OPT_VERSION
  exclSigPool = NULL;
  exclChkSigPool = NULL;
  exclSigPoolSize = 0;
  eachAtomExclSig = NULL;

  fixedAtomsSet = NULL;
  constrainedAtomsSet = NULL;
  #else
  exclusionsByAtom=NULL;
  fullExclusionsByAtom=NULL;
  modExclusionsByAtom=NULL;
  all_exclusions=NULL;
  #endif

  langevinParams=NULL;
  fixedAtomFlags=NULL;

  #ifdef MEM_OPT_VERSION  
  clusterSigs=NULL;
  #else
  cluster=NULL;  
  #endif
  clusterSize=NULL;

  exPressureAtomFlags=NULL;
  rigidBondLengths=NULL;
  consIndexes=NULL;
  consParams=NULL;
  /* BEGIN gf */
  gridfrcIndexes=NULL;
  gridfrcParams=NULL;
  gridfrcGrid=NULL;
  numGridforces=NULL;
  /* END gf */
  stirIndexes=NULL;
  stirParams=NULL;
  movDragIndexes=NULL;
  movDragParams=NULL;
  rotDragIndexes=NULL;
  rotDragParams=NULL;
  consTorqueIndexes=NULL;
  consTorqueParams=NULL;
  consForceIndexes=NULL;
  consForce=NULL;
//fepb
  fepAtomFlags=NULL;
//fepe

  nameArena = new ObjectArena<char>;
  // nameArena->setAlignment(8);
  // arena->setAlignment(32);
  #ifndef MEM_OPT_VERSION
  arena = new ObjectArena<int32>;
  exclArena = new ObjectArena<char>;
  #endif
  // exclArena->setAlignment(32);

  /*  Initialize counts to 0 */
  numAtoms=0;
  numRealBonds=0;
  numBonds=0;
  numAngles=0;
  numDihedrals=0;
  numImpropers=0;
  numTholes=0;
  numAnisos=0;
  numCrossterms=0;
  // JLai
  numLJPair=0;
  // End of JLai
  numDonors=0;
  numAcceptors=0;
  numExclusions=0;

  // DRUDE
  numLonepairs=0;
  numDrudeAtoms=0;
  numLphosts=0;
  numAnisos=0;
  // DRUDE

  numConstraints=0;
  numStirredAtoms=0;
  numMovDrag=0;
  numRotDrag=0;
  numConsTorque=0;
  numConsForce=0;
  numFixedAtoms=0;
  numFixedGroups=0;
  numExPressureAtoms=0;
  numRigidBonds=0;
  numFixedRigidBonds=0;
  numMultipleDihedrals=0;
  numMultipleImpropers=0;
  numCalcBonds=0;
  numCalcAngles=0;
  numCalcDihedrals=0;
  numCalcImpropers=0;
  numCalcTholes=0;
  numCalcAnisos=0;
  numCalcCrossterms=0;
  numCalcExclusions=0;
  numCalcFullExclusions=0;
  // JLai
  numCalcLJPair=0;
  // End of JLai

//fepb
  numFepInitial = 0;
  numFepFinal = 0;
//fepe

  //fields related with pluginIO-based loading molecule structure
  occupancy = NULL;
  bfactor = NULL;

  qmElementArray=0;
  qmDummyElement=0;
  qmGrpSizes=0;
  qmAtomGroup=0;
  qmAtmChrg=0;
  qmAtmIndx=0;
  qmGrpID=0;
  qmGrpChrg=0;
  qmGrpMult=0;
  qmGrpNumBonds=0;
  qmMMBond=0;
  qmGrpBonds=0;
  qmMMBondedIndx=0;
  qmMMChargeTarget=0;
  qmMMNumTargs=0;
  qmDummyBondVal=0;
  qmMeMMindx=0;
  qmMeQMGrp=0;
  qmCustomPCIdxs=0;
  qmCustPCSizes=0;
  qmLSSSize=0;
  qmLSSIdxs=0;
  qmLSSMass=0;
  qmLSSRefIDs=0;
  qmLSSRefMass=0;
  qmLSSRefSize=0;
  qmNumBonds=0;
  
  goInit();
}

/*      END OF FUNCTION initialize */

/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for the Molecule class. */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param)
{
  initialize(simParams,param);
}

/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for the Molecule class from CHARMM/XPLOR files. */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, char *filename, ConfigList *cfgList)
{
  initialize(simParams,param);

#ifdef MEM_OPT_VERSION
  if(simParams->useCompressedPsf)
      read_mol_signatures(filename, param, cfgList);
#else
	read_psf_file(filename, param);	
 //LCPO
  if (simParams->LCPOOn)
    assignLCPOTypes( 0 );
#endif      
 }

/************************************************************************/
/*                                                                      */
/*      FUNCTION Molecule                                               */
/*                                                                      */
/*  This is the constructor for the Molecule class from plugin IO.      */
/*                                                                      */
/************************************************************************/
Molecule::Molecule(SimParameters *simParams, Parameters *param, molfile_plugin_t *pIOHdl, void *pIOFileHdl, int natoms)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Sorry, plugin IO is not supported in the memory optimized version.");
#else
    initialize(simParams, param);
    numAtoms = natoms;
    int optflags = MOLFILE_BADOPTIONS;
    molfile_atom_t *atomarray = (molfile_atom_t *) malloc(natoms*sizeof(molfile_atom_t));
    memset(atomarray, 0, natoms*sizeof(molfile_atom_t));

    //1a. read basic atoms information
    int rc = pIOHdl->read_structure(pIOFileHdl, &optflags, atomarray);
    if (rc != MOLFILE_SUCCESS && rc != MOLFILE_NOSTRUCTUREDATA) {
        free(atomarray);
        NAMD_die("ERROR: plugin failed reading structure data");
    }
    if(optflags == MOLFILE_BADOPTIONS) {
        free(atomarray);
        NAMD_die("ERROR: plugin didn't initialize optional data flags");
    }
    if(optflags & MOLFILE_OCCUPANCY) {
        setOccupancyData(atomarray);
    }
    if(optflags & MOLFILE_BFACTOR) {
        setBFactorData(atomarray);
    }
    //1b. load basic atoms information to the molecule object
    plgLoadAtomBasics(atomarray);    
    free(atomarray);

    //2a. read bonds
    //indices are one-based in read_bonds
    int *from, *to;
    float *bondorder;
    int *bondtype, nbondtypes;
    char **bondtypename;
    if(pIOHdl->read_bonds!=NULL) {
        if(pIOHdl->read_bonds(pIOFileHdl, &numBonds, &from, &to, &bondorder,
                                 &bondtype, &nbondtypes, &bondtypename)){
            NAMD_die("ERROR: failed reading bond information.");
        }
    }    
    //2b. load bonds information to the molecule object
    if(numBonds!=0) {
        plgLoadBonds(from,to);
    }

    //3a. read other bonded structures
    int *plgAngles, *plgDihedrals, *plgImpropers, *plgCterms;
    int ctermcols, ctermrows;
    int *angletypes, numangletypes, *dihedraltypes, numdihedraltypes;
    int *impropertypes, numimpropertypes; 
    char **angletypenames, **dihedraltypenames, **impropertypenames;

    plgAngles=plgDihedrals=plgImpropers=plgCterms=NULL;
    if(pIOHdl->read_angles!=NULL) {
        if(pIOHdl->read_angles(pIOFileHdl,
                  &numAngles, &plgAngles,
                  &angletypes, &numangletypes, &angletypenames,
                  &numDihedrals, &plgDihedrals,
                  &dihedraltypes, &numdihedraltypes, &dihedraltypenames,
                  &numImpropers, &plgImpropers,
                  &impropertypes, &numimpropertypes, &impropertypenames,
                  &numCrossterms, &plgCterms, &ctermcols, &ctermrows)) {
            NAMD_die("ERROR: failed reading angle information.");
        }
    }
    //3b. load other bonded structures to the molecule object
    if(numAngles!=0) plgLoadAngles(plgAngles);
    if(numDihedrals!=0) plgLoadDihedrals(plgDihedrals);
    if(numImpropers!=0) plgLoadImpropers(plgImpropers);
    if(numCrossterms!=0) plgLoadCrossterms(plgCterms);

  numRealBonds = numBonds;
  build_atom_status();
  //LCPO
  if (simParams->LCPOOn)
    assignLCPOTypes( 2 );
#endif
}

/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                  */
/*        FUNCTION Molecule      */
/*                  */
/*  This is the destructor for the class Molecule.  It simply frees */
/*  the memory allocated for each of the arrays used to store the       */
/*  structure information.            */
/*                  */
/************************************************************************/

Molecule::~Molecule()
{
  /*  Check to see if each array was ever allocated.  If it was   */
  /*  then free it            */
  if (atoms != NULL)
    delete [] atoms;

  if (atomNames != NULL)
  {
    // subarrarys allocated from arena - automatically deleted
    delete [] atomNames;
  }
  delete nameArena;

  if (resLookup != NULL)
    delete resLookup;

  // DRUDE: free arrays read from PSF
  if (drudeConsts != NULL) delete [] drudeConsts;
  if (lphosts != NULL) delete [] lphosts;
  if (anisos != NULL) delete [] anisos;
  if (tholes != NULL) delete [] tholes;
  if (lphostIndexes != NULL) delete [] lphostIndexes;
  // DRUDE

  //LCPO
  if (lcpoParamType != NULL) delete [] lcpoParamType;

  #ifdef MEM_OPT_VERSION
  if(eachAtomSig) delete [] eachAtomSig;
  if(atomSigPool) delete [] atomSigPool;
  #else
  if (bonds != NULL)
    delete [] bonds;

  if (angles != NULL)
    delete [] angles;

  if (dihedrals != NULL)
    delete [] dihedrals;

  if (impropers != NULL)
    delete [] impropers;

  if (crossterms != NULL)
    delete [] crossterms;

  if (exclusions != NULL)
    delete [] exclusions;
  #endif

  if (donors != NULL)
    delete [] donors;

  if (acceptors != NULL)
    delete [] acceptors;  

  #ifdef MEM_OPT_VERSION
  if(exclSigPool) delete [] exclSigPool;
  if(exclChkSigPool) delete [] exclChkSigPool;
  if(eachAtomExclSig) delete [] eachAtomExclSig;
  if(fixedAtomsSet) delete fixedAtomsSet;
  if(constrainedAtomsSet) delete constrainedAtomsSet;
  #else
  if (bondsByAtom != NULL)
       delete [] bondsByAtom;
  
  if (anglesByAtom != NULL)
       delete [] anglesByAtom;
  
  if (dihedralsByAtom != NULL)
       delete [] dihedralsByAtom;
  
  if (impropersByAtom != NULL)
       delete [] impropersByAtom;
  
  if (crosstermsByAtom != NULL)
       delete [] crosstermsByAtom;  

  if (exclusionsByAtom != NULL)
       delete [] exclusionsByAtom;
  
  if (fullExclusionsByAtom != NULL)
       delete [] fullExclusionsByAtom;
  
  if (modExclusionsByAtom != NULL)
       delete [] modExclusionsByAtom;
  
  if (all_exclusions != NULL)
       delete [] all_exclusions;

  // JLai
  if (gromacsPairByAtom != NULL)
      delete [] gromacsPairByAtom;
  // End of JLai

  // DRUDE
  if (tholesByAtom != NULL)
       delete [] tholesByAtom;
  if (anisosByAtom != NULL)
       delete [] anisosByAtom;
  // DRUDE
  #endif

  //LCPO
  if (lcpoParamType != NULL)
    delete [] lcpoParamType;

  if (fixedAtomFlags != NULL)
       delete [] fixedAtomFlags;

  if (stirIndexes != NULL)
    delete [] stirIndexes;


  #ifdef MEM_OPT_VERSION
  if(clusterSigs != NULL){      
      delete [] clusterSigs;
  }  
  #else
  if (cluster != NULL)
       delete [] cluster;  
  #endif
  if (clusterSize != NULL)
       delete [] clusterSize;

  if (exPressureAtomFlags != NULL)
       delete [] exPressureAtomFlags;

  if (rigidBondLengths != NULL)
       delete [] rigidBondLengths;

//fepb
  if (fepAtomFlags != NULL)
       delete [] fepAtomFlags;
//fepe

  if (qmAtomGroup != NULL)
       delete [] qmAtomGroup;
  
  if (qmAtmIndx != NULL)
       delete [] qmAtmIndx;
  
  if (qmAtmChrg != NULL)
       delete [] qmAtmChrg;
  
  
  if (qmGrpNumBonds != NULL)
       delete [] qmGrpNumBonds;
  
  if (qmGrpSizes != NULL)
       delete [] qmGrpSizes;
  
  if (qmDummyBondVal != NULL)
       delete [] qmDummyBondVal;
  
  if (qmMMNumTargs != NULL)
       delete [] qmMMNumTargs;
  
  if (qmGrpID != NULL)
       delete [] qmGrpID;
  
  if (qmGrpChrg != NULL)
       delete [] qmGrpChrg;
  
  if (qmGrpMult != NULL)
       delete [] qmGrpMult;
  
  if (qmMeMMindx != NULL)
       delete [] qmMeMMindx;
  
  if (qmMeQMGrp != NULL)
       delete [] qmMeQMGrp;
  
  if (qmLSSSize != NULL)
       delete [] qmLSSSize;
  
  if (qmLSSIdxs != NULL)
       delete [] qmLSSIdxs;
  
  if (qmLSSMass != NULL)
       delete [] qmLSSMass;
  
  if (qmLSSRefSize != NULL)
       delete [] qmLSSRefSize;
  
  if (qmLSSRefIDs != NULL)
       delete [] qmLSSRefIDs;
  
  if (qmLSSRefMass != NULL)
       delete [] qmLSSRefMass;
  
  if (qmMMBond != NULL) {
      for(int grpIndx = 0 ; grpIndx < qmNumBonds; grpIndx++) {
          if (qmMMBond[grpIndx] != NULL)
              delete [] qmMMBond[grpIndx];
      }
      delete [] qmMMBond;
  }
  
  if (qmGrpBonds != NULL) {
      for(int grpIndx = 0 ; grpIndx < qmNumGrps; grpIndx++) {
          if (qmGrpBonds[grpIndx] != NULL)
              delete [] qmGrpBonds[grpIndx];
      }
      delete [] qmGrpBonds;
  }
  
  if (qmMMBondedIndx != NULL) {
      for(int grpIndx = 0 ; grpIndx < qmNumGrps; grpIndx++) {
          if (qmMMBondedIndx[grpIndx] != NULL)
              delete [] qmMMBondedIndx[grpIndx];
      }
      delete [] qmMMBondedIndx;
  }
  
  if (qmMMChargeTarget != NULL) {
      for(int grpIndx = 0 ; grpIndx < qmNumBonds; grpIndx++) {
          if (qmMMChargeTarget[grpIndx] != NULL)
              delete [] qmMMChargeTarget[grpIndx];
      }
      delete [] qmMMChargeTarget;
  }
  
    if (qmElementArray != NULL){
        for(int atmInd = 0 ; atmInd < numAtoms; atmInd++) {
            if (qmElementArray[atmInd] != NULL)
                delete [] qmElementArray[atmInd];
        }
        delete [] qmElementArray;
    }
    
    if (qmDummyElement != NULL){
        for(int atmInd = 0 ; atmInd < numAtoms; atmInd++) {
            if (qmDummyElement[atmInd] != NULL)
                delete [] qmDummyElement[atmInd];
        }
        delete [] qmDummyElement;
    }

    if (qmCustPCSizes != NULL){
        delete [] qmCustPCSizes;
    }
    
    if (qmCustomPCIdxs != NULL){
        delete [] qmCustomPCIdxs;
    }

  #ifndef MEM_OPT_VERSION
  delete arena;
  delete exclArena;
  #endif
}
/*      END OF FUNCTION Molecule      */

#ifndef MEM_OPT_VERSION

//===Non-memory optimized version of functions that read Molecule file===//

/************************************************************************/
/*                  */
/*        FUNCTION read_psf_file      */
/*                  */
/*   INPUTS:                */
/*  fname - Name of the .psf file to read        */
/*  params - pointer to Parameters object to use to obtain          */
/*     parameters for vdWs, bonds, etc.      */
/*                  */
/*  This function reads a .psf file in.  This is where just about   */
/*   all of the structural information for this class comes from.  The  */
/*   .psf file contains descriptions of the atom, bonds, angles,        */
/*   dihedrals, impropers, and exclusions.  The parameter object is     */
/*   used to look up parameters for each of these entities.    */
/*                  */
/************************************************************************/

void Molecule::read_psf_file(char *fname, Parameters *params)
{
  char err_msg[512];  //  Error message for NAMD_die
  char buffer[512];  //  Buffer for file reading
  int i;      //  Loop counter
  int NumTitle;    //  Number of Title lines in .psf file
  FILE *psf_file;    //  pointer to .psf file
  int ret_code;    //  ret_code from NAMD_read_line calls

  /* Try and open the .psf file           */
  if ( (psf_file = Fopen(fname, "r")) == NULL)
  {
    sprintf(err_msg, "UNABLE TO OPEN .psf FILE %s", fname);
    NAMD_die(err_msg);
  }

  /*  Read till we have the first non-blank line of file    */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to see if we dropped out of the loop because of a     */
  /*  read error.  This shouldn't happen unless the file is empty */
  if (ret_code!=0)
  {
    sprintf(err_msg, "EMPTY .psf FILE %s", fname);
    NAMD_die(err_msg);
  }

  /*  The first non-blank line should contain the word "psf".    */
  /*  If we can't find it, die.               */
  if (!NAMD_find_word(buffer, "psf"))
  {
    sprintf(err_msg, "UNABLE TO FIND \"PSF\" STRING IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  // DRUDE: set flag if we discover Drude PSF
  if (NAMD_find_word(buffer, "drude"))
  {
    if ( ! simParams->drudeOn ) {
      iout << iWARN << "Reading PSF supporting DRUDE without "
        "enabling the Drude model in the simulation config file\n" << endi;
    }
    is_drude_psf = 1;
    is_lonepairs_psf = 1;
  }
  else if (simParams->lonepairs) {
    is_lonepairs_psf = 1;
  }
  // DRUDE

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to see if we dropped out of the loop because of a     */
  /*  read error.  This shouldn't happen unless there is nothing  */
  /*  but the PSF line in the file        */
  if (ret_code!=0)
  {
    sprintf(err_msg, "MISSING EVERYTHING BUT PSF FROM %s", fname);
    NAMD_die(err_msg);
  }

  /*  This line should have the word "NTITLE" in it specifying    */
  /*  how many title lines there are        */
  if (!NAMD_find_word(buffer, "NTITLE"))
  {
    sprintf(err_msg,"CAN NOT FIND \"NTITLE\" STRING IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  sscanf(buffer, "%d", &NumTitle);

  /*  Now skip the next NTITLE non-blank lines and then read in the*/
  /*  line which should contain NATOM        */
  i=0;

  while ( ((ret_code=NAMD_read_line(psf_file, buffer)) == 0) && 
    (i<NumTitle) )
  {
    if (!NAMD_blank_string(buffer))
      i++;
  }

  /*  Make sure we didn't exit because of a read error    */
  if (ret_code!=0)
  {
    sprintf(err_msg, "FOUND EOF INSTEAD OF NATOM IN PSF FILE %s", 
       fname);
    NAMD_die(err_msg);
  }

  while (NAMD_blank_string(buffer))
  {
    NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we have the line we want      */
  if (!NAMD_find_word(buffer, "NATOM"))
  {
    sprintf(err_msg, "DIDN'T FIND \"NATOM\" IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  /*  Read in the number of atoms, and then the atoms themselves  */
  sscanf(buffer, "%d", &numAtoms);

  read_atoms(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NBONDS IN PSF");
  }

  /*  Look for the string "NBOND"          */
  if (!NAMD_find_word(buffer, "NBOND"))
  {
    NAMD_die("DID NOT FIND NBOND AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of bonds and then the bonds themselves  */
  sscanf(buffer, "%d", &numBonds);

  if (numBonds)
    read_bonds(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NTHETA IN PSF");
  }

  /*  Look for the string "NTHETA"        */
  if (!NAMD_find_word(buffer, "NTHETA"))
  {
    NAMD_die("DID NOT FIND NTHETA AFTER BOND LIST IN PSF");
  }

  /*  Read in the number of angles and then the angles themselves */
  sscanf(buffer, "%d", &numAngles);

  if (numAngles)
    read_angles(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NPHI IN PSF");
  }

  /*  Look for the string "NPHI"          */
  if (!NAMD_find_word(buffer, "NPHI"))
  {
    NAMD_die("DID NOT FIND NPHI AFTER ANGLE LIST IN PSF");
  }

  /*  Read in the number of dihedrals and then the dihedrals      */
  sscanf(buffer, "%d", &numDihedrals);

  if (numDihedrals)
    read_dihedrals(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NIMPHI IN PSF");
  }

  /*  Look for the string "NIMPHI"        */
  if (!NAMD_find_word(buffer, "NIMPHI"))
  {
    NAMD_die("DID NOT FIND NIMPHI AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of Impropers and then the impropers  */
  sscanf(buffer, "%d", &numImpropers);

  if (numImpropers)
    read_impropers(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NDON IN PSF");
  }

  /*  Look for the string "NDON"        */
  if (!NAMD_find_word(buffer, "NDON"))
  {
    NAMD_die("DID NOT FIND NDON AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of hydrogen bond donors and then the donors */
  sscanf(buffer, "%d", &numDonors);

  if (numDonors)
    read_donors(psf_file);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NACC IN PSF");
  }

  /*  Look for the string "NACC"        */
  if (!NAMD_find_word(buffer, "NACC"))
  {
    NAMD_die("DID NOT FIND NACC AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of hydrogen bond donors and then the donors */
  sscanf(buffer, "%d", &numAcceptors);

  if (numAcceptors)
    read_acceptors(psf_file);

  /*  look for the explicit non-bonded exclusion section.     */
  while (!NAMD_find_word(buffer, "NNB"))
  {
    ret_code = NAMD_read_line(psf_file, buffer);

    if (ret_code != 0)
    {
      NAMD_die("EOF ENCOUNTERED LOOKING FOR NNB IN PSF FILE");
    }
  }

  /*  Read in the number of exclusions and then the exclusions    */
  sscanf(buffer, "%d", &numExclusions);

  if (numExclusions)
    read_exclusions(psf_file);

  // DRUDE: read lone pair hosts and anisotropic terms from PSF
  if (is_lonepairs_psf)
  {
    while (!NAMD_find_word(buffer, "NUMLP"))
    {
      ret_code = NAMD_read_line(psf_file, buffer);
      if (ret_code != 0)
      {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NUMLP IN DRUDE PSF FILE");
      }
    }
    sscanf(buffer, "%d", &numLphosts);
    if (numLphosts) read_lphosts(psf_file);
  }

  if (is_drude_psf)
  {
    while (!NAMD_find_word(buffer, "NUMANISO"))
    {
      ret_code = NAMD_read_line(psf_file, buffer);
      if (ret_code != 0)
      {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NUMANISO IN DRUDE PSF FILE");
      }
    }
    sscanf(buffer, "%d", &numAnisos);
    if (numAnisos) read_anisos(psf_file);
  }
  // DRUDE

  /*  look for the cross-term section.     */
  int crossterms_present = 1;
  while (!NAMD_find_word(buffer, "NCRTERM"))
  {
    ret_code = NAMD_read_line(psf_file, buffer);

    if (ret_code != 0)
    {
      // hit EOF before finding cross-term section
      crossterms_present = 0;
      break;
    }
  }

  if ( crossterms_present) {

    /*  Read in the number of cross-terms and then the cross-terms*/
    sscanf(buffer, "%d", &numCrossterms);

    if (numCrossterms)
      read_crossterms(psf_file, params);

  }

  /*  Close the .psf file.  */
  Fclose(psf_file);

  //  analyze the data and find the status of each atom
  numRealBonds = numBonds;
  build_atom_status();
  return;
}

/************************************************************************/
/*                  */
/*        FUNCTION read_atoms      */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*  params - Parameters object to use for parameters    */
/*                  */
/*  this function reads in the Atoms section of the .psf file.      */
/*   This section consists of numAtoms lines that are of the form:  */
/*     <atom#> <mol> <seg#> <res> <atomname> <atomtype> <charge> <mass> */
/*   Each line is read into the appropriate entry in the atoms array.   */
/*   The parameters object is then used to determine the vdW constants  */
/*   for this atom.              */
/*                  */
/************************************************************************/

void Molecule::read_atoms(FILE *fd, Parameters *params)

{
  char buffer[512];  // Buffer for reading from file
  int atom_number=0;  // Atom number 
  int last_atom_number=0; // Last atom number, used to assure
        // atoms are in order
  char segment_name[11]; // Segment name
  char residue_number[11]; // Residue number
  char residue_name[11];  // Residue name
  char atom_name[11];  // Atom name
  char atom_type[11];  // Atom type
  Real charge;    // Charge for the current atom
  Real mass;    // Mass for the current atom
  int read_count;    // Number of fields read by sscanf

  /*  Allocate the atom arrays          */
  atoms     = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];
  if(simParams->genCompressedPsf) {
      atomSegResids = new AtomSegResInfo[numAtoms];
  }

  // DRUDE: supplement Atom data
  if (is_drude_psf) {
    drudeConsts = new DrudeConst[numAtoms];
  }
  // DRUDE

  hydrogenGroup.resize(0);

  if (atoms == NULL || atomNames == NULL )
  {
    NAMD_die("memory allocation failed in Molecule::read_atoms");
  }

  ResidueLookupElem *tmpResLookup = resLookup;

  /*  Loop and read in numAtoms atom lines.      */
  while (atom_number < numAtoms)
  {
    // Standard PSF format has 8 columns:
    // ATOMNUM  SEGNAME  RESIDUE  RESNAME  ATOMNAME  ATOMTYPE  CHARGE  MASS

    /*  Get the line from the file        */
    NAMD_read_line(fd, buffer);

    /*  If its blank or a comment, skip it      */
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') )
      continue;

    /*  Parse up the line          */
    read_count=sscanf(buffer, "%d %s %s %s %s %s %f %f",
       &atom_number, segment_name, residue_number,
       residue_name, atom_name, atom_type, &charge, &mass);

    /*  Check to make sure we found what we were expecting  */
    if (read_count != 8)
    {
      char err_msg[128];

      sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE IN ATOM LINE %d\nLINE=%s",
         last_atom_number+1, buffer);
      NAMD_die(err_msg);
    }

    // DRUDE: read alpha and thole parameters from atom line
    if (is_drude_psf)
    {
      // Drude model PSF format has 11 columns, the 8 above plus 3 more:
      //   (unknown integer)  ALPHA  THOLE
      // These constants are used for the Thole interactions
      // (dipole interactions occurring between excluded non-bonded terms).

      Real alpha, thole;
      read_count=sscanf(buffer,
//          "%*d %*s %*s %*s %*s %*s %*f %*f %*d %*f %*f %f %f", &alpha, &thole);
                // the two columns preceding alpha and thole will disappear
          "%*d %*s %*s %*s %*s %*s %*f %*f %*d %f %f", &alpha, &thole);
      if (read_count != 2)
      {
        char err_msg[128];

        sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE "
            "IN ATOM LINE %d\nLINE=%s", last_atom_number+1, buffer);
        NAMD_die(err_msg);
      }
      drudeConsts[atom_number-1].alpha = alpha;
      drudeConsts[atom_number-1].thole = thole;
    }
    // DRUDE

    /*  Check if this is in XPLOR format  */
    int atom_type_num;
    if ( sscanf(atom_type, "%d", &atom_type_num) > 0 )
    {
      NAMD_die("Structure (psf) file is either in CHARMM format (with numbers for atoms types, the X-PLOR format using names is required) or the segment name field is empty.");
    }

    /*  Make sure the atoms were in sequence    */
    if (atom_number != last_atom_number+1)
    {
      char err_msg[128];

      sprintf(err_msg, "ATOM NUMBERS OUT OF ORDER AT ATOM #%d OF PSF FILE",
         last_atom_number+1);
      NAMD_die(err_msg);
    }

    last_atom_number++;

    /*  Dynamically allocate strings for atom name, atom    */
    /*  type, etc so that we only allocate as much space    */
    /*  for these strings as we really need      */
    int reslength = strlen(residue_name)+1;
    int namelength = strlen(atom_name)+1;
    int typelength = strlen(atom_type)+1;

    atomNames[atom_number-1].resname = nameArena->getNewArray(reslength);
    atomNames[atom_number-1].atomname = nameArena->getNewArray(namelength);
    atomNames[atom_number-1].atomtype = nameArena->getNewArray(typelength);
  
    if (atomNames[atom_number-1].resname == NULL)
    {
      NAMD_die("memory allocation failed in Molecule::read_atoms");
    }

    /*  Put the values from this atom into the atoms array  */
    strcpy(atomNames[atom_number-1].resname, residue_name);
    strcpy(atomNames[atom_number-1].atomname, atom_name);
    strcpy(atomNames[atom_number-1].atomtype, atom_type);
    atoms[atom_number-1].mass = mass;
    atoms[atom_number-1].charge = charge;
    atoms[atom_number-1].status = UnknownAtom;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append(segment_name, atoi(residue_number), atom_number-1);

    if(atomSegResids) { //for compressing molecule information
        AtomSegResInfo *one = atomSegResids + (atom_number - 1);
        memcpy(one->segname, segment_name, strlen(segment_name)+1);
        one->resid = atoi(residue_number);
    }

    /*  Determine the type of the atom (H or O) */
    if ( simParams->ignoreMass ) {
    } else if (atoms[atom_number-1].mass <= 0.05) {
      atoms[atom_number-1].status |= LonepairAtom;
    } else if (atoms[atom_number-1].mass < 1.0) {
      atoms[atom_number-1].status |= DrudeAtom;
    } else if (atoms[atom_number-1].mass <=3.5) {
      atoms[atom_number-1].status |= HydrogenAtom;
    } else if ((atomNames[atom_number-1].atomname[0] == 'O') && 
         (atoms[atom_number-1].mass >= 14.0) && 
         (atoms[atom_number-1].mass <= 18.0)) {
      atoms[atom_number-1].status |= OxygenAtom;
    }

    /*  Look up the vdw constants for this atom    */
    params->assign_vdw_index(atomNames[atom_number-1].atomtype, 
       &(atoms[atom_number-1]));
        }

  return;
}
/*      END OF FUNCTION read_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION read_bonds        */
/*                  */
/*  read_bonds reads in the bond section of the .psf file.  This    */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are bonded together.  Each atom pair is   */
/*  read in.  Then that parameter object is queried to determine the    */
/*  force constant and rest distance for the bond.      */
/*                  */
/************************************************************************/

void Molecule::read_bonds(FILE *fd, Parameters *params)

{
  int atom_nums[2];  // Atom indexes for the bonded atoms
  char atom1name[11];  // Atom type for atom #1
  char atom2name[11];  // Atom type for atom #2
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
  int origNumBonds = numBonds;   // number of bonds in file header

  /*  Allocate the array to hold the bonds      */
  bonds=new Bond[numBonds];

  if (bonds == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_bonds");
  }

  /*  Loop through and read in all the bonds      */
  while (num_read < numBonds)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "BONDS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "BOND INDEX %d GREATER THAN NATOM %d IN BOND # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }
    }

    /*  Get the atom type for the two atoms.  When we query */
    /*  the parameter object, we need to send the atom type */
    /*  that is alphabetically first as atom 1.    */
    if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
         atomNames[atom_nums[1]].atomtype) < 0)
    {
      strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
      strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    }
    else
    {
      strcpy(atom2name, atomNames[atom_nums[0]].atomtype);
      strcpy(atom1name, atomNames[atom_nums[1]].atomtype);
    }

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(bonds[num_read]);
    b->atom1=atom_nums[0];
    b->atom2=atom_nums[1];

    /*  Query the parameter object for the constants for    */
    /*  this bond            */
    params->assign_bond_index(atom1name, atom2name, b);

    /*  Make sure this isn't a fake bond meant for shake in x-plor.  */
    Real k, x0;
    params->get_bond_params(&k,&x0,b->bond_type);
    if (simParams->lonepairs) {
      // need to retain Lonepair bonds for Drude
      if ( k == 0. && !is_lp(b->atom1) && !is_lp(b->atom2)) --numBonds;
      else ++num_read;
    }
    else {
      if ( k == 0. ) --numBonds;  // fake bond
      else ++num_read;  // real bond
    }
  }

  /*  Tell user about our subterfuge  */
  if ( numBonds != origNumBonds ) {
    iout << iWARN << "Ignored " << origNumBonds - numBonds <<
            " bonds with zero force constants.\n" << endi;
    iout << iWARN <<
	"Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
  }

  return;
}
/*      END OF FUNCTION read_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION read_angles        */
/*                  */
/*   INPUTS:                */
/*  fd - File descriptor for .psf file        */
/*  params - Parameters object to query for parameters    */
/*                  */
/*  read_angles reads the angle parameters from the .psf file.      */
/*   This section of the .psf file consists of a list of triplets of    */
/*   atom indexes.  Each triplet represents three atoms connected via   */
/*   an angle bond.  The parameter object is queried to obtain the      */
/*   constants for each bond.            */
/*                  */
/************************************************************************/

void Molecule::read_angles(FILE *fd, Parameters *params)

{
  int atom_nums[3];  //  Atom numbers for the three atoms
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  register int j;      //  Loop counter
  int num_read=0;    //  Number of angles read so far
  int origNumAngles = numAngles;  // Number of angles in file
  /*  Alloc the array of angles          */
  angles=new Angle[numAngles];

  if (angles == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_angles");
  }

  /*  Loop through and read all the angles      */
  while (num_read < numAngles)
  {
    /*  Loop through the 3 atom indexes in the current angle*/
    for (j=0; j<3; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "ANGLES")-1;

      /*  Check to make sure the atom index doesn't   */
      /*  exceed the Number of Atoms      */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "ANGLES INDEX %d GREATER THAN NATOM %d IN ANGLES # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }
    }

    /*  Place the bond name that is alphabetically first  */
    /*  in the atom1name.  This is OK since the order of    */
    /*  atom1 and atom3 are interchangable.  And to search  */
    /*  the tree of angle parameters, we need the order     */
    /*  to be predictable.          */
    if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
         atomNames[atom_nums[2]].atomtype) < 0)
    {
      strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
      strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
      strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    }
    else
    {
      strcpy(atom1name, atomNames[atom_nums[2]].atomtype);
      strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
      strcpy(atom3name, atomNames[atom_nums[0]].atomtype);
    }

    /*  Assign the three atom indices      */
    angles[num_read].atom1=atom_nums[0];
    angles[num_read].atom2=atom_nums[1];
    angles[num_read].atom3=atom_nums[2];

    /*  Get the constant values for this bond from the  */
    /*  parameter object          */
    params->assign_angle_index(atom1name, atom2name, 
       atom3name, &(angles[num_read]), simParams->alchOn ? -1 : 0);
    if ( angles[num_read].angle_type == -1 ) {
      iout << iWARN << "ALCHEMY MODULE WILL REMOVE ANGLE OR RAISE ERROR\n"
           << endi;
    }

    /*  Make sure this isn't a fake angle meant for shake in x-plor.  */
    Real k, t0, k_ub, r_ub;
    if ( angles[num_read].angle_type == -1 ) { k = -1.;  k_ub = -1.; } else
    params->get_angle_params(&k,&t0,&k_ub,&r_ub,angles[num_read].angle_type);
    if ( k == 0. && k_ub == 0. ) --numAngles;  // fake angle
    else ++num_read;  // real angle
  }

  /*  Tell user about our subterfuge  */
  if ( numAngles != origNumAngles ) {
    iout << iWARN << "Ignored " << origNumAngles - numAngles <<
            " angles with zero force constants.\n" << endi;
  }

  return;
}
/*      END OF FUNCTION read_angles      */

/************************************************************************/
/*                  */
/*        FUNCTION read_dihedrals      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for the .psf file        */
/*  params - pointer to parameter object        */
/*                  */
/*  read_dihedreals reads the dihedral section of the .psf file.    */
/*   This section of the file contains a list of quartets of atom       */
/*   numbers.  Each quartet represents a group of atoms that form a     */
/*   dihedral bond.              */
/*                  */
/************************************************************************/

void Molecule::read_dihedrals(FILE *fd, Parameters *params)
{
  int atom_nums[4];  // The 4 atom indexes
  int last_atom_nums[4];  // Atom numbers from previous bond
  char atom1name[11];  // Atom type for atom 1
  char atom2name[11];  // Atom type for atom 2
  char atom3name[11];  // Atom type for atom 3
  char atom4name[11];  // Atom type for atom 4
  register int j;      // loop counter
  int num_read=0;    // number of dihedrals read so far
  int multiplicity=1;  // multiplicity of the current bond
  Bool duplicate_bond;  // Is this a duplicate of the last bond
  int num_unique=0;   // Number of unique dihedral bonds

  //  Initialize the array used to check for duplicate dihedrals
  for (j=0; j<4; j++)
    last_atom_nums[j] = -1;

  /*  Allocate an array to hold the Dihedrals      */
  dihedrals = new Dihedral[numDihedrals];

  if (dihedrals == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_dihedrals");
  }

  /*  Loop through and read all the dihedrals      */
  while (num_read < numDihedrals)
  {
    duplicate_bond = TRUE;

    /*  Loop through and read the 4 indexes for this bond   */
    for (j=0; j<4; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "DIHEDRALS")-1;

      /*  Check for an atom index that is too large  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "DIHEDRALS INDEX %d GREATER THAN NATOM %d IN DIHEDRALS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      //  Check to see if this atom matches the last bond
      if (atom_nums[j] != last_atom_nums[j])
      {
        duplicate_bond = FALSE;
      }

      last_atom_nums[j] = atom_nums[j];
    }

    /*  Get the atom types for the 4 atoms so we can look  */
    /*  up the constants in the parameter object    */
    strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
    strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

    //  Check to see if this is really a new bond or just
    //  a repeat of the last one
    if (duplicate_bond)
    {
      //  This is a duplicate, so increase the multiplicity
      multiplicity++;

      if (multiplicity == 2)
      {
        numMultipleDihedrals++;
      }
    }
    else
    {
      multiplicity=1;
      num_unique++;
    }

    /*  Assign the atom indexes        */
    dihedrals[num_unique-1].atom1=atom_nums[0];
    dihedrals[num_unique-1].atom2=atom_nums[1];
    dihedrals[num_unique-1].atom3=atom_nums[2];
    dihedrals[num_unique-1].atom4=atom_nums[3];

    /*  Get the constants for this dihedral bond    */
    params->assign_dihedral_index(atom1name, atom2name, 
       atom3name, atom4name, &(dihedrals[num_unique-1]),
       multiplicity, simParams->alchOn ? -1 : 0);
    if ( dihedrals[num_unique-1].dihedral_type == -1 ) {
      iout << iWARN << "ALCHEMY MODULE WILL REMOVE DIHEDRAL OR RAISE ERROR\n"
           << endi;
    }

    num_read++;
  }

  numDihedrals = num_unique;

  return;
}
/*      END OF FUNCTION read_dihedral      */

/************************************************************************/
/*                  */
/*        FUNCTION read_impropers      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*  params - parameter object          */
/*                  */
/*  read_impropers reads the improper section of the .psf file.  */
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.          */
/*                  */
/************************************************************************/

void Molecule::read_impropers(FILE *fd, Parameters *params)
{
  int atom_nums[4];  //  Atom indexes for the 4 atoms
  int last_atom_nums[4];  //  Atom indexes from previous bond
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  char atom4name[11];  //  Atom type for atom 4
  register int j;      //  Loop counter
  int num_read=0;    //  Number of impropers read so far
  int multiplicity=1;  // multiplicity of the current bond
  Bool duplicate_bond;  // Is this a duplicate of the last bond
  int num_unique=0;   // Number of unique dihedral bonds

  //  Initialize the array used to look for duplicate improper
  //  entries.  Set them all to -1 so we know nothing will match
  for (j=0; j<4; j++)
    last_atom_nums[j] = -1;

  /*  Allocate the array to hold the impropers      */
  impropers=new Improper[numImpropers];

  if (impropers == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_impropers");
  }

  /*  Loop through and read all the impropers      */
  while (num_read < numImpropers)
  {
    duplicate_bond = TRUE;

    /*  Loop through the 4 indexes for this improper  */
    for (j=0; j<4; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "IMPROPERS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "IMPROPERS INDEX %d GREATER THAN NATOM %d IN IMPROPERS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      if (atom_nums[j] != last_atom_nums[j])
      {
        duplicate_bond = FALSE;
      }

      last_atom_nums[j] = atom_nums[j];
    }

    /*  Get the atom types so we can look up the parameters */
    strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
    strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

    //  Check to see if this is a duplicate improper
    if (duplicate_bond)
    {
      //  This is a duplicate improper.  So we don't
      //  really count this entry, we just update
      //  the parameters object
      multiplicity++;

      if (multiplicity == 2)
      {
        //  Count the number of multiples.
        numMultipleImpropers++;
      }
    }
    else
    {
      //  Not a duplicate
      multiplicity = 1;
      num_unique++;
    }

    /*  Assign the atom indexes        */
    impropers[num_unique-1].atom1=atom_nums[0];
    impropers[num_unique-1].atom2=atom_nums[1];
    impropers[num_unique-1].atom3=atom_nums[2];
    impropers[num_unique-1].atom4=atom_nums[3];

    /*  Look up the constants for this bond      */
    params->assign_improper_index(atom1name, atom2name, 
       atom3name, atom4name, &(impropers[num_unique-1]),
       multiplicity);

    num_read++;
  }

  //  Now reset the numImpropers value to the number of UNIQUE
  //  impropers.  Sure, we waste a few entries in the improper_array
  //  on the master node, but it is very little space . . .
  numImpropers = num_unique;

  return;
}
/*      END OF FUNCTION read_impropers      */

/************************************************************************/
/*                  */
/*        FUNCTION read_crossterms      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*  params - parameter object          */
/*                  */
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.          */
/*                  */
/************************************************************************/

void Molecule::read_crossterms(FILE *fd, Parameters *params)

{
  int atom_nums[8];  //  Atom indexes for the 4 atoms
  int last_atom_nums[8];  //  Atom indexes from previous bond
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  char atom4name[11];  //  Atom type for atom 4
  char atom5name[11];  //  Atom type for atom 5
  char atom6name[11];  //  Atom type for atom 6
  char atom7name[11];  //  Atom type for atom 7
  char atom8name[11];  //  Atom type for atom 8
  register int j;      //  Loop counter
  int num_read=0;    //  Number of items read so far
  Bool duplicate_bond;  // Is this a duplicate of the last bond

  //  Initialize the array used to look for duplicate crossterm
  //  entries.  Set them all to -1 so we know nothing will match
  for (j=0; j<8; j++)
    last_atom_nums[j] = -1;

  /*  Allocate the array to hold the cross-terms */
  crossterms=new Crossterm[numCrossterms];

  if (crossterms == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_crossterms");
  }

  /*  Loop through and read all the cross-terms      */
  while (num_read < numCrossterms)
  {
    duplicate_bond = TRUE;

    /*  Loop through the 8 indexes for this cross-term */
    for (j=0; j<8; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "CROSS-TERMS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "CROSS-TERM INDEX %d GREATER THAN NATOM %d IN CROSS-TERMS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      if (atom_nums[j] != last_atom_nums[j])
      {
        duplicate_bond = FALSE;
      }

      last_atom_nums[j] = atom_nums[j];
    }

    /*  Get the atom types so we can look up the parameters */
    strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
    strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    strcpy(atom4name, atomNames[atom_nums[3]].atomtype);
    strcpy(atom5name, atomNames[atom_nums[4]].atomtype);
    strcpy(atom6name, atomNames[atom_nums[5]].atomtype);
    strcpy(atom7name, atomNames[atom_nums[6]].atomtype);
    strcpy(atom8name, atomNames[atom_nums[7]].atomtype);

    //  Check to see if this is a duplicate term
    if (duplicate_bond)
    {
      iout << iWARN << "Duplicate cross-term detected.\n" << endi;
    }

    /*  Assign the atom indexes        */
    crossterms[num_read].atom1=atom_nums[0];
    crossterms[num_read].atom2=atom_nums[1];
    crossterms[num_read].atom3=atom_nums[2];
    crossterms[num_read].atom4=atom_nums[3];
    crossterms[num_read].atom5=atom_nums[4];
    crossterms[num_read].atom6=atom_nums[5];
    crossterms[num_read].atom7=atom_nums[6];
    crossterms[num_read].atom8=atom_nums[7];

    /*  Look up the constants for this bond      */
    params->assign_crossterm_index(atom1name, atom2name, 
       atom3name, atom4name, atom5name, atom6name,
       atom7name, atom8name, &(crossterms[num_read]));

    if(!duplicate_bond) num_read++;
  }

  numCrossterms = num_read;

  return;
}
/*      END OF FUNCTION read_impropers      */

/************************************************************************/
/*                  */
/*      FUNCTION read_donors        */
/*                  */
/*  read_donors reads in the bond section of the .psf file.  This   */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*                  */
/*  Donor atoms are the heavy atoms to which hydrogens are bonded.      */
/*  There will always be a donor atom for each donor pair.  However,    */
/*  for a united-atom model there may not be an explicit hydrogen       */
/*  present, in which case the second atom index in the pair will be    */
/*  given as 0 in the PSF (and stored as -1 in this program's internal  */
/*  storage).                                                           */
/************************************************************************/

void Molecule::read_donors(FILE *fd)
{
  int d[2];               // temporary storage of donor atom index
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
  int num_no_hydr=0;      // Number of bonds with no hydrogen given

  /*  Allocate the array to hold the bonds      */
  donors=new Bond[numDonors];

  if (donors == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_donors");
  }

  /*  Loop through and read in all the donors      */
  while (num_read < numDonors)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      d[j]=NAMD_read_int(fd, "DONORS")-1;

      /*  Check to make sure the index isn't too big  */
      if (d[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg,
    "DONOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE",
          d[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      /*  Check if there is a hydrogen given */
      if (d[j] < 0)
                          num_no_hydr++;
    }

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(donors[num_read]);
    b->atom1=d[0];
    b->atom2=d[1];

    num_read++;
  }

  return;
}
/*      END OF FUNCTION read_donors      */


/************************************************************************/
/*                  */
/*      FUNCTION read_acceptors        */
/*                  */
/*  read_acceptors reads in the bond section of the .psf file.      */
/*  This section contains a list of pairs of numbers where each pair is */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*                  */
/*  Acceptor atoms are the heavy atoms to which hydrogens directly      */
/*  orient in a hydrogen bond interaction.  There will always be an     */
/*  acceptor atom for each acceptor pair.  The antecedent atom, to      */
/*  which the acceptor is bound, may not be given in the structure,     */
/*  however, in which case the second atom index in the pair will be    */
/*  given as 0 in the PSF (and stored as -1 in this program's internal  */
/*  storage).                                                           */
/************************************************************************/

void Molecule::read_acceptors(FILE *fd)
{
  int d[2];               // temporary storage of atom index
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
        int num_no_ante=0;      // number of pairs with no antecedent

  /*  Allocate the array to hold the bonds      */
  acceptors=new Bond[numAcceptors];

  if (acceptors == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_acceptors");
  }

  /*  Loop through and read in all the acceptors      */
  while (num_read < numAcceptors)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      d[j]=NAMD_read_int(fd, "ACCEPTORS")-1;

      /*  Check to make sure the index isn't too big  */
      if (d[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "ACCEPTOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE", d[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      /*  Check if there is an antecedent given */
      if (d[j] < 0)
                          num_no_ante++;
    }

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(acceptors[num_read]);
    b->atom1=d[0];
    b->atom2=d[1];

    num_read++;
  }

  return;
}
/*      END OF FUNCTION read_acceptors      */


/************************************************************************/
/*                  */
/*      FUNCTION read_exclusions      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*                  */
/*  read_exclusions reads in the explicit non-bonded exclusions     */
/*  from the .psf file.  This section is a little funky, so hang on.    */
/*  Ok, first there is a list of atom indexes that is NumExclusions     */
/*  long.  These are in some sense the atoms that will be exlcuded.     */
/*  Following this list is a list of NumAtoms length that is a list     */
/*  of indexes into the list of excluded atoms.  So an example.  Suppose*/
/*  we have a 5 atom simulation with 3 explicit exclusions.  The .psf   */
/*  file could look like:            */
/*                  */
/*  3!NNB                */
/*  3 4 5                */
/*  0 1 3 3 3              */
/*                  */
/*  This would mean that atom 1 has no explicit exclusions.  Atom 2     */
/*  has an explicit exclusion with atom 3.  Atom 3 has an explicit      */
/*  exclusion with atoms 4 AND 5.  And atoms 4 and 5 have no explicit   */
/*  exclusions.  Got it!?!  I'm not sure who dreamed this up . . .      */
/*                  */
/************************************************************************/

void Molecule::read_exclusions(FILE *fd)

{
  int *exclusion_atoms;  //  Array of indexes of excluded atoms
  register int num_read=0;    //  Number fo exclusions read in
  int current_index;  //  Current index value
  int last_index;    //  the previous index value
  register int insert_index=0;  //  index of where we are in exlcusions array

  /*  Allocate the array of exclusion structures and the array of */
  /*  exlcuded atom indexes          */
  exclusions      = new Exclusion[numExclusions];
  exclusion_atoms = new int[numExclusions];

  if ( (exclusions == NULL) || (exclusion_atoms == NULL) )
  {
    NAMD_die("memory allocation failed in Molecule::read_exclusions");
  }

  /*  First, read in the excluded atoms list      */
  for (num_read=0; num_read<numExclusions; num_read++)
  {
    /*  Read the atom number from the file. Subtract 1 to   */
    /*  convert the index from the 1 to NumAtoms used in the*/
    /*  file to the  0 to NumAtoms-1 that we need    */
    exclusion_atoms[num_read]=NAMD_read_int(fd, "IMPROPERS")-1;

    /*  Check for an illegal index        */
    if (exclusion_atoms[num_read] >= numAtoms)
    {
      char err_msg[128];

      sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN PSF FILE", exclusion_atoms[num_read]+1, numAtoms, num_read+1);
      NAMD_die(err_msg);
    }
  }

  /*  Now, go through and read the list of NumAtoms pointers into */
  /*  the array that we just read in        */
  last_index=0;

  for (num_read=0; num_read<numAtoms; num_read++)
  {
    /*  Read in the current index value      */
    current_index=NAMD_read_int(fd, "EXCLUSIONS");

    /*  Check for an illegal pointer      */
    if (current_index>numExclusions)
    {
      char err_msg[128];

      sprintf(err_msg, "EXCLUSION INDEX %d LARGER THAN NUMBER OF EXLCUSIONS %d IN PSF FILE, EXCLUSION #%d\n", 
         current_index+1, numExclusions, num_read);
      NAMD_die(err_msg);
    }

    /*  Check to see if it matches the last index.  If so   */
    /*  than this atom has no exclusions.  If not, then     */
    /*  we have to build some exclusions      */
    if (current_index != last_index)
    {
      /*  This atom has some exclusions.  Loop from   */
      /*  the last_index to the current index.  This  */
      /*  will include how ever many exclusions this  */
      /*  atom has          */
      for (insert_index=last_index; 
           insert_index<current_index; insert_index++)
      {
        /*  Assign the two atoms involved.      */
        /*  The first one is our position in    */
        /*  the list, the second is based on    */
        /*  the pointer into the index list     */
        int a1 = num_read;
        int a2 = exclusion_atoms[insert_index];
        if ( a1 < a2 ) {
          exclusions[insert_index].atom1 = a1;
          exclusions[insert_index].atom2 = a2;
        } else if ( a2 < a1 ) {
          exclusions[insert_index].atom1 = a2;
          exclusions[insert_index].atom2 = a1;
        } else {
          char err_msg[128];
          sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN PSF FILE\n", a1+1);
          NAMD_die(err_msg);
        }
      }

      last_index=current_index;
    }
  }

  /*  Free our temporary list of indexes        */
  delete [] exclusion_atoms;

  return;
}
/*      END OF FUNCTION read_exclusions      */

/************************************************************************/
/*      FUNCTION read_exclusions      */
/*                  */
/*   INPUTS:                */
/*   int* atom_i - array of atom i indices    */
/*   int* atom_j - array of atom j indices    */
/*   int num_exclusion - length of array           */
/*                  */
/* JLai August 16th, 2012 */
/************************************************************************/
void Molecule::read_exclusions(int* atom_i, int* atom_j, int num_exclusion)
{
    /*  Allocate the array of exclusion structures and the array of */
  /*  exlcuded atom indexes          */
  exclusions       = new Exclusion[num_exclusion];
  int loop_counter = 0;  
  int a=0;
  int b=0;

  if ( (exclusions == NULL) )
  {
    NAMD_die("memory allocation failed in Molecule::read_exclusions");
  }

  /* The following code only guarantees that exclusion.atom1 is < exclusion.atom2 */
  for (loop_counter = 0; loop_counter < num_exclusion; loop_counter++) {
	
	if ( (atom_i == NULL) || (atom_j == NULL) ) {
	  NAMD_die("null pointer expection in Molecule::read_exclusions");
	}

	a = atom_i[loop_counter];
	b = atom_j[loop_counter];
	if(a < b) {
		exclusions[loop_counter].atom1 = a;
		exclusions[loop_counter].atom2 = b;
	} else {
		exclusions[loop_counter].atom1 = b;
		exclusions[loop_counter].atom2 = a;
	}
	exclusionSet.add(Exclusion(exclusions[loop_counter].atom1,exclusions[loop_counter].atom2));
  }

  if ( ! CkMyPe() ) {
    iout << iINFO << "ADDED " << num_exclusion << " EXPLICIT EXCLUSIONS: THIS VALUE WILL *NOT* BE ADDED TO THE STRUCTURE SUMMARY\n" << endi;
  }

   return;
}
/*      END OF FUNCTION read_exclusions      */

/************************************************************************/
/*                  */
/*        FUNCTION read_lphosts    */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*                  */
/*  this function reads in the lone pair host section of the .psf file. */
/*                  */
void Molecule::read_lphosts(FILE *fd)
{
  char buffer[512];  // Buffer for reading from file
  char lptype[8];
  int numhosts, index, i, read_count;
  Real distance, angle, dihedral;

  lphosts = new Lphost[numLphosts];
  if (lphosts == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_lphosts");
  }
  for (i = 0;  i < numLphosts;  i++)
  {
    NAMD_read_line(fd, buffer);
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') ) continue;
    read_count=sscanf(buffer, "%d %d %6s %f %f %f",
        &numhosts, &index, lptype, &distance, &angle, &dihedral);
    if (read_count != 6 || numhosts != 3 || index != 4*i + 1
        || strcmp(lptype,"F") != 0)
    {
      char err_msg[128];
      sprintf(err_msg, "BAD FORMAT FOR LPHOST LINE %d IN PSF FILE LINE\n"
          "LINE=%s\n", i+1, buffer);
      NAMD_die(err_msg);
    }
    lphosts[i].distance = distance;
    lphosts[i].angle = angle * (M_PI/180);        // convert to radians
    lphosts[i].dihedral = dihedral * (M_PI/180);  // convert to radians
  }
  for (i = 0;  i < numLphosts;  i++) {
    lphosts[i].atom1 = NAMD_read_int(fd, "LPHOSTS")-1;
    lphosts[i].atom2 = NAMD_read_int(fd, "LPHOSTS")-1;
    lphosts[i].atom3 = NAMD_read_int(fd, "LPHOSTS")-1;
    lphosts[i].atom4 = NAMD_read_int(fd, "LPHOSTS")-1;
  }
}
/*      END OF FUNCTION read_lphosts    */

/************************************************************************/
/*                  */
/*        FUNCTION read_anisos     */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*                  */
/*  this function reads in the anisotropic terms section of .psf file. */
/*                  */
void Molecule::read_anisos(FILE *fd)
{
  char buffer[512];  // Buffer for reading from file
  int numhosts, index, i, read_count;
  Real k11, k22, k33;

  anisos = new Aniso[numAnisos];
  if (anisos == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_anisos");
  }
  for (i = 0;  i < numAnisos;  i++)
  {
    NAMD_read_line(fd, buffer);
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') ) continue;
    read_count=sscanf(buffer, "%f %f %f", &k11, &k22, &k33);
    if (read_count != 3)
    {
      char err_msg[128];
      sprintf(err_msg, "BAD FORMAT FOR ANISO LINE %d IN PSF FILE LINE\n"
          "LINE=%s\n", i+1, buffer);
      NAMD_die(err_msg);
    }
    anisos[i].k11 = k11;
    anisos[i].k22 = k22;
    anisos[i].k33 = k33;
  }
  for (i = 0;  i < numAnisos;  i++) {
    anisos[i].atom1 = NAMD_read_int(fd, "ANISOS")-1;
    anisos[i].atom2 = NAMD_read_int(fd, "ANISOS")-1;
    anisos[i].atom3 = NAMD_read_int(fd, "ANISOS")-1;
    anisos[i].atom4 = NAMD_read_int(fd, "ANISOS")-1;
  }
}
/*      END OF FUNCTION read_anisos     */

//LCPO
inline int getLCPOTypeAmber(char atomType[11], int numBonds) {

  //Hydrogen
  if (atomType[0] == 'H' || atomType[0] == 'h') {
    return 0;

  //Carbon
  } else if (atomType[0] == 'C' || atomType[0] == 'c') {
    if (//Sp3 Carbon
      //atomType[1] == 'T')// ||
      strcmp(atomType, "CT" )==0 )
      //strcmp(atomType, "CP1" )==0 ||
      //strcmp(atomType, "CP2" )==0 ||
      //strcmp(atomType, "CP3" )==0 ||
      //strcmp(atomType, "CS"  )==0 )
      {
      if (numBonds == 1)
        return 1;
      else if (numBonds == 2)
        return 2;
      else if (numBonds == 3)
        return 3;
      else if (numBonds == 4)
        return 4;
      else
        return 1;

    } else {//Sp2 or other
      if (numBonds == 2)
        return 5;
      else if (numBonds == 3)
        return 6;
      else
        return 1;
    }

  //Nitrogen
  } else if (atomType[0] == 'N' || atomType[0] == 'n') {
    if ( strcmp(atomType, "N3"  ) == 0 ) { //Sp3 Nitrogen
      if (numBonds == 1)
        return 11;
      else if (numBonds == 2)
        return 12;
      else if (numBonds == 3)
        return 13;
      else
        return 11;

    } else {//SP2 Nitrogen
      if (numBonds == 1)
        return 14;
      else if (numBonds == 2)
        return 15;
      else if (numBonds == 3)
        return 16;
      else
        return 11;
    }

  //Oxygen
  } else if (atomType[0] == 'O' || atomType[0] == 'o') {

    if ( strcmp(atomType, "O" )==0) {//Sp2 Oxygen
      return 9;
    } else if (strcmp(atomType, "O2" )==0) {//Carboxylate Oxygen
      return 10;
    } else { // Sp3 Oxygen
      if (numBonds == 1)
        return 7;
      else if (numBonds == 2)
        return 8;
      else
        return 7;
    }

  //Sulfur
  } else if (atomType[0] == 'S' || atomType[0] == 's') {
    if ( strcmp(atomType, "SH" )==0) { //Sulfur 1 neighbor
      return 17;
    } else {
      return 18;
    }

  //Phosphorus
  } else if (atomType[0] == 'P' || atomType[0] == 'p') {
      if (numBonds == 3)
        return 19;
      else if (numBonds == 4)
        return 20;
      else
        return 19;
  } else if (atomType[0] == 'Z') { // ? just to agree with Amber mdread.f
    return 0;
  } else  if ( strcmp(atomType, "MG" )==0) { //Mg
    return 22;
  } else { // unknown atom type
    return 5;
  }
  return 5;
} // getLCPOTypeAmber

inline int getLCPOTypeCharmm(char atomType[11], int numBonds) {

  //Hydrogen
  if (atomType[0] == 'H') {
    return 0;

  //Carbon
  } else if (atomType[0] == 'C') {
    if (//Sp3 Carbon
      atomType[1] == 'T' ||
      strcmp(atomType, "CP1" )==0 ||
      strcmp(atomType, "CP2" )==0 ||
      strcmp(atomType, "CP3" )==0 ||
      strcmp(atomType, "CS"  )==0 ) {
      if (numBonds == 1)
        return 1;
      else if (numBonds == 2)
        return 2;
      else if (numBonds == 3)
        return 3;
      else if (numBonds == 4)
        return 4;
      else
        return 1;

    } else if (//Sp2
      strcmp(atomType, "C"   )==0 ||
      strcmp(atomType, "CA"  )==0 ||
      strcmp(atomType, "CC"  )==0 ||
      strcmp(atomType, "CD"  )==0 ||
      strcmp(atomType, "CN"  )==0 ||
      strcmp(atomType, "CY"  )==0 ||
      strcmp(atomType, "C3"  )==0 ||
      strcmp(atomType, "CE1" )==0 ||
      strcmp(atomType, "CE2" )==0 ||
      strcmp(atomType, "CST" )==0 ||
      strcmp(atomType, "CAP" )==0 ||
      strcmp(atomType, "COA" )==0 ||
      strcmp(atomType, "CPT" )==0 ||
      strcmp(atomType, "CPH1")==0 ||
      strcmp(atomType, "CPH2")==0
      ) {
      if (numBonds == 2)
        return 5;
      else if (numBonds == 3)
        return 6;
      else
        return 1;
    } else { // other Carbon
        return 1;
    }

  //Nitrogen
  } else if (atomType[0] == 'N') {
    if (//Sp3 Nitrogen
      //strcmp(atomType, "N"   )==0 ||
      //strcmp(atomType, "NH1" )==0 ||
      //strcmp(atomType, "NH2" )==0 ||
      strcmp(atomType, "NH3" )==0 ||
      //strcmp(atomType, "NC2" )==0 ||
      //strcmp(atomType, "NY"  )==0 ||
      strcmp(atomType, "NP"  )==0
      ) {
      if (numBonds == 1)
        return 11;
      else if (numBonds == 2)
        return 12;
      else if (numBonds == 3)
        return 13;
      else
        return 11;

    } else if (//SP2 Nitrogen
      strcmp(atomType, "NY"  )==0 || //
      strcmp(atomType, "NC2" )==0 || //
      strcmp(atomType, "N"   )==0 || //
      strcmp(atomType, "NH1" )==0 || //
      strcmp(atomType, "NH2" )==0 || //
      strcmp(atomType, "NR1" )==0 ||
      strcmp(atomType, "NR2" )==0 ||
      strcmp(atomType, "NR3" )==0 ||
      strcmp(atomType, "NPH" )==0 ||
      strcmp(atomType, "NC"  )==0
      ) {
      if (numBonds == 1)
        return 14;
      else if (numBonds == 2)
        return 15;
      else if (numBonds == 3)
        return 16;
      else
        return 11;
    } else { // other Nitrogen
      return 11;
    }

  //Oxygen
  } else if (atomType[0] == 'O') {
    if (//Sp3 Oxygen
      strcmp(atomType, "OH1" )==0 ||
      strcmp(atomType, "OS"  )==0 ||
      strcmp(atomType, "OC"  )==0 || //
      strcmp(atomType, "OT"  )==0
      ) {
      if (numBonds == 1)
        return 7;
      else if (numBonds == 2)
        return 8;
      else
        return 7;
    } else if ( // Sp2 Oxygen
      strcmp(atomType, "O"   )==0 ||
      strcmp(atomType, "OB"  )==0 ||
      strcmp(atomType, "OST" )==0 ||
      strcmp(atomType, "OCA" )==0 ||
      strcmp(atomType, "OM"  )==0
      ) {
      return 9;
    } else if ( // SP1 Oxygen
      strcmp(atomType, "OC"  )==0
      ) {
      return 10;
    } else { // other Oxygen
      return 7;
    }

  //Sulfur
  } else if (atomType[0] == 'S') {
      if (numBonds == 1)
        return 17;
      else
        return 18;

  //Phosphorus
  } else if (atomType[0] == 'P') {
      if (numBonds == 3)
        return 19;
      else if (numBonds == 4)
        return 20;
      else
        return 19;
  } else { // unknown atom type
    return 5;
  }
  return 5;
} // getLCPOTypeCharmm

//input type is Charmm/Amber/other
//0 - Charmm/Xplor
//1 - Amber
//2 - Plugin
//3 - Gromacs
void Molecule::assignLCPOTypes(int inputType) {
  int *heavyBonds = new int[numAtoms];
  for (int i = 0; i < numAtoms; i++)
    heavyBonds[i] = 0;
  for (int i = 0; i < numBonds; i++ ) {
    Bond curBond = bonds[i];
    int a1 = bonds[i].atom1;
    int a2 = bonds[i].atom2;
    if (atoms[a1].mass > 2.f && atoms[a2].mass > 2.f) {
      heavyBonds[a1]++;
      heavyBonds[a2]++;
    }
  }

  lcpoParamType = new int[numAtoms];

  int warning = 0;
  for (int i = 0; i < numAtoms; i++) {
    //print vdw_type and numbonds

    if (inputType == 1) { // Amber
      lcpoParamType[i] = getLCPOTypeAmber(atomNames[i].atomtype, heavyBonds[i]);
    } else { // Charmm
      lcpoParamType[i] = getLCPOTypeCharmm(atomNames[i].atomtype, heavyBonds[i]);
    }
/*
    CkPrintf("%d MOL: ATOM[%05d] = { %4s %d } : %d\n",
      inputType,
      i+1,
      atomNames[i].atomtype,
      heavyBonds[i],
      lcpoParamType[i]
      );
*/
    if ( atoms[i].mass < 1.5 && lcpoParamType[i] != 0 ) {
      if (atoms[i].status & LonepairAtom) {
        warning |= LonepairAtom;
        lcpoParamType[i] = 0;  // reset LCPO type for LP
      }
      else if (atoms[i].status & DrudeAtom) {
        warning |= DrudeAtom;
        lcpoParamType[i] = 0;  // reset LCPO type for Drude
      }
      else {
        CkPrintf("ERROR in Molecule::assignLCPOTypes(): "
            "Light atom given heavy atom LCPO type.\n");
      }
    }

    //CkPrintf("VDW_TYPE %02d %4s\n", atoms[i].vdw_type, atomNames[i].atomtype);
  } // for atoms

  if (warning & LonepairAtom) {
    iout << iWARN << "LONE PAIRS TO BE IGNORED BY SASA\n" << endi;
  }
  if (warning & DrudeAtom) {
    iout << iWARN << "DRUDE PARTICLES TO BE IGNORED BY SASA\n" << endi;
  }

  delete [] heavyBonds;

} // buildLCPOTable

void Molecule::plgLoadAtomBasics(molfile_atom_t *atomarray){
    atoms = new Atom[numAtoms];
    atomNames = new AtomNameInfo[numAtoms];
    if(simParams->genCompressedPsf) {
        atomSegResids = new AtomSegResInfo[numAtoms];
    }    
    hydrogenGroup.resize(0);

    ResidueLookupElem *tmpResLookup = resLookup;

    for(int i=0; i<numAtoms; i++) {
        int reslength = strlen(atomarray[i].resname)+1;
        int namelength = strlen(atomarray[i].name)+1;
        int typelength = strlen(atomarray[i].type)+1;
        atomNames[i].resname = nameArena->getNewArray(reslength);
        atomNames[i].atomname = nameArena->getNewArray(namelength);
        atomNames[i].atomtype = nameArena->getNewArray(typelength);
        strcpy(atomNames[i].resname, atomarray[i].resname);
        strcpy(atomNames[i].atomname, atomarray[i].name);
        strcpy(atomNames[i].atomtype, atomarray[i].type);

        atoms[i].mass = atomarray[i].mass;
        atoms[i].charge = atomarray[i].charge;
        atoms[i].status = UnknownAtom;

        //add this atom to residue lookup table
        if(tmpResLookup) {
            tmpResLookup = tmpResLookup->append(atomarray[i].segid, atomarray[i].resid, i);
        }

        if(atomSegResids) { //for compressing molecule information
            AtomSegResInfo *one = atomSegResids + i;
            memcpy(one->segname, atomarray[i].segid, strlen(atomarray[i].segid)+1);
            one->resid = atomarray[i].resid;
        }
        //Determine the type of the atom
        if ( simParams->ignoreMass ) {
        }else if(atoms[i].mass <= 0.05) {
            atoms[i].status |= LonepairAtom;
        }else if(atoms[i].mass < 1.0) {
            atoms[i].status |= DrudeAtom;
        }else if(atoms[i].mass <= 3.5) {
            atoms[i].status |= HydrogenAtom;
        }else if((atomNames[i].atomname[0] == 'O') &&
                 (atoms[i].mass>=14.0) && (atoms[i].mass<=18.0)){
            atoms[i].status |= OxygenAtom;
        }
        //Look up the vdw constants for this atom
        params->assign_vdw_index(atomNames[i].atomtype, &atoms[i]);
    }
}

void Molecule::plgLoadBonds(int *from, int *to){
    bonds = new Bond[numBonds];
    char atom1name[11];
    char atom2name[11];
    int realNumBonds = 0;
    for(int i=0; i<numBonds; i++) {
        Bond *thisBond = bonds+realNumBonds;
        thisBond->atom1 = from[i]-1;
        thisBond->atom2 = to[i]-1;
        /* Get the atom type for the two atoms.
         * When we query the parameter object, we
         * need to send the atom type that is alphabetically
         * first as atom 1.
         */
        if(strcasecmp(atomNames[thisBond->atom1].atomtype,
                      atomNames[thisBond->atom2].atomtype)<0) {
            strcpy(atom1name, atomNames[thisBond->atom1].atomtype);
            strcpy(atom2name, atomNames[thisBond->atom2].atomtype);
        }else{
            strcpy(atom2name, atomNames[thisBond->atom1].atomtype);
            strcpy(atom1name, atomNames[thisBond->atom2].atomtype);
        }
        params->assign_bond_index(atom1name, atom2name, thisBond);

        //Make sure this isn't a fake bond meant for shake in x-plor
        Real k, x0;
        params->get_bond_params(&k, &x0, thisBond->bond_type);
        if(simParams->lonepairs) {
            //need to retain Lonepair bonds for Drude
            if(k!=0. || is_lp(thisBond->atom1) || 
               is_lp(thisBond->atom2)) {               
                realNumBonds++;
	    }
        }else{
            if(k != 0.) realNumBonds++;
        }
    }

    if(numBonds != realNumBonds) {
        iout << iWARN << "Ignored" << numBonds-realNumBonds <<
            "bonds with zero force constants.\n" <<endi;
        iout << iWARN << "Will get H-H distance in rigid H20 from H-O-H angle.\n" <<endi;
    }
    numBonds = realNumBonds;
}

void Molecule::plgLoadAngles(int *plgAngles)
{    
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];

    angles=new Angle[numAngles];
    int *atomid = plgAngles;
    int numRealAngles = 0;
    for(int i=0; i<numAngles; i++) {
        Angle *thisAngle = angles+numRealAngles;
        thisAngle->atom1 = atomid[0]-1;
        thisAngle->atom2 = atomid[1]-1;
        thisAngle->atom3 = atomid[2]-1;
        atomid += 3;

        if(strcasecmp(atomNames[thisAngle->atom1].atomtype,
                      atomNames[thisAngle->atom2].atomtype)<0) {
            strcpy(atom1name, atomNames[thisAngle->atom1].atomtype);
            strcpy(atom2name, atomNames[thisAngle->atom2].atomtype);
            strcpy(atom3name, atomNames[thisAngle->atom3].atomtype);
        }else{
            strcpy(atom1name, atomNames[thisAngle->atom3].atomtype);
            strcpy(atom2name, atomNames[thisAngle->atom2].atomtype);
            strcpy(atom3name, atomNames[thisAngle->atom1].atomtype);
        }

        params->assign_angle_index(atom1name, atom2name, atom3name,
				thisAngle, simParams->alchOn ? -1 : 0);
        if ( thisAngle->angle_type == -1 ) {
          iout << iWARN << "ALCHEMY MODULE WILL REMOVE ANGLE OR RAISE ERROR\n"
               << endi;
        }

        Real k, t0, k_ub, r_ub;
        if ( thisAngle->angle_type == -1 ) { k = -1.;  k_ub = -1.; } else
        params->get_angle_params(&k, &t0, &k_ub, &r_ub, thisAngle->angle_type);
        if(k!=0. || k_ub!=0.) numRealAngles++;
    }

    if(numAngles != numRealAngles) {
        iout << iWARN << "Ignored" << numAngles-numRealAngles << 
            " angles with zero force constants.\n" << endi; 
    }
    numAngles = numRealAngles;
}

void Molecule::plgLoadDihedrals(int *plgDihedrals)
{
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];
    char atom4name[11];
    int lastAtomIds[4];
    int multiplicity = 1; //multiplicity of the current bond

    lastAtomIds[0]=lastAtomIds[1]=lastAtomIds[2]=lastAtomIds[3]=-1;
    dihedrals = new Dihedral[numDihedrals];
    int numRealDihedrals = 0;
    int *atomid = plgDihedrals;
    for(int i=0; i<numDihedrals; i++, atomid+=4) {
        Dihedral *thisDihedral = dihedrals + numRealDihedrals;
        Bool duplicate_bond = TRUE;
        for(int j=0; j<4; j++) {
            if(atomid[j] != lastAtomIds[j]) {
                duplicate_bond = FALSE;
            }
            lastAtomIds[j] = atomid[j];
        }

        strcpy(atom1name, atomNames[atomid[0]-1].atomtype);
        strcpy(atom2name, atomNames[atomid[1]-1].atomtype);
        strcpy(atom3name, atomNames[atomid[2]-1].atomtype);
        strcpy(atom4name, atomNames[atomid[3]-1].atomtype);

        if(duplicate_bond) {
            multiplicity++;
            if(multiplicity==2) {
                numMultipleDihedrals++;
            }
        }else{
            multiplicity=1;
            numRealDihedrals++;
        }

        thisDihedral->atom1 = atomid[0]-1;
        thisDihedral->atom2 = atomid[1]-1;
        thisDihedral->atom3 = atomid[2]-1;
        thisDihedral->atom4 = atomid[3]-1;

        params->assign_dihedral_index(atom1name, atom2name,
                                      atom3name, atom4name, thisDihedral,
                                      multiplicity, simParams->alchOn ? -1 : 0);
        if ( thisDihedral->dihedral_type == -1 ) {
          iout << iWARN << "ALCHEMY MODULE WILL REMOVE DIHEDRAL OR RAISE ERROR\n"
               << endi;
        }
    }

    numDihedrals = numRealDihedrals;
}

void Molecule::plgLoadImpropers(int *plgImpropers)
{
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];
    char atom4name[11];
    int lastAtomIds[4];
    int multiplicity = 1; //multiplicity of the current bond

    lastAtomIds[0]=lastAtomIds[1]=lastAtomIds[2]=lastAtomIds[3]=-1;
    impropers = new Improper[numImpropers];
    int numRealImpropers = 0;
    int *atomid = plgImpropers;
    for(int i=0; i<numImpropers; i++, atomid+=4) {
        Improper *thisImproper = impropers + numRealImpropers;
        Bool duplicate_bond = TRUE;
        for(int j=0; j<4; j++) {
            if(atomid[j] != lastAtomIds[j]) {
                duplicate_bond = FALSE;
            }
            lastAtomIds[j] = atomid[j];
        }

        strcpy(atom1name, atomNames[atomid[0]-1].atomtype);
        strcpy(atom2name, atomNames[atomid[1]-1].atomtype);
        strcpy(atom3name, atomNames[atomid[2]-1].atomtype);
        strcpy(atom4name, atomNames[atomid[3]-1].atomtype);

        if(duplicate_bond) {
            multiplicity++;
            if(multiplicity==2) {
                numMultipleImpropers++;
            }
        }else{
            multiplicity=1;
            numRealImpropers++;
        }

        thisImproper->atom1 = atomid[0]-1;
        thisImproper->atom2 = atomid[1]-1;
        thisImproper->atom3 = atomid[2]-1;
        thisImproper->atom4 = atomid[3]-1;

        params->assign_improper_index(atom1name, atom2name,
                                      atom3name, atom4name, thisImproper,
                                      multiplicity);
    }

    numImpropers = numRealImpropers;
}

void Molecule::plgLoadCrossterms(int *plgCterms)
{
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];
    char atom4name[11];
    char atom5name[11];
    char atom6name[11];
    char atom7name[11];
    char atom8name[11];
    int lastAtomIds[8];    

    for(int i=0; i<8; i++)
        lastAtomIds[i]=-1;
    
    crossterms = new Crossterm[numCrossterms];
    int numRealCrossterms = 0;
    int *atomid = plgCterms;
    for(int i=0; i<numCrossterms; i++, atomid+=8) {
        Crossterm *thisCrossterm = crossterms + numRealCrossterms;
        Bool duplicate_bond = TRUE;
        for(int j=0; j<8; j++) {
            if(atomid[j] != lastAtomIds[j]) {
                duplicate_bond = FALSE;
            }
            lastAtomIds[j] = atomid[j];
        }

        strcpy(atom1name, atomNames[atomid[0]-1].atomtype);
        strcpy(atom2name, atomNames[atomid[1]-1].atomtype);
        strcpy(atom3name, atomNames[atomid[2]-1].atomtype);
        strcpy(atom4name, atomNames[atomid[3]-1].atomtype);
        strcpy(atom5name, atomNames[atomid[4]-1].atomtype);
        strcpy(atom6name, atomNames[atomid[5]-1].atomtype);
        strcpy(atom7name, atomNames[atomid[6]-1].atomtype);
        strcpy(atom8name, atomNames[atomid[7]-1].atomtype);

        if(duplicate_bond) {
            iout << iWARN <<"Duplicate cross-term detected.\n" << endi;
        } else
            numRealCrossterms++;

        thisCrossterm->atom1 = atomid[0]-1;
        thisCrossterm->atom2 = atomid[1]-1;
        thisCrossterm->atom3 = atomid[2]-1;
        thisCrossterm->atom4 = atomid[3]-1;
        thisCrossterm->atom5 = atomid[4]-1;
        thisCrossterm->atom6 = atomid[5]-1;
        thisCrossterm->atom7 = atomid[6]-1;
        thisCrossterm->atom8 = atomid[7]-1;

        params->assign_crossterm_index(atom1name, atom2name,
                                       atom3name, atom4name, atom5name,
                                       atom6name, atom7name, atom8name,
                                       thisCrossterm);
    }

    numCrossterms = numRealCrossterms;
}

void Molecule::setOccupancyData(molfile_atom_t *atomarray){
    occupancy = new float[numAtoms];
    for(int i=0; i<numAtoms; i++) {
        occupancy[i] = atomarray[i].occupancy;
    }
}

void Molecule::setBFactorData(molfile_atom_t *atomarray){
    bfactor = new float[numAtoms];
    for(int i=0; i<numAtoms; i++) {
        bfactor[i] = atomarray[i].bfactor;
    }
}

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_lists_by_atom      */
    /*                  */
    /*  This function builds O(NumAtoms) arrays that store the bonds,   */
    /*  angles, dihedrals, and impropers, that each atom is involved in.    */
    /*  This is a space hog, but VERY fast.  This will certainly have to    */
    /*  change to make things scalable in memory, but for now, speed is the */
    /*  thing!                */
    /*                  */
    /************************************************************************/
    void Molecule::build_lists_by_atom()       
    {
       register int i;      //  Loop counter

       register int numFixedAtoms = this->numFixedAtoms;
       // if we want forces on fixed atoms then just pretend
       // there are none for the purposes of this routine;
       if ( simParams->fixedAtomsForces ) numFixedAtoms = 0;

//fepb
//     int numFepInitial = this->numFepInitial;
//     int numFepFinal = this->numFepFinal;
//fepe
       tmpArena = new ObjectArena<int32>;
       bondsWithAtom = new int32 *[numAtoms];
       cluster = new int32 [numAtoms];
       clusterSize = new int32 [numAtoms];

       bondsByAtom = new int32 *[numAtoms];
       anglesByAtom = new int32 *[numAtoms];
       dihedralsByAtom = new int32 *[numAtoms];
       impropersByAtom = new int32 *[numAtoms];
       crosstermsByAtom = new int32 *[numAtoms];

       exclusionsByAtom = new int32 *[numAtoms];
       fullExclusionsByAtom = new int32 *[numAtoms];
       modExclusionsByAtom = new int32 *[numAtoms];

       // JLai
       gromacsPairByAtom = new int32 *[numAtoms];
       // End of JLai

       int32 *byAtomSize = new int32[numAtoms];

       const int pair_self = 
         simParams->pairInteractionOn ? simParams->pairInteractionSelf : 0;

       DebugM(3,"Building bond lists.\n");
    
       //  Build the bond lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       for (i=0; i<numRealBonds; i++)
       {
         byAtomSize[bonds[i].atom1]++;
         byAtomSize[bonds[i].atom2]++;
       }
       for (i=0; i<numAtoms; i++)
       {
         bondsWithAtom[i] = tmpArena->getNewArray(byAtomSize[i]+1);
         bondsWithAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numRealBonds; i++)
       {
         int a1 = bonds[i].atom1;
         int a2 = bonds[i].atom2;
         bondsWithAtom[a1][byAtomSize[a1]++] = i;
         bondsWithAtom[a2][byAtomSize[a2]++] = i;
       }

        
       // Updates all bond, angle, dihedral, improper and crossterm
       // to reflect the QM region (which can't have any of there terms)
       if (simParams->qmForcesOn) {
           
           DebugM(3,"Calculating exclusions for QM simulation.\n");
           build_exclusions();
           
           delete_qm_bonded() ;
           
           DebugM(3,"Re-Building bond lists.\n");
           
           // We re-calculate the bondsWithAtom list for cluster 
           // info calculation below.
           for (i=0; i<numAtoms; i++)
           {
             byAtomSize[i] = 0;
           }
           for (i=0; i<numRealBonds; i++)
           {
             byAtomSize[bonds[i].atom1]++;
             byAtomSize[bonds[i].atom2]++;
           }
           for (i=0; i<numAtoms; i++)
           {
             bondsWithAtom[i][byAtomSize[i]] = -1;
             byAtomSize[i] = 0;
           }
           for (i=0; i<numRealBonds; i++)
           {
             int a1 = bonds[i].atom1;
             int a2 = bonds[i].atom2;
             bondsWithAtom[a1][byAtomSize[a1]++] = i;
             bondsWithAtom[a2][byAtomSize[a2]++] = i;
           }
       }
        
       //  Build cluster information (contiguous clusters)
       for (i=0; i<numAtoms; i++) {
         cluster[i] = i;
       }
       for (i=0; i<numAtoms; i++) {
         int ci = i;
         while ( cluster[ci] != ci ) ci = cluster[ci];
         for ( int32 *b = bondsWithAtom[i]; *b != -1; ++b ) {
           int a = bonds[*b].atom1;
           if ( a == i ) a = bonds[*b].atom2;
           if ( a > i ) {
             int ca = a;
             while ( cluster[ca] != ca ) ca = cluster[ca];
             if ( ca > ci ) cluster[ca] = cluster[ci];
             else cluster[ci] = cluster[ca];
           }
         }
       }
       while ( 1 ) {
         int allok = 1;
         for (i=0; i<numAtoms; i++) {
           int ci = cluster[i];
           if ( cluster[ci] != ci ) {
             allok = 0;
             cluster[i] = cluster[ci];
           }
         }
         if ( allok ) break;
       }
       
       for (i=0; i<numAtoms; i++) {
         clusterSize[i] = 0;
       }       
       for (i=0; i<numAtoms; i++) {           
         clusterSize[cluster[i]] += 1;
       }

/*
       //Getting number of clusters for debugging
       int numClusters=0;
       for(int i=0; i<numAtoms; i++){
           if(clusterSize[i]!=0) numClusters++;
       }
       printf("Num of clusters: %d\n", numClusters);
*/

       //  Build the bond lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcBonds = 0;
       for (i=0; i<numBonds; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
   
         if ( pair_self && fepAtomFlags[bonds[i].atom1] != 1) continue;
         byAtomSize[bonds[i].atom1]++;
         numCalcBonds++;
       }
       for (i=0; i<numAtoms; i++)
       {
         bondsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         bondsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numBonds; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
         if ( pair_self && fepAtomFlags[bonds[i].atom1] != 1) continue;
         int a1 = bonds[i].atom1;
         bondsByAtom[a1][byAtomSize[a1]++] = i;
       }
       for (i=0; i<numBonds; ++i) {
         int a1 = bonds[i].atom1;
         int a2 = bonds[i].atom2;
         int j;
         if ( a1 == a2 ) {
           char buff[512];
           sprintf(buff,"Atom %d is bonded to itself", a1+1);
           NAMD_die(buff);
         }
         for ( j = 0; j < byAtomSize[a1]; ++j ) {
           int b = bondsByAtom[a1][j];
           int ba1 = bonds[b].atom1;
           int ba2 = bonds[b].atom2;
           if ( b != i && ( (ba1==a1 && ba2==a2) || (ba1==a2 && ba2==a1) ) ) {
             char buff[512];
             sprintf(buff,"Duplicate bond from atom %d to atom %d", a1+1, a2+1);
             NAMD_die(buff);
           }
         }
         for ( j = 0; j < byAtomSize[a2]; ++j ) {
           int b = bondsByAtom[a2][j];
           int ba1 = bonds[b].atom1;
           int ba2 = bonds[b].atom2;
           if ( b != i && ( (ba1==a1 && ba2==a2) || (ba1==a2 && ba2==a1) ) ) {
             char buff[512];
             sprintf(buff,"Duplicate bond from atom %d to atom %d", a1+1, a2+1);
             NAMD_die(buff);
           }
         }
       }
       
       DebugM(3,"Building angle lists.\n");
    
       //  Build the angle lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcAngles = 0;
       for (i=0; i<numAngles; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
         if ( pair_self && fepAtomFlags[angles[i].atom1] != 1) continue;
         byAtomSize[angles[i].atom1]++;
         numCalcAngles++;
       }
       for (i=0; i<numAtoms; i++)
       {
         anglesByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         anglesByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numAngles; i++)
       {
         if ( pair_self && fepAtomFlags[angles[i].atom1] != 1) continue;
         if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
         int a1 = angles[i].atom1;
         anglesByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building improper lists.\n");
    
       //  Build the improper lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcImpropers = 0;
       for (i=0; i<numImpropers; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[impropers[i].atom1] != 1) continue;
         byAtomSize[impropers[i].atom1]++;
         numCalcImpropers++;
       }
       for (i=0; i<numAtoms; i++)
       {
         impropersByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         impropersByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numImpropers; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[impropers[i].atom1] != 1) continue;
         int a1 = impropers[i].atom1;
         impropersByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building dihedral lists.\n");
    
       //  Build the dihedral lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcDihedrals = 0;
       for (i=0; i<numDihedrals; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                            && fixedAtomFlags[dihedrals[i].atom2]
                            && fixedAtomFlags[dihedrals[i].atom3]
                            && fixedAtomFlags[dihedrals[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[dihedrals[i].atom1] != 1) continue;
         byAtomSize[dihedrals[i].atom1]++;
         numCalcDihedrals++;
       }
       for (i=0; i<numAtoms; i++)
       {
         dihedralsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         dihedralsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numDihedrals; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                            && fixedAtomFlags[dihedrals[i].atom2]
                            && fixedAtomFlags[dihedrals[i].atom3]
                            && fixedAtomFlags[dihedrals[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[dihedrals[i].atom1] != 1) continue;
         int a1 = dihedrals[i].atom1;
         dihedralsByAtom[a1][byAtomSize[a1]++] = i;
       }
    
       DebugM(3,"Building crossterm lists.\n");
    
       //  Build the crossterm lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcCrossterms = 0;
       for (i=0; i<numCrossterms; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[crossterms[i].atom1]
                            && fixedAtomFlags[crossterms[i].atom2]
                            && fixedAtomFlags[crossterms[i].atom3]
                            && fixedAtomFlags[crossterms[i].atom4]
                            && fixedAtomFlags[crossterms[i].atom5]
                            && fixedAtomFlags[crossterms[i].atom6]
                            && fixedAtomFlags[crossterms[i].atom7]
                            && fixedAtomFlags[crossterms[i].atom8] ) continue;
         if ( pair_self && fepAtomFlags[crossterms[i].atom1] != 1) continue;
         byAtomSize[crossterms[i].atom1]++;
         numCalcCrossterms++;
       }
       for (i=0; i<numAtoms; i++)
       {
         crosstermsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         crosstermsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numCrossterms; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[crossterms[i].atom1]
                            && fixedAtomFlags[crossterms[i].atom2]
                            && fixedAtomFlags[crossterms[i].atom3]
                            && fixedAtomFlags[crossterms[i].atom4]
                            && fixedAtomFlags[crossterms[i].atom5]
                            && fixedAtomFlags[crossterms[i].atom6]
                            && fixedAtomFlags[crossterms[i].atom7]
                            && fixedAtomFlags[crossterms[i].atom8] ) continue;
         if ( pair_self && fepAtomFlags[crossterms[i].atom1] != 1) continue;
         int a1 = crossterms[i].atom1;
         crosstermsByAtom[a1][byAtomSize[a1]++] = i;
       }

       // DRUDE: init lphostIndexes array
       if (simParams->lonepairs) {
         // allocate lone pair host index array only if we need it!
         DebugM(3,"Initializing lone pair host index array.\n");
         lphostIndexes = new int32[numAtoms];
         for (i = 0;  i < numAtoms;  i++) {
           lphostIndexes[i] = -1;
         }
         for (i = 0;  i < numLphosts;  i++) {
           int32 index = lphosts[i].atom1;
           lphostIndexes[index] = i;
         }
       }
       // DRUDE
    
       // JLai
       DebugM(3,"Building gromacsPair lists.\n");
    
       //  Build the gromacsPair lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcLJPair = 0;
       for (i=0; i<numLJPair; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[gromacsPair[i].atom1]
                            && fixedAtomFlags[gromacsPair[i].atom2] ) continue;
         if ( pair_self && fepAtomFlags[gromacsPair[i].atom1] != 1) continue;
         byAtomSize[gromacsPair[i].atom1]++;
         numCalcLJPair++;
       }
       for (i=0; i<numAtoms; i++)
       {
         gromacsPairByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         gromacsPairByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numLJPair; i++)
       {
	   if ( numFixedAtoms && fixedAtomFlags[gromacsPair[i].atom1]
		&& fixedAtomFlags[gromacsPair[i].atom2] ) continue;
	   if ( pair_self && fepAtomFlags[gromacsPair[i].atom1] != 1) continue;
	   int a1 = gromacsPair[i].atom1;
	   gromacsPairByAtom[a1][byAtomSize[a1]++] = i;
	   }

       // End of JLai

       DebugM(3,"Building exclusion data.\n");
    
       //  Build the arrays of exclusions for each atom
       if (! simParams->qmForcesOn)
       build_exclusions();

       //  Remove temporary structures
       delete [] bondsWithAtom;  bondsWithAtom = 0;
       delete tmpArena;  tmpArena = 0;

       if (exclusions != NULL)
      delete [] exclusions;

       // 1-4 exclusions which are also fully excluded were eliminated by hash table
       numTotalExclusions = exclusionSet.size();
       if ( ! CkMyPe() ) {
         iout << iINFO << "ADDED " << (numTotalExclusions - numExclusions) << " IMPLICIT EXCLUSIONS\n" << endi;
       }
       exclusions = new Exclusion[numTotalExclusions];
       UniqueSetIter<Exclusion> exclIter(exclusionSet);
       for ( exclIter=exclIter.begin(),i=0; exclIter != exclIter.end(); exclIter++,i++ )
       {
         exclusions[i] = *exclIter;
       }
       // Free exclusionSet storage
       // exclusionSet.clear(1);
       exclusionSet.clear();

       DebugM(3,"Building exclusion lists.\n");
    
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcExclusions = 0;
       numCalcFullExclusions = 0;
       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         byAtomSize[exclusions[i].atom1]++;
         numCalcExclusions++;
         if ( ! exclusions[i].modified ) numCalcFullExclusions++;
       }

       for (i=0; i<numAtoms; i++)
       {
         exclusionsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         exclusionsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         int a1 = exclusions[i].atom1;
         exclusionsByAtom[a1][byAtomSize[a1]++] = i;
       }

       int32 *byAtomSize2 = new int32[numAtoms];

       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
         byAtomSize2[i] = 0;
       }

       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         if ( exclusions[i].modified ) {
           byAtomSize2[exclusions[i].atom1]++;
           byAtomSize2[exclusions[i].atom2]++;
         } else {
           byAtomSize[exclusions[i].atom1]++;
           byAtomSize[exclusions[i].atom2]++;
         }
       }

       for (i=0; i<numAtoms; i++)
       {
         fullExclusionsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         fullExclusionsByAtom[i][0] = 0;
         modExclusionsByAtom[i] = arena->getNewArray(byAtomSize2[i]+1);
         modExclusionsByAtom[i][0] = 0;
       }

       for (i=0; i<numTotalExclusions; i++)
       {
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         int32 *l1, *l2;
         if ( exclusions[i].modified ) {
           l1 = modExclusionsByAtom[a1];
           l2 = modExclusionsByAtom[a2];
         } else {
           l1 = fullExclusionsByAtom[a1];
           l2 = fullExclusionsByAtom[a2];
         }
         l1[++(*l1)] = a2;
         l2[++(*l2)] = a1;
       }

       if ( ! CkMyPe() && simParams->printExclusions ) {
         for (i=0; i<numAtoms; i++) {
           int32 *lf = fullExclusionsByAtom[i];
           iout << "EXCL " << i << " FULL";
           int nf = *(lf++);
           for ( int j = 0; j < nf; ++j ) {
             iout << " " << *(lf++);
           }
           iout << "\n";
           int32 *lm = modExclusionsByAtom[i];
           iout << "EXCL " << i << " MOD";
           int nm = *(lm++);
           for ( int j = 0; j < nm; ++j ) {
             iout << " " << *(lm++);
           }
           iout << "\n" << endi;
         }
       }

       // DRUDE
       if (is_drude_psf || simParams->drudeOn) {

         // build Thole (screened Coulomb) correction terms;
         // they are constructed implicitly from exclusions

         // free the previous Thole array if already allocated
         if (tholes != NULL) delete[] tholes;
         numTholes = 0;

         // count the number of Thole terms
         for (i = 0;  i < numTotalExclusions;  i++) {
           /* skip over the modified exclusions */
           if (exclusions[i].modified) continue;
           int a1 = exclusions[i].atom1;
           int a2 = exclusions[i].atom2;
           if (a2 < numAtoms-1 && is_drude(a1+1) && is_drude(a2+1)) {
             numTholes++;
           }
         }

         // allocate space for Thole terms
         if (numTholes != 0) tholes = new Thole[numTholes];
         else tholes = NULL;
         int nt = 0;

         Real c = COULOMB*simParams->nonbondedScaling/simParams->dielectric;

         // store Thole terms
         for (i = 0;  i < numTotalExclusions;  i++) {
           /* skip over the modified exclusions */
           if (exclusions[i].modified) continue;
           int a1 = exclusions[i].atom1;
           int a2 = exclusions[i].atom2;
           // exclusions are stored with a1 < a2
           if (a2 < numAtoms-1 && is_drude(a1+1) && is_drude(a2+1)) {
             Real thsum = drudeConsts[a1].thole + drudeConsts[a2].thole;
             Real aprod = drudeConsts[a1].alpha * drudeConsts[a2].alpha;
             // guard against having alpha==0
             Real apower = (aprod <= 0 ? 0 : powf(aprod, -1.f/6));
             tholes[nt].atom1 = a1;
             tholes[nt].atom2 = a1+1;
             tholes[nt].atom3 = a2;
             tholes[nt].atom4 = a2+1;
             tholes[nt].aa = apower * thsum;
             tholes[nt].qq = c * atoms[a1+1].charge * atoms[a2+1].charge;
             nt++;
           }
         }

         // build Thole lists by atom
         DebugM(3, "Building Thole correction term lists.\n");
         tholesByAtom = new int32 *[numAtoms];

         for (i = 0;  i < numAtoms;  i++) {
           byAtomSize[i] = 0;
         }
         numCalcTholes = 0;
         for (i = 0;  i < numTholes;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[tholes[i].atom1]
                              && fixedAtomFlags[tholes[i].atom2]
                              && fixedAtomFlags[tholes[i].atom3]
                              && fixedAtomFlags[tholes[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[tholes[i].atom1] != 1) continue;
           byAtomSize[tholes[i].atom1]++;
           numCalcTholes++;
         }
         for (i = 0;  i < numAtoms;  i++) {
           tholesByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
           tholesByAtom[i][byAtomSize[i]] = -1;
           byAtomSize[i] = 0;
         }
         for (i = 0;  i < numTholes;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[tholes[i].atom1]
                              && fixedAtomFlags[tholes[i].atom2]
                              && fixedAtomFlags[tholes[i].atom3]
                              && fixedAtomFlags[tholes[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[tholes[i].atom1] != 1) continue;
           int a1 = tholes[i].atom1;
           tholesByAtom[a1][byAtomSize[a1]++] = i;
         }

         // build anisotropic lists by atom
         DebugM(3, "Building anisotropic term lists.\n");
         anisosByAtom = new int32 *[numAtoms];

         for (i = 0;  i < numAtoms;  i++) {
           byAtomSize[i] = 0;
         }
         numCalcAnisos = 0;
         for (i = 0;  i < numAnisos;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[anisos[i].atom1]
                              && fixedAtomFlags[anisos[i].atom2]
                              && fixedAtomFlags[anisos[i].atom3]
                              && fixedAtomFlags[anisos[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[anisos[i].atom1] != 1) continue;
           byAtomSize[anisos[i].atom1]++;
           numCalcAnisos++;
         }
         for (i = 0;  i < numAtoms;  i++) {
           anisosByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
           anisosByAtom[i][byAtomSize[i]] = -1;
           byAtomSize[i] = 0;
         }
         for (i = 0;  i < numAnisos;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[anisos[i].atom1]
                              && fixedAtomFlags[anisos[i].atom2]
                              && fixedAtomFlags[anisos[i].atom3]
                              && fixedAtomFlags[anisos[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[anisos[i].atom1] != 1) continue;
           int a1 = anisos[i].atom1;
           anisosByAtom[a1][byAtomSize[a1]++] = i;
         }

       }
       // DRUDE

       delete [] byAtomSize;  byAtomSize = 0;
       delete [] byAtomSize2;  byAtomSize2 = 0;


       //  Allocate an array to hold the exclusions for each atom
       all_exclusions = new ExclusionCheck[numAtoms];

       for (i=0; i<numAtoms; i++)
       {
         all_exclusions[i].min = numAtoms;
         all_exclusions[i].max = -1;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         // first atom should alway have lower number!
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         if ( all_exclusions[a1].min > a2 ) all_exclusions[a1].min = a2;
         if ( all_exclusions[a2].min > a1 ) all_exclusions[a2].min = a1;
         if ( a2 > all_exclusions[a1].max ) all_exclusions[a1].max = a2;
         if ( a1 > all_exclusions[a2].max ) all_exclusions[a2].max = a1;
       }

       // build array of all full exclusions for water etc.
       int maxDenseAllFull = 0;
       int numDenseAllFull = 0;
       for (i=0; i<numAtoms; i++) {
         int iInMiddle = ( i < all_exclusions[i].max &&
                       i > all_exclusions[i].min ) ? 1 : 0;
         int s = all_exclusions[i].max - all_exclusions[i].min + 1;
         if ( s == fullExclusionsByAtom[i][0] + iInMiddle ) {
            if ( s > maxDenseAllFull ) maxDenseAllFull = s;
            all_exclusions[i].flags = (char*)-1;  // shared array
         } else {
            all_exclusions[i].flags = 0;  // individual array
         }
       }
       char *denseFullArray = exclArena->getNewArray(maxDenseAllFull);
       for ( i=0; i<maxDenseAllFull; ++i ) denseFullArray[i] = EXCHCK_FULL;

       int exclmem = maxDenseAllFull;
       int maxExclusionFlags = simParams->maxExclusionFlags;
       for (i=0; i<numAtoms; i++) {
         int s = all_exclusions[i].max - all_exclusions[i].min + 1;
         if ( all_exclusions[i].max != -1 ) {
           if ( all_exclusions[i].flags ) {
             all_exclusions[i].flags = denseFullArray;
             ++numDenseAllFull;
           } else if ( s < maxExclusionFlags ) {
             char *f = all_exclusions[i].flags = exclArena->getNewArray(s);
             for ( int k=0; k<s; ++k ) f[k] = 0;
             exclmem += s;
           } else {
             all_exclusions[i].flags = 0;  // need to build on the fly
           }
         } else {
           all_exclusions[i].flags = (char*)-1; // should never dereference
         }
       }
       if ( 0 ) {
         iout << iINFO << numTotalExclusions << " exclusions consume "
            << exclmem << " bytes.\n" << endi;
         iout << iINFO << numDenseAllFull
            << " atoms sharing one array.\n" << endi;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         if ( exclusions[i].modified ) {
           if ( all_exclusions[a1].flags )
             all_exclusions[a1].flags[a2-all_exclusions[a1].min] = EXCHCK_MOD;
           if ( all_exclusions[a2].flags )
             all_exclusions[a2].flags[a1-all_exclusions[a2].min] = EXCHCK_MOD;
         } else {
           if ( all_exclusions[a1].flags )
             all_exclusions[a1].flags[a2-all_exclusions[a1].min] = EXCHCK_FULL;
           if ( all_exclusions[a2].flags )
             all_exclusions[a2].flags[a1-all_exclusions[a2].min] = EXCHCK_FULL;
         }
       }
    }

    /*    END OF FUNCTION build_lists_by_atom    */

    /****************************************************************/
    /*                */
    /*      FUNCTION build_exclusions    */
    /*                */
    /*  This function builds a list of all the exlcusions       */
    /*  atoms.  These lists include explicit exclusions as well as  */
    /*  exclusions that are calculated based on the bonded structure*/
    /*  and the exclusion flag.  For each pair of atoms that are    */
    /*  excluded, the larger of the 2 atom indexes is stored in the */
    /*  array of the smaller index.  All the arrays are not sorted. */
    /*  Then to determine if two atoms have an exclusion, a linear  */
    /*  search is done on the array of the atom with the smaller    */
    /*  index for the larger index.          */
    /*  If the exclusion policy is set to scaled1-4, there are  */
    /*  actually two lists built.  One contains the pairs of atoms  */
    /*  that are to be exlcuded (i.e., explicit exclusions, 1-2,    */
    /*  and 1-3 interactions) and the other contains just the 1-4   */
    /*  interactions, since they will need to be determined   */
    /*  independantly of the other exclusions.      */
    /*                */
    /****************************************************************/

    void Molecule::build_exclusions()
    {
      register int i;          //  Loop counter
      ExclusionSettings exclude_flag;    //  Exclusion policy

      exclude_flag = simParams->exclude;

      //  Go through the explicit exclusions and add them to the arrays
      for (i=0; i<numExclusions; i++)
      {
        exclusionSet.add(exclusions[i]);
      }

      // If this is AMBER force field, and readExclusions is TRUE,
      // then all the exclusions were read from parm file, and we
      // shouldn't generate any of them.
      if (!simParams->amberOn || !simParams->readExclusions)
      { //  Now calculate the bonded exlcusions based on the exclusion policy
        switch (exclude_flag)
        {
         case NONE:
           break;
         case ONETWO:
           build12excl();
           break;
          case ONETHREE:
            build12excl();
            build13excl();
            break;
          case ONEFOUR:
            build12excl();
            build13excl();
            build14excl(0);
            break;
          case SCALED14:
            build12excl();
            build13excl();
            build14excl(1);
            break;
        }
      }

      stripFepExcl();

      // DRUDE
      if (is_lonepairs_psf) {
        build_inherited_excl(SCALED14 == exclude_flag);
      }
    }
    /*      END OF FUNCTION build_exclusions    */


    // Extend exclusions for the Drude model.  The Drude model is generally
    // used with the 1-3 exclusion policy, although the code below also
    // supports the 1-2 exclusion policy.  The use of light (or massless)
    // pseudo-atoms requires the introduction of extra exclusions.
    //
    // Here is the algorithm for determining Drude model exclusions:
    // (1)  Each Drude particle and each lone pair has a single parent atom.
    //      The parent atom must be a heavy atom.
    // (2)  Each Drude particle and lone pair inherit the exclusions of its
    //      parent atom.
    // (3)  If two heavy atoms are excluded and they both have either a
    //      Drude particle or a lone pair, the these light (or massless)
    //      particles are also excluded from interacting with each other.
    void Molecule::build_inherited_excl(int modified) {
      ExclusionSettings exclude_flag = simParams->exclude;
      int32 *bond1, *bond2, *bond3, *bond4, *bond5;
      int32 i, j, mid1, mid2, mid3, mid4;

      // validate that each Drude or lone pair particle
      // has a unique parent that is a heavy atom
      for (i = 0;  i < numAtoms;  i++) {

        if (!is_drude(i) && !is_lp(i)) continue;
        // make sure that i is either Drude or LP

        // find parent (heavy) atom of particle i
        bond1 = bondsWithAtom[i];

        if (-1 == *bond1) {  // i must have one bond
          char err_msg[512];
          const char *idescrip = (is_drude(i) ? "DRUDE" : "LONE PAIR");
          sprintf(err_msg, "FOUND ISOLATED %s PARTICLE %d", idescrip, i+1);
          NAMD_die(err_msg);
        }
        if (-1 != *(bond1+1)) {  // and only one bond
          char err_msg[512];
          const char *idescrip = (is_drude(i) ? "DRUDE" : "LONE PAIR");
          sprintf(err_msg, "FOUND MULTIPLY LINKED %s PARTICLE %d",
              idescrip, i+1);
          NAMD_die(err_msg);
        }

        // mid1 is parent of particle i
        mid1 = bonds[*bond1].atom1;
        if (mid1 == i) mid1 = bonds[*bond1].atom2;

        // make sure that mid1 is a heavy atom
        if (is_drude(mid1) || is_lp(mid1) || is_hydrogen(mid1)) {
          char err_msg[512];
          const char *idescrip = (is_drude(i) ? "DRUDE" : "LONE PAIR");
          sprintf(err_msg, "PARENT ATOM %d of %s PARTICLE %d "
              "IS NOT HEAVY ATOM", mid1+1, idescrip, i+1);
          NAMD_die(err_msg);
        }

        if (exclude_flag == NONE) {
          // add (i,mid1) as an exclusion
          if (i < mid1) {
            exclusionSet.add(Exclusion(i, mid1));
          }
          else {
            exclusionSet.add(Exclusion(mid1, i));
          }

          // also exclude any Drude particles or LPs bonded to mid1
          bond2 = bondsWithAtom[mid1];
          while (*bond2 != -1) {
            j = bonds[*bond2].atom1;
            if ((is_drude(j) || is_lp(j)) && j != mid1) {
              if      (i < j) exclusionSet.add(Exclusion(i, j));
              else if (j < i) exclusionSet.add(Exclusion(j, i));
            }
            j = bonds[*bond2].atom2;
            if ((is_drude(j) || is_lp(j)) && j != mid1) {
              if      (i < j) exclusionSet.add(Exclusion(i, j));
              else if (j < i) exclusionSet.add(Exclusion(j, i));
            }
            bond2++;
          }
        }
        else {  // if ONETWO or ONETHREE or ONEFOUR or SCALED14

          // find the next link
          bond2 = bondsWithAtom[mid1];

          // loop through all the bonds connected to atom mid1
          while (*bond2 != -1) {
            if (bonds[*bond2].atom1 == mid1) {
              mid2 = bonds[*bond2].atom2;
            }
            else {
              mid2 = bonds[*bond2].atom1;
            }

            // Make sure that we don't double back to where we started from.
            // Doing so causes strange behavior.
            if (mid2 == i) {
              bond2++;
              continue;
            }

            if (exclude_flag == ONETWO) {
              // add (i,mid2) as an exclusion
              if (i < mid2) {
                exclusionSet.add(Exclusion(i, mid2));
              }
              else {
                exclusionSet.add(Exclusion(mid2, i));
              }

              // also exclude any Drude particles or LPs bonded to mid2
              bond3 = bondsWithAtom[mid2];
              while (*bond3 != -1) {
                j = bonds[*bond3].atom1;
                if ((is_drude(j) || is_lp(j)) && j != mid2) {
                  if      (i < j) exclusionSet.add(Exclusion(i, j));
                  else if (j < i) exclusionSet.add(Exclusion(j, i));
                }
                j = bonds[*bond3].atom2;
                if ((is_drude(j) || is_lp(j)) && j != mid2) {
                  if      (i < j) exclusionSet.add(Exclusion(i, j));
                  else if (j < i) exclusionSet.add(Exclusion(j, i));
                }
                bond3++;
              }
            }
            else { // if ONETHREE or ONEFOUR or SCALED14

              // find the next link
              bond3 = bondsWithAtom[mid2];

              // loop through all the bonds connected to mid2
              while (*bond3 != -1) {

                if (bonds[*bond3].atom1 == mid2) {
                  mid3 = bonds[*bond3].atom2;
                }
                else {
                  mid3 = bonds[*bond3].atom1;
                }

                // Make sure we don't double back to where we started.
                // Doing so causes strange behavior.
                if (mid3 == mid1) {
                  bond3++;
                  continue;
                }

                // add (i,mid3) as an exclusion
                if (i < mid3) {
                  exclusionSet.add(Exclusion(i, mid3));
                }
                else if (mid3 < i) {
                  exclusionSet.add(Exclusion(mid3, i));
                }

                if (exclude_flag == ONETHREE) {
                  // also exclude any Drude particles or LPs bonded to mid3
                  bond4 = bondsWithAtom[mid3];
                  while (*bond4 != -1) {
                    j = bonds[*bond4].atom1;
                    if ((is_drude(j) || is_lp(j)) && j != mid3) {
                      if      (i < j) exclusionSet.add(Exclusion(i, j));
                      else if (j < i) exclusionSet.add(Exclusion(j, i));
                    }
                    j = bonds[*bond4].atom2;
                    if ((is_drude(j) || is_lp(j)) && j != mid3) {
                      if      (i < j) exclusionSet.add(Exclusion(i, j));
                      else if (j < i) exclusionSet.add(Exclusion(j, i));
                    }
                    bond4++;
                  }
                }
                else { // if ONEFOUR or SCALED14

                  // find next link
                  bond4 = bondsWithAtom[mid3];

                  // loop through all the bonds connected to mid3
                  while (*bond4 != -1) {

                    if (bonds[*bond4].atom1 == mid3) {
                      mid4 = bonds[*bond4].atom2;
                    }
                    else {
                      mid4 = bonds[*bond4].atom1;
                    }

                    // Make sure we don't double back to where we started.
                    // Doing so causes strange behavior.
                    if (mid4 == mid2) {
                      bond4++;
                      continue;
                    }

                    if (is_drude(mid4) || is_lp(mid4)) {
                      // (i,mid4) is 1-3 excl
                      if (i < mid4) {
                        exclusionSet.add(Exclusion(i, mid4));
                      }
                      else if (mid4 < i) {
                        exclusionSet.add(Exclusion(mid4, i));
                      }
                      bond4++;
                      continue;
                    }

                    // (mid1,mid4) is an existing heavy atom exclusion
                    // if we have modified 1-4 exclusions, make sure
                    // that (mid1,mid4) is modified 1-4 exclusion
                    // rather than something closer due to a ring
                    int modi = modified;
                    if (modified) {
                      int amin = (mid1 < mid4 ? mid1 : mid4);
                      int amax = (mid1 >= mid4 ? mid1 : mid4);
                      Exclusion *pe = exclusionSet.find(Exclusion(amin,amax));
                      if (pe==0) {
                        // since there is not an existing exclusion
                        // between (mid1,mid4), don't inherit!
                        bond4++;
                        continue;
                      }
                      modi = pe->modified;
                    }

                    if (i < mid4) {
                      exclusionSet.add(Exclusion(i, mid4, modi));
                    }
                    else if (mid4 < i) {
                      exclusionSet.add(Exclusion(mid4, i, modi));                    
                    }

                    // also exclude any Drude particles or LPs bonded to mid4
                    // using the "modi" setting of (mid1,mid4) exclusion
                    bond5 = bondsWithAtom[mid4];
                    while (*bond5 != -1) {
                      j = bonds[*bond5].atom1;
                      if ((is_drude(j) || is_lp(j)) && j != mid4) {
                        if      (i<j) exclusionSet.add(Exclusion(i,j,modi));
                        else if (j<i) exclusionSet.add(Exclusion(j,i,modi));
                      }
                      j = bonds[*bond5].atom2;
                      if ((is_drude(j) || is_lp(j)) && j != mid4) {
                        if      (i<j) exclusionSet.add(Exclusion(i,j,modi));
                        else if (j<i) exclusionSet.add(Exclusion(j,i,modi));
                      }
                      bond5++;
                    }
                    ++bond4;
                  } // while bond4

                } // else (if ONEFOUR or SCALED14)

                ++bond3;
              } // while bond3

            } // else (if ONETHREE or ONEFOUR or SCALED14)
           
            ++bond2;
          } // while bond2

        } // else (if ONETWO or ONETHREE or ONEFOUR or SCALED14)

      } // for i
    } 
    // DRUDE


    /************************************************************************/
    /*                  */
    /*      FUNCTION build12excl        */
    /*                  */
    /************************************************************************/

    void Molecule::build12excl(void)       
    {
       int32 *current_val;  //  Current value to check
       register int i;    //  Loop counter to loop through all atoms
       
       //  Loop through all the atoms marking the bonded interactions for each one
       for (i=0; i<numAtoms; i++)
       {
      current_val = bondsWithAtom[i];
       
      //  Loop through all the bonds for this atom
      while (*current_val != -1)
      {
         if (bonds[*current_val].atom1 == i)
         {
      if (i<bonds[*current_val].atom2)
      {
         exclusionSet.add(Exclusion(i,bonds[*current_val].atom2));
      }
         }
         else
         {
      if (i<bonds[*current_val].atom1)
      {
         exclusionSet.add(Exclusion(i,bonds[*current_val].atom1));
      }
         }
    
         ++current_val;
      }
       }
    }
    /*      END OF FUNCTION build12excl      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build13excl        */
    /*                  */
    /************************************************************************/

    void Molecule::build13excl(void)       
    {
       int32 *bond1, *bond2;  //  The two bonds being checked
       int middle_atom;  //  Common third atom
       register int i;    //  Loop counter to loop through all atoms
       
       //  Loop through all the atoms looking at the bonded connections
       //  for each one
       for (i=0; i<numAtoms; i++)
       {
       bond1 = bondsWithAtom[i];
       
       //  Loop through all the bonds directly connect to atom i
       while (*bond1 != -1)
       {
        if (bonds[*bond1].atom1 == i)
        {
          middle_atom=bonds[*bond1].atom2;
        }
        else
        {
          middle_atom=bonds[*bond1].atom1;
        }

        bond2 = bondsWithAtom[middle_atom];

        //  Now loop through all the bonds connect to the
        //  middle atom
        while (*bond2 != -1)
        {
          if (bonds[*bond2].atom1 == middle_atom)
          {
            if (i < bonds[*bond2].atom2)
            {
              exclusionSet.add(Exclusion(i,bonds[*bond2].atom2));
            }
          }
          else
          {
            if (i < bonds[*bond2].atom1)
            {
              exclusionSet.add(Exclusion(i,bonds[*bond2].atom1));
            }
          }

          ++bond2;
        }

        ++bond1;
      }
       }
    }
    /*      END OF FUNCTION build13excl      */

    /************************************************************************/
    /*                  */
    /*        FUNCTION build14excl      */
    /*                  */
    /************************************************************************/


    void Molecule::build14excl(int modified)       
    {
       int32 *bond1, *bond2, *bond3;  //  The two bonds being checked
       int mid1, mid2;    //  Middle atoms
       register int i;      //  Counter to loop through all atoms
       
       //  Loop through all the atoms
       for (i=0; i<numAtoms; i++)
       {  
         if (is_drude(i) || is_lp(i)) continue;  // skip Drude and LP for now

      // Get all the bonds connect directly to atom i
      bond1 = bondsWithAtom[i];
       
      while (*bond1 != -1)
      {
        if (bonds[*bond1].atom1 == i)
        {
          mid1=bonds[*bond1].atom2;
        }
        else
        {
          mid1=bonds[*bond1].atom1;
        }

        bond2 = bondsWithAtom[mid1];

        //  Loop through all the bonds connected to atom mid1
        while (*bond2 != -1)
        {
          if (bonds[*bond2].atom1 == mid1)
          {
            mid2 = bonds[*bond2].atom2;
          }
          else
          {
            mid2 = bonds[*bond2].atom1;
          }

          //  Make sure that we don't double back to where
          //  we started from.  This causes strange behavior.
          //  Trust me, I've been there . . .
          if (mid2 == i)
          {
            ++bond2;
            continue;
          }

          bond3=bondsWithAtom[mid2];

          //  Loop through all the bonds connected to mid2
          while (*bond3 != -1)
          {
            if (bonds[*bond3].atom1 == mid2)
            {
              int j = bonds[*bond3].atom2;
              //  Make sure that we don't double back to where
              //  we started from.  This causes strange behavior.
              //  Trust me, I've been there . . .
              //  I added this!!!  Why wasn't it there before?  -JCP
              if (j != mid1)
              if (i < j && !is_drude(j) && !is_lp(j))  // skip Drude and LP
              {
                 exclusionSet.add(Exclusion(i,j,modified));
              }
            }
            else
            {
              int j = bonds[*bond3].atom1;
              //  Make sure that we don't double back to where
              //  we started from.  This causes strange behavior.
              //  Trust me, I've been there . . .
              //  I added this!!!  Why wasn't it there before?  -JCP
              if (j != mid1)
              if (i < j && !is_drude(j) && !is_lp(j))  // skip Drude and LP
              {
                 exclusionSet.add(Exclusion(i,j,modified));
              }
            }

            ++bond3;
          }

          ++bond2;
        }
    
        ++bond1;
      }
       }
    }
    /*      END OF FUNCTION build14excl      */


    /************************************************************************/
    /*                                                                      */
    /*        FUNCTION stripFepExcl                                         */
    /*                                                                      */
    /************************************************************************/
  void Molecule::stripFepExcl(void)
  {   
    UniqueSet<Exclusion> fepExclusionSet;
    UniqueSetIter<Exclusion> exclIter(exclusionSet);

    if ( simParams->alchOn || simParams->lesOn ) {
       for ( exclIter=exclIter.begin(); exclIter != exclIter.end(); exclIter++ )
       {
         int t1 = get_fep_type(exclIter->atom1);
         int t2 = get_fep_type(exclIter->atom2);
         if ( t1 && t2 && t1 != t2 ) {
           fepExclusionSet.add(*exclIter);
         }
       }
    } else if ( simParams->pairInteractionOn ) {
      for ( exclIter=exclIter.begin(); exclIter != exclIter.end(); exclIter++ )
      {
        int ifep_type = get_fep_type(exclIter->atom1);
        int jfep_type = get_fep_type(exclIter->atom2);
        if ( simParams->pairInteractionSelf ) {
          // for pair-self, both atoms must be in group 1.
          if (ifep_type != 1 || jfep_type != 1) {
            fepExclusionSet.add(*exclIter);
          }
        } else {

          // for pair, must have one from each group.
          if (!(ifep_type == 1 && jfep_type == 2) &&
              !(ifep_type == 2 && jfep_type == 1)) {
            fepExclusionSet.add(*exclIter);
          }
        }
       }
    }

    UniqueSetIter<Exclusion> fepIter(fepExclusionSet);
    for ( fepIter=fepIter.begin(); fepIter != fepIter.end(); fepIter++ )
    {
      exclusionSet.del(*fepIter);
    }
  }
    /*      END OF FUNCTION stripFepExcl      */

#else

//===Memory optimized version of functions that read Molecule file===//
void Molecule::read_mol_signatures(char *fname, Parameters *params, ConfigList *cfgList){
    FILE *psf_file;    //  pointer to .psf file
    int ret_code;    //  ret_code from NAMD_read_line calls
    char buffer[512];

    
    if ( (psf_file = Fopen(fname, "r")) == NULL)
    {
        char err_msg[512];
        sprintf(err_msg, "UNABLE TO OPEN THE COMPRESSED .psf FILE %s", fname);
        NAMD_die(err_msg);
    }

    char strBuf[12];

    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "FORMAT VERSION")) {
        NAMD_die("The compressed psf file format is incorrect, please re-generate!\n");
    }
    float psfVer = 0.0f;
    sscanf(buffer, "FORMAT VERSION: %f\n", &psfVer);
    if(fabs(psfVer - COMPRESSED_PSF_VER)>1e-6) {
        NAMD_die("The compressed psf file format is incorrect, please re-generate!\n");
    }

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NSEGMENTNAMES"))
        NAMD_die("UNABLE TO FIND NSEGMENTNAMES");
    sscanf(buffer, "%d", &segNamePoolSize);
#if 0
    if(segNamePoolSize!=0)
        segNamePool = new char *[segNamePoolSize];
    for(int i=0; i<segNamePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        segNamePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(segNamePool[i], strBuf);
    }
#else
    for(int i=0; i<segNamePoolSize; i++) NAMD_read_line(psf_file, buffer);
#endif
    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NRESIDUENAMES"))
        NAMD_die("UNABLE TO FIND NRESIDUENAMES");
    sscanf(buffer, "%d", &resNamePoolSize);
#if 0
    if(resNamePoolSize!=0)
        resNamePool = new char *[resNamePoolSize];
    for(int i=0; i<resNamePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        resNamePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(resNamePool[i], strBuf);
    }
#else
    for(int i=0; i<resNamePoolSize; i++) NAMD_read_line(psf_file, buffer);
#endif

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NATOMNAMES"))
        NAMD_die("UNABLE TO FIND NATOMNAMES");
    sscanf(buffer, "%d", &atomNamePoolSize);
    if(atomNamePoolSize!=0)
        atomNamePool = new char *[atomNamePoolSize];
    for(int i=0; i<atomNamePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        atomNamePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(atomNamePool[i], strBuf);
    }
    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NATOMTYPES"))
        NAMD_die("UNABLE TO FIND NATOMTYPES");
    sscanf(buffer, "%d", &atomTypePoolSize);
#if 0
    if(atomTypePoolSize!=0)
        atomTypePool = new char *[atomTypePoolSize];
    for(int i=0; i<atomTypePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        atomTypePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(atomTypePool[i], strBuf);
    }
#else
    for(int i=0; i<atomTypePoolSize; i++) NAMD_read_line(psf_file, buffer);
#endif
    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NCHARGES"))
        NAMD_die("UNABLE TO FIND NCHARGES");
    sscanf(buffer, "%d", &chargePoolSize);
    if(chargePoolSize!=0)
        atomChargePool = new Real[chargePoolSize];
    for(int i=0; i<chargePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%f", atomChargePool+i);
    }
    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NMASSES"))
        NAMD_die("UNABLE TO FIND NMASSES");
    sscanf(buffer, "%d", &massPoolSize);
    if(massPoolSize!=0)
        atomMassPool = new Real[massPoolSize];
    for(int i=0; i<massPoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%f", atomMassPool+i);
    }
    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "ATOMSIGS"))
        NAMD_die("UNABLE TO FIND ATOMSIGS");
    sscanf(buffer, "%d", &atomSigPoolSize);
    atomSigPool = new AtomSignature[atomSigPoolSize];
    int typeCnt;
    int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    int tisReal;
    int ttype;
    for(int i=0; i<atomSigPoolSize; i++){
        
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NBONDSIGS"))
            NAMD_die("UNABLE TO FIND NBONDSIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].bondCnt = typeCnt;
            atomSigPool[i].bondSigs = new TupleSignature[typeCnt];
        }
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);
            sscanf(buffer, "%d | %d | %d", &tmp1, &ttype, &tisReal);
            TupleSignature oneSig(1, BOND, (Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            atomSigPool[i].bondSigs[j]=oneSig;
            if(tisReal) numRealBonds++;
        }

        
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NTHETASIGS"))
            NAMD_die("UNABLE TO FIND NTHETASIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].angleCnt = typeCnt;
            atomSigPool[i].angleSigs = new TupleSignature[typeCnt];
        }
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);
            sscanf(buffer, "%d %d | %d | %d", &tmp1, &tmp2, &ttype, &tisReal);
            TupleSignature oneSig(2,ANGLE,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            atomSigPool[i].angleSigs[j] = oneSig;
        }
        
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NPHISIGS"))
            NAMD_die("UNABLE TO FIND NPHISIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].dihedralCnt = typeCnt;
            atomSigPool[i].dihedralSigs = new TupleSignature[typeCnt];
        }
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);
            sscanf(buffer, "%d %d %d | %d | %d", &tmp1, &tmp2, &tmp3, &ttype, &tisReal);
            TupleSignature oneSig(3,DIHEDRAL,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            oneSig.offset[2] = tmp3;
            atomSigPool[i].dihedralSigs[j] = oneSig;
        }
        
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NIMPHISIGS"))
            NAMD_die("UNABLE TO FIND NIMPHISIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].improperCnt = typeCnt;
            atomSigPool[i].improperSigs = new TupleSignature[typeCnt];
        }
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);
            sscanf(buffer, "%d %d %d | %d | %d", &tmp1, &tmp2, &tmp3, &ttype, &tisReal);
            TupleSignature oneSig(3,IMPROPER,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            oneSig.offset[2] = tmp3;
            atomSigPool[i].improperSigs[j] = oneSig;
        }
        
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NCRTERMSIGS"))
            NAMD_die("UNABLE TO FIND NCRTERMSIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].crosstermCnt = typeCnt;
            atomSigPool[i].crosstermSigs = new TupleSignature[typeCnt];
        }
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);
            sscanf(buffer, "%d %d %d %d %d %d %d | %d | %d", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &ttype, &tisReal);
            TupleSignature oneSig(7,CROSSTERM,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            oneSig.offset[2] = tmp3;
            oneSig.offset[3] = tmp4;
            oneSig.offset[4] = tmp5;
            oneSig.offset[5] = tmp6;
            oneSig.offset[6] = tmp7;
            atomSigPool[i].crosstermSigs[j] = oneSig;
        }
    }

    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NEXCLSIGS")){
        NAMD_die("UNABLE TO FIND NEXCLSIGS");
    }
    sscanf(buffer, "%d", &exclSigPoolSize);
    if(exclSigPoolSize>0) exclSigPool = new ExclusionSignature[exclSigPoolSize];
    vector<int> fullExcls;
    vector<int> modExcls;
    for(int i=0; i<exclSigPoolSize; i++){
        int fullExclCnt = NAMD_read_int(psf_file, buffer);
        for(int j=0; j<fullExclCnt; j++)
            fullExcls.push_back(NAMD_read_int(psf_file, buffer));
        int modExclCnt = NAMD_read_int(psf_file, buffer);
        for(int j=0; j<modExclCnt; j++)
            modExcls.push_back(NAMD_read_int(psf_file, buffer));

        
        exclSigPool[i].setOffsets(fullExcls, modExcls);

        fullExcls.clear();
        modExcls.clear();
    }

    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NCLUSTERS")) {
        NAMD_die("UNABLE TO FIND NCLUSTERS");
    }
    sscanf(buffer, "%d", &numClusters);

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NATOM"))
        NAMD_die("UNABLE TO FIND NATOM");
    sscanf(buffer, "%d", &numAtoms);

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NHYDROGENGROUP"))
        NAMD_die("UNABLE TO FIND NHYDROGENGROUP");
    sscanf(buffer, "%d", &numHydrogenGroups);

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "MAXHYDROGENGROUPSIZE"))
        NAMD_die("UNABLE TO FIND MAXHYDROGENGROUPSIZE");
    sscanf(buffer, "%d", &maxHydrogenGroupSize);
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NMIGRATIONGROUP"))
        NAMD_die("UNABLE TO FIND NMIGRATIONGROUP");
    sscanf(buffer, "%d", &numMigrationGroups);
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "MAXMIGRATIONGROUPSIZE"))
        NAMD_die("UNABLE TO FIND MAXMIGRATIONGROUPSIZE");
    sscanf(buffer, "%d", &maxMigrationGroupSize);

    int inputRigidType = -1;
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "RIGIDBONDTYPE"))
      NAMD_die("UNABLE TO FIND RIGIDBONDTYPE");
    sscanf(buffer, "%d", &inputRigidType);
    if(simParams->rigidBonds != RIGID_NONE){
      //check whether the input rigid bond type matches
      if(simParams->rigidBonds != inputRigidType){
        char *tmpstr[]={"RIGID_NONE", "RIGID_ALL", "RIGID_WATER"};
        char errmsg[125];
        sprintf(errmsg, "RIGIDBOND TYPE MISMATCH BETWEEN INPUT (%s) AND CURRENT RUN (%s)", 
                tmpstr[inputRigidType], tmpstr[simParams->rigidBonds]);
        NAMD_die(errmsg);
      }
    }
#if 0
//    int isOccupancyValid, isBFactorValid;
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "OCCUPANCYVALID"))
        NAMD_die("UNABLE TO FIND OCCUPANCYVALID");
    sscanf(buffer, "%d", &isOccupancyValid);
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "TEMPFACTORVALID"))
        NAMD_die("UNABLE TO FIND TEMPFACTORVALID");
    sscanf(buffer, "%d", &isBFactorValid);    
#endif

    //Just reading for the parameters values; extra Bonds, Dihedrals etc.
    //have been taken into account when compressing the molecule object.
    //The actual number of Bonds, Dihedrals etc. will be calculated based
    //on atom signatures.
    if(cfgList && simParams->extraBondsOn)
        build_extra_bonds(params, cfgList->find("extraBondsFile"));

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "DIHEDRALPARAMARRAY"))
        NAMD_die("UNABLE TO FIND DIHEDRALPARAMARRAY");
    for(int i=0; i<params->NumDihedralParams; i++){
        params->dihedral_array[i].multiplicity = NAMD_read_int(psf_file, buffer);
    }
    

    NAMD_read_line(psf_file, buffer); //to read a simple single '\n' line
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "IMPROPERPARAMARRAY"))
        NAMD_die("UNABLE TO FIND IMPROPERPARAMARRAY");
    for(int i=0; i<params->NumImproperParams; i++){
        params->improper_array[i].multiplicity = NAMD_read_int(psf_file, buffer);
    }
    
    Fclose(psf_file);
}

/*
 * The following method is called on every input processors. However, in SMP mode, two 
 * input procs are likely to be inside the same SMP node. Additionally, there's only one 
 * Molecule object per SMP node. Therefore, if there are any assignments to the molecule 
 * object in this function, there's going to be  DATA RACE !!  That's why the calculation of 
 * numTupes(Bonds, Angles, etc) and numExclusions has to be done in ParallelIOMgr, and 
 * then reduce those values.   
 * -Chao Mei
 */
void Molecule::read_binary_atom_info(int fromAtomID, int toAtomID, InputAtomList& inAtoms){
    int numAtomsPar = toAtomID-fromAtomID+1;
    CmiAssert(numAtomsPar > 0);
    CmiAssert(inAtoms.size() == numAtomsPar);
    
    /*
    //all the following vars are not needed as the
    //atom info is loaded into inAtoms.
    atoms = new AtomCstInfo[numAtomsPar];
    atomNames = new AtomNameIdx[numAtomsPar];
    eachAtomMass = new Index[numAtomsPar];
    eachAtomCharge = new Index[numAtomsPar];
    eachAtomSig = new Index[numAtomsPar];
    eachAtomExclSig = new Index[numAtomsPar];
    */
    /*
    atoms = new AtomCstInfo[numAtomsPar];
    atomNames = new AtomNameIdx[numAtomsPar];
    */

    atoms = NULL;
    atomNames = NULL;
    eachAtomMass = NULL;
    eachAtomCharge = NULL;
    eachAtomSig = NULL;
    eachAtomExclSig = NULL;
    clusterSigs = NULL;

    /* HydrogenGroup is not needed here anymore
    hydrogenGroup.resize(numAtomsPar);
    ResidueLookupElem *tmpResLookup = resLookup;
    */
/*
    if (isOccupancyValid) {
        occupancy = new float[numAtomsPar];
    }
    if (isBFactorValid) {
        bfactor = new float[numAtomsPar];
    }
*/
    occupancy = NULL;
    bfactor = NULL;
    
    char *segment_name; 

    //use "fopen" instead of "Fopen" because "stat" call inside Fopen may 
    //fail on some platforms (such as BG/P) for very large file because of 
    //EOVERFLOW, say a 2GB file. -Chao Mei
    FILE *perAtomFile = fopen(simParams->binAtomFile, "rb");
    if (perAtomFile==NULL) {
        char err_msg[512];
        sprintf(err_msg, "UNABLE TO OPEN THE ASSOCIATED PER-ATOM FILE FOR THE COMPRESSED .psf FILE %s", simParams->binAtomFile);
        NAMD_die(err_msg);
    }
    int needFlip = 0;
    int magicNum = COMPRESSED_PSF_MAGICNUM;
    int rMagicNum = COMPRESSED_PSF_MAGICNUM;
    flipNum((char *)&rMagicNum, sizeof(int), 1);
    int fMagicNum;
    fread(&fMagicNum, sizeof(int), 1, perAtomFile);
    if (fMagicNum==magicNum) {
        needFlip = 0;
    } else if (fMagicNum==rMagicNum) {
        needFlip = 1;
    } else {
        char err_msg[512];
        sprintf(err_msg, "THE ASSOCIATED PER-ATOM FILE FOR THE COMPRESSED .psf FILE %s IS CORRUPTED", simParams->binAtomFile);
        NAMD_die(err_msg);
    }

    float verNum =  0.0f;
    fread(&verNum, sizeof(float), 1, perAtomFile);
    if (needFlip) flipNum((char *)&verNum, sizeof(float), 1);
    if (fabs(verNum - COMPRESSED_PSF_VER)>1e-6) {
        char err_msg[512];
        sprintf(err_msg, "THE ASSOCIATED PER-ATOM FILE FOR THE COMPRESSED .psf FILE %s IS INCORRECT, PLEASE RE-GENERATE!\n", simParams->binAtomFile);
        NAMD_die(err_msg);
    }

    int recSize = 0;
    fread(&recSize, sizeof(int), 1, perAtomFile);
    if(needFlip) flipNum((char *)&recSize, sizeof(int), 1);
    if(recSize != sizeof(OutputAtomRecord)){
      char err_msg[512];
      sprintf(err_msg, "THE ASSOCIATED PER-ATOM RECORD SIZE FOR THE COMPRESSED .psf FILE %s IS INCORRECT, PLEASE RE-GENERATE!\n", simParams->binAtomFile);
      NAMD_die(err_msg);
    }
    
    const int BUFELEMS = 32*1024; //32K elems

    //remember to convert to long in case of int overflow!
    int64 startbyte=((int64)fromAtomID)*sizeof(OutputAtomRecord);
#ifdef WIN32
    if ( _fseeki64(perAtomFile,startbyte,SEEK_CUR) )
#else
    if ( fseeko(perAtomFile,startbyte,SEEK_CUR) )
#endif
    {
      char errmsg[512];
      sprintf(errmsg, "Error on seeking binary file %s", simParams->binAtomFile);
      NAMD_err(errmsg);
    }

    //reduce the number of fread calls as file I/O is expensive.
    OutputAtomRecord *elemsBuf = new OutputAtomRecord[BUFELEMS];
    int atomsCnt = numAtomsPar;
    int curIdx=0;
    OutputAtomRecord *oneRec = NULL;
    while(atomsCnt >= BUFELEMS) {
      if ( fread((char *)elemsBuf, sizeof(OutputAtomRecord), BUFELEMS, perAtomFile) != BUFELEMS ) {
        char errmsg[512];
        sprintf(errmsg, "Error on reading binary file %s", simParams->binAtomFile);
        NAMD_err(errmsg);
      }
      oneRec = elemsBuf;
      for(int i=0; i<BUFELEMS; i++, curIdx++, oneRec++) {
        InputAtom *fAtom = &(inAtoms[curIdx]);
        int aid = curIdx+fromAtomID;
        if(needFlip) oneRec->flip();
        load_one_inputatom(aid, oneRec, fAtom);        
      }
      atomsCnt -= BUFELEMS;
    }

    if ( fread(elemsBuf, sizeof(OutputAtomRecord), atomsCnt, perAtomFile) != atomsCnt ) {
      char errmsg[512];
      sprintf(errmsg, "Error on reading binary file %s", simParams->binAtomFile);
      NAMD_err(errmsg);
    }
    oneRec = elemsBuf;    
    for(int i=curIdx; i<numAtomsPar; i++, oneRec++) {
      InputAtom *fAtom = &(inAtoms[i]);
      int aid = i+fromAtomID;
      if(needFlip) oneRec->flip();
      load_one_inputatom(aid,oneRec,fAtom);      
    }

    if ( fclose(perAtomFile) ) {
      char errmsg[512];
      sprintf(errmsg, "Error on closing binary file %s", simParams->binAtomFile);
      NAMD_err(errmsg);
    }

    delete [] elemsBuf;

    //deal with fixed atoms info
    if(simParams->fixedAtomsOn){
        int listIdx=0;
        is_atom_fixed(fromAtomID, &listIdx);
        for(int i=listIdx; i<fixedAtomsSet->size(); i++){
            const AtomSet one = fixedAtomsSet->item(i);
            //set the atoms in this range to be fixed
            int sAtomId = one.aid1>fromAtomID ? one.aid1:fromAtomID;
            int eAtomId = one.aid2>toAtomID? toAtomID:one.aid2;
            for(int j=sAtomId; j<=eAtomId; j++)
                inAtoms[j-fromAtomID].atomFixed = 1;
        }
    }
}

void Molecule::load_one_inputatom(int aid, OutputAtomRecord *one, InputAtom *fAtom){

  char *thisAtomName = NULL;
  fAtom->isValid=true;
  
  //segment_name = segNamePool[sIdx[0]];
  /*
  atomNames[i].resnameIdx = sIdx[1];
  atomNames[i].atomnameIdx = sIdx[2];
  atomNames[i].atomtypeIdx = sIdx[3]; 
  */
  thisAtomName = atomNamePool[one->sSet.atomNameIdx];
         
  fAtom->charge = atomChargePool[one->sSet.chargeIdx];
  fAtom->mass = atomMassPool[one->sSet.massIdx];
  fAtom->sigId = one->iSet.atomSigIdx;
  fAtom->exclId = one->iSet.exclSigIdx;
  fAtom->vdwType = one->sSet.vdw_type;
  
  //atoms[i].vdw_type = sIdx[8];
  
  int residue_number; //for residue number
  residue_number = one->iSet.resID;
  
  /*
  atoms[i].partner = iIdx[2];
  atoms[i].hydrogenList= iIdx[3];
  */
  
  fAtom->id=aid;
  fAtom->atomFixed = 0;
  fAtom->hydList = one->iSet.hydrogenList;
  fAtom->hydrogenGroupSize=one->iSet.atomsInGroup;
  fAtom->GPID=one->iSet.GPID;        
  //fAtom->waterVal=one->waterVal;
  fAtom->migrationGroupSize=one->iSet.atomsInMigrationGroup;
  fAtom->MPID=one->iSet.MPID;
  fAtom->isGP=(fAtom->hydrogenGroupSize ? 1 : 0);
  fAtom->isMP=( fAtom->migrationGroupSize ? 1 : 0 ); 
  
  if(simParams->rigidBonds) {
    fAtom->rigidBondLength = one->fSet.rigidBondLength;
  }else{
    fAtom->rigidBondLength = 0.0;
  }
  
  //Node::Object()->ioMgr->maxAtom=fAtom.id;        
  
  /*if (isOccupancyValid)
      occupancy[i] = tmpf[0];
  if (isBFactorValid)
      bfactor[i] = tmpf[1];*/
  
  /*
  ////TODO: DEAL WITH RESIDUE LOOKUP --CHAOMEI
  if (tmpResLookup) tmpResLookup =
          tmpResLookup->append(segment_name, residue_number, i);
  */
  
  Real thisAtomMass = fAtom->mass;
  
  if ( simParams->ignoreMass ) {
  } else if (thisAtomMass <= 0.05) {
      fAtom->status |= LonepairAtom;
  } else if (thisAtomMass < 1.0) {
      fAtom->status |= DrudeAtom;
  } else if (thisAtomMass <= 3.5) {
      fAtom->status = HydrogenAtom;
  } else if (thisAtomName[0]=='O' &&
             (thisAtomMass >= 14.0) && (thisAtomMass <= 18.0)) {
      fAtom->status = OxygenAtom;
  }
  
  //Move the langevinParam setting which depends on atom's status
  //to the time when each home patch is filled with their atoms
  //(WorkDistrib::fillAtomListForOnePatch so that "langevinParam"
  //could be shared with the "hydVal".
  //--Chao Mei 
}

//Well, the exclusion check signatures could also done on PE0 and
//sent to other processors through send_Molecule/receive_Molecule 
//two procedures.
void Molecule::build_excl_check_signatures(){
   exclChkSigPool = new ExclusionCheck[exclSigPoolSize];
   for(int i=0; i<exclSigPoolSize; i++){
       ExclusionSignature *sig = &exclSigPool[i];
       ExclusionCheck *sigChk = &exclChkSigPool[i];
       if(sig->fullExclCnt){
           if(!sig->modExclCnt){ //only having fullExclusion
               sigChk->min = sig->fullOffset[0];
               sigChk->max = sig->fullOffset[sig->fullExclCnt-1];
           }else{ //have both full and modified exclusion
               int fullMin, fullMax, modMin, modMax;
               
               fullMin = sig->fullOffset[0];
               fullMax = sig->fullOffset[sig->fullExclCnt-1];
           
               modMin = sig->modOffset[0];
               modMax = sig->modOffset[sig->modExclCnt-1];
               
               if(fullMin < modMin)
                   sigChk->min = fullMin;
               else
                   sigChk->min = modMin;
               if(fullMax < modMax)
                   sigChk->max = modMax;
               else
                   sigChk->max = fullMax;
           }        
       }else{
           if(sig->modExclCnt){
               sigChk->min = sig->modOffset[0];
               sigChk->max = sig->modOffset[sig->modExclCnt-1];
           }else{ //both count are 0
               if(CkMyPe()==0)
                   iout << iWARN << "an empty exclusion signature with index "
                     << i << "!\n" << endi;
               continue;
           }
       }           

       sigChk->flags = new char[sigChk->max-sigChk->min+1];
       memset(sigChk->flags, 0, sizeof(char)*(sigChk->max-sigChk->min+1));
       for(int j=0; j<sig->fullExclCnt; j++){
           int dist = sig->fullOffset[j] - sigChk->min;
           sigChk->flags[dist] = EXCHCK_FULL;
       }
       for(int j=0; j<sig->modExclCnt; j++){
           int dist = sig->modOffset[j] - sigChk->min;
           sigChk->flags[dist] = EXCHCK_MOD;
       }
   }
}

/**
 * This function loads the fixed atoms into Molecule object. 
 * The input is a text file, each line specifies a range of 
 * fixed atoms in the format of atomid1[-atomid2]. The atom id
 * starts from 0, not 1! This function should be called only on 
 * the master proc. 
 * -Chao Mei
 */

void Molecule::load_atom_set(StringList *setfile, const char *setname,
	int *numAtomsInSet, AtomSetList **atomsSet) const {
  if(setfile == NULL) {
    char errmsg[128];
    sprintf(errmsg,"The text input file for %s atoms is not found!", setname);
    NAMD_die(errmsg);
  }
  FILE *ifp = fopen(setfile->data, "r");
  
  if(ifp==NULL){
      char errmsg[128];
      sprintf(errmsg, "ERROR IN OPENING %s ATOMS FILE: %s\n", setname, setfile->data);
      NAMD_die(errmsg);
  }

  char oneline[128];
  int numLocalAtoms = 0;
  AtomSet one;
  AtomSetList *localAtomsSet = new AtomSetList();  
  while(1) {
    int ret = NAMD_read_line(ifp, oneline, 128);
    if(ret!=0) break;
    if(NAMD_blank_string(oneline)) continue;
    bool hasDash = false;
    for(int i=0; oneline[i] && i<128; i++){
      if(oneline[i]=='-') {
        hasDash = true;
        break;
      }
    }
    if(hasDash) {
      sscanf(oneline,"%d-%d", &(one.aid1), &(one.aid2));
      if(one.aid1>one.aid2 || one.aid1<0 || one.aid2<0) {
        char errmsg[512];
        sprintf(errmsg, "The input for %s atoms is wrong: %s\n", setname, oneline);
        NAMD_die(errmsg);
      }
      numLocalAtoms += (one.aid2-one.aid1+1);
    }else{
      sscanf(oneline, "%d", &(one.aid1));
      if(one.aid1<0) {
        char errmsg[512];
        sprintf(errmsg, "The input for %s atoms is wrong: %s\n", setname, oneline);
        NAMD_die(errmsg);      
      }
      one.aid2 = one.aid1;
      numLocalAtoms++;
    }
    localAtomsSet->add(one);
  }
  //sort the localAtomsSet for binary search to decide 
  //whether an atom is in the set or not
  std::sort(localAtomsSet->begin(), localAtomsSet->end());  

  *numAtomsInSet = numLocalAtoms;
  *atomsSet = localAtomsSet;
}

void Molecule::load_fixed_atoms(StringList *fixedfile){
  load_atom_set(fixedfile, "FIXED", &numFixedAtoms, &fixedAtomsSet);
}

void Molecule::load_constrained_atoms(StringList *constrainedfile){
  load_atom_set(constrainedfile, "CONSTRAINED", &numConstraints, &constrainedAtomsSet);
}

Bool Molecule::is_atom_in_set(AtomSetList *localAtomsSet, int aid, int *listIdx) const {
  int idx = localAtomsSet->size();
  int rIdx = 0;
  int lIdx = localAtomsSet->size()-1;
  
  while(rIdx <= lIdx){
    int mIdx = (rIdx+lIdx)/2;
    const AtomSet one = localAtomsSet->item(mIdx);

    if(aid < one.aid1){
      //aid could be in [rIdx, mIdx);
      idx = mIdx;
      lIdx = mIdx-1;
    }else if(aid > one.aid1){
      //aid could be inside the atom set "one" or in (mIdx, lIdx];
      if(aid<=one.aid2){
        //found, aid in the atom set "one"
        if(listIdx) *listIdx = mIdx;
        return 1;
      }else{
        rIdx = mIdx+1;
      }
    }else{
      //found, aid is exactly same with one.aid1
      if(listIdx) *listIdx = mIdx;
      return 1;
    }
  }

  //not found
  if(listIdx) *listIdx = idx;
  return 0;
}

#endif


/************************************************************************/
/*                  */
/*      FUNCTION print_atoms        */
/*                  */
/*  print_atoms prints out the list of atoms stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Molecule::print_atoms(Parameters *params)
{
#ifdef MEM_OPT_VERSION
    DebugM(2, "WARNING: this function is not availabe in memory optimized version!\n" << endi);
#else	
  register int i;
  Real sigma;
  Real epsilon;
  Real sigma14;
  Real epsilon14;

  DebugM(2,"ATOM LIST\n" \
      << "******************************************\n" \
                  << "NUM  NAME TYPE RES  MASS    CHARGE CHARGE   FEP-CHARGE"  \
      << "SIGMA   EPSILON SIGMA14 EPSILON14\n" \
        << endi);

  for (i=0; i<numAtoms; i++)
  {
    params->get_vdw_params(&sigma, &epsilon, &sigma14, &epsilon14, 
        atoms[i].vdw_type);

    DebugM(2,i+1 << " " << atomNames[i].atomname  \
              << " " << atomNames[i].atomtype << " " \
              << atomNames[i].resname  << " " << atoms[i].mass  \
        << " " << atoms[i].charge << " " << sigma \
        << " " << epsilon << " " << sigma14 \
        << " " << epsilon14 << "\n" \
        << endi);
  }
#endif  
}
/*      END OF FUNCTION print_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION print_bonds        */
/*                  */
/*  print_bonds prints out the list of bonds stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Molecule::print_bonds(Parameters *params)
{
#ifdef MEM_OPT_VERSION
    DebugM(2, "WARNING: this function is not availabe in memory optimized version!\n" << endi);
#else	
  register int i;
  Real k;
  Real x0;

  DebugM(2,"BOND LIST\n" << "********************************\n" \
      << "ATOM1 ATOM2 TYPE1 TYPE2      k        x0" \
      << endi);

  for (i=0; i<numBonds; i++)
  {
    params->get_bond_params(&k, &x0, bonds[i].bond_type);

    DebugM(2,bonds[i].atom1+1 << " " \
       << bonds[i].atom2+1 << " "   \
       << atomNames[bonds[i].atom1].atomtype << " "  \
       << atomNames[bonds[i].atom2].atomtype << " " << k \
       << " " << x0 << endi);
  }
  
#endif  
}
/*      END OF FUNCTION print_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION print_exclusions      */
/*                  */
/*  print_exlcusions prints out the list of exlcusions stored in    */
/*  this object.  It is inteded mainly for debugging purposes.    */
/*                  */
/************************************************************************/

void Molecule::print_exclusions()
{
#ifdef MEM_OPT_VERSION
    DebugM(2, "WARNING: this function is not availabe in memory optimized version!\n" << endi);
#else
  register int i;

  DebugM(2,"EXPLICIT EXCLUSION LIST\n" \
      << "********************************\n" \
            << "ATOM1 ATOM2 " \
      << endi);

  for (i=0; i<numExclusions; i++)
  {
    DebugM(2,exclusions[i].atom1+1 << "  " \
       << exclusions[i].atom2+1 << endi);
  }
#endif
}
/*      END OF FUNCTION print_exclusions    */

/************************************************************************/
/*                  */
/*      FUNCTION send_Molecule        */
/*                  */
/*  send_Molecule is used by the Master node to distribute the      */
/*   structural information to all the client nodes.  It is NEVER called*/
/*   by the client nodes.              */
/*                  */
/************************************************************************/

void Molecule::send_Molecule(MOStream *msg){
#ifdef MEM_OPT_VERSION
//in the memory optimized version, only the atom signatures are broadcast
//to other Nodes. --Chao Mei

  msg->put(numAtoms);

  msg->put(massPoolSize);
  msg->put(massPoolSize, atomMassPool);

  msg->put(chargePoolSize);
  msg->put(chargePoolSize, atomChargePool); 

  //put atoms' signatures
  msg->put(atomSigPoolSize);
  for(int i=0; i<atomSigPoolSize; i++)
      atomSigPool[i].pack(msg);

  //put atom's exclusion signatures
  msg->put(exclSigPoolSize);
  for(int i=0; i<exclSigPoolSize; i++)
      exclSigPool[i].pack(msg);

  msg->put(numHydrogenGroups);      
  msg->put(maxHydrogenGroupSize);      
  msg->put(numMigrationGroups);      
  msg->put(maxMigrationGroupSize);            
  msg->put(isOccupancyValid);
  msg->put(isBFactorValid);
  
  //put names for atoms
  msg->put(atomNamePoolSize);
  for(int i=0; i<atomNamePoolSize;i++) {
    int len = strlen(atomNamePool[i]);
    msg->put(len);
    msg->put(len*sizeof(char), atomNamePool[i]);
  } 
  
  if(simParams->fixedAtomsOn){
    int numFixedAtomsSet = fixedAtomsSet->size();
    msg->put(numFixedAtoms);
    msg->put(numFixedAtomsSet);
    msg->put(numFixedAtomsSet*sizeof(AtomSet), (char *)(fixedAtomsSet->begin()));
  }

  if (simParams->constraintsOn) {
    int numConstrainedAtomsSet = constrainedAtomsSet->size();
    msg->put(numConstraints);
    msg->put(numConstrainedAtomsSet);
    msg->put(numConstrainedAtomsSet*sizeof(AtomSet), (char *)(constrainedAtomsSet->begin()));
  }
    
#else
  msg->put(numAtoms);
  msg->put(numAtoms*sizeof(Atom), (char*)atoms);
  
  //  Send the bond information
  msg->put(numRealBonds);
  msg->put(numBonds);
 
  if (numBonds)
  {
    msg->put(numBonds*sizeof(Bond), (char*)bonds);
  }

  //  Send the angle information
  msg->put(numAngles);  
  if (numAngles)
  {
    msg->put(numAngles*sizeof(Angle), (char*)angles);
  }  

  //  Send the dihedral information
  msg->put(numDihedrals);
  if (numDihedrals)
  {
    msg->put(numDihedrals*sizeof(Dihedral), (char*)dihedrals);
  }  

  //  Send the improper information
  msg->put(numImpropers);  
  if (numImpropers)
  {
    msg->put(numImpropers*sizeof(Improper), (char*)impropers);
  }

  //  Send the crossterm information
  msg->put(numCrossterms);
  if (numCrossterms)
  {
    msg->put(numCrossterms*sizeof(Crossterm), (char*)crossterms);
  }

  // send the hydrogen bond donor information
  msg->put(numDonors);
  if(numDonors)
  {
    msg->put(numDonors*sizeof(Bond), (char*)donors);
  }

  // send the hydrogen bond acceptor information
  msg->put(numAcceptors);
  if(numAcceptors)
  {
    msg->put(numAcceptors*sizeof(Bond), (char*)acceptors);
  }

  //  Send the exclusion information  
  msg->put(numExclusions);
  if (numExclusions)
  {
    msg->put(numExclusions*sizeof(Exclusion), (char*)exclusions);
  }      
  //  Send the constraint information, if used
  if (simParams->constraintsOn)
  {
     msg->put(numConstraints);
     
     msg->put(numAtoms, consIndexes);
     
     if (numConstraints)
     {
       msg->put(numConstraints*sizeof(ConstraintParams), (char*)consParams);
     }
  }
#endif
  
  /* BEGIN gf */
  // Send the gridforce information, if used
  if (simParams->mgridforceOn)
  {
    DebugM(3, "Sending gridforce info\n" << endi);
    msg->put(numGridforceGrids);
    
    for (int gridnum = 0; gridnum < numGridforceGrids; gridnum++) {
      msg->put(numGridforces[gridnum]);
      msg->put(numAtoms, gridfrcIndexes[gridnum]);
      if (numGridforces[gridnum])
      {
       msg->put(numGridforces[gridnum]*sizeof(GridforceParams), (char*)gridfrcParams[gridnum]);
      }
      GridforceGrid::pack_grid(gridfrcGrid[gridnum], msg);
    }
  }
  /* END gf */
  
  //  Send the stirring information, if used
  if (simParams->stirOn)
  {
     //CkPrintf ("DEBUG: putting numStirredAtoms..\n");
     msg->put(numStirredAtoms);
     //CkPrintf ("DEBUG: putting numAtoms,stirIndexes.. numAtoms=%d\n",numStirredAtoms);
     msg->put(numAtoms, stirIndexes);
     //CkPrintf ("DEBUG: if numStirredAtoms..\n");
     if (numStirredAtoms)
     {
       //CkPrintf ("DEBUG: big put, with (char*)stirParams\n");
       msg->put(numStirredAtoms*sizeof(StirParams), (char*)stirParams);
     }
  }
  
  
  //  Send the moving drag information, if used
  if (simParams->movDragOn) {
     msg->put(numMovDrag);
     msg->put(numAtoms, movDragIndexes);
     if (numMovDrag)
     {
       msg->put(numMovDrag*sizeof(MovDragParams), (char*)movDragParams);
     }
  }
  
  //  Send the rotating drag information, if used
  if (simParams->rotDragOn) {
     msg->put(numRotDrag);
     msg->put(numAtoms, rotDragIndexes);
     if (numRotDrag)
     {
       msg->put(numRotDrag*sizeof(RotDragParams), (char*)rotDragParams);
     }
  }
  
  //  Send the "constant" torque information, if used
  if (simParams->consTorqueOn) {
     msg->put(numConsTorque);
     msg->put(numAtoms, consTorqueIndexes);
     if (numConsTorque)
     {
       msg->put(numConsTorque*sizeof(ConsTorqueParams), (char*)consTorqueParams);
     }
  }
  
  // Send the constant force information, if used
  if (simParams->consForceOn)
  { msg->put(numConsForce);
    msg->put(numAtoms, consForceIndexes);
    if (numConsForce)
      msg->put(numConsForce*sizeof(Vector), (char*)consForce);
  }
  
  if (simParams->excludeFromPressure) {
    msg->put(numExPressureAtoms);
    msg->put(numAtoms, exPressureAtomFlags);
  }
  
#ifndef MEM_OPT_VERSION
  //  Send the langevin parameters, if active
  if (simParams->langevinOn || simParams->tCoupleOn)
  {
    msg->put(numAtoms, langevinParams);
  }
  
  //  Send fixed atoms, if active
  if (simParams->fixedAtomsOn)
  {
    msg->put(numFixedAtoms);
    msg->put(numAtoms, fixedAtomFlags);
  msg->put(numFixedRigidBonds);
  }
  
  if (simParams->qmForcesOn)
  {
    msg->put(numAtoms, qmAtomGroup);
    msg->put(qmNumQMAtoms);
    msg->put(qmNumQMAtoms, qmAtmChrg);
    msg->put(qmNumQMAtoms, qmAtmIndx);
    msg->put(qmNoPC);
    msg->put(qmNumBonds);
    msg->put(qmMeNumBonds);
    msg->put(qmMeNumBonds, qmMeMMindx);
    msg->put(qmMeNumBonds, qmMeQMGrp);
    msg->put(qmPCFreq);
    msg->put(qmNumGrps);
    msg->put(qmNumGrps, qmGrpID);
    msg->put(qmNumGrps, qmCustPCSizes);
    msg->put(qmTotCustPCs);
    msg->put(qmTotCustPCs, qmCustomPCIdxs);
  }
  
  //fepb
  // send fep atom info
  if (simParams->alchOn || simParams->lesOn || simParams->pairInteractionOn) {
    msg->put(numFepInitial);
    msg->put(numFepFinal);
    msg->put(numAtoms*sizeof(char), (char*)fepAtomFlags);
  }
  //fepe

  #ifdef OPENATOM_VERSION
  // needs to be refactored into its own openatom version
  if (simParams->openatomOn ) {
    msg->put(numFepInitial);
    msg->put(numAtoms*sizeof(char), (char*)fepAtomFlags);
  }
  #endif //OPENATOM_VERSION
  
  // DRUDE: send data read from PSF
  msg->put(is_lonepairs_psf);
  if (is_lonepairs_psf) {
    msg->put(numLphosts);
    msg->put(numLphosts*sizeof(Lphost), (char*)lphosts);
  }
  msg->put(is_drude_psf);
  if (is_drude_psf) {
    msg->put(numAtoms*sizeof(DrudeConst), (char*)drudeConsts);
    msg->put(numAnisos);
    msg->put(numAnisos*sizeof(Aniso), (char*)anisos);
  }
  // DRUDE

  //LCPO
  if (simParams->LCPOOn) {
    msg->put(numAtoms, (int*)lcpoParamType);
  }
  
  //Send GromacsPairStuff -- JLai
  if (simParams->goGroPair) {
    msg->put(numLJPair);
    msg->put(numLJPair,indxLJA);
    msg->put(numLJPair,indxLJB);
    msg->put(numLJPair,pairC6);
    msg->put(numLJPair,pairC12);
    msg->put(numLJPair,gromacsPair_type);
    msg->put((numAtoms),pointerToLJBeg);
    msg->put((numAtoms),pointerToLJEnd);
    msg->put(numGaussPair);
    msg->put(numGaussPair,indxGaussA);
    msg->put(numGaussPair,indxGaussB);
    msg->put(numGaussPair,gA);
    msg->put(numGaussPair,gMu1);
    msg->put(numGaussPair,giSigma1);
    msg->put(numGaussPair,gMu2);
    msg->put(numGaussPair,giSigma2);
    msg->put(numGaussPair,gRepulsive);
    msg->put((numAtoms),pointerToGaussBeg);
    msg->put((numAtoms),pointerToGaussEnd);
  }
#endif

  // Broadcast the message to the other nodes
  msg->end();
  delete msg;

#ifdef MEM_OPT_VERSION

  build_excl_check_signatures();

  //set num{Calc}Tuples(Bonds,...,Impropers) to 0
  numBonds = numCalcBonds = 0;
  numAngles = numCalcAngles = 0;
  numDihedrals = numCalcDihedrals = 0;
  numImpropers = numCalcImpropers = 0;
  numCrossterms = numCalcCrossterms = 0;
  numTotalExclusions = numCalcExclusions = numCalcFullExclusions = 0;  
  // JLai
  numLJPair = numCalcLJPair = 0;
  // End of JLai

#else

  //  Now build arrays of indexes into these arrays by atom      
  build_lists_by_atom();

#endif
}
 /*      END OF FUNCTION send_Molecule      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION receive_Molecule      */
    /*                  */
    /*  receive_Molecule is used by all the clients to receive the  */
    /*   structural data sent out by the master node.  It is NEVER called   */
    /*   by the Master node.            */
    /*                  */
    /************************************************************************/

void Molecule::receive_Molecule(MIStream *msg){
  //  Get the atom information
  msg->get(numAtoms);

#ifdef MEM_OPT_VERSION
//in the memory optimized version, only the atom signatures are recved
//from the master Node. --Chao Mei

  msg->get(massPoolSize);
  if(atomMassPool) delete [] atomMassPool;
  atomMassPool = new Real[massPoolSize];
  msg->get(massPoolSize, atomMassPool);

  msg->get(chargePoolSize);
  if(atomChargePool) delete [] atomChargePool;
  atomChargePool = new Real[chargePoolSize];
  msg->get(chargePoolSize, atomChargePool);

  //get atoms' signatures
  msg->get(atomSigPoolSize);
  if(atomSigPool) delete [] atomSigPool;
  atomSigPool = new AtomSignature[atomSigPoolSize];
  for(int i=0; i<atomSigPoolSize; i++)
      atomSigPool[i].unpack(msg);

  //get exclusions' signatures
  msg->get(exclSigPoolSize);
  if(exclSigPool) delete [] exclSigPool;
  exclSigPool = new ExclusionSignature[exclSigPoolSize];
  for(int i=0; i<exclSigPoolSize; i++)
      exclSigPool[i].unpack(msg);
 
  msg->get(numHydrogenGroups);      
  msg->get(maxHydrogenGroupSize);      
  msg->get(numMigrationGroups);      
  msg->get(maxMigrationGroupSize);      
  msg->get(isOccupancyValid);
  msg->get(isBFactorValid);

   //get names for atoms
  msg->get(atomNamePoolSize);
  atomNamePool = new char *[atomNamePoolSize];
  for(int i=0; i<atomNamePoolSize;i++) {
    int len;
    msg->get(len);
    atomNamePool[i] = nameArena->getNewArray(len+1);
    msg->get(len, atomNamePool[i]);
  }
  
  if(simParams->fixedAtomsOn){
    int numFixedAtomsSet;
    msg->get(numFixedAtoms);
    msg->get(numFixedAtomsSet);
    fixedAtomsSet = new AtomSetList(numFixedAtomsSet);
    msg->get(numFixedAtomsSet*sizeof(AtomSet), (char *)(fixedAtomsSet->begin()));
  } 

  if(simParams->constraintsOn){
    int numConstrainedAtomsSet;
    msg->get(numConstraints);
    msg->get(numConstrainedAtomsSet);
    constrainedAtomsSet = new AtomSetList(numConstrainedAtomsSet);
    msg->get(numConstrainedAtomsSet*sizeof(AtomSet), (char *)(constrainedAtomsSet->begin()));
  } 

#else
  delete [] atoms;
  atoms= new Atom[numAtoms];  
  msg->get(numAtoms*sizeof(Atom), (char*)atoms);

  //  Get the bond information
  msg->get(numRealBonds);
  msg->get(numBonds);    
  if (numBonds)
  {
    delete [] bonds;
    bonds=new Bond[numBonds]; 
    msg->get(numBonds*sizeof(Bond), (char*)bonds);
  }  
  
  //  Get the angle information
  msg->get(numAngles);  
  if (numAngles)
  {
    delete [] angles;
    angles=new Angle[numAngles];  
    msg->get(numAngles*sizeof(Angle), (char*)angles);
  }  
  
  //  Get the dihedral information
  msg->get(numDihedrals);    
  if (numDihedrals)
  {
    delete [] dihedrals;
    dihedrals=new Dihedral[numDihedrals];  
    msg->get(numDihedrals*sizeof(Dihedral), (char*)dihedrals);
  }  
  
  //  Get the improper information
  msg->get(numImpropers);
  if (numImpropers)
  {
    delete [] impropers;
    impropers=new Improper[numImpropers];  
    msg->get(numImpropers*sizeof(Improper), (char*)impropers);
  }
  
  //  Get the crossterm information
  msg->get(numCrossterms);
  if (numCrossterms)
  {
    delete [] crossterms;
    crossterms=new Crossterm[numCrossterms];  
    msg->get(numCrossterms*sizeof(Crossterm), (char*)crossterms);
  }
  
  //  Get the hydrogen bond donors
  msg->get(numDonors);  
  if (numDonors)
  {
    delete [] donors;
    donors=new Bond[numDonors];  
    msg->get(numDonors*sizeof(Bond), (char*)donors);
  }
  
  //  Get the hydrogen bond acceptors
  msg->get(numAcceptors);  
  if (numAcceptors)
  {
    delete [] acceptors;
    acceptors=new Bond[numAcceptors];  
    msg->get(numAcceptors*sizeof(Bond), (char*)acceptors);
  }
  
  //  Get the exclusion information 
  msg->get(numExclusions);  
  if (numExclusions)
  {
    delete [] exclusions;
    exclusions=new Exclusion[numExclusions];  
    msg->get(numExclusions*sizeof(Exclusion), (char*)exclusions);
  }
        
      //  Get the constraint information, if they are active
      if (simParams->constraintsOn)
      {
         msg->get(numConstraints);

         delete [] consIndexes;
         consIndexes = new int32[numAtoms];
         
         msg->get(numAtoms, consIndexes);
         
         if (numConstraints)
         {
           delete [] consParams;
           consParams = new ConstraintParams[numConstraints];
      
           msg->get(numConstraints*sizeof(ConstraintParams), (char*)consParams);
         }
      }
#endif

      /* BEGIN gf */
      if (simParams->mgridforceOn)
      {
	 DebugM(3, "Receiving gridforce info\n");
	 
	 msg->get(numGridforceGrids);
	 
	 DebugM(3, "numGridforceGrids = " << numGridforceGrids << "\n");
	 
	 delete [] numGridforces;
	 numGridforces = new int[numGridforceGrids];
	 
	 delete [] gridfrcIndexes;	// Should I be deleting elements of these first?
	 delete [] gridfrcParams;
	 delete [] gridfrcGrid;
	 gridfrcIndexes = new int32*[numGridforceGrids];
	 gridfrcParams = new GridforceParams*[numGridforceGrids];
	 gridfrcGrid = new GridforceGrid*[numGridforceGrids];
	 
	 int grandTotalGrids = 0;
	 for (int gridnum = 0; gridnum < numGridforceGrids; gridnum++) {
	     msg->get(numGridforces[gridnum]);
	     
	     gridfrcIndexes[gridnum] = new int32[numAtoms];
	     msg->get(numAtoms, gridfrcIndexes[gridnum]);
	 
	     if (numGridforces[gridnum])
	     {
		 gridfrcParams[gridnum] = new GridforceParams[numGridforces[gridnum]];
		 msg->get(numGridforces[gridnum]*sizeof(GridforceParams), (char*)gridfrcParams[gridnum]);
	     }
	     
	     gridfrcGrid[gridnum] = GridforceGrid::unpack_grid(gridnum, msg);
	     
	     grandTotalGrids += gridfrcGrid[gridnum]->get_total_grids();
	 }
      }
      /* END gf */
      
      //  Get the stirring information, if stirring is  active
      if (simParams->stirOn)
      {
         msg->get(numStirredAtoms);

         delete [] stirIndexes;
         stirIndexes = new int32[numAtoms];
         
         msg->get(numAtoms, stirIndexes);
         
         if (numStirredAtoms)
         {
           delete [] stirParams;
           stirParams = new StirParams[numStirredAtoms];
      
           msg->get(numStirredAtoms*sizeof(StirParams), (char*)stirParams);
         }
      }
      
      //  Get the moving drag information, if it is active
      if (simParams->movDragOn) {
         msg->get(numMovDrag);
         delete [] movDragIndexes;
         movDragIndexes = new int32[numAtoms];
         msg->get(numAtoms, movDragIndexes);
         if (numMovDrag)
         {
           delete [] movDragParams;
           movDragParams = new MovDragParams[numMovDrag];
           msg->get(numMovDrag*sizeof(MovDragParams), (char*)movDragParams);
         }
      }
      
      //  Get the rotating drag information, if it is active
      if (simParams->rotDragOn) {
         msg->get(numRotDrag);
         delete [] rotDragIndexes;
         rotDragIndexes = new int32[numAtoms];
         msg->get(numAtoms, rotDragIndexes);
         if (numRotDrag)
         {
           delete [] rotDragParams;
           rotDragParams = new RotDragParams[numRotDrag];
           msg->get(numRotDrag*sizeof(RotDragParams), (char*)rotDragParams);
         }
      }
      
      //  Get the "constant" torque information, if it is active
      if (simParams->consTorqueOn) {
         msg->get(numConsTorque);
         delete [] consTorqueIndexes;
         consTorqueIndexes = new int32[numAtoms];
         msg->get(numAtoms, consTorqueIndexes);
         if (numConsTorque)
         {
           delete [] consTorqueParams;
           consTorqueParams = new ConsTorqueParams[numConsTorque];
           msg->get(numConsTorque*sizeof(ConsTorqueParams), (char*)consTorqueParams);
         }
      }
      
      // Get the constant force information, if it's active
      if (simParams->consForceOn)
      { msg->get(numConsForce);
        delete [] consForceIndexes;
        consForceIndexes = new int32[numAtoms];
        msg->get(numAtoms, consForceIndexes);
        if (numConsForce)
        { delete [] consForce;
          consForce = new Vector[numConsForce];
          msg->get(numConsForce*sizeof(Vector), (char*)consForce);
        }
      }

      if (simParams->excludeFromPressure) {
        exPressureAtomFlags = new int32[numAtoms];
        msg->get(numExPressureAtoms);
        msg->get(numAtoms, exPressureAtomFlags);
      }

#ifndef MEM_OPT_VERSION
      //  Get the langevin parameters, if they are active
      if (simParams->langevinOn || simParams->tCoupleOn)
      {
        delete [] langevinParams;
        langevinParams = new Real[numAtoms];

        msg->get(numAtoms, langevinParams);
      }

      //  Get the fixed atoms, if they are active
      if (simParams->fixedAtomsOn)
      {
        delete [] fixedAtomFlags;
        fixedAtomFlags = new int32[numAtoms];

        msg->get(numFixedAtoms);
        msg->get(numAtoms, fixedAtomFlags);
        msg->get(numFixedRigidBonds);
      }

      if (simParams->qmForcesOn)
      {
        if( qmAtomGroup != 0)
            delete [] qmAtomGroup;
        qmAtomGroup = new Real[numAtoms];
        
        msg->get(numAtoms, qmAtomGroup);
        
        msg->get(qmNumQMAtoms);
        
        if( qmAtmChrg != 0)
            delete [] qmAtmChrg;
        qmAtmChrg = new Real[qmNumQMAtoms];
        
        msg->get(qmNumQMAtoms, qmAtmChrg);
        
        if( qmAtmIndx != 0)
            delete [] qmAtmIndx;
        qmAtmIndx = new int[qmNumQMAtoms];
        
        msg->get(qmNumQMAtoms, qmAtmIndx);
        
        msg->get(qmNoPC);
        
        msg->get(qmNumBonds);
        
        msg->get(qmMeNumBonds);
        
        if( qmMeMMindx != 0)
            delete [] qmMeMMindx;
        qmMeMMindx = new int[qmMeNumBonds];
        
        msg->get(qmMeNumBonds, qmMeMMindx);
        
        if( qmMeQMGrp != 0)
            delete [] qmMeQMGrp;
        qmMeQMGrp = new Real[qmMeNumBonds];
        
        msg->get(qmMeNumBonds, qmMeQMGrp);
        
        msg->get(qmPCFreq);
        
        msg->get(qmNumGrps);
        
        if( qmGrpID != 0)
            delete [] qmGrpID;
        qmGrpID = new Real[qmNumGrps];
        msg->get(qmNumGrps, qmGrpID);
        
        if( qmCustPCSizes != 0)
            delete [] qmCustPCSizes;
        qmCustPCSizes = new int[qmNumGrps];
        msg->get(qmNumGrps, qmCustPCSizes);
        
        msg->get(qmTotCustPCs);
        
        if( qmCustomPCIdxs != 0)
            delete [] qmCustomPCIdxs;
        qmCustomPCIdxs = new int[qmTotCustPCs];
        msg->get(qmTotCustPCs, qmCustomPCIdxs);
      }
    
//fepb
      //receive fep atom info
      if (simParams->alchOn || simParams->lesOn || simParams->pairInteractionOn) {
        delete [] fepAtomFlags;
        fepAtomFlags = new unsigned char[numAtoms];

        msg->get(numFepInitial);
        msg->get(numFepFinal);
        msg->get(numAtoms*sizeof(unsigned char), (char*)fepAtomFlags);
      }
//fepe

#ifdef OPENATOM_VERSION
      // This needs to be refactored into its own version
      if (simParams->openatomOn) {
        delete [] fepAtomFlags;
        fepAtomFlags = new unsigned char[numAtoms];

        msg->get(numFepInitial);
        msg->get(numAtoms*sizeof(unsigned char), (char*)fepAtomFlags);
#endif //OPENATOM_VERSION

      // DRUDE: receive data read from PSF
      msg->get(is_lonepairs_psf);
      if (is_lonepairs_psf) {
        msg->get(numLphosts);
        delete[] lphosts;
        lphosts = new Lphost[numLphosts];
        msg->get(numLphosts*sizeof(Lphost), (char*)lphosts);
      }
      msg->get(is_drude_psf);
      if (is_drude_psf) {
        delete[] drudeConsts;
        drudeConsts = new DrudeConst[numAtoms];
        msg->get(numAtoms*sizeof(DrudeConst), (char*)drudeConsts);
        msg->get(numAnisos);
        delete[] anisos;
        anisos = new Aniso[numAnisos];
        msg->get(numAnisos*sizeof(Aniso), (char*)anisos);
      }
      // DRUDE

  //LCPO
  if (simParams->LCPOOn) {
    delete [] lcpoParamType;
    lcpoParamType = new int[numAtoms];
    msg->get(numAtoms, (int*)lcpoParamType);
  }

  //Receive GromacsPairStuff -- JLai

  if (simParams->goGroPair) {
    msg->get(numLJPair);
    delete [] indxLJA;
    indxLJA = new int[numLJPair];
    msg->get(numLJPair,indxLJA);
    delete [] indxLJB;
    indxLJB = new int[numLJPair];
    msg->get(numLJPair,indxLJB);
    delete [] pairC6;
    pairC6 = new Real[numLJPair];    
    msg->get(numLJPair,pairC6);
    delete [] pairC12;
    pairC12 = new Real[numLJPair];
    msg->get(numLJPair,pairC12);
    delete [] gromacsPair_type;
    gromacsPair_type = new int[numLJPair];
    msg->get(numLJPair,gromacsPair_type);
    delete [] pointerToLJBeg;
    pointerToLJBeg = new int[numAtoms];
    msg->get((numAtoms),pointerToLJBeg);
    delete [] pointerToLJEnd;
    pointerToLJEnd = new int[numAtoms];
    msg->get((numAtoms),pointerToLJEnd);
    // JLai
    delete [] gromacsPair;
    gromacsPair = new GromacsPair[numLJPair];
    for(int i=0; i < numLJPair; i++) {
	gromacsPair[i].atom1 = indxLJA[i];
	gromacsPair[i].atom2 = indxLJB[i];
	gromacsPair[i].pairC6  = pairC6[i];
	gromacsPair[i].pairC12 = pairC12[i];
	gromacsPair[i].gromacsPair_type = gromacsPair_type[i];
    }
    //
    msg->get(numGaussPair);
    delete [] indxGaussA;
    indxGaussA = new int[numGaussPair];
    msg->get(numGaussPair,indxGaussA);
    delete [] indxGaussB;
    indxGaussB = new int[numGaussPair];
    msg->get(numGaussPair,indxGaussB);
    delete [] gA;
    gA = new Real[numGaussPair];
    msg->get(numGaussPair,gA);
    delete [] gMu1;
    gMu1 = new Real[numGaussPair];
    msg->get(numGaussPair,gMu1);
    delete [] giSigma1;
    giSigma1 = new Real[numGaussPair];
    msg->get(numGaussPair,giSigma1);
    delete [] gMu2;
    gMu2 = new Real[numGaussPair];
    msg->get(numGaussPair,gMu2);
    delete [] giSigma2;
    giSigma2 = new Real[numGaussPair];
    msg->get(numGaussPair,giSigma2);
    delete [] gRepulsive;
    gRepulsive = new Real[numGaussPair];
    msg->get(numGaussPair,gRepulsive);
    delete [] pointerToGaussBeg;
    pointerToGaussBeg = new int[numAtoms];
    msg->get((numAtoms),pointerToGaussBeg);
    delete [] pointerToGaussEnd;
    pointerToGaussEnd = new int[numAtoms];
    msg->get((numAtoms),pointerToGaussEnd);
    //
  }
#endif

      //  Now free the message 
      delete msg;

#ifdef MEM_OPT_VERSION

      build_excl_check_signatures();

    //set num{Calc}Tuples(Bonds,...,Impropers) to 0
    numBonds = numCalcBonds = 0;
    numAngles = numCalcAngles = 0;
    numDihedrals = numCalcDihedrals = 0;
    numImpropers = numCalcImpropers = 0;
    numCrossterms = numCalcCrossterms = 0;
    numTotalExclusions = numCalcExclusions = numCalcFullExclusions = 0;  
    // JLai
    numLJPair = numCalcLJPair = 0;
    // End of JLai

#else

      //  analyze the data and find the status of each atom
      build_atom_status();
      build_lists_by_atom();      

      
#endif
}
 /*      END OF FUNCTION receive_Molecule    */

/* BEGIN gf */
    /************************************************************************/
    /*                                                                      */
    /*      FUNCTION build_gridforce_params                                 */
    /*                                                                      */
    /*   INPUTS:                                                            */
    /*  gridfrcfile - Value of gridforcefile from config file               */
    /*  gridfrccol - Value of gridforcecol from config file                 */
    /*  gridfrcchrgcol - Value of gridforcechargecol from config file	    */
    /*  potfile - Value of gridforcepotfile from config file                  */
    /*  initial_pdb - PDB object that contains initial positions            */
    /*  cwd - Current working directory                                     */
    /*                                                                      */
    // This function builds all the parameters that are necessary to
    // do gridforcing. This involves looking through a PDB object to
    // determine which atoms are to be gridforced, and what the force
    // multiplier is for each atom.  This information is then stored
    // in the arrays gridfrcIndexes and gridfrcParams.
    /************************************************************************/

void Molecule::build_gridforce_params(StringList *gridfrcfile,
				      StringList *gridfrccol,
				      StringList *gridfrcchrgcol,
				      StringList *potfile,
				      PDB *initial_pdb,
				      char *cwd)
{
    PDB *kPDB;
    register int i;		//  Loop counters
    register int j;
    register int k;

    DebugM(3,  "Entered build_gridforce_params multi...\n");
//     DebugM(3, "\tgridfrcfile = " << gridfrcfile->data << endi);
//     DebugM(3, "\tgridfrccol = " << gridfrccol->data << endi);
    
    MGridforceParams* mgridParams = simParams->mgridforcelist.get_first();
    numGridforceGrids = 0;
    while (mgridParams != NULL) {
	numGridforceGrids++;
	mgridParams = mgridParams->next;
    }
    
    DebugM(3, "numGridforceGrids = " << numGridforceGrids << "\n");
    gridfrcIndexes = new int32*[numGridforceGrids];
    gridfrcParams = new GridforceParams*[numGridforceGrids];
    gridfrcGrid = new GridforceGrid*[numGridforceGrids];
    numGridforces = new int[numGridforceGrids];
    
    int grandTotalGrids = 0;	// including all subgrids
    
    mgridParams = simParams->mgridforcelist.get_first();
    for (int gridnum = 0; gridnum < numGridforceGrids; gridnum++) {
	int current_index=0;	//  Index into values used
	int kcol = 5;		//  Column to look for force constant in
	int qcol = 0;		//  Column for charge (default 0: use electric charge)
	Real kval = 0;		//  Force constant value retreived
	char filename[129];	//  PDB filename
	char potfilename[129];	//  Potential file name
	
	if (mgridParams == NULL) {
	    NAMD_die("Problem with mgridParams!");
	}
	
	// Now load values from mgridforcelist object
	if (mgridParams->gridforceFile == NULL)
	{
	    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, gridforceFile required.");
	    kPDB = initial_pdb;
	}
	else
	{
	    DebugM(4, "mgridParams->gridforceFile = " << mgridParams->gridforceFile << "\n" << endi);
	    
	    if ( (cwd == NULL) || (mgridParams->gridforceFile[0] == '/') )
	    {
		strcpy(filename, mgridParams->gridforceFile);
	    }
	    else
	    {
		strcpy(filename, cwd);
		strcat(filename, mgridParams->gridforceFile);
	    }
	
	    kPDB = new PDB(filename);
	    if ( kPDB == NULL )
	    {
		NAMD_die("Memory allocation failed in Molecule::build_gridforce_params");
	    }
	   
	    if (kPDB->num_atoms() != numAtoms)
	    {
		NAMD_die("Number of atoms in grid force PDB doesn't match coordinate PDB");
	    }
	}

	//  Get the column that the force constant is going to be in.  It
	//  can be in any of the 5 floating point fields in the PDB, according
	//  to what the user wants.  The allowable fields are X, Y, Z, O, or
	//  B which correspond to the 1st, 2nd, ... 5th floating point fields.
	//  The default is the 5th field, which is beta (temperature factor)
	if (mgridParams->gridforceCol == NULL)
	{
	    kcol = 5;
	}
	else
	{
	    if (strcasecmp(mgridParams->gridforceCol, "X") == 0)
	    {
		kcol=1;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "Y") == 0)
	    {
		kcol=2;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "Z") == 0)
	    {
		kcol=3;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "O") == 0)
	    {
		kcol=4;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "B") == 0)
	    {
		kcol=5;
	    }
	    else
	    {
		NAMD_die("gridforcecol must have value of X, Y, Z, O, or B");
	    }
	}
    
	//  Get the column that the charge is going to be in.
        if (mgridParams->gridforceQcol == NULL)
	{
	    qcol = 0;	// Default: don't read charge from file, use electric charge
	}
	else
	{
	    if (strcasecmp(mgridParams->gridforceQcol, "X") == 0)
	    {
		qcol=1;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "Y") == 0)
	    {
		qcol=2;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "Z") == 0)
	    {
		qcol=3;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "O") == 0)
	    {
		qcol=4;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "B") == 0)
	    {
		qcol=5;
	    }
	    else
	    {
		NAMD_die("gridforcechargecol must have value of X, Y, Z, O, or B");
	    }
	}
    
	if (kcol == qcol) {
	    NAMD_die("gridforcecol and gridforcechargecol cannot have same value");
	}

    
	//  Allocate an array that will store an index into the constraint
	//  parameters for each atom.  If the atom is not constrained, its
	//  value will be set to -1 in this array.
	gridfrcIndexes[gridnum] = new int32[numAtoms];
       
	if (gridfrcIndexes[gridnum] == NULL)
	{
	    NAMD_die("memory allocation failed in Molecule::build_gridforce_params()");
	}
	
	//  Loop through all the atoms and find out which ones are constrained
	for (i=0; i<numAtoms; i++)
	{
	    //  Get the k value based on where we were told to find it
	    switch (kcol)
	    {
	    case 1:
		kval = (kPDB->atom(i))->xcoor();
		break;
	    case 2:
		kval = (kPDB->atom(i))->ycoor();
		break;
	    case 3:
		kval = (kPDB->atom(i))->zcoor();
		break;
	    case 4:
		kval = (kPDB->atom(i))->occupancy();
		break;
	    case 5:
		kval = (kPDB->atom(i))->temperaturefactor();
		break;
	    }
	   
	    if (kval > 0.0)
	    {
		//  This atom is constrained
		gridfrcIndexes[gridnum][i] = current_index;
		current_index++;
	    }
	    else
	    {
		//  This atom is not constrained
		gridfrcIndexes[gridnum][i] = -1;
	    }
	}
    
	if (current_index == 0)
	{
	    //  Constraints were turned on, but there weren't really any constrained
	    iout << iWARN << "NO GRIDFORCE ATOMS WERE FOUND, BUT GRIDFORCE IS ON . . .\n" << endi;
	}
	else
	{
	    //  Allocate an array to hold the constraint parameters
	    gridfrcParams[gridnum] = new GridforceParams[current_index];
	    if (gridfrcParams[gridnum] == NULL)
	    {
		NAMD_die("memory allocation failed in Molecule::build_gridforce_params");
	    }
	}
    
	numGridforces[gridnum] = current_index;

	//  Loop through all the atoms and assign the parameters for those
	//  that are constrained
	for (i=0; i<numAtoms; i++)
	{
	    if (gridfrcIndexes[gridnum][i] != -1)
	    {
		//  This atom has grid force, so get the k value again
		switch (kcol)
		{
		case 1:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->xcoor();
		    break;
		case 2:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->ycoor();
		    break;
		case 3:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->zcoor();
		    break;
		case 4:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->occupancy();
		    break;
		case 5:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->temperaturefactor();
		    break;
		}
	    
		//  Also get charge column
		switch (qcol)
		{
		case 0:
#ifdef MEM_OPT_VERSION
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = atomChargePool[eachAtomCharge[i]];
#else
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = atoms[i].charge;
#endif
		    break;
		case 1:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->xcoor();
		    break;
		case 2:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->ycoor();
		    break;
		case 3:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->zcoor();
		    break;
		case 4:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->occupancy();
		    break;
		case 5:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->temperaturefactor();
		    break;
		}
	    }
	}
       
	//  If we had to create new PDB objects, delete them now
	if (mgridParams->gridforceFile != NULL)
	{
	    delete kPDB;
	}
    
	//  Now we fill in our grid information
    
	// Open potential file
	if ( (cwd == NULL) || (mgridParams->gridforceVfile[0] == '/') )
	{
	    strcpy(potfilename, mgridParams->gridforceVfile);
	}
	else
	{
	    strcpy(potfilename, cwd);
	    strcat(potfilename, mgridParams->gridforceVfile);
	}
    
//        iout << iINFO << "Allocating grid " << gridnum
//             << "\n" << endi;
	
	DebugM(3, "allocating GridforceGrid(" << gridnum << ")\n");
	gridfrcGrid[gridnum] = GridforceGrid::new_grid(gridnum, potfilename, simParams, mgridParams);
	
	grandTotalGrids += gridfrcGrid[gridnum]->get_total_grids();
	DebugM(4, "grandTotalGrids = " << grandTotalGrids << "\n" << endi);
	
	// Finally, get next mgridParams pointer
	mgridParams = mgridParams->next;
    }
}
/* END gf */


#endif  // MOLECULE2_C undefined = first object file
#ifdef MOLECULE2_C  // second object file


    /************************************************************************/
    /*                  */
    /*      FUNCTION build_constraint_params    */
    /*                  */
    /*   INPUTS:                */
    /*  consref - Value of consref parameter from config file    */
    /*  conskfile - Value of conskfile from config file      */
    /*  conskcol - Value of conskcol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*  cwd - Current working directory          */
    /*                  */
    /*  This function builds all the parameters that are necessary  */
    /*   to do harmonic constraints.  This involves looking through    */
    /*   one or more PDB objects to determine which atoms are constrained,  */
    /*   and what the force constant and reference position is force each   */
    /*   atom that is constrained.  This information is then stored    */
    /*   in the arrays consIndexes and consParams.        */
    /*                  */
    /************************************************************************/

    void Molecule::build_constraint_params(StringList *consref, 
             StringList *conskfile, 
             StringList *conskcol, 
             PDB *initial_pdb,
             char *cwd)
       
    {
       PDB *refPDB, *kPDB;    //  Pointer to other PDB's if used
       register int i;      //  Loop counter
       int current_index=0;    //  Index into values used
       int kcol = 4;      //  Column to look for force constant in
       Real kval = 0;      //  Force constant value retreived
       char filename[129];    //  PDB filename
       
       //  Get the PDB object that contains the reference positions.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  i.e., constrain
       //  the atoms around their initial position.  This is the most likely
       //  case anyway
       if (consref == NULL)
       {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, consref required.");
    refPDB = initial_pdb;
       }
       else
       {
    if (consref->next != NULL)
    {
       NAMD_die("Multiple definitions of constraint reference file in configruation file");
    }

    if ( (cwd == NULL) || (consref->data[0] == '/') )
    {
         strcpy(filename, consref->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, consref->data);
    }
    
    refPDB = new PDB(filename);
    if ( refPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_constraint_params");
    }
    
    if (refPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in constraint reference PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the PDB to read the force constants from.  Again, if the user
       //  gave us another file name, open that one.  Otherwise, just use
       //  the PDB with the initial coordinates
       if (conskfile == NULL)
       {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, conskfile required.");
    kPDB = initial_pdb;
       }
       else
       {
    if (conskfile->next != NULL)
    {
       NAMD_die("Multiple definitions of constraint constant file in configuration file");
    }

    if ( (consref != NULL) && (strcasecmp(consref->data, conskfile->data) == 0) )
    {
       //  Same PDB used for reference positions and force constants
       kPDB = refPDB; 
    }
    else
    {
      if ( (cwd == NULL) || (conskfile->data[0] == '/') )
      {
        strcpy(filename, conskfile->data);
      }
      else
      {
        strcpy(filename, cwd);
        strcat(filename, conskfile->data);
      }

      kPDB = new PDB(filename);
      if ( kPDB == NULL )
      {
        NAMD_die("Memory allocation failed in Molecule::build_constraint_params");
      }
    
      if (kPDB->num_atoms() != numAtoms)
      {
         NAMD_die("Number of atoms in constraint constant PDB doesn't match coordinate PDB");
      }
    }
       }
       
       //  Get the column that the force constant is going to be in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (conskcol == NULL)
       {
    kcol = 4;
       }
       else
       {
    if (conskcol->next != NULL)
    {
       NAMD_die("Multiple definitions of harmonic constraint column in config file");
    }
    
    if (strcasecmp(conskcol->data, "X") == 0)
    {
       kcol=1;
    }
    else if (strcasecmp(conskcol->data, "Y") == 0)
    {
       kcol=2;
    }
    else if (strcasecmp(conskcol->data, "Z") == 0)
    {
       kcol=3;
    }
    else if (strcasecmp(conskcol->data, "O") == 0)
    {
       kcol=4;
    }
    else if (strcasecmp(conskcol->data, "B") == 0)
    {
       kcol=5;
    }
    else
    {
       NAMD_die("conskcol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate an array that will store an index into the constraint
       //  parameters for each atom.  If the atom is not constrained, its
       //  value will be set to -1 in this array.
       consIndexes = new int32[numAtoms];
       
       if (consIndexes == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_constraint_params()");
       }
       
       //  Loop through all the atoms and find out which ones are constrained
       for (i=0; i<numAtoms; i++)
       {
    //  Get the k value based on where we were told to find it
    switch (kcol)
    {
       case 1:
    kval = (kPDB->atom(i))->xcoor();
    break;
       case 2:
    kval = (kPDB->atom(i))->ycoor();
    break;
       case 3:
    kval = (kPDB->atom(i))->zcoor();
    break;
       case 4:
    kval = (kPDB->atom(i))->occupancy();
    break;
       case 5:
    kval = (kPDB->atom(i))->temperaturefactor();
    break;
    }
    
    if (kval > 0.0)
    {
       //  This atom is constrained
       consIndexes[i] = current_index;
       current_index++;
    }
    else
    {
       //  This atom is not constrained
       consIndexes[i] = -1;
    }
       }
       
       if (current_index == 0)
       {
    //  Constraints were turned on, but there weren't really any constrained
    iout << iWARN << "NO CONSTRAINED ATOMS WERE FOUND, BUT CONSTRAINTS ARE ON . . .\n" << endi;
       }
       else
       {
    //  Allocate an array to hold the constraint parameters
    consParams = new ConstraintParams[current_index];
    
    if (consParams == NULL)
    {
       NAMD_die("memory allocation failed in Molecule::build_constraint_params");
    }
       }
       
       numConstraints = current_index;
       
       //  Loop through all the atoms and assign the parameters for those
       //  that are constrained
       for (i=0; i<numAtoms; i++)
       {
    if (consIndexes[i] != -1)
    {
       //  This atom is constrained, so get the k value again
       switch (kcol)
       {
          case 1:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->xcoor();
       break;
          case 2:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->ycoor();
       break;
          case 3:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->zcoor();
       break;
          case 4:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->occupancy();
       break;
          case 5:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->temperaturefactor();
       break;
       }
       
       //  Get the reference position
       consParams[consIndexes[i]].refPos.x = (refPDB->atom(i))->xcoor();
       consParams[consIndexes[i]].refPos.y = (refPDB->atom(i))->ycoor();
       consParams[consIndexes[i]].refPos.z = (refPDB->atom(i))->zcoor();
    }
       }
       
       //  If we had to create new PDB objects, delete them now
       if (consref != NULL)
       {
    delete refPDB;
       }
       
       if ((conskfile != NULL) &&
     !((consref != NULL) && 
       (strcasecmp(consref->data, conskfile->data) == 0)
      )
    )
       {
    delete kPDB;
       }

    }
    /*      END OF FUNCTION build_constraint_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_movdrag_params  */
/*                  */
/*   INPUTS:        */
/*  movDragFile - value of movDragFile from the config file */
/*  movDragCol - value of movDragCol from the config file */
/*  movDragVelFile - value of movDragVelFile from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do moving drag. This involves looking through one or more    */
/*  PDB objects to determine which atoms are dragged,  and what the */
/*  drag parameters for each atom are. This information is then stored */
/*  in the arrays movDragIndexes and movDragParams. */
/*                  */
/************************************************************************/

void Molecule::build_movdrag_params(StringList *movDragFile, 
				    StringList *movDragCol, 
				    StringList *movDragVelFile, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *vPDB;        //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for drag tag in
  Real dtval = 0;          //  Drag tag value retreived
  char mainfilename[129];  //  main moving drag PDB filename
  char velfilename[129];   //  moving drag velocity PDB filename
  
  //  Get the PDB to read the moving drag tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (movDragFile == NULL) {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, movDragFile required.");
    tPDB = initial_pdb;
    
  } else {

    if (movDragFile->next != NULL) {
      NAMD_die("Multiple definitions of moving drag tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (movDragFile->data[0] == '/') ) {
      strcpy(mainfilename, movDragFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, movDragFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_movdrag_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in moving drag tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom velocities. If no name given, use
  // movDragFile if it is defined. Can NOT use the PDB coordinate
  // file!
  
  if (movDragVelFile == NULL) {
    if (movDragFile == NULL) {
      NAMD_die("Moving drag velocity file can not be same as coordinate PDB file");
    } else {
      if (movDragVelFile->next != NULL) {
	NAMD_die("Multiple definitions of moving drag velocity file in configuration file");
      };
      vPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (movDragVelFile->data[0] == '/') ) {
      strcpy(velfilename, movDragVelFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, movDragVelFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_movdrag_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in moving drag velocity PDB doesn't match coordinate PDB");
    }
  };
  
  
  //  Get the column that the drag tag is going to be in. If
  //  movDragFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If movDragFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (movDragCol == NULL) {
    dtcol = 4;
  } else {
    if (movDragCol->next != NULL) {
      NAMD_die("Multiple definitions of drag column in config file");
    };
    
    if (movDragFile == NULL
	&& strcasecmp(movDragCol->data, "B")
	&& strcasecmp(movDragCol->data, "O")) {
      NAMD_die("Can not read moving drag tags from X, Y, or Z column of the coordinate or velocity file");
    };
    if (!strcasecmp(movDragCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(movDragCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(movDragCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(movDragCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(movDragCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("movDragCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the drag
  //  parameters for each atom.  If the atom is not dragged, its
  //  value will be set to -1 in this array.
  movDragIndexes = new int32[numAtoms];
    if (movDragIndexes == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_movdrag_params()");
  };
  
  //  Loop through all the atoms and find out which ones are dragged
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is dragged
      movDragIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not dragged
      movDragIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    //  Drag was turned on, but there weren't really any dragged
    iout << iWARN << "NO DRAGGED ATOMS WERE FOUND, BUT MOVING DRAG IS ON . . . " << endi;
  } else {
    //  Allocate an array to hold the drag parameters
    movDragParams = new MovDragParams[current_index];
    if (movDragParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_movdrag_params");
    }
  };
  
  numMovDrag = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are dragged
  for (i=0; i<numAtoms; i++) {
    if (movDragIndexes[i] != -1) {
      movDragParams[movDragIndexes[i]].v[0] = (vPDB->atom(i))->xcoor();
      movDragParams[movDragIndexes[i]].v[1] = (vPDB->atom(i))->ycoor();
      movDragParams[movDragIndexes[i]].v[2] = (vPDB->atom(i))->zcoor();
    };
  };
      
  if (movDragFile != NULL) delete tPDB;
  if (movDragVelFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_movdrag_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_rotdrag_params  */
/*                  */
/*   INPUTS:        */
/*  rotDragFile - value of rotDragFile from the config file */
/*  rotDragCol - value of rotDragCol from the config file */
/*  rotDragAxisFile - value of rotDragAxisFile from the config file */
/*  rotDragPivotFile - value of rotDragPivotFile from the config file */
/*  rotDragVelFile - value of rotDragVelFile from the config file */
/*  rotDragVelCol - value of rotDragVelCol from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do moving drag. This involves looking through one or more    */
/*  PDB objects to determine which atoms are dragged,  and what the */
/*  drag parameters for each atom are. This information is then stored */
/*  in the arrays rotDragIndexes and rotDragParams. */
/*                  */
/************************************************************************/

void Molecule::build_rotdrag_params(StringList *rotDragFile, 
				    StringList *rotDragCol, 
				    StringList *rotDragAxisFile, 
				    StringList *rotDragPivotFile, 
				    StringList *rotDragVelFile, 
				    StringList *rotDragVelCol, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *aPDB, *pPDB, *vPDB; //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for drag tag in
  Real dtval = 0;          //  Drag tag value retreived
  int dvcol = 4;           //  Column to look for angular velocity in
  Real dvval = 0;          //  Angular velocity value retreived
  char mainfilename[129];  //  main rotating drag PDB filename
  char axisfilename[129];  //  rotating drag axis PDB filename
  char pivotfilename[129]; //  rotating drag pivot point PDB filename
  char velfilename[129];   //  rotating drag angular velocity PDB filename
  
  //  Get the PDB to read the rotating drag tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (rotDragFile == NULL) {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, rotDragFile required.");
    tPDB = initial_pdb;
    
  } else {

    if (rotDragFile->next != NULL) {
      NAMD_die("Multiple definitions of rotating drag tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (rotDragFile->data[0] == '/') ) {
      strcpy(mainfilename, rotDragFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, rotDragFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom rotation axes. If no name given, use
  // rotDragFile if both it AND rotDragPivotFile are defined. Can NOT
  // use the PDB coordinate file, nor rotDragPivotFile!

  if (rotDragAxisFile == NULL) {
    if (rotDragFile == NULL) {
      NAMD_die("Rotating drag axis file can not be same as coordinate PDB file");
    } else {
      if (rotDragAxisFile->next != NULL) {
	NAMD_die("Multiple definitions of rotating drag axis file in configuration file");
      };
      if (rotDragPivotFile == NULL) {
	NAMD_die("Need to specify at least one of rotDragAxisFile and rotDragPivotFile; they can not be same");
      };
      aPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (rotDragAxisFile->data[0] == '/') ) {
      strcpy(axisfilename, rotDragAxisFile->data);
    } else {
      strcpy(axisfilename, cwd);
      strcat(axisfilename, rotDragAxisFile->data);
    }
    
    aPDB = new PDB(axisfilename);
    if ( aPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (aPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag axis PDB doesn't match coordinate PDB");
    }
  };
  
  // Get the PDB to read atom rotation pivot points. If no name given,
  // use rotDragFile if both it AND rotDragAxisFile are defined. Can
  // NOT use the PDB coordinate file, nor rotDragAxisFile!

  if (rotDragPivotFile == NULL) {
    if (rotDragFile == NULL) {
      NAMD_die("Rotating drag pivot point file can not be same as coordinate PDB file");
    } else {
      if (rotDragPivotFile->next != NULL) {
	NAMD_die("Multiple definitions of rotating drag pivot point file in configuration file");
      };
      if (rotDragAxisFile == NULL) {
	NAMD_die("Need to specify at least one of rotDragAxisFile and rotDragPivotFile; they can not be same");
      };
      pPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (rotDragPivotFile->data[0] == '/') ) {
      strcpy(pivotfilename, rotDragPivotFile->data);
    } else {
      strcpy(pivotfilename, cwd);
      strcat(pivotfilename, rotDragPivotFile->data);
    }
    
    pPDB = new PDB(pivotfilename);
    if ( pPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (pPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag pivot point PDB doesn't match coordinate PDB");
    }
  };
  
  
  // Get the PDB to read atom angular velocities. If no name given,
  // use rotDragFile (or the coordinate PDB file if rotDragFile is not
  // defined).

  if (rotDragVelFile == NULL) {
    vPDB = tPDB;
  } else {
    if (rotDragVelFile->next != NULL) {
      NAMD_die("Multiple definitions of rotating drag velocity file in configuration file");
    };
    
    if ( (cwd == NULL) || (rotDragVelFile->data[0] == '/') ) {
      strcpy(velfilename, rotDragVelFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, rotDragVelFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag velocity PDB doesn't match coordinate PDB");
    }
  };
  
  //  Get the column that the drag tag is going to be in. If
  //  rotDragFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If rotDragFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (rotDragCol == NULL) {
    dtcol = 4;
  } else {
    if (rotDragCol->next != NULL) {
      NAMD_die("Multiple definitions of drag tag column in config file");
    };
    
    if ( rotDragFile == NULL
	 && (!strcasecmp(rotDragCol->data, "X")
	     || !strcasecmp(rotDragCol->data, "Y")
	     || !strcasecmp(rotDragCol->data, "Z"))) {
      NAMD_die("Can not read rotating drag tags from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(rotDragCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(rotDragCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(rotDragCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(rotDragCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(rotDragCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("rotDragCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Get the column that the drag angular velocity is going to be
  //  in. If rotDragVelFile is defined, it can be in any of the 5
  //  floating point fields in the PDB (X, Y, Z, O, or B) which
  //  correspond to the 1st, 2nd, ... 5th floating point fields. If
  //  NEITHER of rotDragVelFile OR rotDragFile is defined, it can
  //  only be O or B fileds. The default is the O (4th) field, which
  //  is the occupancy.

  if (rotDragVelCol == NULL) {
    dvcol = 4;
  } else {
    if (rotDragVelCol->next != NULL) {
      NAMD_die("Multiple definitions of drag angular velocity column in config file");
    };
    
    if (rotDragVelFile == NULL
	&& rotDragFile == NULL
	&& strcasecmp(rotDragCol->data, "B")
	&& strcasecmp(rotDragCol->data, "O")) {
      NAMD_die("Can not read rotating drag angular velocities from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(rotDragVelCol->data, "X")) {
      dvcol=1;
    } else if (!strcasecmp(rotDragVelCol->data, "Y")) {
      dvcol=2;
    } else if (!strcasecmp(rotDragVelCol->data, "Z")) {
      dvcol=3;
    } else if (!strcasecmp(rotDragVelCol->data, "O")) {
      dvcol=4;
    } else if (!strcasecmp(rotDragVelCol->data, "B")) {
      dvcol=5;
    }
    else {
      NAMD_die("rotDragVelCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the drag
  //  parameters for each atom.  If the atom is not dragged, its
  //  value will be set to -1 in this array.
  rotDragIndexes = new int32[numAtoms];
  if (rotDragIndexes == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_rotdrag_params()");
  };
  
  //  Loop through all the atoms and find out which ones are dragged
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is dragged
      rotDragIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not dragged
      rotDragIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    iout << iWARN << "NO DRAGGED ATOMS WERE FOUND, BUT ROTATING DRAG IS ON . . . " << endi;
  } else {
    rotDragParams = new RotDragParams[current_index];
    if (rotDragParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_rotdrag_params");
    }
  };
  
  numRotDrag = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are dragged
  for (i=0; i<numAtoms; i++) {
    if (rotDragIndexes[i] != -1) {
      rotDragParams[rotDragIndexes[i]].a[0] = (aPDB->atom(i))->xcoor();
      rotDragParams[rotDragIndexes[i]].a[1] = (aPDB->atom(i))->ycoor();
      rotDragParams[rotDragIndexes[i]].a[2] = (aPDB->atom(i))->zcoor();
      rotDragParams[rotDragIndexes[i]].p[0] = (pPDB->atom(i))->xcoor();
      rotDragParams[rotDragIndexes[i]].p[1] = (pPDB->atom(i))->ycoor();
      rotDragParams[rotDragIndexes[i]].p[2] = (pPDB->atom(i))->zcoor();
      switch (dvcol) {
      case 1:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->xcoor();
	break;
      case 2:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->ycoor();
	break;
      case 3:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->zcoor();
	break;
      case 4:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->occupancy();
	break;
      case 5:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->temperaturefactor();
	break;
      };
    };
  };
      
  if (rotDragFile != NULL) delete tPDB;
  if (rotDragAxisFile != NULL) delete aPDB;
  if (rotDragPivotFile != NULL) delete pPDB;
  if (rotDragVelFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_rotdrag_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_constorque_params  */
/*                  */
/*   INPUTS:        */
/*  consTorqueFile - value of consTorqueFile from the config file */
/*  consTorqueCol - value of consTorqueCol from the config file */
/*  consTorqueAxisFile - value of consTorqueAxisFile from the config file */
/*  consTorquePivotFile - value of consTorquePivotFile from the config file */
/*  consTorqueValFile - value of consTorqueValFile from the config file */
/*  consTorqueValCol - value of consTorqueValCol from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do "constant" torque. This involves looking through one or more    */
/*  PDB objects to determine which atoms are torqued,  and what the */
/*  torque parameters for each atom are. This information is then stored */
/*  in the arrays consTorqueIndexes and consTorqueParams. */
/*                  */
/************************************************************************/

void Molecule::build_constorque_params(StringList *consTorqueFile, 
				    StringList *consTorqueCol, 
				    StringList *consTorqueAxisFile, 
				    StringList *consTorquePivotFile, 
				    StringList *consTorqueValFile, 
				    StringList *consTorqueValCol, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *aPDB, *pPDB, *vPDB; //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for torque tag in
  Real dtval = 0;          //  Torque tag value retreived
  int dvcol = 4;           //  Column to look for angular velocity in
  Real dvval = 0;          //  Angular velocity value retreived
  char mainfilename[129];  //  main "constant" torque PDB filename
  char axisfilename[129];  //  "constant" torque axis PDB filename
  char pivotfilename[129]; //  "constant" torque pivot point PDB filename
  char velfilename[129];   //  "constant" torque angular velocity PDB filename
  
  //  Get the PDB to read the "constant" torque tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (consTorqueFile == NULL) {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, consTorqueFile required.");
    tPDB = initial_pdb;
    
  } else {

    if (consTorqueFile->next != NULL) {
      NAMD_die("Multiple definitions of \"constant\" torque tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (consTorqueFile->data[0] == '/') ) {
      strcpy(mainfilename, consTorqueFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, consTorqueFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom rotation axes. If no name given, use
  // consTorqueFile if both it AND consTorquePivotFile are defined. Can NOT
  // use the PDB coordinate file, nor consTorquePivotFile!

  if (consTorqueAxisFile == NULL) {
    if (consTorqueFile == NULL) {
      NAMD_die("\"Constant\" torque axis file can not be same as coordinate PDB file");
    } else {
      if (consTorqueAxisFile->next != NULL) {
	NAMD_die("Multiple definitions of \"constant\" torque axis file in configuration file");
      };
      if (consTorquePivotFile == NULL) {
	NAMD_die("Need to specify at least one of consTorqueAxisFile and consTorquePivotFile; they can not be same");
      };
      aPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (consTorqueAxisFile->data[0] == '/') ) {
      strcpy(axisfilename, consTorqueAxisFile->data);
    } else {
      strcpy(axisfilename, cwd);
      strcat(axisfilename, consTorqueAxisFile->data);
    }
    
    aPDB = new PDB(axisfilename);
    if ( aPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (aPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque axis PDB doesn't match coordinate PDB");
    }
  };
  
  // Get the PDB to read atom rotation pivot points. If no name given,
  // use consTorqueFile if both it AND consTorqueAxisFile are defined. Can
  // NOT use the PDB coordinate file, nor consTorqueAxisFile!

  if (consTorquePivotFile == NULL) {
    if (consTorqueFile == NULL) {
      NAMD_die("\"Constant\" torque pivot point file can not be same as coordinate PDB file");
    } else {
      if (consTorquePivotFile->next != NULL) {
	NAMD_die("Multiple definitions of \"constant\" torque pivot point file in configuration file");
      };
      if (consTorqueAxisFile == NULL) {
	NAMD_die("Need to specify at least one of consTorqueAxisFile and consTorquePivotFile; they can not be same");
      };
      pPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (consTorquePivotFile->data[0] == '/') ) {
      strcpy(pivotfilename, consTorquePivotFile->data);
    } else {
      strcpy(pivotfilename, cwd);
      strcat(pivotfilename, consTorquePivotFile->data);
    }
    
    pPDB = new PDB(pivotfilename);
    if ( pPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (pPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque pivot point PDB doesn't match coordinate PDB");
    }
  };
  
  
  // Get the PDB to read atom angular velocities. If no name given,
  // use consTorqueFile (or the coordinate PDB file if consTorqueFile is not
  // defined).

  if (consTorqueValFile == NULL) {
    vPDB = tPDB;
  } else {
    if (consTorqueValFile->next != NULL) {
      NAMD_die("Multiple definitions of \"constant\" torque velocity file in configuration file");
    };
    
    if ( (cwd == NULL) || (consTorqueValFile->data[0] == '/') ) {
      strcpy(velfilename, consTorqueValFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, consTorqueValFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque velocity PDB doesn't match coordinate PDB");
    }
  };
  
  //  Get the column that the torque tag is going to be in. If
  //  consTorqueFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If consTorqueFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (consTorqueCol == NULL) {
    dtcol = 4;
  } else {
    if (consTorqueCol->next != NULL) {
      NAMD_die("Multiple definitions of torque tag column in config file");
    };
    
    if ( consTorqueFile == NULL
	 && (!strcasecmp(consTorqueCol->data, "X")
	     || !strcasecmp(consTorqueCol->data, "Y")
	     || !strcasecmp(consTorqueCol->data, "Z"))) {
      NAMD_die("Can not read \"constant\" torque tags from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(consTorqueCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(consTorqueCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(consTorqueCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(consTorqueCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(consTorqueCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("consTorqueCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Get the column that the torque value is going to be
  //  in. If consTorqueValFile is defined, it can be in any of the 5
  //  floating point fields in the PDB (X, Y, Z, O, or B) which
  //  correspond to the 1st, 2nd, ... 5th floating point fields. If
  //  NEITHER of consTorqueValFile OR consTorqueFile is defined, it can
  //  only be O or B fileds. The default is the O (4th) field, which
  //  is the occupancy.

  if (consTorqueValCol == NULL) {
    dvcol = 4;
  } else {
    if (consTorqueValCol->next != NULL) {
      NAMD_die("Multiple definitions of torque value column in config file");
    };
    
    if (consTorqueValFile == NULL
	&& consTorqueFile == NULL
	&& strcasecmp(consTorqueCol->data, "B")
	&& strcasecmp(consTorqueCol->data, "O")) {
      NAMD_die("Can not read \"constant\" torque values from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(consTorqueValCol->data, "X")) {
      dvcol=1;
    } else if (!strcasecmp(consTorqueValCol->data, "Y")) {
      dvcol=2;
    } else if (!strcasecmp(consTorqueValCol->data, "Z")) {
      dvcol=3;
    } else if (!strcasecmp(consTorqueValCol->data, "O")) {
      dvcol=4;
    } else if (!strcasecmp(consTorqueValCol->data, "B")) {
      dvcol=5;
    }
    else {
      NAMD_die("consTorqueValCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the torque
  //  parameters for each atom.  If the atom is not torqued, its
  //  value will be set to -1 in this array.
  consTorqueIndexes = new int32[numAtoms];
  if (consTorqueIndexes == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_constorque_params()");
  };
  
  //  Loop through all the atoms and find out which ones are torqued
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is torqued
      consTorqueIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not torqued
      consTorqueIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    iout << iWARN << "NO TORQUED ATOMS WERE FOUND, BUT \"CONSTANT\" TORQUE IS ON . . . " << endi;
  } else {
    consTorqueParams = new ConsTorqueParams[current_index];
    if (consTorqueParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_constorque_params");
    }
  };
  
  numConsTorque = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are torqued
  for (i=0; i<numAtoms; i++) {
    if (consTorqueIndexes[i] != -1) {
      consTorqueParams[consTorqueIndexes[i]].a[0] = (aPDB->atom(i))->xcoor();
      consTorqueParams[consTorqueIndexes[i]].a[1] = (aPDB->atom(i))->ycoor();
      consTorqueParams[consTorqueIndexes[i]].a[2] = (aPDB->atom(i))->zcoor();
      consTorqueParams[consTorqueIndexes[i]].p[0] = (pPDB->atom(i))->xcoor();
      consTorqueParams[consTorqueIndexes[i]].p[1] = (pPDB->atom(i))->ycoor();
      consTorqueParams[consTorqueIndexes[i]].p[2] = (pPDB->atom(i))->zcoor();
      switch (dvcol) {
      case 1:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->xcoor();
	break;
      case 2:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->ycoor();
	break;
      case 3:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->zcoor();
	break;
      case 4:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->occupancy();
	break;
      case 5:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->temperaturefactor();
	break;
      };
    };
  };
      
  if (consTorqueFile != NULL) delete tPDB;
  if (consTorqueAxisFile != NULL) delete aPDB;
  if (consTorquePivotFile != NULL) delete pPDB;
  if (consTorqueValFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_constorque_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_constant_forces    */
/*                  */
/*   INPUTS:                */
/*  filename - PDB file containing the constant forces    */
/*                  */
/*  This function reads the constant forces from the PDB file.  */
/*   The force vector to be applied on each atom is determined by:    */
/*     occupancy*(X,Y,Z)   */
/*   Only non-zero forces are stored   */
/*                  */
/************************************************************************/

void Molecule::build_constant_forces(char *filename)
{ int i, index;
  PDB *forcePDB;
  
  if (!filename) {
    // then all forces are zero to begin with; may be changed by
    // the consforceconfig command.
    iout << iWARN << "NO CONSTANT FORCES SPECIFIED, BUT CONSTANT FORCE IS ON . . .\n" << endi;
    consForceIndexes = new int32[numAtoms];
    for (i=0; i<numAtoms; i++) consForceIndexes[i] = -1;
    return;
  }

  if ((forcePDB=new PDB(filename)) == NULL)
    NAMD_die("Memory allocation failed in Molecule::build_constant_forces");
  if (forcePDB->num_atoms() != numAtoms)
    NAMD_die("Number of atoms in constant force PDB doesn't match coordinate PDB");

  //  Allocate an array that will store an index into the constant force
  //  array for each atom.  If the atom has no constant force applied, its
  //  value will be set to -1 in this array.
  consForceIndexes = new int32[numAtoms];
  if (consForceIndexes == NULL)
    NAMD_die("memory allocation failed in Molecule::build_constant_forces()");

  //  Loop through all the atoms and find out which ones have constant force
  numConsForce = 0;
  for (i=0; i<numAtoms; i++)
    if ((forcePDB->atom(i)->xcoor()==0 && forcePDB->atom(i)->ycoor()==0 &&
         forcePDB->atom(i)->zcoor()==0) || forcePDB->atom(i)->occupancy()==0)
      //  This atom has no constant force
      consForceIndexes[i] = -1;
    else
      //  This atom has constant force
      consForceIndexes[i] = numConsForce++;

  if (numConsForce == 0)
    // Constant force was turned on, but there weren't really any non-zero forces
    iout << iWARN << "NO NON-ZERO FORCES WERE FOUND, BUT CONSTANT FORCE IS ON . . .\n" << endi;
  else
  { // Allocate an array to hold the forces
    consForce = new Vector[numConsForce];
    if (consForce == NULL)
      NAMD_die("memory allocation failed in Molecule::build_constant_forces");
    // Loop through all the atoms and assign the forces
    for (i=0; i<numAtoms; i++)
      if ((index=consForceIndexes[i]) != -1)
      { //  This atom has constant force on it
        consForce[index].x = forcePDB->atom(i)->xcoor() * forcePDB->atom(i)->occupancy();
        consForce[index].y = forcePDB->atom(i)->ycoor() * forcePDB->atom(i)->occupancy();
        consForce[index].z = forcePDB->atom(i)->zcoor() * forcePDB->atom(i)->occupancy();
      }
  }

  delete forcePDB;
}
/*      END OF FUNCTION build_constant_forces    */


void Molecule::build_langevin_params(BigReal coupling,
    BigReal drudeCoupling, Bool doHydrogen) {

  //  Allocate the array to hold all the data
  langevinParams = new Real[numAtoms];

  if ( (langevinParams == NULL) )
  {
    NAMD_die("memory allocation failed in Molecule::build_langevin_params()");
  }

  //  Loop through all the atoms and get the b value
  for (int i=0; i<numAtoms; i++)
  {
    BigReal bval = coupling;

    if ( (! doHydrogen) && is_hydrogen(i) ) bval = 0;
    else if ( is_lp(i) ) bval = 0;
    else if ( is_drude(i) ) bval = drudeCoupling;

    //  Assign the b value
    langevinParams[i] = bval;
  }

}

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_langevin_params      */
    /*                  */
    /*   INPUTS:                */
    /*  langfile - Value of langevinfile from config file    */
    /*  langcol - Value of langevincol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*      cwd - Current working directory          */
    /*                  */
    /*  This function builds the array of b values necessary for  */
    /*   Langevin dynamics.  It takes the name of the PDB file and the      */
    /*   column in the PDB file that contains the b values.  It then  */
    /*   builds the array langevinParams for use during the program.  */
    /*                  */
    /************************************************************************/

    void Molecule::build_langevin_params(StringList *langfile, 
           StringList *langcol, 
           PDB *initial_pdb,
           char *cwd)
       
    {
       PDB *bPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  
       if (langfile == NULL)
       {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, langevinFile required.");
    bPDB = initial_pdb;
       }
       else
       {
    if (langfile->next != NULL)
    {
       NAMD_die("Multiple definitions of langvein PDB file in configuration file");
    }

    if ( (cwd == NULL) || (langfile->data[0] == '/') )
    {
         strcpy(filename, langfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, langfile->data);
    }
    
    bPDB = new PDB(filename);
    if ( bPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_langevin_params");
    }
    
    if (bPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in langevin parameter PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the column that the b vaules are in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (langcol == NULL)
       {
    bcol = 4;
       }
       else
       {
    if (langcol->next != NULL)
    {
       NAMD_die("Multiple definitions of langevin parameter column in config file");
    }
    
    if (strcasecmp(langcol->data, "X") == 0)
    {
       bcol=1;
    }
    else if (strcasecmp(langcol->data, "Y") == 0)
    {
       bcol=2;
    }
    else if (strcasecmp(langcol->data, "Z") == 0)
    {
       bcol=3;
    }
    else if (strcasecmp(langcol->data, "O") == 0)
    {
       bcol=4;
    }
    else if (strcasecmp(langcol->data, "B") == 0)
    {
       bcol=5;
    }
    else
    {
       NAMD_die("langevincol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate the array to hold all the data
       langevinParams = new Real[numAtoms];
       
       if ( (langevinParams == NULL) )
       {
    NAMD_die("memory allocation failed in Molecule::build_langevin_params()");
       }

       //  Loop through all the atoms and get the b value
       for (i=0; i<numAtoms; i++)
       {
    //  Get the k value based on where we were told to find it
    switch (bcol)
    {
       case 1:
    bval = (bPDB->atom(i))->xcoor();
    break;
       case 2:
    bval = (bPDB->atom(i))->ycoor();
    break;
       case 3:
    bval = (bPDB->atom(i))->zcoor();
    break;
       case 4:
    bval = (bPDB->atom(i))->occupancy();
    break;
       case 5:
    bval = (bPDB->atom(i))->temperaturefactor();
    break;
    }
    
    //  Assign the b value
    langevinParams[i] = bval;
       }
       
       //  If we had to create a PDB object, delete it now
       if (langfile != NULL)
       {
    delete bPDB;
       }
    }
    /*      END OF FUNCTION build_langevin_params    */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_fixed_atoms      */
    /*                  */
    /*   INPUTS:              */
    /*  fixedfile - Value of langevinfile from config file    */
    /*  fixedcol - Value of langevincol from config file    */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*      cwd - Current working directory        */
    /*                  */
    /*  This function builds the list of fixed atoms.      */
    /*   It takes the name of the PDB file and the      */
    /*   column in the PDB file that contains the flags.  It then  */
    /*   builds the array fixedAtomFlags for use during the program.  */
    /*                  */
    /************************************************************************/

    void Molecule::build_fixed_atoms(StringList *fixedfile, 
           StringList *fixedcol, 
           PDB *initial_pdb,
           char *cwd)
       
{
       PDB *bPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  
       if (fixedfile == NULL)
       {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, fixedAtomsFile required.");
    bPDB = initial_pdb;
       }
       else
       {
    if (fixedfile->next != NULL)
    {
       NAMD_die("Multiple definitions of fixed atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (fixedfile->data[0] == '/') )
    {
         strcpy(filename, fixedfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, fixedfile->data);
    }
    
    bPDB = new PDB(filename);
    if ( bPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_fixed_atoms");
    }
    
    if (bPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the column that the b vaules are in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (fixedcol == NULL)
       {
    bcol = 4;
       }
       else
       {
    if (fixedcol->next != NULL)
    {
       NAMD_die("Multiple definitions of fixed atoms column in config file");
    }
    
    if (strcasecmp(fixedcol->data, "X") == 0)
    {
       bcol=1;
    }
    else if (strcasecmp(fixedcol->data, "Y") == 0)
    {
       bcol=2;
    }
    else if (strcasecmp(fixedcol->data, "Z") == 0)
    {
       bcol=3;
    }
    else if (strcasecmp(fixedcol->data, "O") == 0)
    {
       bcol=4;
    }
    else if (strcasecmp(fixedcol->data, "B") == 0)
    {
       bcol=5;
    }
    else
    {
       NAMD_die("fixedatomscol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate the array to hold all the data
       fixedAtomFlags = new int32[numAtoms];
       
       if (fixedAtomFlags == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_fixed_atoms()");
       }
       
  numFixedAtoms = 0;

       //  Loop through all the atoms and get the b value
       for (i=0; i<numAtoms; i++)
       {
    //  Get the k value based on where we were told to find it
    switch (bcol)
    {
       case 1:
    bval = (bPDB->atom(i))->xcoor();
    break;
       case 2:
    bval = (bPDB->atom(i))->ycoor();
    break;
       case 3:
    bval = (bPDB->atom(i))->zcoor();
    break;
       case 4:
    bval = (bPDB->atom(i))->occupancy();
    break;
       case 5:
    bval = (bPDB->atom(i))->temperaturefactor();
    break;
    }
    
    //  Assign the b value
    if ( bval != 0 ) {
      fixedAtomFlags[i] = 1;
      numFixedAtoms++;
    }
    else {
      fixedAtomFlags[i] = 0;
    }
       }

       //  If we had to create a PDB object, delete it now
       if (fixedfile != NULL)
       {
    delete bPDB;
       }

  // now figure out how we interact with rigidBonds 
  // this is mainly for degree of freedom counting
  if ( numRigidBonds ) {
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    int parentIsFixed = 0;
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) {
	parentIsFixed = fixedAtomFlags[h_i->atomID];
	if ( (rigidBondLengths[h_i->atomID] > 0.)  // water
		&& fixedAtomFlags[h_i[1].atomID]
		&& fixedAtomFlags[h_i[2].atomID] ) {
	  ++numFixedRigidBonds;
	}
      } else {
	if ( (rigidBondLengths[h_i->atomID] > 0.)
		&& fixedAtomFlags[h_i->atomID]
		&& parentIsFixed ) {
	  ++numFixedRigidBonds;
	}
      }
    }
  }

  // how many hydrogen groups are completely fixed
  {
    numFixedGroups = 0;
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    int allFixed = 0;
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) {
        if ( allFixed ) ++numFixedGroups;
        allFixed = 1;
      }
      allFixed = allFixed && fixedAtomFlags[h_i->atomID];
    }
    if ( allFixed ) ++numFixedGroups;
  }

}
    /*      END OF FUNCTION build_fixed_atoms    */




/************************************************************************/
/*                                                                      */
/*      FUNCTION build_stirred_atoms                                    */
/*                                                                      */
/*   INPUTS:                                                            */
/*  stirredfile - Value of stirFilename from config file    */
/*  stirredcol - Value of stircol from config file (but B, O only */
/*                    since we need X, Y, Z!     */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory        */
/*                  */
/*  This function builds the list of fixed atoms.      */
/*   It takes the name of the PDB file and the      */
/*   column in the PDB file that contains the flags.  It then  */
/*   builds the array fixedAtomFlags for use during the program.  */
/*                                                                      */
/************************************************************************/

    void Molecule::build_stirred_atoms(StringList *stirredfile, 
           StringList *stirredcol, 
           PDB *initial_pdb,
           char *cwd)
       
{
       PDB *sPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       int current_index=0;    //  Index into values used
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates. 
       //  use posssible only if this is 'optional' in simulation parameters
       // dangerous, since restarted simulations will be different
       if (stirredfile == NULL)
       {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, stirFilename required.");
    sPDB = initial_pdb;
    // dangerous, since restarted simulations will be different, so warn
        iout << iWARN << "STIRRING USING INITIAL POSITION FILE FOR REFERENCE POSITIONS" << endi;
       }
       else
       {
    if (stirredfile->next != NULL)
    {
       NAMD_die("Multiple definitions of stirred atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (stirredfile->data[0] == '/') )
    {
         strcpy(filename, stirredfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, stirredfile->data);
    }
        
    //CkPrintf ("DEBUG: the stir filename is %s\n",filename);    
    sPDB = new PDB(filename);

    if ( sPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_stirred_atoms");

    }

    if (sPDB->num_atoms() != numAtoms)
      {
	NAMD_die("Number of atoms in stirred atoms PDB doesn't match coordinate PDB");
      }
 
       }

//  Get the column that the b vaules are in.  It
//  can be in any of the 5 floating point fields in the PDB, according
//  to what the user wants.  The allowable fields are X, Y, Z, O, or
//  B which correspond to the 1st, 2nd, ... 5th floating point fields.
//  The default is the 4th field, which is the occupancy

 
      if (stirredcol == NULL)
      {
	 bcol = 4;
      }
      else
    {
      if (stirredcol->next != NULL)
	{
	  NAMD_die("Multiple definitions of stirred atoms column in config file");
	}
      
      if  (strcasecmp(stirredcol->data, "O") == 0)
	{
	  bcol=4;
	}
      else if (strcasecmp(stirredcol->data, "B") == 0)
	{
	  bcol=5;
	}
      else
	{
	  NAMD_die("stirredAtomsCol must have value of O or B");
	}
    }
          
       //  Allocate an array that will store an index into the stir
       //  parameters for each atom.  If the atom is not stirred, its
       //  value will be set to -1 in this array.
       stirIndexes = new int32[numAtoms];
       
       if (stirIndexes == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_stirred_params()");
       }
       
       current_index = 0;
       //  Loop through all the atoms and find out which ones are stirred
       for (i=0; i<numAtoms; i++)
       {



	 //  Get the b value based on where we were told to find it
	 switch (bcol)
	   {

	   case 4:
	     bval = (sPDB->atom(i))->occupancy();
	     break;
	   case 5:
	     bval = (sPDB->atom(i))->temperaturefactor();
	     break;
	   }

	// CkPrintf ("DEBUG--for atom i= %d  bcol= %d    bval= %g   occ= %g   realbval= %g  x= %g numAtoms= %d test= %g\n" ,i        ,bcol, bval, (sPDB->atom(i))->occupancy(), (sPDB->atom(i))->temperaturefactor(),(sPDB->atom(i))->xcoor(), sPDB->num_atoms(), 0.123 );
	 //  Assign the b value
	 if ( bval != 0 )
	 {
	   // This atom is stirred 
	   stirIndexes[i] = current_index;
	   current_index++;
	 }
	 else
	 {
           //This atom is not stirred 
	   stirIndexes[i] = -1;
	 }
       }
	 

    

       
       if (current_index == 0)
       {
    //  Stirring was turned on, but there weren't really any stirred atoms found in file
    iout << iWARN << "NO STIRRED ATOMS WERE FOUND, BUT STIRRING TORQUES ARE ON . . .\n" << endi;
       }
       else
       {
    //  Allocate an array to hold the stirring parameters
    stirParams = new StirParams[current_index];
    
    if (stirParams == NULL)
    {
       NAMD_die("memory allocation failed in Molecule::build_stir_params");
    }
       }
       
       numStirredAtoms = current_index;
       
       //  Loop through all the atoms and assign the parameters for those
       //  that are stirred
       for (i=0; i<numAtoms; i++)
       {
    if (stirIndexes[i] != -1)
    {
       
       //  This atom is stirred, so get the reference position
       stirParams[stirIndexes[i]].refPos.x = (sPDB->atom(i))->xcoor();
       stirParams[stirIndexes[i]].refPos.y = (sPDB->atom(i))->ycoor();
       stirParams[stirIndexes[i]].refPos.z = (sPDB->atom(i))->zcoor();
    }
       }
       
       //  If we had to create a PDB object, delete it now
       if (stirredfile != NULL)
       {
	 delete sPDB;
       }
       

    }

    /*      END OF FUNCTION build_stirred_atoms    */



void Molecule::build_extra_bonds(Parameters *parameters, StringList *file) {
//In the memory optimized version, only the parameters of extraBonds are needed
//to load
  char err_msg[512];
  int a1,a2,a3,a4; float k, ref;
  #ifndef MEM_OPT_VERSION
  ResizeArray<Bond> bonds;
  ResizeArray<Angle> angles;
  ResizeArray<Dihedral> dihedrals;
  ResizeArray<Improper> impropers;
  #endif
  ResizeArray<BondValue> bond_params;
  ResizeArray<AngleValue> angle_params;
  ResizeArray<DihedralValue> dihedral_params;
  ResizeArray<ImproperValue> improper_params;
  ResizeArray<GromacsPairValue> gromacsPair_params;

  if ( ! file ) {
    NAMD_die("NO EXTRA BONDS FILES SPECIFIED");
  }

  for ( ; file; file = file->next ) {  // loop over files
    FILE *f = fopen(file->data,"r");
    if ( ! f ) {
      sprintf(err_msg, "UNABLE TO OPEN EXTRA BONDS FILE %s", file->data);
      NAMD_err(err_msg);
    } else {
      iout << iINFO << "READING EXTRA BONDS FILE " << file->data <<"\n"<<endi;
    }
    
    while ( 1 ) {
      char buffer[512];
      int ret_code;
      do {
        ret_code = NAMD_read_line(f, buffer);
      } while ( (ret_code==0) && (NAMD_blank_string(buffer)) );
      if (ret_code!=0) break;

      char type[512];
      sscanf(buffer,"%s",type);

#define CHECKATOMID(ATOMID) if ( ATOMID < 0 || ATOMID >= numAtoms ) badatom = 1;

      int badline = 0;
      int badatom = 0;
      if ( ! strncasecmp(type,"bond",4) ) {
        if ( sscanf(buffer, "%s %d %d %f %f %s",
	    type, &a1, &a2, &k, &ref, err_msg) != 5 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
        }

        #ifndef MEM_OPT_VERSION              
        Bond tmp;
        tmp.bond_type = parameters->NumBondParams + bonds.size();
        tmp.atom1 = a1;  tmp.atom2 = a2;
        bonds.add(tmp);
        #endif

        BondValue tmpv;
        tmpv.k = k;  tmpv.x0 = ref;
        bond_params.add(tmpv);                
      } else if ( ! strncasecmp(type,"angle",4) ) {
        if ( sscanf(buffer, "%s %d %d %d %f %f %s",
	    type, &a1, &a2, &a3, &k, &ref, err_msg) != 6 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
          CHECKATOMID(a3)
        }
        #ifndef MEM_OPT_VERSION
        Angle tmp;
        tmp.atom1 = a1;  tmp.atom2 = a2;  tmp.atom3 = a3;
        tmp.angle_type = parameters->NumAngleParams + angles.size();
        angles.add(tmp);  
        #endif  

        AngleValue tmpv;
        tmpv.k = k;  tmpv.theta0 = ref / 180. * PI;
        tmpv.k_ub = 0;  tmpv.r_ub = 0;
        angle_params.add(tmpv);      
              
      } else if ( ! strncasecmp(type,"dihedral",4) ) {
        int n = 0;
        int ret = 1 + sscanf(buffer, "%s %d %d %d %d %f %f %s",
	                 type, &a1, &a2, &a3, &a4, &k, &ref, err_msg);
        if ( ret != 8 ) {
          ret = sscanf(buffer, "%s %d %d %d %d %f %d %f %s",
	                 type, &a1, &a2, &a3, &a4, &k, &n, &ref, err_msg);
        }
        if ( ret != 8 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
          CHECKATOMID(a3)
          CHECKATOMID(a4)
        }
        #ifndef MEM_OPT_VERSION
        Dihedral tmp;
        tmp.atom1 = a1;  tmp.atom2 = a2;  tmp.atom3 = a3;  tmp.atom4 = a4;
        tmp.dihedral_type = parameters->NumDihedralParams + dihedrals.size();
        dihedrals.add(tmp);
        #endif

        DihedralValue tmpv;
        tmpv.multiplicity = 1;  tmpv.values[0].n = n;
        tmpv.values[0].k = k;  tmpv.values[0].delta = ref / 180. * PI;
        dihedral_params.add(tmpv);
      } else if ( ! strncasecmp(type,"improper",4) ) {
        int n = 0;
        int ret = 1 + sscanf(buffer, "%s %d %d %d %d %f %f %s",
	                 type, &a1, &a2, &a3, &a4, &k, &ref, err_msg);
        if ( ret != 8 ) {
          ret = sscanf(buffer, "%s %d %d %d %d %f %d %f %s",
	                 type, &a1, &a2, &a3, &a4, &k, &n, &ref, err_msg);
        }
        if ( ret != 8 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
          CHECKATOMID(a3)
          CHECKATOMID(a4)
        }
        #ifndef MEM_OPT_VERSION
        Improper tmp;
        tmp.atom1 = a1;  tmp.atom2 = a2;  tmp.atom3 = a3;  tmp.atom4 = a4;
        tmp.improper_type = parameters->NumImproperParams + impropers.size();
        impropers.add(tmp);  
        #endif

        ImproperValue tmpv;
        tmpv.multiplicity = 1;  tmpv.values[0].n = n;
        tmpv.values[0].k = k;  tmpv.values[0].delta = ref / 180. * PI;
        improper_params.add(tmpv);
      } else if ( ! strncasecmp(type,"#",1) ) {
        continue;  // comment
      } else {
        badline = 1;
      }
#undef CHECKATOMID
      if ( badline ) {
        sprintf(err_msg, "BAD LINE IN EXTRA BONDS FILE %s: %s",
						file->data, buffer);
        NAMD_die(err_msg);
      }
      if ( badatom ) {
        sprintf(err_msg, "BAD ATOM ID IN EXTRA BONDS FILE %s: %s",
						file->data, buffer);
        NAMD_die(err_msg);
      }
    }
    fclose(f);
  }  // loop over files

  // append to parameters and molecule data structures
  int extraNumBonds = bond_params.size();
  if ( extraNumBonds ) {
    iout << iINFO << "READ " << extraNumBonds << " EXTRA BONDS\n" << endi;

    #ifndef MEM_OPT_VERSION
    Bond *newbonds = new Bond[numBonds+extraNumBonds];
    memcpy(newbonds, this->bonds, numBonds*sizeof(Bond));
    memcpy(newbonds+numBonds, bonds.begin(), extraNumBonds*sizeof(Bond));
    delete [] this->bonds;
    this->bonds = newbonds;
    numBonds += extraNumBonds;
    #endif

    BondValue *newbondp = new BondValue[
			parameters->NumBondParams + extraNumBonds];
    memcpy(newbondp, parameters->bond_array,
			parameters->NumBondParams * sizeof(BondValue));
    memcpy(newbondp+parameters->NumBondParams, bond_params.begin(),
			extraNumBonds * sizeof(BondValue));
    delete [] parameters->bond_array;
    parameters->bond_array = newbondp;
    parameters->NumBondParams += extraNumBonds;
  }

  int extraNumAngles = angle_params.size();
  if ( extraNumAngles ) {
    iout << iINFO << "READ " << extraNumAngles << " EXTRA ANGLES\n" << endi;
    #ifndef MEM_OPT_VERSION
    Angle *newangles = new Angle[numAngles+extraNumAngles];
    memcpy(newangles, this->angles, numAngles*sizeof(Angle));
    memcpy(newangles+numAngles, angles.begin(), extraNumAngles*sizeof(Angle));
    delete [] this->angles;
    this->angles = newangles;
    numAngles += extraNumAngles;
    #endif

    AngleValue *newanglep = new AngleValue[
			parameters->NumAngleParams + extraNumAngles];
    memcpy(newanglep, parameters->angle_array,
			parameters->NumAngleParams * sizeof(AngleValue));
    memcpy(newanglep+parameters->NumAngleParams, angle_params.begin(),
			extraNumAngles * sizeof(AngleValue));
    delete [] parameters->angle_array;
    parameters->angle_array = newanglep;
    parameters->NumAngleParams += extraNumAngles;
  }

  int extraNumDihedrals = dihedral_params.size();
  if ( extraNumDihedrals ) {
    iout << iINFO << "READ " << extraNumDihedrals << " EXTRA DIHEDRALS\n" << endi;
    #ifndef MEM_OPT_VERSION
    Dihedral *newdihedrals = new Dihedral[numDihedrals+extraNumDihedrals];
    memcpy(newdihedrals, this->dihedrals, numDihedrals*sizeof(Dihedral));
    memcpy(newdihedrals+numDihedrals, dihedrals.begin(), extraNumDihedrals*sizeof(Dihedral));
    delete [] this->dihedrals;
    this->dihedrals = newdihedrals;
    numDihedrals += extraNumDihedrals;
    #endif

    DihedralValue *newdihedralp = new DihedralValue[
			parameters->NumDihedralParams + extraNumDihedrals];
    memcpy(newdihedralp, parameters->dihedral_array,
			parameters->NumDihedralParams * sizeof(DihedralValue));
    memcpy(newdihedralp+parameters->NumDihedralParams, dihedral_params.begin(),
			extraNumDihedrals * sizeof(DihedralValue));
    delete [] parameters->dihedral_array;
    parameters->dihedral_array = newdihedralp;
    parameters->NumDihedralParams += extraNumDihedrals;
  }

  int extraNumImpropers = improper_params.size();
  if ( extraNumImpropers ) {
    iout << iINFO << "READ " << extraNumImpropers << " EXTRA IMPROPERS\n" << endi;
    #ifndef MEM_OPT_VERSION
    Improper *newimpropers = new Improper[numImpropers+extraNumImpropers];
    memcpy(newimpropers, this->impropers, numImpropers*sizeof(Improper));
    memcpy(newimpropers+numImpropers, impropers.begin(), extraNumImpropers*sizeof(Improper));
    delete [] this->impropers;
    this->impropers = newimpropers;
    numImpropers += extraNumImpropers;
    #endif

    ImproperValue *newimproperp = new ImproperValue[
			parameters->NumImproperParams + extraNumImpropers];
    memcpy(newimproperp, parameters->improper_array,
			parameters->NumImproperParams * sizeof(ImproperValue));
    memcpy(newimproperp+parameters->NumImproperParams, improper_params.begin(),
			extraNumImpropers * sizeof(ImproperValue));
    delete [] parameters->improper_array;
    parameters->improper_array = newimproperp;
    parameters->NumImproperParams += extraNumImpropers;
  }
}// end of Molecule::build_extra_bonds()


//Modifications for alchemical fep
/*
 FUNCTION build_fep_flags

 INPUTS:
 alchfile - Value of alchfile read from config file
 alchcol - Value of alch column, read from config file
 initial_pdb - PDB object that contains the initial positions
  cwd - current working directory

 This function builds the array of state values necessary
 for FEP or TI. It takes the name of the PDB file and column in
 the PDB file that contains the alch flag. It then builds
 the array FepParams for use in the program.

 function doubles up for TI as well 
*/
void Molecule::build_fep_flags(StringList *alchfile, StringList *alchcol,
                              PDB *initial_pdb, char *cwd, 
                              const char *simmethod) {
  PDB *bPDB;  //Pointer to PDB object to use
  int bcol = 5;  //Column that the data is in
  Real bval = 0; //flag from PDB file
  int i;         // loop counter
  char filename[129]; // filename

  // get the pdb object that contains the alch flags.
  // if the user gave another filename, use it, else
  // use the pdb file with the initial coordinates
  if (alchfile == NULL) {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, alchfile required.");
    bPDB = initial_pdb;
    strcpy(filename, "coordinate pdb file (default)");
  }
  else {
    if (alchfile->next != NULL) {
      char *new_err_msg = new char[24 + strlen(simmethod) + 26];
      sprintf(new_err_msg,"Multiple definitions of %sFile in configuration file",simmethod);
      NAMD_die(new_err_msg);
    }
   
    if ((cwd == NULL) || (alchfile->data[0] == '/')) {
      strcpy(filename, alchfile->data);
    }
    else {
      strcpy(filename, cwd);
      strcat(filename, alchfile->data);
    }

    bPDB = new PDB(filename);
    if (bPDB == NULL) {
      NAMD_die("Memory allocation failed in Molecule:build_fep_flags");
    }

    if (bPDB->num_atoms() != numAtoms) {
      char *new_err_msg = new char[19 + strlen(simmethod) + 38];
      sprintf(new_err_msg,"Number of atoms in %sFile PDB does not match coordinate PDB",simmethod);
      NAMD_die(new_err_msg);
    }
  }
   
  // Get the column that the alch flag is in. It can be in any of the 5 
  // floating point fields in the PDB ie X, Y, Z, O or B.
  // The default is 5th field ie the beta field
  if (alchcol == NULL) {
    bcol = 5;
  }
  else {
    if (alchcol->next != NULL) {
      char *new_err_msg = new char[24 + strlen(simmethod) + 35];
      sprintf(new_err_msg,"Multiple definitions of %s parameter column in config file",simmethod);
      NAMD_die(new_err_msg);
    }

    if (strcasecmp(alchcol->data, "X") == 0) {
      bcol = 1;
    }
    else if (strcasecmp(alchcol->data, "Y") == 0) {
      bcol = 2;
    }
    else if (strcasecmp(alchcol->data, "Z") == 0) {
      bcol = 3;
    }
    else if (strcasecmp(alchcol->data, "O") == 0) {
      bcol = 4;
    }
    else if (strcasecmp(alchcol->data, "B") == 0) {
      bcol = 5;
    }
    else {
      NAMD_die("alchcol must have value of X, Y, Z, O or B");
    }
  }

  iout << iINFO << "To read " << simmethod << "data from file: " << filename 
       << "\n" << endi;
  iout << iINFO << "To read " << simmethod << "flag data from column: " << bcol
       << "\n" << endi;
 
  //  Allocate the array to hold all the alch data
  fepAtomFlags = new unsigned char[numAtoms];
       
  if (fepAtomFlags == NULL) {
    NAMD_die("Memory allocation failed in Molecule::build_fep_params()");
  }

  double lesMassFactor = 1.0;
  if ( simParams->lesOn && simParams->lesReduceMass ) {
    lesMassFactor = 1.0 / simParams->lesFactor;
  }

  // loop through all the atoms and get the b value
  for (i = 0; i < numAtoms; i++) {
    // Get the alch flag value
    switch (bcol) {
      case 1:
        bval = (bPDB->atom(i))->xcoor();
        break;
      case 2:
        bval = (bPDB->atom(i))->ycoor();
        break;
      case 3:
        bval = (bPDB->atom(i))->zcoor();
        break;
      case 4:
        bval = (bPDB->atom(i))->occupancy();
        break;
      case 5:
        bval = (bPDB->atom(i))->temperaturefactor();
        break;
    }

    // Assign alch flag value
    if (simParams->lesOn) {
      if ( bval == (int) bval && bval > 0 ) {
        if ( bval > simParams->lesFactor ) 
          NAMD_die("LES flag must be less than or equal to lesFactor.");
        fepAtomFlags[i] = (int) bval;
      #ifdef MEM_OPT_VERSION
        Real newMass = atomMassPool[eachAtomMass[i]]*lesMassFactor;
        eachAtomMass[i] = insert_new_mass(newMass);
      #else
        atoms[i].mass *= lesMassFactor;     
      #endif
        numFepFinal++;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->alchOn) {
      if (bval == 1.0) {
        fepAtomFlags[i] = 1;
        numFepFinal++;
      } else if (bval == -1.0) {
        fepAtomFlags[i] = 2;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->pairInteractionOn) {
      if (bval == simParams->pairInteractionGroup1) {
        fepAtomFlags[i] = 1;
        numFepInitial++;
      } else if (bval == simParams->pairInteractionGroup2) {
        fepAtomFlags[i] = 2;
        numFepFinal++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->pressureProfileAtomTypes > 1) {
      fepAtomFlags[i] = (int) bval;
    } 
#ifdef OPENATOM_VERSION
    // This needs to be refactored into its build_openatom_flags fxn
    if (simParams->openatomOn) {
      if (bval != 0) {
        fepAtomFlags[i] = bval;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    }
#endif //OPENATOM_VERSION
  }

  // if PDB object was created, delete it
  if (alchfile != NULL) {
    delete bPDB;
  }
}
// End of function build_fep_flags
 
   //
   //
   //  FUNCTION delete_alch_bonded
   //
   // FB - Loop over bonds, angles, dihedrals and impropers, drop any that 
   // contain atoms of both partitions 1 and 2
   // 
   // 

#ifndef MEM_OPT_VERSION
void Molecule::delete_alch_bonded(void)  {

  // Bonds
  suspiciousAlchBonds = 0;  // these really shouldn't exist...?
  for (int i = 0; i < numBonds; i++) {
    int part1 = fepAtomFlags[bonds[i].atom1];
    int part2 = fepAtomFlags[bonds[i].atom2];
    if ((part1 == 1 || part2 == 1 ) &&
      (part1 == 2 || part2 == 2 )) {
      //CkPrintf("-----BOND ATOMS %i %i partitions %i %i \n",bonds[i].atom1, bonds[i].atom2, part1, part2);
      suspiciousAlchBonds++;
    }
  }

  // Angles
  Angle *nonalchAngles;
  nonalchAngles = new Angle[numAngles];
  int nonalchAngleCount = 0;
  alchDroppedAngles = 0;
  for (int i = 0; i < numAngles; i++) {
    int part1 = fepAtomFlags[angles[i].atom1];
    int part2 = fepAtomFlags[angles[i].atom2];
    int part3 = fepAtomFlags[angles[i].atom3];
    if ((part1 == 1 || part2 == 1 || part3 == 1) &&
      (part1 == 2 || part2 == 2 || part3 == 2)) {
      //CkPrintf("-----ANGLE ATOMS %i %i %i partitions %i %i %i\n",angles[i].atom1, angles[i].atom2, angles[i].atom3, part1, part2, part3);
      alchDroppedAngles++;
    }
    else {
      if ( angles[i].angle_type == -1 ) {
        char err_msg[128];
        sprintf(err_msg,
            "MISSING PARAMETERS FOR ANGLE %i %i %i PARTITIONS %i %i %i\n",
            angles[i].atom1+1, angles[i].atom2+1, angles[i].atom3+1,
            part1, part2, part3);
        NAMD_die(err_msg);
      }
      nonalchAngles[nonalchAngleCount++] = angles[i];
    }
  }
  numAngles = nonalchAngleCount;
  delete [] angles;
  angles = new Angle[numAngles];
  for (int i = 0; i < nonalchAngleCount; i++) {
    angles[i]=nonalchAngles[i];
  }
  delete [] nonalchAngles;


  // Dihedrals
  Dihedral *nonalchDihedrals;
  nonalchDihedrals = new Dihedral[numDihedrals];
  int nonalchDihedralCount = 0;
  alchDroppedDihedrals = 0;
  for (int i = 0; i < numDihedrals; i++) {
    int part1 = fepAtomFlags[dihedrals[i].atom1];
    int part2 = fepAtomFlags[dihedrals[i].atom2];
    int part3 = fepAtomFlags[dihedrals[i].atom3];
    int part4 = fepAtomFlags[dihedrals[i].atom4];
    if ((part1 == 1 || part2 == 1 || part3 == 1 || part4 == 1) &&
      (part1 == 2 || part2 == 2 || part3 == 2 || part4 == 2)) {
      //CkPrintf("-----i %i DIHEDRAL ATOMS %i %i %i %i partitions %i %i %i %i\n",i,dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4, part1, part2, part3,part4);
      alchDroppedDihedrals++;
    }
    else {
      if ( dihedrals[i].dihedral_type == -1 ) {
        char err_msg[128];
        sprintf(err_msg,
        "MISSING PARAMETERS FOR DIHEDRAL %i %i %i %i PARTITIONS %i %i %i %i\n",
            dihedrals[i].atom1+1, dihedrals[i].atom2+1,
            dihedrals[i].atom3+1, dihedrals[i].atom4+1,
            part1, part2, part3, part4);
        NAMD_die(err_msg);
      }
      nonalchDihedrals[nonalchDihedralCount++] = dihedrals[i];
    }
  }
  numDihedrals = nonalchDihedralCount;
  delete [] dihedrals;
  dihedrals = new Dihedral[numDihedrals];
  for (int i = 0; i < numDihedrals; i++) {
    dihedrals[i]=nonalchDihedrals[i];
  }
  delete [] nonalchDihedrals;

  // Impropers
  Improper *nonalchImpropers;
  nonalchImpropers = new Improper[numImpropers];
  int nonalchImproperCount = 0;
  alchDroppedImpropers = 0;
  for (int i = 0; i < numImpropers; i++) {
    int part1 = fepAtomFlags[impropers[i].atom1];
    int part2 = fepAtomFlags[impropers[i].atom2];
    int part3 = fepAtomFlags[impropers[i].atom3];
    int part4 = fepAtomFlags[impropers[i].atom4];
    if ((part1 == 1 || part2 == 1 || part3 == 1 || part4 == 1) &&
      (part1 == 2 || part2 == 2 || part3 == 2 || part4 == 2)) {
      //CkPrintf("-----i %i IMPROPER ATOMS %i %i %i %i partitions %i %i %i %i\n",i,impropers[i].atom1, impropers[i].atom2, impropers[i].atom3, impropers[i].atom4, part1, part2, part3,part4);
      alchDroppedImpropers++;
    }
    else {
      nonalchImpropers[nonalchImproperCount++] = impropers[i];
    }
  }
  numImpropers = nonalchImproperCount;
  delete [] impropers;
  impropers = new Improper[numImpropers];
  for (int i = 0; i < numImpropers; i++) {
    impropers[i]=nonalchImpropers[i];
  }
  delete [] nonalchImpropers;
  
} // end delete_alch_bonded
#endif  

//fepe



void Molecule::build_exPressure_atoms(StringList *fixedfile, 
   StringList *fixedcol, PDB *initial_pdb, char *cwd) {
       
  PDB *bPDB;      //  Pointer to PDB object to use
  int bcol = 4;      //  Column that data is in
  Real bval = 0;      //  b value from PDB file
  int i;      //  Loop counter
  char filename[129];    //  Filename

  //  Get the PDB object that contains the b values.  If
  //  the user gave another file name, use it.  Otherwise, just use
  //  the PDB file that has the initial coordinates.
  if (fixedfile == NULL) {
    if ( ! initial_pdb ) NAMD_die("Initial PDB file unavailable, excludeFromPressureFile required.");
    bPDB = initial_pdb;
  } else {
    if (fixedfile->next != NULL) {
      NAMD_die("Multiple definitions of excluded pressure atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (fixedfile->data[0] == '/') ) {
         strcpy(filename, fixedfile->data);
    } else {
         strcpy(filename, cwd);
         strcat(filename, fixedfile->data);
    }
    bPDB = new PDB(filename);
    if ( bPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_exPressure_atoms");
    }

    if (bPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in excludedPressure atoms PDB doesn't match coordinate PDB");
    }
  }

  //  Get the column that the b vaules are in.  It
  //  can be in any of the 5 floating point fields in the PDB, according
  //  to what the user wants.  The allowable fields are X, Y, Z, O, or
  //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
  //  The default is the 4th field, which is the occupancy
  if (fixedcol == NULL) {
    bcol = 4;
  } else {
    if (fixedcol->next != NULL) {
      NAMD_die("Multiple definitions of excludedPressure atoms column in config file");
    }

    if (strcasecmp(fixedcol->data, "X") == 0) {
       bcol=1;
    } else if (strcasecmp(fixedcol->data, "Y") == 0) {
       bcol=2;
    } else if (strcasecmp(fixedcol->data, "Z") == 0) {
       bcol=3;
    } else if (strcasecmp(fixedcol->data, "O") == 0) {
       bcol=4;
    } else if (strcasecmp(fixedcol->data, "B") == 0) {
       bcol=5;
    } else {
       NAMD_die("excludedPressureFileCol must have value of X, Y, Z, O, or B");
    }
  }

  //  Allocate the array to hold all the data
  exPressureAtomFlags = new int32[numAtoms];

  if (exPressureAtomFlags == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_fixed_atoms()");
  }

  numExPressureAtoms = 0;

  //  Loop through all the atoms and get the b value
  for (i=0; i<numAtoms; i++) {
    //  Get the k value based on where we were told to find it
    switch (bcol) {
       case 1: bval = (bPDB->atom(i))->xcoor(); break;
       case 2: bval = (bPDB->atom(i))->ycoor(); break;
       case 3: bval = (bPDB->atom(i))->zcoor(); break;
       case 4: bval = (bPDB->atom(i))->occupancy(); break;
       case 5: bval = (bPDB->atom(i))->temperaturefactor(); break;
    }

    //  Assign the b value
    if ( bval != 0 ) {
      exPressureAtomFlags[i] = 1;
      numExPressureAtoms++;
    } else {
      exPressureAtomFlags[i] = 0;
    }
  }
  if (fixedfile != NULL) 
    delete bPDB;

  iout << iINFO << "Got " << numExPressureAtoms << " excluded pressure atoms." 
       << endi;
}


    Bool Molecule::is_lp(int anum) {
      return ((atoms[anum].status & LonepairAtom) != 0);
    }

    Bool Molecule::is_drude(int anum) {
      return ((atoms[anum].status & DrudeAtom) != 0);
    }

    Bool Molecule::is_hydrogen(int anum)
    {
  return ((atoms[anum].status & HydrogenAtom) != 0);
    }

    Bool Molecule::is_oxygen(int anum)
    {
  return ((atoms[anum].status & OxygenAtom) != 0);
    }

    Bool Molecule::is_hydrogenGroupParent(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].isGP);
    }

    Bool Molecule::is_water(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].waterVal == 2);
    }

    int Molecule::get_groupSize(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].atomsInGroup);
    }

    int Molecule::get_mother_atom(int anum)
    {
  // for efficiency reasons, we are not checking if anum is already 
  // hydrogen or not. This function must be called for hydrogens only;
  return atoms[anum].partner;
    }

void Molecule::reloadCharges(float charge[], int n){
  if ( n != numAtoms )
    NAMD_die("Incorrect number of atoms in Molecule::reloadCharges().");

#ifdef MEM_OPT_VERSION
    delete [] atomChargePool;
    vector<Real> tmpCharges;
    for(int i=0; i<numAtoms; i++){
        int foundIdx=-1;
        //naive searching, better to be binary searching but requiring 
        //inserting charges in increasing/decreasing order
        for(int j=0; j<tmpCharges.size();j++){
        	if(tmpCharges[j] == charge[i]){
        	    foundIdx = j;
        	    break;
        	}
        }
        if(foundIdx==-1){
        	tmpCharges.push_back(charge[i]);
        	foundIdx = tmpCharges.size()-1;
        }
        eachAtomCharge[i] = (Index)foundIdx;
    }
    chargePoolSize = tmpCharges.size();
    atomChargePool = new Real[chargePoolSize];
    for(int i=0; i<chargePoolSize; i++)
        atomChargePool[i] = tmpCharges[i];
#else
  for( int i=0; i<n; ++i ) atoms[i].charge = charge[i];
#endif
}

#ifndef MEM_OPT_VERSION	
// go through the molecular structure, analyze the status of each atom,
// and save the data in the Atom structures stored for each atom.  This
// could be built up incrementally while the molecule is being read in,
// but doing it all in one shot allows us to just send the basic info
// over the network and have each node calculate the rest of the data on
// it's own.
void Molecule::build_atom_status(void) {
  register int i;
  int a1, a2, a3;
  int numDrudeWaters = 0;

  // if any atoms have a mass of zero set to 0.001 and warn user
  int numZeroMassAtoms = 0;
  for (i=0; i < numAtoms; i++) {
    if ( atoms[i].mass <= 0. ) {
      if (simParams->watmodel == WAT_TIP4 ||
          simParams->lonepairs) {
        ++numLonepairs;
      } else {
        atoms[i].mass = 0.001;
        ++numZeroMassAtoms;
      }
    }
    else if (atoms[i].mass < 1.) {
      ++numDrudeAtoms;
    }
  }
  // DRUDE: verify number of LPs
  if (simParams->lonepairs && numLonepairs != numLphosts) {
    NAMD_die("must have same number of LP hosts as lone pairs");
  }
  // DRUDE
  if ( ! CkMyPe() ) {
    if (simParams->watmodel == WAT_TIP4 || simParams->lonepairs) {
      iout << iWARN << "CORRECTION OF ZERO MASS ATOMS TURNED OFF "
        "BECAUSE LONE PAIRS ARE USED\n" << endi;
    } else if ( numZeroMassAtoms ) {
      iout << iWARN << "FOUND " << numZeroMassAtoms <<
        " ATOMS WITH ZERO OR NEGATIVE MASSES!  CHANGED TO 0.001\n" << endi;
    }
  }
  // initialize information for each atom (note that the status has
  // already been initialized during the read/receive phase)
  hydrogenGroup.resize(numAtoms);
  HydrogenGroupID *hg = hydrogenGroup.begin();
  for (i=0; i < numAtoms; i++) {
    atoms[i].partner = (-1);
    hg[i].atomID = i;  // currently unsorted
    hg[i].atomsInGroup = 1;  // currently only 1 in group
    hg[i].isGP = 1;  // assume it is a group parent
    hg[i].GPID = i;  // assume it is a group parent
    hg[i].waterVal = 0;  // for group sorting
  }

  // deal with H-H bonds in a sane manner
  // this information will be rewritten later if bonded elsewhere
  int hhbondcount = 0;
  for (i=0; i < numRealBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1) && is_hydrogen(a2)) {
      ++hhbondcount;
      // make H atoms point at each other for now
      atoms[a1].partner = a2;
      atoms[a2].partner = a1;
      hg[a1].atomsInGroup++;
      hg[a1].GPID = a2;
      hg[a2].atomsInGroup++;
      hg[a2].GPID = a1;
    }
  }

  if ( hhbondcount && ! CkMyPe() ) {
    iout << iWARN << "Found " << hhbondcount << " H-H bonds.\n" << endi;
  }

  // find which atom each hydrogen is bound to
  // also determine number of atoms in each group
  for (i=0; i < numRealBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1)) {
      if (is_hydrogen(a2)) continue;
      atoms[a1].partner = a2;
      hg[a2].atomsInGroup++;
      hg[a1].atomsInGroup = 0;
      hg[a1].GPID = a2;
      hg[a1].isGP = 0;
      // check for waters (put them in their own groups: OH or OHH)
      if (is_oxygen(a2))  hg[a2].waterVal++;
    }
    if (is_hydrogen(a2)) {
      atoms[a2].partner = a1;
      hg[a1].atomsInGroup++;
      hg[a2].atomsInGroup = 0;
      hg[a2].GPID = a1;
      hg[a2].isGP = 0;
      // check for waters (put them in their own groups: OH or OHH)
      if (is_oxygen(a1))  hg[a1].waterVal++;
    }

    // If we have TIP4P water, check for lone pairs
    if (simParams->watmodel == WAT_TIP4) {
      if (is_lp(a1)) {
        atoms[a1].partner = a2;
        hg[a2].atomsInGroup++;
        hg[a1].atomsInGroup = 0;
        hg[a1].GPID = a2;
        hg[a1].isGP = 0;
      }
      if (is_lp(a2)) {
        atoms[a2].partner = a1;
        hg[a1].atomsInGroup++;
        hg[a2].atomsInGroup = 0;
        hg[a2].GPID = a1;
        hg[a2].isGP = 0;
      }
    }
    // SWM4 water has lone pair and Drude particles
    else if ( /* simParams->watmodel == WAT_SWM4 */ simParams->lonepairs) {
      if (is_lp(a1) || is_drude(a1)) {
        if (is_hydrogen(a2) || is_lp(a2) || is_drude(a2)) {
          char msg[256];
          sprintf(msg, "%s particle %d is bonded to non-parent atom %d",
              (is_lp(a1) ? "Lone pair" : "Drude"), a1+1, a2+1);
          NAMD_die(msg);
        }
        atoms[a1].partner = a2;
        hg[a2].atomsInGroup++;
        hg[a1].atomsInGroup = 0;
        hg[a1].GPID = a2;
        hg[a1].isGP = 0;
      }
      else if (is_lp(a2) || is_drude(a2)) {
        if (is_hydrogen(a1) || is_lp(a1) || is_drude(a1)) {
          char msg[256];
          sprintf(msg, "%s particle %d is bonded to non-parent atom %d",
              (is_lp(a2) ? "Lone pair" : "Drude"), a2+1, a1+1);
          NAMD_die(msg);
        }
        atoms[a2].partner = a1;
        hg[a1].atomsInGroup++;
        hg[a2].atomsInGroup = 0;
        hg[a2].GPID = a1;
        hg[a2].isGP = 0;
      }
    }

  }

  // check up on our H-H bonds and general sanity check
  int hGPcount = 0;
  for(i=0; i<numAtoms; i++) {
    if ( ! hg[hg[i].GPID].isGP ) {
      char msg[256];
      sprintf(msg, "child atom %d bonded only to child H atoms",i+1);
      NAMD_die(msg);
    }
    if ( hg[i].isGP && is_hydrogen(i) ) {
      if ( hg[i].GPID == i ) continue;  // atomic hydrogen ion
      ++hGPcount;  // molecular hydrogen
      if ( is_hydrogen(hg[i].GPID) && hg[hg[i].GPID].GPID != i ) {
        char msg[256];
        sprintf(msg, "H atom %d bonded only to child H atoms",i+1);
        NAMD_die(msg);
      }
      hg[hg[i].GPID].atomsInGroup = 0;
      hg[hg[i].GPID].isGP = 0;
      hg[i].GPID = i;
      if ( hg[i].atomsInGroup != 2 ) {
        char msg[256];
        sprintf(msg, "H atom %d bonded to multiple H atoms",i+1);
        NAMD_die(msg);
      }
    }
  }
  if ( hGPcount && ! CkMyPe() ) {
    iout << iWARN << "Found " << hGPcount << " H-H molecules.\n" << endi;
  }

  // copy hydrogen groups to migration groups
  for (i=0; i<numAtoms; ++i) {
    if ( hg[i].isGP ) hg[i].GPID = i;  // group parent is its own parent
    else hg[i].waterVal = hg[hg[i].GPID].waterVal;  // copy to children
    hg[i].MPID = hg[i].GPID;
  }

  // determine migration groups based on lone pair hosts
  for (i=0; i<numLphosts; ++i) {
    int a1 = lphosts[i].atom1;
    int a2 = lphosts[i].atom2;
    int a3 = lphosts[i].atom3;
    int a4 = lphosts[i].atom4;
    int m1 = hg[a1].MPID;
    while ( hg[m1].MPID != m1 ) m1 = hg[m1].MPID;
    int m2 = hg[a2].MPID;
    while ( hg[m2].MPID != m2 ) m2 = hg[m2].MPID;
    int m3 = hg[a3].MPID;
    while ( hg[m3].MPID != m3 ) m3 = hg[m3].MPID;
    int m4 = hg[a4].MPID;
    while ( hg[m4].MPID != m4 ) m4 = hg[m4].MPID;
    int mp = m1;
    if ( m2 < mp ) mp = m2;
    if ( m3 < mp ) mp = m3;
    if ( m4 < mp ) mp = m4;
    hg[m1].MPID = mp;
    hg[m2].MPID = mp;
    hg[m3].MPID = mp;
    hg[m4].MPID = mp;
  }
  while ( 1 ) {
    int allok = 1;
    for (i=0; i<numAtoms; ++i) {
      int mp = hg[i].MPID;
      if ( hg[mp].MPID != mp ) {
        allok = 0;
        hg[i].MPID = hg[mp].MPID;
      }
    }
    if ( allok ) break;
  }
  for (i=0; i<numAtoms; ++i) {
    hg[i].isMP = ( hg[i].MPID == i );
    hg[i].atomsInMigrationGroup = 0;
  }
  for (i=0; i<numAtoms; ++i) {
    hg[hg[i].MPID].atomsInMigrationGroup++;
  }

  if ( simParams->splitPatch != SPLIT_PATCH_HYDROGEN ) {
    // every atom its own group
    for (i=0; i<numAtoms; i++) {
      hg[i].isGP = 1;
      hg[i].isMP = 1;
      hg[i].atomsInGroup = 1;
      hg[i].atomsInMigrationGroup = 1;
      hg[i].GPID = i;
      hg[i].MPID = i;
    }
  }

  // count number of groups
  numHydrogenGroups = 0;
  maxHydrogenGroupSize = 0;
  numMigrationGroups = 0;
  maxMigrationGroupSize = 0;
  for(i=0; i<numAtoms; i++)
  {
    if (hg[i].isMP) {
      ++numMigrationGroups;
      int mgs = hg[i].atomsInMigrationGroup;
      if ( mgs > maxMigrationGroupSize ) maxMigrationGroupSize = mgs;
    }
    if (hg[i].isGP) {
      ++numHydrogenGroups;
      int hgs = hg[i].atomsInGroup;
      if ( hgs > maxHydrogenGroupSize ) maxHydrogenGroupSize = hgs;
    }
  }

  hydrogenGroup.sort();

  // sanity checking
  int parentid = -1;
  int hgs = 0;
  for(i=0; i<numAtoms; ++i, --hgs) {
    if ( ! hgs ) {  // expect group parent
      if ( hg[i].isGP ) {
        hgs = hg[i].atomsInGroup;
        parentid = hg[i].atomID;
      } else {
        char buff[512];
        sprintf(buff, "Atom %d has bad hydrogen group size.  "
            "Check for duplicate bonds.", parentid+1);
        NAMD_die(buff);
      }
    } else {  // don't expect group parent
      if ( hg[i].isGP ) {
        char buff[512];
        sprintf(buff, "Atom %d has bad hydrogen group size.  "
            "Check for duplicate bonds.", parentid+1);
        NAMD_die(buff);
      }
    }
  }

  parentid = -1;
  int mgs = 0;
  for(i=0; i<numAtoms; ++i, --mgs) {
    if ( ! mgs ) {  // expect group parent
      if ( hg[i].isMP ) {
        mgs = hg[i].atomsInMigrationGroup;
        parentid = hg[i].atomID;
      } else {
        char buff[512];
        sprintf(buff, "Atom %d has bad migration group size.", parentid+1);
        NAMD_die(buff);
      }
    } else {  // don't expect group parent
      if ( hg[i].isMP ) {
        char buff[512];
        sprintf(buff, "Atom %d has bad migration group size.", parentid+1);
        NAMD_die(buff);
      }
    }
  }


  // finally, add the indexing from atoms[] to hydrogenGroup[]
  for(i=0; i<numAtoms; i++) {
    atoms[hydrogenGroup[i].atomID].hydrogenList = i;
  }

  // check ordering of Drude particles and water
  // count number of Drude waters
  if (simParams->watmodel == WAT_SWM4) {
    for (i = 0;  i < numAtoms;  i++) {
      if (is_water(hg[i].atomID) && hg[i].isGP) {
        if (i > numAtoms-5
            || ! is_drude(hg[i+1].atomID)
            || ! is_lp(hg[i+2].atomID)
            || ! is_hydrogen(hg[i+3].atomID)
            || ! is_hydrogen(hg[i+4].atomID) ) {
          char msg[256];
          sprintf(msg, "Drude water molecule from HydrogenGroup i=%d "
              "starting at atom %d is not sorted\n", i, hg[i].atomID+1);
          NAMD_die(msg);
        }
        numDrudeWaters++;
        i += 4;  // +1 from loop
        continue;
      } // if water
      else if (is_drude(hg[i].atomID)) {
        if (i < 1 || hg[i-1].atomID != hg[i].GPID) {
          char msg[256];
          sprintf(msg, "Drude particle from HydrogenGroup i=%d must "
              "immediately follow its parent atom %d\n", i, hg[i].GPID+1);
          NAMD_die(msg);
        }
      } // else if Drude
#if 0
      else if (is_lp(hg[i].atomID)) {
        char msg[256];
        sprintf(msg, "Drude lonepair from HydrogenGroup i=%d "
            "at particle %d is NOT from water - unsupported\n",
            i, hg[i].atomID+1);
        NAMD_die(msg);
      }
#endif
    } // for numAtoms
  } // if SWM4

  #if 0 
  // debugging code for showing sorted atoms
  if(CkMyPe()==1) {  
  for(i=0; i<numAtoms; i++)
    iout << i << " atomID=" << hydrogenGroup[i].atomID
   << " isGP=" << hydrogenGroup[i].isGP
   << " parent=" << hydrogenGroup[i].GPID
   << " #" << hydrogenGroup[i].atomsInGroup
   << " waterVal=" << hydrogenGroup[i].waterVal
   << " partner=" << atoms[i].partner
   << " hydrogenList=" << atoms[i].hydrogenList
   << "\n" << endi;
  }
  #endif

  // now deal with rigidBonds
  if ( simParams->rigidBonds != RIGID_NONE || simParams->mollyOn ) {
    // temporary variables for use by 4+ site water models
    Real r_oh = -1.0;
    Real r_hh = -1.0;

    delete [] rigidBondLengths;
    rigidBondLengths = new Real[numAtoms];
    if ( ! rigidBondLengths ) {
      NAMD_die("Memory allocation failed in Molecule::build_atom_status()\n");
    }
    for (i=0; i<numAtoms; ++i) rigidBondLengths[i] = 0;
    int mode = simParams->rigidBonds;
    if ( simParams->mollyOn ) mode = RIGID_ALL;

    // add H-mother lengths or 0 if not constrained
    for (i=0; i < numRealBonds; i++) {
      a1 = bonds[i].atom1;
      a2 = bonds[i].atom2;
      Real dum, x0;
      params->get_bond_params(&dum,&x0,bonds[i].bond_type);
      if (is_hydrogen(a2)) { int tmp = a1;  a1 = a2;  a2 = tmp; } // swap
      if (is_hydrogen(a1)) {
        if ( is_hydrogen(a2) ) {  // H-H
          if ( ! is_water(a2) ) {  // H-H but not water
	    rigidBondLengths[a1] = ( mode == RIGID_ALL ? x0 : 0. );
	    rigidBondLengths[a2] = ( mode == RIGID_ALL ? x0 : 0. );
          }
        } else if ( is_water(a2) || mode == RIGID_ALL ) {
	  rigidBondLengths[a1] = x0;
    if (is_water(a2)) r_oh = rigidBondLengths[a1];
	} else {
	  rigidBondLengths[a1] = 0.;
        }
      }
      // Handle lone pairs if they're allowed
      if (simParams->watmodel == WAT_TIP4) {
        if (is_lp(a2)) { int tmp = a1;  a1 = a2;  a2 = tmp; } // swap
        if (is_lp(a1)) {
          if (! is_water(a2) ) {
            // Currently, lonepairs are only allowed on waters,
            // although this may change in the future
            char err_msg[128];
            sprintf(err_msg, "ILLEGAL LONE PAIR AT INDEX %i\n"
                "LONE PAIRS ARE CURRENTLY ALLOWED ONLY ON WATER MOLECULES\n",
                a1);
            NAMD_die(err_msg);
          } else {
            rigidBondLengths[a1] = x0;
            r_om = x0;
          }
        }
      }
      // Handle SWM4 lone pairs
      // (Drude bonds remain flexible)
      if (simParams->watmodel == WAT_SWM4) {
        if (is_lp(a2)) {
          int tmp = a1;  a1 = a2;  a2 = tmp;  // swap
        }
        if (is_lp(a1)) {
          if (is_water(a2)) {
            // do not count bonds to LPs as rigid, do not set rigidBondLengths[]
            r_om = x0;  // for faster position update routine for LP on water
          }
          else if ( ! simParams->drudeOn) {
            // if not using Drude model, lone pairs allowed only on water
            char msg[128];
            sprintf(msg, "ILLEGAL LONE PAIR AT INDEX %d\n"
                "LONE PAIRS ARE CURRENTLY ALLOWED ONLY ON WATER MOLECULES\n",
                a1+1);
            NAMD_die(msg);
          }
        }
      }
    }

    // zero out H-H lengths - water handled below
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) rigidBondLengths[h_i->atomID] = 0.;
    }

    // fill in H-H lengths for water by searching angles - yuck
    for (i=0; i < numAngles; i++) {
      a2 = angles[i].atom2;
      if ( ! is_water(a2) ) continue;
      if ( ! is_oxygen(a2) ) continue;
      a1 = angles[i].atom1;
      if ( ! is_hydrogen(a1) ) continue;
      a3 = angles[i].atom3;
      if ( ! is_hydrogen(a3) ) continue;
      if (is_lp(a2) || is_lp(a1) || is_lp(a3) ||
          is_drude(a2) || is_drude(a1) || is_drude(a3)) continue;
      if ( rigidBondLengths[a1] != rigidBondLengths[a3] ) {
        if (rigidBondLengths[a1] >0.3 && rigidBondLengths[a3] >0.3) {
          printf("length1: %f length2: %f\n", rigidBondLengths[a1], rigidBondLengths[a3]);

          NAMD_die("Asymmetric water molecule found???  This can't be right.\n");
        }
      }
      Real dum, t0;
      params->get_angle_params(&dum,&t0,&dum,&dum,angles[i].angle_type);
      rigidBondLengths[a2] = 2. * rigidBondLengths[a1] * sin(0.5*t0);
      r_hh = rigidBondLengths[a2];
    }

    // fill in H-H lengths for waters that are missing angles
    int numBondWaters = 0;
    int numFailedWaters = 0;

    for (i=0; i < numRealBonds; i++) {
      a1 = bonds[i].atom1;
      a2 = bonds[i].atom2;
      if ( ! is_hydrogen(a1) ) continue;
      if ( ! is_hydrogen(a2) ) continue;
      int ma1 = get_mother_atom(a1);
      int ma2 = get_mother_atom(a2);
      if ( ma1 != ma2 ) continue;
      if ( ! is_water(ma1) ) continue;
      if ( rigidBondLengths[ma1] != 0. ) continue;
      Real dum, x0;
      params->get_bond_params(&dum,&x0,bonds[i].bond_type);
      rigidBondLengths[ma1] = x0;
    }
    
    // We now should be able to set the parameters needed for water lonepairs
    // make sure that there is water in the system
    if ( (simParams->watmodel == WAT_TIP4 && numLonepairs > 0)
        || (simParams->watmodel == WAT_SWM4 && numDrudeWaters > 0)) {
      if (r_oh < 0.0 || r_hh < 0.0) {
        //printf("ERROR: r_oh %f / r_hh %f\n", r_oh, r_hh);
        NAMD_die("Failed to find water bond lengths\n");
      } 
      r_ohc = sqrt(r_oh * r_oh - 0.25 * r_hh * r_hh);
      //printf("final r_om and r_ohc are %f and %f\n", r_om, r_ohc);
    }

    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP && is_water(h_i->atomID) &&
                     rigidBondLengths[h_i->atomID] == 0. ) {
        if ( h_i + 1 == h_e || h_i + 2 == h_e ||
             h_i[1].isGP || h_i[2].isGP || h_i->atomsInGroup != 3 ) {
          NAMD_die("Abnormal water detected.");
        }
        if ( CkNumNodes() > 1 ) {
          NAMD_die("Unable to determine H-H distance for rigid water because structure has neither H-O-H angle nor H-H bond.");
        }
        Bond btmp;
        btmp.atom1 = h_i[1].atomID;
        char atom1name[11];
	strcpy(atom1name,get_atomtype(btmp.atom1));
        btmp.atom2 = h_i[2].atomID;
        char atom2name[11];
	strcpy(atom2name,get_atomtype(btmp.atom2));
        params->assign_bond_index(atom1name,atom2name,&btmp);
        Real k, x0;
	x0 = 0.;
        params->get_bond_params(&k,&x0,btmp.bond_type);
	if ( x0 > 0. ) {
          rigidBondLengths[h_i->atomID] = x0;
	  numBondWaters++;
        } else {
	  numFailedWaters++;
        }
      }
    }
    if ( numBondWaters + numFailedWaters ) {
      iout << iWARN << "Missing angles for " <<
	      ( numBondWaters + numFailedWaters ) << " waters.\n" << endi;
    }
    if ( numBondWaters ) {
      iout << iWARN << "Obtained H-H distance from bond parameters for " <<
	      numBondWaters << " waters.\n" << endi;
      iout << iWARN << "This would not be possible in a multi-process run.\n" << endi;
    }
    if ( numFailedWaters ) {
      iout << iERROR << "Failed to obtain H-H distance from angles or bonds for " <<
	      numFailedWaters << " waters.\n" << endi;
    }

    // in case both molly and rigidBonds are in use make lengths which
    // are molly-only negative and leave lengths which are both alone
    if ( simParams->mollyOn ) {
      mode = simParams->rigidBonds;
      if ( mode == RIGID_NONE ) {
	for (i=0; i<numAtoms; ++i) rigidBondLengths[i] *= -1;
      } else if ( mode == RIGID_WATER ) {
	for (i=0; i<numAtoms; ++i) {
	  if ( ! is_water(i) ) rigidBondLengths[i] *= -1;
	}
      }
    }

    numRigidBonds = 0;
    for (i=0; i<numAtoms; ++i) {
      if ( rigidBondLengths[i] > 0. ) ++numRigidBonds;
    }

  }
  }

/****************************************************************************/
/*  FUNCTION compute_LJcorrection                                           */
/*                                                                          */
/*  Compute the energy and virial tail corrections to the Lennard-Jones     */
/*  potential. The approximation used for heterogenous systems is to compute*/
/*  the average pairwise parameters as in Ref 2.  Additional terms are also */
/*  added in the case of potential or force switching.                      */
/*                                                                          */
/*  REFERENCES                                                              */
/*   1) Allen and Tildesley, Computer Simulation of Liquids, 1991           */
/*   2) Shirts, et al. J Phys Chem B. 2007 111:13052                        */
/****************************************************************************/
void Molecule::compute_LJcorrection() {
  // First, calculate the average A and B coefficients. For TI/FEP, decompose
  // by alchemical group (1 or 2).
  BigReal LJAvgA, LJAvgB, LJAvgA1, LJAvgB1, LJAvgA2, LJAvgB2;

  /*This a shortcut to summing over all atoms since it is faster to count how 
    many atoms are of each LJ type.

    NB: In practice it is easier to double count pairs. That is, we get N*(N-1)
        pairs instead of N*(N-1)/2 and the 2 cancels in the numerator and
        denominator. This affects later corrections to the sums!
  */
  int LJtypecount = params->get_num_vdw_params();
  Real A, B, A14, B14;
  Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
  Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;
  Real *ATable = new Real[LJtypecount*LJtypecount];
  Real *BTable = new Real[LJtypecount*LJtypecount];
  int useGeom = simParams->vdwGeometricSigma;
  // copied from LJTable.C
  for (int i = 0; i < LJtypecount; i++) {
    for (int j = 0; j < LJtypecount; j++) {
      if (params->get_vdw_pair_params(i,j, &A, &B, &A14, &B14)) {
        ATable[i*LJtypecount + j] = A;
        BTable[i*LJtypecount + j] = B;
      }
      else {
        params->get_vdw_params(&sigma_i,&epsilon_i,&sigma_i14,&epsilon_i14,i);
        params->get_vdw_params(&sigma_j,&epsilon_j,&sigma_j14,&epsilon_j14,j);
        BigReal sigma_ij =
          useGeom ? sqrt(sigma_i*sigma_j) : 0.5*(sigma_i+sigma_j);
        BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);
        sigma_ij *= sigma_ij*sigma_ij;
        sigma_ij *= sigma_ij;

        ATable[i*LJtypecount + j] = 4.0*sigma_ij*epsilon_ij*sigma_ij;
        BTable[i*LJtypecount + j] = 4.0*sigma_ij*epsilon_ij;
      }
    }
  }

  int *numAtomsByLjType = new int[LJtypecount];
  for (int i=0; i < LJtypecount; i++) {numAtomsByLjType[i]=0;}
  for (int i=0; i < numAtoms; i++) {numAtomsByLjType[atoms[i].vdw_type]++;}

  BigReal sumOfAs = 0;
  BigReal sumOfBs = 0;
  BigReal count = 0; // needed to avoid overflow
  BigReal npairs;
  for (int i=0; i < LJtypecount; i++) {
    for (int j=0; j < LJtypecount; j++) {
      A = ATable[i*LJtypecount + j];
      B = BTable[i*LJtypecount + j];
      if (!A && !B) continue; // don't count zeroed interactions
      npairs = (numAtomsByLjType[i] - int(i==j))*BigReal(numAtomsByLjType[j]);
      sumOfAs += npairs*A;
      sumOfBs += npairs*B;
      count += npairs;
    }
  }
  delete [] numAtomsByLjType;
  delete [] ATable;
  delete [] BTable;

  /*If alchemical interactions exist, account for interactions that disappear
    at the endpoints. Since alchemical transformations are path independent,
    the intermediate values can be treated fairly arbitrarily.  IMO, the
    easiest thing to do is have the lambda dependent correction be a linear
    interpolation of the endpoint corrections:

    Ecorr(lambda) = lambda*Ecorr(1) + (1-lambda)*Ecorr(0)

    This makes the virial and alchemical derivative very simple also. One
    alternative would be to count "fractional interactions," but that makes
    TI derivatives a bit harder and for no obvious gain.
  */
  if (simParams->alchOn) {
    BigReal sumOfAs1 = sumOfAs;
    BigReal sumOfAs2 = sumOfAs;
    BigReal sumOfBs1 = sumOfBs;
    BigReal sumOfBs2 = sumOfBs;
    BigReal count1 = count;
    BigReal count2 = count;
    int alch_counter = 0;
    for (int i=0; i < numAtoms; ++i) {
      int alchFlagi = (get_fep_type(i) == 2 ? -1 : get_fep_type(i));
      for (int j=i+1; j < numAtoms; ++j) {
        int alchFlagj = (get_fep_type(j) == 2 ? -1 : get_fep_type(j));
        int alchFlagSum = alchFlagi + alchFlagj;

        // Ignore completely non-alchemical pairs.
        if (alchFlagi == 0 && alchFlagj == 0) continue;

        if (params->get_vdw_pair_params(atoms[i].vdw_type, atoms[j].vdw_type,
                                        &A, &B, &A14, &B14)) {
        }
        else {
          params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
                                 &epsilon_i14, atoms[i].vdw_type);
          params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14,
                                 &epsilon_j14, atoms[j].vdw_type);
          BigReal sigma_ij =
            useGeom ? sqrt(sigma_i*sigma_j) : 0.5*(sigma_i+sigma_j);
          BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);

          sigma_ij *= sigma_ij*sigma_ij;
          sigma_ij *= sigma_ij;
          A = 4.0*sigma_ij*epsilon_ij*sigma_ij;
          B = 4.0*sigma_ij*epsilon_ij;
        }
        if (!A && !B) continue; // don't count zeroed interactions

        if ( alchFlagSum > 0 ){ // in group 1, remove from group 2
          sumOfAs2 -= 2*A;
          sumOfBs2 -= 2*B;
          count2 -= 2;
        }
        else if ( alchFlagSum < 0 ){ // in group 2, remove from group 1
          sumOfAs1 -= 2*A;
          sumOfBs1 -= 2*B;
          count1 -= 2;
        }
        else{ // between groups 1 and 2, remove entirely (don't exist!)
          sumOfAs1 -= 2*A;
          sumOfBs1 -= 2*B;
          count1 -= 2;
          sumOfAs2 -= 2*A;
          sumOfBs2 -= 2*B;
          count2 -= 2;
        }
      }
      // This should save _tons_ of time, since the alchemical atoms are almost
      // always at the top of the pdb file.
      if ( alchFlagi == 1 || alchFlagi == -1 ) alch_counter++;
      if ( alch_counter == (numFepInitial + numFepFinal) ) break;
    }
    LJAvgA1 = sumOfAs1 / count1;
    LJAvgB1 = sumOfBs1 / count1;
    LJAvgA2 = sumOfAs2 / count2;
    LJAvgB2 = sumOfBs2 / count2;
    if ( ! CkMyPe() ) {
      iout << iINFO << "LONG-RANGE LJ: APPLYING ANALYTICAL CORRECTIONS TO "
           << "ENERGY AND PRESSURE\n" << endi;
      iout << iINFO << "LONG-RANGE LJ: AVERAGE A0 AND B0 COEFFICIENTS "
           << LJAvgA2 << " AND " << LJAvgB2 << "\n" << endi;
      iout << iINFO << "LONG-RANGE LJ: AVERAGE A1 AND B1 COEFFICIENTS "
           << LJAvgA1 << " AND " << LJAvgB1 << "\n" << endi;
    }

    // Pre-scale by the atom counts, as they differ from when alchemy is off.
    LJAvgA1 *= BigReal(numAtoms - numFepInitial)*(numAtoms - numFepInitial);
    LJAvgB1 *= BigReal(numAtoms - numFepInitial)*(numAtoms - numFepInitial);
    LJAvgA2 *= BigReal(numAtoms - numFepFinal)*(numAtoms - numFepFinal);
    LJAvgB2 *= BigReal(numAtoms - numFepFinal)*(numAtoms - numFepFinal);

    LJAvgA = LJAvgA2;
    LJAvgB = LJAvgB2;
  }
  else{
    LJAvgA1 = LJAvgB1 = LJAvgA2 = LJAvgB2 = 0;
    LJAvgA = sumOfAs / count;
    LJAvgB = sumOfBs / count;

    if ( ! CkMyPe() ) {
      iout << iINFO << "LONG-RANGE LJ: APPLYING ANALYTICAL CORRECTIONS TO "
           << "ENERGY AND PRESSURE\n" << endi;
      iout << iINFO << "LONG-RANGE LJ: AVERAGE A AND B COEFFICIENTS "
           << LJAvgA << " AND " << LJAvgB << "\n" << endi;
    }

    // Pre-scale by the atom counts, as they differ from when alchemy is on.
    LJAvgA *= BigReal(numAtoms)*numAtoms;
    LJAvgB *= BigReal(numAtoms)*numAtoms;
  }

  BigReal rcut = simParams->cutoff;
  BigReal rcut2 = rcut*rcut;
  BigReal rcut3 = rcut*rcut2;
  BigReal rcut4 = rcut2*rcut2;
  BigReal rcut5 = rcut2*rcut3;
  BigReal rcut6 = rcut3*rcut3;
  BigReal rcut9 = rcut5*rcut4;
  BigReal rswitch = simParams->switchingDist;
  BigReal rswitch2 = rswitch*rswitch;
  BigReal rswitch3 = rswitch*rswitch2;
  BigReal rswitch4 = rswitch2*rswitch2;
  BigReal rswitch5 = rswitch2*rswitch3;

  /* Here we tabulate the integrals over the untruncated region. This assumes:

   1.) The energy and virial contribution can be well described by a mean field
       approximation (i.e. a constant).

   2.) The pair distribution function, g(r), is very close to unity on the
       interval (i.e. g(r) = 1 for r > switchdist).

   The mean field integrals are of the form:

   4*N^2*PI*int r^2 U(r) dr : for the energy

   (4/3)*N^2*PI*int r^3 dU(r)/dr dr : for the virial

   NB: An extra factor of 1/2 comes from double counting the number of
       interaction pairs (N*(N-1)/2, approximated as N^2).
  */
  BigReal int_U_gofr_A, int_rF_gofr_A, int_U_gofr_B, int_rF_gofr_B;
  if (simParams->switchingActive) {
    if (!simParams->vdwForceSwitching) {
      int_U_gofr_A = int_rF_gofr_A = (16*PI*(3*rcut4 + 9*rcut3*rswitch
                                      + 11*rcut2*rswitch2 + 9*rcut*rswitch3
                                      + 3*rswitch4)
                                      / (315*rcut5*rswitch5*(rcut + rswitch)
                                         *(rcut + rswitch)*(rcut + rswitch)));
      int_U_gofr_B = int_rF_gofr_B = (-16*PI / (3*(rcut + rswitch)
                                                *(rcut + rswitch)
                                                *(rcut + rswitch)));
    }
    else {
      /* BKR - This only includes the volume dependent portion of the energy
         correction. In order to fully correct to a true Lennard-Jones
         potential one also needs a "core" correction to account for the shift
         inside rswitch; this is a _global_ potential energy shift, regardless
         of volume. Neglecting this correction is equivalent to stating that
         the number of atoms in the whole system is much greater than the
         number within rswitch of any given atom.  Interestingly, this
         approximation gets better as rswitch shrinks whereas the tail
         correction gets worse. The extra term has _zero_ affect on the virial
         since the forces are unchanged.
      */
      BigReal lnr = log(rswitch/rcut);
      int_rF_gofr_A = 16*PI / (9*rswitch3*rcut3*(rcut3 + rswitch3));
      int_rF_gofr_B = 4*PI*lnr / (rcut3 - rswitch3);
      int_U_gofr_A = (2*PI*(5*rcut3 - 3*rswitch3)
                      / (9*rswitch3*rcut6*(rcut3 + rswitch3)));
      int_U_gofr_B = (-2*PI*(rswitch3 - rcut3 - 6*rcut3*lnr)
                      / (3*rcut3*(rcut3 - rswitch3)));
    }
  }
  else {
    int_rF_gofr_A = 8*PI / (9*rcut9);
    int_rF_gofr_B = -4*PI / (3*rcut3);
    int_U_gofr_A = 2*PI / (9*rcut9);
    int_U_gofr_B = -2*PI / (3*rcut3);
  }
  // If TI/FEP is on, these come back with values at alchLambda = 0 and are
  // thus equivalent to alchemical group 2.
  tail_corr_virial = int_rF_gofr_A*LJAvgA + int_rF_gofr_B*LJAvgB;
  tail_corr_ener = int_U_gofr_A*LJAvgA + int_U_gofr_B*LJAvgB;

  tail_corr_dUdl_1 = int_U_gofr_A*LJAvgA1 + int_U_gofr_B*LJAvgB1;
  tail_corr_virial_1 = int_rF_gofr_A*LJAvgA1 + int_rF_gofr_B*LJAvgB1;
}

// Convenience function to simplify lambda scaling.
BigReal Molecule::getEnergyTailCorr(const BigReal alchLambda){
  if (simParams->alchOn) {
    const BigReal vdw_lambda_1 = simParams->getVdwLambda(alchLambda);
    const BigReal vdw_lambda_2 = simParams->getVdwLambda(1-alchLambda);
    // NB: Rather than duplicate variables, dUdl_2 is stored as the energy.
    //     Put another way, dUdl_2 _is_ the energy, if alchLambda = 0.
    return vdw_lambda_1*tail_corr_dUdl_1 + vdw_lambda_2*tail_corr_ener;
  }
  else {
    return tail_corr_ener;
  }
}

// Convenience function to simplify lambda scaling.
BigReal Molecule::getVirialTailCorr(const BigReal alchLambda){
  if (simParams->alchOn) {
    const BigReal vdw_lambda_1 = simParams->getVdwLambda(alchLambda);
    const BigReal vdw_lambda_2 = simParams->getVdwLambda(1-alchLambda);
    // NB: Rather than duplicate variables, virial_2 is stored as the virial.
    //     Put another way, virial_2 _is_ the virial, if alchLambda = 0.
    return vdw_lambda_1*tail_corr_virial_1 + vdw_lambda_2*tail_corr_virial;
  }
  else {
    return tail_corr_virial;
  }
}
#endif

#ifdef MEM_OPT_VERSION
//idx1: atom1's exclusion check signature
//to check whether atom1 and atom2 are excluded from each other
int Molecule::checkExclByIdx(int idx1, int atom1, int atom2) const {

  int amin = exclChkSigPool[idx1].min;
  int amax = exclChkSigPool[idx1].max;
  int dist21 = atom2 - atom1;
  if ( dist21 < amin || dist21 > amax ) return 0;
  else return exclChkSigPool[idx1].flags[dist21-amin];

}
#else
int Molecule::checkexcl(int atom1, int atom2) const {

  int amin = all_exclusions[atom1].min;
  int amax = all_exclusions[atom1].max;
  if ( atom2 < amin || atom2 > amax ) return 0;
  else return all_exclusions[atom1].flags[atom2-amin];

}
#endif


/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for reading AMBER topology data    */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, Ambertoppar *amber_data)
{
  initialize(simParams,param);

  read_parm(amber_data);

#ifndef MEM_OPT_VERSION
  //LCPO
  if (simParams->LCPOOn)
    assignLCPOTypes( 1 );
#endif
}
/*      END OF FUNCTION Molecule      */


/************************************************************************/
/*                  */
/*      FUNCTION read_parm   */
/*                  */
/*   INPUTS:                */
/*  amber_data - AMBER data structure    */
/*                  */
/*  This function copys AMBER topology data to the corresponding data  */
/*   structures      */
/*                  */
/************************************************************************/

void Molecule::read_parm(Ambertoppar *amber_data)
{
#ifdef MEM_OPT_VERSION
    NAMD_die("When reading a compressed file or using the memory-optimized version, amber data is not supported!");    
#else
  int i,j,ntheth,nphih,current_index,a1,a2,
      max,min,index,found;

  if (!amber_data->data_read)
    NAMD_die("No data read from parm file yet!");

  // Copy atom informations
  numAtoms = amber_data->Natom;
  atoms = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];

  if(simParams->genCompressedPsf) {
      atomSegResids = new AtomSegResInfo[numAtoms];
  }

  if (atoms == NULL || atomNames == NULL )
    NAMD_die("memory allocation failed when reading atom information");
  ResidueLookupElem *tmpResLookup = resLookup;
  for (i=0; i<numAtoms; ++i)
  { atomNames[i].resname = nameArena->getNewArray(5);
    atomNames[i].atomname = nameArena->getNewArray(5);
    atomNames[i].atomtype = nameArena->getNewArray(5);
    if (atomNames[i].resname == NULL || atomNames[i].atomname == NULL || atomNames[i].atomtype == NULL)
      NAMD_die("memory allocation failed when reading atom information");
    for (j=0; j<4; ++j)
    { atomNames[i].resname[j] = amber_data->ResNames[amber_data->AtomRes[i]*4+j];
      atomNames[i].atomname[j] = amber_data->AtomNames[i*4+j];
      atomNames[i].atomtype[j] = amber_data->AtomSym[i*4+j];
    }
    atomNames[i].resname[4] = atomNames[i].atomname[4] = atomNames[i].atomtype[4] = '\0';
    strtok(atomNames[i].resname," ");
    strtok(atomNames[i].atomname," ");
    strtok(atomNames[i].atomtype," ");
    atoms[i].mass = amber_data->Masses[i];
    // Divide by 18.2223 to convert to charge in units of the electron charge
    atoms[i].charge = amber_data->Charges[i] / 18.2223;
    atoms[i].vdw_type = amber_data->Iac[i] - 1;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append("MAIN", amber_data->AtomRes[i]+1, i);

    if(atomSegResids) { //for compressing molecule information
        AtomSegResInfo *one = atomSegResids + i;
        memcpy(one->segname, "MAIN", strlen("MAIN")+1);
        one->resid = amber_data->AtomRes[i]+1;
    }
    

    /*  Determine the type of the atom (H or O) */
    atoms[i].status = UnknownAtom; // the default
    if ( simParams->ignoreMass ) {
    } else if (atoms[i].mass <= 0.05) {
      atoms[i].status |= LonepairAtom;
    } else if (atoms[i].mass < 1.0) {
      atoms[i].status |= DrudeAtom;
    } else if (atoms[i].mass <=3.5) {
      atoms[i].status |= HydrogenAtom;
    } else if ((atomNames[i].atomname[0] == 'O') && 
         (atoms[i].mass >= 14.0) && 
         (atoms[i].mass <= 18.0)) {
      atoms[i].status |= OxygenAtom;
    }
  }

// Note: In AMBER, the atom numbers in bond, angle and dihedral arrays are in fact
// (3*(atnum-1)). So we divide them by 3 to get the real indices of atom array. Also
// note that NAMD indexes arrays from 0 to NumAtoms-1.

  // Copy bond information
  // Fake bonds (bonds with 0 force constant) are ignored
  Real k, x0;
  numBonds = 0;
  if (amber_data->Nbonh + amber_data->Nbona > 0)
  { bonds = new Bond[amber_data->Nbonh + amber_data->Nbona];
    if (bonds == NULL || amber_data->Nbona < 0)
      NAMD_die("memory allocation failed when reading bond information");
    // Bonds WITH hydrogen
    for (i=0; i<amber_data->Nbonh; ++i)
    { bonds[numBonds].atom1 = amber_data->BondHAt1[i] / 3;
      bonds[numBonds].atom2 = amber_data->BondHAt2[i] / 3;
      bonds[numBonds].bond_type = amber_data->BondHNum[i] - 1;
      if (bonds[numBonds].atom1>=numAtoms || bonds[numBonds].atom2>=numAtoms ||
          bonds[numBonds].bond_type>=amber_data->Numbnd)
      { char err_msg[128];
        sprintf(err_msg, "BOND (WITH H) # %d OVERFLOW IN PARM FILE", i+1);
        NAMD_die(err_msg);
      }
      params->get_bond_params(&k,&x0,bonds[numBonds].bond_type);
      // if ( k != 0. ) ++numBonds;  // real bond
      ++numBonds;  // keep all bonds in case needed for rigid water
    }
    // Bonds WITHOUT hydrogen
    for (i=amber_data->Nbonh; i<amber_data->Nbonh+amber_data->Nbona; ++i)
    { bonds[numBonds].atom1 = amber_data->BondAt1[i-amber_data->Nbonh] / 3;
      bonds[numBonds].atom2 = amber_data->BondAt2[i-amber_data->Nbonh] / 3;
      bonds[numBonds].bond_type = amber_data->BondNum[i-amber_data->Nbonh] - 1;
      if (bonds[i].atom1>=numAtoms || bonds[i].atom2>=numAtoms ||
          bonds[i].bond_type>=amber_data->Numbnd)
      { char err_msg[128];
        sprintf(err_msg, "BOND (WITHOUT H) # %d OVERFLOW IN PARM FILE", i+1-amber_data->Nbonh);
        NAMD_die(err_msg);
      }
      params->get_bond_params(&k,&x0,bonds[numBonds].bond_type);
      // if ( k != 0. ) ++numBonds;  // real bond
      ++numBonds;  // keep all bonds in case needed for rigid water
    }
  }
  /*  Tell user about our subterfuge  */
  if ( numBonds !=  amber_data->Nbonh + amber_data->Nbona) {
    iout << iWARN << "Ignored " << amber_data->Nbonh + amber_data->Nbona - numBonds <<
            " bonds with zero force constants.\n" << endi;
    iout << iWARN <<
	"Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
  }

  // Copy angle information
  numAngles = amber_data->Ntheth + amber_data->Ntheta;
  if (numAngles > 0)
  { ntheth = amber_data->Ntheth;
    angles = new Angle[numAngles];
    if (angles == NULL || numAngles < ntheth)
      NAMD_die("memory allocation failed when reading angle information");
    // Angles WITH hydrogen
    for (i=0; i<ntheth; ++i)
    { angles[i].atom1 = amber_data->AngleHAt1[i] / 3;
      angles[i].atom2 = amber_data->AngleHAt2[i] / 3;
      angles[i].atom3 = amber_data->AngleHAt3[i] / 3;
      angles[i].angle_type = amber_data->AngleHNum[i] - 1;
      if (angles[i].atom1>=numAtoms || angles[i].atom2>=numAtoms ||
          angles[i].atom3>=numAtoms || angles[i].angle_type>=amber_data->Numang)
      { char err_msg[128];
        sprintf(err_msg, "ANGLE (WITH H) # %d OVERFLOW IN PARM FILE", i+1);
        NAMD_die(err_msg);
      }
    }
    // Angles WITHOUT hydrogen
    for (i=ntheth; i<numAngles; ++i)
    { angles[i].atom1 = amber_data->AngleAt1[i-ntheth] / 3;
      angles[i].atom2 = amber_data->AngleAt2[i-ntheth] / 3;
      angles[i].atom3 = amber_data->AngleAt3[i-ntheth] / 3;
      angles[i].angle_type = amber_data->AngleNum[i-ntheth] - 1;
      if (angles[i].atom1>=numAtoms || angles[i].atom2>=numAtoms ||
          angles[i].atom3>=numAtoms || angles[i].angle_type>=amber_data->Numang)
      { char err_msg[128];
        sprintf(err_msg, "ANGLE (WITHOUT H) # %d OVERFLOW IN PARM FILE", i+1-ntheth);
        NAMD_die(err_msg);
      }
    }
  }

  numExclusions =  0;
  // If readExclusions is TRUE, then we copy exclusions from parm
  // file; otherwise we skip the exclusions here and generate
  // them later in build_exclusions()
  if (simParams->readExclusions)
  { // Copy exclusion information
    // In Amber data structure, Iblo[] is the number of exclusions
    // for each atom; ExclAt[] is the atom index for the excluded atoms.
    exclusions = new Exclusion[amber_data->Nnb];
    if (exclusions == NULL &&  amber_data->Nnb > 0)
      NAMD_die("memory allocation failed when reading exclusion information");
    current_index = 0;
    for (i=0; i<numAtoms; ++i)
      for (j=0; j<amber_data->Iblo[i]; ++j)
      { if (current_index >= amber_data->Nnb)
        { char err_msg[128];
          sprintf(err_msg, "EXCLUSION INDEX EXCEEDS NUMBER OF EXLCUSIONS %d IN AMBER FILE, AT ATOM #%d\n",
            amber_data->Nnb, i+1);
          NAMD_die(err_msg);
        }
        // There's some 0 in the ExclAt[] list, which is strange
        // and redundant. In this case, I simply ignore such entries.
        if (amber_data->ExclAt[current_index] != 0)
        { // Subtract 1 to convert the index from the 1 to NumAtoms
          // used in the file to the 0 to NumAtoms-1 that we need
          a2 = amber_data->ExclAt[current_index] - 1;
          if (a2 < i)
          { // I assume the latter index be larger than the former
            // one, so that the same exclusion won't be double-counted;
            // if not, give error
            char err_msg[128];
            sprintf(err_msg, "Atom #%d has exclusion with atom #%d, in reverse order.", i+1, a2+1);
            NAMD_die(err_msg);
          }
          else if (a2 == i)
          { char err_msg[128];
            sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN AMBER FILE\n", i+1);
              NAMD_die(err_msg);
          }
          else if (a2 >= numAtoms)
          {  char err_msg[128];
             sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN AMBER FILE",
               a2+1, numAtoms, current_index+1);
             NAMD_die(err_msg);
          }
          exclusions[numExclusions].atom1 = i;
          exclusions[numExclusions].atom2 = a2;
          ++numExclusions;
        }
        ++current_index;
      }
    if (current_index < amber_data->Nnb)
    { char err_msg[128];
      sprintf(err_msg, "Num of exclusions recorded (%d) is smaller than what it's supposed to be (%d)",
        current_index,amber_data->Nnb);
      NAMD_die(err_msg);
    }
  }

  // Copy dihedral information
  numDihedrals = amber_data->Nphih + amber_data->Nphia;
  if (numDihedrals > 0)
  { nphih = amber_data->Nphih;
    dihedrals = new Dihedral[numDihedrals];
    if (dihedrals == NULL || numDihedrals < nphih)
      NAMD_die("memory allocation failed when reading dihedral information");
    // Dihedral WITH hydrogen
    for (i=0; i<nphih; ++i)
    { dihedrals[i].atom1 = amber_data->DihHAt1[i] / 3;
      dihedrals[i].atom2 = amber_data->DihHAt2[i] / 3;
      dihedrals[i].atom3 = amber_data->DihHAt3[i] / 3;
      dihedrals[i].atom4 = amber_data->DihHAt4[i] / 3;
      dihedrals[i].dihedral_type = amber_data->DihHNum[i] - 1;
    }
    // Dihedral WITHOUT hydrogen
    for (i=nphih; i<numDihedrals; ++i)
    { dihedrals[i].atom1 = amber_data->DihAt1[i-nphih] / 3;
      dihedrals[i].atom2 = amber_data->DihAt2[i-nphih] / 3;
      dihedrals[i].atom3 = amber_data->DihAt3[i-nphih] / 3;
      dihedrals[i].atom4 = amber_data->DihAt4[i-nphih] / 3;
      dihedrals[i].dihedral_type = amber_data->DihNum[i-nphih] - 1;
    }
  }
  // In AMBER parm file, dihedrals contain 1-4 exclusion infomation:
  // the 1st and 4th atoms have 1-4 nonbond interation. So we should
  // find them in the exclusion array and change their exclusion to
  // 1-4 type. However, there're two exceptions --
  // 1.If the third atom is negative, it means the end group
  //   interactions are to be ignored;
  // 2.If the fourth atom is negative, it means this is an improper.
  // For the above two cases, the actual atom index is the absolute
  // value of the atom number read; and there's no 1-4 interation
  // for these dihedrals.
  // If readExclusions is not TRUE, then we don't worry about
  // exclusions here.
  for (i=0; i<numDihedrals; ++i)
  { if (dihedrals[i].atom3 < 0 || dihedrals[i].atom4 < 0)
    { dihedrals[i].atom3 = abs(dihedrals[i].atom3);
      dihedrals[i].atom4 = abs(dihedrals[i].atom4);
    }
    else if (simParams->readExclusions)
    { if (dihedrals[i].atom1 < dihedrals[i].atom4)
        a1=dihedrals[i].atom1, a2=dihedrals[i].atom4;
      else
        a1=dihedrals[i].atom4, a2=dihedrals[i].atom1;
      // Since in the exclusion array, atom1 is guaranteed to be
      // ordered, we can do a binary serch to find it first.
      found = 0;
      min=0, max=numExclusions-1;
      while (!found && min<=max)
      { index = (min+max)/2;
        if (exclusions[index].atom1 == a1)
          found = 1;
        else if (exclusions[index].atom1 < a1)
          min = index+1;
        else
          max = index-1;
      }
      if (!found)
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
      // After finding atom1, we do a linear serch to find atom2,
      // in both directions.
      for (j=index-1; j>=0 && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; --j);
      if (j<0 || exclusions[j].atom1!=a1)
        for (j=index; j<numExclusions && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; ++j);
      if (j<numExclusions && exclusions[j].atom1==a1)
        exclusions[j].modified = 1;  // Change the exclusion type to 1-4
      else
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
    }
    if (dihedrals[i].atom1>=numAtoms || dihedrals[i].atom2>=numAtoms ||
        dihedrals[i].atom3>=numAtoms || dihedrals[i].atom4>=numAtoms ||
        dihedrals[i].dihedral_type>=amber_data->Nptra)
    { char err_msg[128];
      sprintf(err_msg, "DIHEDRAL # %d OVERFLOW IN PARM FILE", i+1);
      NAMD_die(err_msg);
    }
  }
  
  //  analyze the data and find the status of each atom
  numRealBonds = numBonds;
  build_atom_status();
#endif
}
/*      END OF FUNCTION read_parm    */


/************************************************************************/
/*                                                                      */
/*      FUNCTION Molecule                                               */
/*                                                                      */
/*  This is the constructor for reading GROMACS topology data           */
/*                                                                      */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param,
		   const GromacsTopFile *gromacsTopFile)
{
  initialize(simParams,param);

  read_parm(gromacsTopFile);

#ifndef MEM_OPT_VERSION
  //LCPO
  if (simParams->LCPOOn)
    assignLCPOTypes( 3 );
#endif
}
/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                                                                      */
/*      FUNCTION read_parm                                              */
/*                                                                      */
/*   INPUTS:                                                            */
/*  amber_data - AMBER data structure                                   */
/*                                                                      */
/*  This function copys AMBER topology data to the corresponding data   */
/*   structures                                                         */
/*                                                                      */
/************************************************************************/

void Molecule::read_parm(const GromacsTopFile *gf) {
#ifdef MEM_OPT_VERSION
    NAMD_die("When reading a compressed file or using the memory-optimized version, amber data is not supported!");    
#else
  /*  int i,j,ntheth,nphih,current_index,a1,a2,
      max,min,index,found;*/
  int i;
  
  // Initializes the atom array
  numAtoms = gf->getNumAtoms();
  atoms = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];

  if(simParams->genCompressedPsf) {
      atomSegResids = new AtomSegResInfo[numAtoms];
  }

  if (atoms == NULL || atomNames == NULL )
    NAMD_die("memory allocation failed when reading atom information");
  ResidueLookupElem *tmpResLookup = resLookup;

  // Copy the individual atoms over
  for (i=0; i<numAtoms; ++i) {
    char *resname = nameArena->getNewArray(11);
    char *atomname = nameArena->getNewArray(11);
    char *atomtype = nameArena->getNewArray(11);
    int resnum,typenum;
    Real charge,mass;

    if (resname == NULL || atomname == NULL || atomtype == NULL)
      NAMD_die("memory allocation failed when reading atom information");

    // get the data out of the GROMACS file
    gf->getAtom(i,&resnum,resname,atomname,atomtype,&typenum,
		&charge,&mass);

    atomNames[i].resname = resname;
    atomNames[i].atomname = atomname;
    atomNames[i].atomtype = atomtype;
    atoms[i].mass = mass;
    atoms[i].charge = charge;
    atoms[i].vdw_type = typenum;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append("MAIN", resnum+1, i);

    if(atomSegResids) { //for compressing molecule information
        AtomSegResInfo *one = atomSegResids + i;
        memcpy(one->segname, "MAIN", strlen("MAIN")+1);
        one->resid = resnum+1;
    }

    /*  Determine the type of the atom (H or O) */
    // XXX this cannot be done this way in GROMACS
    // For example, in dppc LO2 appears to be an oxygen.
    // And how do the hydrogens in CH3 etc factor in to this?
    atoms[i].status = UnknownAtom; // the default
    if ( simParams->ignoreMass ) {
    } else if (atoms[i].mass <= 0.05) {
      atoms[i].status |= LonepairAtom;
    } else if (atoms[i].mass < 1.0) {
      atoms[i].status |= DrudeAtom;
    } else if (atoms[i].mass <=3.5) {
      atoms[i].status |= HydrogenAtom;
    } else if ((atomNames[i].atomname[0] == 'O') && 
         (atoms[i].mass >= 14.0) && 
         (atoms[i].mass <= 18.0)) {
      atoms[i].status |= OxygenAtom;
    }
  }

  // Copy bond information
  numBonds = gf->getNumBonds();
  bonds = new Bond[numBonds];
  if (bonds == NULL)
    NAMD_die("memory allocation failed when reading bond information");
  for(i=0;i<numBonds;i++) {
    int type; // to convert the type correctly
    int atom1,atom2;
    gf->getBond(i,&atom1,&atom2,&type);
    bonds[i].atom1 = atom1;
    bonds[i].atom2 = atom2;
    bonds[i].bond_type = (Index)type;
  }

  // Copy angle information
  numAngles = gf->getNumAngles();
  angles = new Angle[numAngles];
  if (angles == NULL)
    NAMD_die("memory allocation failed when reading angle information");
  for(i=0;i<numAngles;i++) {
    int type; // to convert the type correctly
    int atom1,atom2,atom3;
    gf->getAngle(i,&atom1,&atom2,&atom3,&type);

    angles[i].atom1 = atom1;
    angles[i].atom2 = atom2;
    angles[i].atom3 = atom3;

    angles[i].angle_type=type;
  }

  numExclusions =  0;
  exclusions = new Exclusion[numExclusions];

  /*
  // If readExclusions is TRUE, then we copy exclusions from parm
  // file; otherwise we skip the exclusions here and generate
  // them later in build_exclusions()
  if (simParams->readExclusions)
  { // Copy exclusion information
    // In Amber data structure, Iblo[] is the number of exclusions
    // for each atom; ExclAt[] is the atom index for the excluded atoms.
    exclusions = new Exclusion[amber_data->Nnb];
    if (exclusions == NULL &&  amber_data->Nnb > 0)
      NAMD_die("memory allocation failed when reading exclusion information");
    current_index = 0;
    for (i=0; i<numAtoms; ++i)
      for (j=0; j<amber_data->Iblo[i]; ++j)
      { if (current_index >= amber_data->Nnb)
        { char err_msg[128];
          sprintf(err_msg, "EXCLUSION INDEX EXCEEDS NUMBER OF EXLCUSIONS %d IN AMBER FILE, AT ATOM #%d\n",
            amber_data->Nnb, i+1);
          NAMD_die(err_msg);
        }
        // There's some 0 in the ExclAt[] list, which is strange
        // and redundant. In this case, I simply ignore such entries.
        if (amber_data->ExclAt[current_index] != 0)
        { // Subtract 1 to convert the index from the 1 to NumAtoms
          // used in the file to the 0 to NumAtoms-1 that we need
          a2 = amber_data->ExclAt[current_index] - 1;
          if (a2 < i)
          { // I assume the latter index be larger than the former
            // one, so that the same exclusion won't be double-counted;
            // if not, give error
            char err_msg[128];
            sprintf(err_msg, "Atom #%d has exclusion with atom #%d, in reverse order.", i+1, a2+1);
            NAMD_die(err_msg);
          }
          else if (a2 == i)
          { char err_msg[128];
            sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN AMBER FILE\n", i+1);
              NAMD_die(err_msg);
          }
          else if (a2 >= numAtoms)
          {  char err_msg[128];
             sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN AMBER FILE",
               a2+1, numAtoms, current_index+1);
             NAMD_die(err_msg);
          }
          exclusions[numExclusions].atom1 = i;
          exclusions[numExclusions].atom2 = a2;
          ++numExclusions;
        }
        ++current_index;
      }
    if (current_index < amber_data->Nnb)
    { char err_msg[128];
      sprintf(err_msg, "Num of exclusions recorded (%d) is smaller than what it's supposed to be (%d)",
        current_index,amber_data->Nnb);
      NAMD_die(err_msg);
    }
  }
  */

  // Copy dihedral information
  numDihedrals = gf->getNumDihedrals();
  dihedrals = new Dihedral[numDihedrals];
  if (dihedrals == NULL)
    NAMD_die("memory allocation failed when reading dihedral information");
  for(i=0;i<numDihedrals;i++) {
    int type; // to convert the type correctly
    int atom1,atom2,atom3,atom4;
    gf->getDihedral(i,&atom1,&atom2,&atom3,&atom4,&type);
    dihedrals[i].atom1 = atom1;
    dihedrals[i].atom2 = atom2;
    dihedrals[i].atom3 = atom3;
    dihedrals[i].atom4 = atom4;
    dihedrals[i].dihedral_type = type;
  }

#if GROMACS_PAIR
  // JLai modifications on August 16th, 2012
  numPair = gf->getNumPair();
  numLJPair = gf->getNumLJPair();
  //std::cout << "Number of LJ pairs defined: " << numLJPair << "\n";
  indxLJA = new int[numLJPair];
  indxLJB = new int[numLJPair];
  pairC6 = new Real[numLJPair];
  pairC12 = new Real[numLJPair];
  gromacsPair_type = new int[numLJPair];
  const_cast<GromacsTopFile*>(gf)->getPairLJArrays2(indxLJA, indxLJB, pairC6, pairC12);
  gromacsPair = new GromacsPair[numLJPair];
  for(int i=0; i < numLJPair; i++) {
      gromacsPair_type[i] = i;
      gromacsPair[i].atom1 = indxLJA[i];
      gromacsPair[i].atom2 = indxLJB[i];
      gromacsPair[i].pairC6  = pairC6[i];
      gromacsPair[i].pairC12 = pairC12[i];
      //std::cout << "GromacsPairInitialization: " << gromacsPair[i].atom1 << " " << gromacsPair[i].atom2 << " " << gromacsPair[i].pairC6 << " " << gromacsPair[i].pairC12 << "\n";
      gromacsPair[i].gromacsPair_type = gromacsPair_type[i];
  }
  
  pointerToLJBeg = new int[numAtoms];
  pointerToLJEnd = new int[numAtoms];
  int oldIndex = -1;
  for(int i=0; i < numAtoms; i++) {
    pointerToLJBeg[i] = -1;
    pointerToLJEnd[i] = -2;
  }
  for(int i=0; i < numLJPair; i++) {
    if(pointerToLJBeg[indxLJA[i]] == -1) {
      pointerToLJBeg[indxLJA[i]] = i;
      oldIndex = indxLJA[i];
    }
    pointerToLJEnd[oldIndex] = i; 
  }

  // Initialize Gaussian arrays
  numGaussPair = gf->getNumGaussPair();
  indxGaussA = new int[numGaussPair];
  indxGaussB = new int[numGaussPair];
  gA = new Real[numGaussPair];
  gMu1 = new Real[numGaussPair];
  giSigma1 = new Real[numGaussPair];
  gMu2 = new Real[numGaussPair];
  giSigma2 = new Real[numGaussPair];
  gRepulsive = new Real[numGaussPair];
  const_cast<GromacsTopFile*>(gf)->getPairGaussArrays2(indxGaussA, indxGaussB, gA, gMu1, giSigma1, gMu2, giSigma2, gRepulsive);
  
  // Create an array of pointers to index indxGaussA
  pointerToGaussBeg = new int[numAtoms];
  pointerToGaussEnd = new int[numAtoms];
  for(int i=0; i < numAtoms; i++) {
    pointerToGaussBeg[i] = -1;
    pointerToGaussEnd[i] = -2;
  }
  oldIndex = -1;
  for(int i=0; i < numGaussPair; i++) {
    if(pointerToGaussBeg[indxGaussA[i]] == -1) { 
      pointerToGaussBeg[indxGaussA[i]] = i;
      oldIndex = indxGaussA[i];
    }
    pointerToGaussEnd[oldIndex] = i;
  }
  
  iout << iINFO << "Finished reading explicit pair from Gromacs file:\n" << 
    iINFO << "Found a total of: " << numPair << " explicit pairs--of which: " <<
    numLJPair << " are LJ style pairs and " << numGaussPair << 
    " are Gaussian style pairs.\n" << endi; //(Note: A->B is counted twice as A->B and B->A)\n" << endi;
#endif

  // Start of JLai Modifications August 16th, 2012 
#if GROMACS_EXCLUSIONS
  // Initialize exclusion information
  int numExclusions = gf->getNumExclusions();
  int* atom1 = new int[numExclusions];
  int* atom2 = new int[numExclusions];
  for(int j=0; j<numExclusions;j++) {
      atom1[j] = 0;
      atom2[j] = 0;
  }
  // Get exclusion arrays from gf module 
  const_cast<GromacsTopFile*>(gf)->getExclusions(atom1,atom2);
  read_exclusions(atom1,atom2,numExclusions);

  // Dump array 
  delete [] atom1;
  delete [] atom2;
#endif
  /*
  // In AMBER parm file, dihedrals contain 1-4 exclusion infomation:
  // the 1st and 4th atoms have 1-4 nonbond interation. So we should
  // find them in the exclusion array and change their exclusion to
  // 1-4 type. However, there're two exceptions --
  // 1.If the third atom is negative, it means the end group
  //   interactions are to be ignored;
  // 2.If the fourth atom is negative, it means this is an improper.
  // For the above two cases, the actual atom index is the absolute
  // value of the atom number read; and there's no 1-4 interation
  // for these dihedrals.
  // If readExclusions is not TRUE, then we don't worry about
  // exclusions here.
  for (i=0; i<numDihedrals; ++i)
  { if (dihedrals[i].atom3 < 0 || dihedrals[i].atom4 < 0)
    { dihedrals[i].atom3 = abs(dihedrals[i].atom3);
      dihedrals[i].atom4 = abs(dihedrals[i].atom4);
    }
    else if (simParams->readExclusions)
    { if (dihedrals[i].atom1 < dihedrals[i].atom4)
        a1=dihedrals[i].atom1, a2=dihedrals[i].atom4;
      else
        a1=dihedrals[i].atom4, a2=dihedrals[i].atom1;
      // Since in the exclusion array, atom1 is guaranteed to be
      // ordered, we can do a binary serch to find it first.
      found = 0;
      min=0, max=numExclusions-1;
      while (!found && min<=max)
      { index = (min+max)/2;
        if (exclusions[index].atom1 == a1)
          found = 1;
        else if (exclusions[index].atom1 < a1)
          min = index+1;
        else
          max = index-1;
      }
      if (!found)
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
      // After finding atom1, we do a linear serch to find atom2,
      // in both directions.
      for (j=index-1; j>=0 && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; --j);
      if (j<0 || exclusions[j].atom1!=a1)
        for (j=index; j<numExclusions && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; ++j);
      if (j<numExclusions && exclusions[j].atom1==a1)
        exclusions[j].modified = 1;  // Change the exclusion type to 1-4
      else
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
    }
    if (dihedrals[i].atom1>=numAtoms || dihedrals[i].atom2>=numAtoms ||
        dihedrals[i].atom3>=numAtoms || dihedrals[i].atom4>=numAtoms ||
        dihedrals[i].dihedral_type>=amber_data->Nptra)
    { char err_msg[128];
      sprintf(err_msg, "DIHEDRAL # %d OVERFLOW IN PARM FILE", i+1);
      NAMD_die(err_msg);
    }
  }
  */
  //  analyze the data and find the status of each atom
  numRealBonds = numBonds;
  build_atom_status();
#endif
}
/*      END OF FUNCTION read_parm    */

#ifndef MEM_OPT_VERSION
/*
int32 *Molecule::get_bonds_for_atom(int anum){
    NAMD_die("In bonds for atom!");
    return bondsByAtom[anum];
}

Bond *Molecule::get_bond(int bnum){
    NAMD_die("In get_bond!");
    return &bonds[bnum];
}
*/
#endif

#ifdef MEM_OPT_VERSION
//return the index of the new mass in the mass pool
Index Molecule::insert_new_mass(Real newMass){
    //first search
    for(int i=massPoolSize-1; i>=0; i--){
        if(fabs(atomMassPool[i]-newMass)<=1e-6)
            return i;
    }
    //otherwise increase one entry for the new mass
    Real *tmp = new Real[massPoolSize+1];
    tmp[massPoolSize] = newMass;
    memcpy((void *)tmp, (const void *)atomMassPool, sizeof(Real)*massPoolSize);
    delete [] atomMassPool;
    atomMassPool = tmp;
    massPoolSize++;
    return (Index)(massPoolSize-1);
}

void Molecule::addNewExclSigPool(const vector<ExclusionSignature>& newExclSigPool){
    ExclusionSignature *tmpExclSigPool = new ExclusionSignature[exclSigPoolSize+newExclSigPool.size()];
    for(int i=0; i<exclSigPoolSize; i++)
        tmpExclSigPool[i] = exclSigPool[i];
    for(int i=0; i<newExclSigPool.size(); i++)
        tmpExclSigPool[i+exclSigPoolSize] = newExclSigPool[i];

    exclSigPoolSize += newExclSigPool.size();
    exclSigPool = tmpExclSigPool;
}

void TupleSignature::pack(MOStream *msg){
    msg->put((short)tupleType);
    msg->put(numOffset);
    msg->put(numOffset, offset);    
    msg->put(tupleParamType);
    msg->put(isReal);
}

void TupleSignature::unpack(MIStream *msg){
    short ttype;
    msg->get(ttype);
    tupleType = (TupleSigType)ttype;

    msg->get(numOffset);
    delete [] offset;
    offset = new int[numOffset];
    msg->get(numOffset*sizeof(int), (char *)offset);
    
    msg->get(tupleParamType);
    msg->get(isReal);
}

void AtomSignature::pack(MOStream *msg){
    msg->put(bondCnt);
    for(int i=0; i<bondCnt; i++)
        bondSigs[i].pack(msg);

    msg->put(angleCnt);
    for(int i=0; i<angleCnt; i++)
        angleSigs[i].pack(msg);

    msg->put(dihedralCnt);
    for(int i=0; i<dihedralCnt; i++)
        dihedralSigs[i].pack(msg);

    msg->put(improperCnt);
    for(int i=0; i<improperCnt; i++)
        improperSigs[i].pack(msg);

    msg->put(crosstermCnt);
    for(int i=0; i<crosstermCnt; i++)
        crosstermSigs[i].pack(msg);
    
    // JLai
    msg->put(gromacsPairCnt);
    for(int i=0; i<gromacsPairCnt; i++)
	gromacsPairSigs[i].pack(msg);
}

void AtomSignature::unpack(MIStream *msg){
    msg->get(bondCnt);
    delete [] bondSigs;
    if(bondCnt>0){
        bondSigs = new TupleSignature[bondCnt];
        for(int i=0; i<bondCnt; i++)
            bondSigs[i].unpack(msg);
    } else bondSigs = NULL;

    msg->get(angleCnt);
    delete [] angleSigs;
    if(angleCnt>0){
        angleSigs = new TupleSignature[angleCnt];
        for(int i=0; i<angleCnt; i++)
            angleSigs[i].unpack(msg);
    } else angleSigs = NULL;

    msg->get(dihedralCnt);
    delete [] dihedralSigs;
    if(dihedralCnt>0){
        dihedralSigs = new TupleSignature[dihedralCnt];
        for(int i=0; i<dihedralCnt; i++)
            dihedralSigs[i].unpack(msg);
    } else dihedralSigs = NULL;

    msg->get(improperCnt);
    delete [] improperSigs;
    if(improperCnt>0){
        improperSigs = new TupleSignature[improperCnt];
        for(int i=0; i<improperCnt; i++)
            improperSigs[i].unpack(msg);
    } else improperSigs = NULL;

    msg->get(crosstermCnt);
    delete [] crosstermSigs;
    if(crosstermCnt>0){
        crosstermSigs = new TupleSignature[crosstermCnt];
        for(int i=0; i<crosstermCnt; i++)
            crosstermSigs[i].unpack(msg);
    } else crosstermSigs = NULL;

    // JLai

    msg->get(gromacsPairCnt);
    delete [] gromacsPairSigs;
    if(gromacsPairCnt>0){
	gromacsPairSigs = new TupleSignature[gromacsPairCnt];
	for(int i=0; i<gromacsPairCnt; i++)
	    gromacsPairSigs[i].unpack(msg);
    } else gromacsPairSigs = NULL;

    // End of JLai

}

void AtomSignature::removeEmptyTupleSigs(){
    int origTupleCnt;
    int idx;
    TupleSignature *tupleSigs;
    TupleSignature *newTupleSigs;

    //bonds
  {
    origTupleCnt = bondCnt;
    tupleSigs= bondSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            bondCnt--;
    }
    if(bondCnt==0){
        delete [] tupleSigs;
        bondSigs = NULL;
    }else if(bondCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[bondCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        bondSigs = newTupleSigs;
    }
  }

    //angles
  {
    origTupleCnt = angleCnt;
    tupleSigs = angleSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            angleCnt--;
    }
    if(angleCnt==0){
        delete [] tupleSigs;
        angleSigs = NULL;
    }else if(angleCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[angleCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        angleSigs = newTupleSigs;
    }
  }

    //dihedrals
  {
    origTupleCnt = dihedralCnt;
    tupleSigs = dihedralSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            dihedralCnt--;
    }
    if(dihedralCnt==0){
        delete [] tupleSigs;
        dihedralSigs = NULL;
    }else if(dihedralCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[dihedralCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        dihedralSigs = newTupleSigs;        
    }
  }


    //impropers
  {
    origTupleCnt = improperCnt;
    tupleSigs = improperSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            improperCnt--;
    }
    if(improperCnt==0){
        delete [] tupleSigs;
        improperSigs = NULL;
    }else if(improperCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[improperCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        improperSigs = newTupleSigs;
    }    
  }

    //crossterms
  {
    origTupleCnt = crosstermCnt;
    tupleSigs = crosstermSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            crosstermCnt--;
    }
    if(crosstermCnt==0){
        delete [] tupleSigs;
        crosstermSigs = NULL;
    }else if(crosstermCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[crosstermCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        crosstermSigs = newTupleSigs;
    }    
  }

  // JLai
  // gromacs pair force
  {
    origTupleCnt = gromacsPairCnt;
    tupleSigs = gromacsPairSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            gromacsPairCnt--;
    }
    if(gromacsPairCnt==0){
        delete [] tupleSigs;
        gromacsPairSigs = NULL;
    }else if(gromacsPairCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[gromacsPairCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        gromacsPairSigs = newTupleSigs;
    }    
  }

  // End of JLai

}

void ExclusionSignature::removeEmptyOffset(){
	int newCnt=0;
	for(int i=0; i<fullExclCnt; i++){
	    if(fullOffset[i]==0) continue;
	    newCnt++;
	}
    if(newCnt==0){
        fullExclCnt = 0;
        delete [] fullOffset;
        fullOffset = NULL;
    }else if(newCnt!=fullExclCnt){
        int *tmpOffset = new int[newCnt];
    	newCnt=0;
    	for(int i=0; i<fullExclCnt; i++){
    	    if(fullOffset[i]==0) continue;
    	    tmpOffset[newCnt] = fullOffset[i];
    	    newCnt++;
    	}
    	delete [] fullOffset;
    	fullOffset = tmpOffset;
        fullExclCnt = newCnt;
    }
	
	
	newCnt=0;
	for(int i=0; i<modExclCnt; i++){
	    if(modOffset[i]==0) continue;
	    newCnt++;
	}
    if(newCnt==0){
        modExclCnt = 0;
        delete [] modOffset;
        modOffset = NULL;
    }else if(newCnt!=modExclCnt){
        int *tmpOffset = new int[newCnt];
        newCnt=0;
        for(int i=0; i<modExclCnt; i++){
            if(modOffset[i]==0) continue;
            tmpOffset[newCnt] = modOffset[i];
            newCnt++;
        }
        delete [] modOffset;
        modOffset = tmpOffset;
        modExclCnt = newCnt;
    }	
}

//returns the index of the offset. If not found, -1 is returned
//fullOrMod indicates where is the offset found. 0 indicates in
//the full exclusion lists, 1 indicates in the modified exclusion
//lists
int ExclusionSignature::findOffset(int offset, int *fullOrMod){
	//assuming all offsets have been sorted increasingly
	//so that binary search could be used	
	int retidx = -1;
	
	*fullOrMod = 0;	
	int low = 0;
	int high = fullExclCnt-1;
	int mid = (low+high)/2;
	while(low<=high){
		if(offset<fullOffset[mid]){
			high = mid-1;
			mid = (high+low)/2;						
		}else if(offset>fullOffset[mid]){
			low = mid+1;
			mid = (high+low)/2;
		}else{
			retidx = mid;
			break;
		}		
	}
	if(retidx!=-1) return retidx;
	
	*fullOrMod = 1;	
	low = 0;
	high = modExclCnt-1;
	mid = (low+high)/2;
	while(low<=high){
		if(offset<modOffset[mid]){
			high = mid-1;
			mid = (high+low)/2;						
		}else if(offset>modOffset[mid]){
			low = mid+1;
			mid = (high+low)/2;
		}else{
			retidx = mid;
			break;
		}		
	}
	return retidx;	
}

void ExclusionSignature::pack(MOStream *msg){
    msg->put(fullExclCnt);    
    msg->put(fullExclCnt, fullOffset);
    msg->put(modExclCnt);    
    msg->put(modExclCnt, modOffset);
}

void ExclusionSignature::unpack(MIStream *msg){    
    msg->get(fullExclCnt);
    delete [] fullOffset;
    fullOffset = new int[fullExclCnt];
    msg->get(fullExclCnt*sizeof(int), (char *)fullOffset);    
    msg->get(modExclCnt);
    delete [] modOffset;
    modOffset = new int[modExclCnt];
    msg->get(modExclCnt*sizeof(int), (char *)modOffset);    
#ifdef NAMD_CUDA
    buildTuples();
#endif
}
#endif

#endif  // MOLECULE2_C defined = second object file

