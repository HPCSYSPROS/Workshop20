/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "common.h"
#include <vector>


// status elements, used for Atom status 
#define UnknownAtom      0x00
#define HydrogenAtom     0x01
#define OxygenAtom       0x02
#define HBDonorAtom      0x04
#define HBAcceptorAtom   0x08
#define HBAntecedentAtom 0x10
#define HBHydrogenAtom   0x20
#define LonepairAtom     0x40
#define DrudeAtom        0x80


typedef int Index;		//  Used for index into arrays
					//  or parameters

typedef struct atom_name_info
{
	char *resname;
	char *atomname;
	char *atomtype;
} AtomNameInfo;

typedef struct atom_constants
{
	Real mass;
	Real charge;
	Index vdw_type;
	int32 status;	         // flags telling about this atom
	int32 partner;             // connecting atom, for hydrogens
	int32 hydrogenList;	// index of atom in hydrogenGroup list
} Atom;

typedef struct bond
{
	int32 atom1;
	int32 atom2;
	Index bond_type;
} Bond;

typedef struct angle
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	Index angle_type;
} Angle;

typedef struct dihedral
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	Index dihedral_type;
} Dihedral;

typedef struct improper
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	Index improper_type;
} Improper;

typedef struct crossterm
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	int32 atom5;
	int32 atom6;
	int32 atom7;
	int32 atom8;
	Index crossterm_type;
} Crossterm;

typedef struct gromacsPair
{
    int32 atom1;
    int32 atom2;
    Real pairC12;
    Real pairC6;
    Index gromacsPair_type;
} GromacsPair;


// DRUDE: data read from PSF
typedef struct drude_constants  // supplement Atom data
{                        // create array length N of these
  Real alpha;
  Real thole;
} DrudeConst;

typedef struct lphost    // lone pair host
{
  int32 atom1;
  int32 atom2;
  int32 atom3;
  int32 atom4;
  Real distance;
  Real angle;
  Real dihedral;
} Lphost;

typedef struct aniso     // anisotropic term
{
  int32 atom1;
  int32 atom2;
  int32 atom3;
  int32 atom4;
  Real k11;
  Real k22;
  Real k33;
} Aniso;

typedef struct thole    // Thole terms - constructed implicitly from exclusions
{
  int32 atom1;  // first heavy atom
  int32 atom2;  // Drude particle of first heavy atom
  int32 atom3;  // second heavy atom
  int32 atom4;  // Drude particle of second heavy atom
  Real aa;      // constant (alpha_i * alpha_j)^(-6) / (thole_i + thole_j)
  Real qq;      // combined charge of Drudes (C * q_{i+1} * q_{j+1})
} Thole;
// DRUDE


class Exclusion
{
public:
	Exclusion(void) : modified(0) {;}
	Exclusion(int a1, int a2, int mod = 0) :
		atom1(a1), atom2(a2), modified(mod) {;}
	int32 atom1;
	int32 atom2;
	Index modified;
	int hash(void) const
	{
		return atom1 + atom2;
	}
	int operator==(const Exclusion &o) const
	{
		return atom1 == o.atom1 && atom2 == o.atom2;
	}
	int operator<(const Exclusion &o) const
	{
		return
		(
		  ( atom1 < o.atom1 ) ||
		  ( atom1 == o.atom1 && atom2 < o.atom2 )
		);
	}
};


class MOStream;
class MIStream;

enum TupleSigType {BOND=0, ANGLE, DIHEDRAL, IMPROPER, DONOR, ACCEPTOR, CROSSTERM, EXCLUSION};

inline unsigned int circShift(unsigned int h, unsigned int by) {
  const unsigned int intBits=8*sizeof(unsigned int);
  by%=intBits;
  return (h<<by)|(h>>(intBits-by));
}

class TupleSignature{
public:
    //This field indicates which class this tuple belongs to (bond or angle or ...)
    TupleSigType tupleType;
    int numOffset;
    //the relative offset from the atom index
    int *offset;   

    //This type is determined by the Parameter object
    //for Exclusion, this indicates the modified field
    Index tupleParamType;

    //indicate this tuple is specified in the psf file, not from other files.
    //this is added due to the extraBonds feature. All tuples from extraBonds
    //are not real!
    char isReal;

public:
    TupleSignature(){        
        offset = NULL;
        isReal = 1;
    }

    TupleSignature(int n, TupleSigType t, Index paramType, char ir=1){
        tupleType = t;
        numOffset = n;
        offset = new int[n];        
        tupleParamType = paramType;
        isReal = ir;
    }
    TupleSignature(const TupleSignature& oneSig){
        tupleType = oneSig.tupleType;
        numOffset = oneSig.numOffset;
        offset = new int[numOffset];
        setOffsets(oneSig.offset);        
        tupleParamType = oneSig.tupleParamType;
        isReal = oneSig.isReal;
    }
    TupleSignature &operator=(const TupleSignature& oneSig){
        tupleType = oneSig.tupleType;
        numOffset = oneSig.numOffset;
        if(offset) delete [] offset;
        offset = new int[numOffset];
        setOffsets(oneSig.offset);        
        tupleParamType = oneSig.tupleParamType;
        isReal = oneSig.isReal;
        return *this;
    }
    int operator==(const TupleSignature& sig) const{
	    if(tupleType!=sig.tupleType)
	        return 0;

	    if(tupleParamType != sig.tupleParamType)
	        return 0;

        if(isReal != sig.isReal)
            return 0;

		if(numOffset != sig.numOffset) return 0;
		
	    int equalCnt=0;
	    	    
	    for(int i=0; i<numOffset; i++){
	        equalCnt += (offset[i]==sig.offset[i]);      
	    }
	    return equalCnt==numOffset;
	}
    ~TupleSignature(){
        if(offset) delete[] offset;
    }
    void setOffsets(int *offs){
        for(int i=0; i<numOffset; i++)
            offset[i] = offs[i];
	//sort the offset in increasing order
	//based on the input files, this offset is almost sorted increasingly
	//Therefore using insertion sort
	/*switch(numOffset){
	    case 1:
		break;
	    case 2:
		if(offset[0]>offset[1]){
		    int tmp = offset[0];
		    offset[0] = offset[1];
		    offset[1] = tmp;
		}
		break;
	    default: //insertion sort
		for(int ii=1; ii<numOffset; ii++){
		    int val = offset[ii];
		    for(int jj=ii-1; jj>=0; jj--){
			if(offset[jj]<val){
			    offset[jj+1] = offset[jj];
			    offset[jj] = val;
			}else
			    break;
		    }
		}
	}*/
    }
    void setEmpty(){
        delete [] offset;
        offset = NULL;
    }
    int isEmpty(){
        return offset==NULL;
    }
    void output(FILE *ofp){
        for(int i=0; i<numOffset; i++)
            fprintf(ofp, "%d ", offset[i]);
        fprintf(ofp, "| %d | %d\n", tupleParamType, isReal);         
    }
  
    int hash() const {
      unsigned int code = tupleType;
      unsigned int codesz = 8 * sizeof(int);
      unsigned int shift = codesz / numOffset;
      
      if (shift == 0) shift=1;
      unsigned int i;
      for(i=0; i < numOffset; i++) {
        code = circShift(code,shift);
        code ^= offset[i];
      }
      return code;
    }
  
    void pack(MOStream *msg);
    void unpack(MIStream *msg);
};

//represents the signatures for atoms
class AtomSignature{
public:
    int bondCnt;
    int angleCnt;
    int dihedralCnt;
    int improperCnt;
    int crosstermCnt;
    // JLai
    int gromacsPairCnt;

    TupleSignature *bondSigs;
    TupleSignature *angleSigs;
    TupleSignature *dihedralSigs;
    TupleSignature *improperSigs;
    TupleSignature *crosstermSigs;
    // JLai
    TupleSignature *gromacsPairSigs;

    AtomSignature(){
        bondCnt=angleCnt=dihedralCnt=improperCnt=crosstermCnt=0;
	// JLai
	gromacsPairCnt=0;
        bondSigs = NULL;
        angleSigs = NULL;
        dihedralSigs = NULL;
        improperSigs = NULL;
        crosstermSigs = NULL;
	// JLai
	gromacsPairSigs = NULL;
    }
    AtomSignature(const AtomSignature &sig){
        bondSigs = NULL;
        angleSigs = NULL;
        dihedralSigs = NULL;
        improperSigs = NULL;
        crosstermSigs = NULL;
	// JLai
	gromacsPairSigs = NULL;

        bondCnt = sig.bondCnt;
        if(bondCnt>0){            
            bondSigs = new TupleSignature[bondCnt];
            for(int i=0; i<bondCnt; i++)
                bondSigs[i] = sig.bondSigs[i];
        }

        angleCnt = sig.angleCnt;
        if(angleCnt>0){            
            angleSigs = new TupleSignature[angleCnt];
            for(int i=0; i<angleCnt; i++)
                angleSigs[i] = sig.angleSigs[i];
        }
        
        dihedralCnt = sig.dihedralCnt;
        if(dihedralCnt>0){            
            dihedralSigs = new TupleSignature[dihedralCnt];
            for(int i=0; i<dihedralCnt; i++)
                dihedralSigs[i] = sig.dihedralSigs[i];
        }
        
        improperCnt = sig.improperCnt;
        if(improperCnt>0){
            improperSigs = new TupleSignature[improperCnt];
            for(int i=0; i<improperCnt; i++)
                improperSigs[i] = sig.improperSigs[i];
        }      

        crosstermCnt = sig.crosstermCnt;
        if(crosstermCnt>0){
            crosstermSigs = new TupleSignature[crosstermCnt];
            for(int i=0; i<crosstermCnt; i++)
                crosstermSigs[i] = sig.crosstermSigs[i];
        }        

	// JLai
	gromacsPairCnt = sig.gromacsPairCnt;
	if(gromacsPairCnt>0){
	    gromacsPairSigs = new TupleSignature[gromacsPairCnt];
	    for(int i=0; i<gromacsPairCnt; i++) {
		gromacsPairSigs[i] = sig.gromacsPairSigs[i];
	    }
	}
    }
    AtomSignature& operator=(const AtomSignature& sig){        
        bondCnt = sig.bondCnt;
        if(bondSigs) delete [] bondSigs;
        if(bondCnt>0){
            bondSigs = new TupleSignature[bondCnt];
            for(int i=0; i<bondCnt; i++)
                bondSigs[i] = sig.bondSigs[i];
        }else
            bondSigs = NULL;

        angleCnt = sig.angleCnt;
        if(angleSigs) delete [] angleSigs;
        if(angleCnt>0){
            angleSigs = new TupleSignature[angleCnt];
            for(int i=0; i<angleCnt; i++)
                angleSigs[i] = sig.angleSigs[i];
        }else
            angleSigs = NULL;
        
        dihedralCnt = sig.dihedralCnt;
        if(dihedralSigs) delete [] dihedralSigs;
        if(dihedralCnt>0){
            dihedralSigs = new TupleSignature[dihedralCnt];
            for(int i=0; i<dihedralCnt; i++)
                dihedralSigs[i] = sig.dihedralSigs[i];
        }else
            dihedralSigs = NULL;
        
        improperCnt = sig.improperCnt;
        if(improperSigs) delete [] improperSigs;
        if(improperCnt>0){
            improperSigs = new TupleSignature[improperCnt];
            for(int i=0; i<improperCnt; i++)
                improperSigs[i] = sig.improperSigs[i];
        }else
            improperSigs = NULL;      

        crosstermCnt = sig.crosstermCnt;
        if(crosstermSigs) delete [] crosstermSigs;
        if(crosstermCnt>0){
            crosstermSigs = new TupleSignature[crosstermCnt];
            for(int i=0; i<crosstermCnt; i++)
                crosstermSigs[i] = sig.crosstermSigs[i];
        }else
            crosstermSigs = NULL;

	// JLai
	gromacsPairCnt = sig.gromacsPairCnt;
	if(gromacsPairSigs) delete [] gromacsPairSigs;
	if(gromacsPairCnt>0){
	    gromacsPairSigs = new TupleSignature[gromacsPairCnt];
	    for(int i=0; i<gromacsPairCnt; i++)
		gromacsPairSigs[i] = sig.gromacsPairSigs[i];
	}else
	    gromacsPairSigs = NULL;

        return *this;
    }
    int operator==(const AtomSignature& sig) const{
	    if(bondCnt!=sig.bondCnt) return 0;
	    if(angleCnt!=sig.angleCnt) return 0;
	    if(dihedralCnt!=sig.dihedralCnt) return 0;
	    if(improperCnt!=sig.improperCnt) return 0;
	    if(crosstermCnt!=sig.crosstermCnt) return 0;
	    //JLai
	    if(gromacsPairCnt!=sig.gromacsPairCnt) return 0;
	
	#define CMPSIGS(TUPLE) \
	for(int i=0; i<sig.TUPLE##Cnt; i++){ \
	    if(!(TUPLE##Sigs[i]==sig.TUPLE##Sigs[i])) return 0; \
	} \
	
	    CMPSIGS(bond)
	    CMPSIGS(angle)
	    CMPSIGS(dihedral)
	    CMPSIGS(improper)
	    CMPSIGS(crossterm)
	    // JLai
	    CMPSIGS(gromacsPair)
	
	    return 1;
	}
    ~AtomSignature(){
        if(bondSigs) delete[] bondSigs;
        if(angleSigs) delete[] angleSigs;
        if(dihedralSigs) delete[] dihedralSigs;
        if(improperSigs) delete[] improperSigs;
        if(crosstermSigs) delete[] crosstermSigs;
	// JLai
	if(gromacsPairSigs) delete[] gromacsPairSigs;
    }    

    void removeEmptyTupleSigs();
    void pack(MOStream *msg);
    void unpack(MIStream *msg);
};

struct AtomNameIdx{
    Index resnameIdx;
    Index atomnameIdx;
    Index atomtypeIdx;
};
struct AtomCstInfo{
    Index vdw_type;
    int32 status;
    int32 partner;
    int32 hydrogenList;
};

struct ExclusionSignature{
    int fullExclCnt; //1-2, 1-3 exclusion
    int *fullOffset; //should be in increasing order
    int modExclCnt; //1-4 exclusion
    int *modOffset; //should be in increasing order
#ifdef NAMD_CUDA
    int allExclCnt;
    TupleSignature *allTuples;
#endif

    ExclusionSignature(){
    	fullExclCnt = modExclCnt = 0;
    	fullOffset = modOffset = NULL;
#ifdef NAMD_CUDA
    	allExclCnt = 0;
    	allTuples = NULL;
#endif
    }    
    ExclusionSignature(const ExclusionSignature& sig){
        fullOffset = modOffset = NULL;
    	fullExclCnt = sig.fullExclCnt;
        if(fullExclCnt>0){
            fullOffset = new int[fullExclCnt];
            for(int i=0; i<fullExclCnt; i++)
                fullOffset[i] = sig.fullOffset[i];
        }
    	
    	modExclCnt = sig.modExclCnt;
        if(modExclCnt>0){
            modOffset = new int[modExclCnt];
            for(int i=0; i<modExclCnt; i++)
                modOffset[i] = sig.modOffset[i];
        }
#ifdef NAMD_CUDA
    	allTuples = NULL;
    	allExclCnt = sig.allExclCnt;
        if(allExclCnt>0){
            allTuples = new TupleSignature[allExclCnt];
            for(int i=0; i<allExclCnt; i++)
                allTuples[i] = sig.allTuples[i];
        }
#endif
    }
    ~ExclusionSignature(){
    	if(fullOffset) delete [] fullOffset;
    	if(modOffset) delete [] modOffset;
#ifdef NAMD_CUDA
    	if(allTuples) delete [] allTuples;
#endif
    }
    
    ExclusionSignature& operator=(const ExclusionSignature& sig){
        fullExclCnt = sig.fullExclCnt;
        if(fullOffset) delete [] fullOffset;
        if(fullExclCnt>0){
            fullOffset = new int[fullExclCnt];
            for(int i=0; i<fullExclCnt; i++)
                fullOffset[i] = sig.fullOffset[i];
        }else
            fullOffset = NULL;
    
    	modExclCnt = sig.modExclCnt;
        if(modOffset) delete [] modOffset;
        if(modExclCnt>0){
            modOffset = new int[modExclCnt];
            for(int i=0; i<modExclCnt; i++)
                modOffset[i] = sig.modOffset[i];
        }else
            modOffset = NULL;
#ifdef NAMD_CUDA
    	allExclCnt = sig.allExclCnt;
    	if(allTuples) delete [] allTuples;
        if(allExclCnt>0){
            allTuples = new TupleSignature[allExclCnt];
            for(int i=0; i<allExclCnt; i++)
                allTuples[i] = sig.allTuples[i];
        }else
            allTuples = NULL;
#endif

        return *this;
    }
	int operator==(const ExclusionSignature& sig) const{
	    if(fullExclCnt!=sig.fullExclCnt) return 0;
	    if(modExclCnt!=sig.modExclCnt) return 0;
	    
	    for(int i=0; i<fullExclCnt; i++){
			if(fullOffset[i]!=sig.fullOffset[i]) return 0;
	    }
	    for(int i=0; i<modExclCnt; i++){
			if(modOffset[i]!=sig.modOffset[i]) return 0;
	    }
	    return 1;
	}
    //both input should be sorted in increasing order
    void setOffsets(std::vector<int>& fullVec, std::vector<int>& modVec){
    	fullExclCnt = fullVec.size();
    	modExclCnt = modVec.size();
    	if(fullExclCnt>0) {        
            fullOffset = new int[fullExclCnt];
    	    for(int i=0; i<fullExclCnt; i++)
                fullOffset[i] = fullVec[i];
        }

    	if(modExclCnt>0) {        
            modOffset = new int[modExclCnt];
    	    for(int i=0; i<modExclCnt; i++)
                modOffset[i] = modVec[i];	
        }
#ifdef NAMD_CUDA
        buildTuples();
    }
    void buildTuples() {
        delete [] allTuples;
        allTuples = NULL;
    	allExclCnt = 0;
    	for(int i=0; i<fullExclCnt; i++)
            if ( fullOffset[i] > 0 ) ++allExclCnt;
    	for(int i=0; i<modExclCnt; i++)
            if ( modOffset[i] > 0 ) ++allExclCnt;
        if(allExclCnt>0){
            allTuples = new TupleSignature[allExclCnt];
            int j = 0;
            for(int i=0; i<fullExclCnt; i++){
                if ( fullOffset[i] <= 0 ) continue;
                TupleSignature oneSig(1,EXCLUSION,0);
                oneSig.offset[0] = fullOffset[i];
                allTuples[j++] = oneSig;
            }
            for(int i=0; i<modExclCnt; i++){
                if ( modOffset[i] <= 0 ) continue;
                TupleSignature oneSig(1,EXCLUSION,1);
                oneSig.offset[0] = modOffset[i];
                allTuples[j++] = oneSig;
            }
        }
#endif
    }

    int hash() const {
      unsigned int numOffset = fullExclCnt + modExclCnt;
      unsigned int code = 0x12345678;
      unsigned int codesz = 8 * sizeof(int);
      unsigned int shift = codesz / numOffset;
    
      if (shift == 0) shift=1;
      unsigned int i;
      for(i=0; i < fullExclCnt; i++) {
        code = circShift(code,shift);
        code ^= fullOffset[i];
      }
      for(i=0; i < modExclCnt; i++) {
        code = circShift(code,shift);
        code ^= modOffset[i];
      }
      return code;
    }
  
  void removeEmptyOffset();
    //assuming offsets in the signature have been sorted increasingly
    int findOffset(int offset, int  *fullOrMod);
    void pack(MOStream *msg);
    void unpack(MIStream *msg);
};

#endif

