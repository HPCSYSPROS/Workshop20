#include <algorithm>
#include "CompressPsf.h"
#include "strlib.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "InfoStream.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"

using namespace std;

/**
 * Very tricky part: 
 * 
 * Generating the tuple type affects the Parameter object.
 * Particularly, the "multiplicity" field in the TUPLE_array 
 * will be changed in the function assign_TUPLE_index. 
 * TUPLE=dihedral, improper.
 * 
 * And when we read the compressed psf files, assign_TUPLE_index
 * functions will not be called so that the value of "multiplicity"
 * will not be updated. This results in the incorrect energy from bonded
 * interaction.
 * 
 * Therefore, we need to store its value in the compressed psf
 * file. When it comes to read them, we need to update the Parameter object
 * with these values.
 */

//global variables recording compressed psf file information
Molecule *g_mol = NULL;
Parameters *g_param = NULL;
SimParameters *g_simParam = NULL;
ConfigList *g_cfgList = NULL;

struct BasicAtomInfo
{
    Index segNameIdx;
    Index resNameIdx;
    Index atomNameIdx;
    Index atomTypeIdx;
    Index chargeIdx;
    Index massIdx;    
    Index atomSigIdx;
    Index exclSigIdx;
    int resID;
};

void OutputAtomRecord::flip(){
    flipNum((char *)&sSet, sizeof(short), sizeof(sSet)/sizeof(short));
    flipNum((char *)&iSet, sizeof(int), sizeof(iSet)/sizeof(int));
    flipNum((char *)&fSet, sizeof(float), sizeof(fSet)/sizeof(float));
}

struct AtomSigInfo
{
    vector<SigIndex> bondSigIndices;
    vector<SigIndex> angleSigIndices;
    vector<SigIndex> dihedralSigIndices;
    vector<SigIndex> improperSigIndices;
    vector<SigIndex> crosstermSigIndices;

    AtomSigInfo()
    {}
    AtomSigInfo(const AtomSigInfo& sig)
    {
        bondSigIndices.clear();
        for(int i=0; i<sig.bondSigIndices.size(); i++)
            bondSigIndices.push_back(sig.bondSigIndices[i]);

        angleSigIndices.clear();
        for(int i=0; i<sig.angleSigIndices.size(); i++)
            angleSigIndices.push_back(sig.angleSigIndices[i]);

        dihedralSigIndices.clear();
        for(int i=0; i<sig.dihedralSigIndices.size(); i++)
            dihedralSigIndices.push_back(sig.dihedralSigIndices[i]);

        improperSigIndices.clear();
        for(int i=0; i<sig.improperSigIndices.size(); i++)
            improperSigIndices.push_back(sig.improperSigIndices[i]);

        crosstermSigIndices.clear();
        for(int i=0; i<sig.crosstermSigIndices.size(); i++)
            crosstermSigIndices.push_back(sig.crosstermSigIndices[i]);
    }

    ~AtomSigInfo()
    {
        bondSigIndices.clear();
        angleSigIndices.clear();
        dihedralSigIndices.clear();
        improperSigIndices.clear();
        crosstermSigIndices.clear();
    }

    void sortTupleSigIndices()
    {
        sort(bondSigIndices.begin(), bondSigIndices.end());
        sort(angleSigIndices.begin(), angleSigIndices.end());
        sort(dihedralSigIndices.begin(), dihedralSigIndices.end());
        sort(improperSigIndices.begin(), improperSigIndices.end());
        sort(crosstermSigIndices.begin(), crosstermSigIndices.end());
    }
  
    inline CkHashCode hash() const {
      // What's a good hash function for this? Lets make something up
      // Concatenate all the index lists into a list of chars, then hash that
      // string using Charm's string hash function

      // To keep memory allocation cheap, we'll just use a 32-byte buffer
      // and wrap around if we have more sigs
      const int maxlen = 32;
      unsigned char keydata[maxlen+1];
      const int maxchar = 256;
      int i,j;
      for(j=0;j<=maxlen;j++) keydata[j] = 0;
      j=0;
      for(i=0; i<bondSigIndices.size(); i++,j++) {
        keydata[j % maxlen] ^= (bondSigIndices[i] % maxchar);
      }
      for(i=0; i<angleSigIndices.size(); i++,j++) {
        keydata[j % maxlen] ^= (angleSigIndices[i] % maxchar);
      }
      for(i=0; i<dihedralSigIndices.size(); i++,j++) {
        keydata[j % maxlen] ^= (dihedralSigIndices[i] % maxchar);
      }
      for(i=0; i<improperSigIndices.size(); i++,j++) {
        keydata[j % maxlen] ^= (improperSigIndices[i] % maxchar);
      }
      for(i=0; i<crosstermSigIndices.size(); i++,j++) {
        keydata[j % maxlen] ^= (crosstermSigIndices[i] % maxchar);
      }
//      CmiPrintf("Computed hash string len %d,%d\n",j,maxlen);
      if (j > maxlen) j = maxlen;
//      for(i=0; i < j; i++) {
//        if (keydata[i] == 0)
//          keydata[i] = 255;
//        CmiPrintf("key[%d]=%d %p\n",i,keydata[i],keydata);
//      }
      return CkHashFunction_default((const void*)keydata,(size_t)j);
    }
};

int operator==(const AtomSigInfo &s1, const AtomSigInfo& s2)
{
    if(s1.bondSigIndices.size() != s2.bondSigIndices.size())
        return 0;
    if(s1.angleSigIndices.size() != s2.angleSigIndices.size())
        return 0;
    if(s1.dihedralSigIndices.size() != s2.dihedralSigIndices.size())
        return 0;
    if(s1.improperSigIndices.size() != s2.improperSigIndices.size())
        return 0;
    if(s1.crosstermSigIndices.size() != s2.crosstermSigIndices.size())
        return 0;

    int equalCnt;
    equalCnt=0;
    int bondSigCnt = s1.bondSigIndices.size();
    for(int i=0; i<bondSigCnt; i++)
        equalCnt += (s1.bondSigIndices[i]==s2.bondSigIndices[i]);
    if(equalCnt!=bondSigCnt)
        return 0;

    equalCnt=0;
    int angleSigCnt = s1.angleSigIndices.size();
    for(int i=0; i<angleSigCnt; i++)
        equalCnt += (s1.angleSigIndices[i]==s2.angleSigIndices[i]);
    if(equalCnt!=angleSigCnt)
        return 0;

    equalCnt=0;
    int dihedralSigCnt = s1.dihedralSigIndices.size();
    for(int i=0; i<dihedralSigCnt; i++)
        equalCnt += (s1.dihedralSigIndices[i]==s2.dihedralSigIndices[i]);
    if(equalCnt!=dihedralSigCnt)
        return 0;

    equalCnt=0;
    int improperSigCnt = s1.improperSigIndices.size();
    for(int i=0; i<improperSigCnt; i++)
        equalCnt += (s1.improperSigIndices[i]==s2.improperSigIndices[i]);
    if(equalCnt!=improperSigCnt)
        return 0;

    equalCnt=0;
    int crosstermSigCnt = s1.crosstermSigIndices.size();
    for(int i=0; i<crosstermSigCnt; i++)
        equalCnt += (s1.crosstermSigIndices[i]==s2.crosstermSigIndices[i]);
    if(equalCnt!=crosstermSigCnt)
        return 0;

    return 1;
}

struct ExclSigInfo
{
    vector<int> fullExclOffset;
    vector<int> modExclOffset;

    ExclSigInfo()
    {}
    ExclSigInfo(const ExclSigInfo& sig)
    {
        fullExclOffset.clear();
        for(int i=0; i<sig.fullExclOffset.size(); i++)
            fullExclOffset.push_back(sig.fullExclOffset[i]);

        modExclOffset.clear();
        for(int i=0; i<sig.modExclOffset.size(); i++)
            modExclOffset.push_back(sig.modExclOffset[i]);
    }

    ~ExclSigInfo()
    {
        fullExclOffset.clear();
        modExclOffset.clear();
    }

    void sortExclOffset()
    {
        sort(fullExclOffset.begin(), fullExclOffset.end());
        sort(modExclOffset.begin(), modExclOffset.end());
    }

    int hash() const {
      unsigned int code = 0x1234;
      unsigned int codesz = 8 * sizeof(int);
      const unsigned int numFoffset = fullExclOffset.size();
      const unsigned int numMoffset = modExclOffset.size();
      const unsigned int numOffsets = numFoffset + numMoffset;
      
      // No excluded atoms? Just hash to 0
      if (numOffsets == 0)
        return 0;
      
      unsigned int shift = codesz / numOffsets;
      if (shift == 0) shift=1;
      unsigned int i;
      for(i=0; i < numFoffset; i++) {
        code = circShift(code,shift);
        code ^= fullExclOffset[i];
      }
      for(i=0; i < numMoffset; i++) {
        code = circShift(code,shift);
        code ^= modExclOffset[i];
      }
      return code;
    }
};
int operator==(const ExclSigInfo &s1, const ExclSigInfo &s2)
{
    if(s1.fullExclOffset.size()!=s2.fullExclOffset.size())
        return 0;
    if(s1.modExclOffset.size()!=s2.modExclOffset.size())
        return 0;

    for(int i=0; i<s1.fullExclOffset.size(); i++)
    {
        if(s1.fullExclOffset[i] != s2.fullExclOffset[i])
            return 0;
    }

    for(int i=0; i<s1.modExclOffset.size(); i++)
    {
        if(s1.modExclOffset[i] != s2.modExclOffset[i])
            return 0;
    }
    return 1;
}

class HashString : public string {
public:
  int hash() const {
    const char* d = this->c_str();
    int ret=0;
    for (int i=0;d[i]!=0;i++) {
      int shift1=((5*i)%16)+0;
      int shift2=((6*i)%16)+8;
      ret+=((0xa5^d[i])<<shift2)+(d[i]<<shift1);
    }
    return ret;
  }
};

class HashReal {
  Real val;
public:
  HashReal(Real v) : val(v) {}
  operator Real & () { return val; }
  operator const Real & () const { return val; }
  
  int hash() const {
    const char* d = (const char *)&val;
    int ret=0;
    for (int i=0;i < sizeof(Real);i++) {
      int shift1=((5*i)%16)+0;
      int shift2=((6*i)%16)+8;
      ret+=((0xa5^d[i])<<shift2)+(d[i]<<shift1);
    }
    return ret;
  };
};

HashPool<HashString> segNamePool;
HashPool<HashString> resNamePool;
HashPool<HashString> atomNamePool;
HashPool<HashString> atomTypePool;
HashPool<HashReal> chargePool;
HashPool<HashReal> massPool;
HashPool<AtomSigInfo> atomSigPool;
BasicAtomInfo *atomData;

//Recording cluster information after reading all bonds info
int *eachAtomClusterID = NULL;
vector<int> eachClusterSize;
vector<int> eachClusterID;
int g_numClusters = 0;

HashPool<TupleSignature> sigsOfBonds;
HashPool<TupleSignature> sigsOfAngles;
HashPool<TupleSignature> sigsOfDihedrals;
HashPool<TupleSignature> sigsOfImpropers;
HashPool<TupleSignature> sigsOfCrossterms;
AtomSigInfo *eachAtomSigs;

HashPool<ExclSigInfo> sigsOfExclusions;
ExclSigInfo *eachAtomExclSigs;

//Structures for extraBond options
vector<Bond> extraBonds;
vector<Angle> extraAngles;
vector<Dihedral> extraDihedrals;
vector<Improper> extraImpropers;

vector<BondValue> extraBondParams;
vector<AngleValue> extraAngleParams;
vector<DihedralValue> extraDihedralParams;
vector<ImproperValue> extraImproperParams;

int operator==(const BondValue &b1, const BondValue &b2)
{
    return (b1.k==b2.k) && (b1.x0==b2.x0);
}

int operator==(const AngleValue &a1, const AngleValue &a2)
{
    return (a1.k==a2.k) && (a1.k_ub==a2.k_ub) && (a1.r_ub==a2.r_ub) && (a1.theta0==a2.theta0);
}

int operator!=(const FourBodyConsts& f1, const FourBodyConsts& f2)
{
    return (f1.delta!=f2.delta) || (f1.k!=f2.k) || (f1.n!=f2.n);
}

int operator==(const DihedralValue &d1, const DihedralValue &d2)
{
    if(d1.multiplicity != d2.multiplicity)
        return 0;
    for(int i=0; i<MAX_MULTIPLICITY; i++)
    {
        if(d1.values[i] != d2.values[i])
            return 0;
    }
    return 1;
}

int operator==(const ImproperValue &d1, const ImproperValue &d2)
{
    if(d1.multiplicity != d2.multiplicity)
        return 0;
    for(int i=0; i<MAX_MULTIPLICITY; i++)
    {
        if(d1.values[i] != d2.values[i])
            return 0;
    }
    return 1;
}

void loadMolInfo();

void integrateAllAtomSigs();
void outputCompressedFile(FILE *txtOfp, FILE *binOfp);

//reading extraBond's information
void getExtraBonds(StringList *file);

void buildAtomData();
void buildBondData();
void buildAngleData();
void buildDihedralData();
void buildImproperData();
void buildCrosstermData();

void buildExclusionData();

//Functions related with building exclusions
void buildExclusions();
void build12Excls(UniqueSet<Exclusion>&, vector<int> *);
void build13Excls(UniqueSet<Exclusion>&, vector<int> *);
void build14Excls(UniqueSet<Exclusion>&, vector<int> *, int);

//reverse the byte-order of every element starting at "elem"
void flipNum(char *elem, int elemSize, int numElems){
    int mid = elemSize/2;
    char *ptr = elem;
    for(int i=0; i<numElems; i++) {
        for(int j=0; j<mid; j++) {
            char tmp = ptr[j];
            ptr[j] = ptr[elemSize-1-j];
            ptr[elemSize-1-j] = tmp;
        }
        ptr += elemSize;
    }
}

void clearGlobalVectors()
{
    segNamePool.clear();
    resNamePool.clear();
    atomNamePool.clear();
    atomTypePool.clear();
    chargePool.clear();
    massPool.clear();
    sigsOfBonds.clear();
    sigsOfAngles.clear();
    sigsOfDihedrals.clear();
    sigsOfImpropers.clear();
    sigsOfExclusions.clear();

    eachClusterSize.clear();
}

void compress_molecule_info(Molecule *mol, char *psfFileName, Parameters *param, SimParameters *simParam, ConfigList* cfgList)
{
    g_mol = mol;
    g_param = param;
    g_simParam = simParam; //used for building exclusions
    g_cfgList = cfgList; //used for integrating extra bonds

    //read psf files
    //readPsfFile(psfFileName);
    loadMolInfo();

    integrateAllAtomSigs();

    buildExclusions();

    //buildParamData();

    char *outFileName = new char[strlen(psfFileName)+20];
    sprintf(outFileName, "%s.inter", psfFileName);
    //the text file for signatures and other non-per-atom info
    FILE *txtOfp = fopen(outFileName, "w");
    sprintf(outFileName, "%s.inter.bin", psfFileName);
    //the binary file for per-atom info
    FILE *binOfp = fopen(outFileName, "wb");
    delete [] outFileName;

    //output compressed psf file
    outputCompressedFile(txtOfp, binOfp);

    fclose(txtOfp);
    fclose(binOfp);
}

/** Before entering this function, all information about
 *  Molecule has been obtained. The
 *  Molecule::build_atom_status should also be called.
 */
void loadMolInfo()
{
    char err_msg[512];  //  Error message for NAMD_die
    char buffer[512];  //  Buffer for file reading
    int i;      //  Loop counter
    int NumTitle;    //  Number of Title lines in .psf file
    int ret_code;    //  ret_code from NAMD_read_line calls
    FILE *psf_file;
    Parameters *params = g_param;

    buildAtomData();

    //read extra bonds/angles/dihedrals/impropers information first
    //and then integrate them into the following reading procedures
    /*if(g_simParam->extraBondsOn)
    {
        getExtraBonds(g_cfgList->find("extraBondsFile"));
    }*/

    //initialize eachAtomSigs
    eachAtomSigs = new AtomSigInfo[g_mol->numAtoms];    

    buildBondData();
    buildAngleData();
    buildDihedralData();
    buildImproperData();
    buildCrosstermData();

    //  analyze the data and find the status of each atom
    //build_atom_status();
}

void integrateAllAtomSigs()
{
    printf("Bond sigs:  %d\n", (int)sigsOfBonds.size());
    printf("Angle sigs:  %d\n", (int)sigsOfAngles.size());
    printf("Dihedral sigs:  %d\n", (int)sigsOfDihedrals.size());
    printf("Improper sigs:  %d\n", (int)sigsOfImpropers.size());
    printf("Crossterm sigs:  %d\n", (int)sigsOfCrossterms.size());


    for(int i=0; i<g_mol->numAtoms; i++)
    {
        eachAtomSigs[i].sortTupleSigIndices();
        int poolIndex = atomSigPool.lookupCstPool(eachAtomSigs[i]);
        if(poolIndex==-1)
        {
            atomSigPool.push_back(eachAtomSigs[i]);
            poolIndex = atomSigPool.size()-1;
        }
        atomData[i].atomSigIdx = poolIndex;
    }

    printf("Atom's sigs: %d\n", (int)atomSigPool.size());

    delete[] eachAtomSigs;
}

/**
 * Output the compressed psf files. The binary per-atom file 
 * contains two part. The first part is used for the parallel 
 * input containing info such as atom signature ids; the second 
 * part is used for the parallel output containing infor such as 
 * cluster ids and whether an atom is water or not. 
 *  
 * -Chao Mei 
 *  
 */
void outputCompressedFile(FILE *txtOfp, FILE *binOfp)
{
#ifndef MEM_OPT_VERSION
    fprintf(txtOfp, "FORMAT VERSION: %f\n", COMPRESSED_PSF_VER);

    fprintf(txtOfp, "%d !NSEGMENTNAMES\n", segNamePool.size());
    for(int i=0; i<segNamePool.size(); i++)
    {
        fprintf(txtOfp, "%s\n", segNamePool[i].c_str());
    }

    fprintf(txtOfp, "%d !NRESIDUENAMES\n", resNamePool.size());
    for(int i=0; i<resNamePool.size(); i++)
    {
        fprintf(txtOfp, "%s\n", resNamePool[i].c_str());
    }

    fprintf(txtOfp, "%d !NATOMNAMES\n", atomNamePool.size());
    for(int i=0; i<atomNamePool.size(); i++)
    {
        fprintf(txtOfp, "%s\n", atomNamePool[i].c_str());
    }

    fprintf(txtOfp, "%d !NATOMTYPES\n", atomTypePool.size());
    for(int i=0; i<atomTypePool.size(); i++)
    {
        fprintf(txtOfp, "%s\n", atomTypePool[i].c_str());
    }

    fprintf(txtOfp, "%d !NCHARGES\n", chargePool.size());
    for(int i=0; i<chargePool.size(); i++)
    {
        const Real charge = chargePool[i];
        fprintf(txtOfp, "%f\n", charge);
    }

    fprintf(txtOfp, "%d !NMASSES\n", massPool.size());
    for(int i=0; i<massPool.size(); i++)
    {
        const Real mass = massPool[i];
        fprintf(txtOfp, "%f\n", mass);
    }


    fprintf(txtOfp, "%d !NATOMSIGS\n", atomSigPool.size());
    for(int i=0; i<atomSigPool.size(); i++)
    {
        AtomSigInfo& oneAtomSig = atomSigPool[i];
        int oneTypeCnt = oneAtomSig.bondSigIndices.size();
        fprintf(txtOfp, "%d !%sSIGS\n", oneTypeCnt, "NBOND");
        for(int j=0; j<oneTypeCnt; j++)
        {
            SigIndex idx = oneAtomSig.bondSigIndices[j];
            TupleSignature& tSig = sigsOfBonds[idx];
            tSig.output(txtOfp);
        }

        oneTypeCnt = oneAtomSig.angleSigIndices.size();
        fprintf(txtOfp, "%d !%sSIGS\n", oneTypeCnt, "NTHETA");
        for(int j=0; j<oneTypeCnt; j++)
        {
            SigIndex idx = oneAtomSig.angleSigIndices[j];
            TupleSignature& tSig = sigsOfAngles[idx];
            tSig.output(txtOfp);
        }

        oneTypeCnt = oneAtomSig.dihedralSigIndices.size();
        fprintf(txtOfp, "%d !%sSIGS\n", oneTypeCnt, "NPHI");
        for(int j=0; j<oneTypeCnt; j++)
        {
            SigIndex idx = oneAtomSig.dihedralSigIndices[j];
            TupleSignature& tSig = sigsOfDihedrals[idx];
            tSig.output(txtOfp);
        }

        oneTypeCnt = oneAtomSig.improperSigIndices.size();
        fprintf(txtOfp, "%d !%sSIGS\n", oneTypeCnt, "NIMPHI");
        for(int j=0; j<oneTypeCnt; j++)
        {
            SigIndex idx = oneAtomSig.improperSigIndices[j];
            TupleSignature& tSig = sigsOfImpropers[idx];
            tSig.output(txtOfp);
        }

        oneTypeCnt = oneAtomSig.crosstermSigIndices.size();
        fprintf(txtOfp, "%d !%sSIGS\n", oneTypeCnt, "NCRTERM");
        for(int j=0; j<oneTypeCnt; j++)
        {
            SigIndex idx = oneAtomSig.crosstermSigIndices[j];
            TupleSignature& tSig = sigsOfCrossterms[idx];
            tSig.output(txtOfp);
        }
    }

    //2. Output exclusion signatures
    int exclSigCnt = sigsOfExclusions.size();
    fprintf(txtOfp, "%d !NEXCLSIGS\n", exclSigCnt);
    for(int i=0; i<exclSigCnt; i++)
    {
        ExclSigInfo *sig = &sigsOfExclusions[i];
        //first line is for full exclusions (1-2, 1-3) in the format of count offset1 offset2 offset3 ...
        fprintf(txtOfp, "%d", sig->fullExclOffset.size());
        for(int j=0; j<sig->fullExclOffset.size(); j++)
            fprintf(txtOfp, " %d", sig->fullExclOffset[j]);
        fprintf(txtOfp, "\n");

        //second line is for modified exclusions (1-4)
        fprintf(txtOfp, "%d", sig->modExclOffset.size());
        for(int j=0; j<sig->modExclOffset.size(); j++)
            fprintf(txtOfp, " %d", sig->modExclOffset[j]);
        fprintf(txtOfp, "\n");
    }

    //3. Output the cluster information
    fprintf(txtOfp, "%d !NCLUSTERS\n", g_numClusters);

    //4. Output atom info
    fprintf(txtOfp, "%d !NATOM\n", g_mol->numAtoms);
    fprintf(txtOfp, "%d !NHYDROGENGROUP\n", g_mol->numHydrogenGroups);
    fprintf(txtOfp, "%d !MAXHYDROGENGROUPSIZE\n", g_mol->maxHydrogenGroupSize);
    fprintf(txtOfp, "%d !NMIGRATIONGROUP\n", g_mol->numMigrationGroups);
    fprintf(txtOfp, "%d !MAXMIGRATIONGROUPSIZE\n", g_mol->maxMigrationGroupSize);

    //5. Output rigid bond type
    fprintf(txtOfp, "%d !RIGIDBONDTYPE\n", g_simParam->rigidBonds);
#if 0
    const float *atomOccupancy = g_mol->getOccupancyData();
    const float *atomBFactor = g_mol->getBFactorData();
    fprintf(txtOfp, "%d !OCCUPANCYVALID\n", (atomOccupancy==NULL)?0:1);
    fprintf(txtOfp, "%d !TEMPFACTORVALID\n", (atomBFactor==NULL)?0:1);

    float *zeroFloats = NULL;
    if(atomOccupancy==NULL || atomBFactor==NULL) {
        zeroFloats = new float[g_mol->numAtoms];
        memset(zeroFloats, 0, sizeof(float)*g_mol->numAtoms);
        if(atomOccupancy==NULL) atomOccupancy = (const float *)zeroFloats;
        if(atomBFactor==NULL) atomBFactor = (const float *)zeroFloats;
    }
#endif

    Atom *atoms = g_mol->getAtoms(); //need to output its partner and hydrogenList
    HydrogenGroupID *hg = g_mol->hydrogenGroup.begin();

    //First, output magic number
    int magicNum = COMPRESSED_PSF_MAGICNUM;
    fwrite(&magicNum, sizeof(int), 1, binOfp);
    //Second, version number
    float verNum = (float)COMPRESSED_PSF_VER;
    fwrite(&verNum, sizeof(float), 1, binOfp);
    //Third, the per-atom record size
    int recSize = sizeof(OutputAtomRecord);
    fwrite(&recSize, sizeof(int), 1, binOfp);
    //Fourth, each atom info
    OutputAtomRecord oneRec;
    for(int i=0; i<g_mol->numAtoms; i++)
    {                
        oneRec.sSet.segNameIdx = atomData[i].segNameIdx;
        oneRec.sSet.resNameIdx = atomData[i].resNameIdx;
        oneRec.sSet.atomNameIdx = atomData[i].atomNameIdx;
        oneRec.sSet.atomTypeIdx = atomData[i].atomTypeIdx;
        oneRec.sSet.chargeIdx = atomData[i].chargeIdx;
        oneRec.sSet.massIdx = atomData[i].massIdx;
        oneRec.iSet.atomSigIdx = atomData[i].atomSigIdx;

        oneRec.iSet.exclSigIdx = atomData[i].exclSigIdx;        
        oneRec.sSet.vdw_type = atoms[i].vdw_type;
        oneRec.iSet.resID = atomData[i].resID;        
        int hydIdx = atoms[i].hydrogenList;       
        oneRec.iSet.hydrogenList = hydIdx;
        oneRec.iSet.atomsInGroup = hg[hydIdx].atomsInGroup;
        oneRec.iSet.GPID = hg[hydIdx].GPID;
        //oneRec.waterVal = hg[hydIdx].waterVal;
        oneRec.iSet.atomsInMigrationGroup = hg[hydIdx].atomsInMigrationGroup;
        oneRec.iSet.MPID = hg[hydIdx].MPID;
        //oneRec.fSet.occupancy = atomOccupancy[i];
        //oneRec.fSet.bfactor = atomBFactor[i];
        oneRec.fSet.rigidBondLength = g_mol->rigid_bond_length(i);
        fwrite(&oneRec, sizeof(OutputAtomRecord), 1, binOfp);
    }
    //if(zeroFloats) delete[] zeroFloats;
    delete[] atomData;

    //Output the info required for parallel output:
    //1. Cluster ids of each atom
    //Since the cluster info is only when doing output under
    //wrapAll or wrapWater conditions, for the purpose of parallel
    //output, it is reasonable to separate the cluster info from
    //the above per-atom info. So whole the binary per-atom file
    //contains two parts, one for parallel input to create patches
    //and computes; the other for parallel output -Chao Mei
    //Fifth: output the cluster id of each atoms
    fwrite(eachAtomClusterID, sizeof(int), g_mol->numAtoms, binOfp);    
    //2. Whether each atom is water or not
    char *isWater = new char[g_mol->numAtoms];
    for(int i=0; i<g_mol->numAtoms; i++){
        isWater[i] = g_mol->is_water(i);
    }
    fwrite(isWater, sizeof(char), g_mol->numAtoms, binOfp);
    delete [] isWater;
    delete[] atoms;
    g_mol->hydrogenGroup.resize(0);
    delete[] eachAtomClusterID;    

    //Output the parameter new values if extraBonds are present.
    //The parameters are not needed since now extraBonds' parameters will be
    //read again during running the simulation

    //6. Output the "multiplicity" field TUPLE_array of the Parameter object
    fprintf(txtOfp, "!DIHEDRALPARAMARRAY\n");
    for(int i=0; i<g_param->NumDihedralParams; i++)
    {
        fprintf(txtOfp, "%d ", g_param->dihedral_array[i].multiplicity);
    }
    fprintf(txtOfp, "\n");
    fprintf(txtOfp, "!IMPROPERPARAMARRAY\n");
    for(int i=0; i<g_param->NumImproperParams; i++)
    {
        fprintf(txtOfp, "%d ", g_param->improper_array[i].multiplicity);
    }
    fprintf(txtOfp, "\n");
#endif
}

void buildAtomData()
{
#ifndef MEM_OPT_VERSION
    int numAtoms = g_mol->numAtoms;

    //1. parse atom data to build constant pool (atom name, mass, charge etc.)
    atomData = new BasicAtomInfo[numAtoms];
    Atom *atoms = g_mol->getAtoms();
    AtomNameInfo *atomNames = g_mol->getAtomNames();
    AtomSegResInfo *atomSegResids = g_mol->getAtomSegResInfo();

    for(int atomID=0; atomID < numAtoms; atomID++)
    {
        //building constant pool
        int poolIndex;
        HashString fieldName;
        fieldName.assign(atomSegResids[atomID].segname);
        poolIndex = segNamePool.lookupCstPool(fieldName);
        if(poolIndex==-1)
        {
            segNamePool.push_back(fieldName);
            poolIndex = segNamePool.size()-1;
        }
        atomData[atomID].segNameIdx = poolIndex;
        
        atomData[atomID].resID = atomSegResids[atomID].resid;

        fieldName.assign(atomNames[atomID].resname);
        poolIndex = resNamePool.lookupCstPool(fieldName);
        if(poolIndex==-1)
        {
            resNamePool.push_back(fieldName);
            poolIndex = resNamePool.size()-1;
        }
        atomData[atomID].resNameIdx = poolIndex;

        fieldName.assign(atomNames[atomID].atomname);
        poolIndex = atomNamePool.lookupCstPool(fieldName);
        if(poolIndex==-1)
        {
            atomNamePool.push_back(fieldName);
            poolIndex = atomNamePool.size()-1;
        }
        atomData[atomID].atomNameIdx = poolIndex;

        fieldName.assign(atomNames[atomID].atomtype);
        poolIndex = atomTypePool.lookupCstPool(fieldName);
        if(poolIndex==-1)
        {
            atomTypePool.push_back(fieldName);
            poolIndex = atomTypePool.size()-1;
        }
        atomData[atomID].atomTypeIdx = poolIndex;

        poolIndex = chargePool.lookupCstPool(atoms[atomID].charge);
        if(poolIndex==-1)
        {
            chargePool.push_back(atoms[atomID].charge);
            poolIndex = chargePool.size()-1;
        }
        atomData[atomID].chargeIdx = poolIndex;

        poolIndex = massPool.lookupCstPool(atoms[atomID].mass);
        if(poolIndex==-1)
        {
            massPool.push_back(atoms[atomID].mass);
            poolIndex = massPool.size()-1;
        }
        atomData[atomID].massIdx = poolIndex;
    }

    //Free those space to reduce transient memory usage
    //delete [] atoms; (deleted until per-atom info is output)
    delete [] atomNames;
    delete [] atomSegResids;
#endif
}


void buildBondData()
{
#ifndef MEM_OPT_VERSION
    Bond *bonds = g_mol->getAllBonds();    

    //then creating bond's tupleSignature
    for(int i=0; i<g_mol->numBonds; i++)
    {
        Bond *b = bonds+i;
        TupleSignature oneSig(1,BOND,b->bond_type);
        oneSig.offset[0] = b->atom2 - b->atom1;
        oneSig.isReal = (i<g_mol->numRealBonds);

        int poolIndex = sigsOfBonds.lookupCstPool(oneSig);
        int newSig=0;
        if(poolIndex == -1)
        {
            sigsOfBonds.push_back(oneSig);
            poolIndex = (SigIndex)sigsOfBonds.size()-1;
            newSig=1;
        }

      if(!newSig)
        {//check duplicate bonds in the form of (a, b) && (a, b);
          int dupIdx = lookupCstPool(eachAtomSigs[b->atom1].bondSigIndices, (SigIndex)poolIndex);
            if(dupIdx!=-1)
            {
                char err_msg[128];
                sprintf(err_msg, "Duplicate bond %d-%d!", b->atom1+1, b->atom2+1);
                NAMD_die(err_msg);
            }
        }
        eachAtomSigs[b->atom1].bondSigIndices.push_back(poolIndex);
    }

    //check duplicate bonds in the form of (a, b) && (b, a)
    for(int i=0; i<g_mol->numBonds; i++)
    {
        Bond *b=bonds+i;
        int atom2 = b->atom2;
        int thisOffset = atom2 - b->atom1;
        for(int j=0; j<eachAtomSigs[atom2].bondSigIndices.size(); j++)
        {
            SigIndex atom2BondId = eachAtomSigs[atom2].bondSigIndices[j];
            TupleSignature *secSig = &(sigsOfBonds[atom2BondId]);
            if(thisOffset== -(secSig->offset[0]))
            {
                char err_msg[128];
                sprintf(err_msg, "Duplicate bond %d-%d because two atoms are just reversed!", b->atom1+1, atom2+1);
                NAMD_die(err_msg);
            }
        }
    }      

    //building clusters for this simulation system in two steps
    //1. create a list for each atom where each atom in the list is bonded with that atom
    vector<int> *atomListOfBonded = new vector<int>[g_mol->numAtoms];

    for(int i=0; i<g_mol->numRealBonds; i++)
    {
        Bond *b=bonds+i;
        int atom1 = b->atom1;
        int atom2 = b->atom2;
        atomListOfBonded[atom1].push_back(atom2);
        atomListOfBonded[atom2].push_back(atom1);
    }

    delete [] bonds;

    //2. using breadth-first-search to build the clusters. Here, we avoid recursive call
    // because the depth of calls may be of thousands which will blow up the stack, and
    //recursive call is slower than the stack-based BFS.
    //Considering such structure
    //1->1245; 7->1243; 1243->1245
    eachAtomClusterID = new int[g_mol->numAtoms];
    for(int i=0; i<g_mol->numAtoms; i++)
        eachAtomClusterID[i] = -1;

    //It is guaranteed that the clusters found in this way use the
    //smallest atom id of this cluster because each atom is at least
    //connected to one of the atoms in its cluster (by atomListOfBonded
    //constructed above).
    //--Chao Mei
    for(int i=0; i<g_mol->numAtoms; i++)
    {
        int curClusterID=eachAtomClusterID[i];
        //if the atom's cluster id is not -1, the atom has been visited
        if(curClusterID!=-1) continue;

        curClusterID=i;
        deque<int> toVisitAtoms;
	eachAtomClusterID[i] = curClusterID;
        toVisitAtoms.push_back(i);
        while(!toVisitAtoms.empty())
        {
            int visAtomID = toVisitAtoms.front();
            toVisitAtoms.pop_front();
            for(int j=0; j<atomListOfBonded[visAtomID].size(); j++)
            {
                int otherAtom = atomListOfBonded[visAtomID][j];
                if(eachAtomClusterID[otherAtom]!=curClusterID){
                    eachAtomClusterID[otherAtom]=curClusterID;
                    toVisitAtoms.push_back(otherAtom);
		}
            }
        }
    }

#if 0
    //Now the clusterID of each atom should be usually in the non-decreasing
    //order. In other words, the atom ids of a cluster are generally contiguous.
    //If this is the case, the temporary memory usage of output IO during
    //the simulation can be dramatically reduced. So g_isClusterContiguous
    //is used to differentiate the two cases: 

    //1. if the cluster id of atoms is monotonically increasing 
    //(g_isClusterContiguous=1), the size of the cluster can be used as this
    //this cluster's signature (represented by "eachAtomClusterID").
    //The atom who is the first atom (in terms of atom id) in this cluster will 
    //store the cluster size as its signature. The remaining atoms in this 
    //cluster store -1.
    
    //2. if the cluster id of atoms is not monotonically increasing, that is,
    //the atom ids of a cluster are not contiguous (g_isClusterContiguous=0).
    //Then we have to still record each atom's cluster id in variable 
    //"eachAtomClusterID", and we have to change each unique cluster id (valued
    //at the scale of numAtoms) into indexes at the scale of numClusters. For
    //example, the cluster ids of atoms up to this point may be (0...0, 24...24,
    //789...789,...), after re-scaling, cluster ids should be (0...0,1...1,
    //2...2,....) for atoms.
    //We made an ASSUMPTION here: the cluster ids are in the non-decreasing
    //order, and when a cluster discontiguity happens, the cluster id has
    //appeared before! Such ASSUMPTION is to make implementation easy.
    
    int curClusterID;
    int prevClusterID = eachAtomClusterID[0];
    int curClusterSize = 1;
    //step1: segment all atoms according to each cluster
    for(int i=1; i<g_mol->numAtoms; i++){
    	curClusterID = eachAtomClusterID[i];
    	if(curClusterID == prevClusterID){
    		curClusterSize++;
    	}else{
    		eachClusterSize.push_back(curClusterSize);
    		eachClusterID.push_back(prevClusterID);
    		curClusterSize=1;
    	}
    	prevClusterID = curClusterID;	
    }
    //record the info of the last cluster segment
    eachClusterSize.push_back(curClusterSize);
    eachClusterID.push_back(prevClusterID);
    
    //step2: detect contiguity of atoms in each cluster and re-scale
    //cluster id.
    g_isClusterContiguous = 1;    
    int *newClusterIDs = new int[eachClusterID.size()];
    memset(newClusterIDs, 0, sizeof(int)*eachClusterID.size());
    prevClusterID = eachClusterID[0];
    int newCId = 0;    
    newClusterIDs[0] = newCId;
    for(int seg=1; seg<eachClusterID.size(); seg++){
    	curClusterID = eachClusterID[seg];
    	if(curClusterID > prevClusterID){
    		newClusterIDs[seg] = ++newCId;
    		prevClusterID = curClusterID;
    	}else{
    		//non-contiguity happens    		
    		g_isClusterContiguous = 0;
    		
    		//we ASSUME this id appears before
    		//binary search in eachAtomClusterID[0...seg-1).
    		int jl=0, jh=seg-2;
    		int isFound = 0;
    		while(jh>=jl){
    			int mid = (jl+jh)/2;
    			if(curClusterID > eachClusterID[mid])
    				jl = mid+1;
    			else if(curClusterID < eachClusterID[mid])
    				jh = mid-1;
    			else{
    				newClusterIDs[seg] = newClusterIDs[mid];    				
    				isFound = 1;
    				break;
    			}    			
    		}
    		if(!isFound){
    			//ASSUMPTION is wrong and abort
    			char errmsg[300];
    			sprintf(errmsg, "Assumption about building cluster is broken in file %s at line %d\n", __FILE__, __LINE__);
    			NAMD_die(errmsg);
    		}   		    		
    	}    	
    }
    
    //step 3: modify eachAtomClusterID according to g_isClusterContiguous
    //newCId is the id of the last cluster, as id starts from 0, so the
    //total number clusters should be newCId+1
    g_numClusters = newCId+1;
    if(g_isClusterContiguous){
	int aid=0;
        for(int seg=0; seg<eachClusterSize.size(); seg++)
        {
            int curSize = eachClusterSize[seg];
            eachAtomClusterID[aid] = curSize;
            for(int i=aid+1; i<aid+curSize; i++)
                eachAtomClusterID[i] = -1;
            aid += curSize;
        }    	
    }else{
	int aid=0;
        for(int seg=0; seg<eachClusterSize.size(); seg++)
        {
            int curSize = eachClusterSize[seg];            
            for(int i=aid; i<aid+curSize; i++)
                eachAtomClusterID[i] = newClusterIDs[seg];
            aid += curSize;
        }
    }
    free(newClusterIDs);
    eachClusterSize.clear();
    eachClusterID.clear();
#endif 
		
/*
    //check whether cluster is built correctly
    printf("num clusters: %d\n", g_numClusters);
    FILE *checkFile = fopen("cluster.opt", "w");
    for(int i=0; i<g_mol->numAtoms; i++)  fprintf(checkFile, "%d\n", eachAtomClusterID[i]);
    fclose(checkFile);
*/
    
    for(int i=0; i<g_mol->numAtoms; i++)
        atomListOfBonded[i].clear();
    delete [] atomListOfBonded;
#endif
}

void buildAngleData()
{
#ifndef MEM_OPT_VERSION
    Angle *angles = g_mol->getAllAngles();
    //create angles' tupleSignature
    for(int i=0; i<g_mol->numAngles; i++)
    {
        Angle *tuple = angles+i;
        TupleSignature oneSig(2,ANGLE,tuple->angle_type);
        int offset[2];
        offset[0] = tuple->atom2 - tuple->atom1;
        offset[1] = tuple->atom3 - tuple->atom1;
        oneSig.setOffsets(offset);

        int poolIndex = sigsOfAngles.lookupCstPool(oneSig);
        if(poolIndex == -1)
        {
            sigsOfAngles.push_back(oneSig);
            poolIndex = (SigIndex)sigsOfAngles.size()-1;
        }
        eachAtomSigs[tuple->atom1].angleSigIndices.push_back(poolIndex);
    }
    delete [] angles;
#endif
}

void buildDihedralData()
{
#ifndef MEM_OPT_VERSION
    Dihedral *dihedrals = g_mol->getAllDihedrals();    

    //create dihedrals' tupleSignature
    for(int i=0; i<g_mol->numDihedrals; i++)
    {
        Dihedral *tuple = dihedrals+i;
        TupleSignature oneSig(3,DIHEDRAL,tuple->dihedral_type);
        int offset[3];
        offset[0] = tuple->atom2 - tuple->atom1;
        offset[1] = tuple->atom3 - tuple->atom1;
        offset[2] = tuple->atom4 - tuple->atom1;
        oneSig.setOffsets(offset);

        int poolIndex = sigsOfDihedrals.lookupCstPool(oneSig);
        if(poolIndex == -1)
        {
            sigsOfDihedrals.push_back(oneSig);
            poolIndex = (SigIndex)sigsOfDihedrals.size()-1;
        }
        eachAtomSigs[tuple->atom1].dihedralSigIndices.push_back(poolIndex);
    }

    delete[] dihedrals;
#endif
}

void buildImproperData()
{ 
#ifndef MEM_OPT_VERSION
    Improper *impropers=g_mol->getAllImpropers();

    //create improper's tupleSignature
    for(int i=0; i<g_mol->numImpropers; i++)
    {
        Improper *tuple = impropers+i;
        TupleSignature oneSig(3,IMPROPER,tuple->improper_type);
        int offset[3];
        offset[0] = tuple->atom2 - tuple->atom1;
        offset[1] = tuple->atom3 - tuple->atom1;
        offset[2] = tuple->atom4 - tuple->atom1;
        oneSig.setOffsets(offset);

        int poolIndex = sigsOfImpropers.lookupCstPool(oneSig);
        if(poolIndex == -1)
        {
            sigsOfImpropers.push_back(oneSig);
            poolIndex = (SigIndex)sigsOfImpropers.size()-1;
        }
        eachAtomSigs[tuple->atom1].improperSigIndices.push_back(poolIndex);
    }

    delete[] impropers;
#endif
}

void buildCrosstermData()
{
#ifndef MEM_OPT_VERSION
  Crossterm *crossterms = g_mol->getAllCrossterms();
  //create crossterm's tupleSignature
  for(int i=0; i<g_mol->numCrossterms; i++)
  {
    Crossterm *tuple = crossterms+i;
    TupleSignature oneSig(7, CROSSTERM, tuple->crossterm_type);
    int offset[7];
    offset[0] = tuple->atom2 - tuple->atom1;
    offset[1] = tuple->atom3 - tuple->atom1;
    offset[2] = tuple->atom4 - tuple->atom1;
    offset[3] = tuple->atom5 - tuple->atom1;
    offset[4] = tuple->atom6 - tuple->atom1;
    offset[5] = tuple->atom7 - tuple->atom1;
    offset[6] = tuple->atom8 - tuple->atom1;
    oneSig.setOffsets(offset);
   
    int poolIndex = sigsOfCrossterms.lookupCstPool(oneSig);
    if(poolIndex == -1)
    {
      sigsOfCrossterms.push_back(oneSig);
      poolIndex = (SigIndex)sigsOfCrossterms.size()-1;
    }
    eachAtomSigs[tuple->atom1].crosstermSigIndices.push_back(poolIndex);
  }
  
  delete[] crossterms;
#endif
}

void buildExclusions()
{
    //1. Build exclusions: mainly accomplish the function of
    //Molecule::build_exclusions (based on the bonds)
    UniqueSet<Exclusion> allExclusions;

    int exclude_flag; //Exclusion policy
    exclude_flag = g_simParam->exclude;
    //int stripHGroupExclFlag = (simParams->splitPatch == SPLIT_PATCH_HYDROGEN);

    //Commented now since no explicit exclusions are read
    //  Go through the explicit exclusions and add them to the arrays
    //for(i=0; i<numExclusions; i++){
    //	exclusionSet.add(exclusions[i]);
    //}

    // If this is AMBER force field, and readExclusions is TRUE,
    // then all the exclusions were read from parm file, and we
    // shouldn't generate any of them.
    // Comment on stripHGroupExcl:
    // 1. Inside this function, hydrogenGroup is initialized in
    // build_atom_status, therefore, not available when reading psf files
    // 2. this function's main purpose is to reduce memory usage. Since exclusion
    // signatures are used, this function could be overlooked  --Chao Mei

    vector<int> *eachAtomNeighbors = new vector<int>[g_mol->numAtoms];   
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        AtomSigInfo *aSig = &atomSigPool[atomData[atom1].atomSigIdx];
        for(int j=0; j<aSig->bondSigIndices.size(); j++)
        {
            TupleSignature *tSig = &sigsOfBonds[aSig->bondSigIndices[j]];
            if(!tSig->isReal) continue;
            int atom2 = atom1+tSig->offset[0];
            eachAtomNeighbors[atom1].push_back(atom2);
            eachAtomNeighbors[atom2].push_back(atom1);
        }
    }

    if (!g_simParam->amberOn || !g_simParam->readExclusions)
    { //  Now calculate the bonded exlcusions based on the exclusion policy
        switch (exclude_flag)
        {
        case NONE:
            break;
        case ONETWO:
            build12Excls(allExclusions, eachAtomNeighbors);
            break;
        case ONETHREE:
            build12Excls(allExclusions, eachAtomNeighbors);
            build13Excls(allExclusions, eachAtomNeighbors);
            //if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        case ONEFOUR:
            build12Excls(allExclusions, eachAtomNeighbors);
            build13Excls(allExclusions, eachAtomNeighbors);
            build14Excls(allExclusions, eachAtomNeighbors, 0);
            //if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        case SCALED14:
            build12Excls(allExclusions, eachAtomNeighbors);
            build13Excls(allExclusions, eachAtomNeighbors);
            build14Excls(allExclusions, eachAtomNeighbors, 1);
            //if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        }
    }
    //else if (stripHGroupExclFlag && exclude_flag!=NONE && exclude_flag!=ONETWO)
    //  stripHGroupExcl();

    //Commented since atomFepFlags information is not available when reading psf file
    //stripFepExcl();

    for(int i=0; i<g_mol->numAtoms; i++)
        eachAtomNeighbors[i].clear();
    delete [] eachAtomNeighbors;

    //2. Build each atom's list of exclusions
    iout << iINFO << "ADDED " << allExclusions.size() << " IMPLICIT EXCLUSIONS\n" << endi;
    UniqueSetIter<Exclusion> exclIter(allExclusions);
    eachAtomExclSigs = new ExclSigInfo[g_mol->numAtoms];
    for(exclIter=exclIter.begin(); exclIter!=exclIter.end(); exclIter++)
    {
        int atom1 = exclIter->atom1;
        int atom2 = exclIter->atom2;
        int offset21 = atom2-atom1;
        if(exclIter->modified)
        {
            eachAtomExclSigs[atom1].modExclOffset.push_back(offset21);
            eachAtomExclSigs[atom2].modExclOffset.push_back(-offset21);
        }
        else
        {
            eachAtomExclSigs[atom1].fullExclOffset.push_back(offset21);
            eachAtomExclSigs[atom2].fullExclOffset.push_back(-offset21);
        }
    }
    allExclusions.clear();

    //3. Build up exclusion signatures and determine each atom's
    //exclusion signature index
    for(int i=0; i<g_mol->numAtoms; i++)
    {
        eachAtomExclSigs[i].sortExclOffset();
        int poolIndex = sigsOfExclusions.lookupCstPool(eachAtomExclSigs[i]);
        if(poolIndex==-1)
        {
            poolIndex = sigsOfExclusions.size();
            sigsOfExclusions.push_back(eachAtomExclSigs[i]);
        }
        atomData[i].exclSigIdx = poolIndex;
    }
    delete [] eachAtomExclSigs;
    eachAtomExclSigs = NULL;
    printf("Exclusion signatures: %d\n", (int)sigsOfExclusions.size());
}

void build12Excls(UniqueSet<Exclusion>& allExcls, vector<int> *eachAtomNeighbors)
{
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        vector<int> *atom1List = &eachAtomNeighbors[atom1];
        for(int j=0; j<atom1List->size(); j++)
        {
            int atom2 = atom1List->at(j);
            if(atom1<atom2)
                allExcls.add(Exclusion(atom1, atom2));
            else
                allExcls.add(Exclusion(atom2, atom1));
        }
    }
}

void build13Excls(UniqueSet<Exclusion>& allExcls, vector<int> *eachAtomNeighbors)
{
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        vector<int> *atom1List = &eachAtomNeighbors[atom1];
        for(int j=0; j<atom1List->size(); j++)
        {
            int atom2 = atom1List->at(j);
            vector<int> *atom2List = &eachAtomNeighbors[atom2];
            for(int k=0; k<atom2List->size(); k++)
            {
                int atom3 = atom2List->at(k);
                //atom1-atom2, so atom2List contains atom1 which should not be considered
                if(atom3 == atom1)
                    continue;
                if(atom1<atom3)
                    allExcls.add(Exclusion(atom1, atom3));
                else
                    allExcls.add(Exclusion(atom3, atom1));
            }
        }
    }
}

void build14Excls(UniqueSet<Exclusion>& allExcls, vector<int> *eachAtomNeighbors, int modified)
{
    for(int atom1=0; atom1<g_mol->numAtoms; atom1++)
    {
        vector<int> *atom1List = &eachAtomNeighbors[atom1];
        for(int j=0; j<atom1List->size(); j++)
        {
            int atom2 = atom1List->at(j);
            vector<int> *atom2List = &eachAtomNeighbors[atom2];
            for(int k=0; k<atom2List->size(); k++)
            {
                int atom3 = atom2List->at(k);
                //atom1-atom2, so atom2List contains atom1 which should not be considered
                if(atom3 == atom1)
                    continue;
                vector<int> *atom3List = &eachAtomNeighbors[atom3];
                for(int l=0; l<atom3List->size(); l++)
                {
                    int atom4 = atom3List->at(l);
                    //atom1-atom2, so atom2List contains atom1 which should not be considered
                    if(atom4 == atom2 || atom4 == atom1)
                        continue;
                    if(atom1<atom4)
                        allExcls.add(Exclusion(atom1, atom4, modified));
                    else
                        allExcls.add(Exclusion(atom4, atom1, modified));
                }
            }
        }
    }
}

template <class T> void HashPool<T>::dump_tables()
{
  for(int j=0; j < pool.size(); j++) {
    HashPoolAdaptorT<T>* pval = pool[j];
    CmiPrintf("Pool[%d]=%p %p  hash = %d\n",j,pool[j],pval,pval->hash());
  }
  CkHashtableIterator *iter = index_table.iterator();
  void *key,*indx;
  while (iter->hasNext()) {
    indx = iter->next(&key);
    HashPoolAdaptorT<T> *pkey = (HashPoolAdaptorT<T>*)key;
    CmiPrintf("key %p indx %p %d hash=%d\n",key,indx,*((int *)indx),pkey->hash());
  }
}
