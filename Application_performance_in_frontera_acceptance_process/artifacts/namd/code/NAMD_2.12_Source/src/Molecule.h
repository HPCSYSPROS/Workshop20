/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   This class is used to store all of the structural       
   information for a simulation.  It reads in this information 
   from a .psf file, cross checks and obtains some information
   from the Parameters object that is passed in, and then     
   stores all this information for later use. 
*/


#ifndef MOLECULE_H

#define MOLECULE_H

#include "parm.h"

#include "common.h"
#include "NamdTypes.h"
#include "structures.h"
#include "ConfigList.h"
#include "Vector.h"
#include "UniqueSet.h"
#include "Hydrogen.h"
#include "GromacsTopFile.h"
/* BEGIN gf */
#include "GridForceGrid.h"
#include "Tensor.h"
/* END gf */

// Go -- JLai
#define MAX_GO_CHAINS          10
#define MAX_RESTRICTIONS       10
// End of Go

#include "molfile_plugin.h"

#include <vector>

class SimParameters;
class Parameters;
class ConfigList;
class PDB;
class MIStream;
class MOStream;

class ExclElem;
class BondElem;
class AngleElem;
class DihedralElem;
class ImproperElem;
class TholeElem;  // Drude model
class AnisoElem;  // Drude model
class CrosstermElem;
// JLai
class GromacsPairElem;
// End of JLai
class ResidueLookupElem;

struct OutputAtomRecord;

template<class Type> class ObjectArena;

class ExclusionCheck {
public:
  int32 min,max;
  char *flags;

  ExclusionCheck(){
      min=0;
      max=-1;
      flags = NULL;
  }
private:
  ExclusionCheck(const ExclusionCheck& chk){
  }
  ExclusionCheck &operator=(const ExclusionCheck& chk){
    return *this;
  }
};
#define EXCHCK_FULL 1
#define EXCHCK_MOD 2

class ResidueLookupElem
{
public:
  char mySegid[11];
  ResidueLookupElem *next;      // stored as a linked list
  int firstResid;               // valid resid is >= firstResid
  int lastResid;                // valid resid is <= lastResid
  ResizeArray<int> atomIndex;   // 0-based index for first atom in residue

  ResidueLookupElem(void) { next = 0; firstResid = -1; lastResid = -1; }
  ~ResidueLookupElem(void) { delete next; }
  int lookup(const char *segid, int resid, int *begin, int *end) const;
  ResidueLookupElem* append(const char *segid, int resid, int aid);
};

// Ported by JLai -- JE
typedef struct go_val
{
  Real epsilon;      // Epsilon
  int exp_a;         // First exponent for attractive L-J term
  int exp_b;         // Second exponent for attractive L-J term
  int exp_rep;       // Exponent for repulsive L-J term
  Real sigmaRep;     // Sigma for repulsive term
  Real epsilonRep;   // Epsilon for replusive term
  Real cutoff;       // Cutoff distance for Go calculation
  int  restrictions[MAX_RESTRICTIONS];  //  List of residue ID differences to be excluded from Go calculation
} GoValue;

typedef struct go_pair
{
  int goIndxA;       // indexA
  int goIndxB;       // indexB
  double A;          // double A for the LJ pair
  double B;          // double B for the LJ pair
} GoPair;
// End of port - JL

#define QMLSSMODEDIST 1
#define QMLSSMODECOM 2

#define QMFormatORCA 1
#define QMFormatMOPAC 2
#define QMFormatUSR 3

#define QMCHRGNONE 0
#define QMCHRGMULLIKEN 1
#define QMCHRGCHELPG 2

#define QMSCHEMECS 1
#define QMSCHEMERCD 2
#define QMSCHEMEZ1 3
#define QMSCHEMEZ2 4
#define QMSCHEMEZ3 5

//only used for compressing the molecule information
typedef struct seg_resid
{
    char segname[11];
    int resid;
}AtomSegResInfo;

// List maintaining the global atom indicies sorted by helix groups.
class Molecule
{
 private:
typedef struct constraint_params
{
   Real k;    //  Force constant
   Vector refPos; //  Reference position for restraint
} ConstraintParams;



/* BEGIN gf */
typedef struct gridfrc_params
{
    Real k;	// force multiplier
    Charge q;	// charge
} GridforceParams;
/* END gf */


typedef struct stir_params
{
  Real startTheta;
  Vector refPos;  //  Reference position for stirring
} StirParams;

typedef struct movdrag_params
{
  Vector v;            //  Linear velocity (A/step)
} MovDragParams;


typedef struct rotdrag_params
{
   Real v;              //  Angular velocity (deg/step)
   Vector a;            //  Rotation axis
   Vector p;            //  Rotation pivot point
} RotDragParams;

typedef struct constorque_params
{
   Real v;              //  "Torque" value (Kcal/(mol*A^2))
   Vector a;            //  Rotation axis
   Vector p;            //  Rotation pivot point
} ConsTorqueParams;

#ifdef MEM_OPT_VERSION
//indicate a range of atoms from aid1 to aid2
struct AtomSet{
    AtomID aid1;
    AtomID aid2;

    int operator < (const AtomSet &a) const{
		return aid1 < a.aid1;
	}
};
typedef ResizeArray<AtomSet> AtomSetList;

  void load_atom_set(StringList *setfile, const char *setname,
        int *numAtomsInSet, Molecule::AtomSetList **atomsSet) const;

#endif

friend class ExclElem;
friend class BondElem;
friend class AngleElem;
friend class DihedralElem;
friend class ImproperElem;
friend class TholeElem;  // Drude model
friend class AnisoElem;  // Drude model
friend class CrosstermElem;
// JLai
friend class GromacsPairElem;
// End of JLai
friend class WorkDistrib;

private:

#ifndef MEM_OPT_VERSION
	Atom *atoms;    //  Array of atom structures
	ObjectArena<char> *nameArena;
	AtomNameInfo *atomNames;//  Array of atom name info.  Only maintained on node 0 for VMD interface
	//replaced by atom signatures
	Bond *bonds;    //  Array of bond structures
	Angle *angles;    //  Array of angle structures
	Dihedral *dihedrals;  //  Array of dihedral structures
	Improper *impropers;  //  Array of improper structures                          
	Crossterm *crossterms;  //  Array of cross-term structures
        GromacsPair *gromacsPair; //  Array of gromacs-pair structures

	//These will be replaced by exclusion signatures
	Exclusion *exclusions;  //  Array of exclusion structures
	UniqueSet<Exclusion> exclusionSet;  //  Used for building

	int32 *cluster;   //  first atom of connected cluster

	ObjectArena<int32> *tmpArena;
	int32 **bondsWithAtom;  //  List of bonds involving each atom
	ObjectArena<int32> *arena;

	//function is replaced by atom signatures
	int32 **bondsByAtom;  //  List of bonds owned by each atom
	int32 **anglesByAtom;     //  List of angles owned by each atom
	int32 **dihedralsByAtom;  //  List of dihedrals owned by each atom
	int32 **impropersByAtom;  //  List of impropers owned by each atom
	int32 **crosstermsByAtom;  //  List of crossterms owned by each atom

	int32 **exclusionsByAtom; //  List of exclusions owned by each atom
	int32 **fullExclusionsByAtom; //  List of atoms excluded for each atom
	int32 **modExclusionsByAtom; //  List of atoms modified for each atom
// JLai
	int32 **gromacsPairByAtom; // gromacsPair exclusion list which by definition should not have any exclusions (still not sure if it should be empty or zeroed out)
// End of JLai

	ObjectArena<char> *exclArena;
	ExclusionCheck *all_exclusions;

	// DRUDE
	int32 **tholesByAtom;  // list of Thole correction terms owned by each atom
	int32 **anisosByAtom;  // list of anisotropic terms owned by each atom
	// DRUDE

#else
	//Indexing to constant pools to save space
	AtomCstInfo *atoms;
	Index *eachAtomMass; //after initialization, this could be freed (possibly)
	Index *eachAtomCharge; //after initialization, this could be freed (possibly)
	AtomNameIdx *atomNames;
	ObjectArena<char> *nameArena; //the space for names to be allocated  

	//A generally true assumption: all atoms are arranged in the order of clusters. In other words,
	//all atoms for a cluster must appear before/after any atoms in other clusters
	//The first atom in the cluster (which has the lowest atom id) stores the cluster size
	//other atoms in the cluster stores -1
	int32 *clusterSigs;
	int32 numClusters;
	
	AtomSigID *eachAtomSig;
	ExclSigID *eachAtomExclSig;

    AtomSetList *fixedAtomsSet;    
    AtomSetList *constrainedAtomsSet;    
#endif
 
  ResidueLookupElem *resLookup; // find residues by name

  AtomSegResInfo *atomSegResids; //only used for generating compressed molecule info

  Bond *donors;         //  Array of hydrogen bond donor structures
  Bond *acceptors;  //  Array of hydrogen bond acceptor

  // DRUDE
  DrudeConst *drudeConsts;  // supplement Atom data (length of Atom array)
  Thole *tholes;            // Thole (screened Coulomb) interactions
  Aniso *anisos;            // anisotropic terms
  Lphost *lphosts;          // lone pair hosts
  int32 *lphostIndexes;     // index for each LP into lphosts array
  // DRUDE

  int32 *consIndexes; //  Constraint indexes for each atom
  ConstraintParams *consParams;

/* BEGIN gf */
  int32 **gridfrcIndexes;
  GridforceParams **gridfrcParams;
  GridforceGrid **gridfrcGrid;
/* END gf */

        //  Parameters for each atom constrained
  int32 *stirIndexes; //Stirring indexes for each atoms
  StirParams *stirParams;
                          // Paramters for each atoms stirred
  int32 *movDragIndexes;  //  Moving drag indexes for each atom
  MovDragParams *movDragParams;
                                //  Parameters for each atom moving-dragged
  int32 *rotDragIndexes;  //  Rotating drag indexes for each atom
  RotDragParams *rotDragParams;
                                //  Parameters for each atom rotation-dragged

  Real *langevinParams;   //  b values for langevin dynamics
  int32 *fixedAtomFlags;  //  1 for fixed, -1 for fixed group, else 0
  int32 *exPressureAtomFlags; // 1 for excluded, -1 for excluded group.

  //In the memory optimized version: it will be NULL if the general
  //true assumption mentioned above holds. If not, its size is numClusters.
  //In the ordinary version: its size is numAtoms, and indicates the size
  //of connected cluster or 0.
  int32 *clusterSize;                            

	Real *rigidBondLengths;  //  if H, length to parent or 0. or
	//  if not H, length between children or 0.
//fepb
	unsigned char *fepAtomFlags; 
//fepe

  //occupancy and bfactor data from plugin-based IO implementation of loading structures
  float *occupancy;
  float *bfactor;


  // added during the trasition from 1x to 2
  SimParameters *simParams;
  Parameters *params;

private:
	void initialize(SimParameters *, Parameters *param);
	// Sets most data to zero

  //LCPO
  int *lcpoParamType;

#ifndef MEM_OPT_VERSION
	void read_psf_file(char *, Parameters *);
        //  Read in a .psf file given
        //  the filename and the parameter
        //  object to use          
  
  void read_atoms(FILE *, Parameters *);
        //  Read in atom info from .psf
  void read_bonds(FILE *, Parameters *);
        //  Read in bond info from .psf
  void read_angles(FILE *, Parameters *);
        //  Read in angle info from .psf
  void read_dihedrals(FILE *, Parameters *);
        //  Read in dihedral info from .psf
  void read_impropers(FILE *, Parameters *);
        //  Read in improper info from .psf
  void read_crossterms(FILE *, Parameters *);
        //  Read in cross-term info from .psf
  void read_donors(FILE *);
        //  Read in hydrogen bond donors from .psf
  void read_acceptors(FILE *);
        //  Read in hydrogen bond acceptors from .psf
  void read_exclusions(FILE *);
        //  Read in exclusion info from .psf
  // JLai August 16th, 2012 Modifications
  void read_exclusions(int*, int*, int);
  /* Read in exclusion array and sort entries */
  static bool goPairCompare (GoPair, GoPair);
  // End of JLai August 16th, 2012 Modifications


  // DRUDE: PSF reading
  void read_lphosts(FILE *);
        //  Read in lone pair hosts from Drude PSF
  void read_anisos(FILE *);
        //  Read in anisotropic terms from Drude PSF
  // DRUDE

  //LCPO
  //input type is Charmm/Amber/other
  //0 - Charmm/Xplor
  //1 - Amber TODO
  //2 - Plugin TODO
  //3 - Gromacs TODO
  void assignLCPOTypes(int inputType);

  //pluginIO-based loading atoms' structure
  void plgLoadAtomBasics(molfile_atom_t *atomarray);
  void plgLoadBonds(int *from, int *to); //atom index is 1-based in the parameters
  void plgLoadAngles(int *plgAngles);
  void plgLoadDihedrals(int *plgDihedrals);
  void plgLoadImpropers(int *plgImpropers);
  void plgLoadCrossterms(int *plgCterms);

	//  List of all exclusions, including
	//  explicit exclusions and those calculated
	//  from the bonded structure based on the
	//  exclusion policy.  Also includes 1-4's.
  void build_lists_by_atom();
	//  Build the list of structures by atom

  void build12excl(void);
  void build13excl(void);
  void build14excl(int);

  // DRUDE: extend exclusions for Drude and LP
  void build_inherited_excl(int);
  // DRUDE
  void stripFepExcl(void);

  void build_exclusions();
  // analyze the atoms, and determine which are oxygen, hb donors, etc.
  // this is called after a molecule is sent our (or received in)
  void build_atom_status(void);

#else
  //the method to load the signatures of atoms etc. (i.e. reading the file in 
  //text fomrat of the compressed psf file)
  void read_mol_signatures(char *fname, Parameters *params, ConfigList *cfgList=0);	
  void load_one_inputatom(int aid, OutputAtomRecord *one, InputAtom *fAtom);
  void build_excl_check_signatures();  
#endif

	void read_parm(const GromacsTopFile *);  
	  
public:
  // DRUDE
  int is_drude_psf;      // flag for reading Drude PSF
  int is_lonepairs_psf;  // flag for reading lone pairs from PSF
  // DRUDE

  // data for TIP4P
  Real r_om;
  Real r_ohc;

  // data for tail corrections
  BigReal tail_corr_ener;
  BigReal tail_corr_virial;
  BigReal tail_corr_dUdl_1;
  BigReal tail_corr_virial_1;
  void compute_LJcorrection();
  BigReal getEnergyTailCorr(const BigReal);
  BigReal getVirialTailCorr(const BigReal);

  int const * getLcpoParamType() {
    return lcpoParamType;
  }

  BigReal GetAtomAlpha(int i) const {
    return drudeConsts[i].alpha;
  }

#ifdef MEM_OPT_VERSION
  AtomCstInfo *getAtoms() const { return atoms; }
  AtomNameIdx *getAtomNames() const { return atomNames; }
#else
  Atom *getAtoms () const { return atoms; }
  AtomNameInfo *getAtomNames() const { return atomNames; }
#endif

  AtomSegResInfo *getAtomSegResInfo() const { return atomSegResids; }
  
  // return number of fixed atoms, guarded by SimParameters
  int num_fixed_atoms() const {
    // local variables prefixed by s_
    int s_NumFixedAtoms = (simParams->fixedAtomsOn ? numFixedAtoms : 0);
    return s_NumFixedAtoms;  // value is "turned on" SimParameters
  }

  int num_fixed_groups() const {
    // local variables prefixed by s_
    int s_NumFixedAtoms = num_fixed_atoms();
    int s_NumFixedGroups = (s_NumFixedAtoms ? numFixedGroups : 0);
    return s_NumFixedGroups;  // value is "turned on" SimParameters
  }

  int num_group_deg_freedom() const {
    // local variables prefixed by s_
    int s_NumGroupDegFreedom = 3 * numHydrogenGroups;
    int s_NumFixedAtoms = num_fixed_atoms();
    int s_NumFixedGroups = num_fixed_groups();
    if (s_NumFixedGroups) s_NumGroupDegFreedom -= 3 * s_NumFixedGroups;
    if ( ! (s_NumFixedAtoms || numConstraints
          || simParams->comMove || simParams->langevinOn) ) {
      s_NumGroupDegFreedom -= 3;
    }
    return s_NumGroupDegFreedom;
  }

  int num_deg_freedom(int isInitialReport = 0) const {
    // local variables prefixed by s_
    int s_NumDegFreedom = 3 * numAtoms;
    int s_NumFixedAtoms = num_fixed_atoms();
    if (s_NumFixedAtoms) s_NumDegFreedom -= 3 * s_NumFixedAtoms;
    if (numLonepairs) s_NumDegFreedom -= 3 * numLonepairs;
    if ( ! (s_NumFixedAtoms || numConstraints
          || simParams->comMove || simParams->langevinOn) ) {
      s_NumDegFreedom -= 3;
    }
    if ( ! isInitialReport && simParams->pairInteractionOn) {
      //
      // DJH: a kludge?  We want to initially report system degrees of freedom
      //
      // this doesn't attempt to deal with fixed atoms or constraints
      s_NumDegFreedom = 3 * numFepInitial;
    }
    int s_NumFixedRigidBonds = 
      (simParams->fixedAtomsOn ? numFixedRigidBonds : 0);
    if (simParams->watmodel == WAT_TIP4) {
      // numLonepairs is subtracted here because all lonepairs have a rigid bond
      // to oxygen, but all of the LP degrees of freedom are dealt with above
      s_NumDegFreedom -= (numRigidBonds - s_NumFixedRigidBonds - numLonepairs);
    }
    else {
      // Non-polarized systems don't have LPs.
      // For Drude model, bonds that attach LPs are not counted as rigid;
      // LPs have already been subtracted from degrees of freedom.
      s_NumDegFreedom -= (numRigidBonds - s_NumFixedRigidBonds);
    }
    return s_NumDegFreedom;
  }

  int numAtoms;   //  Number of atoms                   

  int numRealBonds;   //  Number of bonds for exclusion determination
  int numBonds;   //  Number of bonds calculated, including extras
  int numAngles;    //  Number of angles
  int numDihedrals; //  Number of dihedrals
  int suspiciousAlchBonds;    //  angles dropped due to crossing FEP partitions
  int alchDroppedAngles;    //  angles dropped due to crossing FEP partitions
  int alchDroppedDihedrals; //  dihedrals dropped due to crossing FEP partitions
  int alchDroppedImpropers; //  impropers dropped due to crossing FEP partitions
  int numImpropers; //  Number of impropers
  int numCrossterms; //  Number of cross-terms
  int numDonors;          //  Number of hydrogen bond donors
  int numAcceptors; //  Number of hydrogen bond acceptors
  int numExclusions;  //  Number of exclusions
  int numTotalExclusions; //  Real Total Number of Exclusions // hack

  // DRUDE
  int numLonepairs; // Number of lone pairs
  int numDrudeAtoms;  // Number of Drude particles
  int numTholes;  // Number of Thole terms
  int numAnisos;  // Number of anisotropic terms
  int numLphosts;  // Number of lone pair hosts
  // DRUDE
  
  int numConstraints; //  Number of atoms constrained
/* BEGIN gf */
  int numGridforceGrids;//  Number of gridforce grids
  int *numGridforces;	//  Number of atoms in gridforce file (array, one per grid)
/* END gf */
  int numMovDrag;         //  Number of atoms moving-dragged
  int numRotDrag;         //  Number of atoms rotating-dragged
  int numConsTorque;  //  Number of atoms "constant"-torqued
  int numFixedAtoms;  //  Number of fixed atoms
  int numStirredAtoms;  //  Number of stirred atoms
  int numExPressureAtoms; //  Number of atoms excluded from pressure
  int numHydrogenGroups;  //  Number of hydrogen groups
  int maxHydrogenGroupSize;  //  Max atoms per hydrogen group
  int numMigrationGroups;  //  Number of migration groups
  int maxMigrationGroupSize;  //  Max atoms per migration group
  int numFixedGroups; //  Number of totally fixed hydrogen groups
  int numRigidBonds;  //  Number of rigid bonds
  int numFixedRigidBonds; //  Number of rigid bonds between fixed atoms
//fepb
  int numFepInitial;  // no. of fep atoms with initial flag
  int numFepFinal;  // no. of fep atoms with final flag
//fepe

  int numConsForce; //  Number of atoms that have constant force applied
  int32 *consForceIndexes;//  Constant force indexes for each atom
  Vector *consForce;  //  Constant force array

  int32 *consTorqueIndexes; //  "Constant" torque indexes for each atom
  ConsTorqueParams *consTorqueParams;
                                //  Parameters for each atom "constant"-torqued

  // The following are needed for error checking because we
  // eliminate bonds, etc. which involve only fixed atoms
  int numCalcBonds; //  Number of bonds requiring calculation
  int numCalcAngles;  //  Number of angles requiring calculation
  int numCalcDihedrals; //  Number of dihedrals requiring calculation
  int numCalcImpropers; //  Number of impropers requiring calculation
  int numCalcCrossterms; //  Number of cross-terms requiring calculation
  int numCalcExclusions;  //  Number of exclusions requiring calculation
  int numCalcFullExclusions;  //  Number of full exclusions requiring calculation

  // DRUDE
  int numCalcTholes;  // Number of Thole correction terms requiring calculation
  int numCalcAnisos;  // Number of anisotropic terms requiring calculation
  // DRUDE

  //  Number of dihedrals with multiple periodicity
  int numMultipleDihedrals; 
  //  Number of impropers with multiple periodicity
  int numMultipleImpropers; 
  // indexes of "atoms" sorted by hydrogen groups
  HydrogenGroup hydrogenGroup;

  // Ported by JLai -- JE - Go
  int numGoAtoms;         //  Number of atoms subject to Go forces -- ported by JLai/ Original by JE
  int32 *atomChainTypes;  //  Go chain type for each atom; from 1 to MAX_GO_CHAINS
  int32 *goSigmaIndices;  //  Indices into goSigmas
  int32 *goResidIndices;  //  Indices into goSigmas
  Real  *goSigmas;        //  Sigma values for Go forces L-J type formula
  bool *goWithinCutoff;   //  Whether the reference atom-atom distance is within the Go cutoff
  Real  *goCoordinates;   //  Coordinates (x,y,z) for Go atoms in the native structure
  int *goResids;          //  Residue ID from PDB
  PDB *goPDB;             //  Pointer to PDB object to use
  // NAMD-Go2 calculation code
  int goNumLJPair;        //  Integer storing the total number of explicit pairs (LJ)
  int *goIndxLJA;         //  Pointer to the array of atom indices for LJ atom A
  int *goIndxLJB;         //  Pointer to the array of atom indices for LJ atom B
  double *goSigmaPairA;  //  Pointer to the array of A LJ parameters
  double *goSigmaPairB;  //  Pointer to the array of B LJ parameters
  int *pointerToGoBeg;    //  Array of pointers to array
  int *pointerToGoEnd;    //  Array of pointers to array
  // Gromacs LJ Pair list calculation code
  int numPair;            //  Integer storing the total number of explicit pairs (LJ + Gaussian)
  int numLJPair;          //  Integer storing the number of explicit LJ pairs
  int numCalcLJPair;         //  Number of explicit LJ pairs requiring calculation
  int *pointerToLJBeg;       //  Array of pointers to array 
  int *pointerToLJEnd;       //  Array of pointers to array B
  int *indxLJA;           //  Pointer to the array of atom indices for LJ atom A
  int *indxLJB;           //  Pointer to the array of atom indices for LJ atom B
  Real *pairC6;           //  Pointer to the array of C6 LJ parameters
  Real *pairC12;          //  Pointer to the array of C12 LJ parameters
  int *gromacsPair_type;   //  
  // Gromacs Gauss Pair list calculation code
  int *pointerToGaussBeg;    //  Array of pointers to array B
  int *pointerToGaussEnd;    //  Array of pointers to array B
  int numGaussPair;       //  Integer storing the number of explicit Gaussian pairs  
  int *indxGaussA;        //  Pointer to the array of atom indices for Gaussian atom A 
  int *indxGaussB;        //  Pointer to the array of atom indices for Gaussian atom B 
  Real *gA;               //  Pointer to the array of force constants to the Gaussian potential
  Real *gMu1;             //  Pointer to the array of mu (shifts Gaussian)
  Real *giSigma1;          //  Pointer to the array of sigma (controls spread of Gaussian)
  Real *gMu2;             //  Pointer to the array of mu (shifts Gaussian 2)
  Real *giSigma2;          //  Pointer to the array of sigma (controls spread of Gaussian 2)
  Real *gRepulsive;       //  Pointer to the a LJ-12 repulsive parameter that adds to the Gaussian

  // GO ENERGY CALCULATION CODE
  BigReal energyNative;    // Holds the energy value of the native structure
  BigReal energyNonnative; // Holds the energy value of the nonnative structure
  // GO ENERGY CALCULATION CODE
  // End of port - JL

  
private:
  /////////////////////////
  // Begin QM
  /////////////////////////
    
  int qmDroppedBonds, qmDroppedAngles, qmDroppedDihedrals;
  int qmDroppedImpropers, qmDroppedCrossterms;
    
  Bool qmReplaceAll;              // Indicates if all forces should be replaced.
  int qmNumQMAtoms;           // Number of QM atoms, from all groups.
  
  // Array of abbreviated element names for all atoms.
  char **qmElementArray;
  // Array of elements for dummy atoms.
  char **qmDummyElement;
  
  // Number of QM atoms in each QM group
  int *qmGrpSizes;
  
  // QM group for each atom of the system. 0 means MM atom.
  Real *qmAtomGroup;
  // Number of different QM regions
  int qmNumGrps ;
  
  // QM Atom charges.
  Real *qmAtmChrg ;
  // global indices of all QM atoms.
  int *qmAtmIndx ;
  
  // Number of QM-MM bonds.
  int qmNumBonds ;
  // IDs of each group, that is, the value in the qmColumn, per group.
  Real *qmGrpID ;
  // The total charge of each QM group.
  Real *qmGrpChrg;
  // The multiplicity of QM groups
  Real *qmGrpMult;
  // Number of QM-MM bonds in each group.
  int *qmGrpNumBonds ;
  // Sequential list of bonds between a MM atom and a QM atom.
  // (with the bonded atoms ins this order: MM first, QM second).
  int **qmMMBond ;
  // List of bond indexes for each QM group.
  int **qmGrpBonds ;
  // List of atom indexes for MM atoms in QM-MM bonds, per group.
  // This is used to check if a point charge is actually an MM atom
  // that need to be ignored, and will be substituted by a dummy atom (for example).
  int **qmMMBondedIndx ;
  // List of indices for MM atoms which will receive fractions of charges from 
  // the MM atom of a QM-MM bond. This is used in the Charge Shift scheme.
  int **qmMMChargeTarget;
  // Number of target MM atoms per QM-MM bond. Targets will receive the partial
  // charge of the MM atom in a QM-MM bond. (in Charge shift scheme)
  int *qmMMNumTargs;
  // List of distances from thr QM atom where the dummy atom will be placed.
  BigReal *qmDummyBondVal;
  // Indicate if point charges will be used in QM calculations.
  Bool qmNoPC;
  
  /////////////////////////
  // These variables are used in case mechanichal embedding is requested (NoPC = TRUE)
  // AND there are QM-MM bonds in the system.
  
  // Indicates the number of items in the arrays for Mechanical embedding with QM-MM bonds.
  int qmMeNumBonds;
  // Indicates the indices of MM atoms which participate in QM-MM bonds.
  int *qmMeMMindx;
  // Indicates the QM group to which the QM-MM bonds belongs.
  Real *qmMeQMGrp;
  /////////////////////////
  
  // Indicates frequency of selection of point charges.
  int qmPCFreq;
  // Custom PC array;
  int *qmCustomPCIdxs;
  // Number of Indexes associated with each QM Group.
  int *qmCustPCSizes;
  // Total number of custom PC charges.
  int qmTotCustPCs;
  // Number of solvent molecules that will be selected and updated by NAMD, per group.
  int *qmLSSSize;
  // Number of atoms per solvent molecule. We need all molecules to have the same size.
  int qmLSSResidueSize;
  // IDs of all atoms belonging to QM solvent molecules, from all groups.
  int *qmLSSIdxs;
  // MAss of each atom in LSS solvent molecules.
  Mass *qmLSSMass;
  // Atom IDs of QM atoms which will be used to select solvent molecules by distance.
  int *qmLSSRefIDs;
  // Mass of each QM atom used as reference for solvent selection.
  Mass *qmLSSRefMass ;
  // Number of atom IDs/Mass per group.
  int *qmLSSRefSize;
  int qmLSSFreq;
  int qmLSSTotalNumAtms;
  // This will map each classical solvent atom with a global index of solvent residues,
  // so we don't depend on segment names and repeated residue indices.
  std::map<int,int> qmClassicSolv ;
  
  /////////////////////////
  // end QM
  /////////////////////////
  
public:
  // QM
  void set_qm_replaceAll(Bool newReplaceAll) { qmReplaceAll = newReplaceAll; } ;
  
  const Real * get_qmAtomGroup() const {return qmAtomGroup; } ;
  Real get_qmAtomGroup(int indx) const {return qmAtomGroup[indx]; } ;
  
  Real *get_qmAtmChrg() {return qmAtmChrg; } ;
  const int *get_qmAtmIndx() {return qmAtmIndx; } ;
  
  int get_numQMAtoms() {return qmNumQMAtoms; } ;
  const char *const * get_qmElements() {return qmElementArray;} ;
  int get_qmNumGrps() const {return qmNumGrps; } ;
  const int * get_qmGrpSizes() {return qmGrpSizes; } ;
  Real * get_qmGrpID() { return qmGrpID; } ;
  Real *get_qmGrpChrg() {return qmGrpChrg; } ;
  Real *get_qmGrpMult() {return qmGrpMult; } ;
  
  int get_qmNumBonds() { return qmNumBonds; } ;
  const BigReal * get_qmDummyBondVal() { return qmDummyBondVal; } ;
  const int * get_qmGrpNumBonds() { return qmGrpNumBonds; } ;
  const int *const * get_qmMMBond() { return qmMMBond; } ;
  const int *const * get_qmGrpBonds() { return qmGrpBonds; } ;
  const int *const * get_qmMMBondedIndx() { return qmMMBondedIndx; } ;
  const char *const * get_qmDummyElement() { return qmDummyElement; } ;
  
  const int *const * get_qmMMChargeTarget() { return qmMMChargeTarget; } ;
  const int * get_qmMMNumTargs() { return qmMMNumTargs; } ;
  
  int* get_qmLSSSize() { return qmLSSSize;} ;
  int* get_qmLSSIdxs() { return qmLSSIdxs;} ;
  Mass *get_qmLSSMass() { return qmLSSMass; } ;
  int get_qmLSSFreq() { return qmLSSFreq; } ;
  int get_qmLSSResSize() { return qmLSSResidueSize; } ;
  int *get_qmLSSRefIDs() { return qmLSSRefIDs ; } ;
  int *get_qmLSSRefSize() { return qmLSSRefSize ; } ;
  Mass *get_qmLSSRefMass() {return qmLSSRefMass; } ;
  std::map<int,int> &get_qmMMSolv() { return qmClassicSolv;} ;
  
  Bool get_qmReplaceAll() {return qmReplaceAll; } ;
  
  Bool get_noPC() { return qmNoPC; } ;
  int get_qmMeNumBonds() { return qmMeNumBonds; } ;
  int *get_qmMeMMindx() { return qmMeMMindx; } ;
  Real *get_qmMeQMGrp() { return qmMeQMGrp; } ;
  
  int get_qmPCFreq() { return qmPCFreq; } ;
  
  int get_qmTotCustPCs() { return qmTotCustPCs; } ;
  int *get_qmCustPCSizes()  { return qmCustPCSizes; } ;
  int *get_qmCustomPCIdxs() { return qmCustomPCIdxs; } ;
  
  void prepare_qm(const char *pdbFileName, Parameters *params, ConfigList *cfgList) ;
  void delete_qm_bonded() ;
  // end QM
  
  
  
  Molecule(SimParameters *, Parameters *param);
  Molecule(SimParameters *, Parameters *param, char *filename, ConfigList *cfgList=NULL);  
  Molecule(SimParameters *simParams, Parameters *param, molfile_plugin_t *pIOHdl, void *pIOFileHdl, int natoms);
  
  Molecule(SimParameters *, Parameters *, Ambertoppar *);
  void read_parm(Ambertoppar *);

  Molecule(SimParameters *, Parameters *, const GromacsTopFile *);

  ~Molecule();    //  Destructor

  void send_Molecule(MOStream *);
        //  send the molecular structure 
        //  from the master to the clients

  void receive_Molecule(MIStream *);
        //  receive the molecular structure
        //  from the master on a client
  
  void build_constraint_params(StringList *, StringList *, StringList *,
             PDB *, char *);
        //  Build the set of harmonic constraint 
        // parameters

/* BEGIN gf */
  void build_gridforce_params(StringList *, StringList *, StringList *, StringList *, PDB *, char *);
        //  Build the set of gridForce-style force pars
/* END gf */

  void build_movdrag_params(StringList *, StringList *, StringList *, 
          PDB *, char *);
        //  Build the set of moving drag pars

  void build_rotdrag_params(StringList *, StringList *, StringList *,
          StringList *, StringList *, StringList *,
          PDB *, char *);
        //  Build the set of rotating drag pars

  void build_constorque_params(StringList *, StringList *, StringList *,
             StringList *, StringList *, StringList *,
             PDB *, char *);
        //  Build the set of "constant" torque pars


  void build_constant_forces(char *);
        //  Build the set of constant forces

  void build_langevin_params(BigReal coupling, BigReal drudeCoupling,
      Bool doHydrogen);
  void build_langevin_params(StringList *, StringList *, PDB *, char *);
        //  Build the set of langevin dynamics parameters

#ifdef MEM_OPT_VERSION
  void load_fixed_atoms(StringList *fixedFile);
  void load_constrained_atoms(StringList *constrainedFile);
#endif

  void build_fixed_atoms(StringList *, StringList *, PDB *, char *);
        //  Determine which atoms are fixed (if any)

  void build_stirred_atoms(StringList *, StringList *, PDB *, char *);
        //  Determine which atoms are stirred (if any)

  void build_extra_bonds(Parameters *parameters, StringList *file);

//fepb
        void build_fep_flags(StringList *, StringList *, PDB *, char *, const char *);
                               // selection of the mutant atoms
        void delete_alch_bonded(void);
//fepe

  void build_exPressure_atoms(StringList *, StringList *, PDB *, char *);
        //  Determine which atoms are excluded from
                                //  pressure (if any)

  // Ported by JLai -- Original JE - Go -- Change the unsigned int to ints
  void print_go_sigmas(); //  Print out Go sigma parameters
  void build_go_sigmas(StringList *, char *);
        //  Determine which atoms have Go forces applied
        //  calculate sigmas from distances between Go atom pairs
  void build_go_sigmas2(StringList *, char *);
        //  Determine which atoms have Go forces applied
        //  calculate sigmas from distances between Go atom pairs
  void build_go_arrays(StringList *, char *);
        //  Determine which atoms have Go forces applied
  BigReal get_gro_force(BigReal, BigReal, BigReal, int, int) const;
  BigReal get_gro_force2(BigReal, BigReal, BigReal, int, int, BigReal *, BigReal *) const;
  BigReal get_go_force(BigReal, int, int, BigReal *, BigReal *) const;
        //  Calculate the go force between a pair of atoms -- Modified to 
        //  output Go energies
  BigReal get_go_force_new(BigReal, int, int, BigReal *, BigReal *) const;
        //  Calculate the go force between a pair of atoms
  BigReal get_go_force2(BigReal, BigReal, BigReal, int, int, BigReal *, BigReal *) const;
        //  Calculate the go force between a pair of atoms
  Bool atoms_1to4(unsigned int, unsigned int);
// End of port -- JL  

  void reloadCharges(float charge[], int n);

        Bool is_lp(int);     // return true if atom is a lone pair
        Bool is_drude(int);     // return true if atom is a Drude particle
        Bool is_hydrogen(int);     // return true if atom is hydrogen
        Bool is_oxygen(int);       // return true if atom is oxygen
  Bool is_hydrogenGroupParent(int); // return true if atom is group parent
  Bool is_water(int);        // return true if atom is part of water 
  int  get_groupSize(int);     // return # atoms in (hydrogen) group
        int get_mother_atom(int);  // return mother atom of a hydrogen

  #ifdef MEM_OPT_VERSION
  //the way to get the cluster size if the atom ids of the cluster are
  //contiguous. The input parameter is the atom id that leads the cluster.
  int get_cluster_size_con(int aid) const { return clusterSigs[aid]; }  
  //the way to get the cluster size if the atoms ids of the cluster are
  //not contiguous. The input parameter is the cluster index.
  int get_cluster_size_uncon(int cIdx) const { return clusterSize[cIdx]; }
  int get_cluster_idx(int aid) const { return clusterSigs[aid]; }
  int get_num_clusters() const { return numClusters; }
  #else
  int get_cluster(int anum) const { return cluster[anum]; }
  int get_clusterSize(int anum) const { return clusterSize[anum]; }
  #endif

#ifndef MEM_OPT_VERSION
  const float *getOccupancyData() { return (const float *)occupancy; }
  void setOccupancyData(molfile_atom_t *atomarray);
  void freeOccupancyData() { delete [] occupancy; occupancy=NULL; }

  const float *getBFactorData() { return (const float *)bfactor; }
  void setBFactorData(molfile_atom_t *atomarray);
  void freeBFactorData() { delete [] bfactor; bfactor=NULL; }
#endif

  //  Get the mass of an atom
  Real atommass(int anum) const
  {
    #ifdef MEM_OPT_VERSION
    return atomMassPool[eachAtomMass[anum]];
    #else
    return(atoms[anum].mass);
    #endif
  }

  //  Get the charge of an atom
  Real atomcharge(int anum) const
  {
    #ifdef MEM_OPT_VERSION
    return atomChargePool[eachAtomCharge[anum]];
    #else
    return(atoms[anum].charge);
    #endif
  }
  
  //  Get the vdw type of an atom
  Index atomvdwtype(int anum) const
  {      
      return(atoms[anum].vdw_type);
  }

  #ifndef MEM_OPT_VERSION
  //  Retrieve a bond structure
  Bond *get_bond(int bnum) const {return (&(bonds[bnum]));}

  //  Retrieve an angle structure
  Angle *get_angle(int anum) const {return (&(angles[anum]));}

  //  Retrieve an improper strutcure
  Improper *get_improper(int inum) const {return (&(impropers[inum]));}

  //  Retrieve a dihedral structure
  Dihedral *get_dihedral(int dnum) const {return (&(dihedrals[dnum]));}

  //  Retrieve a cross-term strutcure
  Crossterm *get_crossterm(int inum) const {return (&(crossterms[inum]));}
  #endif

  // DRUDE: retrieve lphost structure
  Lphost *get_lphost(int atomid) const {
    // don't call unless simParams->drudeOn == TRUE
    // otherwise lphostIndexes array doesn't exist!
    int index = lphostIndexes[atomid];
    return (index != -1 ? &(lphosts[index]) : NULL);
  }
  // DRUDE

  #ifndef MEM_OPT_VERSION
  Bond *getAllBonds() const {return bonds;}
  Angle *getAllAngles() const {return angles;}
  Improper *getAllImpropers() const {return impropers;}
  Dihedral *getAllDihedrals() const {return dihedrals;}
  Crossterm *getAllCrossterms() const {return crossterms;}
  #endif

  // DRUDE: retrieve entire lphosts array
  Lphost *getAllLphosts() const { return lphosts; }
  // DRUDE

  //  Retrieve a hydrogen bond donor structure
  Bond *get_donor(int dnum) const {return (&(donors[dnum]));}  

  //  Retrieve a hydrogen bond acceptor structure
  Bond *get_acceptor(int dnum) const {return (&(acceptors[dnum]));} 

  Bond *getAllDonors() const {return donors;}
  Bond *getAllAcceptors() const {return acceptors;}

  //  Retrieve an exclusion structure
  #ifndef MEM_OPT_VERSION
  Exclusion *get_exclusion(int ex) const {return (&(exclusions[ex]));}
  #endif

  //  Retrieve an atom type
  const char *get_atomtype(int anum) const
  {
    if (atomNames == NULL)
    {
      NAMD_die("Tried to find atom type on node other than node 0");
    }

    #ifdef MEM_OPT_VERSION    
    return atomTypePool[atomNames[anum].atomtypeIdx];
    #else
    return(atomNames[anum].atomtype);
    #endif
  }

  //  Lookup atom id from segment, residue, and name
  int get_atom_from_name(const char *segid, int resid, const char *aname) const;

  //  Lookup number of atoms in residue from segment and residue
  int get_residue_size(const char *segid, int resid) const;

  //  Lookup atom id from segment, residue, and index in residue
  int get_atom_from_index_in_residue(const char *segid, int resid, int index) const;

  
  //  The following routines are used to get the list of bonds
  //  for a given atom.  This is used when creating the bond lists
  //  for the force objects  

  #ifndef MEM_OPT_VERSION
  int32 *get_bonds_for_atom(int anum)
      { return bondsByAtom[anum]; } 
  int32 *get_angles_for_atom(int anum) 
      { return anglesByAtom[anum]; }
  int32 *get_dihedrals_for_atom(int anum) 
      { return dihedralsByAtom[anum]; }
  int32 *get_impropers_for_atom(int anum) 
      { return impropersByAtom[anum]; }  
  int32 *get_crossterms_for_atom(int anum) 
      { return crosstermsByAtom[anum]; }  
  int32 *get_exclusions_for_atom(int anum)
      { return exclusionsByAtom[anum]; }
  const int32 *get_full_exclusions_for_atom(int anum) const
      { return fullExclusionsByAtom[anum]; }
  const int32 *get_mod_exclusions_for_atom(int anum) const
      { return modExclusionsByAtom[anum]; }
  #endif
  
  //  Check for exclusions, either explicit or bonded.
        //  Returns 1 for full, 2 for 1-4 exclusions.
  #ifdef MEM_OPT_VERSION
  int checkExclByIdx(int idx1, int atom1, int atom2) const;
  const ExclusionCheck *get_excl_check_for_idx(int idx) const{      
      return &exclChkSigPool[idx];
  }
  #else
  int checkexcl(int atom1, int atom2) const;

  const ExclusionCheck *get_excl_check_for_atom(int anum) const{      
      return &all_exclusions[anum];             
  }
  #endif

/* BEGIN gf */
  // Return true or false based on whether or not the atom
  // is subject to grid force
  Bool is_atom_gridforced(int atomnum, int gridnum) const
  {
      if (numGridforceGrids)
      {
	  return(gridfrcIndexes[gridnum][atomnum] != -1);
      }
      else
      {
	  return(FALSE);
      }
  }
/* END gf */

#ifndef MEM_OPT_VERSION
  //  Return true or false based on whether the specified atom
  //  is constrained or not.
  Bool is_atom_constrained(int atomnum) const
  {
    if (numConstraints)
    {
      //  Check the index to see if it is constrained
      return(consIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }
#endif

  //  Return true or false based on whether the specified atom
  //  is moving-dragged or not.
  Bool is_atom_movdragged(int atomnum) const
  {
    if (numMovDrag)
    {
      //  Check the index to see if it is constrained
      return(movDragIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }

  //  Return true or false based on whether the specified atom
  //  is rotating-dragged or not.
  Bool is_atom_rotdragged(int atomnum) const
  {
    if (numRotDrag)
    {
      //  Check the index to see if it is constrained
      return(rotDragIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }

  //  Return true or false based on whether the specified atom
  //  is "constant"-torqued or not.
  Bool is_atom_constorqued(int atomnum) const
  {
    if (numConsTorque)
    {
      //  Check the index to see if it is constrained
      return(consTorqueIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }

#ifndef MEM_OPT_VERSION
  //  Get the harmonic constraints for a specific atom
  void get_cons_params(Real &k, Vector &refPos, int atomnum) const
  {
    k = consParams[consIndexes[atomnum]].k;
    refPos = consParams[consIndexes[atomnum]].refPos;
  }
#endif

/* BEGIN gf */
  void get_gridfrc_params(Real &k, Charge &q, int atomnum, int gridnum) const
  {
      k = gridfrcParams[gridnum][gridfrcIndexes[gridnum][atomnum]].k;
      q = gridfrcParams[gridnum][gridfrcIndexes[gridnum][atomnum]].q;
  }
  
  GridforceGrid* get_gridfrc_grid(int gridnum) const
  {
      GridforceGrid *result = NULL;
      if (gridnum >= 0 && gridnum < numGridforceGrids) {
	  result = gridfrcGrid[gridnum];
      }
      return result;
  }
  
  int set_gridfrc_grid(int gridnum, GridforceGrid *grid)
  {
      if (grid && gridnum >= 0 && gridnum < numGridforceGrids) {
	  gridfrcGrid[gridnum] = grid;
	  return 0;
      } else {
	  return -1;
      }
  }
/* END gf */

  Real langevin_param(int atomnum) const
  {
    return(langevinParams ? langevinParams[atomnum] : 0.);
  }

  //  Get the stirring constraints for a specific atom
  void get_stir_refPos(Vector &refPos, int atomnum) const
  {
    refPos = stirParams[stirIndexes[atomnum]].refPos;
  }


  void put_stir_startTheta(Real theta, int atomnum) const
  {
    stirParams[stirIndexes[atomnum]].startTheta = theta;
  }


  Real get_stir_startTheta(int atomnum) const
  {
    return stirParams[stirIndexes[atomnum]].startTheta;
  }
 

  //  Get the moving drag factor for a specific atom
  void get_movdrag_params(Vector &v, int atomnum) const
  {
    v = movDragParams[movDragIndexes[atomnum]].v;
  }

  //  Get the rotating drag pars for a specific atom
  void get_rotdrag_params(BigReal &v, Vector &a, Vector &p, 
        int atomnum) const
  {
    v = rotDragParams[rotDragIndexes[atomnum]].v;
    a = rotDragParams[rotDragIndexes[atomnum]].a;
    p = rotDragParams[rotDragIndexes[atomnum]].p;
  }

  //  Get the "constant" torque pars for a specific atom
  void get_constorque_params(BigReal &v, Vector &a, Vector &p, 
        int atomnum) const
  {
    v = consTorqueParams[consTorqueIndexes[atomnum]].v;
    a = consTorqueParams[consTorqueIndexes[atomnum]].a;
    p = consTorqueParams[consTorqueIndexes[atomnum]].p;
  }

//fepb
  unsigned char get_fep_type(int anum) const
  {
    return(fepAtomFlags[anum]);
  }

  /* BKR - Get the FEP type (i.e. 0, 1, or 2) of a bonded _interaction_ based 
     on the atom indices of the atoms involved (internally converted to FEP 
     types) and the order of the bonded interaction (i.e. 2, 3, or 4). This 
     automatically accounts for whether or not purely alchemical interactions
     are being scaled (e.g. bonds between two atoms of fep_type 1).

     The logic here is admittedly a bit opaque. When the fep_type's are back
     mapped to -1,0,1, we can use the sum of all of the types to determine the 
     type of interaction for all bonded term types:

     0  - only non-alchemical atoms involved
     >0 - _atleast_ one appearing atom
     <0 - _atleast_ one disappearing atom

     If the magnitude of the sum is equal to the interaction order (i.e. 2, 3, 
     or 4), then it is a _purely_ alchemical interaction and it may be 
     desireable to retain it for sampling purposes (note that this adds a 
     constant to the free energy that will almost always differ at the 
     endpoints).  In order to avoid unexpected instability these interactions 
     are retained by default, but can be scaled along with everything else by 
     setting "alchBondDecouple on".

     NB: All pure alchemical interactions beyond order 2 are ALWAYS discarded
     by alchemify. This is good, because higher order interactions would break
     the above logic. For example, if we had an angle term between atoms of 
     types (1,1,-1) the sum would be 1, but this term should receive no scaling
     because it involves groups -1 and 1 but not 0.
  */
  int get_fep_bonded_type(const int *atomID, unsigned int order) const
  {
    int typeSum = 0;
    for ( int i=0; i < order; ++i ) {
      typeSum += (fepAtomFlags[atomID[i]] == 2 ? -1 : fepAtomFlags[atomID[i]]);
    }
    // Increase the cutoff if scaling purely alchemical bonds. 
    // This really only applies when order = 2.
    if ( simParams->alchBondDecouple ) order++;

    if ( typeSum == 0 ) return 0; // Most interactions get caught here.
    else if ( 0 < typeSum && typeSum < order ) return 1;
    else if ( 0 > typeSum && typeSum > -order ) return 2;

    if ( simParams->alchBondDecouple ) {
      // Alchemify should always keep this from bombing, but just in case...
      NAMD_die("Unexpected alchemical bonded interaction!");
    }
  }
//fepe

#ifndef MEM_OPT_VERSION
  Bool is_atom_fixed(int atomnum) const
  {
    return (numFixedAtoms && fixedAtomFlags[atomnum]);
  }
#else
  //Since binary search is more expensive than direct array access,
  //and this function is usually called for consecutive atoms in this context,
  //the *listIdx returns the index to the range of atoms [aid1, aid2]
  //that are fixed. If the atom aid is fixed, then aid1=<aid<=aid2;
  //If the atom aid is not fixed, then aid1 indicates the smallest fixed atom
  //id that is larger than aid; so the listIdx could be equal the size of
  //fixedAtomsSet. --Chao Mei
  Bool is_atom_in_set(AtomSetList *localAtomsSet, int aid, int *listIdx) const;
  inline Bool is_atom_fixed(int aid, int *listIdx=NULL) const {
    return is_atom_in_set(fixedAtomsSet,aid,listIdx);
  }
  inline Bool is_atom_constrained(int aid, int *listIdx=NULL) const {
    return is_atom_in_set(constrainedAtomsSet,aid,listIdx);
  }
#endif
        
  Bool is_atom_stirred(int atomnum) const
  {
    if (numStirredAtoms)
    {
      //  Check the index to see if it is constrained
      return(stirIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }
  

  Bool is_group_fixed(int atomnum) const
  {
    return (numFixedAtoms && (fixedAtomFlags[atomnum] == -1));
  }
  Bool is_atom_exPressure(int atomnum) const
  {
    return (numExPressureAtoms && exPressureAtomFlags[atomnum]);
  }
  // 0 if not rigid or length to parent, for parent refers to H-H length
  // < 0 implies MOLLY but not SHAKE, > 1 implies both if MOLLY is on
  Real rigid_bond_length(int atomnum) const
  {
    return(rigidBondLengths[atomnum]);
  }

  void print_atoms(Parameters *); 
        //  Print out list of atoms
  void print_bonds(Parameters *); 
        //  Print out list of bonds
  void print_exclusions();//  Print out list of exclusions

public:  
  int isOccupancyValid, isBFactorValid;

#ifdef MEM_OPT_VERSION
  //read the per-atom file for the memory optimized version where the file 
  //name already exists the range of atoms to be read are 
  //[fromAtomID, toAtomID].
  void read_binary_atom_info(int fromAtomID, int toAtomID, InputAtomList &inAtoms);

  int getNumCalcExclusions(){return numCalcExclusions;}
  void setNumCalcExclusions(int x){numCalcExclusions= x;}

  Index getEachAtomMass(int i){return eachAtomMass[i];}
  Index getEachAtomCharge(int i){return eachAtomCharge[i];}

  ExclSigID getAtomExclSigId(int aid) const {
      return eachAtomExclSig[aid];
  }

  Real *getAtomMassPool(){return atomMassPool;}
  Real *getAtomChargePool(){return atomChargePool;}
  AtomCstInfo *getAtoms(){return atoms;}  

  int atomSigPoolSize;
  AtomSignature *atomSigPool;

  /* All the following are temporary variables for reading the compressed psf file */
  //declarations for atoms' constant information  
  int segNamePoolSize; //Its value is usually less than 5
  char **segNamePool; //This seems not to be important, but it only occupied very little space.

  int resNamePoolSize;
  char **resNamePool;

  int atomNamePoolSize;
  char **atomNamePool;

  int atomTypePoolSize;
  char **atomTypePool;

  int chargePoolSize;
  Real *atomChargePool;

  int massPoolSize;
  Real *atomMassPool;

  AtomSigID getAtomSigId(int aid) {
      return eachAtomSig[aid]; 
  }

  //Indicates the size of both exclSigPool and exclChkSigPool
  int exclSigPoolSize;
  //this will be deleted after build_lists_by_atom
  ExclusionSignature *exclSigPool;
  //This is the final data structure we want to store  
  ExclusionCheck *exclChkSigPool;

  void addNewExclSigPool(const std::vector<ExclusionSignature>&);  

  void delEachAtomSigs(){    
      //for NAMD-smp version, only one Molecule object is held
      //on each node, therefore, only one deletion operation should
      //be taken on a node, otherwise, there possibly would be some
      //wierd memory problems. The same reason applies to other deletion
      //operations inside the Molecule object.   
      if(CmiMyRank()) return;

      delete [] eachAtomSig;
      delete [] eachAtomExclSig;
      eachAtomSig = NULL;
      eachAtomExclSig = NULL;
  }

  void delChargeSpace(){
      if(CmiMyRank()) return;

      delete [] atomChargePool;
      delete [] eachAtomCharge;
      atomChargePool = NULL;
      eachAtomCharge = NULL;
  }
  
  void delMassSpace(){
      if(CmiMyRank()) return;

      delete [] atomMassPool;
      delete [] eachAtomMass;
      atomMassPool = NULL;
      eachAtomMass = NULL;
  }
  
  void delClusterSigs() {
      if(CmiMyRank()) return;      

      delete [] clusterSigs;
      clusterSigs = NULL;
  }

  void delAtomNames(){
      if(CmiMyRank()) return;
      delete [] atomNamePool;
      delete [] atomNames;
      atomNamePool = NULL;
      atomNames = NULL;
  }

  void delFixedAtoms(){
      if(CmiMyRank()) return;
      delete fixedAtomsSet;
      fixedAtomsSet = NULL;
  }

private:
  Index insert_new_mass(Real newMass);

#endif

// Go stuff
public:

GoValue go_array[MAX_GO_CHAINS*MAX_GO_CHAINS];    //  Array of Go params -- JLai
int go_indices[MAX_GO_CHAINS+1];        //  Indices from chainIDS to go array -- JLai
int NumGoChains;                        //  Number of Go chain types -- JLai

// Declares and initializes Go variables
void goInit();

// Builds explicit Gromacs pairs
void build_gro_pair();

// Builds the initial Go parameters 
void build_go_params(StringList *);

//  Read Go parameter file
void read_go_file(char *);

//  Get Go cutoff for a given chain type pair
Real get_go_cutoff(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].cutoff; };

//  Get Go epsilonRep for a given chain type pair
Real get_go_epsilonRep(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].epsilonRep; };

//  Get Go sigmaRep for a given chain type pair
Real get_go_sigmaRep(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].sigmaRep; };

//  Get Go epsilon for a given chain type pair
Real get_go_epsilon(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].epsilon; };

//  Get Go exp_a for a given chain type pair
int get_go_exp_a(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].exp_a; };

//  Get Go exp_b for a given chain type pair
int get_go_exp_b(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].exp_b; };

//  Get Go exp_rep for a given chain type pair
int get_go_exp_rep(int chain1, int chain2) { return go_array[MAX_GO_CHAINS*chain1 + chain2].exp_rep; };

//  Whether residue IDs with this difference are restricted
Bool go_restricted(int, int, int);

// Prints Go Params
void print_go_params();

void initialize();

void send_GoMolecule(MOStream *);
//  send the molecular structure 
//  from the master to the clients

void receive_GoMolecule(MIStream *);
//  receive the molecular structure
//  from the master on a client
};

#endif

