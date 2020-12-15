#ifndef GROMACSTOPFILE_H
#define GROMACSTOPFILE_H

typedef float Real;

#define NAMESIZE 5
#define LONGNAMESIZE 20

/* A GenericBond represents a bond between two atoms in a template
   that can be used one or more times within a simulation.  This is an
   immutable type.
   Here, atomi and atomj start at 0. */
class GenericBond {
 private:
  int atomi; /* the first atom of the bond */
  int atomj; /* the second atom of the bond */
  int type; /* an index into the bond's parameters, stored elsewhere */

 public:
  /* initializes this to be a bond between atoms <i> and <j> of type
     <type> */
  GenericBond(int i, int j, int theType);

  int getAtomi() const {return atomi;}
  int getAtomj() const {return atomj;}
  int getType() const {return type;}

};

/* A GenericAngle represents a angle between two atoms in a template
   that can be used one or more times within a simulation.  This is an
   immutable type.
   Here, atomi and atomj start at 0. */
class GenericAngle {
 private:
  int atomi; /* the first atom of the angle */
  int atomj; /* the second atom of the angle */
  int atomk; /* the third atom of the angle */
  int type; /* an index into the angle's parameters, stored elsewhere */

 public:
  /* initializes this to be a angle between atoms <i>, <j>, and <k> of
     type <type> */
  GenericAngle(int i, int j, int k, int theType);

  int getAtomi() const {return atomi;}
  int getAtomj() const {return atomj;}
  int getAtomk() const {return atomk;}
  int getType() const {return type;}

};

/* A GenericDihedral represents a dihedral angle between four atoms in
   a template that can be used one or more times within a simulation.
   This is an immutable type.
   Here, atomi and atomj start at 0. */
class GenericDihedral {
 private:
  int atomi; /* the first atom of the angle */
  int atomj; /* the second atom of the angle */
  int atomk; /* the third atom of the angle */
  int atoml; /* the fourth atom of the angle */
  int type; /* an index into the parameters, stored elsewhere */

 public:
  /* initializes this to be a angle between atoms <i>, <j>, <k>, and
     <l> of type <type> */
  GenericDihedral(int i, int j, int k, int l, int theType);

  int getAtomi() const {return atomi;}
  int getAtomj() const {return atomj;}
  int getAtomk() const {return atomk;}
  int getAtoml() const {return atoml;}
  int getType() const {return type;}

};

/* A GenericAtom represents an atom in a template that can be used one
   or more times within a simulation.  This is an immutable type. */
class GenericAtom {
 private:
  char type[NAMESIZE+1];
  int typeNum;
  int resNum;
  char resType[NAMESIZE+1];
  char atomName[NAMESIZE+1];
  Real charge;
  Real mass;

 public: 

  /* Initializes this to be the atom specified by the given
     parameters.  All the strings will be copied, and must be at most
     four characters long. */
  GenericAtom(const char *theType, int theTypeNum, int theResNum,
	      const char *theResType,
	      const char *theAtomName, Real theCharge, Real theMass);

  /* The following just return the parameters of the atom. */
  const char *getType() const { return type; }
  int getTypeNum() const { return typeNum; }
  int getResNum() const {return resNum;}
  const char *getResType() const {return resType;}
  const char *getAtomName() const {return atomName;}
  Real getCharge() const {return charge;}
  Real getMass() const {return mass;}

};

/* A GenericMol represents a complete molecule template.  This will be
   a mutable type, for convenience's sake.  The atoms here are
   numbered starting at ZERO */
class GenericMol {
 private:
  /* the name of this Molecule, so it can be identified later */
  const char *name;

  /* the rest of the list */
  const GenericMol *next;

  /* the list of GenericAtoms */
  ResizeArray <const GenericAtom *> atomList;

  /* the list of GenericBonds */
  ResizeArray <const GenericBond *> bondList;
  
  /* the list of GenericAngles */
  ResizeArray <const GenericAngle *> angleList;
  
  /* the list of GenericDihedrals */
  ResizeArray <const GenericDihedral *> dihedralList;
  
 public:
  /* initializes this to be the molecule with name <name> */
  GenericMol(const char *theName);

  /* getNumAtoms returns the number of atoms stored here. */
  int getNumAtoms() const { return atomList.size(); }
  int getNumBonds() const { return bondList.size(); }
  int getNumAngles() const { return angleList.size(); }
  int getNumDihedrals() const { return dihedralList.size(); }
  int getNumRes() const { return atomList[atomList.size()-1]->getResNum(); }

  /* adds an atom to the end of list.  The residue numbers must start
     at zero and count up in unit steps */
  void addAtom(const char *theType, int theTypeNum, int theResNum,
	       const char *theResType,
	       const char *theAtomName, Real theCharge, Real theMass);

  /* adds a bond/angle/dihedral to the molecule */
  void addBond(int atomi, int atomj, int type);
  void addAngle(int atomi, int atomj, int atomk, int type);
  void addDihedral(int atomi, int atomj, int atomk, int atoml, int type);

  /* gets the name */
  const char *getName() const { return name; }

  /* gets atom <n>, where atom 0 is the first atom in the list. */
  const GenericAtom *getAtom(int n) const;

  /* gets bond <n>, where bond 0 is the first bond in the list. */
  const GenericBond *getBond(int n) const;

  /* gets angle <n>, where angle 0 is the first angle in the list. */
  const GenericAngle *getAngle(int n) const;

  /* gets dihedral <n>, where angle 0 is the first one in the list. */
  const GenericDihedral *getDihedral(int n) const;

};

/* A MolInst represents a number of instances of a molecule.  These
   objects are meant to be assembled (in order) into a list, to
   represent all the instances of molecules in a system. */
class MolInst {
  int num;
  const GenericMol *mol;
 public:
  /* initializes this to represent <theNum> copies of <theMol> */
  MolInst(const GenericMol *theMol, int theNum);

  /* get the number */
  int getNum() const { return num; }

  /* get the total number of X here */
  int getNumAtoms() const;
  int getNumRes() const;
  int getNumBonds() const;
  int getNumAngles() const;
  int getNumDihedrals() const;
  
  /* get the molecule */
  const GenericMol *getMol() const { return mol; }
};

/* An AtomTable stores all the known types of atoms.
   mass: the mass in amu
   charge: the charge in e
   c6: the LJ parameter c6 in kcal/mol A^6
   c12: the LJ parameter c12 in kcal/mol A^12 */
class AtomTable {
 private:
  /* implementation of the extensible array */
  ResizeArray<Real> mArray;
  ResizeArray<Real> qArray;
  ResizeArray<Real> c6Array;
  ResizeArray<Real> c12Array;
  ResizeArray<const char *> typeArray; /* should never be NULL */

 public:
  /* this adds a entry for an atom type to the table */
  void addType(const char *type, Real m, Real q, Real c6, Real c12);

  /* This looks up the atom by type - if it cannot be found, this
     returns -1, otherwise this returns the index of the atomtype */
  int getParams(const char *type, Real *m, Real *q,
		     Real *c6, Real *c12) const;

  /* This looks up the atom type by number, returning it in the string
     <type>, which must be at least 11 characters long. */
  void getType(int num, char *type) const;


  /* returns the number of types in the table */
  int size() const { return mArray.size(); }
};

/* A BondTable stores all the different types of bonds, so that we
   only have to have one copy of each.
   b0: the natural bond length in A.
   kB: the spring constant in kcal/mol/A^2, where E=kx^2 to 1st order.
   funct: 1 for harmonic potential, 2 for fourth-order GROMACS96
   potential. */
class BondTable {
 private:
  /* implementation of the extensible array */
  ResizeArray<Real> b0Array;
  ResizeArray<Real> kBArray;
  ResizeArray<int> functArray;
  ResizeArray<const char *> typeaArray; /* NULL if no type is set */
  ResizeArray<const char *> typebArray; /* NULL if no type is set */

 public:
  /* returns the number of types in the table */
  int size() const { return b0Array.size(); }

  /* this gets the index of a bond in the table (adding an entry if
     none exists). */
  int getIndex(Real b0, Real kB, int funct);

  /* this adds a entry for a bond type to the table, including the two
     atoms involved in the bond */
  void addType(const char *typea, const char *typeb, Real b0,
		   Real kB, int funct);

  /* getParams puts the parameters for bond-type <num> into the
     spaces pointed to by the other arguments. */
  void getParams(int num, Real *b0, Real *kB, int *funct) const;

  /* This version looks up the bond by atom type - the order of the
     types doesn't matter!  If the specified bond/function combination
     cannot be found, it returns -1, otherwise it returns the index of
     the bondtype. */
  int getParams(const char *typea, const char *typeb, int funct,
		     Real *b0, Real *kB) const;
};

/* A AngleTable stores all the different types of angles, so that we
   only have to have one copy of each.  Units:
   th0 - degrees
   kth - kcal/mol/rad^2 (funct 1)  OR  kcal/mol (funct 2)
   funct 1: harmonic angle potential
   funct 2: harmonic G96 cosine potential
   The atoms are in standard geometry order - A--B--C for angle /ABC
*/
class AngleTable {
 private:
  /* implementation of the extensible array */
  ResizeArray<Real> th0Array;
  ResizeArray<Real> kthArray;
  ResizeArray<int> functArray;
  ResizeArray<const char *> typeaArray; /* NULL if no type is set */
  ResizeArray<const char *> typebArray; /* NULL if no type is set */
  ResizeArray<const char *> typecArray; /* NULL if no type is set */

 public:
  /* returns the number of types in the table */
  int size() const { return th0Array.size(); }

  /* this gets the index of a angle in the table (adding an entry if
     none exists).
  */
  int getIndex(Real th0, Real kth, int funct);

  /* this adds a entry for a angle type to the table, including the
     three atoms involved in the angle */
  void addType(const char *typea, const char *typeb,
		    const char *typec, Real th0,
		    Real kth, int funct);

  /* getAngleParams puts the parameters for angle-type <num> into the
     spaces pointed to by the other arguments. 
  */
  void getParams(int num, Real *th0, Real *kth, int *funct) const;

  /* This version looks up the angle by atom type - the direction of
     the types doesn't matter (it can be A--B--C or C--B--A)!  If the
     specified angle/function combination cannot be found, it returns
     -1, otherwise it returns the index of the angletype. */
  int getParams(const char *typea, const char *typeb,
		     const char *typec, int funct,
		     Real *th0, Real *kth) const;
};

/* A DihedralTable stores all the different types of angles, so that we
   only have to have one copy of each.  Units:
   funct 1: proper dihedral
     c0 - thmax (deg)
     c1 - fc (kcal/mol)
     a,b are the two inner atoms
   funct 2: improper dihedral
     c0 - th0 (deg)
     c1 - fc (kcal/mol)
     a,b are the two outer atoms
   funct 3: RB dihedral
     c0-c5 - C0-C5 (kcal/mol)
     a,b are the two inner atoms

   The atoms are in standard geometry order - A--B--C--D for the
   dihedral between the planes containing ABC and BCD.
*/
class DihedralTable {
 private:
  /* implementation of the extensible array */
  ResizeArray<int> functArray;
  ResizeArray<const char *> typeaArray; /* NULL if no type is set */
  ResizeArray<const char *> typebArray; /* NULL if no type is set */
  ResizeArray<int> multArray; /* 0 if funct is not 1 */

  /* cArray is where the actual paramters live - each Real * is a
     pointer to a list of Reals that is at least as long as the
     number of parameters required for the given type. */
  ResizeArray<const Real *> cArray;

 public:
  /* returns the number of types in the table */
  int size() const { return functArray.size(); }

  /* this function gets the index of a dihedral angle in the table
     (adding an entry if none exists).  The required number of
     parameters (see notes on class DihedralTable) must be stored in
     the array <c> */
  int getIndex(const Real *c, int mult, int funct);

  /* this adds a entry for a angle type to the table, including the
     two atoms involved in the angle.  The required number of
     parameters (see notes on class DihedralTable) must be stored in
     the array <c>.  Note that these two angles really refer to either
     atoms A and D or B and C depending on the dihedral type.  */
  void addType(const char *typea, const char *typeb,
	       const Real *c, int mult, int funct);

  /* getParams puts the parameters for angle-type <num> into the
     spaces pointed to by the other arguments.  The required number of
     parameters (see notes on class DihedralTable) will be stored in
     the array <c>, so make sure that c has size >= 6 */
  void getParams(int num, Real *c, int *mult, int *funct) const;

  /* This version looks up the angle by atom type - the direction of
     the types doesn't matter (it can be ABCD or DCBA)!  If the
     specified angle/function combination cannot be found, it returns
     -1, otherwise it returns the index of the dihedraltye.  This is a
     little strange because the different dihedral types use different
     atoms - but since we always will know all four AND the funct when
     we need to look up a dihedral, this will work.  The required
     number of parameters (see notes on class DihedralTable) will be
     stored in the array <c>, so make sure that c is big enough for
     the specified function. */
  int getParams(const char *typea, const char *typeb,
		const char *typec, const char *typed, int funct,
		Real *c, int *mult) const;
};

/* this stores the VDW interaction parameters as a function of the
   atom types. */
class VDWTable {
 private:
  /* atom types a and b */
  ResizeArray<const char *> typeAArray; /* NULL if no type is set */
  ResizeArray<const char *> typeBArray; /* NULL if no type is set */

  /* the VDW parameters */
  ResizeArray<Real> c6Array;
  ResizeArray<Real> c12Array;
  ResizeArray<Real> c6PairArray;  // these two are
  ResizeArray<Real> c12PairArray; // for 1-4 pairs
  
  /* finds the index from the two interacting types, or returns -1 */
  int getIndex(const char *typea, const char *typeb) const;
  
 public:
  /* these add VDW parameters to the list */
  void addType(const char *typea, const char *typeb, Real c6, Real c12);
  void add14Type(const char *typea, const char *typeb,
		 Real c6pair, Real c12pair);

  /* this returns the VDW interaction parameters that were added to
     the list earlier, and returns the index or the parameters (a
     useless number) or -1 if it can't find the specified types. */
  int getParams(const char *typea, const char *typeb,
		 Real *c6, Real *c12, Real *c6pair, Real *c12pair) const;
};

/* this stores the pair interaction parameters as a function of the
   atom types. 
   JLai August 16th, 2012
*/
struct GroLJPair {
  int indxLJA;
  int indxLJB;
  Real c6pair;
  Real c12pair;
};
  
struct GroGaussPair {
  int indxGaussA;
  int indxGaussB;
  Real gA;
  Real gMu1;
  Real giSigma1;
  Real gMu2;
  Real giSigma2;
  Real gRepulsive;
};

class PairTable {
  private:
  // Data structure for LJ Structure-based potential
  ResizeArray<int> indxLJA;
  ResizeArray<int> indxLJB;
  ResizeArray<Real> c6pair;
  ResizeArray<Real> c12pair;
  ResizeArray<GroLJPair> pairlistLJ;

  // Data structure for Gaussian Structure-based potential
  ResizeArray<int> indxGaussA;
  ResizeArray<int> indxGaussB;
  ResizeArray<Real> gA;
  ResizeArray<Real> gMu1;
  ResizeArray<Real> giSigma1;
  ResizeArray<Real> gMu2;
  ResizeArray<Real> giSigma2;
  ResizeArray<Real> gRepulsive;
  ResizeArray<GroGaussPair> pairlistGauss;

 public:
  /* these add pair parameters to the list */
  int addPairLJType2(int typea, int typeb, Real c6, Real c12);
  int addPairGaussType2(int typea, int typeb, Real gA,
			Real gMu1, Real gSigma1);
  int addPairGaussType2(int typea, int typeb, Real gA,
			Real gMu1, Real gSigma1, Real gRepulsive);
  int addPairGaussType2(int typea, int typeb, Real gA, 
			Real gMu1, Real gSigma1, Real gMu2, Real gSigma2, Real gRepulsive);
  void getPairLJArrays2(int *indexA, int *indexB, Real *pairC6, Real *pairC12);
  void getPairGaussArrays2(int *indexA, int *indexB, Real *gaussA, Real *gaussMu1,
			  Real *gaussSigma1, Real *gaussMu2, Real *gaussSigma2,
			   Real *gaussRepulsive);

  static bool GroLJCompare (GroLJPair, GroLJPair);
  static bool GroGaussCompare (GroGaussPair, GroGaussPair);
  
};

/* A GromacsTopFile represents the information stored in a GROMACS
   topology file.  This is an immutable type. */
class GromacsTopFile {
 private:
  /* The system name */
  char *systemName;

  /* 1-4 parameter scaling values */
  Real fudgeLJ, fudgeQQ;

  /* whether or not to generate the 1-4 interactions automatically */
  int genPairs;

  /* a list of molecule templates */
  ResizeArray <GenericMol *> genericMols;

  /* the list of molecule instances */
  ResizeArray <const MolInst *> molInsts;
 
  /* the table of bonds/angles/atoms */
  AtomTable atomTable;
  BondTable bondTable;
  AngleTable angleTable;
  DihedralTable dihedralTable;
  VDWTable vdwTable;
  PairTable pairTable;

 public:

  /* initializes this to represent the data stored in the topology
     file <filename>, or exits on error. */
  GromacsTopFile(char *filename);

  /* getSystemName returns the name of the system */
  char *getSystemName() const { return systemName; }

  /* returns the number of atoms/bonds/angles/dihedrals in the file */
  int getNumAtoms() const;
  int getNumBonds() const;
  int getNumAngles() const;
  int getNumDihedrals() const;

  /* returns the number of atoms/bonds/angles/dihedrals params in the file */
  int getNumAtomParams() const { return atomTable.size(); }
  int getNumBondParams() const { return bondTable.size(); }
  int getNumAngleParams() const { return angleTable.size(); }
  int getNumDihedralParams() const { return dihedralTable.size(); }

  /* getAtom puts the information about atom number <num> into the
     spaces pointed to by the other arguments.  The string buffers
     must be at least 11 characters long.  Atom number 0 is the first
     atom in the list, so that it corresponds with our normal idea of
     atom numbers.
     charge - e
     mass   - amu
  */
  void getAtom(int num, int *residue_number, char *residue_name,
	       char *atom_name, char *atom_type, int *atom_typenum,
	       Real *charge, Real *mass)
    const;

  /* getAtomParams puts the parameters for atom-type <num> into the
     spaces pointed to by the other arguements.  Currently, all anyone
     cares about is the atom type string, so that's all we're
     returning!  This may change later, but for now we get everything
     else with the getAtom() method.  The string buffer must be at
     least 11 characters long. */
  void getAtomParams(int num, char *type) const {
    atomTable.getType(num,type);
  }

  // JLai
  int getNumPair() const;
  int getNumLJPair() const;
  int getNumGaussPair() const;
  int getNumExclusions() const;
  void getExclusions(int *, int *) const;
  // End of JLai

  /* getBond puts the information about bond number <num> into the
     spaces pointed to by the other arguments.  Bond number 0 is the
     first bond in the list.
  */
  void getBond(int num, int *atomi, int *atomj, int *bondtype) const;

  /* getBondParams puts the parameters for bond-type <num> into the
     spaces pointed to by the other arguments. 
     b0 - natural length in A
     kB - spring constant in kcal/mol/A^2 - E=kx^2
     funct - 1 for harmonic, 2 for special fourth-order GROMOS96 thing */
  void getBondParams(int num, Real *b0, Real *kB, int *funct) const;

  /* getAngle puts the information about angle number <num> into the
     spaces pointed to by the other arguments.  Angle number 0 is the
     first angle in the list.
  */
  void getAngle(int num, int *atomi, int *atomj, int *atomk,
		int *angletype) const;

  /* getAngleParams puts the parameters for angle-type <num> into the
     spaces pointed to by the other arguments.
     th0 - natural angle in degrees
     kth - spring constant in kcal/rad - E=kx^2
     funct - 1 for harmonic, 2 for special GROMOS96 thing */
  void getAngleParams(int num, Real *th0, Real *kth, int *funct) const;

  /* getDihedral puts the information about dihedral number <num> into
     the spaces pointed to by the other arguments.  Dihedral number 0
     is the first angle in the list. */
  void getDihedral(int num, int *atomi, int *atomj, int *atomk,
		   int *atoml, int *type) const;

  /* getDihedralParams puts the parameters for dihedral-type <num>
     into the spaces pointed to by the other arguments.  The array of
     Reals <c> must be of size >= 6, since up to 6 parameters will be
     set there.  These parameters are as follows:
   funct 1: proper dihedral
     c0 - thmax (deg)
     c1 - fc (kcal/mol)
     mult is the multiplicity
   funct 2: improper dihedral
     c0 - th0 (deg)
     c1 - fc (kcal/mol/rad2)
   funct 3: RB dihedral
     c0-c5 - C0-C5 (kcal/mol) 
  */
  void getDihedralParams(int num, Real *c, int *mult, int *funct) const;

  // JLai
  void getPairLJArrays2(int *indexA, int *indexB, Real *pairC6, Real *pairC12);
  void getPairGaussArrays2(int *indexA, int *indexB, Real *gaussA, Real *gaussMu1,
  			  Real *gaussSigma1, Real *gaussMu2, Real *gaussSigma2,
  			   Real *gaussRepulsive);

  /* getVDWParams returs the Lennard-Jones bonding parameters for the
     specified two atom types, and the modified bonding parameters
     for 1-4 L-J interactons (c6pair, c12pair).
     Note that unlike the other types of parameters, you don't refer
     to this one by number - any two atom types define a VDWParam */
  void getVDWParams(int typea, int typeb,
		    Real *c6, Real *c12, Real *c6pair, Real *c7) const;
    
};

#endif // GROMACSTOPFILE_H
