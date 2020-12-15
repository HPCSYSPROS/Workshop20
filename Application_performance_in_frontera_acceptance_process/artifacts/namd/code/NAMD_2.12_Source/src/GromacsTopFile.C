#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "common.h"
#include "ResizeArray.h"
#include "InfoStream.h"
#include "GromacsTopFile.h"
#include "InfoStream.h"

#define JOULES_PER_CALORIE 4.184
#define ANGSTROMS_PER_NM 10

/* A GromacsTopFile represents the information stored in a GROMACS
   topolgy file.  This is an immutable type. */

/* XXX warning: this code contains a few algorithms which run in a
   time of O(n^3) or so in the size of the topology file, since I
   haven't bothered to do any sorting.  I don't think it matters much,
   since it still manages to run on the biggest files in less than a
   second, and nothing (memory or time) depends on the total number of
   atoms in the simulation - all that matters is the number that have
   to be individually defined. */

/* GromacsTopFile initializes this to represent the data stored in the
   file <filename>, or exits on error. */
#define LINESIZE 200

/* modes */
#define UNKNOWN 0
#define ATOMS 1
#define MOLECULETYPE 2
#define MOLECULES 3
#define SYSTEM 4
#define BONDS 5
#define BONDTYPES 6
#define ANGLES 7
#define ATOMTYPES 8
#define ANGLETYPES 9
#define DIHEDRALS 10
#define DIHEDRALTYPES 11
#define DEFAULTS 12
#define NONBOND 13
#define PAIRS 14
#define EXCLUSIONS 15
#ifndef CODE_REDUNDANT
#define CODE_REDUNDANT 0
#endif

#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

/* Initialize variables for exclusion calculation */
/* JLai */
ResizeArray<int> exclusions_atom_i;
ResizeArray<int> exclusions_atom_j;
int numExclusion = 0;
int numLJPair = 0;
int numGaussPair = 0;
bool bool_negative_number_warning_flag = false;



GromacsTopFile::GromacsTopFile(char *filename) {
  fudgeLJ = fudgeQQ = 1.0;
  /* open the file */
  FILE *f = fopen(filename,"r");
  char buf[LINESIZE];
  char modename[20];
  int mode;
  if(f==NULL) {
    sprintf(buf,"Error opening file '%s'",filename);
    NAMD_die(buf);
  }

  /* really bad parser XXX probably just works on the files we
     happen to have REWRITE THIS SOON.  It should allow for \- line
     continuations, check for errors in the file, etc. */
  while(fgets(buf,LINESIZE-1,f)) {
    char testchar;
    int i,j;

    /* defaults */
    int nbfunc, combrule;
    char genpairs[20];
    
    /* atom buffers */
    int num, resnum, chargegp, typenum;
    char type[NAMESIZE+1], restype[NAMESIZE+1], atomname[NAMESIZE+1];
    char particletype[NAMESIZE+1];
    float charge, mass, c6, c12, junkf;

    /* moltype buffers */
    int nrexcl;
    char molname[LONGNAMESIZE+1];

    /* molInst buffers */
    int copies;
    MolInst *moleculeinstance;

    /* general atomset buffers */
    int atomi, atomj, atomk, atoml;
    char typea[NAMESIZE+1],typeb[NAMESIZE+1],
      typec[NAMESIZE+1],typed[NAMESIZE+1];
    const char *tmptypea,*tmptypeb,*tmptypec,*tmptyped;
    int funct, index;
    float c0,c1;
    
    /* bonds */
    float b0,kB,th0,kth;

    /* dihedrals */
    float c[6];
    int mult=0;

    /* check for comments */
    if(sscanf(buf," %c",&testchar)==1) {
      if(testchar == ';') continue;
    }
    else { /* this is a blank line */
      continue;
    }

    /* check for a new mode */
    if(sscanf(buf," [ %19[^] ] ]",modename)==1) {
      /* switch the mode */
      if(0==strcmp(modename,"atoms"))              mode = ATOMS;
      else if(0==strcmp(modename,"atomtypes"))	   mode = ATOMTYPES;
      else if(0==strcmp(modename,"moleculetype"))  mode = MOLECULETYPE;
      else if(0==strcmp(modename,"molecules"))	   mode = MOLECULES;
      else if(0==strcmp(modename,"system"))	   mode = SYSTEM;
      else if(0==strcmp(modename,"bonds"))	   mode = BONDS;
      else if(0==strcmp(modename,"bondtypes"))	   mode = BONDTYPES;
      else if(0==strcmp(modename,"angles"))	   mode = ANGLES;
      else if(0==strcmp(modename,"angletypes"))	   mode = ANGLETYPES;
      else if(0==strcmp(modename,"dihedrals"))	   mode = DIHEDRALS;
      else if(0==strcmp(modename,"dihedraltypes")) mode = DIHEDRALTYPES;
      else if(0==strcmp(modename,"defaults"))      mode = DEFAULTS;
      else if(0==strcmp(modename,"nonbond_params")) mode = NONBOND;
      // JLai
      else if(0==strcmp(modename,"pairs")) mode = PAIRS;
      else if(0==strcmp(modename,"exclusions")) mode = EXCLUSIONS;
      else {	
	fprintf(stderr,"Warning: unknown mode %s\n",modename);
	mode = UNKNOWN;
      }

      continue;
    }

    /* now do the appropriate thing based on the current mode */
    switch(mode) {
    case SYSTEM:
      systemName = strdup(buf);
      break;

    case DEFAULTS:
      i = sscanf(buf," %d %d %20s %f %f",
		 &nbfunc,&combrule,genpairs,&fudgeLJ,&fudgeQQ);
      if(i < 3) { // didn't get enough parameters 
	fprintf(stderr,"syntax error in DEFAULTS\n");
	exit(1);
      }
      if(nbfunc != 1) { // I don't know how it works with nbfunc=2
	fprintf(stderr,"Non-bonded function != 1 unsupported in DEFAULTS\n");
	exit(1);
      }
      if(combrule != 1) { // same here
	fprintf(stderr,"Combination rule != 1 unsupported in DEFAULTS\n");
	exit(1);
      }
      if(0==strcmp(genpairs,"yes")) {
	genPairs=1;
	if(i!=5) {
	  fprintf(stderr,"syntax error in DEFAULTS\n");
	  exit(1);
	}
	// else fudgeLJ and fudgeQQ got written automatically
      }
      else genPairs=0;

      break;

    case NONBOND:
      if(5 != sscanf(buf," %5s %5s %d %f %f",
		     typea, typeb, &funct, &c6, &c12)) {
	fprintf(stderr,"Syntax error in NONBOND\n");
	exit(1);
      }
      // convert kJ/mol*nm6 to kcal/mol*A6 and ..12 ..12
      c6 =  c6/JOULES_PER_CALORIE*1E6;
      c12= c12/JOULES_PER_CALORIE*1E12;
      vdwTable.addType(typea,typeb,c6,c12);
      break;

    case BONDS:
      i = sscanf(buf," %d %d %d %f %f",
		 &atomi,&atomj,&funct,&c0,&c1);
      atomi--; // shift right away to a zero-indexing
      atomj--;
      if(i==3) {
       tmptypea = genericMols[genericMols.size()-1]->getAtom(atomi)->getType();
       tmptypeb = genericMols[genericMols.size()-1]->getAtom(atomj)->getType();
	/* find the index and parameters */
	index = bondTable.getParams(tmptypea, tmptypeb, funct, &b0, &kB);
	if(index==-1) {
	  fprintf(stderr,"Required bondtype %s--%s (function %d) not found.\n",
		  tmptypea,tmptypeb,funct);
	  exit(1);
	}
      }
      else if(i==5) {
	/* first set the values of b0 and kB correctly */
	b0 = c0*ANGSTROMS_PER_NM; /* convert nm to A */
	if(funct==1) { /* harmonic potential */
	  /* convert kJ/nm2 to kcal/A2 and use E=kx2 instead of half that. */
	  kB = c1/JOULES_PER_CALORIE/100/2;
	}
	else if(funct==2) { /* special fourth-order potential */
	  /* convert to the normal harmonic constant and kJ/nm2 to kcal/A2 */
	  kB = 2*c1*c0*c0/JOULES_PER_CALORIE/100;
	  kB /= 2; /* use the NAMD system where E=kx2 */
	}
	else {
	  fprintf(stderr,"I don't know what funct=%d means in BONDS\n",funct);
	  exit(1);
	}
	/* look up the index */
	index = bondTable.getIndex(b0,kB,funct);
      }
      else {
	fprintf(stderr,"Syntax error in BONDS\n");
	exit(1);
      }

      genericMols[genericMols.size()-1]->addBond(atomi,atomj,index);
      break;
      
    case BONDTYPES:
      if(5 != sscanf(buf," %5s %5s %d %f %f",
		     typea,typeb,&funct,&c0,&c1)) {
	fprintf(stderr,"Syntax error in BONDTYPES\n");
	exit(1);
      }

      /* first set the values of b0 and kB correctly */
      b0 = c0*10; /* convert nm to A */
      if(funct==1) { /* harmonic potential */
	/* convert kJ/nm2 to kcal/A2 and use E=kx2 instead of half that. */
	kB = c1/JOULES_PER_CALORIE/100/2;
      }
      else if(funct==2) { /* special fourth-order potential */
	/* convert to the normal harmonic constant and kJ/nm2 to kcal/A2 */
	kB = 2*c1*c0*c0/JOULES_PER_CALORIE/100;
	kB /= 2; /* use the NAMD system where E=kx2 */
      }
      else {
 fprintf(stderr,"I don't know what funct=%d means in BONDTYPES\n",funct);
	exit(1);
      }

      bondTable.addType(typea,typeb,b0,kB,funct);
      break;

    case ANGLES:
      i = sscanf(buf," %d %d %d %d %f %f",
		     &atomi,&atomj,&atomk,&funct,&c0,&c1);
      atomi--; // shift right away to a zero-indexing
      atomj--;
      atomk--;
      if(i == 4) { /* we have to look up the last two parameters */
       tmptypea = genericMols[genericMols.size()-1]->getAtom(atomi)->getType();
       tmptypeb = genericMols[genericMols.size()-1]->getAtom(atomj)->getType();
       tmptypec = genericMols[genericMols.size()-1]->getAtom(atomk)->getType();
	/* find the index and parameters */
	index = angleTable.getParams(tmptypea, tmptypeb, tmptypec,
					funct, &th0, &kth);
	if(index==-1) {
	  fprintf(stderr,
		  "Required angletype %s--%s--%s (function %d) not found.\n",
		  tmptypea,tmptypeb,tmptypec,funct);
	  exit(1);
	}
      }
      else if(i == 6) {
	/* first set the values of th0 and kth correctly */
	if(funct == 1) {
	  th0 = c0; /* both are in degrees */
	  kth = c1/JOULES_PER_CALORIE/2; /* convert kJ/rad2 to kcal/rad2 and use E=kx2 */
	}
	else if(funct == 2) {
	  th0 = c0; /* both are in degrees */
	  /* convert G96 kJ to kcal/rad2 and use E=kx2 */
	  kth = sin(th0*PI/180)*sin(th0*PI/180)*c1/JOULES_PER_CALORIE/2;
	}
	else {
	  fprintf(stderr,"I don't know what funct=%d means in ANGLES\n",funct);
	  exit(1);
	}
	/* add the angle type to our table */
	index = angleTable.getIndex(th0,kth,funct);
      }
      else {
	fprintf(stderr,"Syntax error (%d args) in ANGLES: %s\n",i,buf);
	exit(1);
      }

      /* add the angle to our table */
      genericMols[genericMols.size()-1]->addAngle(atomi,atomj,atomk,index);
      break;

    case ANGLETYPES:
      if(6 != sscanf(buf," %5s %5s %5s %d %f %f",
		     typea,typeb,typec,&funct,&c0,&c1)) {
	fprintf(stderr,"Syntax error in ANGLETYPES\n");
	exit(1);
      }
      /* first set the values of th0 and kth correctly */
      if(funct == 1) {
	th0 = c0; /* both are in degrees */
	kth = c1/JOULES_PER_CALORIE/2; /* convert kJ/rad2 to kcal/rad2 and use E=kx2 */
      }
      else if(funct == 2) {
	th0 = c0; /* both are in degrees */
	/* convert G96 kJ to kcal/rad2 and use E=kx2 */
	kth = sin(th0*PI/180)*sin(th0*PI/180)*c1/JOULES_PER_CALORIE/2;
      }
      else {
	fprintf(stderr,"I don't know what funct=%d means in ANGLETYPES\n",
		funct);
	exit(1);
      }
      angleTable.addType(typea,typeb,typec,th0,kth,funct);
      break;

    case DIHEDRALS:
      i = sscanf(buf," %d %d %d %d %d %f %f %f %f %f %f",
		 &atomi,&atomj,&atomk,&atoml,&funct,
		 &c[0],&c[1],&c[2],&c[3],&c[4],&c[5]);
      atomi--; // shift right away to a zero-indexing
      atomj--;
      atomk--;
      atoml--;
      if(i==5) { /* we have to look up the parameters */
       tmptypea = genericMols[genericMols.size()-1]->getAtom(atomi)->getType();
       tmptypeb = genericMols[genericMols.size()-1]->getAtom(atomj)->getType();
       tmptypec = genericMols[genericMols.size()-1]->getAtom(atomk)->getType();
       tmptyped = genericMols[genericMols.size()-1]->getAtom(atoml)->getType();
	/* find the index and parameters */
	index = dihedralTable.getParams(tmptypea, tmptypeb, tmptypec,
					tmptyped, funct, c, &mult);
	if(index==-1) {
	  fprintf(stderr,
	     "Required dihedraltype %s--%s--%s--%s (function %d) not found.\n",
		  tmptypea,tmptypeb,tmptypec,tmptyped,funct);
	  exit(1);
	}
      }
      else if(i==7 || i==8 || i==11) { /* the parameters are given */
	if(funct==1 || funct==2) { /* we should have two parameters */
	  if(i!=7+(funct==1)) {    /* plus a multiplicity for funct==1 */
	    fprintf(stderr,"Must have 7 args for funct=1,2\n");
	    exit(1);
	  }
	  c[0] = c[0]; /* both in deg */
	  if(i==7) {
	      c[1] = c[1]/(2*JOULES_PER_CALORIE); /* convert kJ to kcal and still use E=kx2*/
	  } else if (i==8 || i==11) {
	      c[1] = c[1]/(1*JOULES_PER_CALORIE); /* convert kJ to kcal and still use E=kx2*/
	  } 
	  //c[1] = c[1]/(2*JOULES_PER_CALORIE); /* convert kJ to kcal and still use E=kx2*/
	  /* for funct==1 these are both divided by rad^2 */
	  if(funct==1) {
	    mult=(int)(c[2]+0.5); /* round to nearest integer :p */
	  }
	}
	else if(funct==3) {
	  if(i!=11){
	    fprintf(stderr,"Must have 11 args for funct=3\n");
	    exit(1);
	  }

	  for(j=0;j<6;j++) {
	    c[j] = c[j]/JOULES_PER_CALORIE; /* convert kJ to kcal and E=Cn cos^n */
	  }
	}
	else {
	  fprintf(stderr,
		  "I don't know what funct=%d means in DIHEDRALS\n",funct);
	  exit(1);
	}
	index = dihedralTable.getIndex(c,mult,funct);
      }
      else {
	fprintf(stderr,"Syntax error (%d args) in DIHEDRALS: %s\n",i,buf);
	exit(1);
      }

      /* add the dihedrals to our table */
      genericMols[genericMols.size()-1]->addDihedral(atomi,atomj,atomk,atoml,
						     index);
      break;
    case DIHEDRALTYPES:
      i = sscanf(buf," %5s %5s %d %f %f %f %f %f %f",
		 typea,typeb,&funct,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5]);
      if(funct == 1 || funct == 2) {
	if(i!=5+(funct==1)) { /* 6 for f=2, 5 for f=1 */
	  fprintf(stderr,"Syntax error in DIHEDRALTYPES: %s\n",buf);
	  exit(1);
	}
	c[0] = c[0]; /* both in deg */
	c[1] = c[1]/JOULES_PER_CALORIE; /* convert kJ to kcal and still use E=kx2*/
	/* for funct==1 these are both divided by rad^2 */
	if(funct==1) {
	  mult=(int)(c[2]+0.5); /* round to nearest integer :p */
	}
      }
      else if(funct == 3) {
	if(i!=9) {
	  fprintf(stderr,"Syntax error in DIHEDRALTYPES\n");
	  exit(1);
	}
	for(j=0;j<6;j++) {
	  c[j] = c[j]/JOULES_PER_CALORIE; /* convert kJ to kcal and E=Cn cos^n */
	}
      }
      else {
	fprintf(stderr,"I don't know what funct=%d means in DIHEDRALTYPES\n",
		funct);
	exit(1);
      }
      dihedralTable.addType(typea,typeb,c,mult,funct);
      break;
      
    case ATOMS:
      i = sscanf(buf," %d %5s %d %5s %5s %d %f %f",
		   &num, type, &resnum, restype,
		 atomname, &chargegp, &charge, &mass);
      if(i==7) { /* XXX temporary - I should be able to get more
		    params */
        typenum = atomTable.getParams(type,&mass,&junkf,&junkf,&junkf);
	i=8;
      }
      else {
	if(i!=8) {
	  fprintf(stderr,"Syntax error in ATOMS\n");
	  exit(1);
	}
	// just get the type number
	typenum = atomTable.getParams(type,&junkf,&junkf,&junkf,&junkf);
      }
      genericMols[genericMols.size()-1]->addAtom(type,typenum,resnum,
						 restype,atomname,charge,mass);
      break;

    case ATOMTYPES:
      if(6 != sscanf(buf," %5s %f %f %5s %f %f",type,&mass,&charge,
		     particletype,&c6,&c12)) {
	fprintf(stderr,"Syntax error in ATOMTYPES: %s\n",buf);
	exit(1);
      }
      /* conversions:
	 c6  - kJ/mol nm6  -> kcal/mol A6
	 c12 - kJ/mol nm12 -> kcal/mol A12 */
      atomTable.addType(type,mass,charge,
			c6/(JOULES_PER_CALORIE)*1E6,
			c12/(JOULES_PER_CALORIE)*1E12);
      break;

    case MOLECULETYPE:
      if(2!=sscanf(buf," %20s %d",molname,&nrexcl)) {
	fprintf(stderr,"Syntax error in MOLECULETYPE: %s\n",buf);
	exit(1);
      }

      /* add another generic molecule holder */
      genericMols.add(new GenericMol(molname));
      break;

    case MOLECULES:
      if(2!=sscanf(buf," %20s %d",molname,&copies)) {
	fprintf(stderr,"Syntax error in MOLECULES: %s\n",buf);
	exit(1);
      }
      
      /* search for the specified molecule and make a molInst of it*/
      moleculeinstance = NULL;
      for(i=0;i<genericMols.size();i++) {
	if(0==strcmp(molname,genericMols[i]->getName())) {
	  moleculeinstance = new MolInst(genericMols[i],copies);
	  break;
	}
      }

      if(moleculeinstance==NULL) {
	fprintf(stderr,"Molecule %s undefined.\n",molname);
	exit(1);
      }
      
      /* put it at the end of the list */
      molInsts.add(moleculeinstance);

      break;
    case PAIRS:
      int indexA;
      int indexB;
      int pairFunction;
      Real fA;
      Real fB;
      Real fC;
      Real fD;
      Real fE;
      Real fF;
      
      int countVariables;
      countVariables = sscanf(buf," %d %d %d %f %f %f %f %f %f",&indexA,&indexB,&pairFunction,&fA,&fB,&fC,&fD,&fE,&fF);

      if ((countVariables >= 4 && countVariables >= 10)) {
	fprintf(stderr,"Syntax error in PAIRS: %s\n",buf);
	exit(1);
      }
      
      // Shift the atom indices to be zero-based
      indexA--;
      indexB--;

      // Make sure that the indexA is less than indexB
      /*if ( indexA > indexB ) {
	  int tmpIndex = indexA;
	  indexB = indexA;
	  indexA = tmpIndex;
	  }*/

      if (pairFunction == 1) {
	
	// LJ code
	fA = (fA/JOULES_PER_CALORIE)*1E6;
	fB= (fB/JOULES_PER_CALORIE)*1E12;
	pairTable.addPairLJType2(indexA,indexB,fA,fB);
      } else if (pairFunction == 5){
	
	// Bare Gaussian potential
	fA = (fA/JOULES_PER_CALORIE); //-->gaussA
	fB = (fB*ANGSTROMS_PER_NM); //-->gaussMu1
	if(fC == 0) {
	  char buff[100];
	  sprintf(buff,"GromacsTopFile.C::Attempting to load zero into gaussSigma.  Please check the pair: %s\n",buf);
	  NAMD_die(buff);
	}
	if(fC < 0 && !bool_negative_number_warning_flag) {
	  iout << iWARN << "Attempting to load a negative standard deviation into the gaussSigma.  Taking the absolute value of the standard deviation.";
	  bool_negative_number_warning_flag = true;
	}
	fC = (fC*ANGSTROMS_PER_NM); //-->gaussSigma1
	fC = 1.0/(2 * fC * fC); // Normalizes sigma
	pairTable.addPairGaussType2(indexA,indexB,fA,fB,fC);
      } else if (pairFunction == 6) {
	
	// Combined Gaussian + repulsive r^-12 potential
	fA = (fA/JOULES_PER_CALORIE); //-->gaussA
	fB = (fB*ANGSTROMS_PER_NM); //-->gaussMu1
	if(fC == 0) {
	  char buff[100];
	  sprintf(buff,"GromacsTopFile.C::Attempting to load zero into gaussSigma.  Please check the pair: %s\n",buf);
	  NAMD_die(buff);
	}
	if(fC < 0 && !bool_negative_number_warning_flag) {
	  iout << iWARN << "Attempting to load a negative standard deviation into the gaussSigma.  Taking the absolute value of the standard deviation.";
	  bool_negative_number_warning_flag = true;
	}
	fC = (fC*ANGSTROMS_PER_NM); //-->gaussSigma1
	fC = 1.0/(2 * fC * fC); // Normalizes sigma
	fD = (fD*ANGSTROMS_PER_NM); //-->gaussRepulsive
	pairTable.addPairGaussType2(indexA,indexB,fA,fB,fC,fD);
      } else if (pairFunction == 7) {
	
	// Double well Guassian function
	fA = (fA/JOULES_PER_CALORIE); //-->gaussA
	fB = (fB*ANGSTROMS_PER_NM); //-->gaussMu1
        if(fC == 0 || fE == 0) {
	  char buff[100];
	  sprintf(buff,"GromacsTopFile.C::Attempting to load zero into gaussSigma.  Please check the pair: %s\n",buf);
	  NAMD_die(buff);
	}
	if((fC < 0 || fE < 0)&& !bool_negative_number_warning_flag) {
	  iout << iWARN << "Attempting to load a negative standard deviation into the gaussSigma.  Taking the absolute value of the standard deviation.";
	  bool_negative_number_warning_flag = true;
	}
	fC = (fC*ANGSTROMS_PER_NM); //-->gaussSigma1
	fC = 1.0/(2 * fC * fC); // Normalizes sigma
        fD = (fD*ANGSTROMS_PER_NM); //-->gaussMu2
	fE = (fE*ANGSTROMS_PER_NM); //-->gaussSigma2
	fE = 1.0/(2 * fE * fE); // Normalizes sigma
	fF = (fE*ANGSTROMS_PER_NM); //-->gaussRepulsive
	pairTable.addPairGaussType2(indexA,indexB,fA,fB,fC,fD,fE,fF);
      } else {
	
	// Generic error statement
	fprintf(stderr,"Unknown pairFunction in GromacsTopFile.C under the PAIRS section: %d\n",pairFunction);
      }
      break;
    case EXCLUSIONS:
      /* Start of JLai modifications August 16th, 2012 */
      if(2 != sscanf(buf," %d %d ",&atomi,&atomj)) {
	fprintf(stderr,"Syntax error in EXCLUSIONS: %s\n",buf);
	exit(1);
      }
      
      // Shift the atom indices to be zero-based
      atomi--;
      atomj--;

      /*Load exclusion information into file*/
      exclusions_atom_i.add(atomi);
      exclusions_atom_j.add(atomj);
      numExclusion++;

      /* Reading in exclusions information from file and loading */
      break;
      // End of JLai modifications August 16th, 2012
    }
  }

  fclose(f);
}

/* returns the number of bonds in the file */
int GromacsTopFile::getNumBonds() const {
  int n=0,i;
  for(i=0;i<molInsts.size();i++) {
    n += molInsts[i]->getNum() *
         molInsts[i]->getMol()->getNumBonds();
  }
  return n;
}

/* returns the number of angles in the file */
int GromacsTopFile::getNumAngles() const {
  int n=0,i;
  for(i=0;i<molInsts.size();i++) {
    n += molInsts[i]->getNum() *
         molInsts[i]->getMol()->getNumAngles();
  }
  return n;
}

/* returns the number of dihedral angles in the file */
int GromacsTopFile::getNumDihedrals() const {
  int n=0,i;
  for(i=0;i<molInsts.size();i++) {
    n += molInsts[i]->getNum() *
         molInsts[i]->getMol()->getNumDihedrals();
  }
  return n;
}

/* JLai -- August 16th, 2012 modifications*/
/* returns the number of pair interactions in the file */
int GromacsTopFile::getNumPair() const {
  int numPair = 0;
  numPair = numLJPair + numGaussPair;
  return numPair;
}

int GromacsTopFile::getNumLJPair() const {
  return numLJPair;
}

int GromacsTopFile::getNumGaussPair() const {
  return numGaussPair;
}

/* return the number of exclusions in the file */
int GromacsTopFile::getNumExclusions() const {
  return numExclusion;
}

/* return the list of exclusions from the file */
void GromacsTopFile::getExclusions(int* atomi, int* atomj) const {
    for(int i =0; i < exclusions_atom_i.size(); i++) {
	atomi[i] = exclusions_atom_i[i];
	atomj[i] = exclusions_atom_j[i];
    }
  return;
}
/* End of JLai modifications */

/* getBond puts the information about bond number <num> into the
   spaces pointed to by the other arguments.  Bond number 0 is the
   first bond in the list.
   For atomi and atomj, 1 is the first atom in the list. */
void GromacsTopFile::getBond(int n, int *atomi, int *atomj, int *bondtype) const {
  /* figure out what instance this bond is in */
  int nbonds=0;
  int natoms=0;
  int i;
  for(i=0;i<molInsts.size();i++) {
    int molbonds = molInsts[i]->getMol()->getNumBonds();
    int molatoms = molInsts[i]->getMol()->getNumAtoms();
    int newbonds = molInsts[i]->getNumBonds();
    int newatoms = molInsts[i]->getNumAtoms();
    
    if(nbonds+newbonds-1 >= n) {
      /* the bond is in this MolInst */
      int localnum = (n-nbonds) % molbonds;   /* number within this inst */
      int instnum = (n-nbonds) / molbonds;    /* number of instances before */
      int addatoms = natoms+instnum*molatoms; /* extra atoms to add to atomi */

      const GenericBond *b = molInsts[i]->getMol()->getBond(localnum);

      *atomi = b->getAtomi() + addatoms;
      *atomj = b->getAtomj() + addatoms;
      *bondtype = b->getType();
      break;
    }

    nbonds += newbonds;
    natoms += newatoms;
  }
}

/* getAngle puts the information about angle number <num> into the
   spaces pointed to by the other arguments.  Angle number 0 is the
   first angle in the list.
*/
void GromacsTopFile::getAngle(int n, int *atomi, int *atomj, int *atomk,
			  int *angletype) const {
  /* figure out what instance this angle is in */
  int nangles=0;
  int natoms=0;
  int i;
  for(i=0;i<molInsts.size();i++) {
    int molangles = molInsts[i]->getMol()->getNumAngles();
    int molatoms = molInsts[i]->getMol()->getNumAtoms();
    int newangles = molInsts[i]->getNumAngles();
    int newatoms = molInsts[i]->getNumAtoms();
    
    if(nangles+newangles-1 >= n) {
      /* the angle is in this MolInst */
      int localnum = (n-nangles) % molangles; /* number within this inst */
      int instnum = (n-nangles) / molangles;  /* number of instances before */
      int addatoms = natoms+instnum*molatoms; /* extra atms to add to atmi */

      const GenericAngle *a = molInsts[i]->getMol()->getAngle(localnum);

      *atomi = a->getAtomi() + addatoms;
      *atomj = a->getAtomj() + addatoms;
      *atomk = a->getAtomk() + addatoms;
      *angletype = a->getType();
      break;
    }

    nangles += newangles;
    natoms += newatoms;
  }
}

/* getDihedral puts the information about dihedral number <num> into
   the spaces pointed to by the other arguments.  Dihedral number 0
   is the first angle in the list. */
void GromacsTopFile::getDihedral(int n, int *atomi, int *atomj, int *atomk,
			      int *atoml, int *type) const {
  /* figure out what instance this angle is in */
  int ndihedrals=0;
  int natoms=0;
  int i;
  for(i=0;i<molInsts.size();i++) {
    int moldihedrals = molInsts[i]->getMol()->getNumDihedrals();
    int molatoms = molInsts[i]->getMol()->getNumAtoms();
    int newdihedrals = molInsts[i]->getNumDihedrals();
    int newatoms = molInsts[i]->getNumAtoms();
    
    if(ndihedrals+newdihedrals-1 >= n) {
      /* the dihedral is in this MolInst */
      int localnum = (n-ndihedrals) % moldihedrals; /* number in this inst */
      int instnum = (n-ndihedrals) / moldihedrals; /* number of insts before */
      int addatoms = natoms+instnum*molatoms; /*extra atms to add to atmi*/

      const GenericDihedral *a = molInsts[i]->getMol()->getDihedral(localnum);

      *atomi = a->getAtomi() + addatoms;
      *atomj = a->getAtomj() + addatoms;
      *atomk = a->getAtomk() + addatoms;
      *atoml = a->getAtoml() + addatoms;
      *type = a->getType();
      break;
    }

    ndihedrals += newdihedrals;
    natoms += newatoms;
  }
}

/* getNumAtoms returns the number of atoms stored in the file. */
int GromacsTopFile::getNumAtoms() const {
  int n=0;
  int i;
  for(i=0;i<molInsts.size();i++) {
    n += molInsts[i]->getNum() *
         molInsts[i]->getMol()->getNumAtoms();
  }
  return n;
}

/* getAtom puts the information about atom number <n> into the
   spaces pointed to by the other arguments.  The string buffers
   must be at least 11 characters long. */
void GromacsTopFile::getAtom(int num,
	       int *residue_number, char *residue_name,
	       char *atom_name, char *atom_type, int *atom_typenum,
	       Real *charge, Real *mass) 
  const {
  int natoms=0,n=num; // zero-indexed array
  int resnum=0;
  /* figure out what instance the atom is in, and what residue number
     it has */
  int i;

  for(i=0;i<molInsts.size();i++) {
    int numnew = molInsts[i]->getNumAtoms(); /* atoms in this MolInst */
    int resmol = molInsts[i]->getMol()->getNumRes(); /* # res/mol */
    int newres = molInsts[i]->getNumRes(); /* residues in this MolInst */

    if(natoms+numnew-1 >= n) {
      /* the atom is in this molInst */
      int localnum = (n-natoms) % molInsts[i]->getMol()->getNumAtoms();
      int instnum  = (n-natoms) / molInsts[i]->getMol()->getNumAtoms();
      
      // getAtom is zero-indexed
      const GenericAtom *a = molInsts[i]->getMol()->getAtom(localnum);
      int residue = resnum + resmol*instnum + a->getResNum();

      *residue_number = residue;
      strncpy(residue_name,a->getResType(),11);
      strncpy(atom_name,a->getAtomName(),11);
      strncpy(atom_type,a->getType(),11);
      *charge=a->getCharge();
      *mass=a->getMass();
      *atom_typenum=a->getTypeNum();
      break;
    }

    /* go to the next molInst */
    natoms += numnew;
    resnum += newres;
  }
}

GenericAtom::GenericAtom(const char *theType, int theTypeNum, int theResNum,
			 const char *theResType,const char *theAtomName,
			 Real theCharge, Real theMass) {
  strncpy(type,theType,NAMESIZE+1);
  typeNum = theTypeNum;
  resNum = theResNum;
  strncpy(resType,theResType,NAMESIZE+1);
  strncpy(atomName,theAtomName,NAMESIZE+1);
  charge = theCharge;
  mass = theMass;
}

/* initializes this to be a bond between atoms <i> and <j> of type
   <type>, with <next> pointing to the next bond in the list. */
GenericBond::GenericBond(int i, int j, int theType) {
  atomi=i;
  atomj=j;
  type=theType;
}

/* initializes this to be a angle between atoms <i>, <j>, and <k> of
   type <type>, with <next> pointing to the next angle in the list. */
GenericAngle::GenericAngle(int i, int j, int k, int theType) {
  atomi=i;
  atomj=j;
  atomk=k;
  type=theType;
}

/* initializes this to be a angle between atoms <i>, <j>, <k>, and
   <l> of type <type> */
GenericDihedral::GenericDihedral(int i, int j, int k, int l, int theType) {
  atomi=i;
  atomj=j;
  atomk=k;
  atoml=l;
  type=theType;
}

/* adds a bond to the list */
void GenericMol::addBond(int atomi, int atomj, int type) {
  bondList.add(new GenericBond(atomi, atomj, type));
}

/* adds an angle to the list */
void GenericMol::addAngle(int atomi, int atomj, int atomk, int type) {
  angleList.add(new GenericAngle(atomi, atomj, atomk, type));
}

/* adds a dihedral to the list */
void GenericMol::addDihedral(int atomi, int atomj, int atomk,
			     int atoml, int type) {
  dihedralList.add(new GenericDihedral(atomi, atomj, atomk,
					     atoml, type));
}

/* adds an atom to the list */
void GenericMol::addAtom(const char *theType, int theTypeNum, int theResNum,
			 const char *theResType,const char *theAtomName,
			 Real theCharge, Real theMass) {

  atomList.add(new GenericAtom(theType,theTypeNum,theResNum,theResType,
			   theAtomName,theCharge,theMass));
}

/* initializes this to be the molecule with name <name> */
GenericMol::GenericMol(const char *theName) {
  name = strdup(theName);
}

/* gets a bond */
const GenericBond *GenericMol::getBond(int n) const {
  /* double-check */
  if(n >= bondList.size() || n<0) {
    fprintf(stderr,"Bond index %d out of bounds.\n",n);
    exit(1);
  }
  return bondList[n];
}

/* gets a angle */
const GenericAngle *GenericMol::getAngle(int n) const {
  /* double-check */
  if(n >= angleList.size() || n<0) {
    fprintf(stderr,"Angle index %d out of bounds.\n",n);
    exit(1);
  }
  return angleList[n];
}

/* gets a dihedral */
const GenericDihedral *GenericMol::getDihedral(int n) const {
  /* double-check */
  if(n >= dihedralList.size() || n<0) {
    fprintf(stderr,"Dihedral index %d out of bounds.\n",n);
    exit(1);
  }
  return dihedralList[n];
}

/* gets an atom */
const GenericAtom *GenericMol::getAtom(int n) const {
  /* double-check */
  if(n >= atomList.size() || n<0) {
    fprintf(stderr,"Atom index %d out of bounds for %s.\n",n,name);
    exit(1);
  }
  return atomList[n];
}

/* initializes this to represent <theNum> copies of <theMol> */
MolInst::MolInst(const GenericMol *theMol, int theNum) {
  mol = theMol;
  num = theNum;
}

/* get the total number of various things here */
int MolInst::getNumAtoms() const {
  return getMol()->getNumAtoms() * getNum();
}
int MolInst::getNumBonds() const {
  return getMol()->getNumBonds() * getNum();
}
int MolInst::getNumAngles() const {
  return getMol()->getNumAngles() * getNum();
}
int MolInst::getNumDihedrals() const {
  return getMol()->getNumDihedrals() * getNum();
}
int MolInst::getNumRes() const {
  return getMol()->getNumRes()
    * getNum();
}

/* this gets the index of a bond in the table (adding an entry if
   none exists).
   b0: the natural bond length in A.
   kB: the spring constant in kcal/mol/A^2, where E=kx^2 to 1st order.
   funct: 1 for harmonic potential, 2 for fourth-order GROMACS96
   potential. */
int BondTable::getIndex(float b0, float kB, int funct) {
  /* check to see if it is in the table already */
  int i;
  for(i=0;i<b0Array.size();i++) {
    if(fabs(b0-b0Array[i])<0.00001 &&
       fabs(kB-kBArray[i])<0.00001 &&
       funct == functArray[i]) {
      return i;
    }
  }
  
  /* nope, it wasn't in the table add a new element! */
  b0Array.add(b0);
  kBArray.add(kB);
  functArray.add(funct);
  typeaArray.add(NULL);
  typebArray.add(NULL);
  return b0Array.size()-1;
}

/* this gets the index of a angle in the table (adding an entry if
   none exists).
   b0: the natural angle length in A.
   kB: the spring constant in kcal/mol/A^2, where E=kx^2 to 1st order.
   funct: 1 for harmonic potential, 2 for fourth-order GROMACS96
   potential. */
int AngleTable::getIndex(float th0, float kth, int funct) {
  /* check to see if it is in the table already */
  int i;
  for(i=0;i<th0Array.size();i++) {
    if(fabs(th0-th0Array[i])<0.00001 &&
       fabs(kth-kthArray[i])<0.00001 &&
       funct == functArray[i]) {
      return i;
    }
  }
  
  /* nope, it wasn't in the table add a new element! */
  th0Array.add(th0);
  kthArray.add(kth);
  functArray.add(funct);
  typeaArray.add(NULL);
  typebArray.add(NULL);
  typecArray.add(NULL);
  return th0Array.size()-1;
}

/* this function gets the index of a dihedral angle in the table
   (adding an entry if none exists).  The required number of
   parameters (see notes on class DihedralTable) must be stored in
   the array <c> */
int DihedralTable::getIndex(const float *c, int mult, int funct) {
  /* check to see if it is in the table already */
  int i,j,jmax;
  if(funct==1 || funct==2) { /* for these we only need two params */
    jmax=2;
  } else { /* for RB we need all 6 */
    jmax=6;
  }

  for(i=0;i<cArray.size();i++) {
    int mismatch=0;
    if(mult != multArray[i]) continue;
    if(funct != functArray[i]) continue;

    for(j=0;j<jmax;j++) {
      if(fabs(c[j]-cArray[i][j])>=0.00001) {
	mismatch=1;
	break;
      }
    }
    if(!mismatch) {
      /* all of the parameters matched */
      return i;
    }
  }
  
  /* nope, it wasn't in the table add a new element! */
  addType(NULL,NULL,c,mult,funct);
  return cArray.size()-1;
}

/* getBondParams puts the parameters for bond-type <num> into the
   spaces pointed to by the other arguments. 
   b0 - natural length in A
   kB - spring constant in kcal/A^2 - E=kx^2
   funct - 1 for harmonic, 2 for special fourth-order GROMOS96 thing */
void GromacsTopFile::getBondParams(int num, float *b0, float *kB, int *funct)
  const {
  bondTable.getParams(num,b0,kB,funct);
}

/* getAngleParams puts the parameters for angle-type <num> into the
   spaces pointed to by the other arguments.  */
void GromacsTopFile::getAngleParams(int num, float *th0, float *kth, int *funct)
  const {
  angleTable.getParams(num,th0,kth,funct);
}

void GromacsTopFile::getDihedralParams(int num, float *c, int *mult, int *funct)
  const {
  dihedralTable.getParams(num,c,mult,funct);
}

void GromacsTopFile::getPairLJArrays2(int *indexA, int *indexB, Real *pairC6, Real *pairC12) {
  pairTable.getPairLJArrays2(indexA, indexB, pairC6, pairC12);
}

void GromacsTopFile::getPairGaussArrays2(int *indexA, int *indexB, Real *gaussA, Real *gaussMu1,
				    Real *gaussSigma1, Real *gaussMu2, Real *gaussSigma2,
				    Real *gaussRepulsive) {
  pairTable.getPairGaussArrays2(indexA, indexB, gaussA, gaussMu1, gaussSigma1, gaussMu2, gaussSigma2, gaussRepulsive);
}

/* getParams puts the parameters for bond-type <num> into the
   spaces pointed to by the other arguments. 
   b0 - natural length in A
   kB - spring constant in kcal/A^2 - E=kx^2
   funct - 1 for harmonic, 2 for special fourth-order GROMOS96 thing */
void BondTable::getParams(int num, float *b0, float *kB, int *funct)
  const {
  *b0=b0Array[num];
  *kB=kBArray[num];
  *funct=functArray[num];
}

/* getParams puts the parameters for angle-type <num> into the
   spaces pointed to by the other arguments.  */
void AngleTable::getParams(int num, float *th0, float *kth, int *funct)
  const {
  *th0=th0Array[num];
  *kth=kthArray[num];
  *funct=functArray[num];
}

/* getParams puts the parameters for angle-type <num> into the
   spaces pointed to by the other arguments.  The required number of
   parameters (see notes on class DihedralTable) will be stored in
   the array <c>, so make sure that c has size >= 6 */
void DihedralTable::getParams(int num, float *c, int *mult, int *funct) const {
  int i;

  *funct=functArray[num]; /* first set the function */

  if(*funct==1 || *funct==2) { /* for these we only need two params */
    c[0]=cArray[num][0];
    c[1]=cArray[num][1];
  } else if(*funct==3) { /* for RB we need all 6 */
    for(i=0;i<6;i++) {
      c[i]=cArray[num][i];
    }
  } else { /* error */
    fprintf(stderr,"Bad function number %d - don't know what to do!\n",
	    *funct);
    exit(1);
  }

  if(*funct==1) { /* return the multiplicity */
    *mult=multArray[num];
  }
}

/* this adds a entry for a angle type to the table, including the
   three atoms involved in the angle */
void AngleTable::addType(const char *typea, const char *typeb,
		  const char *typec, float th0, float kth, int funct) {
  typeaArray.add(strdup(typea));
  typebArray.add(strdup(typeb));
  typecArray.add(strdup(typec));
  th0Array.add(th0);
  kthArray.add(kth);
  functArray.add(funct);
}

/* this adds a entry for a angle type to the table, including the
   two atoms involved in the angle.  The required number of
   parameters (see notes on class DihedralTable) must be stored in
   the array <c>.  Note that these two angles really refer to either
   atoms A and D or B and C depending on the dihedral type.  */
void DihedralTable::addType(const char *typea, const char *typeb,
			 const float *c, int mult, int funct) {
  float *cNew;
  int i;

  if(typea != NULL) typeaArray.add(strdup(typea));
  if(typeb != NULL) typebArray.add(strdup(typeb));
  functArray.add(funct);

  if(funct==1) multArray.add(mult);
  else multArray.add(0);

  if(funct==1 || funct==2) { /* for these we only need two params */
    cNew = new float[2];
    cNew[0]=c[0];
    cNew[1]=c[1];
  } else if(funct==3) { /* for RB we need all 6 */
    cNew = new float[6];
    for(i=0;i<6;i++) {
      cNew[i]=c[i];
    }
  } else { /* error */
    fprintf(stderr,"Bad function number %d - don't know what to do!\n",
	    funct);
    exit(1);
  }
  
  /* now push the actual parameters */
  cArray.add(cNew);
}

/* This version looks up the angle by atom type - the direction of
   the types doesn't matter (it can be A--B--C or C--B--A)!  If the
   specified angle/function combination cannot be found, it returns
   -1, otherwise it returns the index of the angletype. */
int AngleTable::getParams(const char *typea, const char *typeb,
		   const char *typec, int funct,
		   float *th0, float *kth) const {
  int i;
  for(i=0;i<th0Array.size();i++) {
    if(typeaArray[i] == NULL || typebArray[i] == NULL
       || typecArray[i] == NULL) continue;
    if( (0==strcmp(typea,typeaArray[i]) &&  /* A--B--C */
	 0==strcmp(typeb,typebArray[i]) &&
	 0==strcmp(typec,typecArray[i]))    /* or */
     || (0==strcmp(typec,typeaArray[i]) &&
	 0==strcmp(typeb,typebArray[i]) &&  /* C--B--A */
	 0==strcmp(typea,typecArray[i]))
     && funct==functArray[i]) {
      *th0 = th0Array[i];
      *kth = kthArray[i];
      return i;
    }
  }
  return -1;
}

/* see the (long) notes in the prototype definition */
int DihedralTable::getParams(const char *typea, const char *typeb,
			     const char *typec, const char *typed,
			     int funct, float *c, int *mult) const {
  int i,j;
  const char *comparea, *compareb;
  
  if(funct == 1 || funct == 3) { /* for these, we use the inner atoms */
    comparea = typeb;
    compareb = typec;
  } else { /* use the outer atoms */
    comparea = typea;
    compareb = typed;
  }

  for(i=0;i<cArray.size();i++) {
    if(typeaArray[i] == NULL || typebArray[i] == NULL)
      continue; /* no atom types specified */
    if(functArray[i] != funct)
      continue; /* wrong function type */

    if( (0==strcmp(comparea,typeaArray[i]) &&  /* A--B */
	 0==strcmp(compareb,typebArray[i]))    /*  or  */
     || (0==strcmp(compareb,typeaArray[i]) &&
	 0==strcmp(comparea,typebArray[i]))    /* B--A */
	  ) {
      if(funct==1 || funct==2) { /* for these we only need two params */
	c[0]=cArray[i][0];
	c[1]=cArray[i][1];
	if(funct==1) {
	  *mult = multArray[i];
	}
      } else if(funct==3) { /* for RB we need all 6 */
	for(j=0;j<6;j++) {
	  c[j]=cArray[i][j];
	}
      }
      return i;
    }
  }
  return -1;
}

/* this adds a entry for a bond type to the table, including the two
   atoms involved in the bond */
void BondTable::addType(const char *typea, const char *typeb,
			    float b0, float kB, int funct) {
  b0Array.add(b0);
  kBArray.add(kB);
  functArray.add(funct);
  typeaArray.add(strdup(typea));
  typebArray.add(strdup(typeb));
}

/* This version looks up the bond by atom type - the order of the
   types doesn't matter! */
int BondTable::getParams(const char *typea, const char *typeb,
			      int funct, float *b0, float *kB) const {
  int i;
  for(i=0;i<b0Array.size();i++) {
    if(typeaArray[i] == NULL || typebArray[i] == NULL) continue;
    if( (0==strcmp(typea,typeaArray[i]) &&
	 0==strcmp(typeb,typebArray[i]))
     || (0==strcmp(typeb,typeaArray[i]) &&
	 0==strcmp(typea,typebArray[i]))
     && funct==functArray[i]) {
      *b0 = b0Array[i];
      *kB = kBArray[i];
      return i;
    }
  }
  return -1;
}

/* this adds a entry for an atom type to the table */
void AtomTable::addType(const char *type, float m, float q,
			    float c6, float c12) {
  typeArray.add(strdup(type));
  mArray.add(m);
  qArray.add(q);
  c6Array.add(c6);
  c12Array.add(c12);
}

/* This looks up the atom type by number, returning it in the string
   <type>, which must be at least 11 characters long. */
void AtomTable::getType(int num, char *type) const {
  if(num>=mArray.size() || num<0) {
    fprintf(stderr,"atomParam index %d out of bounds!\n",num);
    exit(1);
  }
  strncpy(type,typeArray[num],NAMESIZE+1);
}

/* This looks up the atom by type - if it cannot be found, this
   returns -1, otherwise this returns the index of the atomtype (for
   consistency - this is really a useless number.) */
int AtomTable::getParams(const char *type, float *m, float *q,
		  float *c6, float *c12) const {
  int i;
  for(i=0;i<mArray.size();i++) {
    if(typeArray[i]==NULL) {
      fprintf(stderr,"Found NULL atom type in list.\n");
      exit(1);
    }
    if(0==strcmp(typeArray[i],type)) {
      *m = mArray[i];
      *q = qArray[i];
      *c6 = c6Array[i];
      *c12 = c12Array[i];
      return i;
    }
  }
  return -1;
}

/* finds the index from the two interacting types, or returns -1 */
int VDWTable::getIndex(const char *typea, const char *typeb) const {
  int i;
  for(i=0;i<c6Array.size();i++) {
    if((0==strcmp(typea,typeAArray[i]) &&
	0==strcmp(typeb,typeBArray[i])) ||
       (0==strcmp(typeb,typeAArray[i]) &&
	0==strcmp(typea,typeBArray[i]))) {
      return i;
    }
  }
  return -1;
}

/* these add VDW parameters to the list */
void VDWTable::addType(const char *typea, const char *typeb, Real c6, 
		       Real c12) {
  int i;

  /* check to see if the pair is already in the table */
  i = getIndex(typea,typeb);
  if(i != -1) {  /* it was in the table */
    c6Array[i] = c6;
    c12Array[i] = c12;
  }
  else { /* it wasn't in the table - add it! */
    typeAArray.add(strdup(typea));
    typeBArray.add(strdup(typeb));
    c6Array.add(c6);
    c12Array.add(c12);
    c6PairArray.add(0);
    c12PairArray.add(0);
  }
}

void VDWTable::add14Type(const char *typea, const char *typeb,
	       Real c6pair, Real c12pair) {
  int i;

  /* check to see if the pair is already in the table */
  i = getIndex(typea,typeb);
  if(i != -1) {  /* it was in the table */
    c6PairArray[i] = c6pair;
    c12PairArray[i] = c12pair;
  }
  else { /* it wasn't in the table - add it! */
    typeAArray.add(strdup(typea));
    typeBArray.add(strdup(typeb));
    c6PairArray.add(c6pair);
    c12PairArray.add(c12pair);
    c6Array.add(0);
    c12Array.add(0);
  }
}

/* this returns the VDW interaction parameters that were added to
   the list earlier, and returns the index or the parameters (a
   useless number) or -1 if it can't find the specified types. */
int VDWTable::getParams(const char *typea, const char *typeb,
	      Real *c6, Real *c12, Real *c6pair, Real *c12pair) const {
  int i;
  /* check to see if the pair is already in the table */
  i = getIndex(typea,typeb);
  if(i != -1) {  /* it was in the table - return the parameters  */
    *c6 = c6Array[i];
    *c12 = c12Array[i];
    *c6pair = c6PairArray[i];
    *c12pair = c12PairArray[i];
  }
  return i;
}

/* getVDWParams returs the Lennard-Jones bonding parameters for the
   specified two atom types, and the modified bonding parameters
   for 1-4 L-J interactons (c6pair, c12pair). */
void GromacsTopFile::getVDWParams(int numa, int numb,
		  Real *c6, Real *c12, Real *c6pair, Real *c12pair) const {
  int i,ret;
  Real c6a,c6b,c12a,c12b, mjunk, qjunk;
  char typea[11]="",typeb[11]="";

  /* get the string names corresponding to numa and numb */
  getAtomParams(numa,typea);
  getAtomParams(numb,typeb);
  if(typea[0]==0) { /* found a bug in my use of strncpy here once */
    fprintf(stderr,"Failed to get name of atom %d\n",numa);
    exit(1);
  }
  if(typeb[0]==0) {
    fprintf(stderr,"Failed to get name of atom %d\n",numb);
    exit(1);
  }

  /* first try - use the VDW table */
  i = vdwTable.getParams(typea,typeb,c6,c12,c6pair,c12pair);

  if(i==-1) {
    // QQ151069
    /*if ( !genPairs && numa != numb ) {
      iout << iWARN << "VDW table using combining rule for types "
      << typea << " " << typeb << "\n" << endi;
      }*/

    /* fallback - use the individual atom's parameters */
    ret = atomTable.getParams(typea, &mjunk, &qjunk, &c6a, &c12a);
    if(ret==-1) {
      fprintf(stderr,"Couldn't find atom type %s\n",typea);
      exit(1);
    }
    ret = atomTable.getParams(typeb, &mjunk, &qjunk, &c6b, &c12b);
    if(ret==-1) {
      fprintf(stderr,"Couldn't find atom type %s\n",typeb);
      exit(1);
    }
    
    *c6  = (float)(sqrt(c6a * c6b));
    *c12 = (float)(sqrt(c12a*c12b));

    /*  this is the only reasonable option  */
    *c6pair  = *c6  * fudgeLJ;
    *c12pair = *c12 * fudgeLJ;

  }

}

/* Start of JLai modifications August 16th, 2012 */
int PairTable::addPairLJType2(int indexA, int indexB, Real pairC6, Real pairC12) {
  GroLJPair glp;
  glp.indxLJA = indexA;
  glp.indxLJB = indexB;
  glp.c6pair = pairC6;
  glp.c12pair = pairC12;
  pairlistLJ.add(glp);
  numLJPair++;
  return 0;

  // Insert the second copy
  GroLJPair glp2;
  glp2.indxLJA = indexB;
  glp2.indxLJB = indexA;
  glp2.c6pair = pairC6;
  glp2.c12pair = pairC12;
  pairlistLJ.add(glp2);
  numLJPair++;
  return 0;
}

int PairTable::addPairGaussType2(int indexA, int indexB, Real gaussA, Real gaussMu1,
				Real gaussSigma1) {
  return addPairGaussType2(indexA, indexB, gaussA, gaussMu1, gaussSigma1, 0.0, 0.0, 0.0);
}

int PairTable::addPairGaussType2(int indexA, int indexB, Real gaussA, Real gaussMu1,
			     Real gaussSigma1, Real gaussRepulsive) {
  return addPairGaussType2(indexA, indexB, gaussA, gaussMu1, gaussSigma1, 0.0, 0.0, gaussRepulsive);
}

int PairTable::addPairGaussType2(int indexA, int indexB, Real gaussA, Real gaussMu1,
			     Real gaussSigma1, Real gaussMu2, Real gaussSigma2, 
			     Real gaussRepulsive) {
  GroGaussPair ggp;
  ggp.indxGaussA = indexA;
  ggp.indxGaussB = indexB;
  ggp.gA = gaussA;
  ggp.gMu1 = gaussMu1;
  ggp.giSigma1 = gaussSigma1;
  ggp.gMu2 = 0.0;
  ggp.giSigma2 = 0.0;
  ggp.gRepulsive = 0.0;
  pairlistGauss.add(ggp);
  numGaussPair++;
  GroGaussPair ggp2;
  ggp2.indxGaussA = indexB;
  ggp2.indxGaussB = indexA;
  ggp2.gA = gaussA;
  ggp2.gMu1 = gaussMu1;
  ggp2.giSigma1 = gaussSigma1;
  ggp2.gMu2 = 0.0;
  ggp2.giSigma2 = 0.0;
  ggp2.gRepulsive = 0.0;
  numGaussPair++;
  pairlistGauss.add(ggp2);
  return 0;
}

void PairTable::getPairLJArrays2(int *indexA, int *indexB, Real *pairC6, Real *pairC12) {

  std::sort(pairlistLJ.begin(),pairlistLJ.end(),GroLJCompare);
  ResizeArray<GroLJPair>::iterator it;
  for(int i = 0; i < numLJPair; i++) {  
    indexA[i] = pairlistLJ[i].indxLJA;
    indexB[i] = pairlistLJ[i].indxLJB;
    pairC6[i] = pairlistLJ[i].c6pair;
    pairC12[i] = pairlistLJ[i].c12pair;
    }
}

void PairTable::getPairGaussArrays2(int *indexA, int *indexB, Real *gaussA, Real *gaussMu1,
				   Real *gaussSigma1, Real *gaussMu2, Real *gaussSigma2,
				   Real *gaussRepulsive) {
  std::sort(pairlistGauss.begin(),pairlistGauss.end(),GroGaussCompare);
  for(int i = 0; i < numGaussPair; i++) {
    indexA[i] = pairlistGauss[i].indxGaussA;
    indexB[i] = pairlistGauss[i].indxGaussB;
    gaussA[i] = pairlistGauss[i].gA;
    gaussMu1[i] = pairlistGauss[i].gMu1;
    gaussSigma1[i] = pairlistGauss[i].giSigma1;
    gaussMu2[i] = pairlistGauss[i].gMu2;
    gaussSigma2[i] = pairlistGauss[i].giSigma2;
    gaussRepulsive[i] = pairlistGauss[i].gRepulsive;
  }
}

bool PairTable::GroLJCompare (GroLJPair A, GroLJPair B) {
  if(A.indxLJA < B.indxLJA) {
    return true;
  } else if(A.indxLJA == B.indxLJA) {
    return (A.indxLJB < B.indxLJB);
  } 
  return false;
}

bool PairTable::GroGaussCompare (GroGaussPair A, GroGaussPair B) {
  if(A.indxGaussA < B.indxGaussA) {
    return true;
  } else if(A.indxGaussA == B.indxGaussA) {
    return (A.indxGaussB < B.indxGaussB);
  } 
  return false;
}
/* End of JLai Modifications */
