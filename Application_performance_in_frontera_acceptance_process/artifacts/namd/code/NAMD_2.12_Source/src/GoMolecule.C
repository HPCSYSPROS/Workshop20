/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "ResizeArray.h"
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
#ifndef CODE_REDUNDANT
#define CODE_REDUNDANT 0
#endif

// Ported by JLai from NAMD 2.7
/************************************************************************/
/*                                                                      */
/*      FUNCTION goInit                                                 */
/*                                                                      */
/*        This function is only called from the Molecule constructor.   */
/*   It only builds Go specific code into the Molecule object           */
/*                                                                      */
/************************************************************************/
void Molecule::goInit() {
  numGoAtoms=0;
  energyNative=0;
  energyNonnative=0;
  atomChainTypes=NULL;
  goSigmaIndices=NULL;
  goSigmas=NULL;
  goWithinCutoff=NULL;
  goCoordinates=NULL;
  goResids=NULL;
  goPDB=NULL;
}

/************************************************************************/
/*                                                                      */
/*      FUNCTION build_go_params                                        */
/*                                                                      */
/*   INPUTS:                                                            */
/*        fname - name of the parameter file to read                    */
/*                                                                      */
/*        This function reads in multiple Go parameter files            */
/*   from a StringList and exhaustively processes them.                 */
/*                                                                      */
/************************************************************************/
// JE - read in a Go parameter file
void Molecule::build_go_params(StringList *g) {
  iout << iINFO << "Building Go Parameters" << "\n" << endi;
#ifdef MEM_OPT_VERSION
  NAMD_die("Go forces are not supported in memory-optimized builds.");
#else
  build_lists_by_atom();
#endif
  int iterator = 0;
    do
    {
      iout << iINFO << "Reading Go File: " << iterator << "\n" << endi;
      read_go_file(g->data);
      g = g->next;
      iterator++;
    } while ( g != NULL && iterator < 100);    
}

// Ported by JLai from NAMD 2.7
/************************************************************************/
/*                                                                      */
/*      FUNCTION read_go_file                                           */
/*                                                                      */
/*   INPUTS:                                                            */
/*        fname - name of the parameter file to read                    */
/*                                                                      */
/*        This function reads in a Go parameter file and adds the       */ 
/*   parameters from this file to the current group of parameters.      */
/*   The basic logic of the routine is to first find out what type of   */
/*   parameter we have in the file. Then look at each line in turn      */
/*   and call the appropriate routine to add the parameters until we hit*/
/*   a new type of parameter or EOF.                                    */
/*                                                                      */
/************************************************************************/
// JE - read in a Go parameter file
void Molecule::read_go_file(char *fname)

{

  int i;                   // Counter
  int j;                   // Counter
  int  par_type=0;         //  What type of parameter are we currently
                           //  dealing with? (vide infra)
  // JLai -- uncommented
  int  skipline;           //  skip this line?
  int  skipall = 0;        //  skip rest of file;
  char buffer[512];           //  Buffer to store each line of the file
  char first_word[512];           //  First word of the current line
  int read_count = 0;      //  Count of input parameters on a given line
  int chain1 = 0;          //  First chain type for pair interaction
  int chain2 = 0;          //  Second chain type for pair interaction
  int int1;                //  First parameter int
  int int2;                //  Second parameter int
  Real r1;                 //  Parameter Real
  char in2[512];           //  Second parameter word
  GoValue *goValue1 = NULL;    //  current GoValue for loading parameters
  GoValue *goValue2 = NULL;    //  other current GoValue for loading parameters
  Bool sameGoChain = FALSE;    //  whether the current GoValue is within a chain or between chains
  int restrictionCount = 0;    //  number of Go restrictions set for a given chain pair
  FILE *pfile;                 //  File descriptor for the parameter file

  /*  Check to make sure that we haven't previously been told     */
  /*  that all the files were read                                */
  /*if (AllFilesRead)
    {
    NAMD_die("Tried to read another parameter file after being told that all files were read!");
    }*/
  
  /*  Initialize go_indices  */
  for (i=0; i<MAX_GO_CHAINS+1; i++) {
    go_indices[i] = -1;
  }

  /*  Initialize go_array   */
  for (i=0; i<MAX_GO_CHAINS*MAX_GO_CHAINS; i++) {
    go_array[i].epsilon = 1.25;
    go_array[i].exp_a = 12;
    go_array[i].exp_b = 6;
    go_array[i].exp_rep = 12;
    go_array[i].sigmaRep = 2.25;
    go_array[i].epsilonRep = 0.03;
    go_array[i].cutoff = 4.0;
    for (j=0; j<MAX_RESTRICTIONS; j++) {
      go_array[i].restrictions[j] = -1;
    }
  }

  /*  Try and open the file                                        */
  if ( (pfile = fopen(fname, "r")) == NULL)
    {
      char err_msg[256];
      
      sprintf(err_msg, "UNABLE TO OPEN GO PARAMETER FILE %s\n", fname);
      NAMD_die(err_msg);
    }
  
  /*  Keep reading in lines until we hit the EOF                        */
  while (NAMD_read_line(pfile, buffer) != -1)
    {
      /*  Get the first word of the line                        */
      NAMD_find_first_word(buffer, first_word);
      skipline=0;
      
      /*  First, screen out things that we ignore.                   */   
      /*  blank lines, lines that start with '!' or '*', lines that  */
      /*  start with "END".                                          */
      if (!NAMD_blank_string(buffer) &&
	  (strncmp(first_word, "!", 1) != 0) &&
	  (strncmp(first_word, "*", 1) != 0) &&
	  (strncasecmp(first_word, "END", 3) != 0))
	{
	  if ( skipall ) {
	    iout << iWARN << "SKIPPING PART OF GO PARAMETER FILE AFTER RETURN STATEMENT\n" << endi;
	    break;
	  }
	  /*  Now, determine the apropriate parameter type.   */
	  if (strncasecmp(first_word, "chaintypes", 10)==0)
	    {
	      read_count=sscanf(buffer, "%s %d %d\n", first_word, &int1, &int2);
	      if (read_count != 3) {
		char err_msg[512];
		sprintf(err_msg, "UNKNOWN PARAMETER IN GO PARAMETER FILE %s\nLINE=*%s*\nread_count=%d, int1=%d, int2=%d", fname, buffer, read_count, int1, int2);
		NAMD_die(err_msg);
	      }
	      chain1 = int1;
	      chain2 = int2;
	      if (chain1 < 1 || chain1 > MAX_GO_CHAINS ||
		  chain2 < 1 || chain2 > MAX_GO_CHAINS) {
		char err_msg[512];
		sprintf(err_msg, "GO PARAMETER FILE: CHAIN INDEX MUST BE [1-%d] %s\nLINE=*%s*\nread_count=%d, int1=%d, int2=%d", MAX_GO_CHAINS, fname, buffer, read_count, int1, int2);
		NAMD_die(err_msg);
	      }
	      if (go_indices[chain1] == -1) {
		go_indices[chain1] = NumGoChains;
		NumGoChains++;
	      }
	      if (go_indices[chain2] == -1) {
		go_indices[chain2] = NumGoChains;
		NumGoChains++;
	      }
	      if (chain1 == chain2) {
		sameGoChain = TRUE;
	      } else {
		sameGoChain = FALSE;
	      }
	      //goValue = &go_array[(chain1 * MAX_GO_CHAINS) + chain2];
	      goValue1 = &go_array[(chain1*MAX_GO_CHAINS) + chain2];
	      goValue2 = &go_array[(chain2*MAX_GO_CHAINS) + chain1];
#if CODE_REDUNDANT
	      goValue1 = &go_array[(go_indices[chain1]*MAX_GO_CHAINS) + go_indices[chain2]];
	      goValue2 = &go_array[(go_indices[chain2]*MAX_GO_CHAINS) + go_indices[chain1]];
#endif
	      restrictionCount = 0;    //  restrictionCount applies to each chain pair separately
	    }
	  else if (strncasecmp(first_word, "epsilonRep", 10)==0)
	    {
	      read_count=sscanf(buffer, "%s %f\n", first_word, &r1);
	      if (read_count != 2) {}
	      goValue1->epsilonRep = r1;
	      if (!sameGoChain) {
		goValue2->epsilonRep = r1;
	      }
	    }
	  else if (strncasecmp(first_word, "epsilon", 7)==0)
	    {
	      // Read in epsilon
	      read_count=sscanf(buffer, "%s %f\n", first_word, &r1);
	      if (read_count != 2) {}
	      goValue1->epsilon = r1;
	      if (!sameGoChain) {
		goValue2->epsilon = r1;
	      }
	    }
	  else if (strncasecmp(first_word, "exp_a", 5)==0)
	    {
	      read_count=sscanf(buffer, "%s %d\n", first_word, &int1);
	      if (read_count != 2) {}
	      goValue1->exp_a = int1;
	      if (!sameGoChain) {
		goValue2->exp_a = int1;
	      }
	    }
	  else if (strncasecmp(first_word, "exp_b", 5)==0)
	    {
	      read_count=sscanf(buffer, "%s %d\n", first_word, &int1);
	      if (read_count != 2) {}
	      goValue1->exp_b = int1;
	      if (!sameGoChain) {
		goValue2->exp_b = int1;
	      }
	    }
	  else if (strncasecmp(first_word, "exp_rep", 5)==0)
	    {
	      read_count=sscanf(buffer, "%s %d\n", first_word, &int1);
	      if (read_count != 2) {}
	      goValue1->exp_rep = int1;
	      if (!sameGoChain) {
		goValue2->exp_rep = int1;
	      }
	    }
          else if (strncasecmp(first_word, "exp_Rep", 5)==0)
	    {
	      read_count=sscanf(buffer, "%s %d\n", first_word, &int1);
	      if (read_count != 2) {}
	      goValue1->exp_rep = int1;
	      if (!sameGoChain) {
	      goValue2->exp_rep = int1;
	      }
	    }
	  else if (strncasecmp(first_word, "sigmaRep", 8)==0)
	    {
	      read_count=sscanf(buffer, "%s %f\n", first_word, &r1);
	      if (read_count != 2) {}
	      goValue1->sigmaRep = r1;
	      if (!sameGoChain) {
		goValue2->sigmaRep = r1;
	      }
	    }
	  else if (strncasecmp(first_word, "cutoff", 6)==0)
	    {
	      read_count=sscanf(buffer, "%s %f\n", first_word, &r1);
	      if (read_count != 2) {}
	      goValue1->cutoff = r1;
	      if (!sameGoChain) {
		goValue2->cutoff = r1;
	      }
	    }
	  else if (strncasecmp(first_word, "restriction", 10)==0)
	    {
	      read_count=sscanf(buffer, "%s %d\n", first_word, &int1);
	      if (read_count != 2) {}
	      if (int1 < 0) {
		DebugM(3, "ERROR: residue restriction value must be nonnegative.  int1=" << int1 << "\n");
	      }
	      if (!sameGoChain) {
		//goValue2->restrictions[restrictionCount] = int1;
		DebugM(3, "ERROR: residue restrictions should not be defined between two separate chains.  chain1=" << chain1 << ", chain2=" << chain2 << "\n");
	      }
	      else {
		goValue1->restrictions[restrictionCount] = int1;
	      }
	      restrictionCount++;
	    }
	  else if (strncasecmp(first_word, "return", 4)==0)
	    {
	      skipall=8;
	      skipline=1;
	    }        
	  else // if (par_type == 0)
	    {
	      /*  This is an unknown paramter.        */
	      /*  This is BAD                                */
	      char err_msg[512];
	      
	      sprintf(err_msg, "UNKNOWN PARAMETER IN GO PARAMETER FILE %s\nLINE=*%s*",fname, buffer);
	      NAMD_die(err_msg);
	    }
	}
      else
	{
	  skipline=1;
	}
    }
  
  /*  Close the file  */
  fclose(pfile);
  
  return;
}
/*              END OF FUNCTION read_go_file      */

#if CODE_REDUNDANT
/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_cutoff                                          */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go cutoff for two chain types.  The cutoff   */
/*   determines whether the Go force is attractive or repulsive.        */
/*                                                                      */
/************************************************************************/
// JE
Real Molecule::get_go_cutoff(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].cutoff;
  #if CODE_REDUNDANT
  return go_array[MAX_GO_CHAINS*go_indices[chain1] + go_indices[chain2]].cutoff;
  #endif
}
/*             END OF FUNCTION get_go_cutoff      */


/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_epsilonRep                                      */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go epsilonRep value for two chain types.     */
/*   epsilonRep is a factor in the repulsive Go force formula.          */
/*                                                                      */
/************************************************************************/
// JE
Real Molecule::get_go_epsilonRep(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].epsilonRep;
  #if CODE_REDUNDANT
  return go_array[MAX_GO_CHAINS*go_indices[chain1] + go_indices[chain2]].epsilonRep;
  #endif
}
/*             END OF FUNCTION get_go_epsilonRep      */


/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_sigmaRep                                        */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go sigmaRep value for two chain types.       */
/*   sigmaRep is a factor in the repulsive Go force formula.            */
/*                                                                      */
/************************************************************************/
// JE
Real Molecule::get_go_sigmaRep(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].sigmaRep;
}
/*             END OF FUNCTION get_go_sigmaRep        */


/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_epsilon                                         */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go epsilon value for two chain types.        */
/*   epsilon is a factor in the attractive Go force formula.            */
/*                                                                      */
/************************************************************************/
// JE
Real Molecule::get_go_epsilon(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].epsilon;
  #if CODE_REDUNDANT
  return go_array[MAX_GO_CHAINS*go_indices[chain1] + go_indices[chain2]].epsilon;
  #endif
}
/*             END OF FUNCTION get_go_epsilon         */


/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_exp_a                                           */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go exp_a value for two chain types.          */
/*   exp_a is an exponent in the repulsive term of the attractive Go    */
/*   force formula.                                                     */
/*                                                                      */
/************************************************************************/
// JE
int Molecule::get_go_exp_a(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].exp_a;
  #if CODE_REDUNDANT
  return go_array[MAX_GO_CHAINS*go_indices[chain1] + go_indices[chain2]].exp_a;
  #endif
}
/*             END OF FUNCTION get_go_exp_a        */


/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_exp_b                                           */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go exp_b value for two chain types.          */
/*   exp_b is an exponent in the attractive term of the attractive Go   */
/*   force formula.                                                     */
/*                                                                      */
/************************************************************************/
// JE
int Molecule::get_go_exp_b(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].exp_b;
  #if CODE_REDUNDANT
  return go_array[MAX_GO_CHAINS*go_indices[chain1] + go_indices[chain2]].exp_b;
  #endif
}
/*             END OF FUNCTION get_go_exp_b        */


/************************************************************************/
/*                                                                      */
/*      FUNCTION get_go_exp_rep                                         */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*                                                                      */
/*  This function gets the Go exp_rep value for two chain types.        */
/*   exp_b is an exponent in the attractive term of the attractive Go   */
/*   force formula.                                                     */
/*                                                                      */
/************************************************************************/
// JE
int Molecule::get_go_exp_rep(int chain1, int chain2)
{
  return go_array[MAX_GO_CHAINS*chain1 + chain2].exp_rep;
  #if CODE_REDUNDANT
  return go_array[MAX_GO_CHAINS*go_indices[chain1] + go_indices[chain2]].exp_rep;
  #endif
}
/*             END OF FUNCTION get_go_exp_rep      */
#endif

/************************************************************************/
/*                                                                      */
/*      FUNCTION go_restricted                                          */
/*                                                                      */
/*   INPUTS:                                                            */
/*     chain1 - first chain type                                        */
/*     chain2 - second chain type                                       */
/*     rDiff - residue ID difference to check                           */
/*                                                                      */
/*  This function checks whether residues with IDs rDiff apart are      */
/*   restricted from Go force calculation.                              */
/*                                                                      */
/************************************************************************/
// JE
Bool Molecule::go_restricted(int chain1, int chain2, int rDiff)
{
  int i;      //  Loop counter
  for(i=0; i<MAX_RESTRICTIONS; i++) {
    if (go_array[(MAX_GO_CHAINS*chain1) + chain2].restrictions[i]  == rDiff) {
      return TRUE;
    } else if (go_array[(MAX_GO_CHAINS*chain1) + chain2].restrictions[i] == -1) {
      return FALSE;
    }
  }
  return FALSE;
}
/*              END OF FUNCTION go_restricted      */

// Original by JE
/************************************************************************/
/*                  */
/*      FUNCTION print_go_params      */
/*                  */
/*  This is a debugging routine used to print out all the Go  */
/*  parameters                */
/*                  */
/************************************************************************/
void Molecule::print_go_params()
{
  int i;
  int j;
  int index;

  DebugM(3,NumGoChains << " Go PARAMETERS 3\n" \
	 << "*****************************************" << std::endl);

  for (i=0; i<NumGoChains; i++) {
    for (j=0; j<NumGoChains; j++) {
      index = (i * MAX_GO_CHAINS) + j;
      //  Real epsilon;    // Epsilon
      //  Real exp_a;      // First exponent for attractive L-J term
      //  Real exp_b;      // Second exponent for attractive L-J term
      //  Real sigmaRep;   // Sigma for repulsive term
      //  Real epsilonRep; // Epsilon for replusive term
      DebugM(3,"Go index=(" << i << "," << j << ") epsilon=" << go_array[index].epsilon \
	     << " exp_a=" << go_array[index].exp_a << " exp_b=" << go_array[index].exp_b \
	     << " exp_rep=" << go_array[index].exp_rep << " sigmaRep=" \
	     << go_array[index].sigmaRep << " epsilonRep=" << go_array[index].epsilonRep \
	     << " cutoff=" << go_array[index].cutoff << std::endl);
    }
  }

}
// End of port -- JLai

#ifndef MEM_OPT_VERSION
void Molecule::build_go_sigmas(StringList *goCoordFile, 
			       char *cwd)
{
  DebugM(3,"->build_go_sigmas" << std::endl);
  PDB *goPDB;      //  Pointer to PDB object to use
  int bcol = 4;      //  Column that data is in
  int32 chainType = 0;      //  b value from PDB file
  int i;         //  Loop counter
  int j;         //  Loop counter
  int resid1;    //  Residue ID for first atom
  int resid2;    //  Residue ID for second atom
  int residDiff;     //  Difference between resid1 and resid2
  Real sigma;    //  Sigma calculated for a Go atom pair
  Real atomAtomDist;     //  Distance between two atoms
  Real exp_a;            //  First exponent in L-J like formula
  Real exp_b;            //  Second exponent in L-J like formula
  char filename[129];    //  Filename
  
  //JLai
  BigReal nativeEnergy, *native;
  BigReal nonnativeEnergy, *nonnative;
  nativeEnergy = 0;
  nonnativeEnergy = 0;
  native = &nativeEnergy;
  nonnative = &nonnativeEnergy;

  long nativeContacts = 0;
  long nonnativeContacts = 0;

  //  Get the PDB object that contains the Go coordinates.  If
  //  the user gave another file name, use it.  Otherwise, just use
  //  the PDB file that has the initial coordinates.  
  if (goCoordFile == NULL)
    {
      //goPDB = initial_pdb;
      NAMD_die("Error: goCoordFile is NULL - build_go_sigmas");
    }
  else
  {
    if (goCoordFile->next != NULL)
      {
	NAMD_die("Multiple definitions of Go atoms PDB file in configuration file");
      }
    
    if ( (cwd == NULL) || (goCoordFile->data[0] == '/') )
      {
	strcpy(filename, goCoordFile->data);
      }
    else
      {
	strcpy(filename, cwd);
	strcat(filename, goCoordFile->data);
      }
    
    goPDB = new PDB(filename);
    if ( goPDB == NULL )
      {
	NAMD_die("Memory allocation failed in Molecule::build_go_sigmas");
      }
    
    if (goPDB->num_atoms() != numAtoms)
      {
	NAMD_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
      }
  }
  //  Allocate the array to hold the chain types
  atomChainTypes = new int32[numAtoms];
  //  Allocate the array to hold Go atom indices into the sigma array
  goSigmaIndices = new int32[numAtoms];
  
  if (atomChainTypes == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_go_sigmas");
  }
  
  numGoAtoms = 0;
  
  //  Loop through all the atoms and get the Go chain types
  for (i=0; i<numAtoms; i++) {
    //  Get the chainType from the occupancy field
    chainType = (int32)(goPDB->atom(i))->occupancy();
    //  Assign the chainType value
    if ( chainType != 0 ) {
      //DebugM(3,"build_go_sigmas - atom:" << i << ", chainType:" << chainType << std::endl);
      atomChainTypes[i] = chainType;
      goSigmaIndices[i] = numGoAtoms;
      numGoAtoms++;
    }
    else {
      atomChainTypes[i] = 0;
      goSigmaIndices[i] = -1;
    }
    //printf("CT: %d %d %d %d\n",i,numGoAtoms,atomChainTypes[i],goSigmaIndices[i]);
  }

  // Allocate the array to hold the sigma values for each Go atom pair
  goSigmas = new Real[numGoAtoms*numGoAtoms];
  goWithinCutoff = new bool[numGoAtoms*numGoAtoms];
  for (i=0; i<numGoAtoms; i++) {
    for (j=0; j<numGoAtoms; j++) {
      goSigmas[i*numGoAtoms + j] = -1.0;
      goWithinCutoff[i*numGoAtoms + j] = false;
    }
  }
  //  Loop through all atom pairs and calculate sigmas
  DebugM(3,"    numAtoms=" << numAtoms << std::endl);
  for (i=0; i<numAtoms; i++) {
    //DebugM(3,"    i=" << i << std::endl);
    resid1 = (goPDB->atom(i))->residueseq();
    //DebugM(3,"    resid1=" << resid1 << std::endl);
    //if ( goSigmaIndices[i] != -1) {
    //  goSigmas[goSigmaIndices[i]*numGoAtoms + goSigmaIndices[i]] = 0.0;
    //}
     for (j=i+1; j<numAtoms; j++) {
      //DebugM(3,"    j=" << j << std::endl);
      resid2 = (goPDB->atom(j))->residueseq();
      //printf("GSIi %d %d %d\n",i,numAtoms,goSigmaIndices[i]);
      //printf("SAN CHECK: %d\n",goSigmaIndices[37]);
      //printf("GSIj %d %d %d\n",j,numAtoms,goSigmaIndices[j]);
      //printf("ATOMS_1to4 %d\n",!atoms_1to4(i,j));
      //DebugM(3,"    resid2=" << resid2 << std::endl);
      //  if goSigmaIndices aren't defined, don't set anything in goSigmas
      if ( goSigmaIndices[i] != -1 && goSigmaIndices[j] != -1 && !atoms_1to4(i,j) ) {
	//printf("TAKING DIFFERENCE\n");
	residDiff = resid2 - resid1;
	//printf("RESIDDIFF %d\n",residDiff);
	if (residDiff < 0) residDiff = -residDiff;
	//printf("RESIDDIFF2 %d\n",residDiff);
	//  if this is a Go atom pair that is not restricted
	//    calculate sigma
	//  sigmas are initially set to -1.0 if the atom pair fails go_restricted
	//printf("CHECKING LOOPING\n");
	if ( atomChainTypes[i] && atomChainTypes[j] &&
	     !(this->go_restricted(atomChainTypes[i],atomChainTypes[j],residDiff)) ) {
	  atomAtomDist = sqrt(pow((goPDB->atom(i))->xcoor() - (goPDB->atom(j))->xcoor(), 2.0) +
			      pow((goPDB->atom(i))->ycoor() - (goPDB->atom(j))->ycoor(), 2.0) +
			      pow((goPDB->atom(i))->zcoor() - (goPDB->atom(j))->zcoor(), 2.0));
	  if ( atomAtomDist <= this->get_go_cutoff(atomChainTypes[i],atomChainTypes[j]) ) {
	    exp_a = this->get_go_exp_a(atomChainTypes[i],atomChainTypes[j]);
	    exp_b = this->get_go_exp_b(atomChainTypes[i],atomChainTypes[j]);
	    sigma = pow(static_cast<double>(exp_b/exp_a),(1.0/(exp_a-exp_b))) * atomAtomDist;
	    goSigmas[goSigmaIndices[i]*numGoAtoms + goSigmaIndices[j]] = sigma;
	    goSigmas[goSigmaIndices[j]*numGoAtoms + goSigmaIndices[i]] = sigma;
	    goWithinCutoff[goSigmaIndices[i]*numGoAtoms + goSigmaIndices[j]] = true;
	    goWithinCutoff[goSigmaIndices[j]*numGoAtoms + goSigmaIndices[i]] = true;
	    nativeContacts++;
	  } else {
	    goSigmas[goSigmaIndices[i]*numGoAtoms + goSigmaIndices[j]] = 0.0;
	    goSigmas[goSigmaIndices[j]*numGoAtoms + goSigmaIndices[i]] = 0.0;
	    nonnativeContacts++;
	  }
	} else {
	  goSigmas[goSigmaIndices[i]*numGoAtoms + goSigmaIndices[j]] = -1.0;
	  goSigmas[goSigmaIndices[j]*numGoAtoms + goSigmaIndices[i]] = -1.0;
	}
      } 
    }
  }

  iout << iINFO << "Number of UNIQUE    native contacts: " << nativeContacts << "\n" << endi;
  iout << iINFO << "Number of UNIQUE nonnative contacts: " << nonnativeContacts << "\n" << endi;
  
  //  If we had to create a PDB object, delete it now
  if (goCoordFile != NULL) {
    delete goPDB;
  }
  
  return;
}
/*      END OF FUNCTION build_go_sigmas    */

void Molecule::build_go_sigmas2(StringList *goCoordFile, 
			       char *cwd)
{
  DebugM(3,"->build_go_sigmas2" << std::endl);
  PDB *goPDB;      //  Pointer to PDB object to use
  int bcol = 4;      //  Column that data is in
  int32 chainType = 0;      //  b value from PDB file
  int32 residType = 0;      //  resid value from PDB file
  int i;         //  Loop counter
  int j;         //  Loop counter
  int resid1;    //  Residue ID for first atom
  int resid2;    //  Residue ID for second atom
  int residDiff;     //  Difference between resid1 and resid2
  Real sigma;    //  Sigma calculated for a Go atom pair
  Real atomAtomDist;     //  Distance between two atoms
  Real exp_a;            //  First exponent in L-J like formula
  Real exp_b;            //  Second exponent in L-J like formula
  char filename[129];    //  Filename
  
  //JLai
  int numLJPair = 0;
  long nativeContacts = 0;
  long nonnativeContacts = 0;

  //  Get the PDB object that contains the Go coordinates.  If
  //  the user gave another file name, use it.  Otherwise, just use
  //  the PDB file that has the initial coordinates.  
  if (goCoordFile == NULL)
    {
      //goPDB = initial_pdb;
      NAMD_die("Error: goCoordFile is NULL - build_go_sigmas2");
    }
  else
  {
    if (goCoordFile->next != NULL)
      {
	NAMD_die("Multiple definitions of Go atoms PDB file in configuration file");
      }
    
    if ( (cwd == NULL) || (goCoordFile->data[0] == '/') )
      {
	strcpy(filename, goCoordFile->data);
      }
    else
      {
	strcpy(filename, cwd);
	strcat(filename, goCoordFile->data);
      }
    
    goPDB = new PDB(filename);
    if ( goPDB == NULL )
      {
	NAMD_die("Memory allocation failed in Molecule::build_go_sigmas2");
      }
    
    if (goPDB->num_atoms() != numAtoms)
      {
	NAMD_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
      }
  }
  //  Allocate the array to hold the chain types
  atomChainTypes = new int32[numAtoms];
  //  Allocate the array to hold Go atom indices into the sigma array
  goSigmaIndices = new int32[numAtoms];
  //  Allocate the array to hold resid 
  goResidIndices = new int32[numAtoms];

  if (atomChainTypes == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_go_sigmas2");
  }
  
  numGoAtoms = 0;
  
  //  Loop through all the atoms and get the Go chain types
  for (i=0; i<numAtoms; i++) {
    //  Get the chainType from the occupancy field
    chainType = (int32)(goPDB->atom(i))->occupancy();
    //  Get the resid from the resid field
    residType = (int32)(goPDB->atom(i))->residueseq();
    //  Assign the chainType value
    if ( chainType != 0 ) {
      //DebugM(3,"build_go_sigmas2 - atom:" << i << ", chainType:" << chainType << std::endl);
      atomChainTypes[i] = chainType;
      goSigmaIndices[i] = numGoAtoms;
      goResidIndices[i] = residType;
      numGoAtoms++;
    }
    else {
      atomChainTypes[i] = 0;
      goSigmaIndices[i] = -1;
      goResidIndices[i] = -1;
    }
  }

  //Define:
  ResizeArray<GoPair> tmpGoPair;
  
  //  Loop through all atom pairs and calculate sigmas
  // Store sigmas into sorted array
  DebugM(3,"    numAtoms=" << numAtoms << std::endl);
  for (i=0; i<numAtoms; i++) {
    resid1 = (goPDB->atom(i))->residueseq();
     for (j=i+1; j<numAtoms; j++) {
      resid2 = (goPDB->atom(j))->residueseq();
      if ( goSigmaIndices[i] != -1 && goSigmaIndices[j] != -1 && !atoms_1to4(i,j) ) {
	residDiff = resid2 - resid1;
	if (residDiff < 0) residDiff = -residDiff;
	if ( atomChainTypes[i] && atomChainTypes[j] &&
	     !(this->go_restricted(atomChainTypes[i],atomChainTypes[j],residDiff)) ) {
	  atomAtomDist = sqrt(pow((goPDB->atom(i))->xcoor() - (goPDB->atom(j))->xcoor(), 2.0) +
			      pow((goPDB->atom(i))->ycoor() - (goPDB->atom(j))->ycoor(), 2.0) +
			      pow((goPDB->atom(i))->zcoor() - (goPDB->atom(j))->zcoor(), 2.0));
	  if ( atomAtomDist <= this->get_go_cutoff(atomChainTypes[i],atomChainTypes[j]) ) {
	    exp_a = this->get_go_exp_a(atomChainTypes[i],atomChainTypes[j]);
	    exp_b = this->get_go_exp_b(atomChainTypes[i],atomChainTypes[j]);
	    sigma = pow(static_cast<double>(exp_b/exp_a),(1.0/(exp_a-exp_b))) * atomAtomDist;
	    double tmpA = pow(sigma,exp_a);
	    double tmpB = pow(sigma,exp_b);
	    GoPair gp;
	    GoPair gp2;
	    gp.goIndxA = i;
	    gp.goIndxB = j;
	    gp.A = tmpA;
	    gp.B = tmpB;
	    tmpGoPair.add(gp);
	    gp2.goIndxA = j;
	    gp2.goIndxB = i;
	    gp2.A = tmpA;
	    gp2.B = tmpB;
	    tmpGoPair.add(gp2);
	    nativeContacts++;
	  } else {
	    nonnativeContacts++;
	  }
	} 
      } 
    }
  }

  iout << iINFO << "Number of UNIQUE    native contacts: " << nativeContacts << "\n" << endi;
  iout << iINFO << "Number of UNIQUE nonnative contacts: " << nonnativeContacts << "\n" << endi;

  // Copies the resizeArray into a block of continuous memory
  std::sort(tmpGoPair.begin(),tmpGoPair.end(),goPairCompare);
  goNumLJPair = 2*nativeContacts;
  goIndxLJA = new int[goNumLJPair];
  goIndxLJB = new int[goNumLJPair];
  goSigmaPairA = new double[goNumLJPair];
  goSigmaPairB = new double[goNumLJPair];
  for(i=0; i< goNumLJPair; i++) {
    goIndxLJA[i] = tmpGoPair[i].goIndxA;
    goIndxLJB[i] = tmpGoPair[i].goIndxB;
    goSigmaPairA[i] = tmpGoPair[i].A;
    goSigmaPairB[i] = tmpGoPair[i].B;
  }

  pointerToGoBeg = new int[numAtoms];
  pointerToGoEnd = new int[numAtoms];
  int oldIndex = -1;
  for(i=0; i<numAtoms; i++) {
    pointerToGoBeg[i] = -1;
    pointerToGoEnd[i] = -2;
  }
  for(i=0; i < goNumLJPair; i++) {
    if(pointerToGoBeg[goIndxLJA[i]] == -1) {
      pointerToGoBeg[goIndxLJA[i]] = i;
      oldIndex = goIndxLJA[i];
    }
    pointerToGoEnd[oldIndex] = i;
  }

  //  If we had to create a PDB object, delete it now
  if (goCoordFile != NULL) {
    delete goPDB;
  }
    
  return;
}
/*      END OF FUNCTION build_go_sigmas2    */

bool Molecule::goPairCompare(GoPair first, GoPair second) {
  if(first.goIndxA < second.goIndxA) {
    return true;
  } else if(first.goIndxA == second.goIndxA) {
    return (first.goIndxB == second.goIndxB);
  } 
  return false;
}

    /************************************************************************/
    /*                                                                      */
    /*      JE - FUNCTION build_go_arrays                                   */
    /*                                                                      */
    /*   INPUTS:                                                            */
    /*  goCoordFile - Value of Go coordinate file from config file          */
    /*  cwd - Current working directory                                     */
    /*                                                                      */
    /*  This function builds arrays that support sigma calculation for L-J  */
    /* style Go forces.  It takes the name of the PDB file.  It then builds */
    /* an array identifying atoms to which Go forces apply.                 */
    /*                                                                      */
    /************************************************************************/
// JE
void Molecule::build_go_arrays(StringList *goCoordFile, 
			      char *cwd)
{
  DebugM(3,"->build_go_arrays" << std::endl);
  //PDB *goPDB;      //  Pointer to PDB object to use
  int bcol = 4;      //  Column that data is in
  int32 chainType = 0;      //  b value from PDB file
  int i;         //  Loop counter
  int j;         //  Loop counter
  BigReal atomAtomDist;     //  Distance between two atoms -- JLai put back 
  int resid1;    //  Residue ID for first atom
  int resid2;    //  Residue ID for second atom
  int residDiff;     //  Difference between resid1 and resid2
  int goIndex;       //  Index into the goCoordinates array
  int goIndx;        //  Index into the goCoordinates array
  char filename[129];    //  Filename
  
  //JLai
  BigReal nativeEnergy, *native;
  BigReal nonnativeEnergy, *nonnative;
  nativeEnergy = 0;
  nonnativeEnergy = 0;
  native = &nativeEnergy;
  nonnative = &nonnativeEnergy;

  long nativeContacts = 0;
  long nonnativeContacts = 0;

  //  Get the PDB object that contains the Go coordinates.  If
  //  the user gave another file name, use it.  Otherwise, just use
  //  the PDB file that has the initial coordinates.  
  if (goCoordFile == NULL)
    {
      //goPDB = initial_pdb;
      NAMD_die("Error: goCoordFile is NULL - build_go_arrays");
    }
  else
  {
    if (goCoordFile->next != NULL)
      {
	NAMD_die("Multiple definitions of Go atoms PDB file in configuration file");
      }
    
    if ( (cwd == NULL) || (goCoordFile->data[0] == '/') )
      {
	strcpy(filename, goCoordFile->data);
      }
    else
      {
	strcpy(filename, cwd);
	strcat(filename, goCoordFile->data);
      }
    
    goPDB = new PDB(filename);
    if ( goPDB == NULL )
      {
	NAMD_die("goPDB memory allocation failed in Molecule::build_go_arrays");
      }
    
    if (goPDB->num_atoms() != numAtoms)
      {
	NAMD_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
      }
  }
  
  //  Allocate the array to hold Go atom indices into the sigma array
  goSigmaIndices = new int32[numAtoms];
  
  if (goSigmaIndices == NULL) {
    NAMD_die("goSigmaIndices memory allocation failed in Molecule::build_go_arrays");
  }
  
  numGoAtoms = 0;
  
  //  Loop through all the atoms and get the Go chain types
  for (i=0; i<numAtoms; i++) {
    chainType = (int32)(goPDB->atom(i))->occupancy();
    if ( chainType != 0 ) {
      DebugM(3,"build_go_arrays - atom:" << i << std::endl);
      goSigmaIndices[i] = numGoAtoms;
      numGoAtoms++;
    }
    else {
      goSigmaIndices[i] = -1;
    }
  }

  // Allocate the array to hold the sigma values for each Go atom pair
  /**
     goSigmas = new Real[numGoAtoms*numGoAtoms];
     goWithinCutoff = new bool[numGoAtoms*numGoAtoms];
     for (i=0; i<numGoAtoms; i++) {
     for (j=0; j<numGoAtoms; j++) {
     goSigmas[i*numGoAtoms + j] = 0.0;
     goWithinCutoff[i*numGoAtoms + j] = false;
     }
     }
  **/

  //  Allocate the array to hold the chain types
  atomChainTypes = new int32[numGoAtoms];

  if (atomChainTypes == NULL) {
    NAMD_die("atomChainTypes memory allocation failed in Molecule::build_go_arrays");
  }

  // Allocate the array to hold (x,y,z) coordinates for all of the Go atoms
  goCoordinates = new Real[numGoAtoms*3];

  if (goCoordinates == NULL) {
    NAMD_die("goCoordinates memory allocation failed in Molecule::build_go_arrays");
  }

  goResids = new int[numGoAtoms];

  // Allocate the array to hold PDB residu IDs for all of the Go atoms
  if (goResids == NULL) {
    NAMD_die("goResids memory allocation failed in Molecule::build_go_arrays");
  }
  
  for (i=0; i<numAtoms; i++) {
    goIndex = goSigmaIndices[i];
    if (goIndex != -1) {
      //  Assign the chainType value!
      //  Get the chainType from the occupancy field
      atomChainTypes[goIndex] = (int32)(goPDB->atom(i))->occupancy();
      goCoordinates[goIndex*3] = goPDB->atom(i)->xcoor();
      goCoordinates[goIndex*3 + 1] = goPDB->atom(i)->ycoor();
      goCoordinates[goIndex*3 + 2] = goPDB->atom(i)->zcoor();
      goResids[goIndex] = goPDB->atom(i)->residueseq();
    }
  }
      // JLai
  energyNative = 0;
  energyNonnative = 0;
  //printf("INIT ENERGY: (N) %f (NN) %f\n", energyNative, energyNonnative);
  for (i=0; i<numAtoms-1; i++) {
    goIndex = goSigmaIndices[i];
    if (goIndex != -1) {
      for (j=i+1; j<numAtoms;j++) {
	goIndx = goSigmaIndices[j];
	if (goIndx != -1) {
          resid1 = (goPDB->atom(i))->residueseq();
          resid2 = (goPDB->atom(j))->residueseq();
          residDiff = resid2 - resid1;
          if (residDiff < 0) residDiff = -residDiff;
          if (atomChainTypes[goIndex] && atomChainTypes[goIndx] &&
              !(this->go_restricted(atomChainTypes[goIndex],atomChainTypes[goIndx],residDiff)) &&
	      !atoms_1to4(i,j)) {
	    atomAtomDist = sqrt(pow((goPDB->atom(i))->xcoor() - (goPDB->atom(j))->xcoor(), 2.0) +
				pow((goPDB->atom(i))->ycoor() - (goPDB->atom(j))->ycoor(), 2.0) +
				pow((goPDB->atom(i))->zcoor() - (goPDB->atom(j))->zcoor(), 2.0));
            if (atomAtomDist <= this->get_go_cutoff(atomChainTypes[goIndex],atomChainTypes[goIndx]) ) {
              nativeContacts++;
            } else {
              nonnativeContacts++;
            }    
	  }
	}
      }
    }
  }
  iout << iINFO << "Number of UNIQUE    native contacts: " << nativeContacts     << "\n" << endi;
  iout << iINFO << "Number of UNIQUE nonnative contacts: " << nonnativeContacts  << "\n" << endi;
  
  //  If we had to create a PDB object, delete it now
  if (goCoordFile != NULL) {
    delete goPDB;
  }
  
  return;
}
/*      END OF FUNCTION build_go_arrays    */
#endif // #ifndef MEM_OPT_VERSION

/************************************************************************/
/*                                                                      */
/*      FUNCTION print_go_sigmas                                        */
/*                                                                      */
/*  print_go_sigmas prints out goSigmas, the array of sigma parameters  */
/*   used in the L-J type Go force calculations                         */
/*                                                                      */
/************************************************************************/
// JE
void Molecule::print_go_sigmas()
{
  int i;  //  Counter
  int j;  //  Counter
  Real sigma;

  DebugM(3,"GO SIGMA ARRAY\n" \
         << "***************************" << std::endl);
  DebugM(3, "numGoAtoms: " << numGoAtoms << std::endl);

  if (goSigmaIndices == NULL) {
    DebugM(3, "GO SIGMAS HAVE NOT BEEN BUILT" << std::endl);
    return;
  }

  for (i=0; i<numAtoms; i++) {
    for (j=0; j<numAtoms; j++) {
      if ( goSigmaIndices[i] != -1 && goSigmaIndices[j] != -1 ) {
        //DebugM(3, "i: " << i << ", j: " << j << std::endl);
        sigma = goSigmas[goSigmaIndices[i]*numGoAtoms + goSigmaIndices[j]];
        if (sigma > 0.0) {
          DebugM(3, "(" << i << "," << j << ") - +" << sigma << " ");
        }
        else {
          DebugM(3, "(" << i << "," << j << ") - " << sigma << " ");
        }
      } else {
        //DebugM(3, "0 ");
      }
    }
    if ( goSigmaIndices[i] != -1 ) {
      DebugM(3, "-----------" << std::endl);
    }
  }
  return;
}
/*      END OF FUNCTION print_go_sigmas     */

// JLai
BigReal Molecule::get_gro_force2(BigReal x,
				 BigReal y,
				 BigReal z,
				 int atom1,
				 int atom2,
				 BigReal* pairLJEnergy,
				 BigReal* pairGaussEnergy) const
{
  //Initialize return energies to zero
  *pairLJEnergy = 0.0;
  *pairGaussEnergy = 0.0;

  // Linear search for LJ data
  int LJIndex = -1;
  int LJbegin = pointerToLJBeg[atom1];
  int LJend = pointerToLJEnd[atom1];
  for(int i = LJbegin; i <= LJend; i++) {
    if(indxLJB[i] == atom2) {
      LJIndex = i;
      break;
    }
  }

  // Linear search for Gaussian data
  int GaussIndex = -1;
  int Gaussbegin = pointerToGaussBeg[atom1];
  int Gaussend = pointerToGaussEnd[atom1];
  for(int i = Gaussbegin; i <= Gaussend; i++) {
    if(indxGaussB[i] == atom2) {
      GaussIndex = i;
      break;
    }
  }

  if( LJIndex == -1 && GaussIndex == -1) {
    return 0;
  } else {
    // Code to calculate distances because the pair was found in one of the lists
    BigReal r2 = x*x + y*y + z*z;
    BigReal r = sqrt(r2);
    BigReal ri = 1/r;
    BigReal ri6 = (ri*ri*ri*ri*ri*ri);
    BigReal ri12 = ri6*ri6;
    BigReal ri13 = ri12*ri;
    BigReal LJ = 0;
    BigReal Gauss = 0;
    // Code to calculate LJ 
    if (LJIndex != -1) {
      BigReal ri7 = ri6*ri;
      LJ = (12*(pairC12[LJIndex]*ri13) - 6*(pairC6[LJIndex]*ri7));
      *pairLJEnergy = (pairC12[LJIndex]*ri12 - pairC6[LJIndex]*ri6);
      //std::cout << pairC12[LJIndex] << " " << pairC6[LJIndex] << " " << ri13 << " " << ri7 << " " << LJ << " " << r << "\n";
    }
    // Code to calculate Gaussian
    if (GaussIndex != -1) {
      BigReal gr = 12*gRepulsive[GaussIndex]*ri13;
      BigReal r1prime = r - gMu1[GaussIndex];
      BigReal tmp1 = r1prime * r1prime;
      BigReal r2prime = r - gMu2[GaussIndex];
      BigReal tmp2 = r2prime * r2prime;
      BigReal tmp_gauss1 = 0;
      BigReal one_gauss1 = 1;
      BigReal tmp_gauss2 = 0;
      BigReal one_gauss2 = 1;
      if (giSigma1[GaussIndex] != 0) {
	tmp_gauss1 = exp(-tmp1*giSigma1[GaussIndex]);
	one_gauss1 = 1 - tmp_gauss1;
      }
      if (giSigma2[GaussIndex] != 0) {
	tmp_gauss2 = exp(-tmp2*giSigma2[GaussIndex]);
	one_gauss2 = 1 - tmp_gauss2;
      } 
      BigReal A = gA[GaussIndex];
      Gauss = gr*one_gauss1*one_gauss2 - 2*A*tmp_gauss1*one_gauss2*r1prime*giSigma1[GaussIndex] \
	- 2*tmp_gauss1*one_gauss2*r1prime*giSigma1[GaussIndex]*gRepulsive[GaussIndex]*ri12 - 2*A*tmp_gauss2*one_gauss1*r2prime*giSigma2[GaussIndex] \
	- 2*tmp_gauss2*one_gauss1*r2prime*giSigma2[GaussIndex]*gRepulsive[GaussIndex]*ri12;
      *pairGaussEnergy = A*(-1+(one_gauss1)*(one_gauss2)*(1+gRepulsive[GaussIndex]*ri12/A));
    }
    //std::cout << "Net force: " << (LJ + Gauss) << " with ri " << (LJ + Gauss)*ri << "\n";
    return (LJ + Gauss)*ri;
  }
  return 0;
}
// End of get_gro_force2
// JLai

// JE
BigReal Molecule::get_go_force(BigReal r, 
			       int atom1,
			       int atom2,
			       BigReal* goNative,
			       BigReal* goNonnative) const
{
  BigReal goForce = 0.0;
  Real pow1;
  Real pow2;
  //  determine which Go chain pair we are working with
  //DebugM(3,"get_go_force - (" << atom1 << "," << atom2 << ")" << std::endl);
  int32 chain1 = atomChainTypes[atom1];
  int32 chain2 = atomChainTypes[atom2];

  //DebugM(3,"  chain1:" << chain1 << ", chain2:" << chain2 << std::endl);
  if (chain1 == 0 || chain2 == 0)  return 0.0;

  //  retrieve Go cutoff for this chain pair
  //TMP// JLai -- I'm going to replace this with a constant accessor.  This is just a temporary thing
  Real goCutoff = const_cast<Molecule*>(this)->get_go_cutoff(chain1,chain2);
  //DebugM(3,"  goCutoff:" << goCutoff << std::endl);
  if (goCutoff == 0)  return 0.0;
  //  if repulsive then calculate repulsive
  //  sigmas are initially set to -1.0 if the atom pair fails go_restricted
  if (goSigmas[goSigmaIndices[atom1]*numGoAtoms + goSigmaIndices[atom2]] != -1.0) {
    if (!goWithinCutoff[goSigmaIndices[atom1]*numGoAtoms + goSigmaIndices[atom2]]) {
      Real epsilonRep = const_cast<Molecule*>(this)->get_go_epsilonRep(chain1,chain2);
      Real sigmaRep = const_cast<Molecule*>(this)->get_go_sigmaRep(chain1,chain2);
      int exp_rep = const_cast<Molecule*>(this)->get_go_exp_rep(chain1,chain2);
      pow1 = pow(sigmaRep/r,exp_rep);
      goForce = 4*((exp_rep/(r*r)) * epsilonRep * pow1);
      *goNative = 0.0;
      *goNonnative = (4 * epsilonRep * pow1 );
      //DebugM(3,"get_go_force - (" << atom1 << "," << atom2 << ") chain1:" << chain1 << ", chain2:" << chain2 << ", epsilonRep:" << epsilonRep << ", sigmaRep:" << sigmaRep << ", r:" << r << ", goForce:" << goForce << std::endl);
    }
    //  if attractive then calculate attractive
    else {
      int goSigmaIndex1 = goSigmaIndices[atom1];
      int goSigmaIndex2 = goSigmaIndices[atom2];
      if (goSigmaIndex1 != -1 && goSigmaIndex2 != -1) {
	Real epsilon = const_cast<Molecule*>(this)->get_go_epsilon(chain1,chain2);
	int exp_a = const_cast<Molecule*>(this)->get_go_exp_a(chain1,chain2);
	int exp_b = const_cast<Molecule*>(this)->get_go_exp_b(chain1,chain2);
	Real sigma_ij = goSigmas[goSigmaIndices[atom1]*numGoAtoms + goSigmaIndices[atom2]];
	// Positive gradient of potential, not negative gradient of potential
        pow1 = pow(sigma_ij/r,exp_a);
        pow2 = pow(sigma_ij/r,exp_b);
	goForce = ((4/(r*r)) * epsilon * (exp_a * pow1 - exp_b * pow2));
	//DebugM(3,"get_go_force - (" << atom1 << "," << atom2 << ") chain1:" << chain1 << ", chain2:" << chain2 << ", sigma_ij:" << sigma_ij << ", r:" << r << ", goForce:" << goForce << std::endl);
	*goNative = (4 * epsilon * ( pow1 -  pow2 ) );
        *goNonnative = 0.0;
      }
    }
  }
  //DebugM(3,"goForce:" << goForce << std::endl);
  return goForce;
}
/*      END OF FUNCTION get_go_force_old       */


    /************************************************************************/
    /*                                                                      */
    /*      JE - FUNCTION get_go_force_new                                  */
    /*                                                                      */
    /*   INPUTS:                                                            */
    /*  r - distance between the two atoms                                  */
    /*  atom1 - the ID of the first atom                                    */
    /*  atom2 - the ID of the second atom                                   */
    /*                                                                      */
    /*  This function calculates the Go force between two atoms.  If the    */
    /*   atoms do not have Go parameters or sigmas, 0 is returned.          */
    /*                                                                      */
    /************************************************************************/
// JE
BigReal Molecule::get_go_force_new(BigReal r,
				   int atom1,
				   int atom2,
				   BigReal* goNative,
				   BigReal* goNonnative) const
{
  int resid1;
  int resid2;
  int residDiff;
  Real xcoorDiff;
  Real ycoorDiff;
  Real zcoorDiff;
  Real atomAtomDist;
  Real exp_a;
  Real exp_b;
  Real sigma_ij;
  Real epsilon;
  Real epsilonRep;
  Real sigmaRep;
  Real expRep;
  Real pow1;
  Real pow2;
  
  BigReal goForce = 0.0;
  *goNative = 0;
  *goNonnative = 0;

  //  determine which Go chain pair we are working with
  DebugM(3,"get_go_force - (" << atom1 << "," << atom2 << ")" << std::endl);
  int goIndex1 = goSigmaIndices[atom1];
  int goIndex2 = goSigmaIndices[atom2];

  int32 chain1 = atomChainTypes[goIndex1];
  int32 chain2 = atomChainTypes[goIndex2];

  DebugM(3,"  chain1:" << chain1 << ", chain2:" << chain2 << std::endl);
  if (chain1 == 0 || chain2 == 0)  return 0.0;

  //  retrieve Go cutoff for this chain pair
  Real goCutoff = const_cast<Molecule*>(this)->get_go_cutoff(chain1,chain2);
  DebugM(3,"  goCutoff:" << goCutoff << std::endl);
  if (goCutoff == 0)  return 0.0;

  //  sigmas are initially set to -1.0 if the atom pair fails go_restricted
  //  no goSigmas array anymore
  //Real sigma_ij = goSigmas[goSigmaIndices[atom1]*numGoAtoms + goSigmaIndices[atom2]];

  // XXX - used to be a condition for the following if
  //if the atoms are within 4 of each other
  //if ( !atoms_1to4(atom1,atom2) ) {

  //  if goSigmaIndices aren't defined, don't calculate forces
  if ( goIndex1 != -1 && goIndex2 != -1 ) {
    resid1 = goResids[goIndex1];
    resid2 = goResids[goIndex2];
    residDiff = resid2 - resid1;
    if (residDiff < 0) residDiff = -residDiff;
    //  if this is a Go atom pair that is not restricted
    if ( !(const_cast<Molecule*>(this)->go_restricted(chain1,chain2,residDiff)) ) {
      xcoorDiff = goCoordinates[goIndex1*3] - goCoordinates[goIndex2*3];
      ycoorDiff = goCoordinates[goIndex1*3 + 1] - goCoordinates[goIndex2*3 + 1];
      zcoorDiff = goCoordinates[goIndex1*3 + 2] - goCoordinates[goIndex2*3 + 2];
      atomAtomDist = sqrt(xcoorDiff*xcoorDiff + ycoorDiff*ycoorDiff + zcoorDiff*zcoorDiff);
      
      //  if attractive then calculate attractive
      if ( atomAtomDist <= const_cast<Molecule*>(this)->get_go_cutoff(chain1,chain2) ) {
	exp_a = const_cast<Molecule*>(this)->get_go_exp_a(chain1,chain2);
	exp_b = const_cast<Molecule*>(this)->get_go_exp_b(chain1,chain2);
	sigma_ij = pow(static_cast<double>(exp_b/exp_a),(1.0/(exp_a-exp_b))) * atomAtomDist;
	
	// [JLai] print out atoms involved in native contacts
	// printf("ATOM1: %d C1: %d ATOM2: %d C2: %d\n", atom1,chain1,atom2,chain2);

	epsilon = const_cast<Molecule*>(this)->get_go_epsilon(chain1,chain2);
	pow1 = pow(sigma_ij/r,static_cast<double>(exp_a));
	pow2 = pow(sigma_ij/r,static_cast<double>(exp_b));
	//goForce = ((4/r) * epsilon * (exp_a * pow(sigma_ij/r,exp_a) - exp_b * pow(sigma_ij/r,exp_b)));
	goForce = ((4/(r*r)) * epsilon * (exp_a * pow1 - exp_b * pow2));
	DebugM(3,"get_go_force - (" << atom1 << "," << atom2 << ") chain1:" << chain1 << ", chain2:" << chain2 << ", exp_a:" << exp_a << ", exp_b:" << exp_b << ", sigma_ij:" << sigma_ij << ", r:" << r << ", goForce:" << goForce << std::endl);
	//goEnergy = (4 * epsilon * ( pow(sigma_ij/r,exp_a) -  pow(sigma_ij/r,exp_b) ) ); // JLai I changed some of the expressions
	*goNative = (4 * epsilon * ( pow1 -  pow2 ) ); 
	*goNonnative = 0;
      }
      
      //  if repulsive then calculate repulsive
      else {
	epsilonRep = const_cast<Molecule*>(this)->get_go_epsilonRep(chain1,chain2);
	sigmaRep = const_cast<Molecule*>(this)->get_go_sigmaRep(chain1,chain2);
	expRep = const_cast<Molecule*>(this)->get_go_exp_rep(chain1,chain2);
	pow1 = pow(sigmaRep/r,(BigReal)expRep);
	//goForce = ((12.0/r) * epsilonRep * pow(sigmaRep/r,12.0));
	goForce = (4*(expRep/(r*r)) * epsilonRep * pow1);
	DebugM(3,"get_go_force - (" << atom1 << "," << atom2 << ") chain1:" << chain1 << ", chain2:" << chain2 << ", epsilonRep:" << epsilonRep << ", sigmaRep:" << sigmaRep << ", r:" << r << ", goForce:" << goForce << std::endl);
	//goEnergy = (4 * epsilonRep * pow(sigmaRep/r,12.0)); // JLai I changed some of the expressions
	*goNonnative = (4 * epsilonRep * pow1); 
	*goNative = 0;
      }
    }
  }
  
  //DebugM(3,"goForce:" << goForce << std::endl);
  return goForce;
}
/*      END OF FUNCTION get_go_force_new   */


    /************************************************************************/
    /*                                                                      */
    /*   JLai - FUNCTION get_go_force2                                      */
    /*                                                                      */
    /*   INPUTS:                                                            */
    /*  x - the x distance between p_i and p_j                              */
    /*  y - the y distance between p_i and p_j                              */
    /*  z - the z distance between p_i and p_j                              */
    /*  atom1 - the ID of the second atom                                   */
    /*  atom2 - the ID of the second atom                                   */
    /*                                                                      */
    /*  This function returns the force between two input atoms give their  */
/*  distance and their atom indices.                                        */
    /*                                                                      */
    /************************************************************************/
// JLai
BigReal Molecule::get_go_force2(BigReal x,
				BigReal y,
				BigReal z,
				int atom1,
				int atom2,
				BigReal *goNative,
				BigReal *goNonnative) const
{
 
  // Check to see if restricted.  If so, escape function early
  int32 chain1 = atomChainTypes[atom1];
  int32 chain2 = atomChainTypes[atom2];

  if(chain1 == 0 || chain2 == 0) return 0.0;
  Molecule *mol = const_cast<Molecule*>(this);
  Real goCutoff = mol->get_go_cutoff(chain1,chain2);
  if(goCutoff == 0) return 0.0;

  int resid1 = goResidIndices[atom1];
  int resid2 = goResidIndices[atom2];
  int residDiff = abs(resid1 - resid2);
  if((mol->go_restricted(chain1,chain2,residDiff))) {
    return 0.0;
  }

  int LJIndex = -1;
  int LJbegin = pointerToGoBeg[atom1];
  int LJend = pointerToGoEnd[atom1];
  for(int i = LJbegin; i <= LJend; i++) {
    if(goIndxLJB[i] == atom2) {
      LJIndex = i;
    }
  }
  
  BigReal r2 = x*x + y*y + z*z;
  BigReal r = sqrt(r2);

  if (LJIndex == -1) {
    int exp_rep = const_cast<Molecule*>(this)->get_go_exp_rep(chain1,chain2);
    BigReal epsilonRep = const_cast<Molecule*>(this)->get_go_epsilonRep(chain1, chain2);
    BigReal sigmaRep = const_cast<Molecule*>(this)->get_go_sigmaRep(chain1, chain2);
    double sigmaRepPow = pow(sigmaRep,exp_rep);
    BigReal LJ = (4*epsilonRep*exp_rep*sigmaRepPow*pow(r,-(exp_rep+1)));
    *goNative = 0;
    *goNonnative = (4*epsilonRep*sigmaRepPow*pow(r,-exp_rep));
    //*goNonnative = (4*epsilonRep * pow(sigmaRep/r,exp_rep));
    return (LJ/r);
  } else {
    // Code to calculate distances because the pair was found in one of the lists
    int exp_a = const_cast<Molecule*>(this)->get_go_exp_a(chain1,chain2);
    int exp_b = const_cast<Molecule*>(this)->get_go_exp_b(chain1,chain2);
    // We want the force, so we have to take the n+1 derivative
    BigReal powA = pow(r,-(exp_a + 1));
    BigReal powB = pow(r,-(exp_b + 1));
    BigReal powaa = pow(r,-exp_a);
    BigReal powbb = pow(r,-exp_b);
    BigReal epsilon = const_cast<Molecule*>(this)->get_go_epsilon(chain1,chain2);
    BigReal LJ = 4 * epsilon * (exp_a*goSigmaPairA[LJIndex]*powA - exp_b*goSigmaPairB[LJIndex]*powB);
    *goNative =  4 * epsilon * (goSigmaPairA[LJIndex]*powaa - goSigmaPairB[LJIndex]*powbb);
    *goNonnative = 0;
    return (LJ/r);
  }
  return 0;
}
// JLai
/*      END OF FUNCTION get_go_force2   */

#ifndef MEM_OPT_VERSION
    /************************************************************************/
    /*                                                                      */
    /*      JE - FUNCTION atoms_1to4                                        */
    /*                                                                      */
    /*   INPUTS:                                                            */
    /*  atom1 - the ID of the first atom                                  */
    /*  atom2 - the ID of the second atom                                 */
    /*                                                                      */
    /*  This function tells whether or not the two input atoms are within   */
    /*   1 to 4 bonds away from each other: bonds, angles, or dihedrals.    */
    /*                                                                      */
    /************************************************************************/
// JE
Bool Molecule::atoms_1to4(unsigned int atom1,
			  unsigned int atom2)
{
  int bondNum;   //  Bonds in bonded list
  int angleNum;  //  Angles in angle list
  int dihedralNum;   //  Dihedrals in dihedral list
  int *bonds;
  int *angles;
  int *dihedrals;
  Bond *bond;     //  Temporary bond structure
  Angle *angle;   //  Temporary angle structure
  Dihedral *dihedral; //  Temporary dihedral structure

  DebugM(2,"atoms_1to4(" << atom1 << "," << atom2 << ")" << std::endl);

  bonds = get_bonds_for_atom(atom1);
  bondNum = *bonds;
  while(bondNum != -1) {
    bond = get_bond(bondNum);
    DebugM(2,"bond  atom1:" << bond->atom1 << ", atom2:" << bond->atom2 << std::endl);
    if (atom2 == bond->atom1 || atom2 == bond->atom2) {
      return TRUE;
    }
    bondNum = *(++bonds);
  }

  bonds = get_bonds_for_atom(atom2);
  bondNum = *bonds;
  while(bondNum != -1) {
    bond = get_bond(bondNum);
    DebugM(2,"bond  atom1:" << bond->atom1 << ", atom2:" << bond->atom2 << std::endl);
    if (atom1 == bond->atom1 || atom1 == bond->atom2) {
      return TRUE;
    }
    bondNum = *(++bonds);
  }

  angles = get_angles_for_atom(atom1);
  angleNum = *angles;
  while(angleNum != -1) {
    angle = get_angle(angleNum);
    DebugM(2,"angle  atom1:" << angle->atom1 << ", atom2:" << angle->atom2 << ", atom3:" << angle->atom3 << std::endl);
    if (atom2 == angle->atom1 || atom2 == angle->atom2 || atom2 == angle->atom3) {
      return TRUE;
    }
    angleNum = *(++angles);
  }

  angles = get_angles_for_atom(atom2);
  angleNum = *angles;
  while(angleNum != -1) {
    angle = get_angle(angleNum);
    DebugM(2,"angle  atom1:" << angle->atom1 << ", atom2:" << angle->atom2 << ", atom3:" << angle->atom3 << std::endl);
    if (atom1 == angle->atom1 || atom1 == angle->atom2 || atom1 == angle->atom3) {
      return TRUE;
    }
    angleNum = *(++angles);
  }

  dihedrals = get_dihedrals_for_atom(atom1);
  dihedralNum = *dihedrals;
  while(dihedralNum != -1) {
    dihedral = get_dihedral(dihedralNum);
    DebugM(2,"dihedral  atom1:" << dihedral->atom1 << ", atom2:" << dihedral->atom2 << ", atom3:" << dihedral->atom3 << ", atom4:" << dihedral->atom4 << std::endl);
    if (atom2 == dihedral->atom1 || atom2 == dihedral->atom2 \
	|| atom2 == dihedral->atom3 || atom2 == dihedral->atom4) {
      return TRUE;
    }
    dihedralNum = *(++dihedrals);
  }

  dihedrals = get_dihedrals_for_atom(atom2);
  dihedralNum = *dihedrals;
  while(dihedralNum != -1) {
    dihedral = get_dihedral(dihedralNum);
    DebugM(2,"dihedral  atom1:" << dihedral->atom1 << ", atom2:" << dihedral->atom2 << ", atom3:" << dihedral->atom3 << ", atom4:" << dihedral->atom4 << std::endl);
    if (atom1 == dihedral->atom1 || atom1 == dihedral->atom2 \
	|| atom1 == dihedral->atom3 || atom1 == dihedral->atom4) {
      return TRUE;
    }
    dihedralNum = *(++dihedrals);
  }
  
  return FALSE;
}
/*      END OF FUNCTION atoms_1to4       */
#endif // #ifndef MEM_OPT_VERSION

//JLai
/************************************************************************/
/*                  */
/*      FUNCTION send_GoMolecule        */
/*                  */
/*  send_Molecule is used by the Master node to distribute the      */
/*   Go information to all the client nodes.  It is NEVER called*/
/*   by the client nodes.              */
/*                  */
/************************************************************************/
void Molecule::send_GoMolecule(MOStream *msg) {
  Real *a1, *a2, *a3, *a4;
  int *i1, *i2, *i3, *i4;
  int maxGoChainsSqr = MAX_GO_CHAINS*MAX_GO_CHAINS;  // JE JLai Go code
  msg->put(NumGoChains);
  
  if (NumGoChains) {
    //      int go_indices[MAX_GO_CHAINS+1];        //  Indices from chainIDs to go_array      
    //      GoValue go_array[MAX_GO_CHAINS*MAX_GO_CHAINS];   //  Array of Go params
    msg->put(MAX_GO_CHAINS+1,go_indices);

    a1 = new Real[maxGoChainsSqr];
    a2 = new Real[maxGoChainsSqr];
    a3 = new Real[maxGoChainsSqr];
    a4 = new Real[maxGoChainsSqr];
    i1 = new int[maxGoChainsSqr];
    i2 = new int[maxGoChainsSqr];
    i3 = new int[maxGoChainsSqr];
    i4 = new int[maxGoChainsSqr*MAX_RESTRICTIONS];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (a4 == NULL) || 
         (i1 == NULL) || (i2 == NULL) || (i3 == NULL) || (i4 == NULL) )
    {
      NAMD_die("memory allocation failed in Molecules::send_Molecules");
    }

    for (int i=0; i<maxGoChainsSqr; i++) {
      a1[i] = go_array[i].epsilon;
      a2[i] = go_array[i].sigmaRep;
      a3[i] = go_array[i].epsilonRep;
      a4[i] = go_array[i].cutoff;
      i1[i] = go_array[i].exp_a;
      i2[i] = go_array[i].exp_b;
      i3[i] = go_array[i].exp_rep;
      for (int j=0; j<MAX_RESTRICTIONS; j++) {
	i4[i*MAX_RESTRICTIONS + j] = go_array[i].restrictions[j];
      }
    }

    msg->put(maxGoChainsSqr, a1);
    msg->put(maxGoChainsSqr, a2);
    msg->put(maxGoChainsSqr, a3);
    msg->put(maxGoChainsSqr, a4);
    msg->put(maxGoChainsSqr, i1);
    msg->put(maxGoChainsSqr, i2);
    msg->put(maxGoChainsSqr, i3);
    msg->put(maxGoChainsSqr*MAX_RESTRICTIONS, i4);

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
    delete [] i1;
    delete [] i2;
    delete [] i3;
    delete [] i4;
  }

  //Ported JLai
  if (simParams->goForcesOn) {
    switch(simParams->goMethod) {
    case 1:
      msg->put(numGoAtoms);
      msg->put(numAtoms, goSigmaIndices);
      msg->put(numGoAtoms, atomChainTypes);
      msg->put(numGoAtoms*numGoAtoms, goSigmas);
      msg->put(numGoAtoms*numGoAtoms*sizeof(bool), (char*)goWithinCutoff);
      // printf("Molecule.C sending atomChainTypes %d %d \n", numGoAtoms, atomChainTypes);
      break;
    case 2: //GSS
      msg->put(numGoAtoms);
      msg->put(numAtoms,pointerToGoBeg);
      msg->put(numAtoms,pointerToGoEnd);
      msg->put(numAtoms,goSigmaIndices);
      msg->put(numAtoms,goResidIndices);
      msg->put(numGoAtoms,atomChainTypes);
      msg->put(goNumLJPair);
      msg->put(goNumLJPair,goIndxLJA);
      msg->put(goNumLJPair,goIndxLJB);
      msg->put(goNumLJPair,goSigmaPairA);
      msg->put(goNumLJPair,goSigmaPairB);
      break;
    case 3:
      msg->put(numGoAtoms);
      msg->put(numAtoms, goSigmaIndices);
      msg->put(numGoAtoms, atomChainTypes);
      //msg->put(numGoAtoms*numGoAtoms, goSigmas);
      //msg->put(numGoAtoms*numGoAtoms*sizeof(bool), (char*)goWithinCutoff);
      msg->put(numGoAtoms*3, goCoordinates);
      msg->put(numGoAtoms, goResids);
      break;
    }
  } 

  msg->end();
  delete msg;
}
/*      END OF FUNCTION send_GoMolecule     */

// JLai
/************************************************************************/
/*                  */
/*      FUNCTION receive_Molecule      */
/*                  */
/*  receive_Molecule is used by all the clients to receive the  */
/*   Go structural data sent out by the master node.  It is NEVER called   */
/*   by the Master node.            */
/*                  */
/************************************************************************/
void Molecule::receive_GoMolecule(MIStream *msg) {
      // Ported by JLai -- Original by JE
      // JE - receive Go info
      Real *a1, *a2, *a3, *a4;
      int *i1, *i2, *i3, *i4;
      int maxGoChainsSqr = MAX_GO_CHAINS*MAX_GO_CHAINS;  // JE JLai Go code
      msg->get(NumGoChains);
      
      if (NumGoChains) {
	//go_indices = new int[MAX_GO_CHAINS+1];
	//go_array = new GoValue[MAX_GO_CHAINS*MAX_GO_CHAINS];
	
	//      int go_indices[MAX_GO_CHAINS+1];        //  Indices from chainIDs to go_array      
	//      GoValue go_array[MAX_GO_CHAINS*MAX_GO_CHAINS];   //  Array of Go params
	msg->get(MAX_GO_CHAINS+1,go_indices);
	
	a1 = new Real[maxGoChainsSqr];
	a2 = new Real[maxGoChainsSqr];
	a3 = new Real[maxGoChainsSqr];
	a4 = new Real[maxGoChainsSqr];
	i1 = new int[maxGoChainsSqr];
	i2 = new int[maxGoChainsSqr];
	i3 = new int[maxGoChainsSqr];
	i4 = new int[maxGoChainsSqr*MAX_RESTRICTIONS];
	
	if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (a4 == NULL) || 
	     (i1 == NULL) || (i2 == NULL) || (i3 == NULL) || (i4 == NULL) )
	  {
	    NAMD_die("memory allocation failed in Molecule::send_Molecule");
	  }
	
	msg->get(maxGoChainsSqr, a1);
	msg->get(maxGoChainsSqr, a2);
	msg->get(maxGoChainsSqr, a3);
	msg->get(maxGoChainsSqr, a4);
	msg->get(maxGoChainsSqr, i1);
	msg->get(maxGoChainsSqr, i2);
	msg->get(maxGoChainsSqr, i3);
	msg->get(maxGoChainsSqr*MAX_RESTRICTIONS, i4);
	
	for (int i=0; i<maxGoChainsSqr; i++) {
	  go_array[i].epsilon = a1[i];
	  go_array[i].sigmaRep = a2[i];
	  go_array[i].epsilonRep = a3[i];
	  go_array[i].cutoff = a4[i];
	  go_array[i].exp_a = i1[i];
	  go_array[i].exp_b = i2[i];
	  go_array[i].exp_rep = i3[i];
	  for (int j=0; j<MAX_RESTRICTIONS; j++) {
	    go_array[i].restrictions[j] = i4[i*MAX_RESTRICTIONS + j];
	  }
	}
	
	delete [] a1;
	delete [] a2;
	delete [] a3;
	delete [] a4;
	delete [] i1;
	delete [] i2;
	delete [] i3;
	delete [] i4;

	//msg->get(MAX_GO_CHAINS*MAX_GO_CHAINS, (char*)go_array);
	
	/*DebugM(3,"NumGoChains:" << NumGoChains << std::endl);
	  for (int ii=0; ii<MAX_GO_CHAINS; ii++) {
	  for (int jj=0; jj<MAX_GO_CHAINS; jj++) {
	  DebugM(3,"go_array[" << ii << "][" << jj << "]:" << go_array[ii][jj] << std::endl);
	  }
	  }
	  for (int ii=0; ii<MAX_GO_CHAINS+1; ii++) {
	  DebugM(3,"go_indices[" << ii << "]:" << go_indices[ii] << std::endl);
	  }*/
      }

      if (simParams->goForcesOn) {
	switch(simParams->goMethod) {
	case 1:
          msg->get(numGoAtoms);
	  //printf("Deleting goSigmaIndiciesA\n");
          delete [] goSigmaIndices;
          goSigmaIndices = new int32[numAtoms];
	  //printf("Deleting atomChainTypesA\n");
          delete [] atomChainTypes;
          atomChainTypes = new int32[numGoAtoms];
	  //printf("Deleting goSigmasA\n");
          delete [] goSigmas;
          goSigmas = new Real[numGoAtoms*numGoAtoms];
	  //printf("Deleting goWithinCutoffA\n"); 
          delete [] goWithinCutoff;
          goWithinCutoff = new bool[numGoAtoms*numGoAtoms];
          msg->get(numAtoms, goSigmaIndices);
          msg->get(numGoAtoms, atomChainTypes);
          msg->get(numGoAtoms*numGoAtoms, goSigmas);
          msg->get(numGoAtoms*numGoAtoms*sizeof(bool), (char*)goWithinCutoff);
          break;	  
	case 2: //GSR
	  msg->get(numGoAtoms);
	  delete [] pointerToGoBeg;
	  pointerToGoBeg = new int[numAtoms];
	  msg->get(numAtoms,pointerToGoBeg);
	  delete [] pointerToGoEnd;
	  pointerToGoEnd = new int[numAtoms];
	  msg->get(numAtoms,pointerToGoEnd);
	  delete [] goSigmaIndices;
	  goSigmaIndices = new int32[numAtoms];
	  msg->get(numAtoms,goSigmaIndices);
	  delete [] goResidIndices;
	  goResidIndices = new int32[numAtoms];
	  msg->get(numAtoms,goResidIndices);	  
	  delete [] atomChainTypes;
	  atomChainTypes = new int32[numGoAtoms];
	  msg->get(numGoAtoms,atomChainTypes);
	  msg->get(goNumLJPair);
	  delete [] goIndxLJA;
	  goIndxLJA = new int[goNumLJPair];
	  msg->get(goNumLJPair,goIndxLJA);
	  delete [] goIndxLJB;
	  goIndxLJB = new int[goNumLJPair];
	  msg->get(goNumLJPair,goIndxLJB);
	  delete [] goSigmaPairA;
	  goSigmaPairA = new double[goNumLJPair];
	  msg->get(goNumLJPair,goSigmaPairA);
	  delete [] goSigmaPairB;
	  goSigmaPairB = new double[goNumLJPair];
	  msg->get(goNumLJPair,goSigmaPairB);  
	  break;
	case 3:
	  msg->get(numGoAtoms);
	  //printf("Deleting goSigmaIndiciesB\n");
	  delete [] goSigmaIndices;
	  goSigmaIndices = new int32[numAtoms];
	  //printf("Deleting atomChainTypesB\n");
	  delete [] atomChainTypes;
	  atomChainTypes = new int32[numGoAtoms];
	  //delete [] goSigmas;
	  //goSigmas = new Real[numGoAtoms*numGoAtoms];
	  //delete [] goWithinCutoff;
	  //goWithinCutoff = new bool[numGoAtoms*numGoAtoms];
	  //printf("Deleting goCoordinatesB\n");
	  delete [] goCoordinates;
	  goCoordinates = new Real[numGoAtoms*3];
	  //printf("Deleting goResidsB\n");
	  delete [] goResids;
	  goResids = new int[numGoAtoms];
	  msg->get(numAtoms, goSigmaIndices);
	  msg->get(numGoAtoms, atomChainTypes);
	  //msg->get(numGoAtoms*numGoAtoms, goSigmas);
	  //msg->get(numGoAtoms*numGoAtoms*sizeof(bool), (char*)goWithinCutoff);
	  msg->get(numGoAtoms*3, goCoordinates);
	  msg->get(numGoAtoms, goResids);
	  break;
	}
      }

      delete msg;

}
/*      END OF FUNCTION receive_GoMolecule     */
