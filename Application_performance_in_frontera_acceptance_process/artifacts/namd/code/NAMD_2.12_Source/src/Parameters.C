/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   The class Parameters is used to hold all of the parameters read
   in from the parameter files.  The class provides a routine to read in
   parameter files (as many parameter files as desired can be read in) and
   a series of routines that allow the parameters that have been read in
   to be queried.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#ifndef WIN32
#include <strings.h>
#endif
#include "InfoStream.h"
#include <charm++.h>
#include "Parameters.h"
#include "Communicate.h"
#include "ConfigList.h"
//****** BEGIN CHARMM/XPLOR type changes
#include "SimParameters.h"
//****** END CHARMM/XPLOR type changes

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#define INDEX(ncols,i,j)  ((i)*ncols + (j))

#define ENABLETABLES

static char** table_types;

//  struct bond_params is used to form a binary tree of bond parameters.
//  The two atom names are used to determine the order of the nodes in the
//  tree.  atom1name should ALWAYS be lexically before atom2name

struct bond_params
{
  char atom1name[11];
  char atom2name[11];
  Real forceconstant;
  Real distance;
  Index index;
  struct bond_params *left;
  struct bond_params *right;
};

//  struct angle_params is used to form a binary tree of bond parameters.
//  The three atom names are used to determine the order of the nodes in
//  the tree.  atom1name should ALWAYS be lexically before atom3name

struct angle_params
{
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  Real forceconstant;
  int normal;
  Real angle;
  Real k_ub;
  Real r_ub;
  Index index;
  struct angle_params *left;
  struct angle_params *right;
};

//  struct dihedral_params is used to form a linked list of the dihedral
//  parameters.  The linked list is arranged in such a way that any
//  bonds with wildcards are at the end of the list so that a linear
//  search can be done but we will still find exact matches before
//  wildcard matches

struct dihedral_params
{
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  char atom4name[11];
  char atom1wild;
  char atom2wild;
  char atom3wild;
  char atom4wild;
  int multiplicity;
  FourBodyConsts values[MAX_MULTIPLICITY];
  Index index;
  struct dihedral_params *next;
  dihedral_params() { memset(this, 0, sizeof(dihedral_params)); }
};

//  struct improper_params is used to form a linked list of the improper
//  parameters.  The linked list is arranged in such a way that any
//  bonds with wildcards are at the end of the list so that a linear
//  search can be done but we will still find exact matches before
//  wildcard matches

struct improper_params
{
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  char atom4name[11];
  int multiplicity;
  FourBodyConsts values[MAX_MULTIPLICITY];
  Index index;
  struct improper_params *next;
};

struct crossterm_params
{
  crossterm_params(int dim) : dimension(dim) {
    values = new double[dimension*dimension];
  }
  ~crossterm_params() {
    delete [] values;
  }
  char atom1name[11];
  char atom2name[11];
  char atom3name[11];
  char atom4name[11];
  char atom5name[11];
  char atom6name[11];
  char atom7name[11];
  char atom8name[11];
  int dimension;  // usually 24
  double *values;  // dimension * dimension data
  Index index;
  struct crossterm_params *next;
};

//  struct vdw_params is used to form a binary serach tree of the
//  vdw paramters for a single atom.

struct vdw_params
{
  char atomname[11];
  Real sigma;
  Real epsilon;
  Real sigma14;
  Real epsilon14;
  Index index;
  struct vdw_params *left;
  struct vdw_params *right;
};

//  struct vdw_pair_params is used to form a linked list of the
//  vdw parameters for a pair of atoms

struct vdw_pair_params
{
  char atom1name[11];
  char atom2name[11];
  Real A;
  Real B;
  Real A14;
  Real B14;
  struct vdw_pair_params *next;
};

struct table_pair_params
{
  char atom1name[11];
  char atom2name[11];
  int K;
  struct table_pair_params *next;
};

struct nbthole_pair_params
{
  char atom1name[11];
  char atom2name[11];
  Real alphai;
  Real alphaj;
  Real tholeij;
  Index index;
  struct nbthole_pair_params *next;
};

Parameters::Parameters() {
  initialize();
}

void Parameters::initialize() {

  paramType = -1;

  /*  Set all the pointers to NULL        */
  atomTypeNames=NULL;
  bondp=NULL;
  anglep=NULL;
  improperp=NULL;
  dihedralp=NULL;
  crosstermp=NULL;
  vdwp=NULL;
  vdw_pairp=NULL;
  nbthole_pairp=NULL;
  table_pairp=NULL;
  bond_array=NULL;
  angle_array=NULL;
  dihedral_array=NULL;
  improper_array=NULL;
  crossterm_array=NULL;
  // JLai
  gromacsPair_array=NULL;
  // End of JLai
  vdw_array=NULL;
  vdw_pair_tree=NULL;
  nbthole_pair_tree=NULL;
  tab_pair_tree=NULL;
  maxDihedralMults=NULL;
  maxImproperMults=NULL;
  table_ener = NULL;

  /*  Set all the counts to 0          */
  NumBondParams=0;
  NumAngleParams=0;
  NumDihedralParams=0;
  NumImproperParams=0;
  NumCrosstermParams=0;
  // JLai
  NumGromacsPairParams=0;
  // End of JLai
  NumVdwParams=0;
  NumVdwPairParams=0;
  NumNbtholePairParams=0;
  NumTablePairParams=0;
  NumCosAngles=0;
  numenerentries=0;
}

/************************************************************************/
/*                  */
/*      FUNCTION Parameters        */
/*                  */
/*  This is the constructor for the class.  It simlpy sets the      */
/*  pointers to the list and trees to NULL and the count of all the     */
/*  parameters to 0.              */
/*  The type (format) of the input parameters (Xplor,Charmm) is set here. */
/*                  */
/************************************************************************/

Parameters::Parameters(SimParameters *simParams, StringList *f)
{
  initialize();

  //// get current parameter format
  if (simParams->paraTypeXplorOn)
  {
    paramType = paraXplor;
  }
  else if (simParams->paraTypeCharmmOn)
  {
    paramType = paraCharmm;
  }
  //****** END CHARMM/XPLOR type changes
  //Test for cos-based angles
  if (simParams->cosAngles) {
    cosAngles = true;
  } else {
    cosAngles = false;
  }

  if (simParams->tabulatedEnergies) {
	  CkPrintf("Working on tables\n");
	  read_ener_table(simParams);
  }

  //****** BEGIN CHARMM/XPLOR type changes
  /* Set up AllFilesRead flag to FALSE.  Once all of the files    */
  /* have been read in, then this will be set to true and the     */
  /* arrays of parameters will be set up        */
  AllFilesRead = FALSE;

  if (NULL != f) 
  {
    do
    {
      //****** BEGIN CHARMM/XPLOR type changes
      if (paramType == paraXplor)
      {
        read_parameter_file(f->data);
      }
      else if (paramType == paraCharmm)
      {
        read_charmm_parameter_file(f->data);
      }
      //****** END CHARMM/XPLOR type changes
      f = f->next;
    } while ( f != NULL );

    done_reading_files();
  }

}
/*      END OF FUNCTION Parameters      */

/************************************************************************/
/*                  */
/*      FUNCTION ~Parameters        */
/*                  */
/*  This is the destructor for this class.  It basically just       */
/*  frees all of the memory allocated for the parameters.    */
/*                  */
/************************************************************************/

Parameters::~Parameters()

{
        if (atomTypeNames)
          delete [] atomTypeNames;

  if (bondp != NULL)
    free_bond_tree(bondp);

  if (anglep != NULL)
    free_angle_tree(anglep);

  if (dihedralp != NULL)
    free_dihedral_list(dihedralp);

  if (improperp != NULL)
    free_improper_list(improperp);

  if (crosstermp != NULL)
    free_crossterm_list(crosstermp);

  if (vdwp != NULL)
    free_vdw_tree(vdwp);

  if (vdw_pairp != NULL)
    free_vdw_pair_list();

  if (nbthole_pairp != NULL)
    free_nbthole_pair_list();

  if (bond_array != NULL)
    delete [] bond_array;

  if (angle_array != NULL)
    delete [] angle_array;

  if (dihedral_array != NULL)
    delete [] dihedral_array;

  if (improper_array != NULL)
    delete [] improper_array;

  if (crossterm_array != NULL)
    delete [] crossterm_array;

  // JLai
  if (gromacsPair_array != NULL)
    delete [] gromacsPair_array;
  // End of JLai

  if (vdw_array != NULL)
    delete [] vdw_array;
  
  if (tab_pair_tree != NULL)
    free_table_pair_tree(tab_pair_tree);

  if (vdw_pair_tree != NULL)
    free_vdw_pair_tree(vdw_pair_tree);

  if (nbthole_pair_tree != NULL)
    free_nbthole_pair_tree(nbthole_pair_tree);

  if (maxDihedralMults != NULL)
    delete [] maxDihedralMults;

  if (maxImproperMults != NULL)
    delete [] maxImproperMults;

  for( int i = 0; i < error_msgs.size(); ++i ) {
    delete [] error_msgs[i];
  }
  error_msgs.resize(0);
}
/*      END OF FUNCTION ~Parameters      */

/************************************************************************/
/*                  */
/*      FUNCTION read_paramter_file      */
/*                  */
/*   INPUTS:                */
/*  fname - name of the parameter file to read      */
/*                  */
/*  This function reads in a parameter file and adds the parameters */
/*   from this file to the current group of parameters.  The basic      */
/*   logic of the routine is to read in a line from the file, looks at  */
/*   the first word of the line to determine what kind of parameter we  */
/*   have, and then call the appropriate routine to add the parameter   */
/*   to the parameters that we have.          */
/*                  */
/************************************************************************/

void Parameters::read_parameter_file(char *fname)

{
  char buffer[512];  //  Buffer to store each line of the file
  char first_word[512];  //  First word of the current line
  FILE *pfile;    //  File descriptor for the parameter file

  /*  Check to make sure that we haven't previously been told     */
  /*  that all the files were read        */
  if (AllFilesRead)
  {
    NAMD_die("Tried to read another parameter file after being told that all files were read!");
  }

  /*  Try and open the file          */
  if ( (pfile = Fopen(fname, "r")) == NULL)
  {
    char err_msg[256];

    sprintf(err_msg, "UNABLE TO OPEN XPLOR PARAMETER FILE %s\n", fname);
    NAMD_die(err_msg);
  }

  /*  Keep reading in lines until we hit the EOF      */
  while (NAMD_read_line(pfile, buffer) != -1)
  {
    /*  Get the first word of the line      */
    NAMD_find_first_word(buffer, first_word);

    /*  First, screen out things that we ignore, such as    */
    /*  blank lines, lines that start with '!', lines that  */
    /*  start with "REMARK", lines that start with set",    */
    /*  and most of the HBOND parameters which include      */
    /*  AEXP, REXP, HAEX, AAEX, but not the HBOND statement */
    /*  which is parsed.                                    */
    if ((buffer[0] != '!') && 
        !NAMD_blank_string(buffer) &&
        (strncasecmp(first_word, "REMARK", 6) != 0) &&
        (strcasecmp(first_word, "set")!=0) &&
        (strncasecmp(first_word, "AEXP", 4) != 0) &&
        (strncasecmp(first_word, "REXP", 4) != 0) &&
        (strncasecmp(first_word, "HAEX", 4) != 0) &&
        (strncasecmp(first_word, "AAEX", 4) != 0) &&
        (strncasecmp(first_word, "NBOND", 5) != 0) &&
        (strncasecmp(first_word, "CUTNB", 5) != 0) &&
        (strncasecmp(first_word, "END", 3) != 0) &&
        (strncasecmp(first_word, "CTONN", 5) != 0) &&
        (strncasecmp(first_word, "EPS", 3) != 0) &&
        (strncasecmp(first_word, "VSWI", 4) != 0) &&
        (strncasecmp(first_word, "NBXM", 4) != 0) &&
        (strncasecmp(first_word, "INHI", 4) != 0) )
    {
      /*  Now, call the appropriate function based    */
      /*  on the type of parameter we have    */
      if (strncasecmp(first_word, "bond", 4)==0)
      {
        add_bond_param(buffer);
        NumBondParams++;
      }
      else if (strncasecmp(first_word, "angl", 4)==0)
      {
        add_angle_param(buffer);
        NumAngleParams++;
      }
      else if (strncasecmp(first_word, "dihe", 4)==0)
      {
        add_dihedral_param(buffer, pfile);
        NumDihedralParams++;
      }
      else if (strncasecmp(first_word, "impr", 4)==0)
      {
        add_improper_param(buffer, pfile);
        NumImproperParams++;
      }
      else if (strncasecmp(first_word, "nonb", 4)==0)
      {
        add_vdw_param(buffer);
        NumVdwParams++; 
      }
      else if (strncasecmp(first_word, "nbfi", 4)==0)
      {
        add_vdw_pair_param(buffer);
        NumVdwPairParams++; 
      }
      else if (strncasecmp(first_word, "nbta", 4)==0)
      {

        if (table_ener == NULL) {
          continue;
        }

        add_table_pair_param(buffer);
        NumTablePairParams++; 
      }
      else if (strncasecmp(first_word, "hbon", 4)==0)
      {
        add_hb_pair_param(buffer);
      }
      else
      {
        /*  This is an unknown paramter.        */
        /*  This is BAD        */
        char err_msg[512];

        sprintf(err_msg, "UNKNOWN PARAMETER IN XPLOR PARAMETER FILE %s\nLINE=*%s*",
           fname, buffer);
        NAMD_die(err_msg);
      }
    }
  }

  /*  Close the file            */
  Fclose(pfile);

  return;
}
/*      END OF FUNCTION read_paramter_file    */

//****** BEGIN CHARMM/XPLOR type changes
/************************************************************************/
/*                                                                        */
/*                        FUNCTION read_charmm_paramter_file                */
/*                                                                        */
/*   INPUTS:                                                                */
/*        fname - name of the parameter file to read                        */
/*                                                                        */
/*        This function reads in a CAHRMM parameter file and adds the     */ 
/*   parameters from this file to the current group of parameters.      */
/*   The basic logic of the routine is to first find out what type of   */
/*   parameter we have in the file. Then look at each line in turn      */
/*   and call the appropriate routine to add the parameters until we hit*/
/*   a new type of parameter or EOF.                                    */
/*                                                                        */
/************************************************************************/

void Parameters::read_charmm_parameter_file(char *fname)

{
  int  par_type=0;         //  What type of parameter are we currently
                           //  dealing with? (vide infra)
  int  skipline;           //  skip this line?
  int  skipall = 0;        //  skip rest of file;
  char buffer[512];           //  Buffer to store each line of the file
  char first_word[512];           //  First word of the current line
  FILE *pfile;                   //  File descriptor for the parameter file

  /*  Check to make sure that we haven't previously been told     */
  /*  that all the files were read                                */
  if (AllFilesRead)
  {
    NAMD_die("Tried to read another parameter file after being told that all files were read!");
  }

  /*  Try and open the file                                        */
  if ( (pfile = fopen(fname, "r")) == NULL)
  {
    char err_msg[256];

    sprintf(err_msg, "UNABLE TO OPEN CHARMM PARAMETER FILE %s\n", fname);
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
	iout << iWARN << "SKIPPING PART OF PARAMETER FILE AFTER RETURN STATEMENT\n" << endi;
        break;
      }
      /*  Now, determine the apropriate parameter type.   */
      if (strncasecmp(first_word, "bond", 4)==0)
      {
        par_type=1; skipline=1;
      }
      else if (strncasecmp(first_word, "angl", 4)==0)
      {
        par_type=2; skipline=1;
      }
      else if (strncasecmp(first_word, "thet", 4)==0)
      {
        par_type=2; skipline=1;
      }
      else if (strncasecmp(first_word, "dihe", 4)==0)
      {
        par_type=3; skipline=1;
      }
      else if (strncasecmp(first_word, "phi", 3)==0)
      {
        par_type=3; skipline=1;
      }
      else if (strncasecmp(first_word, "impr", 4)==0)
      {
        par_type=4; skipline=1;
      }
      else if (strncasecmp(first_word, "imph", 4)==0)
      {
        par_type=4; skipline=1;
      }
      else if (strncasecmp(first_word, "nbon", 4)==0)
      {
        par_type=5; skipline=1;
      }
      else if (strncasecmp(first_word, "nonb", 4)==0)
      {
        par_type=5; skipline=1;
      }
      else if (strncasecmp(first_word, "nbfi", 4)==0)
      {
        par_type=6; skipline=1;
      }
      else if (strncasecmp(first_word, "hbon", 4)==0)
      {
        par_type=7; skipline=1;
      }
      else if (strncasecmp(first_word, "cmap", 4)==0)
      {
        par_type=8; skipline=1;
      }
      else if (strncasecmp(first_word, "nbta", 4)==0)
      {
        par_type=9; skipline=1;
      }
      else if (strncasecmp(first_word, "thol", 4)==0)
      {
        par_type=10; skipline=1;
      }
      else if (strncasecmp(first_word, "atom", 4)==0)
      {
        par_type=11; skipline=1;
      }
      else if (strncasecmp(first_word, "ioformat", 8)==0)
      {
        skipline=1;
      }
      else if (strncasecmp(first_word, "read", 4)==0)
      {
        skip_stream_read(buffer, pfile);  skipline=1;
      }
      else if (strncasecmp(first_word, "return", 4)==0)
      {
        skipall=1;  skipline=1;
      }
      else if ((strncasecmp(first_word, "nbxm", 4) == 0) ||
               (strncasecmp(first_word, "grou", 4) == 0) ||
               (strncasecmp(first_word, "cdie", 4) == 0) ||
               (strncasecmp(first_word, "shif", 4) == 0) ||
               (strncasecmp(first_word, "vgro", 4) == 0) ||
               (strncasecmp(first_word, "vdis", 4) == 0) ||
               (strncasecmp(first_word, "vswi", 4) == 0) ||
               (strncasecmp(first_word, "cutn", 4) == 0) ||
               (strncasecmp(first_word, "ctof", 4) == 0) ||
               (strncasecmp(first_word, "cton", 4) == 0) ||
               (strncasecmp(first_word, "eps", 3) == 0) ||
               (strncasecmp(first_word, "e14f", 4) == 0) ||
               (strncasecmp(first_word, "wmin", 4) == 0) ||
               (strncasecmp(first_word, "aexp", 4) == 0) ||
               (strncasecmp(first_word, "rexp", 4) == 0) ||
               (strncasecmp(first_word, "haex", 4) == 0) ||
               (strncasecmp(first_word, "aaex", 4) == 0) ||
               (strncasecmp(first_word, "noac", 4) == 0) ||
               (strncasecmp(first_word, "hbno", 4) == 0) ||
               (strncasecmp(first_word, "cuth", 4) == 0) ||
               (strncasecmp(first_word, "ctof", 4) == 0) ||
               (strncasecmp(first_word, "cton", 4) == 0) ||
               (strncasecmp(first_word, "cuth", 4) == 0) ||
               (strncasecmp(first_word, "ctof", 4) == 0) ||
               (strncasecmp(first_word, "cton", 4) == 0) ) 
      {
        if ((par_type != 5) && (par_type != 6) && (par_type != 7) && (par_type != 9))
        {
          char err_msg[512];

          sprintf(err_msg, "ERROR IN CHARMM PARAMETER FILE %s\nLINE=*%s*",fname, buffer);
          NAMD_die(err_msg);
        }
        else 
        {
          skipline = 1;
        }
      }        
      else if (par_type == 0)
      {
        /*  This is an unknown paramter.        */
        /*  This is BAD                                */
        char err_msg[512];

        sprintf(err_msg, "UNKNOWN PARAMETER IN CHARMM PARAMETER FILE %s\nLINE=*%s*",fname, buffer);
        NAMD_die(err_msg);
      }
    }
    else
    {
      skipline=1;
    }

    if ( (par_type != 0) && (!skipline) )
    {
      /*  Now, call the appropriate function based    */
      /*  on the type of parameter we have                */
      /*  I know, this should really be a switch ...  */
      if (par_type == 1)
      {
        add_bond_param(buffer);
        NumBondParams++;
      }
      else if (par_type == 2)
      {
        add_angle_param(buffer);
        NumAngleParams++;
      }
      else if (par_type == 3)
      {
        add_dihedral_param(buffer, pfile);
        NumDihedralParams++;
      }
      else if (par_type == 4)
      {
        add_improper_param(buffer, pfile);
        NumImproperParams++;
      }
      else if (par_type == 5)
      {
        add_vdw_param(buffer);
        NumVdwParams++;
      }
      else if (par_type == 6)
      {
        add_vdw_pair_param(buffer);
        NumVdwPairParams++; 
      }
      else if (par_type == 7)
      {
        add_hb_pair_param(buffer);                  
      }
      else if (par_type == 8)
      {
        add_crossterm_param(buffer, pfile);                  
        NumCrosstermParams++;
      }
      else if (par_type == 9)
      {

        if (table_ener == NULL) {
          continue;
        }

        add_table_pair_param(buffer);                  
        NumTablePairParams++;
      }
      else if (par_type == 10)
      {
        add_nbthole_pair_param(buffer);
        NumNbtholePairParams++;
      }
      else if (par_type == 11)
      {
        if ( strncasecmp(first_word, "mass", 4) != 0 ) {
          char err_msg[512];
          sprintf(err_msg, "UNKNOWN PARAMETER IN CHARMM PARAMETER FILE %s\nLINE=*%s*",fname, buffer);
          NAMD_die(err_msg);
        }
      }
      else
      {
        /*  This really should not occour!      */
        /*  This is an internal error.          */
        /*  This is VERY BAD                        */
        char err_msg[512];

        sprintf(err_msg, "INTERNAL ERROR IN CHARMM PARAMETER FILE %s\nLINE=*%s*",fname, buffer);
        NAMD_die(err_msg);
      }
    }
  }

  /*  Close the file                                                */
  fclose(pfile);

  return;
}
/*                        END OF FUNCTION read_charmm_paramter_file                */
//****** END CHARMM/XPLOR type changes


void Parameters::skip_stream_read(char *buf, FILE *fd) {

  char buffer[513];
  char first_word[513];
  char s1[128];
  char s2[128];
  int rval = sscanf(buf, "%s %s", s1, s2);
  if (rval != 2) {
        char err_msg[512];
        sprintf(err_msg, "BAD FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*", buf);
        NAMD_die(err_msg);
  }
  if ( ! strncasecmp(s2,"PARA",4) ) return;  // read parameters
  
  iout << iINFO << "SKIPPING " << s2 << " SECTION IN STREAM FILE\n" << endi;

  while (NAMD_read_line(fd, buffer) != -1)
  {
    // read until we find "END"
    NAMD_find_first_word(buffer, first_word);
    if (!NAMD_blank_string(buffer) &&
        (strncmp(first_word, "!", 1) != 0) &&
         (strncmp(first_word, "*", 1) != 0) &&
         (strncasecmp(first_word, "END", 3) == 0) ) return;
  }

}


/************************************************************************/
/*                  */
/*      FUNCTION add_bond_param        */
/*                  */
/*   INPUTS:                */
/*  buf - Line from parameter file containing bond parameters  */
/*                  */
/*  This function adds a new bond paramter to the binary tree of    */
/*   angle paramters that we have.  If a duplicate is found, a warning  */
/*   message is printed and the new parameters are used.    */
/*                  */
/************************************************************************/

void Parameters::add_bond_param(char *buf)

{
  char atom1name[11];    //  Atom type for atom 1
  char atom2name[11];    //  Atom type for atom 2
  Real forceconstant;    //  Force constant for bond
  Real distance;      //  Rest distance for bond
  int read_count;      //  Count from sscanf
  struct bond_params *new_node;  //  New node in tree

  //****** BEGIN CHARMM/XPLOR type changes
  /*  Use sscanf to parse up the input line      */
  if (paramType == paraXplor)
  {
    /* read XPLOR format */
    read_count=sscanf(buf, "%*s %s %s %f %f\n", atom1name, atom2name, 
       &forceconstant, &distance);
  }
  else if (paramType == paraCharmm)
  {
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %s %f %f\n", atom1name, atom2name, 
       &forceconstant, &distance);
  }
  //****** END CHARMM/XPLOR type changes

  /*  Check to make sure we found everything we expeceted    */
  if (read_count != 4)
  {
    char err_msg[512];

    if (paramType == paraXplor)
    {
      sprintf(err_msg, "BAD BOND FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buf);
    }
    else if (paramType == paraCharmm)
    {
      sprintf(err_msg, "BAD BOND FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*\n", buf);
    }
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new bond_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_bond_param\n");
  }

  /*  Order the atoms so that the atom that comes alphabetically  */
  /*  first is atom 1.  Since the bond is symmetric, it doesn't   */
  /*  matter physically which atom is first.  And this allows the */
  /*  search of the binary tree to be done in a logical manner    */
  if (strcasecmp(atom1name, atom2name) < 0)
  {
    strcpy(new_node->atom1name, atom1name);
    strcpy(new_node->atom2name, atom2name);
  }
  else
  {
    strcpy(new_node->atom2name, atom1name);
    strcpy(new_node->atom1name, atom2name);
  }

  /*  Assign force constant and distance        */
  new_node->forceconstant = forceconstant;
  new_node->distance = distance;

  /*  Set pointers to null          */
  new_node->left = NULL;
  new_node->right = NULL;

  /*  Make call to recursive call to actually add the node to the */
  /*  tree              */
  bondp=add_to_bond_tree(new_node, bondp);

  return;
}
/*      END OF FUNCTION add_bond_param      */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_bond_tree      */
/*                  */
/*   INPUTS:                */
/*  new_node - Node to add to the tree        */
/*  tree - tree to add the node to          */
/*                  */
/*   OUTPUTS:                */
/*  ths function returns a pointer to the new tree with the node    */
/*   added to it.  Most of the time it will be the same pointer as was  */
/*   passed in, but not if the current tree is empty.      */
/*                  */
/*  this is a receursive function that adds a node to the binary    */
/*   tree used to store bond parameters.        */
/*                  */
/************************************************************************/

struct bond_params *Parameters::add_to_bond_tree(struct bond_params *new_node,
             struct bond_params *tree)

{
  int compare_code;  //  Results from strcasecmp

  /*  If the tree is currently empty, then the new tree consists  */
  /*  only of the new node          */
  if (tree == NULL)
    return(new_node);

  /*  Compare the atom1 name from the new node and the head of    */
  /*  the tree              */
  compare_code = strcasecmp(new_node->atom1name, tree->atom1name);

  /*  Check to see if they are the same        */
  if (compare_code == 0)
  {
    /*  The atom 1 names are the same, compare atom 2  */
    compare_code = strcasecmp(new_node->atom2name, tree->atom2name);

    /*  If atom 1 AND atom 2 are the same, we have a duplicate */
    if (compare_code == 0)
    {
      /*  We have a duplicate.  So print out a warning*/
      /*  message.  Then assign the new values to the */
      /*  tree and free the new_node      */
      //****** BEGIN CHARMM/XPLOR type changes
      /* we do not care about identical replacement */
      if ((tree->forceconstant != new_node->forceconstant) || 
          (tree->distance != new_node->distance))
      {
        iout << "\n" << iWARN << "DUPLICATE BOND ENTRY FOR "
          << new_node->atom1name << "-"
          << new_node->atom2name
          << "\nPREVIOUS VALUES  k=" << tree->forceconstant
          << "  x0=" << tree->distance
          << "\n   USING VALUES  k=" << new_node->forceconstant
          << "  x0=" << new_node->distance
          << "\n" << endi;

        tree->forceconstant=new_node->forceconstant;
        tree->distance=new_node->distance;
      }
      //****** END CHARMM/XPLOR type changes

      delete new_node;

      return(tree);
    }
  }

  /*  We don't have a duplicate, so if the new value is less      */
  /*  than the head of the tree, add it to the left child,   */
  /*  otherwise add it to the right child        */
  if (compare_code < 0)
  {
    tree->left = add_to_bond_tree(new_node, tree->left);
  }
  else
  {
    tree->right = add_to_bond_tree(new_node, tree->right);
  }

  return(tree);
}
/*    END OF FUNCTION add_to_bond_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION add_angle_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with angle parameters    */
/*                  */
/*  this function adds an angle parameter.  It parses up the input  */
/*   line and then adds it to the binary tree used to store the angle   */
/*   parameters.              */
/*                  */
/************************************************************************/

void Parameters::add_angle_param(char *buf)

{
  char atom1name[11];    // Type for atom 1
  char atom2name[11];    // Type for atom 2
  char atom3name[11];    // Type for atom 3
  char norm[4]="xxx";
  Real forceconstant;    // Force constant
  Real angle;      // Theta 0
  Real k_ub;      // Urey-Bradley force constant
  Real r_ub;      // Urey-Bradley distance
  int read_count;      // count from sscanf
  struct angle_params *new_node;  // new node in tree

  //****** BEGIN CHARMM/XPLOR type changes
  /*  parse up the input line with sscanf        */
  if (paramType == paraXplor)
  {
    /* read XPLOR format */
    read_count=sscanf(buf, "%*s %s %s %s %f %f UB %f %f\n", 
       atom1name, atom2name, atom3name, &forceconstant, &angle,
       &k_ub, &r_ub);
  }
  else if ((paramType == paraCharmm) && cosAngles) {
    read_count=sscanf(buf, "%s %s %s %f %f %3s %f %f\n", 
       atom1name, atom2name, atom3name, &forceconstant, &angle, norm,
       &k_ub, &r_ub);
//    printf("%s\n", buf);
//    printf("Data: %s %s %s %f %f %s %f %f\n", atom1name, atom2name, atom3name, forceconstant, angle, norm, k_ub, r_ub);
  }  
  else if (paramType == paraCharmm)
  {
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %s %s %f %f %f %f\n", 
       atom1name, atom2name, atom3name, &forceconstant, &angle,
       &k_ub, &r_ub);
//    printf("%s\n", buf);
//    printf("Data: %s %s %s %f %f\n", atom1name, atom2name, atom3name, forceconstant, angle);
  }
  //****** END CHARMM/XPLOR type changes

  /*  Check to make sure we got what we expected      */
  if ( (read_count != 5) && (read_count != 7) && !(cosAngles && read_count == 6))
  {
    char err_msg[512];

    if (paramType == paraXplor)
    {
      sprintf(err_msg, "BAD ANGLE FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buf);
    }
    else if (paramType == paraCharmm)
    {
      sprintf(err_msg, "BAD ANGLE FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*\n", buf);
    }
    NAMD_die(err_msg);
  }

  /*  Allocate the new node          */
  new_node = new angle_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_angle_param");
  }

  /*  As with the bond, we want the atom type is comes first  */
  /*  alphbetically first between atom 1 and atom 3 to be in      */
  /*  atom 1 so that we can search the tree reliably.    */
  if (strcasecmp(atom1name, atom3name) < 0)
  {
    strcpy(new_node->atom1name, atom1name);
    strcpy(new_node->atom2name, atom2name);
    strcpy(new_node->atom3name, atom3name);
  }
  else
  {
    strcpy(new_node->atom3name, atom1name);
    strcpy(new_node->atom2name, atom2name);
    strcpy(new_node->atom1name, atom3name);
  }

  /*  Assign the constants and pointer values      */
  new_node->forceconstant = forceconstant;
  new_node->angle = angle;

  if (cosAngles) {
    if (strcasecmp("cos",norm)==0) {
//      iout << "Info: Using cos mode for angle " << buf << endl;
      NumCosAngles++;
      new_node->normal = 0;
    } else {
//      iout << "Info: Using x^2 mode for angle " << buf << endl;
      new_node->normal = 1;
    }
  } else {
    new_node->normal = 1;
  }

  if (read_count == 7)
  {
    //  Urey-Bradley constants
    if (new_node->normal == 0) {
      NAMD_die("ERROR: Urey-Bradley angles can't be used with cosine-based terms\n");
    }
    new_node->k_ub = k_ub;
    new_node->r_ub = r_ub;
  }
  else
  {
    new_node->k_ub = 0.0;
    new_node->r_ub = 0.0;
  }

  new_node->left = NULL;
  new_node->right = NULL;

  /*  Insert it into the tree          */
  anglep = add_to_angle_tree(new_node, anglep);

  return;
}
/*      END OF FUNCTION add_angle_param      */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_angle_tree      */
/*                  */
/*   INPUTS:                */
/*  new_node - new node to add to the angle tree      */
/*  tree - tree to add the node to          */
/*                  */
/*   OUTPUTS:                */
/*  the function returns a pointer to the new tree with the node    */
/*   added.  Most of the time, this will be the same as the value passed*/
/*   in, but not in the case where the tree is empty.      */
/*                  */
/*  this is a recursive function that adds an angle parameter  */
/*   to the binary tree storing the angle parameters.  If a duplicate   */
/*   is found, a warning message is printed, the current values in the  */
/*   tree are replaced with the new values, and the new node is free'd  */
/*                  */
/************************************************************************/

struct angle_params *Parameters::add_to_angle_tree(struct angle_params *new_node,
             struct angle_params *tree)

{
  int compare_code;  //  Return code from strcasecmp

  /*  If the tree is empty, then the new_node is the tree    */
  if (tree == NULL)
    return(new_node);

  /*  Compare atom 1 from the new node and the head of the tree   */
  compare_code = strcasecmp(new_node->atom1name, tree->atom1name);

  if (compare_code == 0)
  {
    /*  Atom 1 is the same, compare atom 2      */
    compare_code = strcasecmp(new_node->atom2name, tree->atom2name);

    if (compare_code == 0)
    {
      /*  Atoms 1 & 2 are the same, compare atom 3  */
      compare_code = strcasecmp(new_node->atom3name, 
            tree->atom3name);

      if (compare_code == 0)
      {
        /*  All three atoms were the same, this */
        /*  is a duplicate.  Print a warning    */
        /*  message, replace the current values,*/
        /*  and free the new node    */
        //****** BEGIN CHARMM/XPLOR type changes
        /* we do not care about identical replacement */
        if ((tree->forceconstant != new_node->forceconstant) ||
            (tree->angle != new_node->angle) ||
            (tree->k_ub != new_node->k_ub) ||
            (tree->r_ub != new_node->r_ub) || (tree->normal != new_node->normal))
        {
          iout << "\n" << iWARN << "DUPLICATE ANGLE ENTRY FOR "
            << new_node->atom1name << "-"
            << new_node->atom2name << "-"
            << new_node->atom3name
            << "\nPREVIOUS VALUES  k="
            << tree->forceconstant << "  theta0="
            << tree->angle << " k_ub="
            << tree->k_ub << " r_ub="
            << tree->r_ub
            << "\n   USING VALUES  k="
            << new_node->forceconstant << "  theta0="
            << new_node->angle << " k_ub="
            << new_node->k_ub << " r_ub=" << new_node->r_ub 
            << "\n" << endi;

          tree->forceconstant=new_node->forceconstant;
          tree->angle=new_node->angle;
          tree->k_ub=new_node->k_ub;
          tree->r_ub=new_node->r_ub;
          tree->normal=new_node->normal;
        }
        //****** END CHARMM/XPLOR type changes

        delete new_node;

        return(tree);
      }
    }
  }

  /*  Didn't find a duplicate, so if the new_node is smaller  */
  /*  than the current head, add it to the left child.  Otherwise */
  /*  add it to the right child.          */
  if (compare_code < 0)
  {
    tree->left = add_to_angle_tree(new_node, tree->left);
  }
  else
  {
    tree->right = add_to_angle_tree(new_node, tree->right);
  }

  return(tree);
}
/*      END OF FUNCTION add_to_angle_tree    */

/************************************************************************/
/*                  */
/*      FUNCTION add_dihedral_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with dihedral parameters    */
/*                  */
/*  this function adds an dihedral parameter.  It parses up the     */
/*   input line and then adds it to the binary tree used to store the   */
/*   dihedral parameters.            */
/*                  */
/************************************************************************/

void Parameters::add_dihedral_param(char *buf, FILE *fd)

{
  char atom1name[11];       //  Type of atom 1
  char atom2name[11];       //  Type of atom 2
  char atom3name[11];       //  Type of atom 3
  char atom4name[11];       //  Type of atom 4
  Real forceconstant;       //  Force constant
  int periodicity;       //  Periodicity
  Real phase_shift;       //  Phase shift
  int read_count;         //  Count from sscanf
  struct dihedral_params *new_node;  //  New node
  int multiplicity;       //  Multiplicity for bonds
  int i;           //  Loop counter
  char buffer[513];       //  Buffer for new line
  int ret_code;         //  Return code

  //****** BEGIN CHARMM/XPLOR type changes
  /*  Parse up the input line using sscanf      */
  if (paramType == paraXplor)
  {
    /* read XPLOR format */
    read_count=sscanf(buf, "%*s %s %s %s %s MULTIPLE= %d %f %d %f\n", 
       atom1name, atom2name, atom3name, atom4name, &multiplicity,
       &forceconstant, &periodicity, &phase_shift);
  }
  else if (paramType == paraCharmm)
  {
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %s %s %s %f %d %f\n", 
       atom1name, atom2name, atom3name, atom4name,
       &forceconstant, &periodicity, &phase_shift);
    multiplicity=1; 
  }

  if ( (read_count != 4) && (read_count != 8) && (paramType == paraXplor) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD DIHEDRAL FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }
  else if ( (read_count != 7) && (paramType == paraCharmm) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD DIHEDRAL FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }

  if ( (read_count == 4) && (paramType == paraXplor) )
  //****** END CHARMM/XPLOR type changes
  {
    read_count=sscanf(buf, "%*s %*s %*s %*s %*s %f %d %f\n", 
          &forceconstant, &periodicity, &phase_shift);

    /*  Check to make sure we got what we expected    */
    if (read_count != 3)
    {
      char err_msg[512];

      sprintf(err_msg, "BAD DIHEDRAL FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buf);
      NAMD_die(err_msg);
    }

    multiplicity = 1;
  }

  if (multiplicity > MAX_MULTIPLICITY)
  {
    char err_msg[181];

    sprintf(err_msg, "Multiple dihedral with multiplicity of %d greater than max of %d",
       multiplicity, MAX_MULTIPLICITY);
    NAMD_die(err_msg);
  }

  /*  Allocate new node            */
  new_node = new dihedral_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_dihedral_param\n");
  }

  /*  Assign all of the values for this node.  Notice that since  */
  /*  the dihedrals and impropers are implemented with a linked   */
  /*  list rather than a binary tree, we don't really care about  */
  /*  the order of the atoms any more        */
  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);
  strcpy(new_node->atom3name, atom3name);
  strcpy(new_node->atom4name, atom4name);
  new_node->atom1wild = ! strcasecmp(atom1name, "X");
  new_node->atom2wild = ! strcasecmp(atom2name, "X");
  new_node->atom3wild = ! strcasecmp(atom3name, "X");
  new_node->atom4wild = ! strcasecmp(atom4name, "X");
  new_node->multiplicity = multiplicity;
  if (paramType == paraXplor && periodicity != 0) phase_shift *= -1;
  new_node->values[0].k = forceconstant;
  new_node->values[0].n = periodicity;
  new_node->values[0].delta = phase_shift;

  new_node->next = NULL;

  //  If the multiplicity is greater than 1, then read in other parameters
  if (multiplicity > 1)
  {
    for (i=1; i<multiplicity; i++)
    {
      ret_code = NAMD_read_line(fd, buffer);

      //  Get rid of comments at the end of a line
      if (ret_code == 0)
      {
        NAMD_remove_comment(buffer);
      }

      //  Keep reading lines until we get one that isn't blank
      while ( (ret_code == 0) && (NAMD_blank_string(buffer)) )
      {
        ret_code = NAMD_read_line(fd, buffer);
      }

      if (ret_code != 0)
      {
        NAMD_die("EOF encoutner in middle of multiple dihedral");
      }

      read_count=sscanf(buffer, "%f %d %f\n", 
            &forceconstant, &periodicity, &phase_shift);

      if (read_count != 3)
      {
        char err_msg[512];

        sprintf(err_msg, "BAD MULTIPLE FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buffer);
        NAMD_die(err_msg);
      }

      if (paramType == paraXplor && periodicity != 0) phase_shift *= -1;
      new_node->values[i].k = forceconstant;
      new_node->values[i].n = periodicity;
      new_node->values[i].delta = phase_shift;
    }
  }

  //****** BEGIN CHARMM/XPLOR type changes
  /*  Add this node to the list          */
  if (paramType == paraXplor)
  {
    add_to_dihedral_list(new_node); // XPLOR
  }
  else if (paramType == paraCharmm)
  {
    add_to_charmm_dihedral_list(new_node); // CHARMM
  }
 //****** END CHARMM/XPLOR type changes

  return;
}
/*      END OF FUNCTION add_dihedral_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_dihedral_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node that is to be added to dihedral_list    */
/*                  */
/*  this function adds a new dihedral parameter to the linked list  */
/*   of dihedral parameters.  First, it checks for duplicates.  If a    */
/*   duplicate is found, a warning message is printed, the old values   */
/*   are replaced with the new values, and the new node is freed.  If   */
/*   Otherwise, the node is added to the list.  This list is arranged   */
/*   so that bods with wildcards are placed at the tail of the list.    */
/*   This will guarantee that if we just do a linear search, we will    */
/*   always find an exact match before a wildcard match.    */
/*                  */
/************************************************************************/

void Parameters::add_to_dihedral_list(
        struct dihedral_params *new_node)

{
  static struct dihedral_params *ptr;   //  position within list
  static struct dihedral_params *tail;  //  Pointer to the end of 
                //  the list so we can add
                //  entries to the end of the
                //  list in constant time
  int i;              //  Loop counter

  /*  If the list is currently empty, then the new node is the list*/
  if (dihedralp == NULL)
  {
    dihedralp=new_node;
    tail=new_node;

    return;
  }

  /*  The list isn't empty, so check for a duplicate    */
  ptr=dihedralp;

  while (ptr != NULL)
  {
    if ( ( (strcasecmp(new_node->atom1name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom4name, ptr->atom4name) == 0) ) ||
         ( (strcasecmp(new_node->atom4name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom1name, ptr->atom4name) == 0) ) )
    {
      /*  Found a duplicate        */
      //****** BEGIN CHARMM/XPLOR type changes
      /* we do not care about identical replacement */
      int echoWarn=0;  // echo warning messages ?

      if (ptr->multiplicity != new_node->multiplicity) {echoWarn=1;}
      
      if (!echoWarn)
      {
        for (i=0; i<ptr->multiplicity; i++)
        {
          if (ptr->values[i].k != new_node->values[i].k) {echoWarn=1; break;}
          if (ptr->values[i].n != new_node->values[i].n) {echoWarn=1; break;}
          if (ptr->values[i].delta != new_node->values[i].delta) {echoWarn=1; break;}
        }
      }

      if (echoWarn)
      {
        iout << "\n" << iWARN << "DUPLICATE DIHEDRAL ENTRY FOR "
          << ptr->atom1name << "-"
          << ptr->atom2name << "-"
          << ptr->atom3name << "-"
          << ptr->atom4name
          << "\nPREVIOUS VALUES MULTIPLICITY " << ptr->multiplicity << "\n";
        
        for (i=0; i<ptr->multiplicity; i++)
        {
          iout     << "  k=" << ptr->values[i].k
                   << "  n=" << ptr->values[i].n
                   << "  delta=" << ptr->values[i].delta;
        }

        iout << "\nUSING VALUES MULTIPLICITY " << new_node->multiplicity << "\n";

        for (i=0; i<new_node->multiplicity; i++)
        {
          iout <<     "  k=" << new_node->values[i].k
                   << "  n=" << new_node->values[i].n
                   << "  delta=" << new_node->values[i].delta;
        }

        iout << endi;

        ptr->multiplicity = new_node->multiplicity;

        for (i=0; i<new_node->multiplicity; i++)
        {
          ptr->values[i].k = new_node->values[i].k;
          ptr->values[i].n = new_node->values[i].n;
          ptr->values[i].delta = new_node->values[i].delta;
        }

      }
      //****** END CHARMM/XPLOR type changes

      delete new_node;

      return;
    }

    ptr=ptr->next;
  }

  /*  Check to see if we have any wildcards.  Since specific  */
  /*  entries are to take precedence, we'll put anything without  */
  /*  wildcards at the begining of the list and anything with     */
  /*  wildcards at the end of the list.  Then, we can just do a   */
  /*  linear search for a bond and be guaranteed to have specific */
  /*  entries take precendence over over wildcards          */
  if ( new_node->atom1wild ||
       new_node->atom2wild ||
       new_node->atom3wild ||
       new_node->atom4wild )
  {
    /*  add to the end of the list        */
    tail->next=new_node;
    tail=new_node;

    return;
  }
  else
  {
    /*  add to the head of the list        */
    new_node->next=dihedralp;
    dihedralp=new_node;

    return;
  }

}
/*    END OF FUNCTION add_to_dihedral_list      */

//****** BEGIN CHARMM/XPLOR type changes
/************************************************************************/
/*                                                                        */
/*                        FUNCTION add_to_charmm_dihedral_list                */
/*                                                                        */
/*   INPUTS:                                                                */
/*        new_node - node that is to be added to dihedral_list                */
/*                                                                        */
/*        this function adds a new dihedral parameter to the linked list  */
/*   of dihedral parameters in CHARMM format.                           */
/*   First, it checks for duplicates.  If a duplicate is found, a       */
/*   warning message is printed. If the periodicity is the same as of   */
/*   a previous dihedral the old values are replaced with the new       */
/*   values, otherwise, the dihedral is added and the multiplicity is   */
/*   increased.                                                         */
/*   Otherwise, the node is added to the list.  This list is arranged   */
/*   so that bonds with wildcards are placed at the tail of the list.   */
/*   This will guarantee that if we just do a linear search, we will    */
/*   always find an exact match before a wildcard match.                */
/*                                                                        */
/************************************************************************/

void Parameters::add_to_charmm_dihedral_list(
                                struct dihedral_params *new_node)

{
        static struct dihedral_params *ptr;   //  position within list
        static struct dihedral_params *tail;  //  Pointer to the end of 
                                              //  the list so we can add
                                              //  entries to the end of the
                                              //  list in constant time
        int i;                                      //  Loop counter
        int replace;                          //  replace values?

        // keep track of the last dihedral param read to avoid spurious
        // error messages.
        static struct dihedral_params last_dihedral; 

        /*  If the list is currently empty, then the new node is the list*/
        if (dihedralp == NULL)
        {
                dihedralp=new_node;
                tail=new_node;
                memcpy(&last_dihedral, new_node, sizeof(dihedral_params));

                return;
        }

        /*  The list isn't empty, so check for a duplicate                */
        ptr=dihedralp;

        while (ptr != NULL)
        {
                int same_as_last = 0;
                if (  ( (strcasecmp(new_node->atom1name, ptr->atom1name) == 0) &&
                       (strcasecmp(new_node->atom2name, ptr->atom2name) == 0) &&
                       (strcasecmp(new_node->atom3name, ptr->atom3name) == 0) &&
                       (strcasecmp(new_node->atom4name, ptr->atom4name) == 0) ) ||
                     ( (strcasecmp(new_node->atom4name, ptr->atom1name) == 0) &&
                       (strcasecmp(new_node->atom3name, ptr->atom2name) == 0) &&
                       (strcasecmp(new_node->atom2name, ptr->atom3name) == 0) &&
                       (strcasecmp(new_node->atom1name, ptr->atom4name) == 0) )
                       )
                {
                        /*  Found a duplicate                                */
                        
                        // check for same_as_last.  Note: don't believe the echoWarn crap; it controls
                        // not just whether we print warning messages, but whether we actually change
                        // values or not!  

                        if ( ( !strcmp(ptr->atom1name, last_dihedral.atom1name) && 
                               !strcmp(ptr->atom2name, last_dihedral.atom2name) &&
                               !strcmp(ptr->atom3name, last_dihedral.atom3name) &&
                               !strcmp(ptr->atom4name, last_dihedral.atom4name)))
                          same_as_last = 1;

                        //****** BEGIN CHARMM/XPLOR type changes
                        /* we do not care about identical replacement */
                        int echoWarn=1;  // echo warning messages ?

                        // ptr->multiplicity will always be >= new_node->multiplicity
                        for (i=0; i<ptr->multiplicity; i++)
                        {
                          if ((ptr->values[i].k == new_node->values[0].k) && 
                              (ptr->values[i].n == new_node->values[0].n) &&
                              (ptr->values[i].delta == new_node->values[0].delta)) 
                          {
                            // found an identical replacement
                            echoWarn=0; 
                            break;
                          }

                        }
                  
                        if (echoWarn)
                        {
                          if (!same_as_last) {
                            iout << "\n" << iWARN << "DUPLICATE DIHEDRAL ENTRY FOR "
                                 << ptr->atom1name << "-"
                                 << ptr->atom2name << "-"
                                 << ptr->atom3name << "-"
                                 << ptr->atom4name
                                 << "\nPREVIOUS VALUES MULTIPLICITY: " << ptr->multiplicity << "\n";
                          }
                          replace=0;
                          
                          for (i=0; i<ptr->multiplicity; i++)
                          {
                            if (!same_as_last) {
                              iout << "  k=" << ptr->values[i].k
                                   << "  n=" << ptr->values[i].n
                                   << "  delta=" << ptr->values[i].delta << "\n";
                            }
                            if (ptr->values[i].n == new_node->values[0].n)
                            {
                              iout << iWARN << "IDENTICAL PERIODICITY! REPLACING OLD VALUES BY: \n";
                              ptr->values[i].k = new_node->values[0].k;
                              ptr->values[i].delta = new_node->values[0].delta;
                              iout << "  k=" << ptr->values[i].k
                                   << "  n=" << ptr->values[i].n
                                   << "  delta=" << ptr->values[i].delta<< "\n";
                              replace=1;
                              break;
                            }
                          }

                          if (!replace)
                          {
                            ptr->multiplicity += 1;

                            if (ptr->multiplicity > MAX_MULTIPLICITY)
                            {
                              char err_msg[181];

                              sprintf(err_msg, "Multiple dihedral with multiplicity of %d greater than max of %d",
                                      ptr->multiplicity, MAX_MULTIPLICITY);
                              NAMD_die(err_msg);
                            }
                            if (!same_as_last) 
                              iout << "INCREASING MULTIPLICITY TO: " << ptr->multiplicity << "\n";

                            i= ptr->multiplicity - 1; 
                            ptr->values[i].k = new_node->values[0].k;
                            ptr->values[i].n = new_node->values[0].n;
                            ptr->values[i].delta = new_node->values[0].delta;

                            if (!same_as_last) 
                              iout << "  k=" << ptr->values[i].k
                                   << "  n=" << ptr->values[i].n
                                   << "  delta=" << ptr->values[i].delta<< "\n";
                          }
                        
                          iout << endi;
                        } 
                        //****** END CHARMM/XPLOR type changes

                        memcpy(&last_dihedral, new_node, sizeof(dihedral_params));
                        delete new_node;

                        return;
                }

                ptr=ptr->next;
        }

        /*  CHARMM and XPLOR wildcards for dihedrals are luckily the same */
        /*  Check to see if we have any wildcards.  Since specific        */
        /*  entries are to take precedence, we'll put anything without  */
        /*  wildcards at the begining of the list and anything with     */
        /*  wildcards at the end of the list.  Then, we can just do a   */
        /*  linear search for a bond and be guaranteed to have specific */
        /*  entries take precendence over over wildcards                */
        if ( new_node->atom1wild ||
             new_node->atom2wild ||
             new_node->atom3wild ||
             new_node->atom4wild )
        {
                /*  add to the end of the list                                */
                tail->next=new_node;
                tail=new_node;

                memcpy(&last_dihedral, new_node, sizeof(dihedral_params));
                return;
        }
        else
        {
                /*  add to the head of the list                                */
                new_node->next=dihedralp;
                dihedralp=new_node;

                memcpy(&last_dihedral, new_node, sizeof(dihedral_params));
                return;
        }

}
/*                END OF FUNCTION add_to_charmm_dihedral_list                */
//****** END CHARMM/XPLOR type changes

/************************************************************************/
/*                  */
/*      FUNCTION add_improper_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with improper parameters    */
/*                  */
/*  this function adds an improper parameter.  It parses up the     */
/*   input line and then adds it to the binary tree used to store the   */
/*   improper parameters.            */
/*                  */
/************************************************************************/

void Parameters::add_improper_param(char *buf, FILE *fd)

{
  char atom1name[11];       //  Atom 1 type
  char atom2name[11];       //  Atom 2 type
  char atom3name[11];       //  Atom 3 type
  char atom4name[11];       //  Atom 4 type
  Real forceconstant;       //  Force constant 
  int periodicity;       //  Periodicity
  Real phase_shift;       //  Phase shift
  int read_count;         //  Count from sscanf
  struct improper_params *new_node;  //  New node
  int multiplicity;       //  Multiplicity for bonds
  int i;           //  Loop counter
  char buffer[513];       //  Buffer for new line
  int ret_code;         //  Return code

  //****** BEGIN CHARMM/XPLOR type changes
  /*  Parse up the line with sscanf                                */
  if (paramType == paraXplor)
  {
    /* read XPLOR format */
    read_count=sscanf(buf, "%*s %s %s %s %s MULTIPLE= %d %f %d %f\n", 
       atom1name, atom2name, atom3name, atom4name, &multiplicity, 
       &forceconstant, &periodicity, &phase_shift);
  }
  else if (paramType == paraCharmm)
  {
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %s %s %s %f %d %f\n", 
       atom1name, atom2name, atom3name, atom4name,  
       &forceconstant, &periodicity, &phase_shift); 
    multiplicity=1;      
  }

  if ( (read_count != 4) && (read_count != 8) && (paramType == paraXplor) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD IMPROPER FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }
  else if ( (read_count != 7) && (paramType == paraCharmm) )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD IMPROPER FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }

  if ( (read_count == 4) && (paramType == paraXplor) )
  //****** END CHARMM/XPLOR type changes
  {
    read_count=sscanf(buf, "%*s %*s %*s %*s %*s %f %d %f\n", 
          &forceconstant, &periodicity, &phase_shift);

    /*  Check to make sure we got what we expected    */
    if (read_count != 3)
    {
      char err_msg[512];

      sprintf(err_msg, "BAD IMPROPER FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buf);
      NAMD_die(err_msg);
    }

    multiplicity = 1;
  }

  if (multiplicity > MAX_MULTIPLICITY)
  {
    char err_msg[181];

    sprintf(err_msg, "Multiple improper with multiplicity of %d greater than max of %d",
       multiplicity, MAX_MULTIPLICITY);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new improper_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_improper_param");
  }

  /*  Assign the values for this bond.  As with the dihedrals,    */
  /*  the atom order doesn't matter        */
  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);
  strcpy(new_node->atom3name, atom3name);
  strcpy(new_node->atom4name, atom4name);
  new_node->multiplicity = multiplicity;
  if (paramType == paraXplor && periodicity != 0) phase_shift *= -1;
  new_node->values[0].k = forceconstant;
  new_node->values[0].n = periodicity;
  new_node->values[0].delta = phase_shift;

  new_node->next = NULL;

  //  Check to see if this improper has multiple values
  if (multiplicity > 1)
  {
    //  Loop through and read the other values
    for (i=1; i<multiplicity; i++)
    {
      ret_code = NAMD_read_line(fd, buffer);

      //  Strip off comments at the end of the line
      if (ret_code == 0)
      {
        NAMD_remove_comment(buffer);
      }

      //  Skip blank lines
      while ( (ret_code == 0) && (NAMD_blank_string(buffer)) )
      {
        ret_code = NAMD_read_line(fd, buffer);
      }

      if (ret_code != 0)
      {
        NAMD_die("EOF encoutner in middle of multiple improper");
      }

      //  Get the values from the line
      read_count=sscanf(buffer, "%f %d %f\n", 
            &forceconstant, &periodicity, &phase_shift);

      if (read_count != 3)
      {
        char err_msg[512];

        sprintf(err_msg, "BAD MULTIPLE FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buffer);
        NAMD_die(err_msg);
      }

      if (paramType == paraXplor && periodicity != 0) phase_shift *= -1;
      new_node->values[i].k = forceconstant;
      new_node->values[i].n = periodicity;
      new_node->values[i].delta = phase_shift;
    }
  }

  /*  Add the paramter to the list        */
  add_to_improper_list(new_node);  // works for both XPLOR & CHARMM

  return;
}
/*      END OF FUNCTION add_improper_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_improper_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node that is to be added to imporper_list    */
/*                  */
/*  this function adds a new dihedral parameter to the linked list  */
/*   of improper parameters.  First, it checks for duplicates.  If a    */
/*   duplicate is found, a warning message is printed, the old values   */
/*   are replaced with the new values, and the new node is freed.  If   */
/*   Otherwise, the node is added to the list.  This list is arranged   */
/*   so that bods with wildcards are placed at the tail of the list.    */
/*   This will guarantee that if we just do a linear search, we will    */
/*   always find an exact match before a wildcard match.    */
/*                  */
/************************************************************************/

void Parameters::add_to_improper_list(struct improper_params *new_node)

{
  int i;              //  Loop counter
  static struct improper_params *ptr;   //  position within list
  static struct improper_params *tail;  //  Pointer to the end of 
                //  the list so we can add
                //  entries to the end of the
                //  list in constant time

  /*  If the list is currently empty, then the new node is the list*/
  if (improperp == NULL)
  {
    improperp=new_node;
    tail=new_node;

    return;
  }

  /*  The list isn't empty, so check for a duplicate    */
  ptr=improperp;

  while (ptr != NULL)
  {
    if ( ( (strcasecmp(new_node->atom1name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom4name, ptr->atom4name) == 0) ) ||
         ( (strcasecmp(new_node->atom4name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom1name, ptr->atom4name) == 0) ) )
    {
      /*  Found a duplicate        */
      //****** BEGIN CHARMM/XPLOR type changes
      /* we do not care about identical replacement */
      int echoWarn=0;  // echo warning messages ?

      if (ptr->multiplicity != new_node->multiplicity) {echoWarn=1;}
      
      if (!echoWarn)
      {
        for (i=0; i<ptr->multiplicity; i++)
        {
          if (ptr->values[i].k != new_node->values[i].k) {echoWarn=1; break;}
          if (ptr->values[i].n != new_node->values[i].n) {echoWarn=1; break;}
          if (ptr->values[i].delta != new_node->values[i].delta) {echoWarn=1; break;}
        }
      }

      if (echoWarn)
      {
        iout << "\n" << iWARN << "DUPLICATE IMPROPER DIHEDRAL ENTRY FOR "
          << ptr->atom1name << "-"
          << ptr->atom2name << "-"
          << ptr->atom3name << "-"
          << ptr->atom4name
          << "\nPREVIOUS VALUES MULTIPLICITY " << ptr->multiplicity << "\n";
        
        for (i=0; i<ptr->multiplicity; i++)
        {
          iout <<     "  k=" << ptr->values[i].k
                   << "  n=" << ptr->values[i].n
                   << "  delta=" << ptr->values[i].delta;
        }

        iout << "\n" << "USING VALUES MULTIPLICITY " << new_node->multiplicity << "\n";

        for (i=0; i<new_node->multiplicity; i++)
        {
          iout <<     "  k=" << new_node->values[i].k
                   << "  n=" << new_node->values[i].n
                   << "  delta=" << new_node->values[i].delta;
        }

        iout << endi;

        ptr->multiplicity = new_node->multiplicity;

        for (i=0; i<new_node->multiplicity; i++)
        {
          ptr->values[i].k = new_node->values[i].k;
          ptr->values[i].n = new_node->values[i].n;
          ptr->values[i].delta = new_node->values[i].delta;
        }
      }
      //****** END CHARMM/XPLOR type changes

      delete new_node;

      return;
    }

    ptr=ptr->next;
  }

  /*  Check to see if we have any wildcards.  Since specific  */
  /*  entries are to take precedence, we'll put anything without  */
  /*  wildcards at the begining of the list and anything with     */
  /*  wildcards at the end of the list.  Then, we can just do a   */
  /*  linear search for a bond and be guaranteed to have specific */
  /*  entries take precendence over over wildcards          */
  if ( (strcasecmp(new_node->atom1name, "X") == 0) ||
       (strcasecmp(new_node->atom2name, "X") == 0) ||
       (strcasecmp(new_node->atom3name, "X") == 0) ||
       (strcasecmp(new_node->atom4name, "X") == 0) )
  {
    /*  add to the end of the list        */
    tail->next=new_node;
    tail=new_node;

    return;
  }
  else
  {
    /*  add to the head of the list        */
    new_node->next=improperp;
    improperp=new_node;

    return;
  }
}
/*    END OF FUNCTION add_to_improper_list      */

/************************************************************************/
/*                  */
/*      FUNCTION add_crossterm_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line from paramter file with crossterm parameters    */
/*                  */
/*  this function adds an crossterm parameter.  It parses up the     */
/*   input line and then adds it to the binary tree used to store the   */
/*   crossterm parameters.            */
/*                  */
/************************************************************************/

void Parameters::add_crossterm_param(char *buf, FILE *fd)

{
  char atom1name[11];       //  Atom 1 type
  char atom2name[11];       //  Atom 2 type
  char atom3name[11];       //  Atom 3 type
  char atom4name[11];       //  Atom 4 type
  char atom5name[11];       //  Atom 1 type
  char atom6name[11];       //  Atom 2 type
  char atom7name[11];       //  Atom 3 type
  char atom8name[11];       //  Atom 4 type
  int dimension;
  int read_count;         //  Count from sscanf
  struct crossterm_params *new_node;  //  New node
  char buffer[513];       //  Buffer for new line
  int ret_code;         //  Return code

  /* read CHARMM format */
  read_count=sscanf(buf, "%s %s %s %s %s %s %s %s %d\n", 
     atom1name, atom2name, atom3name, atom4name,  
     atom5name, atom6name, atom7name, atom8name,  
     &dimension); 

  if ( (read_count != 9) || dimension < 1 || dimension > 1000 )
  {
    char err_msg[512];

    sprintf(err_msg, "BAD CMAP FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new crossterm_params(dimension);

  /*  Assign the values for this bond.  As with the dihedrals,    */
  /*  the atom order doesn't matter        */
  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);
  strcpy(new_node->atom3name, atom3name);
  strcpy(new_node->atom4name, atom4name);
  strcpy(new_node->atom5name, atom5name);
  strcpy(new_node->atom6name, atom6name);
  strcpy(new_node->atom7name, atom7name);
  strcpy(new_node->atom8name, atom8name);

  new_node->next = NULL;

  int nterms = dimension * dimension;
  int nread = 0;

  //  Loop through and read the other values
  while ( nread < nterms ) {
    ret_code = NAMD_read_line(fd, buffer);

    //  Strip off comments at the end of the line
    if (ret_code == 0) {
      NAMD_remove_comment(buffer);
    }

    //  Skip blank lines
    while ( (ret_code == 0) && (NAMD_blank_string(buffer)) ) {
      ret_code = NAMD_read_line(fd, buffer);
      if (ret_code == 0) {
        NAMD_remove_comment(buffer);
      }
    }

    if (ret_code != 0) {
      NAMD_die("EOF encoutner in middle of CMAP");
    }

    //  Get the values from the line
    read_count=sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
	new_node->values + nread,
	new_node->values + nread+1,
	new_node->values + nread+2,
	new_node->values + nread+3,
	new_node->values + nread+4,
	new_node->values + nread+5,
	new_node->values + nread+6,
	new_node->values + nread+7);

    nread += read_count;

    if (read_count == 0 || nread > nterms) {
      char err_msg[512];

      sprintf(err_msg, "BAD CMAP FORMAT IN PARAMETER FILE\nLINE=*%s*\n", buffer);
      NAMD_die(err_msg);
    }
  }

  /*  Add the paramter to the list        */
  add_to_crossterm_list(new_node);

  return;
}
/*      END OF FUNCTION add_crossterm_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_crossterm_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node that is to be added to crossterm_list    */
/*                  */
/*  this function adds a new crossterm parameter to the linked list  */
/*   of crossterm parameters.  First, it checks for duplicates.  If a    */
/*   duplicate is found, a warning message is printed, the old values   */
/*   are replaced with the new values, and the new node is freed.  If   */
/*   Otherwise, the node is added to the list.  This list is arranged   */
/*   so that bods with wildcards are placed at the tail of the list.    */
/*   This will guarantee that if we just do a linear search, we will    */
/*   always find an exact match before a wildcard match.    */
/*                  */
/************************************************************************/

void Parameters::add_to_crossterm_list(struct crossterm_params *new_node)

{
  int i;              //  Loop counter
  static struct crossterm_params *ptr;   //  position within list
  static struct crossterm_params *tail;  //  Pointer to the end of 
                //  the list so we can add
                //  entries to the end of the
                //  list in constant time

  /*  If the list is currently empty, then the new node is the list*/
  if (crosstermp == NULL)
  {
    crosstermp=new_node;
    tail=new_node;

    return;
  }

  /*  The list isn't empty, so check for a duplicate    */
  ptr=crosstermp;

  while (ptr != NULL)
  {
    if ( ( (strcasecmp(new_node->atom1name, ptr->atom1name) == 0) &&
           (strcasecmp(new_node->atom2name, ptr->atom2name) == 0) &&
           (strcasecmp(new_node->atom3name, ptr->atom3name) == 0) &&
           (strcasecmp(new_node->atom4name, ptr->atom4name) == 0) &&
           (strcasecmp(new_node->atom5name, ptr->atom5name) == 0) &&
           (strcasecmp(new_node->atom6name, ptr->atom6name) == 0) &&
           (strcasecmp(new_node->atom7name, ptr->atom7name) == 0) &&
           (strcasecmp(new_node->atom8name, ptr->atom8name) == 0) ) )
    {
      /*  Found a duplicate        */
      /* we do not care about identical replacement */
      int echoWarn=0;  // echo warning messages ?

      if (ptr->dimension != new_node->dimension) {echoWarn=1;}
      
      if (!echoWarn)
      {
        int nvals = ptr->dimension * ptr->dimension;
        for (i=0; i<nvals; i++)
        {
          if (ptr->values[i] != new_node->values[i]) {echoWarn=1; break;}
        }
      }

      if (echoWarn)
      {
        iout << "\n" << iWARN << "DUPLICATE CMAP ENTRY FOR "
          << ptr->atom1name << "-"
          << ptr->atom2name << "-"
          << ptr->atom3name << "-"
          << ptr->atom4name << " "
          << ptr->atom5name << "-"
          << ptr->atom6name << "-"
          << ptr->atom7name << "-"
          << ptr->atom8name << ", USING NEW VALUES\n";
        
        iout << endi;

        ptr->dimension = new_node->dimension;

        BigReal *tmpvalues = ptr->values;
        ptr->values = new_node->values;
        new_node->values = tmpvalues;
      }

      delete new_node;

      return;
    }

    ptr=ptr->next;
  }

  /*  Check to see if we have any wildcards.  Since specific  */
  /*  entries are to take precedence, we'll put anything without  */
  /*  wildcards at the begining of the list and anything with     */
  /*  wildcards at the end of the list.  Then, we can just do a   */
  /*  linear search for a bond and be guaranteed to have specific */
  /*  entries take precendence over over wildcards          */
  if ( (strcasecmp(new_node->atom1name, "X") == 0) ||
       (strcasecmp(new_node->atom2name, "X") == 0) ||
       (strcasecmp(new_node->atom3name, "X") == 0) ||
       (strcasecmp(new_node->atom4name, "X") == 0) ||
       (strcasecmp(new_node->atom5name, "X") == 0) ||
       (strcasecmp(new_node->atom6name, "X") == 0) ||
       (strcasecmp(new_node->atom7name, "X") == 0) ||
       (strcasecmp(new_node->atom8name, "X") == 0) )
  {
    /*  add to the end of the list        */
    tail->next=new_node;
    tail=new_node;

    return;
  }
  else
  {
    /*  add to the head of the list        */
    new_node->next=crosstermp;
    crosstermp=new_node;

    return;
  }
}
/*    END OF FUNCTION add_to_crossterm_list      */

/************************************************************************/
/*                  */
/*      FUNCTION add_vdw_param        */
/*                  */
/*  INPUTS:                */
/*  buf - line containing the vdw information      */
/*                  */
/*  add_vdw_param adds a vdw parameter for an atom to the current   */
/*  binary tree of values.            */
/*                  */
/************************************************************************/

void Parameters::add_vdw_param(char *buf)

{
  char atomname[11];    //  atom type of paramter
  Real sigma;      //  sigma value for this atom
  Real epsilon;      //  epsilon value for this atom
  Real sigma14;      //  sigma value for 1-4 interactions
  Real epsilon14;      //  epsilon value for 1-4 interactions
  Real sqrt26;         //  2^(1/6)
  int read_count;      //  count returned by sscanf
  struct vdw_params *new_node;  //  new node for tree

  //****** BEGIN CHARMM/XPLOR type changes
  /*  Parse up the line with sscanf        */
  if (paramType == paraXplor)
  {
    /* read XPLOR format */
    read_count=sscanf(buf, "%*s %s %f %f %f %f\n", atomname, 
       &epsilon, &sigma, &epsilon14, &sigma14);
  }
  else if (paramType == paraCharmm)
  {
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %*f %f %f %*f %f %f\n", atomname, 
       &epsilon, &sigma, &epsilon14, &sigma14);
  }

  /*  Check to make sure we got what we expected      */
  if ((read_count != 5) && (paramType == paraXplor))
  {
    char err_msg[512];

    sprintf(err_msg, "BAD vdW FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }
  else if ( ((read_count != 5) && (read_count != 3)) && (paramType == paraCharmm))
  {
    char err_msg[512];

    sprintf(err_msg, "BAD vdW FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*\n", buf);
    NAMD_die(err_msg);
  }

  if (paramType == paraCharmm)
  {
    // convert CHARMM to XPLOR format
    epsilon*=-1.;
    sqrt26=pow(2.,(1./6.));
    sigma=2.*sigma/sqrt26; 

    if (read_count == 3)
    {
      epsilon14=epsilon;
      sigma14=sigma;
    }
    else
    {
      epsilon14*=-1.;
      sigma14=2.*sigma14/sqrt26; 
    }
  }
  //****** END CHARMM/XPLOR type changes

  if ( epsilon < 0. || epsilon14 < 0. ) {
    iout << iWARN << "Ignoring VDW parameter with negative epsilon:\n"
        << buf << "\n" << endi;
    return;
  }

  /*  Allocate a new node            */
  new_node = new vdw_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_vdw_param");
  }

  /*  Assign the values to the new node        */
  strcpy(new_node->atomname, atomname);
  new_node->sigma = sigma;
  new_node->sigma14 = sigma14;
  new_node->epsilon = epsilon;
  new_node->epsilon14 = epsilon14;

  new_node->left = NULL;
  new_node->right = NULL;

  /*  Add the new node into the tree        */
  vdwp=add_to_vdw_tree(new_node, vdwp);

  return;
}
/*      END OF FUNCTION add_vdw_param      */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_vdw_tree      */
/*                  */
/*   INPUTS:                */
/*  new_node - node to add to tree          */
/*  tree - tree to add the node to          */
/*                  */
/*   OUTPUTS:                */
/*  the function returns a pointer to the tree with the node added  */
/*                  */
/*  this function adds a vdw to the binary tree containing the      */
/*   parameters.              */
/*                  */
/************************************************************************/

struct vdw_params *Parameters::add_to_vdw_tree(struct vdw_params *new_node,
             struct vdw_params *tree)

{
  int compare_code;  //  Return code from strcasecmp

  /*  If the tree is currently empty, the new node is the tree    */
  if (tree == NULL)
    return(new_node);

  compare_code = strcasecmp(new_node->atomname, tree->atomname);

  /*  Check to see if we have a duplicate        */
  if (compare_code==0)
  {
    /*  We have a duplicate.  So print out a warning   */
    /*  message, copy the new values into the current node  */
    /*  of the tree, and then free the new_node    */
    if ((tree->sigma != new_node->sigma) || 
        (tree->epsilon != new_node->epsilon) ||
        (tree->sigma14 != new_node->sigma14) ||
        (tree->epsilon14 != new_node->epsilon14))
    {
      iout << iWARN << "DUPLICATE vdW ENTRY FOR " << tree->atomname
        << "\nPREVIOUS VALUES  sigma=" << tree->sigma
        << " epsilon=" << tree->epsilon
        << " sigma14=" << tree->sigma14
        << " epsilon14=" << tree->epsilon14
        << "\n   USING VALUES  sigma=" << new_node->sigma
        << " epsilon=" << new_node->epsilon
        << " sigma14=" << new_node->sigma14
        << " epsilon14=" << new_node->epsilon14
        << "\n" << endi;

      tree->sigma=new_node->sigma;
      tree->epsilon=new_node->epsilon;
      tree->sigma14=new_node->sigma14;
      tree->epsilon14=new_node->epsilon14;
    }

    delete new_node;

    return(tree);
  }

  /*  Otherwise, if the new node is less than the head of    */
  /*  the tree, add it to the left child, and if it is greater  */
  /*  add it to the right child          */
  if (compare_code < 0)
  {
    tree->left = add_to_vdw_tree(new_node, tree->left);
  }
  else
  {
    tree->right = add_to_vdw_tree(new_node, tree->right);
  }

  return(tree);
}
/*      END OF FUNCTION add_to_vdw_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION add_table_pair_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line containing the table_pair information      */
/*                  */
/*  this function adds a table_pair parameter to the current          */
/*   parameters.              */
/*                  */
/************************************************************************/

void Parameters::add_table_pair_param(char *buf)

{
  char atom1name[11];      //  Atom 1 name
  char atom2name[11];      //  Atom 2 name
  char tabletype[11];      // Name of pair interaction
  int K;           // Table entry type for pair
  int read_count;        //  count from sscanf
  struct table_pair_params *new_node;  //  new node

  /*  Parse up the input line using sscanf      */
  read_count=sscanf(buf, "%s %s %s\n", atom1name, 
  atom2name, tabletype);

  /*  Check to make sure we got what we expected      */
  if ((read_count != 3))
  {
    char err_msg[512];

    sprintf(err_msg, "BAD TABLE PAIR FORMAT IN PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }

  /*  Allocate a new node            */
  new_node = new table_pair_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_table_pair_param\n");
  }

  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);

  /* Find the proper table type for this pairing */
  K = get_int_table_type(tabletype);
//  printf("Looking for type %s; got %i\n", tabletype, K);
  if (K < 0) {
    char err_msg[512];
    sprintf(err_msg, "Couldn't find table parameters for table interaction %s!\n", tabletype);
    NAMD_die(err_msg);
  }

  /*  Assign values to this node          */
  new_node->K = K;

  new_node->next = NULL;

  /*  Add this node to the tree          */
//  printf("Adding pair parameter with index %i\n", K);
  add_to_table_pair_list(new_node);

  return;
}
/*      END OF FUNCTION add_vdw_par_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_vdw_pair_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line containing the vdw_pair information      */
/*                  */
/*  this function adds a vdw_pair parameter to the current          */
/*   parameters.              */
/*                  */
/************************************************************************/

void Parameters::add_vdw_pair_param(char *buf)

{
  char atom1name[11];      //  Atom 1 name
  char atom2name[11];      //  Atom 2 name
  Real A;          //  A value for pair
  Real B;          //  B value for pair
  Real A14;        //  A value for 1-4 ints
  Real B14;        //  B value for 1-4 ints
  int read_count;        //  count from sscanf
  struct vdw_pair_params *new_node;  //  new node

  /*  Parse up the input line using sscanf      */
  if (paramType == paraXplor)
  {
    /* read XPLOR format */
    read_count=sscanf(buf, "%*s %s %s %f %f %f %f\n", atom1name, 
       atom2name, &A, &B, &A14, &B14);
  }
  else if (paramType == paraCharmm)
  {
    Real well, rmin, well14, rmin14;
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %s %f %f %f %f\n", atom1name, 
       atom2name, &well, &rmin, &well14, &rmin14);
    if ( read_count == 4 ) { well14 = well; rmin14 = rmin; }
    A = -1. * well * pow(rmin, 12.);
    B = -2. * well * pow(rmin, 6.);
    A14 = -1. * well14 * pow(rmin14, 12.);
    B14 = -2. * well14 * pow(rmin14, 6.);
  }

  /*  Check to make sure we got what we expected      */
  if ((read_count != 6) && (paramType == paraXplor))
  {
    char err_msg[512];

    sprintf(err_msg, "BAD vdW PAIR FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }
  if ((read_count != 4) && (read_count != 6) && (paramType == paraCharmm))
  {
    char err_msg[512];

    sprintf(err_msg, "BAD vdW PAIR FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }


  /*  Allocate a new node            */
  new_node = new vdw_pair_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::add_vdw_pair_param\n");
  }

  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);

  /*  Assign values to this node          */
  new_node->A = A;
  new_node->A14 = A14;
  new_node->B = B;
  new_node->B14 = B14;

  new_node->next = NULL;

  /*  Add this node to the tree          */
  add_to_vdw_pair_list(new_node);

  return;
}
/*      END OF FUNCTION add_vdw_par_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_nbthole_pair_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line containing the nbthole_pair information      */
/*                  */
/*  this function adds a nbthole_pair parameter to the current          */
/*   parameters.              */
/*                  */
/************************************************************************/

void Parameters::add_nbthole_pair_param(char *buf)

{
  char atom1name[11];      //  Atom 1 name
  char atom2name[11];      //  Atom 2 name
  Real tholeij;            //  nonbonded thole pair thole
  int read_count;        //  count from sscanf
  struct nbthole_pair_params *new_node;  //  new node

  /*  Parse up the input line using sscanf      */
  if (paramType == paraCharmm)
  {
    /* read CHARMM format */
    read_count=sscanf(buf, "%s %s %f\n", atom1name,
       atom2name, &tholeij);
  }

  /*  Check to make sure we got what we expected      */
  if ((read_count != 3) && (paramType == paraCharmm))
  {
    char err_msg[512];

    sprintf(err_msg, "BAD NBTHOLE PAIR FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*", buf);
    NAMD_die(err_msg);
  }


  /*  Allocate a new node            */
  new_node = new nbthole_pair_params;

  if (new_node == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::nbthole_pair_param\n");
  }

  strcpy(new_node->atom1name, atom1name);
  strcpy(new_node->atom2name, atom2name);

  /*  Assign values to this node          */
  new_node->tholeij = tholeij;

  new_node->next = NULL;

  /*  Add this node to the tree          */
  add_to_nbthole_pair_list(new_node);

  return;
}
/*      END OF FUNCTION add_nbthole_par_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_hb_pair_param      */
/*                  */
/*   INPUTS:                */
/*  buf - line containing the hydrogen bond information    */
/*                  */
/*  this function adds data for a hydrogen bond interaction pair    */
/*   to the hbondParams object.                                         */
/*                  */
/************************************************************************/

void Parameters::add_hb_pair_param(char *buf)

{
#if 0
  char a1n[11];      //  Atom 1 name
  char a2n[11];      //  Atom 2 name
  Real A, B;      //  A, B value for pair

  //****** BEGIN CHARMM/XPLOR type changes
  //// luckily the format and units are the same CHARMM is just missing the HBON marker
  /*  Parse up the input line using sscanf      */
  if (paramType == paraXplor) {
    if (sscanf(buf, "%*s %s %s %f %f\n", a1n, a2n, &A, &B) != 4) {
      char err_msg[512];
      sprintf(err_msg, "BAD HBOND PAIR FORMAT IN XPLOR PARAMETER FILE\nLINE=*%s*", buf);
      NAMD_die(err_msg);
    }
  }
  else if (paramType == paraCharmm) {
    if (sscanf(buf, "%s %s %f %f\n", a1n, a2n, &A, &B) != 4) {
      char err_msg[512];
      sprintf(err_msg, "BAD HBOND PAIR FORMAT IN CHARMM PARAMETER FILE\nLINE=*%s*", buf);
      NAMD_die(err_msg);
    }
  }
  //****** END CHARMM/XPLOR type changes

  /*  add data */
  if (hbondParams.add_hbond_pair(a1n, a2n, A, B) == FALSE) {
    iout << "\n" << iWARN << "Duplicate HBOND parameters for types " << a1n
    << " and " << a2n << " found; using latest values." << "\n" << endi;
  }
#endif
}
/*      END OF FUNCTION add_hb_par_param    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_table_pair_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node to be added to list        */
/*                  */
/*  This function adds a link to the end of the table_pair_list list  */
/*                  */
/************************************************************************/

void Parameters::add_to_table_pair_list(struct table_pair_params *new_node)

{
     static struct table_pair_params *tail=NULL;
  struct table_pair_params *ptr;
  int compare_code;
  

  //  If the list was empty, then just make the new node the list
  if (table_pairp == NULL)
  {
     table_pairp = new_node;
     tail = new_node;
     return;
  }
  
  ptr = table_pairp;

  //  Now check the list to see if we have a duplicate entry
  while (ptr!=NULL)
  {
      /*  Compare atom 1            */
      compare_code = strncasecmp(new_node->atom1name, ptr->atom1name, 5);
      
      if (compare_code == 0)
      {
    /*  Atom 1 is the same, compare atom 2      */
    compare_code = strcasecmp(new_node->atom2name, ptr->atom2name);

    if (compare_code==0)
    {
      /*  Found a duplicate.  Print out a warning   */
      /*  message, assign the values to the current   */
      /*  node in the tree, and then free the new_node*/
      iout << iWARN << "DUPLICATE TABLE PAIR ENTRY FOR "
        << new_node->atom1name << "-"
        << new_node->atom2name
        << "\n" << endi;

        ptr->K=new_node->K;

      delete new_node;

      return;
    }
      }
      
      ptr = ptr->next;
  }

  //  We didn't find a duplicate, so add this node to the end
  //  of the list
  tail->next = new_node;
  tail = new_node;
}
/*      END OF FUNCTION add_to_vdw_pair_list    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_vdw_pair_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node to be added to list        */
/*                  */
/*  This function adds a link to the end of the vdw_pair_list list  */
/*                  */
/************************************************************************/

void Parameters::add_to_vdw_pair_list(struct vdw_pair_params *new_node)

{
     static struct vdw_pair_params *tail=NULL;
  struct vdw_pair_params *ptr;
  int compare_code;
  

  //  If the list was empty, then just make the new node the list
  if (vdw_pairp == NULL)
  {
     vdw_pairp = new_node;
     tail = new_node;
     return;
  }
  
  ptr = vdw_pairp;

  //  Now check the list to see if we have a duplicate entry
  while (ptr!=NULL)
  {
      /*  Compare atom 1            */
      compare_code = strcasecmp(new_node->atom1name, ptr->atom1name);
      
      if (compare_code == 0)
      {
    /*  Atom 1 is the same, compare atom 2      */
    compare_code = strcasecmp(new_node->atom2name, ptr->atom2name);

    if (compare_code==0)
    {
      /*  Found a duplicate.  Print out a warning   */
      /*  message, assign the values to the current   */
      /*  node in the tree, and then free the new_node*/
      if ((ptr->A != new_node->A) ||
          (ptr->B != new_node->B) ||
          (ptr->A14 != new_node->A14) ||
          (ptr->B14 != new_node->B14))
      {
        iout << iWARN << "DUPLICATE vdW PAIR ENTRY FOR "
          << new_node->atom1name << "-"
          << new_node->atom2name
          << "\nPREVIOUS VALUES  A=" << ptr->A
          << " B=" << ptr->B
          << " A14=" << ptr->A14
          << " B14" << ptr->B14
          << "\n   USING VALUES  A=" << new_node->A
          << " B=" << new_node->B
          << " A14=" << new_node->A14
          << " B14" << new_node->B14
          << "\n" << endi;

        ptr->A=new_node->A;
        ptr->B=new_node->B;
        ptr->A14=new_node->A14;
        ptr->B14=new_node->B14;
      }

      delete new_node;

      return;
    }
      }
      
      ptr = ptr->next;
  }

  //  We didn't find a duplicate, so add this node to the end
  //  of the list
  tail->next = new_node;
  tail = new_node;
}
/*      END OF FUNCTION add_to_vdw_pair_list    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_nbthole_pair_list      */
/*                  */
/*   INPUTS:                */
/*  new_node - node to be added to list        */
/*                  */
/*  This function adds a link to the end of the nbthole_pair_list list  */
/*                  */
/************************************************************************/

void Parameters::add_to_nbthole_pair_list(struct nbthole_pair_params *new_node)

{
     static struct nbthole_pair_params *tail=NULL;
  struct nbthole_pair_params *ptr;
  int compare_code;


  //  If the list was empty, then just make the new node the list
  if (nbthole_pairp == NULL)
  {
     nbthole_pairp = new_node;
     tail = new_node;
     return;
  }

  ptr = nbthole_pairp;

  tail->next = new_node;
  tail = new_node;
}
/*      END OF FUNCTION add_to_nbthole_pair_list    */

/************************************************************************/
/*                  */
/*      FUNCTION done_reading_files      */
/*                  */
/*  This function is used to signal the Parameters object that all  */
/*  of the parameter files have been read.  Once the object knows this, */
/*  it can set un indexes for all the parameters and transfer the values*/
/*  to linear arrays.  This will allow constant time access from this   */
/*  point on.                */
/*                  */
/************************************************************************/

void Parameters::done_reading_files()

{
  AllFilesRead = TRUE;

  //  Allocate space for all of the arrays
  if (NumBondParams)
  {
    bond_array = new BondValue[NumBondParams];

    if (bond_array == NULL)
    {
      NAMD_die("memory allocation of bond_array failed!");
    }
  }

  if (NumAngleParams)
  {
    angle_array = new AngleValue[NumAngleParams];

    if (angle_array == NULL)
    {
      NAMD_die("memory allocation of angle_array failed!");
    }
  }

  if (NumDihedralParams)
  {
    dihedral_array = new DihedralValue[NumDihedralParams];

    if (dihedral_array == NULL)
    {
      NAMD_die("memory allocation of dihedral_array failed!");
    }
    memset(dihedral_array, 0, NumDihedralParams*sizeof(DihedralValue));
  }

  if (NumImproperParams)
  {
    improper_array = new ImproperValue[NumImproperParams];

    if (improper_array == NULL)
    {
      NAMD_die("memory allocation of improper_array failed!");
    }
    memset(improper_array, 0, NumImproperParams*sizeof(ImproperValue));
  }

  if (NumCrosstermParams)
  {
    crossterm_array = new CrosstermValue[NumCrosstermParams];
    memset(crossterm_array, 0, NumCrosstermParams*sizeof(CrosstermValue));
  }

  // JLai
  if (NumGromacsPairParams)
  {
    gromacsPair_array = new GromacsPairValue[NumGromacsPairParams];
    memset(gromacsPair_array, 0, NumGromacsPairParams*sizeof(GromacsPairValue)); 
  }
  // End of JLai

  if (NumVdwParams)
  {
          atomTypeNames = new char[NumVdwParams*(MAX_ATOMTYPE_CHARS+1)];
    vdw_array = new VdwValue[NumVdwParams];
    
    if (vdw_array == NULL)
    {
      NAMD_die("memory allocation of vdw_array failed!");
    }
  }
  if (NumNbtholePairParams)
  {
    nbthole_array = new NbtholePairValue[NumNbtholePairParams];

    if(nbthole_array == NULL)
    {
      NAMD_die("memory allocation of nbthole_array failed!");
    }
  }
  //  Assign indexes to each of the parameters and populate the
  //  arrays using the binary trees and linked lists that we have
  //  already read in
  index_bonds(bondp, 0);
  index_angles(anglep, 0);
  NumVdwParamsAssigned = index_vdw(vdwp, 0);
  index_dihedrals();
  index_impropers();
  index_crossterms();
  convert_nbthole_pairs();
  //  Convert the vdw pairs
  convert_vdw_pairs();
  convert_table_pairs();
}
/*      END OF FUNCTION done_reading_files    */

/************************************************************************/
/*                  */
/*      FUNCTION index_bonds        */
/*                  */
/*   INPUTS:                */
/*  tree - The tree that is to be indexed        */
/*  index - index to start with          */
/*                  */
/*  This is a recursive routine that will traverse the binary tree  */
/*   of bond parameters, assigning an index to each one, and copying    */
/*   the data from the binary tree to the array that will be used from  */
/*   here on.                */
/*                  */
/************************************************************************/

Index Parameters::index_bonds(struct bond_params *tree, Index index)

{
  //  Tree is empty, do nothing
  if (tree==NULL)
    return(index);

  //  If I have a left subtree, index it first
  if (tree->left != NULL)
  {
    index=index_bonds(tree->left, index);
  }

  //  Now assign an index to top node and populate array
  tree->index = index;
  bond_array[index].k = tree->forceconstant;
  bond_array[index].x0 = tree->distance;
  index++;

  //  If I have a right subtree, index it
  if (tree->right != NULL)
  {
    index=index_bonds(tree->right, index);
  }

  return(index);
}
/*      END OF FUNCTION index_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION index_angles        */
/*                  */
/*   INPUTS:                */
/*  tree - The tree that is to be indexed        */
/*  index - index to start with          */
/*                  */
/*  This is a recursive routine that will traverse the binary tree  */
/*   of angle parameters, assigning an index to each one, and copying   */
/*   the data from the binary tree to the array that will be used from  */
/*   here on.                */
/*                  */
/************************************************************************/

Index Parameters::index_angles(struct angle_params *tree, Index index)

{
  //  Tree is empty, do nothing
  if (tree==NULL)
    return(index);

  //  If I have a left subtree, index it first
  if (tree->left != NULL)
  {
    index=index_angles(tree->left, index);
  }

  //  Now assign an index to top node and populate array
  tree->index = index;

  angle_array[index].k = tree->forceconstant;
  angle_array[index].k_ub = tree->k_ub;
  angle_array[index].r_ub = tree->r_ub;
  angle_array[index].normal = tree->normal;

  //  Convert the angle to radians before storing it
  angle_array[index].theta0 = (tree->angle*PI)/180.0;
  index++;

  //  If I have a right subtree, index it
  if (tree->right != NULL)
  {
    index=index_angles(tree->right, index);
  }

  return(index);
}
/*      END OF FUNCTION index_angles      */

/************************************************************************/
/*                  */
/*      FUNCTION index_dihedrals      */
/*                  */
/*  This function walks down the linked list of dihedral parameters */
/*  and assigns an index to each one.  It also copies the data from this*/
/*  linked list to the arrays that will be used from here on out  */
/*                  */
/************************************************************************/

void Parameters::index_dihedrals()

{
  struct dihedral_params *ptr;  //  Current location in list
  Index index=0;      //  Current index value
  int i;        //  Loop counter

  //  Allocate an array to hold the multiplicity present in the
  //  parameter file for each bond.  This will be used to check
  //  the multiplicities that are detected in the psf file

  //  This is kind of ugly, but necessary because of the way that
  //  X-PLOR psf files deal with Charmm22 parameters.  The way
  //  that multiple periodicities are specified is by having
  //  the bonds appear multiple times in the psf file.  This even
  //  if a bond type has multiple parameters defined, they
  //  will be used if the bond appears multiple times in the
  //  psf file.  So we need to store the number of parameters
  //  we have to make sure the psf file doesn't ask for more
  //  parameters than we really have, and we also need to track
  //  how many times the bond appears in the psf file so that
  //  we can decide how many parameters to actually use.
  //  This is different for CHARMM parameter files as stated below!
  maxDihedralMults = new int[NumDihedralParams];

  if (maxDihedralMults == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::index_dihedrals()");
  }
  
  //  Start at the head
  ptr = dihedralp;

  while (ptr != NULL)
  {
    //  Copy data to array and assign index

    //  Save the multiplicity in another array
    maxDihedralMults[index] = ptr->multiplicity;


    //****** BEGIN CHARMM/XPLOR type changes
    if (paramType == paraXplor)
    {
      //  Assign the multiplicity in the actual structure a bogus value
      //  that we will update in assign_dihedral_index
      dihedral_array[index].multiplicity = -1;
    }
    else if (paramType == paraCharmm)
    {
      // In a CHARMM psf file each dihedral will be only listed once
      // even if it has multiple terms. There is no point in comparing
      // to the psf information
      dihedral_array[index].multiplicity = ptr->multiplicity;
    } 
    //****** END CHARMM/XPLOR type changes

    for (i=0; i<ptr->multiplicity; i++)
    {
      dihedral_array[index].values[i].k = ptr->values[i].k;
      dihedral_array[index].values[i].n = ptr->values[i].n;

      //  Convert the angle to radians before storing it
      dihedral_array[index].values[i].delta = ptr->values[i].delta*PI/180.0;
    }

    ptr->index = index;

    index++;
    ptr=ptr->next;
  }
}
/*      END OF FUNCTION index_dihedrals      */

/************************************************************************/
/*                  */
/*      FUNCTION index_impropers      */
/*                  */
/*  This function walks down the linked list of improper parameters */
/*  and assigns an index to each one.  It also copies the data from this*/
/*  linked list to the arrays that will be used from here on out  */
/*                  */
/************************************************************************/

void Parameters::index_impropers()

{
  struct improper_params *ptr;  //  Current place in list
  Index index=0;      //  Current index value
  int i;        //  Loop counter

  //  Allocate an array to hold the multiplicity present in the
  //  parameter file for each bond.  This will be used to check
  //  the multiplicities that are detected in the psf file

  //  This is kind of ugly, but necessary because of the way that
  //  X-PLOR psf files deal with Charmm22 parameters.  The way
  //  that multiple periodicities are specified is by having
  //  the bonds appear multiple times in the psf file.  This even
  //  if a bond type has multiple parameters defined, they
  //  will be used if the bond appears multiple times in the
  //  psf file.  So we need to store the number of parameters
  //  we have to make sure the psf file doesn't ask for more
  //  parameters than we really have, and we also need to track
  //  how many times the bond appears in the psf file so that
  //  we can decide how many parameters to actually use.
  maxImproperMults = new int[NumImproperParams];

  if (maxImproperMults == NULL)
  {
    NAMD_die("memory allocation failed in Parameters::index_impropers()");
  }
  
  //  Start at the head
  ptr = improperp;

  while (ptr != NULL)
  {
    //  Copy data to array and assign index

    //  Save the multiplicity in another array
    maxImproperMults[index] = ptr->multiplicity;

    //  Assign the multiplicity in the actual structure a bogus value
    //  that we will update in assign_dihedral_index
    improper_array[index].multiplicity = -1;

    for (i=0; i<ptr->multiplicity; i++)
    {
      improper_array[index].values[i].k = ptr->values[i].k;
      improper_array[index].values[i].n = ptr->values[i].n;

      //  Convert the angle to radians before storing it
      improper_array[index].values[i].delta = ptr->values[i].delta*PI/180.0;
    }

    ptr->index=index;

    index++;
    ptr=ptr->next;
  }
}
/*      END OF FUNCTION index_impropers      */


/************************************************************************/
/*                  */
/*      FUNCTION index_crossterms      */
/*                  */
/*  This function walks down the linked list of crossterm parameters */
/*  and assigns an index to each one.  It also copies the data from this*/
/*  linked list to the arrays that will be used from here on out  */
/*                  */
/************************************************************************/

void crossterm_setup(CrosstermData *);

void Parameters::index_crossterms()

{
  struct crossterm_params *ptr;  //  Current place in list
  Index index=0;      //  Current index value
  int i,j,k;        //  Loop counter

  //  Start at the head
  ptr = crosstermp;

  while (ptr != NULL)
  {
    //  Copy data to array and assign index

    int N = CrosstermValue::dim - 1;

    if ( ptr->dimension != N ) {
      NAMD_die("Sorry, only CMAP dimension of 24 is supported");
    }

    k = 0;
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        crossterm_array[index].c[i][j].d00 = ptr->values[k];
        ++k;
      }
    }
    for (i=0; i<N; i++) {
        crossterm_array[index].c[i][N].d00 = 
				crossterm_array[index].c[i][0].d00;
        crossterm_array[index].c[N][i].d00 = 
				crossterm_array[index].c[0][i].d00;
    }
    crossterm_array[index].c[N][N].d00 = 
				crossterm_array[index].c[0][0].d00;

    crossterm_setup(&crossterm_array[index].c[0][0]);

    ptr->index=index;

    index++;
    ptr=ptr->next;
  }
}
/*      END OF FUNCTION index_crossterms      */

/************************************************************************/
/*                  */
/*      FUNCTION index_vdw        */
/*                  */
/*   INPUTS:                */
/*  tree - The tree that is to be indexed        */
/*  index - index to start with          */
/*                  */
/*  This is a recursive routine that will traverse the binary tree  */
/*   of vdw parameters, assigning an index to each one, and copying     */
/*   the data from the binary tree to the array that will be used from  */
/*   here on.                */
/*                  */
/************************************************************************/

Index Parameters::index_vdw(struct vdw_params *tree, Index index)

{
  //  If the tree is empty, do nothing
  if (tree==NULL)
    return(index);

  //  If I have a left subtree, populate it first
  if (tree->left != NULL)
  {
    index=index_vdw(tree->left, index);
  }

  //  Assign the index and copy the data to the array
  tree->index = index;

  vdw_array[index].sigma = tree->sigma;
  vdw_array[index].epsilon = tree->epsilon;
  vdw_array[index].sigma14 = tree->sigma14;
  vdw_array[index].epsilon14 = tree->epsilon14;

  char *nameloc = atom_type_name(index);
  strncpy(nameloc, tree->atomname, MAX_ATOMTYPE_CHARS);
  nameloc[MAX_ATOMTYPE_CHARS] = '\0';

//  iout << iWARN << "Parameters: Stored name for type " << index << ": '";
//      iout << iWARN << nameloc << "'" << "\n" << endi;

  index++;

  //  If I have a right subtree, index it
  if (tree->right != NULL)
  {
    index=index_vdw(tree->right, index);
  }

  return(index);
}
/*      END OF FUNCTION index_vdw      */

/************************************************************************/
/*                  */
/*      FUNCTION assign_vdw_index      */
/*                  */
/*   INPUTS:                */
/*  atomtype - atom type to find          */
/*  atom_ptr - pointer to the atom structure to find vdw paramters  */
/*       for              */
/*                  */
/*   OUTPUTS:                */
/*  the vdw_index field of the atom structure is populated    */
/*                  */
/*  This function searches the binary tree of vdw parameters so     */
/*   that an index can be assigned to this atom.  If the parameter is   */
/*   is found, then the index is assigned.  If the parameter is not     */
/*   found, then NAMD terminates.          */
/*                  */
/************************************************************************/
#ifdef MEM_OPT_VERSION
void Parameters::assign_vdw_index(char *atomtype, AtomCstInfo *atom_ptr)
#else    
void Parameters::assign_vdw_index(char *atomtype, Atom *atom_ptr)
#endif
{
  struct vdw_params *ptr;    //  Current position in trees
  int found=0;      //  Flag 1->found match
  int comp_code;      //  return code from strcasecmp

  /*  Check to make sure the files have all been read    */
  if (!AllFilesRead)
  {
    NAMD_die("Tried to assign vdw index before all parameter files were read");
  }

  /*  Start at the top            */
  ptr=vdwp;

  /*  While we haven't found a match, and we haven't reached      */
  /*  the bottom of the tree, compare the atom passed in with     */
  /*  the current value and decide if we have a match, or if not, */
  /*  which way to go            */
  while (!found && (ptr!=NULL))
  {
    comp_code = strcasecmp(atomtype, ptr->atomname);

    if (comp_code == 0)
    {
      /*  Found a match!        */
      atom_ptr->vdw_type=ptr->index;
      found=1;
    }
    else if (comp_code < 0)
    {
      /*  Go to the left        */
      ptr=ptr->left;
    }
    else
    {
      /*  Go to the right        */
      ptr=ptr->right;
    }
  }

  //****** BEGIN CHARMM/XPLOR type changes
  if (!found)
  {
    // since CHARMM allows wildcards "*" in vdw typenames
    // we have to look again if necessary, this way, if
    // we already had an exact match, this is never executed
    size_t windx;                      //  wildcard index

    /*  Start again at the top                                */
    ptr=vdwp;
  
     while (!found && (ptr!=NULL))
     {
  
       // get index of wildcard wildcard, get index
       windx= strcspn(ptr->atomname,"*"); 
       if (windx == strlen(ptr->atomname))
       {
         // there is no wildcard here
         comp_code = strcasecmp(atomtype, ptr->atomname);   
       }
       else
       {
         comp_code = strncasecmp(atomtype, ptr->atomname, windx); 
       }  

       if (comp_code == 0)
       {
         /*  Found a match!                                */
         atom_ptr->vdw_type=ptr->index;
         found=1;
         char errbuf[100];
         sprintf(errbuf,"VDW TYPE NAME %s MATCHES PARAMETER TYPE NAME %s",
                        atomtype, ptr->atomname);
         int i;
         for(i=0; i<error_msgs.size(); i++) {
           if ( strcmp(errbuf,error_msgs[i]) == 0 ) break;
         }
         if ( i == error_msgs.size() ) {
           char *newbuf = new char[strlen(errbuf)+1];
           strcpy(newbuf,errbuf);
           error_msgs.add(newbuf);
           iout << iWARN << newbuf << "\n" << endi;
         }
       }
       else if (comp_code < 0)
       {
          /*  Go to the left                                */
                ptr=ptr->left;
       }
       else
       {
         /*  Go to the right                                */
                ptr=ptr->right;
       }
     
     }
                
  }
  //****** END CHARMM/XPLOR type changes

  /*  Make sure we found it          */
  if (!found)
  {
    char err_msg[100];

    sprintf(err_msg, "DIDN'T FIND vdW PARAMETER FOR ATOM TYPE %s",
       atomtype);
    NAMD_die(err_msg);
  }

  return;
}
/*      END OF FUNCTION assign_vdw_index    */

/************************************************************************
 * FUNCTION get_table_pair_params
 *
 * Inputs:
 * atom1 - atom type for atom 1
 * atom2 - atom type for atom 2
 * K - an integer value for the table type to populate
 *
 * Outputs:
 * If a match is found, K is populated and 1 is returned. Otherwise,
 * 0 is returned.
 *
 * This function finds the proper type index for tabulated nonbonded 
 * interactions between two atoms. If no such interactions are found,
 * the atoms are assumed to interact through standard VDW potentials.
 * 
 ************************************************************************/

int Parameters::get_table_pair_params(Index ind1, Index ind2, int *K) {
  IndexedTablePair *ptr;
  Index temp;
  int found=FALSE;

  ptr=tab_pair_tree;
  //
  //  We need the smaller type in ind1, so if it isn't already that 
  //  way, switch them        */
  if (ind1 > ind2)
  {
    temp = ind1;
    ind1 = ind2;
    ind2 = temp;
  }

  /*  While we haven't found a match and we're not at the end  */
  /*  of the tree, compare the bond passed in with the tree  */
  while (!found && (ptr!=NULL))
  {
//    printf("Comparing %i with %i and %i with %i\n", ind1, ptr->ind1, ind2, ptr->ind2);
    if ( (ind1 == ptr->ind1) && (ind2 == ptr->ind2) )
    {
       found = TRUE;
    }
    else if ( (ind1 < ptr->ind1) || 
        ( (ind1==ptr->ind1) && (ind2 < ptr->ind2) ) )
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  If we found a match, assign the values      */
  if (found)
  {
    *K = ptr->K;
    return(TRUE);
  }
  else
  {
    return(FALSE);
  }
}
/*      END OF FUNCTION get_table_pair_params    */

/************************************************************************/
/*                  */
/*      FUNCTION get_vdw_pair_params      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  A - A value to populate            */
/*  B - B value to populate            */
/*  A14 - A 1-4 value to populate          */
/*  B14 - B 1-4 value to populate          */
/*                  */
/*   OUTPUTS:                */
/*  If a match is found, A, B, A14, and B14 are all populated and a */
/*   1 is returned.  Otherwise, a 0 is returned.      */
/*                    */
/*  This function finds a set of vdw_pair paramters.  It is given   */
/*   the two types of atoms involved.  This is the only paramter for    */
/*   which a match is NOT guaranteed.  There will only be a match if    */
/*   there are specific van der waals parameters for the two atom types */
/*   involved.                */
/*                  */
/************************************************************************/

int Parameters::get_vdw_pair_params(Index ind1, Index ind2, Real *A, 
        Real *B, Real *A14, Real *B14)

{
  IndexedVdwPair *ptr;    //  Current location in tree
  Index temp;      //  Temporary value for swithcing
          // values
  int found=FALSE;    //  Flag 1-> found a match

  ptr=vdw_pair_tree;

  //  We need the smaller type in ind1, so if it isn't already that 
  //  way, switch them        */
  if (ind1 > ind2)
  {
    temp = ind1;
    ind1 = ind2;
    ind2 = temp;
  }

  /*  While we haven't found a match and we're not at the end  */
  /*  of the tree, compare the bond passed in with the tree  */
  while (!found && (ptr!=NULL))
  {
    if ( (ind1 == ptr->ind1) && (ind2 == ptr->ind2) )
    {
       found = TRUE;
    }
    else if ( (ind1 < ptr->ind1) || 
        ( (ind1==ptr->ind1) && (ind2 < ptr->ind2) ) )
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  If we found a match, assign the values      */
  if (found)
  {
    *A = ptr->A;
    *B = ptr->B;
    *A14 = ptr->A14;
    *B14 = ptr->B14;

    return(TRUE);
  }
  else
  {
    return(FALSE);
  }
}
/*      END OF FUNCTION get_vdw_pair_params    */


/************************************************************************/
/*                  */
/*        FUNCTION assign_bond_index    */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  bond_ptr - pointer to bond structure to populate    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by bond_ptr is populated    */
/*                  */
/*  This function finds a bond in the binary tree of bond values    */
/*   and assigns its index.  If the bond is found, than the bond_type   */
/*   field of the bond structure is populated.  If the parameter is     */
/*   not found, NAMD will terminate.          */
/*                  */
/************************************************************************/

void Parameters::assign_bond_index(char *atom1, char *atom2, Bond *bond_ptr)

{
  struct bond_params *ptr;  //  Current location in tree
  int found=0;      //  Flag 1-> found a match
  int cmp_code;      //  return code from strcasecmp
  char tmp_name[15];    //  Temporary atom name

  /*  Check to make sure the files have all been read    */
  if (!AllFilesRead)
  {
    NAMD_die("Tried to assign bond index before all parameter files were read");
  }

  /*  We need atom1 < atom2, so if that's not the way they        */
  /*  were passed, flip them          */
  if (strcasecmp(atom1, atom2) > 0)
  {
    strcpy(tmp_name, atom1);
    strcpy(atom1, atom2);
    strcpy(atom2, tmp_name);
  }

  /*  Start at the top            */
  ptr=bondp;

  /*  While we haven't found a match and we're not at the end  */
  /*  of the tree, compare the bond passed in with the tree  */
  while (!found && (ptr!=NULL))
  {
    cmp_code=strcasecmp(atom1, ptr->atom1name);

    if (cmp_code == 0)
    {
      cmp_code=strcasecmp(atom2, ptr->atom2name);
    }

    if (cmp_code == 0)
    {
      /*  Found a match        */
      found=1;
      bond_ptr->bond_type = ptr->index;
    }
    else if (cmp_code < 0)
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  Check to see if we found anything        */
  if (!found)
  {
    if ((strcmp(atom1, "DRUD")==0 || strcmp(atom2, "DRUD")==0)
        && (strcmp(atom1, "X")!=0 && strcmp(atom2, "X")!=0)) {
      /* try a wildcard DRUD X match for this Drude bond */
      char a1[8] = "DRUD", a2[8] = "X";
      return assign_bond_index(a1, a2, bond_ptr);  /* recursive call */
    }
    else {
      char err_msg[128];

      sprintf(err_msg, "UNABLE TO FIND BOND PARAMETERS FOR %s %s (ATOMS %i %i)",
         atom1, atom2, bond_ptr->atom1+1, bond_ptr->atom2+1);
      NAMD_die(err_msg);
    }
  }

  return;
}
/*      END OF FUNCTION assign_bond_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_angle_index      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  atom3 - atom type for atom 3          */
/*  angle_ptr - pointer to angle structure to populate    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by angle_ptr is populated    */
/*                  */
/*  This function assigns an angle index to a specific angle.  */
/*   It searches the binary tree of angle parameters for the appropriate*/
/*   values.  If they are found, the index is assigned.  If they are    */
/*   not found, then NAMD will terminate.        */
/*                  */
/************************************************************************/

void Parameters::assign_angle_index(char *atom1, char *atom2, char*atom3,
          Angle *angle_ptr, int notFoundIndex)

{
  struct angle_params *ptr;  //  Current position in tree
  int comp_val;      //  value from strcasecmp
  int found=0;      //  flag 1->found a match
  char tmp_name[15];    //  Temporary atom name

  /*  Check to make sure the files have all been read    */
  if (!AllFilesRead)
  {
    NAMD_die("Tried to assign angle index before all parameter files were read");
  }

  /*  We need atom1 < atom3.  If that was not what we were   */
  /*  passed, switch them            */
  if (strcasecmp(atom1, atom3) > 0)
  {
    strcpy(tmp_name, atom1);
    strcpy(atom1, atom3);
    strcpy(atom3, tmp_name);
  }

  /*  Start at the top            */
  ptr=anglep;

  /*  While we don't have a match and we haven't reached the  */
  /*  bottom of the tree, compare values        */
  while (!found && (ptr != NULL))
  {
    comp_val = strcasecmp(atom1, ptr->atom1name);

    if (comp_val == 0)
    {
      /*  Atom 1 matches, so compare atom 2    */
      comp_val = strcasecmp(atom2, ptr->atom2name);
      
      if (comp_val == 0)
      {
        /*  Atoms 1&2 match, try atom 3    */
        comp_val = strcasecmp(atom3, ptr->atom3name);
      }
    }

    if (comp_val == 0)
    {
      /*  Found a match        */
      found = 1;
      angle_ptr->angle_type = ptr->index;
    }
    else if (comp_val < 0)
    {
      /*  Go left          */
      ptr=ptr->left;
    }
    else
    {
      /*  Go right          */
      ptr=ptr->right;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "UNABLE TO FIND ANGLE PARAMETERS FOR %s %s %s (ATOMS %i %i %i)",
       atom1, atom2, atom3, angle_ptr->atom1+1, angle_ptr->atom2+1, angle_ptr->atom3+1);

    if ( notFoundIndex ) {
      angle_ptr->angle_type = notFoundIndex;
      iout << iWARN << err_msg << "\n" << endi;
      return;
    } else NAMD_die(err_msg);
  }

  return;
}
/*      END OF FUNCTION assign_angle_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_dihedral_index      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  atom3 - atom type for atom 3          */
/*  atom4 - atom type for atom 4          */
/*  dihedral_ptr - pointer to dihedral structure to populate  */
/*  multiplicity - Multiplicity to assign to this bond    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by dihedral_ptr is populated    */
/*                  */
/*  This function searchs the linked list of dihedral parameters for*/
/*   a given bond.  If a match is found, the dihedral type is assigned. */
/*   If no match is found, NAMD terminates        */
/*                  */
/************************************************************************/

void Parameters::assign_dihedral_index(char *atom1, char *atom2, char *atom3,
        char *atom4, Dihedral *dihedral_ptr,
        int multiplicity, int notFoundIndex)

{
  struct dihedral_params *ptr;  //  Current position in list
  int found=0;      //  Flag 1->found a match

  /*  Start at the begining of the list        */
  ptr=dihedralp;

  /*  While we haven't found a match and we haven't reached       */
  /*  the end of the list, keep looking        */
  while (!found && (ptr!=NULL))
  {
    /*  Do a linear search through the linked list of   */
    /*  dihedral parameters.  Since the list is arranged    */
    /*  with wildcard paramters at the end of the list, we  */
    /*  can simply do a linear search and be guaranteed that*/
    /*  we will find exact matches before wildcard matches. */
    /*  Also, we must check for an exact match, and a match */
    /*  in reverse, since they are really the same          */
    /*  physically.            */
    if ( ( ptr->atom1wild || (strcasecmp(ptr->atom1name, atom1)==0) ) && 
         ( ptr->atom2wild || (strcasecmp(ptr->atom2name, atom2)==0) ) &&
         ( ptr->atom3wild || (strcasecmp(ptr->atom3name, atom3)==0) ) &&
         ( ptr->atom4wild || (strcasecmp(ptr->atom4name, atom4)==0) ) ) 
    {
      /*  Found an exact match      */
      found=1;
    }
    else if ( ( ptr->atom4wild || (strcasecmp(ptr->atom4name, atom1)==0) ) &&
              ( ptr->atom3wild || (strcasecmp(ptr->atom3name, atom2)==0) ) &&
              ( ptr->atom2wild || (strcasecmp(ptr->atom2name, atom3)==0) ) &&
              ( ptr->atom1wild || (strcasecmp(ptr->atom1name, atom4)==0) ) )
    {
      /*  Found a reverse match      */
      found=1;
    }
    else
    {
      /*  Didn't find a match, go to the next node  */
      ptr=ptr->next;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "UNABLE TO FIND DIHEDRAL PARAMETERS FOR %s %s %s %s (ATOMS %i %i %i %i)",
       atom1, atom2, atom3, atom4, dihedral_ptr->atom1+1, dihedral_ptr->atom2+1,
       dihedral_ptr->atom3+1, dihedral_ptr->atom4+1);
    
    if ( notFoundIndex ) {
      dihedral_ptr->dihedral_type = notFoundIndex;
      iout << iWARN << err_msg << "\n" << endi;
      return;
    } else NAMD_die(err_msg);
  }

  //  Check to make sure the number of multiples specified in the psf
  //  file doesn't exceed the number of parameters in the parameter
  //  files
  if (multiplicity > maxDihedralMults[ptr->index])
  {
    char err_msg[257];

    sprintf(err_msg, "Multiplicity of Paramters for diehedral bond %s %s %s %s of %d exceeded", atom1, atom2, atom3, atom4, maxDihedralMults[ptr->index]);
    NAMD_die(err_msg);
  }

  //  If the multiplicity from the current bond is larger than that
  //  seen in the past, increase the multiplicity for this bond
  if (multiplicity > dihedral_array[ptr->index].multiplicity)
  {
    dihedral_array[ptr->index].multiplicity = multiplicity;
  }

  dihedral_ptr->dihedral_type = ptr->index;

  return;
}
/*      END OF FUNCTION assign_dihedral_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_improper_index      */
/*                  */
/*   INPUTS:                */
/*  atom1 - atom type for atom 1          */
/*  atom2 - atom type for atom 2          */
/*  atom3 - atom type for atom 3          */
/*  atom4 - atom type for atom 4          */
/*  improper_ptr - pointer to improper structure to populate  */
/*   multiplicity - Multiplicity to assign to this bond    */
/*                  */
/*   OUTPUTS:                */
/*  the structure pointed to by improper_ptr is populated    */
/*                  */
/*  This function searchs the linked list of improper parameters for*/
/*   a given bond.  If a match is found, the improper_type is assigned. */
/*   If no match is found, NAMD will terminate.        */
/*                  */
/************************************************************************/

void Parameters::assign_improper_index(char *atom1, char *atom2, char *atom3,
        char *atom4, Improper *improper_ptr,
        int multiplicity)

{
  struct improper_params *ptr;  //  Current position in list
  int found=0;      //  Flag 1->found a match

  /*  Start at the head of the list        */
  ptr=improperp;

  /*  While we haven't fuond a match and haven't reached the end  */
  /*  of the list, keep looking          */
  while (!found && (ptr!=NULL))
  {
    /*  Do a linear search through the linked list of   */
    /*  improper parameters.  Since the list is arranged    */
    /*  with wildcard paramters at the end of the list, we  */
    /*  can simply do a linear search and be guaranteed that*/
    /*  we will find exact matches before wildcard matches. */
    /*  Also, we must check for an exact match, and a match */
    /*  in reverse, since they are really the same          */
    /*  physically.            */
    if ( ( (strcasecmp(ptr->atom1name, atom1)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom2)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom3)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom4name, atom4)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) )
    {
      /*  Found an exact match      */
      found=1;
    }
    else if ( ( (strcasecmp(ptr->atom4name, atom1)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom2)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom3)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom1name, atom4)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) )
    {
      /*  Found a reverse match      */
      found=1;
    }
    else
    {
      /*  Didn't find a match, go to the next node  */
      ptr=ptr->next;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "UNABLE TO FIND IMPROPER PARAMETERS FOR %s %s %s %s (ATOMS %i %i %i %i)",
       atom1, atom2, atom3, atom4, improper_ptr->atom1+1, improper_ptr->atom2+1,
       improper_ptr->atom3+1, improper_ptr->atom4+1);
    
    NAMD_die(err_msg);
  }

  //  Check to make sure the number of multiples specified in the psf
  //  file doesn't exceed the number of parameters in the parameter
  //  files
  if (multiplicity > maxImproperMults[ptr->index])
  {
    char err_msg[257];

    sprintf(err_msg, "Multiplicity of Paramters for improper bond %s %s %s %s of %d exceeded", atom1, atom2, atom3, atom4, maxImproperMults[ptr->index]);
    NAMD_die(err_msg);
  }

  //  If the multiplicity from the current bond is larger than that
  //  seen in the past, increase the multiplicity for this bond
  if (multiplicity > improper_array[ptr->index].multiplicity)
  {
    improper_array[ptr->index].multiplicity = multiplicity;
  }

  /*  Assign the constants          */
  improper_ptr->improper_type = ptr->index;

  return;
}
/*      END OF FUNCTION assign_improper_index    */

/************************************************************************/
/*                  */
/*      FUNCTION assign_crossterm_index      */
/*                  */
/************************************************************************/

void Parameters::assign_crossterm_index(char *atom1, char *atom2, char *atom3,
        char *atom4, char *atom5, char *atom6, char *atom7,
        char *atom8, Crossterm *crossterm_ptr)
{
  struct crossterm_params *ptr;  //  Current position in list
  int found=0;      //  Flag 1->found a match

  /*  Start at the head of the list        */
  ptr=crosstermp;

  /*  While we haven't fuond a match and haven't reached the end  */
  /*  of the list, keep looking          */
  while (!found && (ptr!=NULL))
  {
    /*  Do a linear search through the linked list of   */
    /*  crossterm parameters.  Since the list is arranged    */
    /*  with wildcard paramters at the end of the list, we  */
    /*  can simply do a linear search and be guaranteed that*/
    /*  we will find exact matches before wildcard matches. */
    /*  Also, we must check for an exact match, and a match */
    /*  in reverse, since they are really the same          */
    /*  physically.            */
    if ( ( (strcasecmp(ptr->atom1name, atom1)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom2)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom3)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom4name, atom4)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) )
    {
      /*  Found an exact match      */
      found=1;
    }
    else if ( ( (strcasecmp(ptr->atom4name, atom1)==0) || 
           (strcasecmp(ptr->atom4name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom3name, atom2)==0) || 
           (strcasecmp(ptr->atom3name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom2name, atom3)==0) || 
           (strcasecmp(ptr->atom2name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom1name, atom4)==0) || 
           (strcasecmp(ptr->atom1name, "X")==0) ) )
    {
      /*  Found a reverse match      */
      found=1;
    }
    if ( ! found ) {
      /*  Didn't find a match, go to the next node  */
      ptr=ptr->next;
      continue;
    }
    found = 0;
    if ( ( (strcasecmp(ptr->atom5name, atom5)==0) || 
           (strcasecmp(ptr->atom5name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom6name, atom6)==0) || 
           (strcasecmp(ptr->atom6name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom7name, atom7)==0) || 
           (strcasecmp(ptr->atom7name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom8name, atom8)==0) || 
           (strcasecmp(ptr->atom8name, "X")==0) ) )
    {
      /*  Found an exact match      */
      found=1;
    }
    else if ( ( (strcasecmp(ptr->atom8name, atom5)==0) || 
           (strcasecmp(ptr->atom8name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom7name, atom6)==0) || 
           (strcasecmp(ptr->atom7name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom6name, atom7)==0) || 
           (strcasecmp(ptr->atom6name, "X")==0) ) &&
       ( (strcasecmp(ptr->atom5name, atom8)==0) || 
           (strcasecmp(ptr->atom5name, "X")==0) ) )
    {
      /*  Found a reverse match      */
      found=1;
    }
    if ( ! found ) {
      /*  Didn't find a match, go to the next node  */
      ptr=ptr->next;
    }
  }

  /*  Make sure we found a match          */
  if (!found)
  {
    char err_msg[128];

    sprintf(err_msg, "UNABLE TO FIND CROSSTERM PARAMETERS FOR %s  %s  %s  %s  %s  %s  %s  %s\n"
      "(ATOMS %i %i %i %i %i %i %i %i)",
      atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8,
      crossterm_ptr->atom1+1, crossterm_ptr->atom1+2, crossterm_ptr->atom1+3, crossterm_ptr->atom4+1,
      crossterm_ptr->atom1+5, crossterm_ptr->atom1+6, crossterm_ptr->atom1+7, crossterm_ptr->atom8+1);
    
    NAMD_die(err_msg);
  }

  /*  Assign the constants          */
  crossterm_ptr->crossterm_type = ptr->index;

  return;
}
/*      END OF FUNCTION assign_improper_index    */

/************************************************************************/
/*                  */
/*      FUNCTION free_bond_tree        */
/*                  */
/*   INPUTS:                */
/*  bond_ptr - pointer to bond tree to free        */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a bond paramter tree.  It makes recursive calls to   */
/*   free the left an right subtress, and then frees the head.  It is   */
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_bond_tree(struct bond_params *bond_ptr)

{
  if (bond_ptr->left != NULL)
  {
    free_bond_tree(bond_ptr->left);
  }

  if (bond_ptr->right != NULL)
  {
    free_bond_tree(bond_ptr->right);
  }

  delete bond_ptr;

  return;
}
/*      END OF FUNCTION free_bond_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION free_angle_tree      */
/*                  */
/*   INPUTS:                */
/*  angle_ptr - pointer to angle tree to free      */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a angle paramter tree.  It makes recursive calls to  */
/*   free the left an right subtress, and then frees the head.  It is   */
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_angle_tree(struct angle_params *angle_ptr)

{
  if (angle_ptr->left != NULL)
  {
    free_angle_tree(angle_ptr->left);
  }

  if (angle_ptr->right != NULL)
  {
    free_angle_tree(angle_ptr->right);
  }

  delete angle_ptr;

  return;
}
/*      END OF FUNCTION free_angle_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION free_dihedral_list      */
/*                  */
/*   INPUTS:                */
/*  dih_ptr - pointer to the list to free        */
/*                  */
/*  this function frees a linked list of dihedral parameters.  It   */
/*   is only called by the destructor.          */
/*                  */
/************************************************************************/

void Parameters::free_dihedral_list(struct dihedral_params *dih_ptr)

{
  struct dihedral_params *ptr;  //  Current position in list
  struct dihedral_params *next; //  Next position in list

  ptr=dih_ptr;

  while (ptr != NULL)
  {
    next=ptr->next;
    delete ptr;
    ptr=next;
  }

  return;
}
/*      END OF FUNCTION free_dihedral_list    */

/************************************************************************/
/*                  */
/*      FUNCTION free_improper_list      */
/*                  */
/*   INPUTS:                */
/*  imp_ptr - pointer to the list to free        */
/*                  */
/*  this function frees a linked list of improper parameters.  It   */
/*   is only called by the destructor.          */
/*                  */
/************************************************************************/

void Parameters::free_improper_list(struct improper_params *imp_ptr)

{
  struct improper_params *ptr;  //  Current position in list
  struct improper_params *next; //  Next position in list

  ptr=imp_ptr;

  while (ptr != NULL)
  {
    next=ptr->next;
    delete ptr;
    ptr=next;
  }

  return;
}
/*      END OF FUNCTION free_improper_list    */
    
/************************************************************************/
/*                  */
/*      FUNCTION free_crossterm_list      */
/*                  */
/*   INPUTS:                */
/*  imp_ptr - pointer to the list to free        */
/*                  */
/*  this function frees a linked list of crossterm parameters.  It   */
/*   is only called by the destructor.          */
/*                  */
/************************************************************************/

void Parameters::free_crossterm_list(struct crossterm_params *imp_ptr)

{
  struct crossterm_params *ptr;  //  Current position in list
  struct crossterm_params *next; //  Next position in list

  ptr=imp_ptr;

  while (ptr != NULL)
  {
    next=ptr->next;
    delete ptr;
    ptr=next;
  }

  return;
}
/*      END OF FUNCTION free_crossterm_list    */
    
/************************************************************************/
/*                  */
/*      FUNCTION free_vdw_tree        */
/*                  */
/*   INPUTS:                */
/*  vdw_ptr - pointer to vdw tree to free        */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a vdw paramter tree.  It makes recursive calls to    */
/*   free the left an right subtress, and then frees the head.  It is   */
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_vdw_tree(struct vdw_params *vdw_ptr)

{
  if (vdw_ptr->left != NULL)
  {
    free_vdw_tree(vdw_ptr->left);
  }

  if (vdw_ptr->right != NULL)
  {
    free_vdw_tree(vdw_ptr->right);
  }

  delete vdw_ptr;

  return;
}
/*      END OF FUNCTION free_vdw_tree      */

/************************************************************************/
/*                  */
/*      FUNCTION free_vdw_pair_list      */
/*                  */
/*  This function frees the vdw_pair_list        */
/*                  */
/************************************************************************/

void Parameters::free_vdw_pair_list()
{
   struct vdw_pair_params *ptr, *next;
   
   ptr=vdw_pairp;
   
   while (ptr != NULL)
   {
      next = ptr->next;
      
      delete ptr;
      
      ptr = next;
   }
   
   vdw_pairp = NULL;
}
/*      END OF FUNCTION free_vdw_pair_list    */

/************************************************************************/
/*                  */
/*      FUNCTION free_nbthole_pair_list      */
/*                  */
/*  This function frees the nbthole_pair_list        */
/*                  */
/************************************************************************/

void Parameters::free_nbthole_pair_list()
{
   struct nbthole_pair_params *ptr, *next;

   ptr=nbthole_pairp;

   while (ptr != NULL)
   {
      next = ptr->next;

      delete ptr;

      ptr = next;
   }

   nbthole_pairp = NULL;
}
/*      END OF FUNCTION free_vdw_pair_list    */

/************************************************************************
 * FUNCTION free_table_pair_tree
 *
 * Free a table_pair_tree given a pointer to its head
 * **********************************************************************/

void Parameters::free_table_pair_tree(IndexedTablePair *table_pair_ptr) {
  if (table_pair_ptr->left != NULL)
  {
    free_table_pair_tree(table_pair_ptr->left);
  }

  if (table_pair_ptr->right != NULL)
  {
    free_table_pair_tree(table_pair_ptr->right);
  }

  delete table_pair_ptr;

  return;
}


/************************************************************************/
/*                  */
/*      FUNCTION free_vdw_pair_tree      */
/*                  */
/*   INPUTS:                */
/*  vdw_pair_ptr - pointer to vdw_pair tree to free      */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a vdw_pair paramter tree.  It makes recursive calls  */
/*   to free the left an right subtress, and then frees the head.  It is*/
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_vdw_pair_tree(IndexedVdwPair *vdw_pair_ptr)

{
  if (vdw_pair_ptr->left != NULL)
  {
    free_vdw_pair_tree(vdw_pair_ptr->left);
  }

  if (vdw_pair_ptr->right != NULL)
  {
    free_vdw_pair_tree(vdw_pair_ptr->right);
  }

  delete vdw_pair_ptr;

  return;
}
/*      END OF FUNCTION free_vdw_pair_tree    */

/************************************************************************/
/*                  */
/*      FUNCTION free_nbthole_pair_tree      */
/*                  */
/*   INPUTS:                */
/*  nbthole_pair_ptr - pointer to nbthole_pair tree to free      */
/*                  */
/*  this is a recursive function that is used to free the memory    */
/*   allocated for a nbthole_pair paramter tree.  It makes recursive calls  */
/*   to free the left an right subtress, and then frees the head.  It is*/
/*   only called by the destructor          */
/*                  */
/************************************************************************/

void Parameters::free_nbthole_pair_tree(IndexedNbtholePair *nbthole_pair_ptr)

{
  if (nbthole_pair_ptr->left != NULL)
  {
    free_nbthole_pair_tree(nbthole_pair_ptr->left);
  }

  if (nbthole_pair_ptr->right != NULL)
  {
    free_nbthole_pair_tree(nbthole_pair_ptr->right);
  }

  delete nbthole_pair_ptr;

  return;
}
/*      END OF FUNCTION free_nbthole_pair_tree    */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_bond_params      */
/*                  */
/*   INPUTS:                */
/*  tree - the bond binary tree to traverse        */
/*                  */
/*  This is a recursive call used for debugging purposes that  */
/*   prints out all the bond paramters in the bond parameter binary  */
/*   search tree. It is only called by print_bond_params    */
/*                  */
/************************************************************************/

void Parameters::traverse_bond_params(struct bond_params *tree)

{
  if (tree==NULL)
    return;

  if (tree->left != NULL)
  {
    traverse_bond_params(tree->left);
  }

  DebugM(3,"BOND " <<  tree->atom1name << "  " << tree->atom2name \
      << " index=" << tree->index << " k=" << tree->forceconstant \
      << " x0=" << tree->distance);

  if (tree->right != NULL)
  {
    traverse_bond_params(tree->right);
  }
}
/*      END OF FUNCTION traverse_bond_params    */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_angle_params      */
/*                  */
/*   INPUTS:                */
/*  tree - the angle binary tree to traverse      */
/*                  */
/*  This is a recursive call used for debugging purposes that  */
/*   prints out all the angle paramters in the angle parameter binary  */
/*   search tree. It is only called by print_angle_params    */
/*                  */
/************************************************************************/

void Parameters::traverse_angle_params(struct angle_params *tree)

{
  if (tree==NULL)
    return;

  if (tree->left != NULL)
  {
    traverse_angle_params(tree->left);
  }
  DebugM(3,"ANGLE " << tree->atom1name << "  " << tree->atom2name \
      << "  " << tree->atom3name << " index=" << tree->index \
      << " k=" << tree->forceconstant << " theta0=" << tree->angle \
      );

  if (tree->right != NULL)
  {
    traverse_angle_params(tree->right);
  }
}
/*      END OF FUNCTION traverse_angle_params    */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_dihedral_params    */
/*                  */
/*   INPUTS:                */
/*  list - the dihedral linked list to traverse      */
/*                  */
/*  This is a call used for debugging purposes that prints out all  */
/*   the bond paramters in the dihedral parameter linked list. It is    */
/*   only called by print_dihedral_params.        */
/*                  */
/************************************************************************/

void Parameters::traverse_dihedral_params(struct dihedral_params *list)

{
  int i;

  while (list != NULL)
  {
    DebugM(3,"DIHEDRAL  " << list->atom1name << "  " \
        << list->atom2name << "  " << list->atom3name \
        << "  " << list->atom4name << " index=" \
        << list->index \
        << " multiplicity=" << list->multiplicity << "\n");
        
    for (i=0; i<list->multiplicity; i++)
    {
      DebugM(3,"k=" << list->values[i].k \
          << " n=" << list->values[i].n  \
          << " delta=" << list->values[i].delta);
    }

    list=list->next;
  }
}
/*      END OF FUNCTION traverse_dihedral_params  */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_improper_params    */
/*                  */
/*   INPUTS:                */
/*  list - the improper linked list to traverse      */
/*                  */
/*  This is a call used for debugging purposes that prints out all  */
/*   the improper paramters in the improper parameter linked list. It is*/
/*   only called by print_improper_params.        */
/*                  */
/************************************************************************/

void Parameters::traverse_improper_params(struct improper_params *list)

{
  int i;

  while (list != NULL)
  {
    DebugM(3,"Improper  " << list->atom1name << "  " \
        << list->atom2name << "  " << list->atom3name \
        << "  " << list->atom4name << " index="  \
        << list->index  \
        << " multiplicity=" << list->multiplicity << "\n");

    for (i=0; i<list->multiplicity; i++)
    {
       DebugM(3,"k=" << list->values[i].k \
           << " n=" << list->values[i].n \
           << " delta=" << list->values[i].delta);
    }

    list=list->next;
  }
}
/*      END OF FUNCTION traverse_improper_params  */


/************************************************************************/
/*                  */
/*      FUNCTION traverse_vdw_params      */
/*                  */
/*   INPUTS:                */
/*  tree - the vw binary tree to traverse        */
/*                  */
/*  This is a recursive call used for debugging purposes that  */
/*   prints out all the vdw paramters in the vdw parameter binary  */
/*   search tree. It is only called by print_vdw_params      */
/*                  */
/************************************************************************/

void Parameters::traverse_vdw_params(struct vdw_params *tree)

{
  if (tree==NULL)
    return;

  if (tree->left != NULL)
  {
    traverse_vdw_params(tree->left);
  }

  DebugM(3,"vdW " << tree->atomname << " index=" << tree->index \
      << " sigma=" << tree->sigma << " epsilon=" << \
      tree->epsilon << " sigma 1-4=" << tree->sigma14 \
      << " epsilon 1-4=" << tree->epsilon14 << std::endl);

  if (tree->right != NULL)
  {
    traverse_vdw_params(tree->right);
  }
}
/*      END OF FUNCTION traverse_vdw_params    */


/************************************************************************/
/*                  */
/*      FUNCTION traverse_vdw_pair_params    */
/*                  */
/*   INPUTS:                */
/*  list - the vdw_pair list to traverse        */
/*                  */
/*  This call simply prints out the vdw_pair list      */
/*                  */
/************************************************************************/

void Parameters::traverse_vdw_pair_params(struct vdw_pair_params *list)

{
  if (list==NULL)
    return;

  DebugM(3,"vdW PAIR  " << list->atom1name << "  "  \
      << list->atom2name << " A=" << list->A \
      << " B=" << list->B << " A 1-4=" \
      << list->A14 << " B 1-4=" << list->B14 \
      );

  traverse_vdw_pair_params(list->next);
}
/*      END OF FUNCTION traverse_vdw_pair_params  */

/************************************************************************/
/*                  */
/*      FUNCTION traverse_nbthole_pair_params    */
/*                  */
/*   INPUTS:                */
/*  list - the nbthole_pair list to traverse        */
/*                  */
/*  This call simply prints out the nbthole_pair list      */
/*                  */
/************************************************************************/

void Parameters::traverse_nbthole_pair_params(struct nbthole_pair_params *list)

{
  if (list==NULL)
    return;

  DebugM(3,"NBTHOLE PAIR  " << list->atom1name << "  "  \
      << list->atom2name << " tholeij =" << list->tholeij \
      );

  traverse_nbthole_pair_params(list->next);
}
/*      END OF FUNCTION traverse_nbthole_pair_params  */

/************************************************************************/
/*                  */
/*      FUNCTION print_bond_params      */
/*                  */
/*  This is a debugging routine used to print out all the bond  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_bond_params()
{
  DebugM(3,NumBondParams << " BOND PARAMETERS\n" \
      << "*****************************************"  \
      );

  traverse_bond_params(bondp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_angle_params      */
/*                  */
/*  This is a debugging routine used to print out all the angle  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_angle_params()
{
  DebugM(3,NumAngleParams << " ANGLE PARAMETERS\n"
      << "*****************************************" );
  traverse_angle_params(anglep);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_dihedral_params      */
/*                  */
/*  This is a debugging routine used to print out all the dihedral  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_dihedral_params()
{
  DebugM(3,NumDihedralParams << " DIHEDRAL PARAMETERS\n" \
      << "*****************************************" );

  traverse_dihedral_params(dihedralp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_improper_params      */
/*                  */
/*  This is a debugging routine used to print out all the improper  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_improper_params()
{
  DebugM(3,NumImproperParams << " IMPROPER PARAMETERS\n" \
      << "*****************************************" );

  traverse_improper_params(improperp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_vdw_params      */
/*                  */
/*  This is a debugging routine used to print out all the vdw  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_vdw_params()
{
  DebugM(3,NumVdwParams << " vdW PARAMETERS\n" \
      << "*****************************************" );

  traverse_vdw_params(vdwp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_vdw_pair_params      */
/*                  */
/*  This is a debugging routine used to print out all the vdw_pair  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_vdw_pair_params()
{
  DebugM(3,NumVdwPairParams << " vdW PAIR PARAMETERS\n" \
      << "*****************************************" );

  traverse_vdw_pair_params(vdw_pairp);
}

/************************************************************************/
/*                  */
/*      FUNCTION print_nbthole_pair_params      */
/*                  */
/*  This is a debugging routine used to print out all the nbthole_pair  */
/*  parameters                */
/*                  */
/************************************************************************/

void Parameters::print_nbthole_pair_params()
{
  DebugM(3,NumNbtholePairParams << " NBTHOLE PAIR PARAMETERS\n" \
      << "*****************************************" );

  traverse_nbthole_pair_params(nbthole_pairp);
} 

/************************************************************************/
/*                  */
/*      FUNCTION print_param_summary      */
/*                  */
/*  This function just prints out a brief summary of the paramters  */
/*  that have been read in.  It is intended for debugging purposes  */
/*                  */
/************************************************************************/

void Parameters::print_param_summary()
{
  iout << iINFO << "SUMMARY OF PARAMETERS:\n" 
       << iINFO << NumBondParams << " BONDS\n" 
       << iINFO << NumAngleParams << " ANGLES\n" << endi;
       if (cosAngles) {
         iout << iINFO << "  " << (NumAngleParams - NumCosAngles) << " HARMONIC\n"
              << iINFO << "  " << NumCosAngles << " COSINE-BASED\n" << endi;
       }
  iout << iINFO << NumDihedralParams << " DIHEDRAL\n"
       << iINFO << NumImproperParams << " IMPROPER\n"
       << iINFO << NumCrosstermParams << " CROSSTERM\n"
       << iINFO << NumVdwParams << " VDW\n"
       << iINFO << NumVdwPairParams << " VDW_PAIRS\n"
       << iINFO << NumNbtholePairParams << " NBTHOLE_PAIRS\n" << endi;
}


/************************************************************************/
/*                  */
/*      FUNCTION done_reading_structure      */
/*                  */
/*  This function is used to tell the Parameters object that the    */
/*  structure has been read in.  This is so that the Parameters object  */
/*  can now release the binary trees and linked lists that it was using */
/*  to search for parameters based on the atom type.  From this point   */
/*  on, only the arrays of parameter data will be used.  If this object */
/*  resides on any node BUT the master node, it will never even have    */
/*  these trees and lists.  For the master node, this just frees up     */
/*  some memory for better uses.          */
/*                  */
/************************************************************************/

void Parameters::done_reading_structure()

{
  if (bondp != NULL)
    free_bond_tree(bondp);

  if (anglep != NULL)
    free_angle_tree(anglep);

  if (dihedralp != NULL)
    free_dihedral_list(dihedralp);

  if (improperp != NULL)
    free_improper_list(improperp);

  if (crosstermp != NULL)
    free_crossterm_list(crosstermp);

  if (vdwp != NULL)
    free_vdw_tree(vdwp);

  //  Free the arrays used to track multiplicity for dihedrals
  //  and impropers
  if (maxDihedralMults != NULL)
    delete [] maxDihedralMults;

  if (maxImproperMults != NULL)
    delete [] maxImproperMults;

  bondp=NULL;
  anglep=NULL;
  dihedralp=NULL;
  improperp=NULL;
  crosstermp=NULL;
  vdwp=NULL;
  maxImproperMults=NULL;
  maxDihedralMults=NULL;
}
/*      END OF FUNCTION done_reading_structure    */

/************************************************************************/
/*                  */
/*      FUNCTION send_Parameters      */
/*                  */
/*  This function is used by the master node to broadcast the       */
/*   structure Parameters to all the other nodes.        */
/*                  */
/************************************************************************/

void Parameters::send_Parameters(MOStream *msg)
{
  Real *a1, *a2, *a3, *a4;        //  Temporary arrays for sending messages
  int *i1, *i2, *i3;      //  Temporary int array
  int i, j;      //  Loop counters
  Real **kvals;      //  Force constant values for dihedrals and impropers
  int **nvals;      //  Periodicity values for  dihedrals and impropers
  Real **deltavals;    //  Phase shift values for  dihedrals and impropers
  BigReal *pairC6, *pairC12; // JLai
  /*MOStream *msg=comm_obj->newOutputStream(ALLBUTME, STATICPARAMSTAG, BUFSIZE);
  if ( msg == NULL )
  {
    NAMD_die("memory allocation failed in Parameters::send_Parameters");
  }*/

  //  Send the bond parameters
  msg->put(NumBondParams);

  if (NumBondParams)
  {
    a1 = new Real[NumBondParams];
    a2 = new Real[NumBondParams];

    if ( (a1 == NULL) || (a2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<NumBondParams; i++)
    {
      a1[i] = bond_array[i].k;
      a2[i] = bond_array[i].x0;
    }

    msg->put(NumBondParams, a1)->put(NumBondParams, a2);

    delete [] a1;
    delete [] a2;
  }

  //  Send the angle parameters
  msg->put(NumAngleParams);

  if (NumAngleParams)
  {
    a1 = new Real[NumAngleParams];
    a2 = new Real[NumAngleParams];
    a3 = new Real[NumAngleParams];
    a4 = new Real[NumAngleParams];
    i1 = new int[NumAngleParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) ||
         (a4 == NULL) || (i1 == NULL))
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<NumAngleParams; i++)
    {
      a1[i] = angle_array[i].k;
      a2[i] = angle_array[i].theta0;
      a3[i] = angle_array[i].k_ub;
      a4[i] = angle_array[i].r_ub;
      i1[i] = angle_array[i].normal;
    }

    msg->put(NumAngleParams, a1)->put(NumAngleParams, a2);
    msg->put(NumAngleParams, a3)->put(NumAngleParams, a4);
    msg->put(NumAngleParams, i1);

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
    delete [] i1;
  }

  //  Send the dihedral parameters
  msg->put(NumDihedralParams);

  if (NumDihedralParams)
  {
    i1 = new int[NumDihedralParams];
    kvals = new Real *[MAX_MULTIPLICITY];
    nvals = new int *[MAX_MULTIPLICITY];
    deltavals = new Real *[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumDihedralParams];
      nvals[i] = new int[NumDihedralParams];
      deltavals[i] = new Real[NumDihedralParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    for (i=0; i<NumDihedralParams; i++)
    {
      i1[i] = dihedral_array[i].multiplicity;

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        kvals[j][i] = dihedral_array[i].values[j].k;
        nvals[j][i] = dihedral_array[i].values[j].n;
        deltavals[j][i] = dihedral_array[i].values[j].delta;
      }
    }

    msg->put(NumDihedralParams, i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->put(NumDihedralParams, kvals[i]);
      msg->put(NumDihedralParams, nvals[i]);
      msg->put(NumDihedralParams, deltavals[i]);

      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Send the improper parameters
  msg->put(NumImproperParams);

  if (NumImproperParams)
  {
    i1 = new int[NumImproperParams];
    kvals = new Real *[MAX_MULTIPLICITY];
    nvals = new int *[MAX_MULTIPLICITY];
    deltavals = new Real *[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumImproperParams];
      nvals[i] = new int[NumImproperParams];
      deltavals[i] = new Real[NumImproperParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    for (i=0; i<NumImproperParams; i++)
    {
      i1[i] = improper_array[i].multiplicity;

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        kvals[j][i] = improper_array[i].values[j].k;
        nvals[j][i] = improper_array[i].values[j].n;
        deltavals[j][i] = improper_array[i].values[j].delta;
      }
    }

    msg->put(NumImproperParams, i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->put(NumImproperParams, kvals[i]);
      msg->put(NumImproperParams, nvals[i]);
      msg->put(NumImproperParams, deltavals[i]);

      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Send the crossterm parameters
  msg->put(NumCrosstermParams);

  if (NumCrosstermParams)
  {
    for (i=0; i<NumCrosstermParams; ++i) {
      int nvals = CrosstermValue::dim * CrosstermValue::dim * 2 * 2;
      msg->put(nvals,&crossterm_array[i].c[0][0].d00);
    }
  }
  //  Send the GromacsPairs parameters
  // JLai
  msg->put(NumGromacsPairParams);
 
  if(NumGromacsPairParams) 
  {
      pairC6 = new BigReal[NumGromacsPairParams];
      pairC12 = new BigReal[NumGromacsPairParams];
      if ( (pairC6 == NULL) || (pairC12 == NULL) ) {
	  NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }

      for (i=0; i<NumGromacsPairParams; i++) {
	  pairC6[i] = gromacsPair_array[i].pairC6;
	  pairC12[i] = gromacsPair_array[i].pairC12;
      }

      msg->put(NumGromacsPairParams,pairC6);
      msg->put(NumGromacsPairParams,pairC12);

      delete [] pairC6;
      delete [] pairC12;
  }
  // End of JLai
  
  //
  //Send the energy table parameters
  msg->put(numenerentries);

  if (numenerentries) {
	  /*
    b1 = new Real[numenerentries];
    if (b1 == NULL) 
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    */
    
    msg->put(numenerentries, table_ener);
  }

  //  Send the vdw parameters
  msg->put(NumVdwParams);
  msg->put(NumVdwParamsAssigned);

  if (NumVdwParams)
  {
    a1 = new Real[NumVdwParams];
    a2 = new Real[NumVdwParams];
    a3 = new Real[NumVdwParams];
    a4 = new Real[NumVdwParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<NumVdwParams; i++)
    {
      a1[i] = vdw_array[i].sigma;
      a2[i] = vdw_array[i].epsilon;
      a3[i] = vdw_array[i].sigma14;
      a4[i] = vdw_array[i].epsilon14;
    }

    msg->put(NumVdwParams * (MAX_ATOMTYPE_CHARS+1), atomTypeNames);
    msg->put(NumVdwParams, a1);
    msg->put(NumVdwParams, a2);
    msg->put(NumVdwParams, a3);
    msg->put(NumVdwParams, a4);

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }
  
  //  Send the vdw pair parameters
  msg->put(NumVdwPairParams);
  
  if (NumVdwPairParams)
  {
    a1 = new Real[NumVdwPairParams];
    a2 = new Real[NumVdwPairParams];
    a3 = new Real[NumVdwPairParams];
    a4 = new Real[NumVdwPairParams];
    i1 = new int[NumVdwPairParams];
    i2 = new int[NumVdwPairParams];    

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (a4 == NULL) || 
         (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    
    vdw_pair_to_arrays(i1, i2, a1, a2, a3, a4, 0, vdw_pair_tree);
    
    msg->put(NumVdwPairParams, i1)->put(NumVdwPairParams, i2);
    msg->put(NumVdwPairParams, a1);
    msg->put(NumVdwPairParams, a2)->put(NumVdwPairParams, a3);
    msg->put(NumVdwPairParams, a4);
  }

  //  Send the nbthole pair parameters
  msg->put(NumNbtholePairParams);

  if (NumNbtholePairParams)
  {
    a1 = new Real[NumNbtholePairParams];
    a2 = new Real[NumNbtholePairParams];
    a3 = new Real[NumNbtholePairParams];
    i1 = new int[NumNbtholePairParams];
    i2 = new int[NumNbtholePairParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    nbthole_pair_to_arrays(i1, i2, a1, a2, a3, 0, nbthole_pair_tree);

   for (i=0; i<NumNbtholePairParams; i++)
   {
    nbthole_array[i].ind1 = i1[i];
    nbthole_array[i].ind2 = i2[i];
    nbthole_array[i].alphai = a1[i];
    nbthole_array[i].alphaj = a2[i];
    nbthole_array[i].tholeij = a3[i];
   }

    msg->put(NumNbtholePairParams, i1)->put(NumNbtholePairParams, i2);
    msg->put(NumNbtholePairParams, a1);
    msg->put(NumNbtholePairParams, a2)->put(NumNbtholePairParams, a3);
  }
  
  //  Send the table pair parameters
  //printf("Pairs: %i\n", NumTablePairParams);
  msg->put(NumTablePairParams);
  
  if (NumTablePairParams)
  {
    i1 = new int[NumTablePairParams];
    i2 = new int[NumTablePairParams];    
    i3 = new int[NumTablePairParams];

    if ( (i3 == NULL) || (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    
    table_pair_to_arrays(i1, i2, i3, 0, tab_pair_tree);
    
    msg->put(NumTablePairParams, i1)->put(NumTablePairParams, i2);
    msg->put(NumTablePairParams, i3);
  }

  // send the hydrogen bond parameters
  // hbondParams.create_message(msg);
  msg->end();
  delete msg;
}

/************************************************************************/
/*                  */
/*      FUNCTION receive_Parameters      */
/*                  */
/*  This function is used by all the client processes to receive    */
/*   the structure parameters from the master node.      */
/*                  */
/************************************************************************/

void Parameters::receive_Parameters(MIStream *msg)

{
  int i, j;      //  Loop counters
  Real *a1, *a2, *a3, *a4;  //  Temporary arrays to get data from message in
  int *i1, *i2, *i3;      //  Temporary int array to get data from message in
  IndexedVdwPair *new_node;  //  New node for vdw pair param tree
  IndexedTablePair *new_tab_node;
  Real **kvals;      //  Force constant values for dihedrals and impropers
  int **nvals;      //  Periodicity values for dihedrals and impropers
  Real **deltavals;    //  Phase shift values for dihedrals and impropers
  BigReal *pairC6, *pairC12;  // JLai

  //  Get the bonded parameters
  msg->get(NumBondParams);

  if (NumBondParams)
  {
    bond_array = new BondValue[NumBondParams];
    a1 = new Real[NumBondParams];
    a2 = new Real[NumBondParams];

    if ( (bond_array == NULL) || (a1 == NULL) || (a2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumBondParams, a1);
    msg->get(NumBondParams, a2);

    for (i=0; i<NumBondParams; i++)
    {
      bond_array[i].k = a1[i];
      bond_array[i].x0 = a2[i];
    }

    delete [] a1;
    delete [] a2;
  }

  //  Get the angle parameters
  msg->get(NumAngleParams);

  if (NumAngleParams)
  {
    angle_array = new AngleValue[NumAngleParams];
    a1 = new Real[NumAngleParams];
    a2 = new Real[NumAngleParams];
    a3 = new Real[NumAngleParams];
    a4 = new Real[NumAngleParams];
    i1 = new int[NumAngleParams];

    if ( (angle_array == NULL) || (a1 == NULL) || (a2 == NULL) ||
         (a3 == NULL) || (a4 == NULL) || (i1 == NULL))
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumAngleParams, a1);
    msg->get(NumAngleParams, a2);
    msg->get(NumAngleParams, a3);
    msg->get(NumAngleParams, a4);
    msg->get(NumAngleParams, i1);

    for (i=0; i<NumAngleParams; i++)
    {
      angle_array[i].k = a1[i];
      angle_array[i].theta0 = a2[i];
      angle_array[i].k_ub = a3[i];
      angle_array[i].r_ub = a4[i];
      angle_array[i].normal = i1[i];
    }

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
    delete [] i1;
  }

  //  Get the dihedral parameters
  msg->get(NumDihedralParams);

  if (NumDihedralParams)
  {
    dihedral_array = new DihedralValue[NumDihedralParams];

    i1 = new int[NumDihedralParams];
    kvals = new Real *[MAX_MULTIPLICITY];
    nvals = new int *[MAX_MULTIPLICITY];
    deltavals = new Real *[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) || (dihedral_array == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumDihedralParams];
      nvals[i] = new int[NumDihedralParams];
      deltavals[i] = new Real[NumDihedralParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    msg->get(NumDihedralParams, i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->get(NumDihedralParams, kvals[i]);
      msg->get(NumDihedralParams, nvals[i]);
      msg->get(NumDihedralParams, deltavals[i]);
    }

    for (i=0; i<NumDihedralParams; i++)
    {
      dihedral_array[i].multiplicity = i1[i];

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        dihedral_array[i].values[j].k = kvals[j][i];
        dihedral_array[i].values[j].n = nvals[j][i];
        dihedral_array[i].values[j].delta = deltavals[j][i];
      }
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Get the improper parameters
  msg->get(NumImproperParams);

  if (NumImproperParams)
  {
    improper_array = new ImproperValue[NumImproperParams];
    i1 = new int[NumImproperParams];
    kvals = new Real *[MAX_MULTIPLICITY];
    nvals = new int *[MAX_MULTIPLICITY];
    deltavals = new Real *[MAX_MULTIPLICITY];

    if ( (i1 == NULL) || (kvals == NULL) || (nvals == NULL) || 
         (deltavals == NULL) || (improper_array==NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      kvals[i] = new Real[NumImproperParams];
      nvals[i] = new int[NumImproperParams];
      deltavals[i] = new Real[NumImproperParams];

      if ( (kvals[i] == NULL) || (nvals[i] == NULL) || (deltavals[i] == NULL) )
      {
        NAMD_die("memory allocation failed in Parameters::send_Parameters");
      }
    }

    msg->get(NumImproperParams,i1);

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      msg->get(NumImproperParams,kvals[i]);
      msg->get(NumImproperParams,nvals[i]);
      msg->get(NumImproperParams,deltavals[i]);
    }

    for (i=0; i<NumImproperParams; i++)
    {
      improper_array[i].multiplicity = i1[i];

      for (j=0; j<MAX_MULTIPLICITY; j++)
      {
        improper_array[i].values[j].k = kvals[j][i];
        improper_array[i].values[j].n = nvals[j][i];
        improper_array[i].values[j].delta = deltavals[j][i];
      }
    }

    for (i=0; i<MAX_MULTIPLICITY; i++)
    {
      delete [] kvals[i];
      delete [] nvals[i];
      delete [] deltavals[i];
    }

    delete [] i1;
    delete [] kvals;
    delete [] nvals;
    delete [] deltavals;
  }

  //  Get the crossterm parameters
  msg->get(NumCrosstermParams);

  if (NumCrosstermParams)
  {
    crossterm_array = new CrosstermValue[NumCrosstermParams];

    for (i=0; i<NumCrosstermParams; ++i) {
      int nvals = CrosstermValue::dim * CrosstermValue::dim * 2 * 2;
      msg->get(nvals,&crossterm_array[i].c[0][0].d00);
    }
  }
  
  // Get GromacsPairs parameters
  // JLai
  msg->get(NumGromacsPairParams);
  
  if (NumGromacsPairParams)
  {
      gromacsPair_array = new GromacsPairValue[NumGromacsPairParams];
      pairC6 = new BigReal[NumGromacsPairParams];
      pairC12 = new BigReal[NumGromacsPairParams];
      
      if ( (pairC6 == NULL) || (pairC12 == NULL) ) {
	  NAMD_die("memory allocation failed in Parameters::receive_Parameters");
      }
    
      msg->get(NumGromacsPairParams,pairC6);
      msg->get(NumGromacsPairParams,pairC12);
      
      for (i=0; i<NumGromacsPairParams; ++i) {
	  gromacsPair_array[i].pairC6 = pairC6[i];
	  gromacsPair_array[i].pairC12 = pairC12[i];
      }

      delete [] pairC6;
      delete [] pairC12;
  }
  // JLai

  //Get the energy table
  msg->get(numenerentries);
  if (numenerentries > 0) {
    //fprintf(stdout, "Getting tables\n");
    //fprintf(infofile, "Trying to allocate table\n");
    table_ener = new BigReal[numenerentries];
    //fprintf(infofile, "Finished table allocation\n");
    if (table_ener==NULL)
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(numenerentries, table_ener);
  }

  //  Get the vdw parameters
  msg->get(NumVdwParams);
  msg->get(NumVdwParamsAssigned);

  if (NumVdwParams)
  {
          atomTypeNames = new char[NumVdwParams*(MAX_ATOMTYPE_CHARS+1)];
    vdw_array = new VdwValue[NumVdwParams];
    a1 = new Real[NumVdwParams];
    a2 = new Real[NumVdwParams];
    a3 = new Real[NumVdwParams];
    a4 = new Real[NumVdwParams];

    if ( (vdw_array==NULL) || (a1==NULL) || (a2==NULL) || (a3==NULL)
             || (a4==NULL) || (atomTypeNames==NULL))
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumVdwParams * (MAX_ATOMTYPE_CHARS+1), atomTypeNames);
    msg->get(NumVdwParams, a1);
    msg->get(NumVdwParams, a2);
    msg->get(NumVdwParams, a3);
    msg->get(NumVdwParams, a4);

    for (i=0; i<NumVdwParams; i++)
    {
      vdw_array[i].sigma = a1[i];
      vdw_array[i].epsilon = a2[i];
      vdw_array[i].sigma14 = a3[i];
      vdw_array[i].epsilon14 = a4[i];
    }

    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }
  
  //  Get the vdw_pair_parameters
  msg->get(NumVdwPairParams);
  
  if (NumVdwPairParams)
  {
    a1 = new Real[NumVdwPairParams];
    a2 = new Real[NumVdwPairParams];
    a3 = new Real[NumVdwPairParams];
    a4 = new Real[NumVdwPairParams];
    i1 = new int[NumVdwPairParams];
    i2 = new int[NumVdwPairParams];

    if ( (a1 == NULL) || (a2 == NULL) || (a3 == NULL) || (a4 == NULL) || 
         (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    
    msg->get(NumVdwPairParams, i1);
    msg->get(NumVdwPairParams, i2);
    msg->get(NumVdwPairParams, a1);
    msg->get(NumVdwPairParams, a2);
    msg->get(NumVdwPairParams, a3);
    msg->get(NumVdwPairParams, a4);
    
    for (i=0; i<NumVdwPairParams; i++)
    {
      new_node = (IndexedVdwPair *) malloc(sizeof(IndexedVdwPair));
      
      if (new_node == NULL)
      {
         NAMD_die("memory allocation failed in Parameters::receive_Parameters");
      }
      
      new_node->ind1 = i1[i];
      new_node->ind2 = i2[i];
      new_node->A = a1[i];
      new_node->A14 = a2[i];
      new_node->B = a3[i];
      new_node->B14 = a4[i];
      new_node->left = NULL;
      new_node->right = NULL;
      
      vdw_pair_tree = add_to_indexed_vdw_pairs(new_node, vdw_pair_tree);
    }
    
    delete [] i1;
    delete [] i2;
    delete [] a1;
    delete [] a2;
    delete [] a3;
    delete [] a4;
  }
 
  //  Get the nbthole_pair_parameters
  msg->get(NumNbtholePairParams); 
    
  if (NumNbtholePairParams)
  {
    nbthole_array = new NbtholePairValue[NumNbtholePairParams];
    a1 = new Real[NumNbtholePairParams];
    a2 = new Real[NumNbtholePairParams];
    a3 = new Real[NumNbtholePairParams];
    i1 = new int[NumNbtholePairParams];
    i2 = new int[NumNbtholePairParams];

    if ( (nbthole_array == NULL) || (a1 == NULL) || (a2 == NULL) || (a3 == NULL)
         || (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::receive_Parameters");
    }

    msg->get(NumNbtholePairParams, i1);
    msg->get(NumNbtholePairParams, i2);
    msg->get(NumNbtholePairParams, a1);
    msg->get(NumNbtholePairParams, a2);
    msg->get(NumNbtholePairParams, a3);

    for (i=0; i<NumNbtholePairParams; i++)
    {

      nbthole_array[i].ind1 = i1[i];
      nbthole_array[i].ind2 = i2[i];
      nbthole_array[i].alphai = a1[i];
      nbthole_array[i].alphaj = a2[i];
      nbthole_array[i].tholeij = a3[i];

    }

    delete [] i1;
    delete [] i2;
    delete [] a1;
    delete [] a2;
    delete [] a3;
  }

  //  Get the table_pair_parameters
  msg->get(NumTablePairParams);
  
  if (NumTablePairParams)
  {
    i1 = new int[NumTablePairParams];
    i2 = new int[NumTablePairParams];
    i3 = new int[NumTablePairParams];

    if ( (i3 == NULL) || (i1 == NULL) || (i2 == NULL) )
    {
      NAMD_die("memory allocation failed in Parameters::send_Parameters");
    }
    
    msg->get(NumTablePairParams, i1);
    msg->get(NumTablePairParams, i2);
    msg->get(NumTablePairParams, i3);
    
    for (i=0; i<NumTablePairParams; i++)
    {
      new_tab_node = (IndexedTablePair *) malloc(sizeof(IndexedTablePair));
      
      if (new_tab_node == NULL)
      {
         NAMD_die("memory allocation failed in Parameters::receive_Parameters");
      }
      
//      printf("Adding new table pair with ind1 %i ind2 %i k %i\n", i1[i], i2[i],i3[i]);
      new_tab_node->ind1 = i1[i];
      new_tab_node->ind2 = i2[i];
      new_tab_node->K = i3[i];
      new_tab_node->left = NULL;
      new_tab_node->right = NULL;
      
      tab_pair_tree = add_to_indexed_table_pairs(new_tab_node, tab_pair_tree);
    }
    
    delete [] i1;
    delete [] i2;
    delete [] i3;
  }
  
  // receive the hydrogen bond parameters
  // hbondParams.receive_message(msg);

  AllFilesRead = TRUE;

  delete msg;
}
/*      END OF FUNCTION receive_Parameters    */

/************************************************************************/
/*                  */
/*      FUNCTION convert_vdw_pairs      */
/*                  */
/*  This function converts the linked list of vdw_pairs indexed by  */
/*  atom name into a binary search tree of parameters stored by vdw     */
/*  type index.  This tree is what will be used for real when searching */
/*  for parameters during the simulation.        */
/*                  */
/************************************************************************/

void Parameters::convert_vdw_pairs()
   
{
   #ifdef MEM_OPT_VERSION
   AtomCstInfo atom_struct;
   #else
   Atom atom_struct;    //  Dummy structure for getting indexes
   #endif
   Index index1, index2;  //  Indexes for the two atoms
   IndexedVdwPair *new_node;  //  New node for tree
   struct vdw_pair_params *ptr, *next;  //  Pointers for traversing list
   
   ptr = vdw_pairp;
   
   //  Go down then entire list and insert each node into the 
   //  binary search tree
   while (ptr != NULL)
   {
      new_node = (IndexedVdwPair *) malloc(sizeof(IndexedVdwPair));
      
      if (new_node == NULL)
      {
   NAMD_die("memory allocation failed in Parameters::convert_vdw_pairs");
      }
      
      //  Get the vdw indexes for the two atoms.  This is kind of a hack
      //  using the goofy Atom structure, but hey, it works
      assign_vdw_index(ptr->atom1name, &atom_struct);
      index1 = atom_struct.vdw_type;
      assign_vdw_index(ptr->atom2name, &atom_struct);
      index2 = atom_struct.vdw_type;
      
      if (index1 > index2)
      {
   new_node->ind1 = index2;
   new_node->ind2 = index1;
      }
      else
      {
   new_node->ind1 = index1;
   new_node->ind2 = index2;
      }
           
      new_node->A = ptr->A;
      new_node->B = ptr->B;
      new_node->A14 = ptr->A14;
      new_node->B14 = ptr->B14;
      
      new_node->left = NULL;
      new_node->right = NULL;
      
      //  Add it to the tree
      vdw_pair_tree = add_to_indexed_vdw_pairs(new_node, vdw_pair_tree);
      
      //  Free the current node and move to the next
      next = ptr->next;
      
      delete ptr;
      
      ptr = next;
   }
   
   vdw_pairp = NULL;

}
/*      END OF FUNCTION convert_vdw_pairs    */

/************************************************************************/
/*                  */
/*      FUNCTION convert_nbthole_pairs      */
/*                  */
/*  This function converts the linked list of nbthole_pairs indexed by  */
/*  atom name into a binary search tree of parameters stored by vdw     */
/*  type index.  This tree is what will be used for real when searching */
/*  for parameters during the simulation.        */
/*                  */
/************************************************************************/

void Parameters::convert_nbthole_pairs()

{
   #ifdef MEM_OPT_VERSION
   AtomCstInfo atom_struct;
   #else
   Atom atom_struct;    //  Dummy structure for getting indexes
   #endif
   Index index1, index2;  //  Indexes for the two atoms
   IndexedNbtholePair *new_node;  //  New node for tree
   struct nbthole_pair_params *ptr, *next;  //  Pointers for traversing list

   ptr = nbthole_pairp;

   //  Go down then entire list and insert each node into the
   //  binary search tree
   while (ptr != NULL)
   {
      new_node = (IndexedNbtholePair *) malloc(sizeof(IndexedNbtholePair));

      if (new_node == NULL)
      {
   NAMD_die("memory allocation failed in Parameters::convert_nbthole_pairs");
      }

      //  Get the vdw indexes for the two atoms.  This is kind of a hack
      //  using the goofy Atom structure, but hey, it works
      assign_vdw_index(ptr->atom1name, &atom_struct);
      index1 = atom_struct.vdw_type;
      assign_vdw_index(ptr->atom2name, &atom_struct);
      index2 = atom_struct.vdw_type;

      if (index1 > index2)
      {
   new_node->ind1 = index2;
   new_node->ind2 = index1;
      }
      else
      {
   new_node->ind1 = index1;
   new_node->ind2 = index2;
      }

      new_node->alphai = ptr->alphai;
      new_node->alphaj = ptr->alphaj;
      new_node->tholeij = ptr->tholeij;

      new_node->left = NULL;
      new_node->right = NULL;

      //  Add it to the tree
      nbthole_pair_tree = add_to_indexed_nbthole_pairs(new_node, nbthole_pair_tree);

      //  Free the current node and move to the next
      next = ptr->next;

      delete ptr;

      ptr = next;
   }

   nbthole_pairp = NULL;

}
/*      END OF FUNCTION convert_nbthole_pairs    */

/************************************************************************/
/*                  */
/*      FUNCTION convert_table_pairs      */
/*                  */
/*  This function converts the linked list of table_pairs indexed by  */
/*  atom name into a binary search tree of parameters stored by table     */
/*  type index.  This tree is what will be used for real when searching */
/*  for parameters during the simulation.        */
/*                  */
/************************************************************************/

void Parameters::convert_table_pairs()
   
{
   #ifdef MEM_OPT_VERSION
   AtomCstInfo atom_struct;
   #else
   Atom atom_struct;    //  Dummy structure for getting indexes
   #endif
   Index index1, index2;  //  Indexes for the two atoms
   IndexedTablePair *new_node;  //  New node for tree
   struct table_pair_params *ptr, *next;  //  Pointers for traversing list
   
   ptr = table_pairp;
   
   //  Go down then entire list and insert each node into the 
   //  binary search tree
   while (ptr != NULL)
   {
      new_node = (IndexedTablePair *) malloc(sizeof(IndexedTablePair));
      
      if (new_node == NULL)
      {
   NAMD_die("memory allocation failed in Parameters::convert_table_pairs");
      }
      
      //  Get the vdw indexes for the two atoms.  This is kind of a hack
      //  using the goofy Atom structure, but hey, it works
      assign_vdw_index(ptr->atom1name, &atom_struct);
      index1 = atom_struct.vdw_type;
      assign_vdw_index(ptr->atom2name, &atom_struct);
      index2 = atom_struct.vdw_type;

//      printf("Adding new table pair with index1 %i, index2 %i, k %i\n", index1, index2, ptr->K);
      
      if (index1 > index2)
      {
   new_node->ind1 = index2;
   new_node->ind2 = index1;
      }
      else
      {
   new_node->ind1 = index1;
   new_node->ind2 = index2;
      }
           
      new_node->K = ptr->K;

      new_node->left = NULL;
      new_node->right = NULL;
      
      //  Add it to the tree
      tab_pair_tree = add_to_indexed_table_pairs(new_node, tab_pair_tree);
      
      //  Free the current node and move to the next
      next = ptr->next;
      
      delete ptr;
      
      ptr = next;
   }
   
   table_pairp = NULL;

}
/*      END OF FUNCTION convert_table_pairs    */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_indexed_table_pairs    */
/*                  */
/*   INPUTS:                */
/*  new_node - new node to be added to the tree      */
/*  tree - tree to add the node to          */
/*                  */
/*  This is a recursive function that adds a node to the    */
/*   binary search tree of table_pair parameters        */
/*                  */
/************************************************************************/

IndexedTablePair *Parameters::add_to_indexed_table_pairs(IndexedTablePair *new_node,
                 IndexedTablePair *tree)
   
{
   if (tree == NULL)
      return(new_node);
   
   if ( (new_node->ind1 < tree->ind1) || 
        ((new_node->ind1 == tree->ind1) && (new_node->ind2 < tree->ind2)) )
   {
      tree->left = add_to_indexed_table_pairs(new_node, tree->left);
   }
   else
   {
      tree->right = add_to_indexed_table_pairs(new_node, tree->right);
   }
   
   return(tree);
}
/*      END OF FUNCTION add_to_indexed_table_pairs  */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_indexed_vdw_pairs    */
/*                  */
/*   INPUTS:                */
/*  new_node - new node to be added to the tree      */
/*  tree - tree to add the node to          */
/*                  */
/*  This is a recursive function that adds a node to the    */
/*   binary search tree of vdw_pair parameters        */
/*                  */
/************************************************************************/

IndexedVdwPair *Parameters::add_to_indexed_vdw_pairs(IndexedVdwPair *new_node,
                 IndexedVdwPair *tree)
   
{
   if (tree == NULL)
      return(new_node);
   
   if ( (new_node->ind1 < tree->ind1) || 
        ((new_node->ind1 == tree->ind1) && (new_node->ind2 < tree->ind2)) )
   {
      tree->left = add_to_indexed_vdw_pairs(new_node, tree->left);
   }
   else
   {
      tree->right = add_to_indexed_vdw_pairs(new_node, tree->right);
   }
   
   return(tree);
}
/*      END OF FUNCTION add_to_indexed_vdw_pairs  */

/************************************************************************/
/*                  */
/*      FUNCTION add_to_indexed_nbthole_pairs    */
/*                  */
/*   INPUTS:                */
/*  new_node - new node to be added to the tree      */
/*  tree - tree to add the node to          */
/*                  */ 
/*  This is a recursive function that adds a node to the    */
/*   binary search tree of nbthole_pair parameters        */
/*                  */
/************************************************************************/

IndexedNbtholePair *Parameters::add_to_indexed_nbthole_pairs(IndexedNbtholePair *new_node,
                 IndexedNbtholePair *tree)

{
   if (tree == NULL)
      return(new_node);

   if ( (new_node->ind1 < tree->ind1) ||
        ((new_node->ind1 == tree->ind1) && (new_node->ind2 < tree->ind2)) )
   {
      tree->left = add_to_indexed_nbthole_pairs(new_node, tree->left);
   }
   else
   {
      tree->right = add_to_indexed_nbthole_pairs(new_node, tree->right);
   }

   return(tree);
}
/*      END OF FUNCTION add_to_indexed_nbthole_pairs  */

/************************************************************************/
/*                  */
/*      FUNCTION vdw_pair_to_arrays      */
/*                  */
/*   INPUTS:                */
/*  ind1_array - Array of index 1 values        */
/*  ind2_array - Array of index 2 values        */
/*  A - Array of A values            */
/*  A14 - Array of A 1-4 values          */
/*  B - Array of B values            */
/*  B14 - Array of B 1-4 values          */
/*  arr_index - current position in arrays        */
/*  tree - tree to traverse            */
/*                  */
/*  This is a recursive function that places all the entries of     */
/*   the tree passed in into arrays of values.  This is done so that    */
/*   the parameters can be sent from the master node to the other       */
/*   nodes.                */
/*                  */
/************************************************************************/

int Parameters::vdw_pair_to_arrays(int *ind1_array, int *ind2_array,
          Real *A, Real *A14,
          Real *B, Real *B14,
          int arr_index, IndexedVdwPair *tree)
      
{
   if (tree == NULL)
      return(arr_index);
   
   ind1_array[arr_index] = tree->ind1;
   ind2_array[arr_index] = tree->ind2;
   A[arr_index] = tree->A;
   A14[arr_index] = tree->A14;
   B[arr_index] = tree->B;
   B14[arr_index] = tree->B14;
   
   arr_index++;
   
   arr_index = vdw_pair_to_arrays(ind1_array, ind2_array, A, A14, B, B14,
          arr_index, tree->left);
   arr_index = vdw_pair_to_arrays(ind1_array, ind2_array, A, A14, B, B14,
          arr_index, tree->right);
   
   return(arr_index);
}
/*      END OF FUNCTION vdw_pair_to_arrays    */

/************************************************************************/
/*                  */
/*      FUNCTION nbthole_pair_to_arrays      */
/*                  */
/*   INPUTS:                */
/*  ind1_array - Array of index 1 values        */
/*  ind2_array - Array of index 2 values        */
/*  tholeij - Array of tholeij values            */
/*  arr_index - current position in arrays        */
/*  tree - tree to traverse            */
/*                  */
/*  This is a recursive function that places all the entries of     */
/*   the tree passed in into arrays of values.  This is done so that    */
/*   the parameters can be sent from the master node to the other       */
/*   nodes.                */
/*                  */
/************************************************************************/

int Parameters::nbthole_pair_to_arrays(int *ind1_array, int *ind2_array,
          Real *alphai, Real *alphaj, Real *tholeij,
          int arr_index, IndexedNbtholePair *tree)

{
   if (tree == NULL)
      return(arr_index);

   ind1_array[arr_index] = tree->ind1;
   ind2_array[arr_index] = tree->ind2;
   alphai[arr_index] = tree->alphai;
   alphaj[arr_index] = tree->alphaj;
   tholeij[arr_index] = tree->tholeij;

   arr_index++;

   arr_index = nbthole_pair_to_arrays(ind1_array, ind2_array, alphai,
          alphaj, tholeij, arr_index, tree->left);
   arr_index = nbthole_pair_to_arrays(ind1_array, ind2_array, alphai,
          alphaj, tholeij, arr_index, tree->right);

   return(arr_index);
}
/*      END OF FUNCTION nbthole_pair_to_arrays    */

/************************************************************************/
/*                  */
/*      FUNCTION table_pair_to_arrays      */
/*                  */
/*   INPUTS:                */
/*  ind1_array - Array of index 1 values        */
/*  ind2_array - Array of index 2 values        */
/*  K - Array of K values            */
/*                  */
/*  This is a recursive function that places all the entries of     */
/*   the tree passed in into arrays of values.  This is done so that    */
/*   the parameters can be sent from the master node to the other       */
/*   nodes.                */
/*                  */
/************************************************************************/

int Parameters::table_pair_to_arrays(int *ind1_array, int *ind2_array,
          int *K,
          int arr_index, IndexedTablePair *tree)
      
{
   if (tree == NULL)
      return(arr_index);
   
   ind1_array[arr_index] = tree->ind1;
   ind2_array[arr_index] = tree->ind2;
   K[arr_index] = tree->K;

   arr_index++;
   
   arr_index = table_pair_to_arrays(ind1_array, ind2_array, K,
          arr_index, tree->left);
   arr_index = table_pair_to_arrays(ind1_array, ind2_array, K,
          arr_index, tree->right);
   
   return(arr_index);
}
/*      END OF FUNCTION table_pair_to_arrays    */

/************************************************************************/
/*                  */
/*      FUNCTION Parameters        */
/*                  */
/*  This is the constructor for reading AMBER parameter data    */
/*                  */
/************************************************************************/

Parameters::Parameters(Ambertoppar *amber_data, BigReal vdw14)
{
  initialize();

  // Read in parm parameters
  read_parm(amber_data,vdw14);
}
/*      END OF FUNCTION Parameters      */


/************************************************************************/
/*                  */
/*      FUNCTION read_parm   */
/*                  */
/*   INPUTS:                */
/*  amber_data - AMBER data structure    */
/*                  */
/*  This function copys AMBER parameter data to the corresponding data  */
/*   structures      */
/*                  */
/************************************************************************/

void Parameters::read_parm(Ambertoppar *amber_data, BigReal vdw14)
{ 
  int i,j,ico,numtype,mul;
  IndexedVdwPair *new_node;

  if (!amber_data->data_read)
    NAMD_die("No data read from parm file yet!");

  // Copy bond parameters
  NumBondParams = amber_data->Numbnd;
  if (NumBondParams)
  { bond_array = new BondValue[NumBondParams];
    if (bond_array == NULL)
      NAMD_die("memory allocation of bond_array failed!");
  }
  for (i=0; i<NumBondParams; ++i)
  { bond_array[i].k = amber_data->Rk[i];
    bond_array[i].x0 = amber_data->Req[i];
  }

  // Copy angle parameters
  NumAngleParams = amber_data->Numang;
  if (NumAngleParams)
  { angle_array = new AngleValue[NumAngleParams];
    if (angle_array == NULL)
      NAMD_die("memory allocation of angle_array failed!");
  }
  for (i=0; i<NumAngleParams; ++i)
  { angle_array[i].k = amber_data->Tk[i];
    angle_array[i].theta0 = amber_data->Teq[i];
    // Amber has no k_ub and r_ub for angle parameters, so they're set to 0
    angle_array[i].k_ub = angle_array[i].r_ub = 0;
    // All angles are harmonic
    angle_array[i].normal = 1;
  }

  // Copy dihedral parameters
  // Note: If the periodicity is negative, it means the following
  //  entry is another term in a multitermed dihedral, and the
  //  absolute value is the true periodicity; in this case the
  //  following entry in "dihedral_array" should be left empty,
  //  NOT be skipped, in order to keep the next dihedral's index
  //  unchanged.
  NumDihedralParams = amber_data->Nptra;
  if (NumDihedralParams)
  { dihedral_array = new DihedralValue[amber_data->Nptra];
    if (dihedral_array == NULL)
      NAMD_die("memory allocation of dihedral_array failed!");
  }
  mul = 0;
  for (i=0; i<NumDihedralParams; ++i)
  { dihedral_array[i-mul].values[mul].k = amber_data->Pk[i];
    dihedral_array[i-mul].values[mul].n = int(fabs(amber_data->Pn[i])+0.5);
    if (dihedral_array[i-mul].values[mul].n == 0)
    { char err_msg[128];
      sprintf(err_msg, "The periodicity of dihedral # %d is zero!", i+1);
      NAMD_die(err_msg);
    }
    dihedral_array[i-mul].values[mul].delta = amber_data->Phase[i];
    // If the periodicity is positive, it means the following
    // entry is a new dihedral term.
    if (amber_data->Pn[i] > 0)
    { dihedral_array[i-mul].multiplicity = mul+1;
      mul = 0;
    }
    else if (++mul >= MAX_MULTIPLICITY)
    { char err_msg[181];
      sprintf(err_msg, "Multiple dihedral with multiplicity of %d greater than max of %d",
         mul+1, MAX_MULTIPLICITY);
      NAMD_die(err_msg);
    }
  }
  if (mul > 0)
    NAMD_die("Negative periodicity in the last dihedral entry!");

  // Copy VDW parameters: AMBER explicitly gives VDW parameters between every
  // 2 atom types, so the data are copied to vdw_pair_tree
  // In AMBER, all 1-4 VDW interactions are divided by factor vdw14
  NumVdwParamsAssigned = numtype = amber_data->Ntypes; // Number of atom types
  NumVdwPairParams = numtype * (numtype+1) / 2;
  for (i=0; i<numtype; ++i)
    for (j=i; j<numtype; ++j)
    { new_node = new IndexedVdwPair;
      if (new_node == NULL)
        NAMD_die("memory allocation of vdw_pair_tree failed!");
      new_node->ind1 = i;
      new_node->ind2 = j;
      new_node->left = new_node->right = NULL;
      // ico is the index of interaction between atom type i and j into
      // the parameter arrays. If it's positive, the interaction is
      // 6-12 VDW; otherwise it's 10-12 H-bond interaction. NAMD doesn't
      // have 10-12 term, so if ico is negative, then the 10-12
      // coefficients must be 0, otherwise die.
      ico = amber_data->Cno[numtype*i+j];
      if (ico>0)
      { new_node->A = amber_data->Cn1[ico-1];
        new_node->A14 = new_node->A / vdw14;
        new_node->B = amber_data->Cn2[ico-1];
        new_node->B14 = new_node->B / vdw14;
      }
      else if (amber_data->HB12[abs(ico)-1]==0.0 && amber_data->HB6[abs(ico)-1]==0.0)
      { new_node->A = new_node->A14 = new_node->B = new_node->B14 = 0.0;
        iout << iWARN << "Encounter 10-12 H-bond term\n";
      }
      else
        NAMD_die("Encounter non-zero 10-12 H-bond term!");
      // Add the node to the binary tree
      vdw_pair_tree = add_to_indexed_vdw_pairs(new_node, vdw_pair_tree);
    }
}
/*      END OF FUNCTION read_parm    */

/************************************************************************/
/*                                                                      */
/*      FUNCTION Parameters                                             */
/*                                                                      */
/*  This is the constructor for reading GROMACS parameter data          */
/*                                                                      */
/************************************************************************/

Parameters::Parameters(const GromacsTopFile *gf, Bool min)
{
  initialize();

  // Read in parm parameters
  read_parm(gf,min);
}
/*      END OF FUNCTION Parameters      */


/************************************************************************/
/*                                                                      */
/*      FUNCTION read_parm                                              */
/*                                                                      */
/*   INPUTS:                                                            */
/*  gf - GROMACS topology file                                          */
/*                                                                      */
/* This function copys GROMACS parameter data to the corresponding data */
/*   structures                                                         */
/*                                                                      */
/************************************************************************/

void Parameters::read_parm(const GromacsTopFile *gf, Bool min)
{ 
  int numtype;
  IndexedVdwPair *new_node;
  int i,j,funct;
  Real test1,test2;

  // Allocate space for all of the arrays first
  NumBondParams = gf->getNumBondParams();
  NumAngleParams = gf->getNumAngleParams();
  NumDihedralParams = gf->getNumDihedralParams();
  if (NumBondParams) {
    bond_array = new BondValue[NumBondParams];
    if (bond_array == NULL)
      NAMD_die("memory allocation of bond_array failed!");
  }
  if (NumDihedralParams) {
    dihedral_array = new DihedralValue[NumDihedralParams];
    if (dihedral_array == NULL)
      NAMD_die("memory allocation of dihedral_array failed!");
  }
  if (NumAngleParams) {
    angle_array = new AngleValue[NumAngleParams];
    if (angle_array == NULL)
      NAMD_die("memory allocation of angle_array failed!");
  }

  // Copy bond parameters
  // XXX Warning: we are discarding the GROMACS function type - since
  // NAMD does not let us choose between different spring models.
  for (i=0;i<NumBondParams;i++) {
    Real x0,k;
    gf->getBondParams(i,
                      &x0, // the bond length
                      &k,  // the spring constant
                      &funct);           // discarded
    bond_array[i].x0 = x0;
    bond_array[i].k = k;
  }

  // Copy angle parameters
  // XXX Warning: we are discarding the GROMACS function type here
  // too.
  for (i=0;i<NumAngleParams;i++) {
    Real theta0,k;
    gf->getAngleParams(i,
                       &theta0, // the angle size
                       &k,      // the spring constant
                       &funct);                // discarded
    angle_array[i].theta0 = theta0*PI/180;
    angle_array[i].k = k;
    // Gromacs has no Urey-Bradley angle parameters, so they're set to 0
    angle_array[i].k_ub = angle_array[i].r_ub = 0;
    angle_array[i].normal = 1;
  }

  // Copy dihedral parameters
  // Here we use the function type (carefully)
  for (i=0; i<NumDihedralParams; ++i) { // iterate over all dihedral angles
    Real c[6];
    int num_periods; // number of periods in one full rotation
    int funct;

    gf->getDihedralParams(i,c,&num_periods,&funct); // get the parameters
    
    switch(funct) {
    case 1: ////////// it is a proper dihedral
      dihedral_array[i].values[0].delta = c[0]*PI/180; // the phase shift
      dihedral_array[i].values[0].k = c[1]; // the spring constant
      dihedral_array[i].values[0].n = num_periods; // the periodicity
      dihedral_array[i].multiplicity = 1;
      break;
    case 2: ////////// it is an improper dihedral
      dihedral_array[i].values[0].delta = c[0]*PI/180; // the phase shift
      dihedral_array[i].values[0].k = c[1]; // the spring constant
      dihedral_array[i].values[0].n = 0; // 'periodicity'=0 for impropers
      dihedral_array[i].multiplicity = 1;
      break;
    case 3: ////////// it is a Ryckaert-Bellemans dihedral

      // first a quick check to make sure this is legal
      if(MAX_MULTIPLICITY < 5)
        NAMD_die("I can't do RB dihedrals with MAX_MULTIPLICITY < 5");
      dihedral_array[i].multiplicity = 5;

      // Next we negate every other term, since GROMACS does this
      // silly thing with psi = 180 - phi:
      for(j=0;j<6;j++) {
        if(j%2 == 1) c[j] = -c[j];
      }

      // Now fill up all the terms.  Each is k(1 + cos(n*x - delta))
      // so we first let delta = 0 and let n range from 1-6:
      for(j=0;j<5;j++) dihedral_array[i].values[j].delta = 0;
      for(j=0;j<5;j++) dihedral_array[i].values[j].n = j+1;

      // and now we have a sum of kn(1+cos(nx))
      // Gromacs RB gives you a sum of cn(cos(x))^n, so we have to
      // use trigonometry to compute the kn:
      dihedral_array[i].values[0].k =    1*c[1] + 3/4.*c[3] + 10/16.*c[5];
      dihedral_array[i].values[1].k =       1/2.*c[2] + 4/8.*c[4]        ;
      dihedral_array[i].values[2].k =             1/4.*c[3] +  5/16.*c[5];
      dihedral_array[i].values[3].k =                   1/8.*c[4]        ;
      dihedral_array[i].values[4].k =                          1/16.*c[5];

      // That was a lot of math, so we had better check it:
      // The constant term (which is missing) is c0 + 1/2 c2 + 3/8 c4
      test1 = 0;
      for(j=5;j>=0;j--) { // sum (..(c5 cos x + c4) cos x + c3)..) + c1
        test1 *= cos(0.5);
        test1 += c[j];
      }

      test2 = c[0]+1/2.*c[2]+3/8.*c[4];
      for(j=0;j<5;j++) { // sum k1 cos x + k2 cos 2x + ... 
        test2 += dihedral_array[i].values[j].k * cos((j+1)*0.5);
      }

      if(fabs(test1-test2) > 0.0001)
        NAMD_die("Internal error: failed to handle RB dihedrals");
      
      // Turn this on to have a look at the values if you *still*
      // don't believe that they are right!

      /*      iout << iINFO << "k: ";
              for(j=0;j<5;j++)
              iout  << dihedral_array[i].values[j].k << " ";
              iout << "\n" << endi;
              
              iout << iINFO << "c: ";
              for(j=0;j<6;j++)
              iout  << c[j] << " ";
              iout << "\n" << endi;*/
      
      break;
    default:
      NAMD_die("unknown dihedral type found");
    }
  }

  // Copy VDW parameters.

  Bool warned=false; // warned the user about extra LJ term yet?

  NumVdwParamsAssigned = numtype = gf->getNumAtomParams(); // # of atom types
  NumVdwPairParams = numtype * (numtype+1) / 2;
  for (i=0; i<numtype; i++) {
    for (j=i; j<numtype; j++) {

      // set up a node to put one set of VDW parameters in
      new_node = new IndexedVdwPair;
      if (new_node == NULL)
        NAMD_die("memory allocation of vdw_pair_tree failed!");
      new_node->ind1 = i;
      new_node->ind2 = j;
      new_node->left = new_node->right = NULL;

      gf->getVDWParams(i,j, &(new_node->B), &(new_node->A),
                       &(new_node->B14), &(new_node->A14));

      /* If we have any VDW radii equal to zero, atoms can just sit on
         each other during minimization.  So, we'll give a minimum of
         1.0 kcal*A^12 to the LJ-repulsion when we are minimizing.
         But a warning should be displayed to the user... */
      if(min && ( fabs(new_node->A) < 1.0 )) {
        new_node->A = 1.0;
        if(!warned) {
          iout << iWARN <<
            "Adding small LJ repulsion term to some atoms.\n" << endi;
          warned=true;
        }
      }

      vdw_pair_tree = add_to_indexed_vdw_pairs(new_node, vdw_pair_tree);
    }
  }

  // JLai
  // Allocate space for all of the GromacsPair arrays first
  int numPair, numLJPair, numGaussPair;
  Real *pairC6,*pairC12; // constants to define LJ potential
  int *atom1,*atom2; // atom indices for LJ code
  atom1 = 0;
  atom2 = 0;
  pairC6 = 0;
  pairC12 = 0;
  numPair = gf->getNumPair();
  NumGromacsPairParams = numLJPair = gf->getNumLJPair();
  if (numLJPair) {
      atom1   = new int[numLJPair];
      atom2   = new int[numLJPair];
      pairC6  = new Real[numLJPair];
      pairC12 = new Real[numLJPair];
      gromacsPair_array = new GromacsPairValue[numLJPair];
  }
 
  // Copy GromacsPair data into gromacsPair array structures
  const_cast<GromacsTopFile*>(gf)->getPairLJArrays2(atom1, atom2, pairC6, pairC12);
  for (i=0;i<numLJPair;i++) {
      gromacsPair_array[i].pairC6  = pairC6[i];
      gromacsPair_array[i].pairC12 = pairC12[i];
  }
  delete [] atom1;
  delete [] atom2;
  delete [] pairC6;
  delete [] pairC12;
}
/*      END OF FUNCTION read_parm    */

/************************************************************************/
/*                                                                      */
/*      FUNCTION read_ener_table                                          */
/*                                                                      */
/*   INPUTS:                                                            */
/*  simParams -- Simulation Parameters   */
/*                                                                      */
/*   This function reads energy tables from a file and places them into                                                       */
/*   memory.                                                                      */
/************************************************************************/

void Parameters::read_ener_table(SimParameters *simParams) {
	char* table_file = simParams->tabulatedEnergiesFile;
  char* interp_type = simParams->tableInterpType;
	FILE* enertable;
	enertable = fopen(table_file, "r");

	if (enertable == NULL) {
		NAMD_die("ERROR: Could not open energy table file!\n");
	}

	char tableline[121];
  char* newtypename;
  int numtypes;
	int atom1;
	int atom2;
	int distbin;
  int readcount;
	Real dist;
	Real energy;
	Real force;
	Real table_spacing;
	Real maxdist;

/* First find the header line and set the variables we need */
	iout << "Beginning table read\n" << endi;
	while(fgets(tableline,120,enertable)) {
		if (strncmp(tableline,"#",1)==0) {
			continue;
		}
    readcount = sscanf(tableline, "%i %f %f", &numtypes, &table_spacing, &maxdist);
    if (readcount != 3) {
      NAMD_die("ERROR: Couldn't parse table header line\n");
    }
    break;
  }

  if (maxdist < simParams->cutoff) {
    NAMD_die("Tabulated energies must at least extend to the cutoff distance\n");
  }

	if (maxdist > simParams->cutoff) {
		iout << "Info: Reducing max table size to match nonbond cutoff\n";
		maxdist = ceil(simParams->cutoff);
	}

/* Now allocate memory for the table; we know what we should be getting */
	numenerentries = 2 * numtypes * int (mynearbyint(maxdist/table_spacing) + 1);
	//Set up a full energy lookup table from a file
	//Allocate the table; layout is atom1 x atom2 x distance energy force
	fprintf(stdout, "Table has %i entries\n",numenerentries);
	//iout << "Allocating original energy table\n" << endi;
	if ( table_ener ) delete [] table_ener;
	table_ener = new BigReal[numenerentries];
  if ( table_types ) delete [] table_types;
  table_types = new char*[numtypes];

/* Finally, start reading the table */
  int numtypesread = 0; //Number of types read so far

	while(fgets(tableline,120,enertable)) {
		if (strncmp(tableline,"#",1)==0) {
			continue;
		}
    if (strncmp(tableline,"TYPE",4)==0) {
      // Read a new type
      newtypename = new char[5];
      int readcount = sscanf(tableline, "%*s %s", newtypename);
      if (readcount != 1) {
        NAMD_die("Couldn't read interaction type from TYPE line\n");
      }
//      printf("Setting table_types[%i] to %s\n", numtypesread, newtypename);
      table_types[numtypesread] = newtypename;

      if (numtypesread == numtypes) {
        NAMD_die("Error: Number of types in table doesn't match table header\n");
      }

      // Read the current energy type with the proper interpolation
      if (!strncasecmp(interp_type, "linear", 6)) {
        if (read_energy_type(enertable, numtypesread, table_ener, table_spacing, maxdist) != 0) {
          char err_msg[512];
          sprintf(err_msg, "Failed to read table energy (linear) type %s\n", newtypename);
          NAMD_die(err_msg);
        }
      } else if (!strncasecmp(interp_type, "cubic", 5)) {
        if (read_energy_type_bothcubspline(enertable, numtypesread, table_ener, table_spacing, maxdist) != 0) {
          char err_msg[512];
          sprintf(err_msg, "Failed to read table energy (cubic) type %s\n", newtypename);
          NAMD_die(err_msg);
        }
      } else {
        NAMD_die("Table interpolation type must be linear or cubic\n");
      }

      numtypesread++;
      continue;
    }
    //char err_msg[512];
    //sprintf(err_msg, "Unrecognized line in energy table file: %s\n", tableline);
    //NAMD_die(err_msg);
  }

  // See if we got what we expected
  if (numtypesread != numtypes) {
    char err_msg[512];
    sprintf(err_msg, "ERROR: Expected %i tabulated energy types but got %i\n", numtypes, numtypesread);
    NAMD_die(err_msg);
  }

  // Move the necessary information into simParams
  simParams->tableNumTypes = numtypes;
  simParams->tableSpacing = table_spacing;
  simParams->tableMaxDist = maxdist;
  tablenumtypes = numtypes;

  /*
xxxxxx
	int numtypes = simParams->tableNumTypes;

	//This parameter controls the table spacing
	BigReal table_spacing = 0.01;
	BigReal maxdist = 20.0;
	if (maxdist > simParams->cutoff) {
		iout << "Info: Reducing max table size to match nonbond cutoff\n";
		maxdist = ceil(simParams->cutoff);
	}

	numenerentries = (numtypes + 1) * numtypes * int (ceil(maxdist/table_spacing));
//	fprintf(stdout, "Table arithmetic: maxdist %f, table_spacing %f, numtypes %i, numentries %i\n", maxdist, table_spacing, numtypes, numenerentries);
	columnsize = 2 * mynearbyint(maxdist/table_spacing);
	rowsize = numtypes * columnsize;
	//Set up a full energy lookup table from a file
	//Allocate the table; layout is atom1 x atom2 x distance energy force
	fprintf(stdout, "Table has %i entries\n",numenerentries);
	//iout << "Allocating original energy table\n" << endi;
	if ( table_ener ) delete [] table_ener;
	table_ener = new Real[numenerentries];
	//
	//Set sentinel values for uninitialized table energies
	for (int i=0 ; i<numenerentries ; i++) {
		table_ener[i] = 1337.0;
	}
	Real compval = 1337.0;

	//    iout << "Entering table reading\n" << endi;
	//iout << "Finished allocating table\n" << endi;

	//Counter to make sure we read the right number of energy entries
	int numentries = 0;

	//Now, start reading from the file
	char* table_file = simParams->tabulatedEnergiesFile;
	FILE* enertable;
	enertable = fopen(table_file, "r");

	if (enertable == NULL) {
		NAMD_die("ERROR: Could not open energy table file!\n");
	}

	char tableline[121];
	int atom1;
	int atom2;
	int distbin;
	Real dist;
	Real energy;
	Real force;

	iout << "Beginning table read\n" << endi;
	//Read the table entries
	while(fgets(tableline,120,enertable)) {
//		IOut << "Processing line " << tableline << "\n" << endi;
		if (strncmp(tableline,"#",1)==0) {
			continue;
		}


		sscanf(tableline, "%i %i %f %f %f\n", &atom1, &atom2, &dist, &energy, &force);
		distbin = int(mynearbyint(dist/table_spacing));
//		fprintf(stdout, "Working on atoms %i and %i at distance %f\n",atom1,atom2,dist);
		if ((atom2 > atom1) || (distbin > int(mynearbyint(maxdist/table_spacing)))) {
//			fprintf(stdout,"Rejected\n");
//			fprintf(stdout, "Error: Throwing out energy line beyond bounds\n");
	//		fprintf(stdout, "Bad line: Atom 1: %i Atom 2: %i Distance Bin: %i Max Distance Bin: %i \n", atom1, atom2, distbin, int(mynearbyint(maxdist/table_spacing)));
		} else {
			//The magic formula for the number of columns from previous rows
			//in the triangular matrix is (2ni+i-i**2)/2
			//Here n is the number of types, and i is atom2
//			fprintf(stdout, "Input: atom1 %f atom2 %f columnsize %f \n", float(atom1), float(atom2), float(columnsize));
//			fprintf(stdout, "Testing arithmetic: Part1: %i Part2: %i Part3: %i Total: %i\n", columnsize*((2*numtypes*atom2 + atom2 - atom2*atom2)/2), columnsize*(atom1-atom2), 2*distbin, columnsize*((2*numtypes*atom2 + atom2 - atom2*atom2)/2) + columnsize*(atom1-atom2) + 2*distbin - 2);
			int eneraddress = columnsize*((2*numtypes*atom2 + atom2 - atom2*atom2)/2) + columnsize*(atom1-atom2) + 2*distbin - 2;
			int forceaddress = eneraddress + 1;
//				fprintf(stdout, "Tableloc: %i Atom 1: %i Atom 2: %i Distance Bin: %i Energy: %f Force: %f\n", eneraddress, atom1, atom2, distbin, table_ener[eneraddress], table_ener[forceaddress]);
		fflush(stdout);
//			fprintf(stdout, "Checking for dupes: Looking at: %f %f \n", table_ener[eneraddress], table_ener[forceaddress]);
			if ((table_ener[eneraddress] == compval && table_ener[forceaddress] == compval)) {
				numentries++;
				table_ener[eneraddress] = energy;
				table_ener[forceaddress] = force;
//				fprintf(stdout, "Tableloc: %i Atom 1: %i Atom 2: %i Distance Bin: %i Energy: %f Force: %f\n", eneraddress, atom1, atom2, distbin, table_ener[eneraddress], table_ener[forceaddress]);
				//table_ener[rowsize*atom2 + columnsize*atom1 + 2*distbin] = energy;
				//table_ener[rowsize*atom2 + columnsize*atom1 + 2*distbin + 1] = force;
//				fflush(stdout);
			} else {
//				fprintf(stdout,"Rejecting duplicate entry\n");
			}
		}
		//      iout << "Done with line\n"<< endi;
	}

	//    int should = numtypes * numtypes * (maxdist/table_spacing);
	//    cout << "Entries: " << numentries << " should be " << should << "\n" << endi;
//	int exptypes = ceil((numtypes+1) * numtypes * (maxdist/table_spacing));
//fprintf(stdout,"numtypes: %i maxdist: %f table_spacing: %f exptypes: %i\n",numtypes,maxdist,table_spacing);
	if (numentries != int(numenerentries/2)) {
		fprintf(stdout,"ERROR: Expected %i entries but got %i\n",numenerentries/2, numentries);
		NAMD_die("Number of energy table entries did not match expected value\n");
	}
	//      iout << "Done with table\n"<< endi;
	fprintf(stdout, "Numenerentries: %i\n",numenerentries/2);
  */
} /* END of function read_ener_table */ 

/**************************************************************************
 * FUNCTION read_energy_type_bothcubspline
 *
 * Read a single type block from an energy table file, using cubic spline interpolation
 * Unlike _cubspline, the forces are interpolated separately
 *
 * Inputs:
 *  enertable - File stream positioned at the start of the type block
 *  typeindex  - integer index of current type
 *  table_ener - pointer to array to be filled with table entries
 *  table_spacing - Spacing between table points (A)
 *  maxdist - Longest distance needed in table
 *
 * Return values:
 *  0 on normal exit
 *  1 if not enough entries were present to fill out the table
 *
 *  Note: enertable should be left positioned immediately BEFORE the next
 *  TYPE block starts
 *  **********************************************************************/

int Parameters::read_energy_type_bothcubspline(FILE* enertable, const int typeindex, BigReal* table_ener, const float table_spacing, const float maxdist) {

  char tableline[120];
  int i,j;

  /* Last values read from table */
  BigReal readdist;
  BigReal readener;
  BigReal readforce;
  BigReal spacing;
//  BigReal readforce;
  BigReal lastdist;
//  BigReal lastener;
//  BigReal lastforce;
//  readdist = -1.0;
//  readener = 0.0;
//  readforce = 0.0;

  /* Create arrays for holding the input values */
  std::vector<BigReal>  dists;
  std::vector<BigReal> enervalues;
  std::vector<BigReal> forcevalues;
  int numentries = 0;


  /* Keep track of where in the table we are */
  BigReal currdist;
  int distbin;
  currdist = 0.0;
  lastdist = -1.0;
  distbin = 0;

  // Read all of the values first -- we'll interpolate later
	while(fgets(tableline,120,enertable) && distbin <= (int) (mynearbyint(maxdist / table_spacing) + 1)) {
		if (strncmp(tableline,"#",1)==0) {
			continue;
		}
    if (strncmp(tableline,"TYPE",4)==0) {
      fseek(enertable, -1 * (long) strlen(tableline), SEEK_CUR); 
      break;
    }

    // Read an energy line from the table
    int readcount = sscanf(tableline, "%lf %lf %lf", &readdist, &readener, &readforce);

    //printf("Read an energy line: %g %g %g\n", readdist, readener, readforce);
    if (readcount != 3) {
      char err_msg[512];
      sprintf(err_msg, "ERROR: Failed to parse table line %s!\n", tableline);
      NAMD_die(err_msg);
    }

    //Sanity check the current entry
    if (readdist < lastdist) {
      NAMD_die("ERROR: Encountered badly ordered entries in energy table!\n");
    }

    lastdist = readdist;
    dists.push_back(readdist);
    enervalues.push_back(readener);
    forcevalues.push_back(readforce);
    numentries++;
  }

  // Check the spacing and make sure it is uniform
  if (dists[0] != 0.0) {
    NAMD_die("ERROR: First data point for energy table must be at r=0\n");
  }
  spacing = dists[1] - dists[0];
  for (i=2; i<(numentries - 1); i++) {
    BigReal myspacing;
    myspacing = dists[i] - dists[i-1];
    if (fabs(myspacing - spacing) > 0.00001) {
      printf("Bad spacing in table: %f should be %f (between distances %f and %f)\n", myspacing, spacing, dists[i-1], dists[i]);
      NAMD_die("ERROR: Nonuniform table encountered on cubic interpolation. Use a uniformly spaced table or switch to linear interpolation.\n");
    }
  }

// Perform cubic spline interpolation to get the energies and forces

  /* allocate spline coefficient matrix */
  // xe and xf / be and bf for energies and forces, respectively
  double* m = new double[numentries*numentries];
  double* xe = new double[numentries];
  double* xf = new double[numentries];
  double* be = new double[numentries];
  double* bf = new double[numentries];

  be[0] = 3 * (enervalues[1] - enervalues[0]);
  for (i=1; i < (numentries - 1); i++) {
//    printf("Control point %i at %f\n", i, enervalues[i]);
    be[i] = 3 * (enervalues[i+1] - enervalues[i-1]);
//    printf("be is %f\n", be[i]);
  }
  be[numentries - 1] = 3 * (enervalues[numentries - 1] - enervalues[numentries - 2]);

  bf[0] = 3 * (forcevalues[1] - forcevalues[0]);
  for (i=1; i < (numentries - 1); i++) {
//    printf("Control point %i at %f\n", i, forcevalues[i]);
    bf[i] = 3 * (forcevalues[i+1] - forcevalues[i-1]);
//    printf("bf is %f\n", bf[i]);
  }
  bf[numentries - 1] = 3 * (forcevalues[numentries - 1] - forcevalues[numentries - 2]);

  memset(m,0,numentries*numentries*sizeof(double));

  /* initialize spline coefficient matrix */
  m[0] = 2;
  for (i = 1;  i < numentries;  i++) {
    m[INDEX(numentries,i-1,i)] = 1;
    m[INDEX(numentries,i,i-1)] = 1;
    m[INDEX(numentries,i,i)] = 4;
  }
  m[INDEX(numentries,numentries-1,numentries-1)] = 2;

  /* Now we need to solve the equation M X = b for X */
  // Do this simultaneously for energy and force -- ONLY because we use the same matrix

  //Put m in upper triangular form and apply corresponding operations to b
  for (i=0; i<numentries; i++) {
    // zero the ith column in all rows below i
    const BigReal baseval = m[INDEX(numentries,i,i)];
    for (j=i+1; j<numentries; j++) {
      const BigReal myval = m[INDEX(numentries,j,i)];
      const BigReal factor = myval / baseval;

      for (int k=i; k<numentries; k++) {
        const BigReal subval = m[INDEX(numentries,i,k)];
        m[INDEX(numentries,j,k)] -= (factor * subval);
      }

      be[j] -= (factor * be[i]);
      bf[j] -= (factor * bf[i]);

    }
  }

  //Now work our way diagonally up from the bottom right to assign values to X
  for (i=numentries-1; i>=0; i--) {

    //Subtract the effects of previous columns
    for (j=i+1; j<numentries; j++) {
      be[i] -= ( m[INDEX(numentries,i,j)] * xe[j] );
      bf[i] -= ( m[INDEX(numentries,i,j)] * xf[j] );
    }

    xe[i] = be[i] / m[INDEX(numentries,i,i)];
    xf[i] = bf[i] / m[INDEX(numentries,i,i)];

  }

  // We now have the coefficient information we need to make the table
  // Now interpolate on each interval we want

  distbin = 0;
  int entriesperseg = (int) ceil(spacing / table_spacing);
  int table_prefix = 2 * typeindex * (int) (mynearbyint(maxdist / table_spacing) + 1);

  for (i=0; i<numentries-1; i++) {
    BigReal Ae,Be,Ce,De;
    BigReal Af,Bf,Cf,Df;
    currdist = dists[i];

//    printf("Interpolating on interval %i\n", i);

    // Set up the coefficients for this section
    Ae = enervalues[i];
    Be = xe[i];
    Ce = 3 * (enervalues[i+1] - enervalues[i]) - (2 * xe[i]) - xe[i+1];
    De = 2 * (enervalues[i] - enervalues[i+1]) + xe[i] + xe[i+1];

    Af = forcevalues[i];
    Bf = xf[i];
    Cf = 3 * (forcevalues[i+1] - forcevalues[i]) - (2 * xf[i]) - xf[i+1];
    Df = 2 * (forcevalues[i] - forcevalues[i+1]) + xf[i] + xf[i+1];

    // Go over the region of interest and fill in the table
    for (j=0; j<entriesperseg; j++) {
      const BigReal mydist = currdist + (j * table_spacing);
      const BigReal mydistfrac = (float) j / (entriesperseg - 1);
      distbin = (int) mynearbyint(mydist / table_spacing);
      if (distbin >= (int) mynearbyint(maxdist / table_spacing)) break;
      BigReal energy;
      BigReal force;

      energy = Ae + (Be * mydistfrac) + (Ce * mydistfrac * mydistfrac) + (De * mydistfrac * mydistfrac * mydistfrac);
      force = Af + (Bf * mydistfrac) + (Cf * mydistfrac * mydistfrac) + (Df * mydistfrac * mydistfrac * mydistfrac);

//      printf("Adding energy/force entry %f / %f in bins %i / %i for distance %f (%f)\n", energy, force, (table_prefix + 2 * distbin), (table_prefix + 2 * distbin + 1), mydist, mydistfrac);
      table_ener[table_prefix + 2 * distbin] = energy;
      table_ener[table_prefix + 2 * distbin + 1] = force;
      distbin++;
    }
    if (currdist >= maxdist) break;
  }

  //The procedure above leaves out the last entry -- add it explicitly
  distbin = (int) mynearbyint(maxdist / table_spacing);
//  printf("Adding energy/force entry %f / %f in bins %i / %i\n", enervalues[numentries - 1], 0.0, (table_prefix + 2 * distbin), (table_prefix + 2 * distbin + 1));
  table_ener[table_prefix + 2 * distbin] = enervalues[numentries - 1];
  table_ener[table_prefix + 2 * distbin + 1] = 0.0;
  distbin++;


  // Clean up and make sure everything worked ok
  delete m;
  delete xe;
  delete xf;
  delete be;
  delete bf;
  distbin--;
  printf("Testing: %i vs %i (from %f / %f)\n", distbin, (int) (mynearbyint(maxdist / table_spacing)), maxdist, table_spacing);
  if (distbin != (int) (mynearbyint(maxdist / table_spacing))) return 1;
  return 0;
} /* end read_energy_type_bothcubspline */

/**************************************************************************
 * FUNCTION read_energy_type_cubspline
 *
 * Read a single type block from an energy table file, using cubic spline interpolation
 *
 * Inputs:
 *  enertable - File stream positioned at the start of the type block
 *  typeindex  - integer index of current type
 *  table_ener - pointer to array to be filled with table entries
 *  table_spacing - Spacing between table points (A)
 *  maxdist - Longest distance needed in table
 *
 * Return values:
 *  0 on normal exit
 *  1 if not enough entries were present to fill out the table
 *
 *  Note: enertable should be left positioned immediately BEFORE the next
 *  TYPE block starts
 *  **********************************************************************/

int Parameters::read_energy_type_cubspline(FILE* enertable, const int typeindex, BigReal* table_ener, const float table_spacing, const float maxdist) {

  char tableline[120];
  int i,j;

  /* Last values read from table */
  BigReal readdist;
  BigReal readener;
  BigReal spacing;
//  BigReal readforce;
  BigReal lastdist;
//  BigReal lastener;
//  BigReal lastforce;
//  readdist = -1.0;
//  readener = 0.0;
//  readforce = 0.0;

  /* Create arrays for holding the input values */
  std::vector<BigReal>  dists;
  std::vector<BigReal> enervalues;
  int numentries = 0;


  /* Keep track of where in the table we are */
  BigReal currdist;
  int distbin;
  currdist = 0.0;
  lastdist = -1.0;
  distbin = 0;

  // Read all of the values first -- we'll interpolate later
	while(fgets(tableline,120,enertable) && distbin <= (int) (mynearbyint(maxdist / table_spacing) + 1)) {
		if (strncmp(tableline,"#",1)==0) {
			continue;
		}
    if (strncmp(tableline,"TYPE",4)==0) {
      fseek(enertable, -1 * (long) strlen(tableline), SEEK_CUR); 
      break;
    }

    // Read an energy line from the table
    int readcount = sscanf(tableline, "%lf %lf", &readdist, &readener);

   // printf("Read an energy line: %f %f %f\n", readdist, readener, readforce);
    if (readcount != 2) {
      char err_msg[512];
      sprintf(err_msg, "ERROR: Failed to parse table line %s!\n", tableline);
      NAMD_die(err_msg);
    }

    //Sanity check the current entry
    if (readdist < lastdist) {
      NAMD_die("ERROR: Encountered badly ordered entries in energy table!\n");
    }

    lastdist = readdist;
    dists.push_back(readdist);
    enervalues.push_back(readener);
    numentries++;
  }

  // Check the spacing and make sure it is uniform
  if (dists[0] != 0.0) {
    NAMD_die("ERROR: First data point for energy table must be at r=0\n");
  }
  spacing = dists[1] - dists[0];
  for (i=2; i<(numentries - 1); i++) {
    BigReal myspacing;
    myspacing = dists[i] - dists[i-1];
    if (fabs(myspacing - spacing) > 0.00001) {
      printf("Bad spacing in table: %f should be %f (between distances %f and %f)\n", myspacing, spacing, dists[i-1], dists[i]);
      NAMD_die("ERROR: Nonuniform table encountered on cubic interpolation. Use a uniformly spaced table or switch to linear interpolation.\n");
    }
  }

// Perform cubic spline interpolation to get the energies and forces

  /* allocate spline coefficient matrix */
  double* m = new double[numentries*numentries];
  double* x = new double[numentries];
  double* b = new double[numentries];

  b[0] = 3 * (enervalues[1] - enervalues[0]);
  for (i=1; i < (numentries - 1); i++) {
    printf("Control point %i at %f\n", i, enervalues[i]);
    b[i] = 3 * (enervalues[i+1] - enervalues[i-1]);
    printf("b is %f\n", b[i]);
  }
  b[numentries - 1] = 3 * (enervalues[numentries - 1] - enervalues[numentries - 2]);

  memset(m,0,numentries*numentries*sizeof(double));

  /* initialize spline coefficient matrix */
  m[0] = 2;
  for (i = 1;  i < numentries;  i++) {
    m[INDEX(numentries,i-1,i)] = 1;
    m[INDEX(numentries,i,i-1)] = 1;
    m[INDEX(numentries,i,i)] = 4;
  }
  m[INDEX(numentries,numentries-1,numentries-1)] = 2;

  /* Now we need to solve the equation M X = b for X */

  printf("Solving the matrix equation: \n");

  for (i=0; i<numentries; i++) {
    printf(" ( ");
    for (j=0; j<numentries; j++) {
      printf(" %6.3f,", m[INDEX(numentries, i, j)]);
    }
    printf(" )  ( D%-3i )  =  ( %6.3f )\n", i, b[i]);
  }

  //Put m in upper triangular form and apply corresponding operations to b
  for (i=0; i<numentries; i++) {
    // zero the ith column in all rows below i
    const BigReal baseval = m[INDEX(numentries,i,i)];
    for (j=i+1; j<numentries; j++) {
      const BigReal myval = m[INDEX(numentries,j,i)];
      const BigReal factor = myval / baseval;

      for (int k=i; k<numentries; k++) {
        const BigReal subval = m[INDEX(numentries,i,k)];
        m[INDEX(numentries,j,k)] -= (factor * subval);
      }

      b[j] -= (factor * b[i]);

    }
  }

  printf(" In upper diagonal form, equation is:\n");
  for (i=0; i<numentries; i++) {
    printf(" ( ");
    for (j=0; j<numentries; j++) {
      printf(" %6.3f,", m[INDEX(numentries, i, j)]);
    }
    printf(" )  ( D%-3i )  =  ( %6.3f )\n", i, b[i]);
  }

  //Now work our way diagonally up from the bottom right to assign values to X
  for (i=numentries-1; i>=0; i--) {

    //Subtract the effects of previous columns
    for (j=i+1; j<numentries; j++) {
      b[i] -= ( m[INDEX(numentries,i,j)] * x[j] );
    }

    x[i] = b[i] / m[INDEX(numentries,i,i)];

  }

  printf(" Solution vector is:\n\t(");
  for (i=0; i<numentries; i++) {
    printf(" %6.3f ", x[i]);
  }
  printf(" ) \n");

  // We now have the coefficient information we need to make the table
  // Now interpolate on each interval we want

  distbin = 0;
  int entriesperseg = (int) ceil(spacing / table_spacing);
  int table_prefix = 2 * typeindex * (int) (mynearbyint(maxdist / table_spacing) + 1);

  for (i=0; i<numentries-1; i++) {
    BigReal A,B,C,D;
    currdist = dists[i];

    printf("Interpolating on interval %i\n", i);

    // Set up the coefficients for this section
    A = enervalues[i];
    B = x[i];
    C = 3 * (enervalues[i+1] - enervalues[i]) - (2 * x[i]) - x[i+1];
    D = 2 * (enervalues[i] - enervalues[i+1]) + x[i] + x[i+1];

    printf("Coefficients for this interval: %f %f %f %f\n", A, B, C, D);

    // Go over the region of interest and fill in the table
    for (j=0; j<entriesperseg; j++) {
      const BigReal mydist = currdist + (j * table_spacing);
      const BigReal mydistfrac = (float) j / (entriesperseg - 1);
      distbin = (int) mynearbyint(mydist / table_spacing);
      if (distbin >= (int) mynearbyint(maxdist / table_spacing)) break;
      BigReal energy;
      BigReal force;

      energy = A + (B * mydistfrac) + (C * mydistfrac * mydistfrac) + (D * mydistfrac * mydistfrac * mydistfrac);
      force = B + (2 * C * mydistfrac) + (3 * D * mydistfrac * mydistfrac);
      // Multiply force by 1 / (interval length)
      force *= (1.0 / spacing);

      printf("Adding energy/force entry %f / %f in bins %i / %i for distance %f (%f)\n", energy, force, (table_prefix + 2 * distbin), (table_prefix + 2 * distbin + 1), mydist, mydistfrac);
      table_ener[table_prefix + 2 * distbin] = energy;
      table_ener[table_prefix + 2 * distbin + 1] = force;
      distbin++;
    }
    if (currdist >= maxdist) break;
  }

  //The procedure above leaves out the last entry -- add it explicitly
  distbin = (int) mynearbyint(maxdist / table_spacing);
  printf("Adding energy/force entry %f / %f in bins %i / %i\n", enervalues[numentries - 1], 0.0, (table_prefix + 2 * distbin), (table_prefix + 2 * distbin + 1));
  table_ener[table_prefix + 2 * distbin] = enervalues[numentries - 1];
  table_ener[table_prefix + 2 * distbin + 1] = 0.0;
  distbin++;


  // Clean up and make sure everything worked ok
  delete m;
  delete x;
  delete b;
  distbin--;
  printf("Testing: %i vs %i (from %f / %f)\n", distbin, (int) (mynearbyint(maxdist / table_spacing)), maxdist, table_spacing);
  if (distbin != (int) (mynearbyint(maxdist / table_spacing))) return 1;
  return 0;
} /* end read_energy_type_cubspline */

/**************************************************************************
 * FUNCTION read_energy_type
 *
 * Read a single type block from an energy table file
 *
 * Inputs:
 *  enertable - File stream positioned at the start of the type block
 *  typeindex  - integer index of current type
 *  table_ener - pointer to array to be filled with table entries
 *  table_spacing - Spacing between table points (A)
 *  maxdist - Longest distance needed in table
 *
 * Return values:
 *  0 on normal exit
 *  1 if not enough entries were present to fill out the table
 *
 *  Note: enertable should be left positioned immediately BEFORE the next
 *  TYPE block starts
 *  **********************************************************************/

int Parameters::read_energy_type(FILE* enertable, const int typeindex, BigReal* table_ener, const float table_spacing, const float maxdist) {

  char tableline[120];

  /* Last values read from table */
  BigReal readdist;
  BigReal readener;
  BigReal readforce;
  BigReal lastdist;
  BigReal lastener;
  BigReal lastforce;
  readdist = -1.0;
  readener = 0.0;
  readforce = 0.0;

  /* Keep track of where in the table we are */
  float currdist;
  int distbin;
  currdist = -1.0;
  distbin = -1;

	while(fgets(tableline,120,enertable) && distbin <= (int) (mynearbyint(maxdist / table_spacing) + 1)) {
    printf("At distance %f + %f vs. %f\n", currdist, table_spacing, maxdist);
		if (strncmp(tableline,"#",1)==0) {
			continue;
		}
    if (strncmp(tableline,"TYPE",4)==0) {
      break;
    }

    // Read an energy line from the table
    lastdist = readdist;
    lastener = readener;
    lastforce = readforce;
    int readcount = sscanf(tableline, "%lf %lf %lf", &readdist, &readener, &readforce);
    if (distbin == -1) {
      if (readdist != 0.0) {
        NAMD_die("ERROR: Energy/force tables must start at d=0.0\n");
      } else {
        distbin = 0;
        continue;
      }
    }
   // printf("Read an energy line: %f %f %f\n", readdist, readener, readforce);
    if (readcount != 3) {
      char err_msg[512];
      sprintf(err_msg, "ERROR: Failed to parse table line %s!\n", tableline);
      NAMD_die(err_msg);
    }

    //Sanity check the current entry
    if (readdist < lastdist) {
      NAMD_die("ERROR: Encountered badly ordered entries in energy table!\n");
    }

    currdist = lastdist;

    while (currdist <= readdist && distbin <= (int) (mynearbyint(maxdist / table_spacing))) {
      distbin = (int) (mynearbyint(currdist / table_spacing));
      int table_loc = 2 * (distbin + (typeindex * (1 + (int) mynearbyint(maxdist / table_spacing))));
      printf("Doing interpolation for energy between %f %f and %f %f: Dist %f\n", readener, readdist, lastener, lastdist, currdist);
      table_ener[table_loc] = interp_lin(readener, lastener, readdist, lastdist, currdist);
      table_ener[table_loc + 1] = interp_lin(readforce, lastforce, readdist, lastdist, currdist);
      printf("Adding energy/force entry: %f %f in distbin %i (distance %f) to address %i/%i\n", table_ener[table_loc], table_ener[table_loc + 1], distbin, currdist, table_loc, table_loc + 1);
      currdist += table_spacing;
      distbin++;
    }
  }

  // Go back one line, since we should now be into the next TYPE block
  fseek(enertable, -1 * (long) strlen(tableline), SEEK_CUR); 

  // Clean up and make sure everything worked ok
  distbin--;
  printf("Testing: %i vs %i (from %f / %f)\n", distbin, (int) (mynearbyint(maxdist / table_spacing)), maxdist, table_spacing);
  if (distbin != (int) (mynearbyint(maxdist / table_spacing))) return 1;
  return 0;
}

/*********************************************************************
 * FUNCTION interp_lin
 *
 * Perform a linear interpolation to fill in energy/force tables
 * This should be replaced in the near future with a better interpolation
 *
 * Input:
 *  val1,val2 --  Y Values at the endpoints of the segments we interpolate on
 *  end1,end2 --  X coordinates of the corresponding endpoints
 *  currdist -- Distance we want the value at
 *  ** It is assumed that end2 > end1 **
 *
 * Output: Returns a floating point value at the requested point
 * ********************************************************************/

BigReal Parameters::interp_lin(BigReal val1, BigReal val2, BigReal end1, BigReal end2, BigReal currdist) {

  BigReal m; //slope of line
  BigReal val; // Value at desired point
  
  m = (val2 - val1) / (end2 - end1);

  val = ((currdist-end1) * m + val1);
  return val;
}

/*************************************************************************
 * FUNCTION get_int_table_type
 *
 * Find and return the integer index of a table type given its name
 *
 * Input:
 *  tabletype -- char array containing the name of the type to be looked up
 *
 * Output:
 * Returns an integer index < the total number of types, or -1 if the type could
 * not be found
 * ************************************************************************/

int Parameters::get_int_table_type(char* tabletype) {
  for (int i=0; i<tablenumtypes; i++) {
    if (!strncmp(tabletype, table_types[i], 5)) {
      return i;
    }
  }

  return -1;
}


