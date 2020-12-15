/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   PDB Class
     Given a PDB file name, read in all the data.
*/

#ifndef PDB_H
#define PDB_H

// These are added to the global namespace:
//   whatever is in PDBData.h
//   the typedef PDBAtomList (a singly linked list of PDBAtom *
//   the class PDB

#include "parm.h"
#include "ResizeArray.h"
#include "GromacsTopFile.h"

#include "PDBData.h"
#include "Vector.h"
#include "Lattice.h"
#include "molfile_plugin.h"

typedef PDBAtom *PDBAtomPtr ;
typedef struct PAL {
  PDBAtom *data;
  struct PAL *next;
} PDBAtomList;
  
class PDB {
  private:
    PDBAtomList *atomListHead, *atomListTail;

#ifdef MEM_OPT_VERSION
    char *altlocArray;
    PDBCoreData *atomArray;
#else
    PDBAtom **atomArray;
#endif
      // this doesn't create a copy 
    void add_atom_element(PDBAtom *newAtom); 
    int atomCount;
    
    ScaledPosition smin, smax;  // extreme edges of the molecular system

    void find_extremes_helper(
        SortableResizeArray<BigReal> &coor,
        BigReal &min, BigReal &max, Vector rec, BigReal frac
        );

  public:
    PDB(const char *pdbfilename);   // read in PDB from a file

#ifdef MEM_OPT_VERSION
    //read in PDB from a file and eliminate the temporary memory usage of pdb atom list    
    PDB(const char *pdbfilename, int expectedNumAtoms); 
#endif

    //Constructor for plugin IO based way of loading atoms' structure
    PDB(molfile_plugin_t *pIOHdl, void *pIOFileHdl, int numAtoms, const float *occupancy, const float *bfactor);

    PDB(const char *, Ambertoppar *);  // read AMBER coordinate file

    /* This constructor initializes the PDB data using a Gromacs
       coordinate file, generating an error message if the file
       can't be parsed or if its contents don't jive with what is in
       the topo file <topology>. */
    PDB(const char *filename, const GromacsTopFile *topology);

    ~PDB(void);               // clear everything
    void write(const char *outfilename, const char *commentline=NULL); // write the coordinates to a file
       // the following deals only with ATOMs and HETATMs
    int num_atoms( void);

#ifdef MEM_OPT_VERSION
    PDBCoreData *atom(int place);           
    char alternatelocation(int place) { return altlocArray[place]; }
    
    void delPDBCoreData() { delete [] atomArray; atomArray=NULL; } 
#else
    PDBAtom *atom(int place); // get the nth atom in the PDB file
#endif    
         // return linked list containing all atoms
    PDBAtomList *atoms(void ) { return atomListHead; }  

#if 0
	// Find the extreme edges of the molecule
    void find_extremes(BigReal *min, BigReal *max, Vector rec,
                                                  BigReal frac=1.0) const;
#else
    // Find the extreme edges of molecule in scaled coordinates,
    // where "frac" sets bounds based on a fraction of the atoms.
    void find_extremes(const Lattice &, BigReal frac=1.0);

    // Obtain results after find_extremes().
    void get_extremes(ScaledPosition &xmin, ScaledPosition &xmax) const {
      xmin = smin;  xmax = smax;
    }
#endif

    void set_all_positions(Vector *);	//  Reset all the positions in PDB

    void get_all_positions(Vector *);	//  Get all positions in PDB

    void get_position_for_atom(Vector *, int); //Get the position for an atom
};

#endif // PDB_H

