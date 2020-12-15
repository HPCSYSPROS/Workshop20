/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   PDB Class
     Code to implement the PDB class.  This reads in a bunch of
   PDBData records from a file, given the filename.  See PDB.h
   for a bit more information.
*/

#include <stdio.h>
#include <stdlib.h>
#ifndef WIN32
#include <strings.h>
#endif
#include "common.h"
#include "PDB.h"
#include "SortableResizeArray.h"


#ifdef MEM_OPT_VERSION
//read in PDB from a file and eliminate the temporary memory usage of pdb atom list    
PDB::PDB(const char *pdbfilename, int expectedNumAtoms){
  FILE *infile;
  char buf[160];

  altlocArray = 0;
  atomArray = new PDBCoreData[expectedNumAtoms];

  atomCount = 0;
  atomListHead = atomListTail = NULL;
  infile = Fopen(pdbfilename, "r");
  if (! infile) {
     char s[500];
     sprintf(s, "Cannot open file '%s' for input in PDB::PDB.", pdbfilename);
     NAMD_err(s);
  }
  
    // loop through each line in the file and get a record
  while ( fgets(buf, 150, infile) ) {
   PDBData *newelement;
   char *s;
   for (s=buf; *s && *s!='\n'; s++)  // chop off the '\n'
    ;
   *s = 0;
   if ( s == (buf + 149) ) {
     char s[500];
     sprintf( s, "Line too long in pdbfile %s:\n%s\n", pdbfilename, buf);
     NAMD_die(s);
   }
   *(s+1) = 0;  // just to be on the safe side

     // I now have a string; make a PDBData element out of it
   newelement = new_PDBData(buf);
   if (!newelement) {
      NAMD_die("Could not allocate PDBData.\n");
   }
     // I only know how to deal with ATOM and HETATM types; and
     //  I want to throw away the unknown data types; so
   if (newelement -> type() != PDBData::ATOM && 
           newelement -> type() != PDBData::HETATM) {
       delete newelement;
   } else {
       PDBAtom *thisAtom = (PDBAtom *)newelement;
       const BigReal *coor = thisAtom->coordinates();
       atomArray[atomCount].coor[0] = coor[0];
       atomArray[atomCount].coor[1] = coor[1];
       atomArray[atomCount].coor[2] = coor[2];
       atomArray[atomCount].myoccupancy = thisAtom->occupancy();
       atomArray[atomCount].tempfactor = thisAtom->temperaturefactor();
       atomCount++;
       delete newelement;
       if(atomCount > expectedNumAtoms) NAMD_die("Number of pdb and psf atoms are not the same!");       
   }
  }  // input while loop
 Fclose(infile);
}
#endif

PDB::PDB(molfile_plugin_t *pIOHdl, void *pIOFileHdl, int numAtoms, const float *occupancy, const float *bfactor){
#ifdef MEM_OPT_VERSION
  NAMD_die("Sorry, plugin IO is not supported in the memory optimized version.");
#else
  molfile_timestep_t ts;
  float *atomcoords;
  memset(&ts, 0, sizeof(molfile_timestep_t));

  /* set defaults for unit cell information */
  ts.A = ts.B = ts.C = 0.0f;
  ts.alpha = ts.beta = ts.gamma = 90.0f; 

  atomCount = numAtoms;

  atomcoords = (float *) malloc(3*numAtoms*sizeof(float));
  memset(atomcoords, 0, 3*numAtoms*sizeof(float));
  ts.coords = atomcoords;
  
  if (pIOHdl->read_next_timestep(pIOFileHdl, numAtoms, &ts)) {
    free(atomcoords);    
    pIOHdl->close_file_read(pIOFileHdl);
    NAMD_die("ERROR: failed reading atom coordinates");
  }

  //load coordinates to PDB object
  //Note: the PDBAtom structure is very redundant in this case
  PDBAtom *tmpAtoms = new PDBAtom[numAtoms];
  atomArray = new PDBAtomPtr[numAtoms];
  BigReal tmpCoords[3];
  for(int i=0; i<numAtoms; i++) {
      PDBAtom *thisAtom = tmpAtoms+i;
      atomArray[i] = thisAtom;

      tmpCoords[0] = atomcoords[0];
      tmpCoords[1] = atomcoords[1];
      tmpCoords[2] = atomcoords[2];
      atomcoords += 3; //set to the next atom
      thisAtom->coordinates(tmpCoords);
  }
  free(ts.coords);

  if(occupancy) {
      for(int i=0; i<numAtoms; i++) {
          PDBAtom *thisAtom = tmpAtoms+i;
          thisAtom->occupancy((BigReal)occupancy[i]);          
      }
  }
  if(bfactor) {
      for(int i=0; i<numAtoms; i++) {
          PDBAtom *thisAtom = tmpAtoms+i;       
          thisAtom->temperaturefactor((BigReal)bfactor[i]);       
      }
  }
#endif
}


// read in a file and stick all the elements on the appropriate list
PDB::PDB( const char *pdbfilename) {
  FILE *infile;
  char buf[160];

  atomCount = 0;
  atomListHead = atomListTail = NULL;
  infile = Fopen(pdbfilename, "r");
  if (! infile) {
     char s[500];
     sprintf(s, "Cannot open file '%s' for input in PDB::PDB.", pdbfilename);
     NAMD_err(s);
  }
  
    // loop through each line in the file and get a record
  while ( fgets(buf, 150, infile) ) {
   PDBData *newelement;
   char *s;
   for (s=buf; *s && *s!='\n'; s++)  // chop off the '\n'
    ;
   *s = 0;
   if ( s == (buf + 149) ) {
     char s[500];
     sprintf( s, "Line too long in pdbfile %s:\n%s\n", pdbfilename, buf);
     NAMD_die(s);
   }
   *(s+1) = 0;  // just to be on the safe side

     // I now have a string; make a PDBData element out of it
   newelement = new_PDBData(buf);
   if (!newelement) {
      NAMD_die("Could not allocate PDBData.\n");
   }
     // I only know how to deal with ATOM and HETATM types; and
     //  I want to throw away the unknown data types; so
   if (newelement -> type() != PDBData::ATOM && 
           newelement -> type() != PDBData::HETATM) {
       delete newelement;
   } else {
       add_atom_element( (PDBAtom *) newelement);
   }
  }  // input while loop
 Fclose(infile);
 
 // now I have a linked list, and I know the size.  However,
 // that's got pretty slow access time for ramdom fetches, so
 // I'll turn it into an array
 {

#ifdef MEM_OPT_VERSION
//pruning the PDBData structure into PDBCoreData by only leaving X/Y/Z, occupancy and temperaturefactor
     altlocArray = new char[atomCount];
     atomArray = new PDBCoreData[atomCount];
     if(atomArray == NULL)
         NAMD_die("memory allocation failed in PDB::PDB");

     PDBAtomList *tmp = atomListHead;
     for(int i=0; tmp!=NULL; tmp=tmp->next, i++){
         const BigReal *coor = tmp->data->coordinates();
         atomArray[i].coor[0] = coor[0];
         atomArray[i].coor[1] = coor[1];
         atomArray[i].coor[2] = coor[2];
         atomArray[i].myoccupancy = tmp->data->occupancy();
         atomArray[i].tempfactor = tmp->data->temperaturefactor();
         altlocArray[i] = tmp->data->alternatelocation()[0];
     }

     //free the PDBAtomList
     tmp = atomListHead;
     while(tmp!=NULL){
         PDBAtomList *nextTmp = tmp->next;
         delete tmp->data;
         delete tmp;
         tmp = nextTmp;
     }
     atomListHead = atomListTail = NULL;
#else
  atomArray = new PDBAtomPtr[atomCount];
  if ( atomArray == NULL )
  {
    NAMD_die("memory allocation failed in PDB::PDB");
  }
  PDBAtomList *tmp = atomListHead;
  int i=0;                              // just need to copy the pointers
  for (i=0, tmp = atomListHead; tmp != NULL; tmp = tmp -> next, i++) {
    atomArray[i] = tmp -> data;
  }
  // now delete the linked list (w/o deleting the data)
  PDBAtomList *tmp2;
  for (tmp2 = tmp = atomListHead; tmp != NULL; tmp = tmp2) {
    tmp2 = tmp->next;
    delete tmp;
  }
  atomListHead = atomListTail = NULL;
#endif 
 }  // everything converted
 
}

//  Destructor - delete all the data pointed to by the array
//   and then delete the array
PDB::~PDB( void )
{
#ifndef MEM_OPT_VERSION
	int i;
	if ( atomArray[atomCount-1] == atomArray[0] + (atomCount-1) ) {
	  delete [] atomArray[0];
	} else {
	  for (i=atomCount-1; i>=0; i--)
	    delete atomArray[i];
        }
#else
	delete [] altlocArray;
#endif
	delete [] atomArray;
	atomArray = NULL;
	atomCount = 0;
}

// print the PDB file out to a given file name
void PDB::write(const char *outfilename, const char *commentline)
{
	int i;
	char s[200];
	FILE *outfile;
	if ((outfile = fopen(outfilename, "w")) == NULL) {
	   sprintf(s, "Cannot open file '%s' in PDB::write.", outfilename);
	   NAMD_err(s);
	}

	if (commentline != NULL)
	{
		sprintf(s, "REMARK  %s\n", commentline);
		if (fputs(s, outfile) == EOF)
		{
			NAMD_err("EOF in PDB::write writing the comment line - file system full?");
		}
	}

	for (i=0; i<atomCount; i++){ // I only contain ATOM/HETATM records
#ifdef MEM_OPT_VERSION    
      atomArray[i].sprint(s, PDBData::COLUMNS);  
#else
	  atomArray[i]->sprint(s, PDBData::COLUMNS);
#endif
	  if ( (fputs(s, outfile)    == EOF) || 
	       (fputs("\n", outfile) == EOF)    ) {
	    sprintf(s, "EOF in PDB::write line %d - file system full?", i);
	    NAMD_err(s);
	  }
	}
	if (fputs("END\n", outfile) == EOF) {
	   NAMD_err("EOF in PDB::write while printing 'END' -- file system full?");
	}
	if (fclose(outfile) == EOF) {
	   NAMD_err("EOF in PDB::write while closing -- file system full?");
	}
	  
}

// store the info on the linked list
void PDB::add_atom_element( PDBAtom *newAtom)
{
  PDBAtomList *tmp = new PDBAtomList;
  if ( tmp == NULL )
  {
    NAMD_die("memory allocation failed in PDB::add_atom_element");
  }
  tmp -> data = newAtom;
  tmp -> next = NULL;
  
  if (atomListHead == NULL) {        // make the list
    atomListHead = atomListTail = tmp;
  } else {
    atomListTail -> next = tmp;       // add to the tail
    atomListTail = tmp;
  }
  atomCount++;
}


// return the number of atoms found
int PDB::num_atoms( void)
{
  return atomCount;
}


// Reset all the atom positions.  This is used in preparation for
// output in cases like the restart files, etc.
void PDB::set_all_positions(Vector *pos)
{
#ifdef MEM_OPT_VERSION
    for(int i=0; i<atomCount; i++){
        atomArray[i].coor[0] = pos[i].x;
        atomArray[i].coor[1] = pos[i].y;
        atomArray[i].coor[2] = pos[i].z;
    }
#else
	int i;
	PDBAtomPtr *atomptr;

	for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++)
	{
		(*atomptr)->xcoor(pos[i].x);
		(*atomptr)->ycoor(pos[i].y);
		(*atomptr)->zcoor(pos[i].z);
	}
#endif
}

void PDB::get_position_for_atom(Vector *pos, int aid){
#ifdef MEM_OPT_VERSION
    pos->x = atomArray[aid].coor[0];
    pos->y = atomArray[aid].coor[1];
    pos->z = atomArray[aid].coor[2];
#else
    PDBAtomPtr *atomptr = &atomArray[aid];
    pos->x = (*atomptr)->xcoor();
    pos->y = (*atomptr)->ycoor();
    pos->z = (*atomptr)->zcoor();
#endif
}

//  Get all the atom positions into a list of Vectors
void PDB::get_all_positions(Vector *pos)
{
#ifdef MEM_OPT_VERSION
    for(int i=0; i<atomCount; i++){
        pos[i].x = atomArray[i].coor[0];
        pos[i].y = atomArray[i].coor[1];
        pos[i].z = atomArray[i].coor[2];        
    }
#else
	int i;
	PDBAtomPtr *atomptr;

	for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++)
	{
		pos[i].x = (*atomptr)->xcoor();
		pos[i].y = (*atomptr)->ycoor();
		pos[i].z = (*atomptr)->zcoor();
	}
#endif
}

//  given an index, return that atom
#ifdef MEM_OPT_VERSION
PDBCoreData *PDB::atom(int place){
    if(place<0 || place>=atomCount) return NULL;
    return &atomArray[place];
}
#else
PDBAtom *PDB::atom(int place)
{
  if (place <0 || place >= atomCount)
    return NULL;
  return atomArray[place];
}
#endif


// find the lowest and highest bounds based on a fraction of the atoms
//void PDB::find_extremes(BigReal *min, BigReal *max, Vector rec, BigReal frac) const
void PDB::find_extremes_helper(
    SortableResizeArray<BigReal> &coor,
    BigReal &min, BigReal &max, Vector rec, BigReal frac
    )
{
    SortableResizeArray<BigReal>::iterator c_i = coor.begin();
#ifdef MEM_OPT_VERSION
    PDBCoreData *atomptr = atomArray;
    for(int i=0; i<atomCount; i++, atomptr++){
        Vector pos(atomptr->xcoor(),atomptr->ycoor(),atomptr->zcoor());
        c_i[i] = rec*pos;
    }
#else
    PDBAtomPtr *atomptr = atomArray;
    for (int i=0; i<atomCount; ++i, ++atomptr) {
      PDBAtom *atom = *atomptr;
      Vector pos(atom->xcoor(),atom->ycoor(),atom->zcoor());
      c_i[i] = rec*pos;
    }
#endif
    coor.sort();
    int ilow = (int)((1.0 - frac) * atomCount);
    if ( ilow < 0 ) ilow = 0;
    if ( ilow > atomCount/2 ) ilow = atomCount/2;
    int ihigh = atomCount - ilow - 1;
    BigReal span = coor[ihigh] - coor[ilow];
    BigReal extension = (1.0 - frac) * span / (2.0 * frac - 1.0);
    max = coor[ihigh] + extension;
    min = coor[ilow] - extension;
}

// Find the extreme edges of molecule in scaled coordinates,
// where "frac" sets bounds based on a fraction of the atoms.
void PDB::find_extremes(const Lattice &lattice, BigReal frac) {
  if (atomCount == 0) {
    smin = smax = 0;
  }
  else if (frac < 1.0) {
    // for now use the previous "sort the array" approach
    // for solving "frac"th largest and smallest selection problems
    SortableResizeArray<BigReal> coor;
    coor.resize(atomCount);  // allocate array space once
    find_extremes_helper(coor, smin.x, smax.x, lattice.a_r(), frac);
    find_extremes_helper(coor, smin.y, smax.y, lattice.b_r(), frac);
    find_extremes_helper(coor, smin.z, smax.z, lattice.c_r(), frac);
  }
  else {
    // finding absolute min and max does not require sorting
#ifdef MEM_OPT_VERSION
    PDBCoreData *atomptr = atomArray;
    Vector p(atomptr->xcoor(),atomptr->ycoor(),atomptr->zcoor());
#else
    PDBAtomPtr *atomptr = atomArray;
    PDBAtom *atom = *atomptr;
    Vector p(atom->xcoor(),atom->ycoor(),atom->zcoor());
#endif
    Vector s(lattice.a_r()*p, lattice.b_r()*p, lattice.c_r()*p);
    smin = smax = s;
    atomptr++;
    for(int i=1; i<atomCount; i++, atomptr++){
#ifdef MEM_OPT_VERSION
      p = Vector(atomptr->xcoor(),atomptr->ycoor(),atomptr->zcoor());
#else
      atom = *atomptr;
      p = Vector (atom->xcoor(),atom->ycoor(),atom->zcoor());
#endif
      s = Vector(lattice.a_r()*p, lattice.b_r()*p, lattice.c_r()*p);
      if      (smin.x > s.x) smin.x = s.x;
      else if (smax.x < s.x) smax.x = s.x;
      if      (smin.y > s.y) smin.y = s.y;
      else if (smax.y < s.y) smax.y = s.y;
      if      (smin.z > s.z) smin.z = s.z;
      else if (smax.z < s.z) smax.z = s.z;
    }
  }
  // shift using the origin
  BigReal origin_shift;
  origin_shift = lattice.a_r() * lattice.origin();
  smin.x -= origin_shift;
  smax.x -= origin_shift;
  origin_shift = lattice.b_r() * lattice.origin();
  smin.y -= origin_shift;
  smax.y -= origin_shift;
  origin_shift = lattice.c_r() * lattice.origin();
  smin.z -= origin_shift;
  smax.z -= origin_shift;
}

//#define TEST_PDB_CLASS
#ifdef TEST_PDB_CLASS

main()
{
 PDB *pdb = new PDB("pti.pdb");
 if ( atomArray == NULL )
 {
   NAMD_die("memory allocation failed in main of test PDB class");
 }
 ilist = pdb->find_atom_name("CA");
 if (!ilist)
   printf("None found.\n");
 else {
   int i;
   char s[200];
   for (i=0; i<ilist->num(); i++) {
     pdb->atom((*ilist)[i]) -> sprint(s);
     printf("%s\n", s);
   }
   delete ilist;
 } // if ilist
 
 printf("Now, search through space.\n");
 
 ilist = pdb->find_atoms_in_region(4.38, 19.5, 3.0,
                                   4.40, 20.0, 3.2 );
 if (!ilist)
   printf("None found.\n");
 else {
   int i;
   char s[200];
   printf("%d found\n", ilist -> num());
   for (i=0; i<ilist->num(); i++) {
     pdb->atom((*ilist)[i]) -> sprint(s);
     printf("%s\n", s);
   }
   delete ilist;
 }
}
#endif // TEST_PDB_CLASS



// This function was borrowed from VMD code in "ReadPARM.C".
// It skips to a new line.
static int readtoeoln(FILE *f) {
  char c;

  /* skip to eoln */
  while((c = getc(f)) != '\n') {
    if (c == EOF) 
      return -1;
  }

  return 0;
}  



// read in an AMBER coordinate file and populate the PDB structure
PDB::PDB( const char *filename, Ambertoppar *amber_data)
{ int i,j,k;
  BigReal coor[3];
  char buf[13],resname[5],atomname[5];
  FILE *infile;
  PDBAtom *pdb;

  if ((infile=Fopen(filename, "r")) == NULL)
    NAMD_err("Can't open AMBER coordinate file!");

  readtoeoln(infile);  // Skip the first line (title)

  fscanf(infile,"%d",&atomCount);  // Read in num of atoms
  if (atomCount != amber_data->Natom)
    NAMD_die("Num of atoms in coordinate file is different from that in parm file!");
  readtoeoln(infile);

#ifdef MEM_OPT_VERSION
  altlocArray = 0;
  atomArray = new PDBCoreData[atomCount];
#else
  atomArray = new PDBAtomPtr[atomCount];
#endif

  if ( atomArray == NULL )
  {
    NAMD_die("memory allocation failed in PDB::PDB");
  }
  
  // Read in the coordinates, which are in the format of 6F12.7
  // Other fields are copied from "amber_data"
  for (i=0; i<atomCount; ++i)
  { // Read x,y,z coordinates
    for (j=0; j<3; ++j)
    { for (k=0; k<12; ++k)
      { buf[k]=getc(infile);
        if (buf[k]=='\n' || buf[k]=='\0' || buf[k]==EOF)
          NAMD_die("Error reading AMBER coordinate file!");
      }
      buf[12] = '\0';
      coor[j] = atof(buf);
    }
    if (i%2 == 1)
      readtoeoln(infile);
#ifdef MEM_OPT_VERSION
    atomArray[i].coor[0] = coor[0];
    atomArray[i].coor[1] = coor[1];
    atomArray[i].coor[2] = coor[2];
    atomArray[i].myoccupancy = PDBAtom::default_occupancy;
    atomArray[i].tempfactor = PDBAtom::default_temperaturefactor;
#else
    // Copy name, resname and resid from "amber_data"
    for (j=0; j<4; ++j)
    { resname[j] = amber_data->ResNames[amber_data->AtomRes[i]*4+j];
      atomname[j] = amber_data->AtomNames[i*4+j];
    }
    resname[4] = atomname[4] = '\0';
    // Create a new PDB record, and fill in its entries
    pdb = new PDBAtomRecord("");
    pdb->name(atomname);
    pdb->residuename(resname);
    pdb->serialnumber(i+1);
    pdb->residueseq(amber_data->AtomRes[i]+1);
    pdb->coordinates(coor);
    atomArray[i] = pdb;  // Include the new record into the array
#endif
  }
}

#define LINESIZE 100

/* This constructor initializes the PDB data using a Gromacs
   coordinate file, generating an error message if the file
   can't be parsed or if its contents don't jive with what is in
   the topo file <topology>. */
PDB::PDB(const char *filename, const GromacsTopFile *topology) {
  int i;
  char buf[LINESIZE];
  FILE *infile;
  
  /* open up the coordinate file */
  infile=Fopen(filename, "r");
  if (infile == NULL)
    NAMD_err("Can't open GROMACS coordinate file!");

  fgets(buf,LINESIZE-1,infile); // get the title
  // if(strcmp(buf,topology->getSystemName()) != 0)
  //   NAMD_die("System names in topology and coordinate files differ.");

  fgets(buf,LINESIZE-1,infile); // get the number of atoms
  sscanf(buf,"%d",&atomCount);
  if (atomCount != topology->getNumAtoms())
    NAMD_die("Num of atoms in coordinate file is different from that in topology file!");

  /* read in the atoms */
#ifdef MEM_OPT_VERSION
  altlocArray = 0;
  atomArray = new PDBCoreData[atomCount];
#else
  atomArray = new PDBAtomPtr[atomCount];
#endif
  if ( atomArray == NULL )
    NAMD_die("memory allocation failed in PDB::PDB");

#ifdef MEM_OPT_VERSION
  for(i=0; i<atomCount; i++){
      fgets(buf,LINESIZE-1,infile); // get a line
      char *buf2 = buf+20; // skip three fields to get to the coordinates
      BigReal coor[3];
      if(3 != sscanf(buf2,"%lf%lf%lf", &coor[0],&coor[1],&coor[2]))
          NAMD_die("Couldn't get three coordinates from file.");

      coor[0] *= 10; // convert to angstroms from nanometers
      coor[1] *= 10;
      coor[2] *= 10;

      atomArray[i].coor[0] = coor[0];
      atomArray[i].coor[1] = coor[1];
      atomArray[i].coor[2] = coor[2];
      atomArray[i].myoccupancy = PDBAtom::default_occupancy;
      atomArray[i].tempfactor = PDBAtom::default_temperaturefactor;
  }
#else
  for (i=0;i<atomCount;i++) {
    char *buf2, resname[11], atomname[11], atmtype[11];
    int resnum, typenum;
    Real charge,mass;
    BigReal coor[3];
    PDBAtom *pdb = new PDBAtomRecord("");  
    
    fgets(buf,LINESIZE-1,infile); // get a line
    buf2 = buf+20; // skip three fields to get to the coordinates
    if(3 != sscanf(buf2,"%lf%lf%lf",
		   &coor[0],&coor[1],&coor[2]))
      NAMD_die("Couldn't get three coordinates from file.");
    topology->getAtom(i,&resnum,resname,
		      atomname,atmtype,&typenum,&charge,&mass);
    coor[0] *= 10; // convert to angstroms from nanometers
    coor[1] *= 10;
    coor[2] *= 10;
    
    pdb->name(atomname);
    pdb->residuename(resname);
    pdb->serialnumber(i+1);
    pdb->residueseq(resnum+1);
    pdb->coordinates(coor);
    
    atomArray[i] = pdb;  // Include the new record into the array
  }
#endif
}

#ifdef MEM_OPT_VERSION
void PDBCoreData::sprint( char *outstr, PDBData::PDBFormatStyle usestyle){
    if(usestyle == PDBData::COLUMNS){
        //sprint_columns copied from PDBAtom::sprint_columns
        for(int i=0; i<79; i++) outstr[i] = 32;
        outstr[79] = 0;
        PDBData::sprintcol(outstr, PDBAtom::SX, PDBAtom::LCOOR, PDBAtom::LCOORPREC, xcoor());
        PDBData::sprintcol(outstr, PDBAtom::SY, PDBAtom::LCOOR, PDBAtom::LCOORPREC, ycoor());
        PDBData::sprintcol(outstr, PDBAtom::SZ, PDBAtom::LCOOR, PDBAtom::LCOORPREC, zcoor());
        PDBData::sprintcol(outstr, PDBAtom::SOCC, PDBAtom::LOCC, PDBAtom::LOCCPREC, occupancy());
        PDBData::sprintcol(outstr, PDBAtom::STEMPF, PDBAtom::LTEMPF, PDBAtom::LTEMPFPREC, temperaturefactor());
    }else{
        //sprint_fields
        char tmpstr[50];        
        if (xcoor() == PDBAtom::default_coor)
            sprintf(tmpstr, " #");
           else
            sprintf(tmpstr, " %*.*f", PDBAtom::LCOOR, PDBAtom::LCOORPREC, xcoor());
        strcat(outstr, tmpstr);
        if (ycoor() == PDBAtom::default_coor)
            sprintf(tmpstr, " #");
           else
            sprintf(tmpstr, " %*.*f", PDBAtom::LCOOR, PDBAtom::LCOORPREC, ycoor());
        strcat(outstr, tmpstr);
        if (zcoor() == PDBAtom::default_coor)
            sprintf(tmpstr, " #");
           else
            sprintf(tmpstr, " %*.*f", PDBAtom::LCOOR, PDBAtom::LCOORPREC, zcoor());
        strcat(outstr, tmpstr);
      
        sprintf(tmpstr, " %*.*f",  PDBAtom::LOCC, PDBAtom::LOCCPREC, occupancy());
        strcat(outstr, tmpstr);
      
        sprintf(tmpstr, " %*.*f", PDBAtom::LTEMPF, PDBAtom::LTEMPFPREC, temperaturefactor());
        strcat(outstr, tmpstr);        
    }
}
#endif

