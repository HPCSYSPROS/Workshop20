/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Read and parse a line of data from a PDB record.  There are many
   different types of PDB records.  This version reads only the ATOM and
   HETATM records and makes all fields accessible via the appropriate
   member function.  In NAMD, this will be called only by the PDB class,
   which reads PDB files, and a PDB writer class.
*/

#ifndef _PDBREADER_H_
#define _PDBREADER_H_

// These are added to the global namespace
//   the class PDBData
//   the class PDBUnknown, derived from PDBData
//   the class PDBAtom, derived from PDBData - contains ATOM and HETATM records
//   the class PDBAtomRecord, derived from PDBAtom - contains only ATOM records
//   the class PDBHetAtm, derived from PDBAtom - contains only HETATM records
//   the function new_PDBData, which creates the right pdb class given a string

#ifndef WIN32
#include <strings.h>
#endif
#include "common.h"
#include <string.h>


class PDBData {  // at the basic level, a PDB record only knows its type
   public:

  // These data types come from the Protein Data Bank format
  // description located via anon. ftp to pdb.pdb.bnl.gov
  // in the file /pub/format.desc.ps
  //  In addition, I define an extra type, UNKNOWN.  If I can't
  // figure out what's going on, I store the complete string
  // and return it when asked.
   enum PDBType {HEADER, OBSLTE, COMPND, SOURCE, EXPDTA, AUTHOR,
     REVDAT, SPRSDE, JRNL, REMARK, SEQRES, FTNOTE, HET, FORMUL,
     HELIX, SHEET, TURN, SSBOND, SITE, CRYST1, ORIGX, SCALE,
     MTRIX, TVECT, MODEL, ATOM, HETATM, SIGATM, ANISOU, SIGUIJ,
     TER, ENDMDL, CONECT, MASTER, END, UNKNOWN};

     enum PDBFormatStyle { COLUMNS, FIELDS };  // used to specify if the



    static const char *PDBNames[UNKNOWN+1]; // string descriptors for each field
    
       // output should be based on columns (FORTRAN style) or
       // fields (C/ awk style).
 // See, there are two different types of formats that this program
 // understands, one is the basic PDB w/ or w/o the XLPOR extension - these
 // are the column based versions.  The other is my own
 // field based version - each data element is seperated by a blank
 // and, if the element is empty, a pound sign ('#') is put in its place.
 // This type of record is denoted by a '#' in the first non-blank
 // character (hence, it is the first non-blank character of the first
 // field.  Basically, I'm a unix/ C/ awk/ yacc ... freak, I like field
 // based data rather than column based data.
 
 
   private:
    PDBType mytype;

#ifdef MEM_OPT_VERSION
//for the sake of pruning PDBData
public:
#else
   protected:
#endif
        // some parsing routines to get info from a line of text
    static void scan( const char *data, int len, int start, int size, 
                         int *ans, int defalt);
    static void scan( const char *data, int len, int start, int size,
                          BigReal *ans, BigReal defalt);
    static void scan( const char *data, int len, int start, int size,
                         char *ans);
    static void field( const char *data, int fld, char *result);
        // some routine to print to a specific column and width
    static void sprintcol( char *s, int start, int len, const char *val);
    static void sprintcol( char *s, int start, int len, int val);
    static void sprintcol( char *s, int start, int len, int prec, BigReal val);
    
   public:
     PDBData(PDBType newtype) {
       mytype = newtype;
     }
     virtual ~PDBData( void) {
     }
     PDBType type( void) {
       return mytype;
     }
                    // I know nothing, so I'll fake it and hope it works
     virtual void sprint( char *s, PDBFormatStyle usestyle = COLUMNS) {
       if (usestyle == COLUMNS)     // get rid of warning
         strcpy(s, "REMARK     (undefined remark - this is a bug)");
        else
         strcpy(s, "REMARK     (undefined remark - this is a bug)");
     }
};

//////******* the "UNKNOWN" class *****//////
class PDBUnknown : public PDBData {
  private:
    char *mystr;
  public:
   PDBUnknown(const char *data): PDBData(PDBData::UNKNOWN) {
     mystr = new char[strlen(data)+1];
     if ( mystr == NULL )
     {
       NAMD_die("memory allocation failed in PDBUnknown::PDBUnknown");
     }
     strcpy(mystr, data);
   }
   virtual ~PDBUnknown( void) {
     delete [] mystr;
   }  
   void sprint(char *s, PDBFormatStyle usestyle) {
     strcpy(s, mystr);
     if (usestyle == PDBData::COLUMNS)   // they are the same, but I won't
       strcpy( s, mystr);                //   get the stupid warning during
      else                               //   compilation
       strcpy( s, mystr);
   }
};


////************* routines used for ATOM and HETATM **********/////
class PDBAtom : public PDBData {
public:
    //extract them out from PDBAtom for the sake of pruning PDBData
      // starting location for each record element
    enum Start {STYPE=1,SSERIAL=7, SNAME=13, SALT=17, SRESNAME=18, SCHAIN=22, 
                SRESSEQ=23, SINSERT=27, SX=31, SY=39, SZ=47,
                SOCC=55, STEMPF=61, SFOOT=68, SSEGNAME=73, SELEMENT=77};
      // length of each element, the PREC is the number of digits
      // in the output after the decimal
// NOTE: The PDB says the length of the residue name is only 3 characters
//  whereas XPLOR allows 4 character names.  We choose 4 for compatability
//  with both systems (since we never change the length, we you give us is
//  what we use)
    enum Length {LTYPE=6, LSERIAL=5, LNAME=4, LALT=1, LRESNAME=4, LCHAIN=1, 
                 LRESSEQ=4, LINSERT=1, LCOOR=8,
                 LCOORPREC=3, LOCC=6, LOCCPREC=2, LTEMPF=6, 
                 LTEMPFPREC=2, LFOOT=3, LSEGNAME=4, LELEMENT=2};

  public:
    static const int default_serial;         // some default values
    static const int default_residueseq;     // these are set in the .C file
    static const BigReal default_coor;
    static const BigReal default_occupancy;
    static const BigReal default_temperaturefactor;
    static const int no_footnote;

  private:
    int myserialnumber;                 // atom serial number
    char myname[LNAME+1];               // atom name
    char myalternatelocation[LALT+1];   // alternamte location identifier
    char myresiduename[LNAME+1];        // residue name
    char mychain[LCHAIN+1];             // chain indentifier
    int myresidueseq;                   // residue seq. no.
    char myinsertioncode[LINSERT+1];    // code for insertions of residues
    BigReal mycoor[3];                     // X, Y, and Z orthogonal A coordinates
    BigReal myoccupancy;                   // occupancy
    BigReal mytemperaturefactor;           // temperature factor
    int myfootnote;                     // footnote number
    char mysegmentname[LSEGNAME+1];     // XPLOR-type segment name
    char myelement[LELEMENT+1];         // element

    void parse_field_data( const char *data);
    void parse_column_data( const char *data);
    void sprint_columns( char *outstr);
    void sprint_fields( char *outstr);

  protected:
    enum PDBPossibleAtoms {USE_ATOM = ATOM, USE_HETATM = HETATM};
    PDBAtom( const char *data,
           PDBPossibleAtoms whichatom);// parses a line from the PDB data file
    //PDBAtom( void);        // makes a generic atom

  public:
    PDBAtom( void);        // makes a generic atom
    virtual ~PDBAtom( void);
    void parse( const char *s);  // reset to new input values
    void  sprint( char *s, PDBFormatStyle usestyle = COLUMNS);// write to string
    int serialnumber( void);
    void serialnumber( int newserialnumber);
    
    const char*name( void);
    void name( const char *newname);
    
    const char*alternatelocation( void);
    void alternatelocation( const char *newalternatelocation);
    
    const char*residuename( void);
    void residuename( const char *newresiduename);
    
    const char*chain( void);
    void chain( const char *newchain);
    
    int residueseq( void);
    void residueseq( int newresidueseq);
    
    const char*insertioncode( void);
    void insertioncode( const char *newinsertioncode);
    
    BigReal xcoor( void);
    void xcoor( BigReal newxcoor);
    BigReal ycoor( void);
    void ycoor( BigReal newycoor); 
    BigReal zcoor( void);
    void zcoor( BigReal newzcoor);
    
    const BigReal *coordinates( void);
    void coordinates(const BigReal *newcoordinates);
    
    BigReal occupancy( void);
    void occupancy( BigReal newoccupancy);

    BigReal temperaturefactor( void);
    void temperaturefactor( BigReal newtemperaturefactor);

    int footnote( void);
    void footnote( int newfootnote);
    
      // this is not part of the PDB format but is used by XPLOR instead of
      // the chain identifier (see XPLOR 3.1 manual, p 104)
    const char*segmentname( void);
    void segmentname( const char *newsegmentname);

    const char* element( void);
    void element( const char *newelement);
};

// The two sub-classes of PDB Atom
class PDBAtomRecord : public PDBAtom{
   public:
     PDBAtomRecord( const char *data ) :
          PDBAtom( data, PDBAtom::USE_ATOM) {
     }
     virtual ~PDBAtomRecord( void) {
     }
};

class PDBHetatm : public PDBAtom {
  public:
    PDBHetatm( const char *data) :
         PDBAtom( data, PDBAtom::USE_HETATM) {
    }
    virtual ~PDBHetatm( void) {
    }
};


#ifdef MEM_OPT_VERSION
struct PDBCoreData{
    BigReal coor[3];              // X, Y, and Z orthogonal A coordinates
    BigReal myoccupancy;          // occupancy
    BigReal tempfactor;           // temperature factor

    //These functions are added for fewer changes in the src code when
    //pruning the PDBData
    BigReal occupancy() { return myoccupancy;   }
    BigReal xcoor() { return coor[0];  }
    BigReal ycoor() { return coor[1];  }
    BigReal zcoor() { return coor[2];  }
    BigReal temperaturefactor() { return tempfactor; }

    void  sprint( char *s, PDBData::PDBFormatStyle usestyle = PDBData::COLUMNS);// write to string
};
#endif

////********* Wrap up everything in one function call ***********//////
// somehow I need the base class to figure out which derived class
// to use to parse.   Since I don't know how to do that, I'll
// fake it with this.  Give it a string and it will create the
// correct PDB data type.
PDBData *new_PDBData(const char *data);  // nasty


#endif

