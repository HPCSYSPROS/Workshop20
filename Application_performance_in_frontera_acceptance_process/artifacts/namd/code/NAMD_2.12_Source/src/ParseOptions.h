/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Simplifies the interaction with the ConfigList.  This allows you to
   define a set of dependencies in the configuration file, along with
   default values, ranges, and units.  Then, given the ConfigList, this
   checks the list for valid entries and reads those values, warns for
   unneeded options, and errors for stuff it couldn't understand, or
   were out of range or in the wrong units, ...
*/

#ifndef PARSEOPTIONS_H

#define PARSEOPTIONS_H

#include "common.h"
#include "Vector.h"

class StringList;
class ConfigList;

enum Range { FREE_RANGE, POSITIVE, NOT_NEGATIVE , NEGATIVE, NOT_POSITIVE };
// const char *string(Range r); // for printing the range
enum Units { N_UNIT, N_FSEC, N_NSEC, N_SEC, N_MIN, N_HOUR, N_ANGSTROM, N_NANOMETER, N_METER, 
             N_KCAL, N_KJOULE, N_EV, N_KELVIN, N_UNITS_UNDEFINED};
// const char *string(Units u); // for printing the units
BigReal convert(Units to, Units from);//return 0 if stupid (like METER to SEC)

#define PARSE_FLOAT (BigReal *) NULL
#define PARSE_BIGREAL (BigReal *) NULL
#define PARSE_VECTOR (Vector *) NULL
#define PARSE_INT (int *) NULL
#define PARSE_BOOL (int *) NULL
#define PARSE_STRING (char *) NULL
#define PARSE_ANYTHING (StringList **) NULL
#define PARSE_MULTIPLES (StringList **) NULL, TRUE

class ParseOptions {
 public:
   
   class DataElement {
    public:
      enum data_types {UNDEF, FLOAT, VECTOR, UINT, INT, BOOL, STRINGLIST, STRING};
      int index;
      char *name;        // name of the variable
      int is_optional;
      int is_defined;
      char *parent;     // I need this for "forward" definitions
      DataElement *parent_ptr;
      char *error_message;
      data_types type;
      int has_default;
      int many_allowed;  // only used by StringList; otherwise assumed FALSE
      union {
	 BigReal fdef;     // the default data value
	 int idef;
	 unsigned int uidef;
	 // there is no default element for the StringList or String
      };
      Vector vdef;        // seperate since it has its own constructor
      union {
	 BigReal *fptr;     // the element to set
	 int *iptr;
	 unsigned int *uiptr;
	 Vector *vptr;
	 StringList **slptr; // for string lists
	 char *sptr;     // for strings
      };
      union {           // the actual data
	 BigReal fdata;
	 int idata;
	 int uidata;
	 StringList *sldata;  // this points directly to the ConfigList, so
      };                     // it had better not be deallocated.
      Vector vdata;
      Range range;
      Units units;
	 
    private:
      void init(const char *newname, const char *newparent, int optional,
		const char *err);
    public:
      ~DataElement(void);
      // FLOAT
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, BigReal *ptr, BigReal defalt);
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, BigReal *ptr);
      // VECTOR
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, Vector *ptr, Vector defalt);
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, Vector *ptr);
      
      // INT and BOOL
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, int *ptr, int defalt);
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, int *ptr);

      //  UNSIGNED INT
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, unsigned int *ptr, unsigned int defalt);
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, unsigned int *ptr);

      // STRINGLIST
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, StringList **ptr, int many_allowed = FALSE);
      // STRING
      DataElement(const char *newname, const char *newparent, int optional,
		  const char *err, char *ptr);
   };
 private:
   ConfigList const *configList;
   DataElement **data_array;
   int array_size;
   int array_max_size;
   void add_element(DataElement *el);
   int make_dependencies(DataElement *el);

   // indicates if the string is TRUE (1), FALSE (0), or unknown (-1)
   // if it can't figure it out
   int atoBool(const char *s);
   // convert from a string to Units, returns UNITS_UNDEFINED if error
   Units atoUnits(const char *s);

   // tells if this node has children (TRUE for yes)
   Bool is_parent_node(DataElement *el);

 public:
   ParseOptions(void);
   ~ParseOptions(void);
   // returns 1 if everything okay, or 0 if already defined
   // the msg is printed to namdErr if it was required and not
   // defined and no default exists
   int require(const char *newname, const char *parent, const char *msg,
	       BigReal *ptr, BigReal defalt);
   int require(const char *newname, const char *parent, const char *msg,
	       BigReal *ptr);
   
   int require(const char *newname, const char *parent, const char *msg,
	       Vector *ptr, Vector defalt);
   int require(const char *newname, const char *parent, const char *msg,
	       Vector *ptr);
   
   int require(const char *newname, const char *parent, const char *msg,
	       int *ptr, int defalt);
   int require(const char *newname, const char *parent, const char *msg,
	       int *ptr);

   int require(const char *newname, const char *parent, const char *msg,
	       unsigned int *ptr, unsigned int defalt);
   int require(const char *newname, const char *parent, const char *msg,
	       unsigned int *ptr);

   // a "boolean" version, which can convert a string to 1 if true
   // and 0 if false;  everywhere else this looks like an int
   // (now, if we have a real Boolean type, we wouldn't need this 'B')
   int requireB(const char *newname, const char *parent, const char *msg,
		int *ptr, int defalt);
   int requireB(const char *newname, const char *parent, const char *msg,
		int *ptr);
   // for the StringList; there is no default version
   int require(const char *newname, const char *parent, const char *msg,
	       StringList **ptr = NULL, int many_allowed = FALSE);
   // Strings don't have a default version, either
   int require(const char *newname, const char *parent, const char *msg,
	       char *ptr);
   
   int optional(const char *newname, const char *parent, const char *msg,
		BigReal *ptr, BigReal defalt);
   int optional(const char *newname, const char *parent, const char *msg,
		BigReal *ptr);

   int optional(const char *newname, const char *parent, const char *msg,
		Vector *ptr, Vector defalt);
   int optional(const char *newname, const char *parent, const char *msg,
		Vector *ptr);

   int optional(const char *newname, const char *parent, const char *msg,
		int *ptr, int defalt);
   int optional(const char *newname, const char *parent, const char *msg,
		int *ptr);

   int optional(const char *newname, const char *parent, const char *msg,
		unsigned int *ptr, unsigned int defalt);
   int optional(const char *newname, const char *parent, const char *msg,
		unsigned int *ptr);

   int optionalB(const char *newname, const char *parent, const char *msg,
		 int *ptr, int defalt);
   int optionalB(const char *newname, const char *parent, const char *msg,
		 int *ptr);

   int optional(const char *newname, const char *parent, const char *msg,
		StringList **ptr = NULL, int many_allowed = FALSE);
   int optional(const char *newname, const char *parent, const char *msg,
		char *ptr);
   

   // get the range of the given variable
   Range range(const char *name);
   // set the range of the given variable
   void range(const char *name, Range newrange);
 private:
   // find the children of the given element; used by check_consistency
   int check_children(int idx, int *flg);

   // read a string into the appropriate data type; do units if need be
   Bool scan_float(DataElement *el, const char *s);
   Bool scan_int(DataElement *el, const char *s);
   Bool scan_uint(DataElement *el, const char *s);
   Bool scan_bool(DataElement *el, const char *s);
   Bool scan_vector(DataElement *el, const char *s);

   // check if the range is correct and, if not NULL, set the variable
   // (at this point, the new value already exists in the data element)
   Bool set_float(DataElement *el);
   void set_vector(DataElement *el);
   Bool set_int(DataElement *el);
   Bool set_uint(DataElement *el);
   void set_bool(DataElement *el);
   void set_stringlist(DataElement *el);  //  also checks for 'many_allowed'
   void set_string(DataElement *el);
 public:
   // make sure the options were defined properly; return TRUE okay, FALSE not
   Bool check_consistency(void);
   // returns TRUE if everything set okay or FALSE if not
   Bool set(const ConfigList& configlist);

   // find the specified element.  Returns NULL if not found;
 private:
   DataElement *internal_find(const char *name);
 public:
   // special accessor for ScriptTcl
   char* getfromptr(const char* name, char *outbuf);
   int istruefromptr(const char* name);
   int issetfromptr(const char* name);
   // get the specified element.  Returns FALSE if not found or undefined;
   // prints warning if had to do a type conversion
   Bool get(const char* name, int *val);
   Bool get(const char* name, BigReal *val);
   Bool get(const char* name, Vector *val);
   Bool get(const char* name, StringList **val);
   Bool get(const char* name, char *val, int n=0);

   // number of elements for the given name
   int num(const char* name);
   // does the given element exist?
   Bool defined(const char *name);
   Bool exists(const char *name);

   // get/ set the units and scale for the given variable
   Bool units(const char *name, Units units);
   Bool units(const char *name, Units *units);
};

#endif

