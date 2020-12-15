/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Defines a set of dependencies, units, ranges, defaults, and info
   messages which simplify the parsing of a ConfigList
*/

#include <stdlib.h>
#include <string.h>
#include "ParseOptions.h"
#include "ConfigList.h"
#include "InfoStream.h"

#include "strlib.h"		//  For strcasecmp and strncasecmp

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif

// given the range, return the string
const char *rstring(Range r)
{
   switch (r) {
    case FREE_RANGE: return "unconstrained";
    case POSITIVE: return "positive";
    case NOT_NEGATIVE: return "non-negative";
    case NEGATIVE: return "negative";
    case NOT_POSITIVE: return "non-positive";
    default: return "error in rstring(Range )";
   }
}
   
static const char *unit_string_array[N_UNITS_UNDEFINED+1] =  {
  "", "fs", "ns", "sec", "min", "hr", "A", "nm", "m", 
  "kcal", "kJ", "eV", "K", "undefined units"
};
// given the units, return the string
const char *ustring(Units u) 
{
   return unit_string_array[u];
}
// little function so I can do for loops nicely
// static, so it never leaves the .C
Units static next(Units u) {
  switch (u) {
   case N_UNIT: return N_FSEC;
   case N_FSEC: return N_NSEC;
   case N_NSEC: return N_SEC;
   case N_SEC: return N_MIN;
   case N_MIN: return N_HOUR;
   case N_HOUR:  return N_ANGSTROM;
   case N_ANGSTROM: return N_NANOMETER;
   case N_NANOMETER: return N_METER;
   case N_METER: return N_KCAL;
   case N_KCAL: return N_KJOULE;
   case N_KJOULE: return N_EV;
   case N_EV: return N_KELVIN;
   case N_KELVIN: return N_UNITS_UNDEFINED;
   default: return N_UNITS_UNDEFINED;
  }
}

// convert from a string to Units
Units ParseOptions::atoUnits(const char *s) {
   Units u;
   for (u=N_UNIT; u!=N_UNITS_UNDEFINED; u = ::next(u)) {
      if (!strcasecmp(unit_string_array[u], s)) return u;
   }
   if (!strcasecmp(s, "Angstrom")) return N_ANGSTROM;
   if (!strcasecmp(s, "kcal/mol")) return N_KCAL;
   if (!strcasecmp(s, "kJ/mol")) return N_KJOULE;
   return N_UNITS_UNDEFINED;
}


// returns the scaling factor "sf" such that from*sf == to
// returns 0 if there is no conversion (like from METER to SEC)
static BigReal scaling_factors[N_UNITS_UNDEFINED+1] = {
   1, 1, 1000, 1E15, 60E15, 3600E15, 1, 10, 1E10, 1, 1/4.1855, 1/23.052,
   1, 0
};

// find the convertion factor "cf" such that to = cf * from
// or return 0 of there is an error
BigReal convert(Units to, Units from)
{
// cout << "Converting from " << string(from) << " to " << string(to) << std::endl;
// cout << scaling_factors[from] << " <--> " << scaling_factors[to] << std::endl;
   if (from == N_UNIT && to == N_UNIT) { return 1.0; }
   if ((from == N_NSEC || from == N_FSEC || from == N_SEC || from == N_MIN || 
        from == N_HOUR) &&
       (to   == N_NSEC || to   == N_FSEC || to   == N_SEC || to   == N_MIN || 
        to   == N_HOUR)) {
      return scaling_factors[from]/scaling_factors[to];
   }
   if ((from == N_METER || from == N_NANOMETER || from == N_ANGSTROM) &&
       (to   == N_METER || to   == N_NANOMETER || to   == N_ANGSTROM)) {
      return scaling_factors[from]/scaling_factors[to];
   }
   if ((from == N_KCAL || from == N_KJOULE || from == N_EV) &&
       (to   == N_KCAL || to   == N_KJOULE || to   == N_EV)) {
      return scaling_factors[from]/scaling_factors[to];
   }
   if (from == N_KELVIN && to == N_KELVIN) {
      return scaling_factors[from]/scaling_factors[to];
   }
   return 0.0;
}

static char *Strdup(const char *newname)
{
  char *tmp = new char[strlen(newname)+1];
  strcpy(tmp, newname);
  return tmp;
}

// Initialize a DataElement; this is called by all the constructors
void ParseOptions::DataElement::init(const char *newname,
     const char *newparent, int optional, const char *err) {
   name = Strdup(newname);
   index = -1;
   is_optional = optional;
   is_defined = FALSE;
   parent = Strdup(newparent);
   parent_ptr = NULL;
   error_message = Strdup(err);
   type = UNDEF;
   has_default = FALSE;
   many_allowed = FALSE;
   idef = 0;
   iptr = NULL;
   range = FREE_RANGE;
   units = N_UNIT;
}

#define dataelement_cons_macro_default(Mtype, MType, Mptr, Mdef) \
ParseOptions::DataElement::DataElement(const char *newname,      \
     const char *newparent, int optional, const char *err,       \
     Mtype *ptr, Mtype defalt)                                   \
{                                                                \
   init(newname, newparent, optional, err);                      \
   type = MType;                                                 \
   Mptr = ptr;                                                   \
   Mdef = defalt;                                                \
   has_default = TRUE;                                           \
   if ( ptr ) *ptr = defalt;                                     \
}

#define dataelement_cons_macro(Mtype, MType, Mptr) \
ParseOptions::DataElement::DataElement(const char *newname,      \
     const char *newparent, int optional, const char *err,       \
     Mtype *ptr)                                                 \
{                                                                \
   init(newname, newparent, optional, err);                      \
   type = MType;                                                 \
   Mptr = ptr;                                                   \
}

dataelement_cons_macro_default(BigReal, FLOAT, fptr, fdef);
dataelement_cons_macro_default(int, INT, iptr, idef);
dataelement_cons_macro_default(unsigned int, UINT, uiptr, uidef);
dataelement_cons_macro_default(Vector, VECTOR, vptr, vdef);

//dataelement_cons_macro_default(int, BOOL, iptr, idef);
dataelement_cons_macro(BigReal, FLOAT, fptr);
dataelement_cons_macro(Vector, VECTOR, vptr);
dataelement_cons_macro(int, INT, iptr);
dataelement_cons_macro(unsigned int, UINT, uiptr);
//dataelement_cons_macro(int, BOOL, iptr);
dataelement_cons_macro(char, STRING, sptr);


// This is a "STRINGLIST", which has NO default element
ParseOptions::DataElement::DataElement(const char *newname,
      const char *newparent, int optional, const char *err, StringList **ptr,
      int many)
{
   init(newname, newparent, optional, err);
   type = DataElement::STRINGLIST;
   slptr = ptr;
   has_default = FALSE;
   many_allowed = many;
}

// free up what needs to be freed
ParseOptions::DataElement::~DataElement(void) {
   if (name) delete[] name;
   if (parent) delete[] parent;
   if (error_message) delete[] error_message;
}
///////////////////////////////////////////////////// ParseOptions

// Initialize the data array to 20 element, and
// start it off with the "main" DataElement
ParseOptions::ParseOptions(void) {
   configList = NULL;
   array_size = 0;
   array_max_size = 20;
   data_array = new DataElement*[array_max_size];
   DataElement *tmp = new DataElement("main", "main", TRUE,
				      "Error in ParseOptions",
				      (int *) NULL,  0);
   tmp->type = DataElement::UNDEF;
   add_element(tmp);
}

// delete the data array
ParseOptions::~ParseOptions(void) {
   for (int i=0; i<array_size; i++) {
      delete data_array[i];
   }
   delete [] data_array;
}
   
// add the new element to the array,
void ParseOptions::add_element(DataElement *el) {
   if (array_size == array_max_size) {   // grow the array, if need be
      array_max_size += 30;
      DataElement **tmp = new DataElement*[array_max_size]; 
      memcpy(tmp, data_array, array_size * sizeof(DataElement *)); // copy
      delete [] data_array;
      data_array = tmp;
   }
   el->index = array_size;       // append the new element to the array
   data_array[array_size++] = el;
}

// update so that the elements are properly dependent
// returns 1 if everything set okay, 0 if no
int ParseOptions::make_dependencies(DataElement *el) {
   int i;
   // check if it is dependent on itself
   if (!strcasecmp(el->name, el->parent)) {
      return FALSE;
   }
   // check that there is no element with this name
   for (i=0; i<array_size; i++) {
      if (!strcasecmp(data_array[i]->name, el->name) &&
	  el != data_array[i]) {
	 return FALSE;
      }
   }

   // check if el's parent has already been inserted
   for (i=0; i<array_size; i++) {
      if (!strcasecmp(data_array[i]->name, el->parent)) {
	 el->parent_ptr = data_array[i];
	 break;
      }
   }
   
   // check if el is the parent of any already inserted
   // element
   for (i=0; i<array_size; i++) {
      if (!data_array[i]->parent_ptr) {
	 if (!strcasecmp(data_array[i]->parent,
			 el->name)) {
	    data_array[i]->parent_ptr = el;
	 }
      }
   }
   return TRUE;
}
   

/// routines to add dependencies to the array
#define parse_input_macro_default(fctnname, type, optional)           \
int ParseOptions::fctnname(const char *parent, const char *newname,   \
			  const char *msg, type *ptr, type defalt)    \
{                                                                     \
   DataElement *tmp = new DataElement(newname, parent, optional, msg, \
				      ptr, defalt);                   \
   if (!make_dependencies(tmp)) {                                     \
      iout << iERROR << "ParseOption '" << newname << "' already exists" << "\n" << endi; \
      return FALSE;                                                   \
   }                                                                  \
   add_element(tmp);                                                  \
   return TRUE;                                                       \
}
#define parse_input_macro(fctnname, type, optional)                   \
int ParseOptions::fctnname(const char *parent, const char *newname,   \
			  const char *msg, type *ptr)                 \
{                                                                     \
   DataElement *tmp = new DataElement(newname, parent, optional, msg, \
				      ptr);                           \
   if (!make_dependencies(tmp)) {                                     \
      iout << iERROR << "ParseOption '" << newname << "' already exists" << "\n" << endi; \
      return FALSE;                                                   \
   }                                                                  \
   add_element(tmp);                                                  \
   return TRUE;                                                       \
}
#define parse_input_macro_default_b(fctnname, type, optional, extra)  \
int ParseOptions::fctnname(const char *parent, const char *newname,   \
			  const char *msg, type *ptr, type defalt)    \
{                                                                     \
   DataElement *tmp = new DataElement(newname, parent, optional, msg, \
				      ptr, defalt);                   \
   if (!make_dependencies(tmp)) {                                     \
      iout << iERROR << "ParseOption '" << newname << "' already exists" << "\n" << endi; \
      return FALSE;                                                   \
   }                                                                  \
   add_element(tmp);                                                  \
   extra;                                                             \
   return TRUE;                                                       \
}
#define parse_input_macro_b(fctnname, type, optional, extra)                 \
int ParseOptions::fctnname(const char *parent, const char *newname,   \
			  const char *msg, type *ptr)                 \
{                                                                     \
   DataElement *tmp = new DataElement(newname, parent, optional, msg, \
				      ptr);                           \
   if (!make_dependencies(tmp)) {                                     \
      iout << iERROR << "ParseOption '" << newname << "' already exists" << "\n" << endi; \
      return FALSE;                                                   \
   }                                                                  \
   add_element(tmp);                                                  \
   extra;                                                             \
   return TRUE;                                                       \
}

parse_input_macro(require, BigReal, FALSE);  // the ; is there to look pretty
parse_input_macro(require, Vector, FALSE);
parse_input_macro(require, int, FALSE);
parse_input_macro(require, unsigned int, FALSE);
parse_input_macro_b(requireB, int, FALSE, tmp->type = DataElement::BOOL);
parse_input_macro(require, char, FALSE);

parse_input_macro(optional, BigReal, TRUE);
parse_input_macro(optional, Vector, TRUE);
parse_input_macro(optional, int, TRUE);
parse_input_macro(optional, unsigned int, TRUE);
parse_input_macro_b(optionalB, int, TRUE, tmp->type = DataElement::BOOL);
parse_input_macro(optional, char, TRUE);

parse_input_macro_default(require, BigReal, FALSE);
parse_input_macro_default(require, Vector, FALSE);
parse_input_macro_default(require, int, FALSE);
parse_input_macro_default(require, unsigned int, FALSE);
parse_input_macro_default_b(requireB, int, FALSE, tmp->type=DataElement::BOOL);

parse_input_macro_default(optional, BigReal, TRUE);
parse_input_macro_default(optional, Vector, TRUE);
parse_input_macro_default(optional, int, TRUE);
parse_input_macro_default(optional, unsigned int, TRUE);
parse_input_macro_default_b(optionalB, int, TRUE, tmp->type=DataElement::BOOL);

#define parse_stringlist_macro(fctn, xxx) \
int ParseOptions::fctn(const char *parent, const char *newname, \
			  const char *msg, StringList **ptr, int many_allowed)\
{                                                                            \
   DataElement *tmp = new DataElement(newname, parent, xxx, msg,             \
				      ptr, many_allowed);                    \
   if (!make_dependencies(tmp)) {                                            \
      iout << iERROR << "ParseOption '" << newname << "' already exists" << "\n" << endi;\
      return FALSE;                                                          \
   }                                                                         \
   add_element(tmp);                                                         \
   return TRUE;                                                              \
}
parse_stringlist_macro(require, FALSE);
parse_stringlist_macro(optional, TRUE);


// find all the children of the given element; returns 0
// if a loop was detected.  Children are marked with a 1 in the
// appropriate location in the flag array
int ParseOptions::check_children(int idx, int *flgs)
{
   if (flgs[idx]) { // oops, been here before
      return 0;
   }
   flgs[idx] = 1;
   for (int i=0; i<array_size; i++) {
      if (data_array[i]->parent_ptr == data_array[idx]) {
	 if (!check_children(i, flgs)) {
	    return 0;
	 }
      }
   }
   return 1;
}
	 
// see if there are elements which have no parent (except main)
// or elements inaccessible via main
// returns a 0 if there was an error
Bool ParseOptions::check_consistency(void) {
   int i;
   // check for lack of parent
   {
      int has_error = FALSE;
      for(i=1; i<array_size; i++) {
	 if (!data_array[i]->parent_ptr) {
	    // missing a parent
	    iout << iERROR << "Configuration element '" << data_array[i]->name
		    << "' defined, but the parent element" << "\n" << endi;
	    iout << iERROR << "  '" << data_array[i]->parent << "' is nowhere "
		    << "to be found" << "\n" << endi;
	    has_error = TRUE;
	 }
      }
      if (has_error) return 0;
   }

   // check for loop constructs in the "main" heirarchy
   int *arr = new int[array_size];
   for (i=0; i<array_size; i++) {  // initialize it
      arr[i] = 0;
   }
   if (!check_children(0, arr)) {
      // a loop was found
      iout << iERROR << "Loop found in ParseOptions data" << "\n" << endi;
      delete [] arr;
      return 0;
   }

   // check for elements inaccessible to "main"
   {
      int has_error = FALSE;
      for (i=1; i<array_size; i++) {
	 if (arr[i] == 0) {
	    // found an inaccesible element
	    if (has_error == FALSE) { // first time, so print message
	       iout << iERROR 
		   << "Found data in ParseOptions which are inaccessible "
		   << "to" << "\n" << endi;
	       iout << iERROR 
		  << "the main data hierarchy.  Errors in:" << "\n" << endi;
	       has_error = TRUE;
	    }
	    iout << iERROR << "   '" << data_array[i]->name << "' depends on '"
		    << data_array[i]->parent << "'" << "\n" << endi;
	 }
      }
      if (has_error) {
	 delete [] arr;
	 return 0;
      }
   }
   // looks like everything went well
   delete [] arr;
   return 1;
}

// convert from a string to Bool; returns 1(TRUE) 0(FALSE) or -1(if unknown)
int ParseOptions::atoBool(const char *s)
{
   if (!strcasecmp(s, "on")) return 1;
   if (!strcasecmp(s, "off")) return 0;
   if (!strcasecmp(s, "true")) return 1;
   if (!strcasecmp(s, "false")) return 0;
   if (!strcasecmp(s, "yes")) return 1;
   if (!strcasecmp(s, "no")) return 0;
   if (!strcasecmp(s, "1")) return 1;
   if (!strcasecmp(s, "0")) return 0;
   return -1;
}

// A "Boolean" can indicate 1 of 2 things, either a simple
// yes/ no flag for a piece of data, or a on/ off for a
// set of parameters.  The distinction is if the option
// has children.  In "set", if the option is defined and
// is false and has children, then it is made undefined,
// so that other elements won't be parsed
Bool ParseOptions::is_parent_node(DataElement *el) {
   for (int i=1; i<array_size; i++) {
      if (data_array[i]->parent_ptr == el) {
	 return 1;
      }
   }
   return 0;
}

// given a string, read in the float and (perhaps) units,
// do the needed conversions, check for errors, etc...
Bool ParseOptions::scan_float(DataElement *data, const char *s)
{
   double input_value;  // I do this since I don't know if "BigReal" is
   BigReal fval;        // float or double, and I need to know for the
   char units_str[80];  // sscanf
   char tmp_str[80];
   int count = sscanf(s, "%lf%s%s", &input_value, units_str, tmp_str);
   if (count > 1 && units_str[0] == '.' &&
       input_value == (double)(long int)input_value) {  // for final . on Mac
     long int input_long;
     count = sscanf(s, "%ld.%s%s", &input_long, units_str, tmp_str);
     if ( count < 1 || input_value != (double)input_long ) {
	 iout << iERROR << "Could not parse option '"
		 << data->name << " = " << s << "'\n" << endi;
	 return FALSE;
     }
   }
   fval = input_value;
   if (count == 1) {     // no units given, so simply apply the number
      data->fdata = fval;
      return TRUE;
   }
   if (count == 2) {     // number and units
      Units u = atoUnits(units_str);
      if (u == N_UNITS_UNDEFINED) {
	 iout << iERROR << "Could not understand units '" << units_str
		 << "' in option '" << data->name << " = " << s
                 << "\n" << endi;
	 return FALSE;
      }
      BigReal scale = convert(data->units, u);
      if (scale == 0) {
	 iout << iERROR << "Could not translate from units '" << ustring(u) 
	         << "' to '" << ustring(data->units) << "' for option '"
                 << data->name << "'" << "\n" << endi;
	 return FALSE;
      }
//      cout << "fval == " << fval << "  scale == " << scale << std::endl;
      data->fdata = fval * scale;
      return TRUE;
   }
   if (count <=0) {  // not enough
      iout << iERROR << "Expecting value and optional units for option '"
              << data->name << "'" << "\n" << endi;
   }
   if (count > 2) {  // too many
      iout << iERROR << "Too much information given to '" << data -> name 
      << " = " <<  s << "'" << "\n" << endi;
      iout << iERROR << "  - expecting a value and optional units" << "\n" << endi;
   }
   return FALSE;
}

// no units are allowed (yet?) for a vector
Bool ParseOptions::scan_vector(DataElement *data, const char *s)
{
   Vector v;
   if (!v.set(s)) {
      iout << iERROR << "Could not translate the value '" << s << "'" << "\n" << endi;
      iout << iERROR << "  into a Vector for the option '" << data->name << "'"
              << "\n" << endi;
      return FALSE;
   }
   data->vdata = v;
   return TRUE;
}

// read an int from string.  No units are supported for ints,
// though they could be.  (eg, 1 KB, ...)
Bool ParseOptions::scan_int(DataElement *data, const char *s)
{
   int ival;
   char units_str[80];
   char tmp_str[80];
   int count = sscanf(s, "%d%s%s", &ival, units_str, tmp_str);
   if (count == 1) {
      data->idata = ival;
      return TRUE;
   }
   iout << iERROR << "Expecting only a number for '" << data->name
	   << "' input, got: " << s << "\n" << endi;
   return FALSE;
}

// read an unsigned int from string.  No units are supported for ints,
// though they could be.  (eg, 1 KB, ...)
Bool ParseOptions::scan_uint(DataElement *data, const char *s)
{
   unsigned int ival;
   char units_str[80];
   char tmp_str[80];
   int count = sscanf(s, "%u%s%s", &ival, units_str, tmp_str);
   if (count == 1) {
      data->idata = ival;
      return TRUE;
   }
   iout << iERROR << "Expecting only a number for '" << data->name
	   << "' input, got: " << s << "\n" << endi;
   return FALSE;
}
   
Bool ParseOptions::scan_bool(DataElement *data, const char *s)
{
   int tmp = atoBool(s);  // convert the boolean

   if (tmp == -1) 
   {
  	iout << iERROR << "ParseOptions can't understand '" << s << "' for the "
	        << "\n" << endi;
  	iout << iERROR << " Boolean variable '" << data->name << "'"  << "\n" << endi;

  	data->idata = FALSE;

        return FALSE;
   }

   data->idata = tmp;  // set the value, if understood

   return TRUE;
}

// check the range and, if valid, set the variable and return 1
// if bad range, print error message and return 0
#define set_macro(type, field, fieldptr)                 \
int ParseOptions::set_##type(DataElement *el)            \
{                                                        \
   if (el->range == FREE_RANGE ||                        \
       (el->range == POSITIVE && el->field > 0) ||       \
       (el->range == NOT_NEGATIVE && el->field >= 0) ||  \
       (el->range == NEGATIVE && el->field < 0) ||       \
       (el->range == NOT_POSITIVE && el->field <= 0)) {  \
      if (el->fieldptr) *(el->fieldptr) = el->field;     \
      return 1;                                          \
   }                                                     \
   iout << iERROR << "'" << el->name << "' was set to " << el->field << " but it " \
	   << "should be " << rstring(el->range)          \
	   << "\n" << endi;      \
   return 0;                \
}
set_macro(float, fdata, fptr);
set_macro(int, idata, iptr);
set_macro(uint, uidata, uiptr);

// for elements without ranges
#define simple_set_macro(type, field, fieldptr)   \
void ParseOptions::set_##type(DataElement *el)     \
{                                                   \
   if (el->fieldptr) *(el->fieldptr) = el->field;    \
}

simple_set_macro(bool, idata, iptr); 
simple_set_macro(vector, vdata, vptr);
simple_set_macro(stringlist, sldata, slptr);
// simple_set_macro(string, sldata->data, sptr);

void ParseOptions::set_string(DataElement *el)   
{                                                  
   if (el->sptr) strcpy(el->sptr, el->sldata->data);  
}

// set the variables, given the contents of the ConfigList
// return FALSE if there was an error
Bool ParseOptions::set(const ConfigList& clist)
{
   // the algorithm is easy, though it looks scary
   int cont = TRUE, i;  // do some initialization
   StringList *slptr;
   DataElement *data;
   int has_error = FALSE;
   int *checked = new int[array_size];
   configList = &clist;

   // I make here a list of element I have already checked, starting
   // at the head.  I check only children of those that have already
   // been defined and add it to the checked list
   for (i=0; i<array_size; i++) 
   {
      checked[i] = FALSE;
   }

   //   make "main" 'defined' (at this point, nothing else should be defined)
   data_array[0]->is_defined = TRUE;
   checked[0] = TRUE;  // and make "main" checked
   
   // while there is still data which hasn't been defined
   while (cont) 
   {
      cont = FALSE;
      for (i=1; i<array_size; i++) 
      {  // check each element
	 data = data_array[i];

	 // find unchecked data which has a parent which was checked
	 // and defined
	 if (!checked[data->index] &&
	     checked[data-> parent_ptr -> index] &&
	     data -> parent_ptr -> is_defined) 
	 {
	    cont = TRUE;
	    checked[data->index] = TRUE;  // so I don't check again

	    // check to see if data is available in the StringList
	    slptr = clist.find(data->name);

	    if (slptr != NULL) 
	    {  // it is 

	       // most data types allow only 1 item, so check if that is
	       // a problem.  (some StringLists, like 'parameters', allow
	       // multiple strings)
	       if (!data->many_allowed && slptr->next != NULL) 
	       {
		  iout << iERROR << "Multiple definitions of '" << data->name << "'" << "\n" << endi;
  		  iout << iERROR << "  in the configuration file are not allowed" << "\n" << endi;
  		  has_error = TRUE;
               }

	       data->is_defined = TRUE;

	       // set the appropriate data field
	       if (data->type == DataElement::FLOAT) 
	       {
		  if (!scan_float(data, slptr->data)) 
			has_error = TRUE;
	       } 
	       else if (data->type == DataElement::VECTOR) 
	       {
		  if (!scan_vector(data, slptr->data)) 
			has_error = TRUE;
	       } 
	       else if (data->type == DataElement::INT) 
	       {
		  if (!scan_int(data, slptr->data)) 
			has_error = TRUE;
	       } 
	       else if (data->type == DataElement::UINT) 
	       {
		  if (!scan_uint(data, slptr->data)) 
			has_error = TRUE;
	       } 
	       else if (data->type == DataElement::BOOL) 
	       {
		  if (!scan_bool(data, slptr->data)) 
			has_error = TRUE;
	       } 
	       else if (data->type == DataElement::STRINGLIST ||
			  data->type == DataElement::STRING ) 
	       {
		  data->sldata = slptr;
	       } 
	       else 
	       {
		  iout << iERROR << "Unknown ParseOption data type " << (int)(data->type) << " for "
			  << "variable " << data->name << "\n" << endi;
  		  has_error = TRUE;
	       }
	    } 
	    else 
	    {  // no definition; is there a default?
	       if (data->has_default) 
	       {
		  data->is_defined = TRUE;

		  if (data->type == DataElement::FLOAT) 
		  {
		     data->fdata = data->fdef;
		  } 
		  else if (data->type == DataElement::VECTOR) 
		  {
		     data->vdata = data->vdef;
		  } 
		  else if (data->type == DataElement::UINT) 
		  {
		     data->uidata = data->uidef;
		  } 
		  else if (data->type == DataElement::INT ||
			     data->type == DataElement::BOOL) 
		  {
		     data->idata = data->idef;
		  } 
		  else 
		  {
   			iout << iERROR << "Unknown ParseOption data type " << (int)(data->type) << " for "
				<< "variable " << data->name << "\n" << endi;
  		    has_error = TRUE;
		  }
	       }
	    }

	    // at this point we should have gotten data from the file or the defaults,
  	    // or it hasn't yet been defined.  If it still isn't defined, check
            // to see if it is optional
	    if (!data->is_defined) 
	    { // if still not defined,
	       if (!data->is_optional) 
	       { // it is it required
		  has_error = TRUE;
		  iout << iERROR << "'" << data->name << "' is a required configuration option" << "\n" << endi;

  		  // print some helpful information if this isn't a "main" option
                  if (data->parent_ptr != data_array[0]) 
		  {
			iout << iERROR << "  when '" << data->parent_ptr -> name << "' is set" << "\n" << endi;
                  }  // printed parent info

  		  iout << iERROR << data->name << " defines:   " << data->error_message << "\n" << endi;

	       }  // printed error message
	    }  
	    else 
	    { // there was a definition, so assign to the variable
	       if (data->type ==  DataElement::FLOAT) 
	       {
		  if (!set_float(data)) 
			has_error = TRUE;
	       } 
	       else if ( data -> type == DataElement::VECTOR) 
	       {
		  set_vector(data);
	       } 
	       else if ( data -> type == DataElement::INT) 
	       {
		  if (!set_int(data)) 
			has_error = TRUE;
	       } 
	       else if ( data -> type == DataElement::UINT) 
	       {
		  if (!set_uint(data)) 
			has_error = TRUE;
	       } 
	       else if ( data -> type == DataElement::BOOL) 
	       {
		  set_bool(data);

		  if (is_parent_node(data)) 
   		  {
      			// this makes the boolean variable undefined if it is 'off'
      			// _and_ it is a parent; this makes it agree with namd's
      			// configuration option semantics
      			data->is_defined = data->idata;
   		  }
	       } 
	       else if ( data -> type == DataElement::STRINGLIST) 
	       {
		  set_stringlist(data);
	       } 
	       else if ( data -> type == DataElement::STRING) 
	       {
		  set_string(data);
	       } 
	       else 
	       {
		  // already printed the error message
	       }
	    }
	 }  // end of checking the available variables
      } // end of pass through the list
   } // end of finding data in the ConfigList


   // and now print the warning messages

   // first, find elements which are in the configuration file and are
   // valid, but which were not needed (ie, the checked array wasn''t set)
   {
      int flg = 0;

      for (int i=1; i<array_size; i++) 
      {
	 if (!checked[i]) { // wasn't needed
	    data = data_array[i];
	    if (clist.find(data->name)) {
	       if (flg == 0) {
		  flg = 1;
                  iout << iWARN 
		    << "The following variables were set in the\n";
		  iout << iWARN 
		    << "configuration file but will be ignored:\n" << endi;
	       }
	       iout << iWARN << "   " << data->name;
               if (data->parent_ptr != data_array[0]) {
	         iout << " (" << data->parent_ptr->name << ")";
	       }
	       iout << "\n" << endi;
	    }
	 }
      }
   }
   // and now look for names which are in the config list but which
   // are not in the parseoptions list
   {
      int flg = 0;
      ConfigList::ConfigListNode const *ptr;
      for (ptr = clist.head(); ptr != NULL; ptr = ptr -> next) {
	 if (!exists(ptr -> name)) {
	    if (flg == 0) {
	       flg = 1;
	       has_error = TRUE;
               iout << iERROR
		  << "The following variables were set in the\n";
               iout << iERROR
		  << "configuration file but are NOT VALID\n" << endi;
	    }
	    iout << iERROR << "   " << ptr -> name << "\n" << endi;
	 }
      }
   }
      
   delete [] checked;
   return !has_error;
}

// simple search though the list; return NULL if not found
ParseOptions::DataElement *ParseOptions::internal_find(const char* name)
{
   for (int i=1; i<array_size; i++) {
      if (!strcasecmp(name, data_array[i]->name)) {
	 return data_array[i];
      }
   }
   return NULL;
}

// returns 1 if the given name exists; 0 if no
Bool ParseOptions::exists(const char *name)
{
   if (!name) return 0;
   if (internal_find(name)) {
      return 1;
   }
   return 0;
}
// returns 1 if the element name exists and is defined
Bool ParseOptions::defined(const char *name)
{
   if (!name) return FALSE;
   DataElement *tmp = internal_find(name);
   if (!tmp) return FALSE;
   if (tmp->is_defined) {
      return TRUE;
   }
   return FALSE;
}

#ifdef NAMD_TCL
#define PRINT_DOUBLE(BUF,VAL) Tcl_PrintDouble(0,VAL,BUF)

static void PRINT_VECTOR(char *buf, Vector val) {
  PRINT_DOUBLE(buf, val.x);
  buf += strlen(buf); buf[0] = ' '; ++buf;
  PRINT_DOUBLE(buf, val.y);
  buf += strlen(buf); buf[0] = ' '; ++buf;
  PRINT_DOUBLE(buf, val.z);
  buf += strlen(buf); buf[0] = ' '; buf[1] = 0;
}
#endif

//// special accessor for ScriptTcl - uses pointer value if available
char* ParseOptions::getfromptr(const char* name, char *outbuf) {
#ifdef NAMD_TCL
   if ( ! name ) NAMD_bug("ParseOptions::getfromptr called with null name");
   if ( ! outbuf ) NAMD_bug("ParseOptions::getfromptr called with null outbuf");
   DataElement *el = internal_find(name);
   if ( el == NULL ) return 0;
   switch (el->type) {
    case DataElement::FLOAT :
      if ( el->fptr ) PRINT_DOUBLE(outbuf, *(el->fptr));
      else PRINT_DOUBLE(outbuf, el->fdata);
      return outbuf;
    case DataElement::INT:
    case DataElement::BOOL:
      if ( el->iptr ) sprintf(outbuf,"%d", *(el->iptr));
      else sprintf(outbuf,"%d", el->idata);
      return outbuf;
    case DataElement::STRINGLIST :
      if ( el->slptr ) return (*(el->slptr))->data;
      else if ( el->sldata ) return el->sldata->data;
      else return 0;
    case DataElement::STRING :
      if ( el->sptr ) return el->sptr;
      else if ( el->sldata ) return el->sldata->data;
      else return 0;
    case DataElement::VECTOR :
      if ( el->vptr ) PRINT_VECTOR(outbuf, *(el->vptr));
      else PRINT_VECTOR(outbuf, el->vdata);
      return outbuf;
    default:
      iout << iERROR 
	 << "Unknown data type " << (int)(el->type) << " for '" << name << "'"
	 << "\n" << endi;
   }
#endif
   return 0;
}

//// special accessor for ScriptTcl - uses pointer value if available
int ParseOptions::istruefromptr(const char* name) {
   if ( ! name ) NAMD_bug("ParseOptions::getfromptr called with null name");
   DataElement *el = internal_find(name);
   if ( el == NULL ) return -1;
   if ( el->type != DataElement::BOOL ) return -2;
   if ( el->iptr ) return ((*(el->iptr)) ? 1 : 0);
   if ( ! el->is_defined ) return -3;  // ignores defaults
   return (el->idata ? 1 : 0);
}

//// special accessor for ScriptTcl - uses pointer value if available
int ParseOptions::issetfromptr(const char* name) {
   if ( ! name ) NAMD_bug("ParseOptions::getfromptr called with null name");
   DataElement *el = internal_find(name);
   if ( el == NULL ) return -1;
   return (el->is_defined ? 1 : 0);
}

//// accessors
// get the value associated with the given name
// if there was a type conversion, print a warning message
// if the element wasn't found, return a FALSE, else return a TRUE
Bool ParseOptions::get(const char* name, int *val) {
   if (!val) return FALSE;
   DataElement *el = internal_find(name);
   if (el == NULL || !el->is_defined) {
      return FALSE;
   }
   switch (el->type) {
    case DataElement::FLOAT :
      iout << iWARN 
	 << "ParseOptions doing a conversion from float to int for '"
	 << name << "'" << "\n" << endi;
      *val = (int) el->fdata;
      return TRUE;
    case DataElement::INT:
    case DataElement::BOOL:
      *val = el->idata;
      return TRUE;
    case DataElement::STRINGLIST :
    case DataElement::STRING :
      iout << iWARN 
	 << "ParseOptions doing a conversion from StringList[0] to int "
	 << "for '" << name << "'" << "\n" << endi;
      *val = atoi(el->sldata->data);
      return TRUE;
    case DataElement::VECTOR :
       iout << iERROR 
	  << "ParseOptions cannot convert from Vector to int for '"
          << name << "'" << "\n" << endi;
       return FALSE;
    default:
      iout << iERROR 
	 << "Unknown data type " << (int)(el->type) << " for '" << name << "'"
	 << "\n" << endi;
   }
   return FALSE;
}

Bool ParseOptions::get(const char* name, BigReal *val) {
   if (!val) return FALSE;
   DataElement *el = internal_find(name);
   if (el == NULL || !el->is_defined) {
      return FALSE;
   }
   switch (el -> type) {
    case DataElement::FLOAT: 
      *val =  el->fdata;
      return TRUE;
    case DataElement::INT:
      iout << iWARN 
	 << "ParseOptions doing a conversion from int to float '"
	 << name << "'" << "\n" << endi;
      *val = (BigReal) el->idata;
      return TRUE;
    case DataElement::BOOL:
      iout << iWARN 
	 << "ParseOptions doing a conversion from boolean to float for '"
	 << name << "'" << "\n" << endi;
      *val = (BigReal) el->idata;
      return TRUE;
    case DataElement::STRING:
    case DataElement::STRINGLIST:
     iout << iWARN 
	<< "ParseOptions doing a conversion from StringList[0] to float "
        << "for '" << name << "'" << "\n" << endi;
      *val = atof(el->sldata->data);
      return TRUE;
    case DataElement::VECTOR :
       iout << iERROR 
	  << "ParseOptions cannot convert from Vector to float for '"
          << name << "'" << "\n" << endi;
       return FALSE;
    default:
      iout << iERROR 
	 << "Unknown data type " << (int)(el->type) << " for '" << name << "'"
	 << "\n" << endi;
   }
   return FALSE;
}
Bool ParseOptions::get(const char *name, Vector *val) {
   if (!val) return FALSE;
   DataElement *el = internal_find(name);
   if (el == NULL || !el->is_defined) {
      return FALSE;
   }
   switch (el -> type) {
    case DataElement::FLOAT:
      iout << iERROR 
	 << "ParseOptions cannot convert from float to Vector for '"
         << name << "'" << "\n" << endi;
      return FALSE;
    case DataElement::INT:
      iout << iERROR 
	 << "ParseOptions cannot convert from int to Vector for '"
         << name << "'" << "\n" << endi;
      return FALSE;
    case DataElement::STRING:
    case DataElement::STRINGLIST:{
      iout << iWARN 
	 << "ParseOptions doing a conversion from StringList[0] to "
         << "Vector for '" << name << "'" << "\n" << endi;
      Vector v;
      if (!v.set(el->sldata->data)) {
	 iout << iERROR 
	    << "Could not convert '" << el->sldata->data
	    << "' to a Vector";
         return FALSE;
      }
      *val = v;
      return TRUE;
     }
     case DataElement::VECTOR : 
      *val = el->vdata;
      return TRUE;
    default:
      iout << iERROR 
	 << "Unknown data type " << (int)(el->type) << " for '" << name << "'"
	 << "\n" << endi;

   }
   return FALSE;
}
Bool ParseOptions::get(const char *name, StringList **val) {
   if (!val) return FALSE;
   DataElement *el = internal_find(name); // first check it is internally valid
   if (el == NULL || !el->is_defined) {
      return FALSE;
   }
   // then simply ask the configList itself for the answer
   // (while I do keep the information internally, in sldata, why bother?
   if (!configList) { return FALSE; }
   *val = configList->find(name);
   if (!*val) { return FALSE; }  // paranoia, I know...
   return TRUE;
}

// get the nth element of the StringList
Bool ParseOptions::get(const char *name, char *val, int n)
{
   if (!val || n<0) {return FALSE;}
   StringList *tmp;
   if (!get(name, &tmp)) {val[0]=STRINGNULL; return FALSE; }
   int i=n;
   while (i>0 && tmp) { // search for the nth element
      tmp=tmp->next;
      i--;
   }
   if (tmp) {  // if it was long enough, return it
      strcpy(val, tmp->data);
      return TRUE;
   }
   val[0] = STRINGNULL;
   return FALSE;
}


int ParseOptions::num(const char *name)
{
   DataElement *el = internal_find(name);
   if (!el || !el ->is_defined) {
      return 0;
   }
   if (!el->many_allowed) {
      return 1;
   }
   StringList *tmp;
   if (!get(name, &tmp)) { return 0; }
   int i=0;
   while (tmp) {
      i++;
      tmp=tmp->next;
   }
   return i;
}
   

// get or set the range for the given variable
void ParseOptions::range(const char *name, Range newrange)
{
   DataElement *el = internal_find(name);
   if (!el) {
      iout << iERROR 
	 << "Trying to set the range of undefined variable '"
	 << name << "'" << "\n" << endi;
      return;
   }
   el->range = newrange;
   
}
Range ParseOptions::range(const char *name)
{
   DataElement *el = internal_find(name);
   if (!el) {
      iout << iERROR 
	 << "Trying to get the range of undefined variable '"
	 << name << "'" << "\n" << endi;
      return FREE_RANGE;
   }
   return el->range;
}

// get / set the units value
Bool ParseOptions::units(const char *name, Units units)  // get
{
  DataElement *tmp = internal_find(name);
  if (!tmp) {
     iout << iERROR 
	<< name << " not found so units not set" << "\n" << endi;
     return FALSE;
  }
  if ((tmp -> type == DataElement::INT && units != N_UNIT) ||
      (tmp -> type != DataElement::INT && 
       tmp -> type != DataElement::FLOAT)) {
     iout << iERROR 
	<< "Cannot set units '" << ustring(units) << "' for option '"
        << name << "'; wrong data type" << "\n" << endi;
     return FALSE;
  }
  tmp -> units = units;
  return TRUE;
}

Bool ParseOptions::units(const char *name, Units *units) // set
{
   DataElement *tmp = internal_find(name);
   *units = N_UNIT;
   if (!tmp) {
      iout << iERROR 
	 << "'" << name << "' doesn't exist so cannot get its units"
         << "\n" << endi;
      return FALSE;
   }
   if (tmp -> type != DataElement::INT && 
       tmp -> type != DataElement::FLOAT) {
      iout << iERROR 
	 << "Can only get units for FLOAT and INT variables, and '"
	 << name << "' isn't one of those" << "\n" << endi;
      return FALSE;
   }
   *units = tmp->units;
   return TRUE;
}

