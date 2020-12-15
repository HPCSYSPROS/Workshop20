/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
     Read in a configuration file of the form:
         keyword = information\n
   -or-  keyword information\n
   and produces a database which can return a linked list of strings (char *)
   to all the information fields associated with that keyword.
   
      A "word" is a seqeunce of characters that are not white space (see
   isspace(3C) ).  The equals sign ('=') is optional (though if there is more
   more than one equals sign, then the 2nd is not ignored).  The "information"
   field may contain more than one word.  White space is ignored except that
   white space between multiple words in the information field is maintained.
   Everything on the line at and beyond a pound sign ('#') is ignored.  Hence
   a data file can be:
     fullname = George   Washington # the first president of the U.S.
     fullname = Martha Washington   # his second wife
   Once that is read in, all data associated with "name" can be retreived as
    StringList *strList = configFile.find("fullname");
    for (StringList *tmp=strList; tmp!=NULL; tmp = tmp -> next)
        std::cout << tmp->data << '\n';
   Note:
     The returned StringList * is NOT new'ed.  Do NOT free it.
     Keywords are case INsensitive
*/

#include <string.h> // for strncpy, strcasecmp
#include <ctype.h>   // for isspace
#include <stdio.h>   // Yes, I do use stdio
#include "InfoStream.h"
#include "ConfigList.h"
#include "common.h"

#include "strlib.h"	// for strncpy, strcasecmp

// given a key word, find the element of theList that matches
// This is a private function.
ConfigList::ConfigListNode *ConfigList::find_key_word(const char *key) const 
{
 ConfigListNode *tmp;
 for (tmp=theList; tmp!= NULL; tmp = tmp -> next) {
   if (!strcasecmp(tmp->name, key))
     return tmp;
 }
 return NULL;
}

//  I have a new name/ data pair.  If the name does not already exist,
// make a new ListData element and put it on the HEAD of theList.  If
// name does already exist, make a new ListNode and add the name to the end
//   String 1 starts at s1 and is of length len1
//      "   2    "   "  s2  "   "  "   "    len2
void ConfigList::add_element(const char *s1, int len1, const char *s2, int len2)
{
    char *temps = new char[len1 + 1];  // what is the name?
    if ( temps == NULL )
    {
      NAMD_die("new failed in ConfigList::add_element");
    }
//printf("%s %s\n",s1,s2);
#ifndef NAMD_NO_STDOUT_FLUSH
    fflush(stdout);
#endif
    strncpy(temps, s1, len1);
    temps[len1] = 0;                   //       terminate the string
    ConfigListNode *tmpList;
                                       
    tmpList = find_key_word( temps);  // go through the list
    if (tmpList == NULL )  {          // if not found
       tmpList = new ConfigListNode( theList, temps, NULL);// create a new node
       if ( tmpList == NULL )
       {
	 NAMD_die("new failed in ConfigList::add_element");
       }
                                    // and stick it on the head of the list
       theList = tmpList;           // note that I can continue to use tmpList
    }

    if (len1 < len2) {                  // if string is smaller, don't realloc
      delete [] temps;
      temps = new char[len2 + 1];
      if ( temps == NULL )
      {
	NAMD_die("new failed in ConfigList::add_element");
      }
    }
    strncpy(temps, s2, len2);           // get the new string
    temps[len2] = 0;

    StringList *newStrList = new StringList(temps);  // new element
    if ( newStrList == NULL )
    {
      NAMD_die("new failed in ConfigList::add_element");
    }
    newStrList -> next = NULL;
    if (tmpList -> data == NULL) {       // no previous definition
      tmpList -> data = newStrList;      //    so this is the new head
    } else {
      StringList *tmpStrList = tmpList -> data;   // else,
      while (tmpStrList -> next != NULL)          // find the end of the list
        tmpStrList = tmpStrList -> next;
      tmpStrList -> next = newStrList;            // and stick it on the tail
    }
    
    delete [] temps;
}


// The scripting interface will call our add_element routine.

ConfigList::ConfigList(void)
{
  isokay = TRUE;
  theList = NULL;
}

// Used to implement "source" for file parser below.

struct FileStack {
  FILE *file;
  int linenumber;
  char *filename;
  FileStack *next;
};

// open, read, parse, and close the file
// make a linked list of "AssocList"
// The file is parsed as (I think):
// [w]?([W])*[w]?[=]?[w]?([W][wW])?[w]
// where w == white space; W is non-whitespace
// [x]? says that there can be 0 or more characters of type x.
// [x]* says there must be at least 1 character of type x
// the terms in () are the two that are used.  Everything after a '#'
// sign is excluded.  I'm not sure how to specify that exclusion in
// my regular expression

// a new one:  if the first character of the data is '{'
//   then I append new keyword values for each line until I get to
//   a line with the first non-blank character as a '}'
ConfigList::ConfigList(const char *filename_in)
{
  char *filename = new char[strlen(filename_in)+1];
  strcpy(filename,filename_in);
  FileStack *fileStack = 0;
  FILE *infile;
  
  isokay = FALSE;
  theList = NULL;
  int linenumber = 0;   // keep track of line numbers for searching out errors

  if (!strcmp(filename,"-")) {  // should the input be from stdin?
    infile = stdin;
  } else {
    if ( (infile = Fopen(filename, "r")) == NULL ) {
        iout << iWARN << "Unable to open configuration file '" 
                 << filename << "'.\n" << endi;
        isokay = FALSE;
        return;
    }
  }
  isokay = TRUE;         // file is now open

     // so read and parse it
  char buf[1000]; // give myself lots of space
  char *namestart, *nameend, *datastart, *dataend;
  char *s;
  int spacecount;
  int fileok;
  while ((fileok = ! ! fgets(buf, 999, infile)) || fileStack) {
    if ( fileStack && ! fileok ) { // done with "source"
      delete [] filename;
      filename = fileStack->filename;
      linenumber = fileStack->linenumber;
      infile = fileStack->file;
      FileStack *delStack = fileStack;
      fileStack = fileStack->next;
      delete delStack;
      continue;
    }
    linenumber ++;        
    namestart = nameend = datastart = dataend = NULL;
    spacecount = 0;
    
    for (s = buf; *s && *s!='\n'; s++) {    // get to the end of the line
       if (*s == '#')                       // found a comment, so break
          break;
       if ( !isspace(*s) )    // dataend will always be the last non-blank char
          dataend = s;
       if ( !isspace(*s) && !namestart)     // found first character of name
          {namestart = s; continue; }
       if ( (isspace(*s)  || *s == '=') &&  // found last character of name
                 namestart && !nameend)
          nameend = s - 1;
       if ( !isspace(*s) && !datastart &&   // found the next char. after name
                 nameend)
          if (*s == '=' && spacecount == 0) // an equals is allowed
             {spacecount++; continue; }     // but only once
            else
             {datastart = s; continue; }    // otherwise, use it
    }
    if (*s == '\n')          // cut out the newline at the end
      *s = 0;
     else
      if (*s != '#') {       // then there was an overflow
        iout << iWARN << "Line " << linenumber << " of configuration file "
                 << filename << " contains more than 999 characters."
                 << "  Excess characters will be ignored.\n" << endi;
      } else {
        *s = 0;  // delete the '#' character
      }

// I will also ignore line that I can't understand
// If there is any text on the line (as compared to a blank or commented)
//   line, then I will say that there is a problem.
    if (!namestart || !nameend || !datastart || !dataend) {
      if (!namestart && datastart || namestart && !datastart) {// was some data
        iout << iWARN << "Couldn't parse line " << linenumber << " in "
                 << "configuration file " << filename << ".  The line was: "
                 << buf << "\n" << endi;
      }
      continue;  // which ever the case, go to the next line
    }

   if ( ! strncmp(namestart, "source", nameend-namestart+1) )  {
     // see if the the name is "source"

     // store the old file data
     FileStack *newStack = new FileStack;
     newStack->filename = filename;
     newStack->linenumber = linenumber;
     newStack->file = infile;
     newStack->next = fileStack;
     fileStack = newStack;

     // copy the filename
     char *cpychar = new char[dataend-datastart+2];
     strcpy(cpychar,datastart);
     filename = cpychar;

     // open the sourced file
     if ( (infile = Fopen(filename, "r")) == NULL ) {
        iout << iWARN << "Unable to open file '" 
                 << filename << "' sourced by '"
		<< fileStack->filename << "' at line "
		<< fileStack->linenumber << ".\n" << endi;
        isokay = FALSE;
        return;
     }
     iout << iINFO << "Sourcing " << filename << "\n" << endi;
     isokay = TRUE;         // file is now open
     linenumber = 0;
     
   } else if (datastart[0] == '{') {
     // check if the data begins with a '{'
     // will remove initial '{' and final '}' to match Tcl.
     std::ostringstream alldata;
     char newdata[1000];
     int found_end = 0;
     ++datastart;  // remove initial '{'
     int open_brace_count = 1;
     char *newline = datastart;
     strcat(newline,"\n");  // put back newline that was removed above

     while ( 1 ) {
       int i;
       int escape_next = 0;
       for (i=0; i<1000; i++) {
	 if (! newline[i]) {
	   break;
         }
	 if (escape_next) {
	   escape_next = 0;
	   continue;
	 }
	 if (newline[i] == '\\' && ! escape_next) {
	   escape_next = 1;
	 }
	 if (newline[i] == '{' && ! escape_next) {
	   ++open_brace_count;
	 }
	 if (newline[i] == '}' && ! escape_next) {
	   if ( ( found_end = ! --open_brace_count ) ) {
	     newline[i] = 0;
	     break;
	   }
	 }
       }
       alldata << newline;
       if (found_end) break;
       newline = newdata;
       if ( ! fgets(newdata, 999, infile) ) break;
       linenumber ++;
     }
     std::string alldatastr = alldata.str();
     add_element(namestart, nameend-namestart+1, alldatastr.c_str(), alldata.str().length());
     // delete string?
     if (!found_end) {
       *(nameend+1) = 0;
       sprintf(newdata, "configuration file ended early while parsing line "
	       "%d of keyword structure %s", linenumber, namestart);
       NAMD_die(newdata);
     }
   } else {
     // now I can add the new values to the linked list
     add_element( namestart, nameend - namestart + 1, datastart,
                  dataend - datastart + 1 );
   }
  } // while I can still get data with fgets
  
  if (strcmp(filename,"-")) {  // close the input file if not stdin
    Fclose(infile);
  }
  delete [] filename;
}

// destructor for the class - just delete a linked list
ConfigList::~ConfigList( void) 
{
  ConfigListNode *curr=theList, *next=NULL;

  while ( curr!=NULL ) 
  {
     next = curr -> next;
     delete curr;     // the nasties for deleted the linked list of string
                      // has already been defined in the typedef struct
                      // for ConfigListNode
     curr = next;
  }

  theList = NULL;
} // all done


// Given "name", find all the elements of the ConfigList that
//  have have that "name"
StringList *ConfigList::find(const char *name) const
{
  ConfigListNode *tmpList;
  tmpList = find_key_word(name);
  if (tmpList != NULL)
    return tmpList -> data;
  return NULL;
}

//#define TEST_CONFIGLIST_C
#ifdef TEST_CONFIGLIST_C
#include <iostream.h>

main()
{
  ConfigList dat("namd.conf"); // an example file
  StringList *strings;
  int i;
  
  if (!dat.okay())
    NAMD_die("Can't continue - cannot get info from file.");
  std::cout << "Searching for: 'fullname'\n";
  strings = dat.find("fullname");
  if (!strings) {
    NAMD_die("Couldn't find fullname.\n");
  }
  i=0;
  while (strings) {   // write all the data associated with 'fullname'
    std::cout << i << " = " << strings->data << '\n';
    strings = strings -> next;
    i++;
  }
}

#endif

