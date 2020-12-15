/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Read in a configuration file of the form:
         keyword = information\n
   -or-  keyword information\n
   -or-  keyword = {\n line 0\n line 1\n ... line n\n}
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

// This header introduces two names to the global name space
// They are:
//    StringList -- a linked list of char *.  It has an initilizer and
//      destructor to handle the char * by alloc'ing new space
//    ConfigList -- its constructor opens a file containing lines of the
//      format "keyword = data" and creates a listing of all the data
//      for each keyword.  This data can be retreived with "find(char *)"

#ifndef CONFIGLIST_H
#define CONFIGLIST_H
#include "common.h"
#include <string.h>

class StringList {
 public:
  char *data;
  StringList *next;
  StringList(char *newdata) {  // take a string, and copy it
     data = new char[strlen(newdata)+1];
     if ( data == NULL )
     {
       NAMD_die("new failed in struct StringList");
     }
     strcpy( data, newdata);
     next = NULL;
  }
  void set(const char *newdata) {  // take a string, and copy it
    delete [] data;
    data = new char[strlen(newdata)+1];
    if ( data == NULL )
    {
      NAMD_die("new failed in struct StringList");
    }
    strcpy( data, newdata);
  }
  ~StringList( void) {  // just clear out my info
    delete [] data;
    data = NULL;
    next = NULL;
  }
};

class ConfigList {
  public:
    class ConfigListNode 
    {
    private:
    public:
      char *name;
      StringList *data;
      ConfigListNode *next;
      ConfigListNode( ConfigListNode *newnext, char *newname, 
                                          StringList *newdata) {
        name = new char[strlen(newname)+1];  // create space for the name
	if ( name == NULL )
	{
	  NAMD_die("new failed in ConfigListNode::ConfigListNode");
	}
        strcpy((char *) name, newname);      // and copy the new name
        data = newdata;
        next = newnext;
      }
      ~ConfigListNode( void) 
      {
        delete [] name;                  // delete the new'ed name
        name = NULL;
        next = NULL;
        StringList *curr=data, *next_local=NULL;  // go through the string list

        while ( curr!=NULL ) 
	{
          next_local = curr->next;           // and delete each element
          delete curr;
	  curr = next_local;
        }
      }
    };
 private:
    ConfigListNode *theList;
       // copy the information into a String, as appropriate
       // this is really a "push"
    ConfigListNode *find_key_word( const char *keyword) const;
    Bool isokay;
  public:
    ConfigList(void);
    void add_element( const char *s1, int len1, const char *s2, int len2);
    ConfigList( const char *filename);
    Bool okay( void) { return isokay; }
    ~ConfigList( void);
    StringList *find( const char *name) const;   //search for values by name
       // NOTE: do not change or delete this information.  It is not new'ed
       //   and any changed you make will be permanent.
       
    ConfigListNode *head( void) const { return theList;  } // return everything
       // NOTE:  you _REALLY_ not not want to change the information
       //   you get from this (unless you really want to)
};


#endif // CONFIGLIST_H

