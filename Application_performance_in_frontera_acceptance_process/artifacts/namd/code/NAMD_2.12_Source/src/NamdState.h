/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef _NAMDSTATE_H
#define _NAMDSTATE_H

#include "Lattice.h"

class Molecule;
class SimParameters;
class Parameters;
class ConfigList;
class PDB;
class Controller;

// Everything needed to specify a simulation is in this object
// For the moment it is really only a structure.  Eventually
// I hope it encapsulates a namd state.  
class NamdState {
  friend class Namd; 
  friend class Node;
  friend class Controller;
  private:
    Molecule *molecule;
    Parameters *parameters;
    SimParameters *simParameters;
    ConfigList *configList;
    PDB *pdb;
    Controller *controller;
    Lattice lattice;
    //char *currentdir;
    std::string callback_labelstring, callback_valuestring;

#ifdef MEM_OPT_VERSION
    void checkMemOptCompatibility();
#endif

public:
    NamdState(void);
    ~NamdState() {}
    int configFileInit(char *);
    friend class ScriptTcl;
    int configListInit(ConfigList *);
    int loadStructure(const char *, const char *, int);
    int status();
    void useController(Controller *controllerPtr);
    void runController(void);
};

#endif /* _NAMDSTATE_H */
 
