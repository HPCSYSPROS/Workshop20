/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"

#include "GlobalMaster.h"
#include "GlobalMasterEasy.h"

#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include <stdio.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"


GlobalMasterEasy::GlobalMasterEasy(const char *the_config_name) {
  DebugM(1,"Here\n");

  molecule = Node::Object()->molecule;
  int len = strlen(the_config_name);
  configName = new char[len+1];
  strcpy(configName,the_config_name);

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  initialize();
}

GlobalMasterEasy::~GlobalMasterEasy() {
  delete [] configName;
  delete reduction;
}

int GlobalMasterEasy::
        getAtomID(const char *segid, int resid, const char *aname)
{
  return molecule->get_atom_from_name(segid,resid,aname);
}

int GlobalMasterEasy::
        getNumAtoms(const char* segid, int resid) // 0 on error
{
  return molecule->get_residue_size(segid,resid);
}

int GlobalMasterEasy::
        getAtomID(const char *segid, int resid, int index)
{
  return molecule->get_atom_from_index_in_residue(segid,resid,index);
}

double GlobalMasterEasy::getMass(int atomid)
{
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1.;  // failure
  return molecule->atommass(atomid);
}


int GlobalMasterEasy::requestAtom(int atomid)
{
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  modifyRequestedAtoms().add(atomid);
  return 0;  // success
}

int GlobalMasterEasy::getPosition(int atomid, Position &position)
{
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    if ( *a_i == atomid ) {
      position = *p_i;
      return 0;  // success
    }
  }
  return -1;  // failure
}

int GlobalMasterEasy::addForce(int atomid, Force force)
{
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  modifyForcedAtoms().add(atomid);
  modifyAppliedForces().add(force);
  return 0;  // success
}

void GlobalMasterEasy::addEnergy(BigReal e) {
  energy += e;
}

void GlobalMasterEasy::initialize() {
  DebugM(1,"Here\n");
  DebugM(4,"Initializing master\n");

  // Get the script for subclasses
  StringList *script = Node::Object()->configList->find(configName);

  config = new char[1];
  config[0] = '\0';

  for ( ; script; script = script->next) {
    if ( strstr(script->data,"\n") ) {
      size_t add_len = strlen(script->data);
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 2];
      strcpy(new_config,config);
      strcat(new_config,script->data);
      strcat(new_config,"\n");  // just to be safe
      delete [] config;
      config = new_config;
    } else {
      FILE *infile = fopen(script->data,"r");
      if ( ! infile ) {
        char errmsg[256];
        sprintf(errmsg,"Error trying to read file %s!\n",script->data);
        NAMD_die(errmsg);
      }
      fseek(infile,0,SEEK_END);
      size_t add_len = ftell(infile);
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 3];
      strcpy(new_config,config);
      delete [] config;
      config = new_config;
      new_config += config_len;
      rewind(infile);
      fread(new_config,sizeof(char),add_len,infile);
      new_config += add_len;
      new_config[0] = '\n';
      new_config[1] = '\0';
      fclose(infile);
    }
  }

}

void GlobalMasterEasy::calculate() {
  DebugM(1,"Here\n");

  /* zero out the forces */
  modifyForcedAtoms().resize(0);
  modifyAppliedForces().resize(0);

  /* XXX is this necessary? */
  modifyGroupForces().resize(getGroupMassEnd() - getGroupMassBegin());
  modifyGroupForces().setall(Vector(0,0,0));
  energy = 0.0;

  // Build results here
  easy_calc();

  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  reduction->submit();
}

void GlobalMasterEasy::easy_init(const char *) {
  CkPrintf("Default GlobalMasterEasy::easy_init() called.\n");
}

void GlobalMasterEasy::easy_calc() {
  CkPrintf("Default GlobalMasterEasy::easy_calc() called.\n");
}
