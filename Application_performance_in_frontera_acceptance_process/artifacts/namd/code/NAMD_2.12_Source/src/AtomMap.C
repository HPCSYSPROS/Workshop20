/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Tracks location of Atoms on node.  Singleton.
*/

#include "ProcessorPrivate.h"
#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

#ifdef MEM_OPT_VERSION
#define MAXBITS 20
#define MAXNUMATOMS (1<<MAXBITS)
#endif


void AtomMapper::registerIDsCompAtomExt(const CompAtomExt *begin, const CompAtomExt *end) {
  if ( mapped ) return;
  mapped = 1;
#ifdef MEM_OPT_VERSION
  if ( ! map->onlyUseTbl ) {
    int n = end - begin;
    entries.resize(n);
    AtomMapEntry *e = entries.begin();
    for ( int i=0; i<n; ++i, ++e ) {
      e->pid = pid;
      e->index = i;
      AtomID aid = begin[i].id;
      short aid_upper = aid >> MAXBITS;
      e->aid_upper = aid_upper;
      int aid_hash = aid & (MAXNUMATOMS-1);
      AtomMapEntry **me = map->entries + aid_hash;
      while ( *me && (*me)->aid_upper < aid_upper ) me = &((*me)->next);
      e->next = *me;
      *me = e;
    }
  } else
#endif
  if ( map->registerIDsCompAtomExt(pid, begin, end) ) NAMD_bug("atom map failed");
}


void AtomMapper::registerIDsFullAtom(const FullAtom *begin, const FullAtom *end) {
  if ( mapped ) return;
  mapped = 1;
#ifdef MEM_OPT_VERSION
  if ( ! map->onlyUseTbl ) {
    int n = end - begin;
    entries.resize(n);
    AtomMapEntry *e = entries.begin();
    for ( int i=0; i<n; ++i, ++e ) {
      e->pid = pid;
      e->index = i;
      AtomID aid = begin[i].id;
      short aid_upper = aid >> MAXBITS;
      e->aid_upper = aid_upper;
      int aid_hash = aid & (MAXNUMATOMS-1);
      AtomMapEntry **me = map->entries + aid_hash;
      while ( *me && (*me)->aid_upper < aid_upper ) me = &((*me)->next);
      e->next = *me;
      *me = e;
    }
  } else
#endif
  if ( map->registerIDsFullAtom(pid, begin, end) ) NAMD_bug("atom map failed");
}


void AtomMapper::unregisterIDsCompAtomExt(const CompAtomExt *begin, const CompAtomExt *end) {
  if ( ! mapped ) return;
  mapped = 0;
#ifdef MEM_OPT_VERSION
  if ( ! map->onlyUseTbl ) {
    int n = end - begin;
    if ( entries.size() != n ) {
      CkPrintf("AtomMapper entries.size() %d != %d\n", entries.size(), n);
      NAMD_bug("AtomMapper::unregisterIDsCompAtomExt size mismatch");
    }
    AtomMapEntry *e = entries.begin();
    for ( int i=0; i<n; ++i, ++e ) {
      AtomID aid = begin[i].id;
      int aid_hash = aid & (MAXNUMATOMS-1);
      AtomMapEntry **me = map->entries + aid_hash;
      while ( *me != e ) me = &((*me)->next);
      *me = e->next;
    }
  } else
#endif
  if ( map->unregisterIDsCompAtomExt(pid, begin, end) ) NAMD_bug("atom map failed");
}


void AtomMapper::unregisterIDsFullAtom(const FullAtom *begin, const FullAtom *end) {
  if ( ! mapped ) return;
  mapped = 0;
#ifdef MEM_OPT_VERSION
  if ( ! map->onlyUseTbl ) {
    int n = end - begin;
    if ( entries.size() != n ) {
      CkPrintf("AtomMapper entries.size() %d != %d\n", entries.size(), n);
      NAMD_bug("AtomMapper::unregisterIDsFullAtom size mismatch");
    }
    AtomMapEntry *e = entries.begin();
    for ( int i=0; i<n; ++i, ++e ) {
      AtomID aid = begin[i].id;
      int aid_hash = aid & (MAXNUMATOMS-1);
      AtomMapEntry **me = map->entries + aid_hash;
      while ( *me != e ) me = &((*me)->next);
      *me = e->next;
    }
  } else
#endif
  if ( map->unregisterIDsFullAtom(pid, begin, end) ) NAMD_bug("atom map failed");
}


// Singleton method
AtomMap *AtomMap::Instance() {
  if (CkpvAccess(AtomMap_instance) == 0) {
    CkpvAccess(AtomMap_instance) = new AtomMap;	// this is never deleted!
  }
  return CkpvAccess(AtomMap_instance);
}

//----------------------------------------------------------------------
AtomMap::AtomMap(void)
{
  localIDTable = NULL;
  tableSz = 0;

#ifdef MEM_OPT_VERSION
  entries = NULL;
  onlyUseTbl = false;
#endif
}

void
AtomMap::checkMap(void)
{ }
  

//----------------------------------------------------------------------
AtomMap::~AtomMap(void)
{
  delete [] localIDTable;  // Delete on a NULL pointer should be ok

#ifdef MEM_OPT_VERSION
  delete [] entries;
#endif
}

//----------------------------------------------------------------------
// Creates fixed size table
void AtomMap::allocateMap(int nAtomIds)
{
#ifdef MEM_OPT_VERSION
  if ( nAtomIds > MAXNUMATOMS ) {
    entries = new AtomMapEntry*[MAXNUMATOMS];
    memset(entries,0,MAXNUMATOMS*sizeof(AtomMapEntry*));
    return;
  } // else use non-memopt strategy
  onlyUseTbl = true;	
#endif
  if ( nAtomIds <= tableSz ) return;
  LocalID *oldTable = localIDTable;
  localIDTable = new LocalID[nAtomIds];
  for(int i=0; i < tableSz; i++)
    localIDTable[i] = oldTable[i];
  for(int i=tableSz; i < nAtomIds; i++)
    localIDTable[i].pid = localIDTable[i].index = notUsed;
  delete [] oldTable;
  tableSz = nAtomIds;
}

//
int AtomMap::unregisterIDsCompAtomExt(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const CompAtomExt *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	if (localIDTable[ali].pid == pid) {
	    localIDTable[ali].pid = notUsed;
	    localIDTable[ali].index = notUsed;
	}
    }
    return 0;
  }
}

//----------------------------------------------------------------------
int AtomMap::unregisterIDsFullAtom(PatchID pid, const FullAtom *begin, const FullAtom *end)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const FullAtom *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	if (localIDTable[ali].pid == pid) {
	    localIDTable[ali].pid = notUsed;
	    localIDTable[ali].index = notUsed;
	}
    }
    return 0;
  }
}

//It's possible to register the same atom for a new patch before it is moved 
//from the old patch on the same processor!

//----------------------------------------------------------------------
int AtomMap::registerIDsCompAtomExt(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const CompAtomExt *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	localIDTable[ali].pid = pid;
	localIDTable[ali].index = a - begin;
    }
    return 0;
  }
}

//----------------------------------------------------------------------
int AtomMap::registerIDsFullAtom(PatchID pid, const FullAtom *begin, const FullAtom *end)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const FullAtom *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	localIDTable[ali].pid = pid;
	localIDTable[ali].index = a - begin;
    }
    return 0;
  }
}


#ifdef MEM_OPT_VERSION
LocalID AtomMap::localID(AtomID id)
{
	if(onlyUseTbl){
		return localIDTable[id];
	}else{

      short aid_upper = id >> MAXBITS;
      int aid_hash = id & (MAXNUMATOMS-1);
      AtomMapEntry *me = entries[aid_hash];
      while ( me && me->aid_upper < aid_upper ) me = me->next;
      LocalID rval;
      if ( me && me->aid_upper == aid_upper ) {
        rval.pid = me->pid;
        rval.index = me->index;
      } else {
        rval.pid = notUsed;
        rval.index = notUsed;
      }
      return rval;
	}
}
#endif

