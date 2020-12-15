/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Uses simple hash table type unique elements set 
   layout allows for faster iteration due to clustering of data
*/

#ifndef USETRAW_H
#define USETRAW_H

#include <new>

template <class Elem> class EntryGlob;
template <class Elem> class UniqueSetRaw;

// Utility class for UniqueSetRaw
// Each entry hold an Obj<Elem>, used flag and next pointer

template <class Elem> class Entry {

  friend class UniqueSetRaw<Elem>;
  friend class EntryGlob<Elem>;

  private:
  
    Entry<Elem> *_next;
    long used;
    Entry(void) { _next = 0; used = 0; }
  
  public:
  
    Elem obj;

    ~Entry(void) { clear(); }

    void clear(void) { 
      if (used) obj.~Elem();
      used = 0;
    }

    Entry<Elem> *next(void) const { return _next; }

    int isUsed(void) const { return used; }
  
    static void clearChain(Entry<Elem> *e) {
      while(e) {
        Entry<Elem> *tmp = e->next();
        e->clear();
        e = tmp;
      }
    }

    // remove this entry from link list
    Entry<Elem> *removeSelf(void) {
      Entry<Elem> *e = this->next();
      clear();
      return (e);
    }

    //
    Entry<Elem> &operator=(const Elem &e) {
      if (!used)
        new(&obj) Elem;
      obj = e;
      used = 1;
      return (*this);
    }

    //
    Entry<Elem> *setNext(Entry<Elem> *e) {
      this->_next = e;
      return this;
    }
};

template <class Elem> class UniqueSetRaw;
template <class Elem> class UniqueSetIter;

// Utility class for UniqueSetRaw
template <class Elem> class EntryGlob {

  friend class UniqueSetIter<Elem>;
  friend class UniqueSetRaw<Elem>;

  public:
  
    ~EntryGlob(void) { clear(); }

    void clear(void) { delete[] glob; }

    EntryGlob<Elem> *next(void) { return _next; }

    EntryGlob<Elem> *setNext(EntryGlob<Elem>* e) { 
      this->_next = e; return this; 
    }

    static void clearChain(EntryGlob<Elem> *e) {
      while(e) { EntryGlob<Elem> *tmp = e->next(); delete e; e = tmp; }
    }

  private:

    inline EntryGlob(int size);
  
    EntryGlob<Elem> *_next;
    Entry<Elem> *glob;
};
  
template<class Elem>
inline EntryGlob<Elem>::EntryGlob(int size) {
  _next = 0;
  Entry<Elem>* entry = glob = new Entry<Elem>[size];
  for (; entry != glob+size-1; entry++) {
    entry->setNext(entry+1);
  }
  entry->setNext(0);
}


//----------------------------------------------------------------------
// UniqueSetRaw container class definition
// Uses hash table type lookup
// main user methods: add(), del(), find()
//----------------------------------------------------------------------

template <class Elem> class UniqueSet;
template <class Elem> class UniqueSetIter;

template <class Elem> class UniqueSetRaw {

  friend class UniqueSet<Elem>;
  friend class UniqueSetIter<Elem>;

  public:

    // Various Constructors
    inline UniqueSetRaw(int size=0);

    // Hopefully, this is never used
    UniqueSetRaw(const UniqueSetRaw<Elem> &us) { copy(us); }

    ~UniqueSetRaw(void) { cleanUp(); }

    UniqueSetRaw<Elem> & operator=(const UniqueSetRaw<Elem> &us) {
      copy(us);
      return *this;
    }  

    void rehash(void) { rehash(numElem); }

    // make a new table of appropriate size
    void rehash(int size) {
      // Store away old pointers and Entry<Elem>(s)
      int oldTableLength= tableLength;
      Entry<Elem> **oldTable = table;
    
      // recreate the table
      tableLength = findSize(numElem > size ? numElem : size);
      if (tableLength == oldTableLength) 
        return;
      numElem = 0;
      table = new Entry<Elem>*[tableLength];
      growSize = tableLength * 4;
    
      // 0 table
      Entry<Elem> **t = table; Entry<Elem> **e = table+tableLength;
      for (; t != e; *t++ = 0);
    
      // go thru old table and chains of Entry<Elem>(s)
      for (int i=0; i<oldTableLength; i++) {
        Entry<Elem> *e = oldTable[i];
        while(e) {
          Entry<Elem> *tmp = e->next();
          remapEntry(e);
          numElem++;
          e = tmp;
        }
      }
      delete[] oldTable;
    } 

  
    void setSafeNoDestruction(int flag) { isSafeNoDestruction = flag; }
  
    // add element to our hash table!
    int add(const Elem &elem) {
      int tableSlot = elem.hash() % tableLength;
      int doadd = 1;
      for (Entry<Elem> *e = table[tableSlot]; e; e=e->next()) {
        if (e->obj == elem) { doadd = 0; break; }
      }
      if (doadd) {
        Entry<Elem> *entry = nextFree();
        *entry = elem; // sets used flag
        table[tableSlot] = entry->setNext(table[tableSlot]);
        if ( ++numElem > growSize)  
	  rehash();
      }
      return(doadd);
    }
  
    int load(const Elem &elem) {
      Entry<Elem> *entry = nextFree();
      *entry = elem; // sets used flag
      table[0] = entry->setNext(table[0]);
      ++numElem;
      return(1);
    }
  
    // find element
    Elem *find(const Elem &elem) {
      register Entry<Elem> *e = table[elem.hash() % tableLength];
      while(e) {
        if (elem == e->obj ) {
	  return (&(e->obj));
        }
        e = e->next();
      }
      return (0);
    }
  
    int del(const Elem &elem) {
      int tableSlot = elem.hash() % tableLength;
      int dodel = 0;
    
      Entry<Elem> *prev = 0;
      Entry<Elem> *e = table[tableSlot];
      while(e) {
        if (e->obj == elem) {
          if (0 == prev) {
	    table[tableSlot] = e->removeSelf();
          } else {
	    prev->setNext(e->removeSelf());
          } 
          addFree(e); dodel = 1; numElem--;
          break;
        }
        prev = e;
        e = e->next();
      }
      return dodel;
    }

  
    int size(void) const { return(numElem); }

    void clear(int n=0) {
      if (!n) n = size();
      cleanUp();
      init(n); // have to have something currently for a table
    }

#ifdef DEBUG
    void status(void) {
      int count = 0;
      for (int i=0; i < tableLength; i++) {
        Entry<Elem> *e = table[i];
        cout << "Table entry [" << i << "]" << std::endl;
        while (e) {
          if (e->isUsed())  { 
	      count++;
	      e->obj.status();
          } else {
	      cout << "Entry is not used" << std::endl;
          }
          e = e->next();
        }
      }
      cout << "===== COUNT = " << count << std::endl;
    }
#endif
  
  protected:

    int refCount;
    // used by iterator and free list
    EntryGlob<Elem> *globHead;
  
  private:

    // basic hash table pointers
    Entry<Elem> **table;
    int tableLength;
    int growSize;
    // number of elements in table
    int numElem;
    // allocation parameters
    int growable;
    int globSize;
    int isSafeNoDestruction;

    // Utilities

    int findSize(int size) {
      size /= 2;
      int rval;
      if ( ( rval = 11 ) <= size )
       if ( ( rval = 31 ) <= size )
        if ( ( rval = 101 ) <= size )
         if ( ( rval = 307 ) <= size )
          if ( ( rval = 1009 ) <= size )
           if ( ( rval = 3001 ) <= size )
            if ( ( rval = 10007 ) <= size )
             if ( ( rval = 30011 ) <= size )
              if ( ( rval = 100003 ) <= size )
               if ( ( rval = 300007 ) <= size )
                if ( ( rval = 1000003 ) <= size )
                 if ( ( rval = 3000017 ) <= size )
                  if ( ( rval = 10000019 ) <= size )
                   if ( ( rval = 30000023 ) <= size )
                    if ( ( rval = 100000007 ) <= size )
                     if ( ( rval = 300000007 ) <= size )
                      rval = 1000000007;
      return rval;
    }

    void cleanUp(void) {
      for (int i=0; i<tableLength; i++)
        Entry<Elem>::clearChain(table[i]);
      tableLength = 0;
      growSize = 0;
      delete[] table;
      table = 0;
      freeListClean();
      numElem = 0;
    }

    void copy(const UniqueSetRaw<Elem> &us) {
      cleanUp();
      isSafeNoDestruction = us.isSafeNoDestruction;
      tableLength = findSize(us.numElem);
      growSize = tableLength * 4;
      table = new Entry<Elem>*[tableLength];
    
      for (int i=0; i<us.tableLength; i++) {
        Entry<Elem> *e = us.table[i];
        while (e) {
          add(e->obj);
          e = e->next();
        }
      }
    }

    void init(int size) {
      tableLength = findSize(size);
      growSize = tableLength * 4;
      table = new Entry<Elem>*[tableLength];
    
      Entry<Elem> **t = table; Entry<Elem> **e = table+tableLength;
      for (; t != e; *t++ = 0);
    
      freeListInit();
    }

    // Only to be used by rehash
    void remapEntry(Entry<Elem> *e) {
      int tableSlot = e->obj.hash() % tableLength;
      table[tableSlot] = e->setNext(table[tableSlot]);
    }
  
    // free list maintenance
    Entry<Elem> *head;

    void freeListInit(void) {
      head = 0; globHead = 0;
    }

    void freeListClean(void) {
      EntryGlob<Elem> *tmp = globHead;
      while( (globHead = tmp) ) {
        tmp = globHead->next();
        delete globHead; // delete's glob
      } 
    }

    Entry<Elem> *nextFree(void) {
      if (!head) {
         EntryGlob<Elem> *tmp = new EntryGlob<Elem>(globSize);
         head = tmp->glob;
         globHead = tmp->setNext(globHead);
      }
      Entry<Elem> *entry = head;
      head = entry->next();
      return(entry);
    }

    void addFree(Entry<Elem> *e) {
      head = e->setNext(head);
      e->used = 0;
    }
  
};

template <class  Elem>
inline UniqueSetRaw<Elem>::UniqueSetRaw(int size) : 
  table(0), tableLength(0), numElem(0), growable(1) {
  const int minGlobSize = 32;
  isSafeNoDestruction = 1;
  growSize = 0;
  globSize = minGlobSize;
  init(size);
}

#endif
