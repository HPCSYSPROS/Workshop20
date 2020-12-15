#ifndef COMPRESSPSF_H
#define COMPRESSPSF_H

#include "structures.h"
#include <string>
#include <deque>
#include <ckhashtable.h>

#define COMPRESSED_PSF_VER 1.72

//used to detemine big-endian or little-endian for 
//the per-atom binary file
#define COMPRESSED_PSF_MAGICNUM 1234

class Molecule;
class Parameters;
class SimParameters;
class ConfigList;

//if the compiler supports anonymous struct, then sSet, iSet
//and fSet could be omitted for simplicity. -Chao Mei
struct OutputAtomRecord{
  struct shortVals{
    short segNameIdx;
    short resNameIdx;
    short atomNameIdx;
    short atomTypeIdx;
    short chargeIdx;
    short massIdx;    
    short vdw_type;
  }sSet;

  struct integerVals{
    int atomSigIdx;
    int exclSigIdx;
    int resID;      
    int hydrogenList;
    int atomsInGroup;
    int GPID;      
    int atomsInMigrationGroup;
    int MPID;    
  }iSet;

  struct floatVals{
    Real rigidBondLength;
  }fSet;
           
  void flip();
};

void compress_molecule_info(Molecule *mol, char *psfFileName, Parameters *param, SimParameters *simParam, ConfigList* cfgList);

void flipNum(char *elem, int elemSize, int numElems);

template <typename T>
int lookupCstPool(const std::vector<T>& pool, const T& val)
{
    for(int i=0; i<pool.size(); i++)
    {
        if(pool[i]==val)
            return i;
    }
    return -1;
}

// Adapt the function prototype for ckhashtable keys to the 
// NAMD XXXSig classes
template <class T> class HashPoolAdaptorT {
  T val;
public:
  HashPoolAdaptorT<T>(const T &v):val(v) { }
  /**added to allow pup to do Key k while unPacking*/
  HashPoolAdaptorT<T>(){}
  operator T & () { return val; }
  operator const T & () const { return val; }
  
  inline CkHashCode hash(void) const { 
    const int hash=val.hash();
    return CkHashFunction_int(&hash,sizeof(int));
  }

  static CkHashCode staticHash(const void *k,size_t) {
    return ((HashPoolAdaptorT<T> *)k)->hash();
  }
  
  inline int compare(const HashPoolAdaptorT<T> &t) const {
    return val==t.val;
  }
  
  static int staticCompare(const void *a,const void *b,size_t) {
    return ((HashPoolAdaptorT<T> *)a)->compare(*(HashPoolAdaptorT<T> *)b);
  }
  
  inline T& getVal() { return val; }
  
  // PUPer not tested
  void pup(PUP::er &p){
    p | *val;
  }
};

template <typename T>
class HashPool {
public:
  ~HashPool() {
    clear();
  }
  
  void clear() {
    // Delete the pool entries
    for ( int i=0; i < pool.size(); i++)
      delete pool[i];
    // Clear the pool and hash table
    pool.clear();
    index_table.empty();
  }
  
  int lookupCstPool(const T& val)
  {
    HashPoolAdaptorT<T> hx(val);
    // Ugly: Can't store zeros in the table, so add 1 to indices on insert and
    // subtract 1 to get real index when retrieving
    int loc = index_table.get(hx) - 1;
#if 0
    int i;
    for(i=0; i < pool.size(); i++)
    {
      if (pool[i]->getVal() == val)  {
        if (i != loc) {
          CmiPrintf("Get[%d] returned %d, actual is %d\n",hx.hash(),loc,i);
          dump_tables();
          loc = i;
        }
        break;
      }
    }
    if (loc != -1 && i == pool.size()) {
      CmiPrintf("Get returned %d, actual not found\n",loc);
    }
#endif
    return loc;
  }
  
  void push_back(const T& x)
  {
    // Behave like a STL vector, but also store the indexing info in
    // the hashtable index
    
    // Add to vector
    HashPoolAdaptorT<T>* val = new HashPoolAdaptorT<T>(x);
    pool.push_back(val);
    // Also add to hash table. Make sure pointer doesn't change
    int* index = &(index_table.put(*val));
    // Ugly: Can't store zeros in the table, so add one to all entries
    *index = pool.size(); // pool.size()  - 1 + 1
    
    //  CmiPrintf("Adding hx=%p hash=%d index[%p]=%d\n",&val,val->hash(),index,*index);
    //  dump_tables();  
  }
  
  void dump_tables(); 
  
  T& operator[](int i) const { return pool[i]->getVal(); }
  
  int size() { return pool.size(); }

private:
  CkHashtableT<HashPoolAdaptorT<T>,int> index_table;
  std::vector<HashPoolAdaptorT<T>*> pool;
};

#endif
