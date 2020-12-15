
#ifndef   __PME_FFT_MAP_H__
#define   __PME_FFT_MAP_H__

#include <charm++.h>
#include <PatchMap.h>
#include <fftlib.h>

CProxy_OptPmePencilMapZ  global_map_z;
CProxy_OptPmePencilMapY  global_map_y;
CProxy_OptPmePencilMapX  global_map_x;

struct PmeFFTInfo {
  int  xBlocks;   //FFT grid dimensions
  int  yBlocks;
  int  zBlocks;    
};

static  inline void initializePmeMap(PmeFFTInfo                    _info,
				     SortableResizeArray<int>    & xprocs,
				     SortableResizeArray<int>    & yprocs,
				     SortableResizeArray<int>    & zprocs) {
    
  // decide which pes to use by bit reversal and patch use
  int i;
  int ncpus = CkNumPes();

  int *basenodes = new int [ncpus];
  memset (basenodes, 0, sizeof(int) * ncpus);
  PatchMap *pmap = PatchMap::Object();
  for (int p = 0; p < pmap->numPatches(); p++)
    basenodes[pmap->basenode(p)] = 1;
  
  // find next highest power of two
  int npow2 = 1;  int nbits = 0;
  while ( npow2 < ncpus ) { npow2 *= 2; nbits += 1; }
  
  // build bit reversal sequence
  SortableResizeArray<int> patches, nopatches, pmeprocs, baseprocs;
  i = 0;
  for ( int icpu=0; icpu<ncpus; ++icpu ) {
    int ri;
    for ( ri = ncpus; ri >= ncpus; ++i ) {
      ri = 0;
      int pow2 = 1;
      int rpow2 = npow2 / 2;
      for ( int j=0; j<nbits; ++j ) {
	ri += rpow2 * ( ( i / pow2 ) % 2 );
	pow2 *= 2;  rpow2 /= 2;
      }
    }
    // seq[icpu] = ri;
    if ( ri ) { // keep 0 for special case
      if ( pmap->numPatchesOnNode(ri) ) 
	patches.add(ri);
      else if (basenodes[ri]) 
	baseprocs.add(ri);
      else nopatches.add(ri);
    }
  }   

  delete [] basenodes;
  
  // only use zero if it eliminates overloading or has patches
  int useZero = 0;
  int npens = _info.xBlocks*_info.yBlocks;
  if ( npens % ncpus == 0 ) useZero = 1;
  if ( npens == nopatches.size() + 1 ) useZero = 1;
  npens += _info.xBlocks*_info.zBlocks;
  if ( npens % ncpus == 0 ) useZero = 1;
  if ( npens == nopatches.size() + 1 ) useZero = 1;
  npens += _info.yBlocks*_info.zBlocks;
  if ( npens % ncpus == 0 ) useZero = 1;
  if ( npens == nopatches.size() + 1 ) useZero = 1;
  
  // add nopatches then patches in reversed order
  for ( i=nopatches.size()-1; i>=0; --i ) pmeprocs.add(nopatches[i]);
  for ( i=baseprocs.size()-1; i>=0; --i ) pmeprocs.add(baseprocs[i]);

  if ( useZero && ! pmap->numPatchesOnNode(0) ) pmeprocs.add(0);
  for ( i=patches.size()-1; i>=0; --i ) pmeprocs.add(patches[i]);
  if ( pmap->numPatchesOnNode(0) ) pmeprocs.add(0);

  int pe = 0;
  int npes = pmeprocs.size();
  int nxzpes = _info.xBlocks * _info.yBlocks;
  if (nxzpes < _info.yBlocks*_info.zBlocks)
    nxzpes = _info.yBlocks*_info.zBlocks;
  
  zprocs.resize (_info.xBlocks * _info.yBlocks);
  for ( i=0; i<_info.xBlocks * _info.yBlocks; ++i, ++pe ) zprocs[i] = pmeprocs[pe%npes];
  zprocs.sort();

  pe = nxzpes; 
  yprocs.resize(_info.xBlocks*_info.zBlocks);
  for ( i=0; i<_info.xBlocks*_info.zBlocks; ++i, ++pe ) yprocs[i] = pmeprocs[pe%npes];
  yprocs.sort();
  
  xprocs.resize(_info.yBlocks*_info.zBlocks);
  //for ( i=0; i<_info.yBlocks*_info.zBlocks; ++i, ++pe ) xprocs[i] = pmeprocs[pe%npes];
  for ( i=0, pe=0; i<_info.yBlocks*_info.zBlocks; ++i, ++pe ) xprocs[i] = pmeprocs[pe%npes];
  xprocs.sort();
}


class OptPmePencilMapX : public CBase_OptPmePencilMapX
{
  PmeFFTInfo   _info;
  int        * _mapcache;
  bool         _initialized;

 public:
  OptPmePencilMapX(int xblock, int yblock, int zblock) {
    _initialized = false;
    _info.xBlocks = xblock;
    _info.yBlocks = yblock;
    _info.zBlocks = zblock;    
    global_map_x = thisProxy;   
  }
    
  inline void initialize () {
    _initialized = true;
    _mapcache = (int *) malloc(_info.yBlocks * _info.zBlocks * sizeof(int));
    
    SortableResizeArray<int>    xprocs;
    SortableResizeArray<int>    yprocs;
    SortableResizeArray<int>    zprocs;
    
    initializePmeMap (_info, xprocs, yprocs, zprocs);
    
    for (int y = 0; y < _info.yBlocks; y++) {
      for (int z = 0; z < _info.zBlocks; z ++) {
	int index = z + y * _info.zBlocks;
	int pe = xprocs[index];
	_mapcache[index] = pe;
        
	if(CkMyRank() == 0) 
	  pencilPMEProcessors[pe] = 1;
      }
    }
  }

  OptPmePencilMapX(CkMigrateMessage *m){}

  int procNum(int foo, const CkArrayIndex &idx) {
    if (!_initialized) initialize();

    CkArrayIndex3D idx3d = *(CkArrayIndex3D *) &idx;
    int index = idx3d.index[2] + idx3d.index[1] * _info.zBlocks;
    
    return _mapcache[index];
  }
};


class OptPmePencilMapY : public CBase_OptPmePencilMapY
{
  PmeFFTInfo   _info;
  int        * _mapcache;
  bool         _initialized;

 public:
  OptPmePencilMapY(int xblock, int yblock, int zblock) {
    _initialized = false;
    _info.xBlocks = xblock;
    _info.yBlocks = yblock;
    _info.zBlocks = zblock;    
    global_map_y = thisProxy;
  }

  inline void initialize() {
    _initialized = true;
    _mapcache = (int *) malloc(_info.xBlocks * _info.zBlocks * sizeof(int)); 
    
    SortableResizeArray<int>    xprocs;
    SortableResizeArray<int>    yprocs;
    SortableResizeArray<int>    zprocs;
    
    initializePmeMap (_info, xprocs, yprocs, zprocs);
    
    for (int x = 0; x < _info.xBlocks; x ++) {
      for (int z = 0; z < _info.zBlocks; z++) {
	int index = z + x * _info.zBlocks;
	int pe = yprocs[index];
	_mapcache [index] = pe;

	if (CkMyPe() == 0) 
	  pencilPMEProcessors[pe] = 1;
      }
    }
  }
  
  OptPmePencilMapY(CkMigrateMessage *m){}

  int procNum(int foo, const CkArrayIndex &idx) {
    if (!_initialized) initialize();

    CkArrayIndex3D idx3d = *(CkArrayIndex3D *) &idx;
    int index = idx3d.index[2] + idx3d.index[0] * _info.zBlocks;
    return _mapcache [index];
  }
};

class OptPmePencilMapZ : public CBase_OptPmePencilMapZ
{
  PmeFFTInfo   _info;
  int        * _mapcache;
  bool         _initialized;

 public:
  OptPmePencilMapZ(int xblock, int yblock, int zblock) {
    _initialized = false;
    _info.xBlocks = xblock;
    _info.yBlocks = yblock;
    _info.zBlocks = zblock;    
    global_map_z = thisProxy;
  }
  
  inline void initialize() {
    _initialized = true;
    _mapcache = (int *) malloc(_info.xBlocks * _info.yBlocks * sizeof(int)); 

    SortableResizeArray<int>    xprocs;
    SortableResizeArray<int>    yprocs;
    SortableResizeArray<int>    zprocs;
    
    initializePmeMap (_info, xprocs, yprocs, zprocs);

    for (int x = 0; x < _info.xBlocks; x++) {
      for (int y = 0; y < _info.yBlocks; y ++) {	
	int index = y + x * _info.yBlocks;
	int pe = zprocs[index];
	_mapcache[index] = pe;

	if (CkMyPe() == 0)
	  pencilPMEProcessors[pe] = 1;
      }
    }
  }
  
  OptPmePencilMapZ(CkMigrateMessage *m){}
  
  int procNum(int foo, const CkArrayIndex &idx) {
    if (!_initialized) initialize();
    
    CkArrayIndex3D idx3d = *(CkArrayIndex3D *) &idx;
    int index = idx3d.index[1] + idx3d.index[0] * _info.yBlocks;
    return _mapcache[index];
  }
};


#endif
