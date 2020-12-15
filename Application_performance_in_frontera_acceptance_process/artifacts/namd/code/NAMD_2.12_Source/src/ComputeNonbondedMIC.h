#include "ComputeNonbondedUtil.h"
#include "ComputeHomeTuples.h"
#include "ComputeNonbondedMICKernel.h"
#include <set>

class ComputeMgr;

class ComputeNonbondedMICKernel;

class float4;
class double4;

int mic_device_pe();

bool mic_device_shared_with_pe(int pe);

class ComputeNonbondedMIC : public Compute, private ComputeNonbondedUtil {
  public:

  struct compute_record {
    ComputeID c;
    PatchID pid[2];
    Vector offset;
    int isSelf;
    int part;
    int numParts;
  };

  struct patch_record {
    int localIndex;
    int localStart;
    int numAtoms;
    int numFreeAtoms;
    int refCount;
    int isLocal;
    int hostPe;
    PatchID patchID;
    Patch *p;
    Box<Patch,CompAtom> *positionBox;
    Box<Patch,Results> *forceBox;
    Box<Patch,Real>   *intRadBox; //5 GBIS Boxes
    Box<Patch,GBReal> *psiSumBox;
    Box<Patch,Real>   *bornRadBox;
    Box<Patch,GBReal> *dEdaSumBox;
    Box<Patch,Real>   *dHdrPrefixBox;
    CompAtom *x;
    CompAtomExt *xExt;
    Results *r;
    Force *f;
    Real   *intRad; //5 GBIS arrays
    GBReal *psiSum;
    Real   *bornRad;
    GBReal *dEdaSum;
    Real   *dHdrPrefix;

    patch_record() { refCount = 0; }
  };


    ComputeNonbondedMIC(ComputeID c, ComputeMgr *mgr,
		ComputeNonbondedMIC *m = 0, int idx = -1);
    ~ComputeNonbondedMIC();

    void atomUpdate();
    void doWork();
    int noWork();
    void skip();

    void recvYieldDevice(int pe);
    LocalWorkMsg *localWorkMsg2;

    int workStarted;
    Lattice lattice;
    int doSlow, doEnergy;
    int step;
    int finishWork();  // returns true when finished, false to continue

    static void bind_lj_table(int deviceNum);
    static void bind_force_table(int deviceNum);
    static void bind_constants(int deviceNum);
    static void bind_exclusions(int deviceNum);

    void build_exclusions();

    void requirePatch(int pid);
    void assignPatches();
    void registerPatches();
    ResizeArray<int> activePatches, localActivePatches, remoteActivePatches;
    ResizeArray<int> hostedPatches, localHostedPatches, remoteHostedPatches;
    ResizeArray<patch_record> patchRecords;
    ResizeArray<compute_record> computeRecords;
    ResizeArray<compute_record> localComputeRecords, remoteComputeRecords;

    int num_atom_records;
    int num_local_atom_records;
    int num_remote_atom_records;
    int num_force_records;

    int num_local_compute_records;
    int num_remote_compute_records;
    int num_local_patch_records;
    int num_remote_patch_records;

    double4 *forces;
    double4 *slow_forces;
    GBReal *psiSumH;
    GBReal *dEdaSumH;

    PatchMap *patchMap;
    AtomMap *atomMap;
    SubmitReduction *reduction;

    ComputeNonbondedMICKernel *kernel;

    ComputeNonbondedMIC *master;
    int masterPe;
    int slaveIndex;
    ComputeNonbondedMIC **slaves;
    int *slavePes;
    int numSlaves;

    int atomsChanged;
    int computesChanged;

    int pairlistsValid;
    float pairlistTolerance;
    int usePairlists;
    int savePairlists;
    float plcutoff2;

    #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
      int micDevice;
      int2 *exclusionsByAtom_ptr;
      std::set<void*> *atomSubmitSignals;
      virtual void patchReady(PatchID patchID, int doneMigration, int seq);
    #endif

    // DMK - NOTE : Moved some variables from global scope in ComputeNonbondedMIC to here
    // DMK - NOTE : For the following members, non-MIC builds will not have the types defined, but
    //   this class should also go unused.  So leave these members out of non-MIC builds.  May
    //   need to move them back to "__thread" variables or find a cleaner solution.
    #if defined(NAMD_MIC)
      ResizeArray<patch_pair> patch_pairs;
      ResizeArray<force_list> force_lists;
      atom * atoms;
      atom_param * atom_params;
    #endif

    int exclContrib;

    // DMK - DEBUG
    int timestep;
};


