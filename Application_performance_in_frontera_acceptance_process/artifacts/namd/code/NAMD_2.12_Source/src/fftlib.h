
#ifndef   __FFT_LIB_H__
#define   __FFT_LIB_H__

#include <cmidirectmanytomany.h>

#define  MANY_TO_MANY_START  10
#define  MANY_TO_MANY_SETUP  2

#define  PHASE_YF  0   // Recv Y-Forward
#define  PHASE_XF  1   // Recv X-Forward
#define  PHASE_YB  2   // Recv Y-Backward
#define  PHASE_ZB  3   // Recv Z-Backward
#define  PHASE_GR  4   // Send Grid
#define  PHASE_UG  5   // Recv Ungrid

extern int many_to_many_start;

class OptPmeGridMsg : public CMessage_OptPmeGridMsg {
public:
  int sourceNode;
  int xstart;
  int xlen;
  int ystart;
  int ylen;
  int zstart;
  int zlen;
  int patchID;
  float *qgrid;
};  //32 byte header on 32 bit architectures


class OptPmeFFTMsg : public CMessage_OptPmeFFTMsg {
public:
  int sourceNode;
  int nx;
  float *qgrid;
};  //12 byte header on 32 bit architectures


class OptPmeDummyMsg : public CMessage_OptPmeDummyMsg {
 public:
  int   to_pe;
};

// use this idiom since messages don't have copy constructors
struct OptPmePencilInitMsgData {
  PmeGrid grid;
  int xBlocks, yBlocks, zBlocks;
  CProxy_OptPmeXPencil xPencil;
  CProxy_OptPmeYPencil yPencil;
  CProxy_OptPmeZPencil zPencil;
  CProxy_OptPmeMgr     pmeProxy;
  CkCallback           cb_energy;
  bool                 constant_pressure;
};


class OptPmePencilInitMsg : public CMessage_OptPmePencilInitMsg {
public:
  OptPmePencilInitMsg(OptPmePencilInitMsgData &d) { data = d; }
  OptPmePencilInitMsgData data;
};

struct CkCallbackWrapper {
  CkCallback    cb;
  void        * msg;
  void        * array;
};


template <class T> class OptPmePencil : public T {
public:
  OptPmePencil() {
    data = 0;
    work = 0;
    send_order = 0;
  }
  ~OptPmePencil() {
    delete [] data;
    delete [] work;
    delete [] send_order;
  }
  void base_init(OptPmePencilInitMsg *msg) {
    initdata = msg->data;
  }
  void order_init(int nBlocks) {
    send_order = new int[nBlocks];
    for ( int i=0; i<nBlocks; ++i ) send_order[i] = i;
    Random rand(CkMyPe());
    rand.reorder(send_order,nBlocks);
  }
  OptPmePencilInitMsgData initdata;
  Lattice lattice;
  PmeReduction evir;
  int imsg;  // used in sdag code
  int _iter; // used in sdag code
  float *data;
  float *many_to_many_data; //data in a differnt format
  int   *many_to_many_nb;
  float *work;
  int *send_order;
  void *handle;
  bool single_pencil;
};


class OptPmeZPencil : public OptPmePencil<CBase_OptPmeZPencil> {
public:
    OptPmeZPencil_SDAG_CODE
    OptPmeZPencil() { __sdag_init(); setMigratable(false); }
    OptPmeZPencil(CkMigrateMessage *) { __sdag_init();  setMigratable (false); }
    void fft_init();
    void recv_grid(const OptPmeGridMsg *);
    void many_to_many_recv_grid();
    void forward_fft();
    void send_trans();
    void many_to_many_send_trans();
    void recv_untrans(const OptPmeFFTMsg *);
    void many_to_many_recv_untrans();
    void backward_fft();
    void send_ungrid(OptPmeGridMsg *);
    void many_to_many_send_ungrid ();
private:
    ResizeArray<OptPmeGridMsg *> grid_msgs;
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_plan  forward_plan, backward_plan;
#else
    rfftwnd_plan forward_plan, backward_plan;
#endif
#endif
    int nx, ny;
    PatchGridElem  * m2m_recv_grid;
    float          * many_to_many_gr_data;
    CkCallbackWrapper  cbw_recvgrid;
    CkCallbackWrapper  cbw_recvuntrans;

    void initialize_manytomany ();
};

class OptPmeYPencil : public OptPmePencil<CBase_OptPmeYPencil> {
public:
    OptPmeYPencil_SDAG_CODE
    OptPmeYPencil() { __sdag_init(); setMigratable(false); }
    OptPmeYPencil(CkMigrateMessage *) { __sdag_init(); }
    void fft_init();
    void recv_trans(const OptPmeFFTMsg *);
    void forward_fft();
    void send_trans();
    void many_to_many_send(int phase);
    void recv_untrans(const OptPmeFFTMsg *);
    void many_to_many_recv_trans();
    void backward_fft();
    void send_untrans();
    void many_to_many_recv_untrans();	
private:
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_plan  forward_plan, backward_plan;
#else
    fftw_plan forward_plan, backward_plan;
#endif
#endif
    int nx, nz;
    CkCallbackWrapper  cbw_recvtrans;
    CkCallbackWrapper  cbw_recvuntrans;

    void initialize_manytomany ();
};

class OptPmeXPencil : public OptPmePencil<CBase_OptPmeXPencil> {
public:
    OptPmeXPencil_SDAG_CODE
    OptPmeXPencil() { __sdag_init();  myKSpace = 0; setMigratable(false); }
    OptPmeXPencil(CkMigrateMessage *) { __sdag_init(); }
    void fft_init();
    void recv_trans(const OptPmeFFTMsg *);
    void many_to_many_recv_trans();	
    void forward_fft();
    void pme_kspace();
    void backward_fft();
    void send_untrans();
    void many_to_many_send_untrans();
    void submit_evir();
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_plan  forward_plan, backward_plan;
#else
    fftw_plan forward_plan, backward_plan;
#endif
#endif
    int ny, nz;
    PmeKSpace *myKSpace;
    CkCallbackWrapper  cbw_recvtrans;
    bool               constant_pressure;

    SubmitReduction *reduction;

    void initialize_manytomany ();
};

#endif
