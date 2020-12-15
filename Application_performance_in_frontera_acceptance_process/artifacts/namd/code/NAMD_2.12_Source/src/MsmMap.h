/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MSMMAP_H
#define MSMMAP_H

// SSE and AVX vector intrinsics and memory alignment macros
#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
#include <emmintrin.h>  // SSE2
#if defined(__INTEL_COMPILER)
#define __align(X) __declspec(align(X) )
#elif defined(__PGI)
#define __align(X)  __attribute__((aligned(X) ))
#define MISSING_mm_cvtsd_f64
#elif defined(__GNUC__)
#define __align(X)  __attribute__((aligned(X) ))
#if (__GNUC__ < 4)
#define MISSING_mm_cvtsd_f64
#endif
#else
#define __align(X) __declspec(align(X) )
#endif
#endif

// migration of MSM computes not currently enabled
#define MSM_MIGRATION
#undef MSM_MIGRATION

// for debugging MSM migration
#define DEBUG_MSM_MIGRATE
#undef DEBUG_MSM_MIGRATE

#define MSM_MAX_BLOCK_SIZE 8
#define MSM_MAX_BLOCK_VOLUME \
  (MSM_MAX_BLOCK_SIZE * MSM_MAX_BLOCK_SIZE * MSM_MAX_BLOCK_SIZE)

#define MSM_C1VECTOR_MAX_BLOCK_SIZE (MSM_MAX_BLOCK_SIZE / 2)
#define MSM_C1VECTOR_MAX_BLOCK_VOLUME \
  (MSM_C1VECTOR_MAX_BLOCK_SIZE * \
   MSM_C1VECTOR_MAX_BLOCK_SIZE * \
   MSM_C1VECTOR_MAX_BLOCK_SIZE)

#define DEBUG_MSM
#undef DEBUG_MSM

#define DEBUG_MSM_VERBOSE
#undef DEBUG_MSM_VERBOSE

#define DEBUG_MSM_GRID
#undef DEBUG_MSM_GRID

// assert macro
#undef ASSERT
#ifdef DEBUG_MSM
#define ASSERT(expr) \
  do { \
    if ( !(expr) ) { \
      char msg[100]; \
      snprintf(msg, sizeof(msg), "ASSERT: \"%s\" " \
          "(%s, %d)\n", #expr, __FILE__, __LINE__); \
      NAMD_die(msg); \
    } \
  } while (0)
#else
#define ASSERT(expr)
#endif 


// employ mixed precision
// (but allow easy change to all double precision for comparison)
typedef float Float;
typedef double Double;


  ///////////////////////////////////////////////////////////////////////////
  //
  // Vector and matrix elements for C1 Hermite interpolation.
  //
  ///////////////////////////////////////////////////////////////////////////

  enum { C1_VECTOR_SIZE = 8, C1_MATRIX_SIZE = 8*8 };

  struct C1Vector {
    Float velem[C1_VECTOR_SIZE];
    C1Vector(Float r=0) { set(r); }
    void set(Float r) {
      for (int n=0;  n < C1_VECTOR_SIZE;  n++)  velem[n] = r;
    }
    C1Vector& operator+=(const C1Vector& v) {
      for (int n=0;  n < C1_VECTOR_SIZE;  n++)  velem[n] += v.velem[n];
      return(*this);
    }
    friend Float operator*(const C1Vector& u, const C1Vector& v) {
      Float r=0;
      for (int n=0;  n < C1_VECTOR_SIZE;  n++)  r += u.velem[n] * v.velem[n];
      return r;
    }
    friend C1Vector operator+(const C1Vector& u, const C1Vector& v) {
      C1Vector w;
      for (int n=0;  n < C1_VECTOR_SIZE;  n++) {
        w.velem[n] = u.velem[n] + v.velem[n];
      }
      return w;
    }
  };

  struct C1Matrix {
    Float melem[C1_MATRIX_SIZE];
    C1Matrix(Float r=0) { set(r); }
    void set(Float r) {
      for (int n=0;  n < C1_MATRIX_SIZE;  n++)  melem[n] = 0;
    }
    friend C1Vector operator*(const C1Matrix& m, const C1Vector& u) {
      C1Vector v;

      // XXX not tested yet
#if 1 && (defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE))
      // Hand-coded SSE2 vectorization
      // This loop requires that the single-precision input arrays be 
      // aligned on 16-byte boundaries, such that array[index % 4 == 0] 
      // can be safely accessed with aligned load/store operations
      for (int k=0, j=0;  j < C1_VECTOR_SIZE;  j++) {
        __m128 melem4 = _mm_load_ps(&m.melem[k]);
        __m128 uelem4 = _mm_load_ps(&u.velem[0]);
        __m128 tmp4 = _mm_mul_ps(melem4, uelem4); 
        melem4 = _mm_load_ps(&m.melem[k+4]);
        uelem4 = _mm_load_ps(&u.velem[4]);
        tmp4 = _mm_add_ps(tmp4, _mm_mul_ps(melem4, uelem4)); 

        // do a 4-element reduction and accumulate result
        __m128 sum4 = tmp4;
        sum4 = _mm_shuffle_ps(sum4, sum4, _MM_SHUFFLE(2, 3, 0, 1));
        sum4 = _mm_add_ps(sum4, tmp4);
        tmp4 = sum4;
        sum4 = _mm_shuffle_ps(sum4, sum4, _MM_SHUFFLE(1, 0, 3, 2));
        sum4 = _mm_add_ps(sum4, tmp4);

        // all 4 elements are now set to the sum
        float sum;
        _mm_store_ss(&sum, sum4); // store lowest element
        v.velem[j] += sum;
        k+=8;
      }
#elif 0 && (defined(__AVX__) && ! defined(NAMD_DISABLE_SSE))
      // Hand-coded AVX vectorization
      // This loop requires that the single-precision input arrays be 
      // aligned on 32-byte boundaries, such that array[index % 8 == 0] 
      // can be safely accessed with aligned load/store operations
      for (int k=0, j=0;  j < C1_VECTOR_SIZE;  j++) {
        __m256 melem8 = _mm256_load_ps(&m.melem[k]);
        __m256 uelem8 = _mm256_load_ps(&u.velem[0]);
        __m256 tmp8 = _mm256_mul_ps(melem8, uelem8); 

        // XXX this still needs to be rewritten a bit for AVX
        // do an 8-element reduction and accumulate result
        __m256 sum8 = tmp8;
        sum8 = _mm256_hadd_ps(sum8, sum8);
        sum8 = _mm256_hadd_ps(sum8, sum8);
        tmp8 = sum8;
        tmp8 = _mm256_permute2f128_ps(tmp8, tmp8, 1);
        sum8 = _mm256_hadd_ps(tmp8, sum8);

        // all 8 elements are now set to the sum
        float sum;
        _mm_store_ss(&sum, sum8); // store lowest element
        v.velem[j] += sum;
        k+=8;
      }
#else
#if defined(__INTEL_COMPILER)
#pragma vector always
#endif
      for (int k=0, j=0;  j < C1_VECTOR_SIZE;  j++) {
        for (int i = 0;  i < C1_VECTOR_SIZE;  i++, k++) {
          v.velem[j] += m.melem[k] * u.velem[i];
        }
      }
#endif
      return v;
    }
  };

  // index vector based on mixed partial derivatives in x,y,z
  enum { D000=0, D100, D010, D001, D110, D101, D011, D111 };

  // index matrix using 2 vector indexes, row-major column ordering,
  // defining partial derivatives of g(xj,yj,zj,xi,yi,zi)
#define C1INDEX(drj,dri)  ((drj)*C1_VECTOR_SIZE + (dri))


namespace msm {

  ///////////////////////////////////////////////////////////////////////////
  //
  // Resizable Array class
  //
  ///////////////////////////////////////////////////////////////////////////

  template <class T> class Array;

  template <class T>
  void swap(Array<T>& s, Array<T>& t);

  template <class T>
  class Array {
    public:
      Array() : abuffer(0), alen(0), amax(0) { }
      Array(int n) : abuffer(0), alen(0), amax(0) { resize(n); }
      Array(const Array& a) : abuffer(0), alen(0), amax(0) { copy(a); }
      ~Array() { setmax(0); }
      Array& operator=(const Array& a) {
        if (this != &a) copy(a);  // don't allow self-assignment
        return(*this);
      }
      int len() const { return alen; }
      int max() const { return amax; }
      const T& operator[](int i) const {
#ifdef DEBUG_MSM
        return elem(i);
#else
        return abuffer[i];
#endif
      }
      const T& elem(int i) const {
        if (i < 0 || i >= alen) {
          char msg[100];
          snprintf(msg, sizeof(msg), "Array index:  alen=%d, i=%d\n", alen, i);
          NAMD_die(msg);
        }
        return abuffer[i];
      }
      T& operator[](int i) {
#ifdef DEBUG_MSM
        return elem(i);
#else
        return abuffer[i];
#endif
      }
      T& elem(int i) {
        if (i < 0 || i >= alen) {
          char msg[100];
          snprintf(msg, sizeof(msg), "Array index:  alen=%d, i=%d\n", alen, i);
          NAMD_die(msg);
        }
        return abuffer[i];
      }
      void append(const T& t) {
        if (alen==amax) setmax(2*amax+1);
        abuffer[alen++] = t;
      }
      void resize(int n) {
        if (n > amax) setmax(n);
        alen = n;
      }
      void setmax(int m);
      const T *buffer() const { return abuffer; }
      T *buffer() { return abuffer; }
      const T *buffer(int& n) const { n = alen; return abuffer; }
      T *buffer(int& n) { n = alen; return abuffer; }
      void reset(const T& t) {
        for (int n = 0;  n < alen;  n++)  abuffer[n] = t;
      }
      friend void swap<T>(Array&, Array&);
#ifdef DEBUG_MSM
      void print(const char *s=0) const {
        if (s) printf("PRINTING DATA FOR ARRAY \"%s\":\n", s);
        printf("abuffer=%p\n  alen=%d  amax=%d\n",
            (void *) abuffer, alen, amax);
      }
#endif
    protected:
      T *abuffer;
      int alen, amax;
      void copy(const Array& a);
  };

  template <class T>
  void Array<T>::setmax(int m) {
    if (m == amax) return;
    else if (m > 0) {
      T *newbuffer = new T[m];
      if ( ! newbuffer) {
        char msg[100];
        snprintf(msg, sizeof(msg),
            "Can't allocate %lu KB for msm::Array\n",
            (unsigned long)(m * sizeof(T) / 1024));
        NAMD_die(msg);
      }
      if (alen > m) alen = m;  // new buffer is shorter than old buffer
      for (int i = 0;  i < alen;  i++) {
        newbuffer[i] = abuffer[i];
      }
      delete[] abuffer;
      abuffer = newbuffer;
      amax = m;
    }
    else {  // consider m == 0
      delete[] abuffer;
      abuffer = 0;
      alen = 0;
      amax = 0;
    }
  }

  template <class T>
  void Array<T>::copy(const Array<T>& a) {
    setmax(a.amax);
    alen = a.alen;
    for (int i = 0;  i < alen;  i++) {
      abuffer[i] = a.abuffer[i];
    }
  }

  // swap arrays without duplicating memory buffer
  template <class T>
  void swap(Array<T>& s, Array<T>& t) {
    T *tmpbuffer = s.abuffer;  s.abuffer = t.abuffer;  t.abuffer = tmpbuffer;
    tmpbuffer = 0;
    int tmpn = s.alen;  s.alen = t.alen;  t.alen = tmpn;
    tmpn = s.amax;  s.amax = t.amax;  t.amax = tmpn;
    tmpn = s.astate;  s.astate = t.astate;  t.astate = tmpn;
  }


  ///////////////////////////////////////////////////////////////////////////
  //
  // Priority queue for static load balancing of work
  //
  ///////////////////////////////////////////////////////////////////////////

  // smallest value has priority
  // T must have partial ordering operator<=() defined along with assignment
  template <class T>
  class PriorityQueue {
    public:
      PriorityQueue(int nelems=0) {
        if (nelems > 0)  init(nelems);
      }
      void init(int nelems) {
        a.resize(nelems);  // pre-allocate space
        a.resize(0);       // nothing stored yet (does not free memory)
      }
      void clear() {
        a.resize(0);
      }
      void insert(const T& t) {
        a.append(t);
        upheap();
      }
      void remove(T& t) {
        int last = a.len() - 1;
        if (last < 0) return;
        t = a[0];
        if (last > 0) a[0] = a[last];
        a.resize(last);  // remove last element from array
        downheap();
      }
    private:
      // bubble up last element to a correct position
      void upheap() {
        int n = a.len() - 1;
        while (n > 0) {
          int parent = (n-1) / 2;
          if (a[parent] <= a[n]) break;
          T tmp = a[parent];
          a[parent] = a[n];
          a[n] = tmp;
          n = parent;
        }
      }
      // trickle down first element to a correct position
      void downheap() {
        int n = 0;
        int len = a.len();
        int left = 2*n + 1;
        int right = left + 1;
        while (left < len) {
          if (right < len && a[right] <= a[left]) {
            if (a[n] <= a[right]) break;
            T tmp = a[right];
            a[right] = a[n];
            a[n] = tmp;
            n = right;
          }
          else {
            if (a[n] <= a[left]) break;
            T tmp = a[left];
            a[left] = a[n];
            a[n] = tmp;
            n = left;
          }
          left = 2*n + 1;
          right = left + 1;
        }
      }
      Array<T> a;
  };


  ///////////////////////////////////////////////////////////////////////////
  //
  // Grid is 3D lattice of grid points with user-definable index ranges.
  //
  ///////////////////////////////////////////////////////////////////////////

  // 3-integer vector, used for indexing from a 3D grid
  struct Ivec {
    int i, j, k;
    Ivec(int n=0) : i(n), j(n), k(n) { }
    Ivec(int ni, int nj, int nk) : i(ni), j(nj), k(nk) { }
    int operator==(const Ivec& n) { return(i==n.i && j==n.j && k==n.k); }
#ifdef MSM_MIGRATION
    void pup(PUP::er& p) {
      p|i, p|j, p|k;
    }
#endif
  };

  // index range for 3D lattice of grid points
  class IndexRange {
    public:
      IndexRange() : nlower(0), nextent(0) { }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        ASSERT(pni >= 0 && pnj >= 0 && pnk >= 0);
        nlower = Ivec(pia, pja, pka);
        nextent = Ivec(pni, pnj, pnk);
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        set(pia, pib-pia+1, pja, pjb-pja+1, pka, pkb-pka+1);
      }
      int ia() const { return nlower.i; }
      int ib() const { return nlower.i + nextent.i - 1; }
      int ja() const { return nlower.j; }
      int jb() const { return nlower.j + nextent.j - 1; }
      int ka() const { return nlower.k; }
      int kb() const { return nlower.k + nextent.k - 1; }
      int ni() const { return nextent.i; }
      int nj() const { return nextent.j; }
      int nk() const { return nextent.k; }
      int nn() const { return nextent.i * nextent.j * nextent.k; }
      Ivec lower() const { return nlower; }
      Ivec extent() const { return nextent; }
      int operator<=(const IndexRange& n) {
        // true if this IndexRange fits inside n
        return ( ia() >= n.ia() && ib() <= n.ib() &&
                 ja() >= n.ja() && jb() <= n.jb() &&
                 ka() >= n.ka() && kb() <= n.kb() );
      }
#ifdef MSM_MIGRATION
      void pup(PUP::er& p) {
        p|nlower, p|nextent;
      }
#endif
    protected:
      Ivec nlower;   // index for lowest corner of rectangular lattice
      Ivec nextent;  // extent of lattice along each dimension
  };

  // storage and indexing for 3D lattice of grid points
  // with fixed buffer storage no larger than size of block
  template <class T> class Grid;

  template <class T, int N>
  class GridFixed : public IndexRange {
    friend class Grid<T>;
    public:
      GridFixed() { }
      void init(const IndexRange& n) {
        nlower = n.lower();
        nextent = n.extent();
        ASSERT(nextent.i * nextent.j * nextent.k <= N);
      }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        IndexRange::set(pia, pni, pja, pnj, pka, pnk);
        ASSERT(nextent.i * nextent.j * nextent.k <= N);
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        IndexRange::setbounds(pia, pib, pja, pjb, pka, pkb);
        ASSERT(nextent.i * nextent.j * nextent.k <= N);
      }
      const T& operator()(int i, int j, int k) const {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      const T& operator()(const Ivec& n) const {
        return this->operator()(n.i, n.j, n.k);
      }
      const T& elem(int i, int j, int k) const {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      T& operator()(int i, int j, int k) {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      T& operator()(const Ivec& n) {
        return this->operator()(n.i, n.j, n.k);
      }
      T& elem(int i, int j, int k) {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      int flatindex(int i, int j, int k) const {
        return ((k-ka())*nj() + (j-ja()))*ni() + (i-ia());
      }
      const T *buffer() const { return gdata; }
      T *buffer() { return gdata; }

      // use to zero out grid
      void reset(const T& t) {
        int len = nn();
        for (int n = 0;  n < len;  n++) { gdata[n] = t; }
      }

      // use to modify the indexing by changing lower corner
      void updateLower(const Ivec& n) { nlower = n; }

      // accumulate another grid into this grid
      // the grid to be added must fit within this grid's index range
      GridFixed<T,N>& operator+=(const GridFixed<T,N>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        const T *gbuf = g.gdata.buffer();
        T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              buf[i + ijkoff] += gbuf[index];
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      void extract(GridFixed<T,N>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        T *gbuf = g.gdata.buffer();
        const T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              gbuf[index] = buf[i + ijkoff];
            }
          }
        }
      }

    private:
      T gdata[N];
  };

  // storage and indexing for 3D lattice of grid points
  template <class T>
  class Grid : public IndexRange {
    public:
      Grid() { }
      void init(const IndexRange& n) {
        nlower = n.lower();
        nextent = n.extent();
        gdata.resize(nn());
      }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        IndexRange::set(pia, pni, pja, pnj, pka, pnk);
        gdata.resize(nn());
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        IndexRange::setbounds(pia, pib, pja, pjb, pka, pkb);
        gdata.resize(nn());
      }
      void resize(int n) { // reserve space but don't set grid indexing
        gdata.resize(n);
      }
      const T& operator()(int i, int j, int k) const {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      const T& operator()(const Ivec& n) const {
        return this->operator()(n.i, n.j, n.k);
      }
      const T& elem(int i, int j, int k) const {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      T& operator()(int i, int j, int k) {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      T& operator()(const Ivec& n) {
        return this->operator()(n.i, n.j, n.k);
      }
      T& elem(int i, int j, int k) {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      int flatindex(int i, int j, int k) const {
        return ((k-ka())*nj() + (j-ja()))*ni() + (i-ia());
      }
      const Array<T>& data() const { return gdata; }
      Array<T>& data() { return gdata; }

      // use to zero out grid
      void reset(const T& t) {
        T *buf = gdata.buffer();
        int len = nn();
        for (int n = 0;  n < len;  n++) { buf[n] = t; }
      }

      // use to modify the indexing by changing lower corner
      void updateLower(const Ivec& n) { nlower = n; }

      // accumulate another grid into this grid
      // the grid to be added must fit within this grid's index range
      Grid<T>& operator+=(const Grid<T>& g) {
#if 1
        ASSERT(IndexRange(g) <= IndexRange(*this));
#else
        if ( ! (IndexRange(g) <= IndexRange(*this)) ) {
          Grid<T> tmp = *this;
          // expand myself to hold sum
          int ia = nlower.i;
          if (ia > g.nlower.i) ia = g.nlower.i;
          int ja = nlower.j;
          if (ja > g.nlower.j) ja = g.nlower.j;
          int ka = nlower.k;
          if (ka > g.nlower.k) ka = g.nlower.k;
          int ib1 = nlower.i + nextent.i;
          int gib1 = g.nlower.i + g.nextent.i;
          if (ib1 < gib1) ib1 = gib1;
          int jb1 = nlower.j + nextent.j;
          int gjb1 = g.nlower.j + g.nextent.j;
          if (jb1 < gjb1) jb1 = gjb1;
          int kb1 = nlower.k + nextent.k;
          int gkb1 = g.nlower.k + g.nextent.k;
          if (kb1 < gkb1) kb1 = gkb1;
          setbounds(ia, ib1-1, ja, jb1-1, ka, kb1-1);
          reset(0);  // make sure constructor for T accepts "0" as its zero
          // now copy "tmp" grid elements into my expanded self
          int index = 0;
          int ni = nextent.i;
          int nij = nextent.i * nextent.j;
          int koff = (tmp.nlower.k - nlower.k) * nij
            + (tmp.nlower.j - nlower.j) * ni + (tmp.nlower.i - nlower.i);
          const T *gbuf = tmp.gdata.buffer();
          T *buf = gdata.buffer();
          for (int k = 0;  k < tmp.nextent.k;  k++) {
            int jkoff = k * nij + koff;
            for (int j = 0;  j < tmp.nextent.j;  j++) {
              int ijkoff = j * ni + jkoff;
              for (int i = 0;  i < tmp.nextent.i;  i++, index++) {
                buf[i + ijkoff] = gbuf[index];
              }
            }
          }
        }
#endif
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        const T *gbuf = g.gdata.buffer();
        T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              buf[i + ijkoff] += gbuf[index];
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      void extract(Grid<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        T *gbuf = g.gdata.buffer();
        const T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              gbuf[index] = buf[i + ijkoff];
            }
          }
        }
      }

      // accumulate a fixed size grid into this grid
      // the grid to be added must fit within this grid's index range
      template <int N>
      Grid<T>& operator+=(const GridFixed<T,N>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        const T *gbuf = g.buffer();
        T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              buf[i + ijkoff] += gbuf[index];
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      template <int N>
      void extract(GridFixed<T,N>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        T *gbuf = g.buffer();
        const T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              gbuf[index] = buf[i + ijkoff];
            }
          }
        }
      }

    private:
      Array<T> gdata;
  }; // Grid


  ///////////////////////////////////////////////////////////////////////////
  //
  // Map object 
  //
  ///////////////////////////////////////////////////////////////////////////

  // index a block from the MSM grid hierarchy
  struct BlockIndex {
    int level;
    Ivec n;
    BlockIndex() : level(0), n(0) { }
    BlockIndex(int ll, const Ivec& nn) : level(ll), n(nn) { }
#ifdef MSM_MIGRATION
    void pup(PUP::er& p) {
      p|level, p|n;
    }
#endif
  };

  // for uppermost levels of hierarchy
  // fold out image charges along periodic boundaries
  // to fill up desired block size
  struct FoldFactor {
    int active;   // is some numrep dimension > 1?
    Ivec numrep;  // number of replications along each dimension
    FoldFactor() : active(0), numrep(1) { }
    FoldFactor(int i, int j, int k) { set(i,j,k); }
    void set(int i, int j, int k) {
      if (i <= 0) i = 1;
      if (j <= 0) j = 1;
      if (k <= 0) k = 1;
      if (i > 1 || j > 1 || k > 1) active = 1;
      numrep = Ivec(i, j, k);
    }
  };

  // sending part of an extended grid calculation to another block
  struct BlockSend {
    BlockIndex nblock;       // relative block index
    IndexRange nrange;       // relative grid index range
    BlockIndex nblock_wrap;  // true block index
    IndexRange nrange_wrap;  // true grid index range
    void reset() {
      nblock = BlockIndex();
      nrange = IndexRange();
      nblock_wrap = BlockIndex();
      nrange_wrap = IndexRange();
    } // reset
#ifdef MSM_MIGRATION
    void pup(PUP::er& p) {
      p|nblock, p|nrange, p|nblock_wrap, p|nrange_wrap;
    }
#endif
  };

  struct PatchSend {
    IndexRange nrange;         // true grid index range from my block
    IndexRange nrange_unwrap;  // relative grid index range for patch
    int patchID;               // patch ID
    void reset() {
      nrange = IndexRange();
      nrange_unwrap = IndexRange();
      patchID = -1;
    } // reset
  };

  // one PatchDiagram for each patch
  // maintain a Grid of PatchDiagram, indexed by patch ID
  struct PatchDiagram {
    IndexRange nrange;       // shows subset of MSM h-grid covering this patch
    Array<BlockSend> send;   // array of blocks to which this patch sends
    int numRecvs;            // number of blocks from which this patch receives
    void reset() {
      nrange = IndexRange();
      send.resize(0);
      numRecvs = 0;
    } // reset
  };

  // one BlockDiagram for each block of each level of each MSM grid
  // maintain a Grid of BlockDiagram for each level
  struct BlockDiagram {
    IndexRange nrange;            // subset of MSM grid for this block
    IndexRange nrangeCutoff;      // expanded subgrid for cutoff calculation
    IndexRange nrangeRestricted;  // (level+1) subgrid for restriction
    IndexRange nrangeProlongated; // (level-1) subgrid for prolongation
    Array<BlockSend> sendUp;      // send up charge to blocks on (level+1)
    Array<BlockSend> sendAcross;  // send across potential to blocks on (level)
    Array<int> indexGridCutoff;   // index of MsmGridCutoff chares to calculate
                                  // each charge -> potential block interaction
    Array<int> recvGridCutoff;    // index of MsmGridCutoff chares contributing
                                  // back into my potential block
    Array<BlockSend> sendDown;    // send down potential to blocks on (level-1)
    Array<PatchSend> sendPatch;   // send my (level=0) potential block to patch
    int numRecvsCharge;           // number of expected receives of charge
    int numRecvsPotential;        // number of expected receives of potential

    void reset() {
      nrange = IndexRange();
      nrangeCutoff = IndexRange();
      nrangeRestricted = IndexRange();
      nrangeProlongated = IndexRange();
      sendUp.resize(0);
      sendAcross.resize(0);
      sendDown.resize(0);
      sendPatch.resize(0);
      numRecvsCharge = 0;
      numRecvsPotential = 0;
    } // reset
  };


  struct Map {
    Array<IndexRange> gridrange;  // dimensions for each MSM grid level

    Array<Grid<Float> > gc;       // grid constant weights for each level
    Array<Grid<Float> > gvc;      // virial grid weights for each level
    Grid<Float> grespro;          // restriction / prolongation nonzero stencil
                                  // requires correct IndexOffset array

    Array<Grid<C1Matrix> > gc_c1hermite;    // grid constant weights C1 Hermite
    Array<Grid<C1Matrix> > gvc_c1hermite;   // virial grid weights C1 Hermite
    Array<Grid<C1Matrix> > gres_c1hermite;  // restriction stencil C1 Hermite
    Array<Grid<C1Matrix> > gpro_c1hermite;  // prolongation stencil C1 Hermite
                                            // requires index offsets

    Array<PatchDiagram> patchList;
    Array<Grid<BlockDiagram> > blockLevel;

    int ispx, ispy, ispz;         // is periodic in x, y, z?

    Array<int> bsx, bsy, bsz;     // block size in x, y, z for each level

    Array<FoldFactor> foldfactor; // for uppermost grid levels
      // replicate periodic dimensions in order to fill up block size

    // clip index to grid level, using periodicity flags
    Ivec clipIndexToLevel(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      Ivec pn(n);
      if ( ! ispx) {
        int a = gridrange[level].ia();
        int b = gridrange[level].ib();
        if (pn.i < a) pn.i = a;
        if (pn.i > b) pn.i = b;
      }
      if ( ! ispy) {
        int a = gridrange[level].ja();
        int b = gridrange[level].jb();
        if (pn.j < a) pn.j = a;
        if (pn.j > b) pn.j = b;
      }
      if ( ! ispz) {
        int a = gridrange[level].ka();
        int b = gridrange[level].kb();
        if (pn.k < a) pn.k = a;
        if (pn.k > b) pn.k = b;
      }
      return pn;
    }

    // determine relative (unwrapped) block index for the given grid index
    BlockIndex blockOfGridIndex(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      BlockIndex bn;
      // we want floor((i - ia) / bsx), etc.
      // modify case i < ia to avoid integer division of negative numbers
      int d = n.i - gridrange[level].ia();
      bn.n.i = (d >= 0 ? d / bsx[level] : -((-d+bsx[level]-1) / bsx[level]));
      d = n.j - gridrange[level].ja();
      bn.n.j = (d >= 0 ? d / bsy[level] : -((-d+bsy[level]-1) / bsy[level]));
      d = n.k - gridrange[level].ka();
      bn.n.k = (d >= 0 ? d / bsz[level] : -((-d+bsz[level]-1) / bsz[level]));
      bn.level = level;
      return bn;
    }

    // determine relative (unwrapped) block index for the given grid index
    // for unfolded replication of image charges
    BlockIndex blockOfGridIndexFold(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      BlockIndex bn;
      int bsi = foldfactor[level].numrep.i * bsx[level];
      int bsj = foldfactor[level].numrep.j * bsy[level];
      int bsk = foldfactor[level].numrep.k * bsz[level];
      // we want floor((i - ia) / bsx), etc.
      // modify case i < ia to avoid integer division of negative numbers
      int d = n.i - gridrange[level].ia();
      bn.n.i = (d >= 0 ? d / bsi : -((-d+bsi-1) / bsi));
      d = n.j - gridrange[level].ja();
      bn.n.j = (d >= 0 ? d / bsj : -((-d+bsj-1) / bsj));
      d = n.k - gridrange[level].ka();
      bn.n.k = (d >= 0 ? d / bsk : -((-d+bsk-1) / bsk));
      bn.level = level;
      return bn;
    }

    // find the natural index range of the given relative block number
    IndexRange indexRangeOfBlock(const BlockIndex& nb) const {
      ASSERT(nb.level >= 0 && nb.level < gridrange.len());
      IndexRange nr;
      int ia = nb.n.i * bsx[nb.level] + gridrange[nb.level].ia();
      int ja = nb.n.j * bsy[nb.level] + gridrange[nb.level].ja();
      int ka = nb.n.k * bsz[nb.level] + gridrange[nb.level].ka();
      nr.set(ia, bsx[nb.level], ja, bsy[nb.level], ka, bsz[nb.level]);
      return nr;
    }

    // find the natural index range of the given relative block number
    // for unfolded replication of image charges
    IndexRange indexRangeOfBlockFold(const BlockIndex& nb) const {
      ASSERT(nb.level >= 0 && nb.level < gridrange.len());
      int bsi = foldfactor[nb.level].numrep.i * bsx[nb.level];
      int bsj = foldfactor[nb.level].numrep.j * bsy[nb.level];
      int bsk = foldfactor[nb.level].numrep.k * bsz[nb.level];
      IndexRange nr;
      int ia = nb.n.i * bsi + gridrange[nb.level].ia();
      int ja = nb.n.j * bsj + gridrange[nb.level].ja();
      int ka = nb.n.k * bsk + gridrange[nb.level].ka();
      nr.set(ia, bsi, ja, bsj, ka, bsk);
      return nr;
    }

    // clip the natural block index range to not exceed the given index range
    IndexRange clipBlockToIndexRange(const BlockIndex& nb,
        const IndexRange& nrange) const {
      IndexRange nr = indexRangeOfBlock(nb);
      int nia = nrange.ia();
      int nib = nrange.ib();
      int nja = nrange.ja();
      int njb = nrange.jb();
      int nka = nrange.ka();
      int nkb = nrange.kb();
      int ia = nr.ia();
      if (ia < nia) ia = nia;
      int ib = nr.ib();
      if (ib > nib) ib = nib;
      int ja = nr.ja();
      if (ja < nja) ja = nja;
      int jb = nr.jb();
      if (jb > njb) jb = njb;
      int ka = nr.ka();
      if (ka < nka) ka = nka;
      int kb = nr.kb();
      if (kb > nkb) kb = nkb;
      nr.setbounds(ia, ib, ja, jb, ka, kb);
      return nr;
    }

    // clip the natural block index range to not exceed the given index range
    // for unfolded replication of image charges
    IndexRange clipBlockToIndexRangeFold(const BlockIndex& nb,
        const IndexRange& nrange) const {
      IndexRange nr = indexRangeOfBlockFold(nb);
      int nia = nrange.ia();
      int nib = nrange.ib();
      int nja = nrange.ja();
      int njb = nrange.jb();
      int nka = nrange.ka();
      int nkb = nrange.kb();
      int ia = nr.ia();
      if (ia < nia) ia = nia;
      int ib = nr.ib();
      if (ib > nib) ib = nib;
      int ja = nr.ja();
      if (ja < nja) ja = nja;
      int jb = nr.jb();
      if (jb > njb) jb = njb;
      int ka = nr.ka();
      if (ka < nka) ka = nka;
      int kb = nr.kb();
      if (kb > nkb) kb = nkb;
      nr.setbounds(ia, ib, ja, jb, ka, kb);
      return nr;
    }

    // set the nblock_wrap and nrange_wrap fields based on periodicity
    void wrapBlockSend(BlockSend& bs) const {
      BlockIndex nb = bs.nblock;
      IndexRange nr = bs.nrange;
      int level = bs.nblock.level;
      ASSERT(level >= 0 && level < blockLevel.len());
      int ni = blockLevel[level].ni();
      int nj = blockLevel[level].nj();
      int nk = blockLevel[level].nk();
      int di=0, dj=0, dk=0;
      if (ispx) {
        while (nb.n.i < 0) {
          nb.n.i += ni;
          di += ni * bsx[level];
        }
        while (nb.n.i >= ni) {
          nb.n.i -= ni;
          di -= ni * bsx[level];
        }
      }
      if (ispy) {
        while (nb.n.j < 0) {
          nb.n.j += nj;
          dj += nj * bsy[level];
        }
        while (nb.n.j >= nj) {
          nb.n.j -= nj;
          dj -= nj * bsy[level];
        }
      }
      if (ispz) {
        while (nb.n.k < 0) {
          nb.n.k += nk;
          dk += nk * bsz[level];
        }
        while (nb.n.k >= nk) {
          nb.n.k -= nk;
          dk -= nk * bsz[level];
        }
      }
      int ia = nr.ia();
      int ib = nr.ib();
      int ja = nr.ja();
      int jb = nr.jb();
      int ka = nr.ka();
      int kb = nr.kb();
      nr.setbounds(ia + di, ib + di, ja + dj, jb + dj, ka + dk, kb + dk);
      bs.nblock_wrap = nb;
      bs.nrange_wrap = nr;
    }

    // set the nblock_wrap and nrange_wrap fields based on periodicity
    // for unfolded replication of image charges
    void wrapBlockSendFold(BlockSend& bs) const {
      BlockIndex nb = bs.nblock;
      IndexRange nr = bs.nrange;
      int level = bs.nblock.level;
      ASSERT(level >= 0 && level < blockLevel.len());
      int foldi = foldfactor[level].numrep.i;
      int foldj = foldfactor[level].numrep.j;
      int foldk = foldfactor[level].numrep.k;
      int ni = blockLevel[level].ni();
      int nj = blockLevel[level].nj();
      int nk = blockLevel[level].nk();
      int bsi = foldi * bsx[level];
      int bsj = foldj * bsy[level];
      int bsk = foldk * bsz[level];
      int di=0, dj=0, dk=0;
      if (ispx) {
        while (nb.n.i < 0) {
          nb.n.i += ni;
          di += ni * bsi;
        }
        while (nb.n.i >= ni) {
          nb.n.i -= ni;
          di -= ni * bsi;
        }
      }
      if (ispy) {
        while (nb.n.j < 0) {
          nb.n.j += nj;
          dj += nj * bsj;
        }
        while (nb.n.j >= nj) {
          nb.n.j -= nj;
          dj -= nj * bsj;
        }
      }
      if (ispz) {
        while (nb.n.k < 0) {
          nb.n.k += nk;
          dk += nk * bsk;
        }
        while (nb.n.k >= nk) {
          nb.n.k -= nk;
          dk -= nk * bsk;
        }
      }
      int ia = nr.ia();
      int ib = nr.ib();
      int ja = nr.ja();
      int jb = nr.jb();
      int ka = nr.ka();
      int kb = nr.kb();
      nr.setbounds(ia + di, ib + di, ja + dj, jb + dj, ka + dk, kb + dk);
      bs.nblock_wrap = nb;
      bs.nrange_wrap = nr;
    }

    void wrapBlockIndex(BlockIndex& bn) const {
      int level = bn.level;
      int ni = blockLevel[level].ni();
      int nj = blockLevel[level].nj();
      int nk = blockLevel[level].nk();
      if (ispx) {
        while (bn.n.i < 0) {
          bn.n.i += ni;
        }
        while (bn.n.i >= ni) {
          bn.n.i -= ni;
        }
      }
      if (ispy) {
        while (bn.n.j < 0) {
          bn.n.j += nj;
        }
        while (bn.n.j >= nj) {
          bn.n.j -= nj;
        }
      }
      if (ispz) {
        while (bn.n.k < 0) {
          bn.n.k += nk;
        }
        while (bn.n.k >= nk) {
          bn.n.k -= nk;
        }
      }
    }

  }; // Map


  struct AtomCoord {
    Position position;
    Real charge;
    int id;
  };

  typedef Array<AtomCoord> AtomCoordArray;
  typedef Array<Force> ForceArray;

  struct PatchData;
  typedef Array<PatchData *> PatchPtrArray;

} // namespace msm

#endif // MSMMAP_H
