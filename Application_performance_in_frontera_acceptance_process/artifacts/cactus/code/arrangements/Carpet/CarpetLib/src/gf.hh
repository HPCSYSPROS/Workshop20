#ifndef GF_HH
#define GF_HH

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "bbox.hh"
#include "bboxset.hh"
#include "data.hh"
#include "defs.hh"
#include "dh.hh"
#include "ggf.hh"
#include "th.hh"
#include "vect.hh"

using namespace std;

// A real grid function
template <typename T> class gf : public ggf {

public:
  // Constructors
  gf(const int varindex, const operator_type transport_operator, th &t, dh &d,
     const int prolongation_order_time, const int vectorlength,
     const int vectorindex, ggf *const vectorleader);

  // Destructors
  virtual ~gf();

  // Helpers

  virtual gdata *typed_data(int tl, int rl, int lc, int ml) const {
    data<T> *const vl =
        this->vectorleader
            ? ((gf<T> *)this->vectorleader)->typed_data_pointer(tl, rl, lc, ml)
            : NULL;
    return new data<T>(this->varindex, h.refcent, this->transport_operator,
                       this->vectorlength, this->vectorindex, vl);
  }

  virtual gdata *new_typed_data() const {
    return new data<T>(this->varindex, h.refcent, this->transport_operator, 1,
                       0, NULL);
  }

  // Access to the data

  data<T> const *typed_data_pointer(int tl, int rl, int lc, int ml) const {
    assert(rl >= 0 and rl < h.reflevels());
    assert(lc >= 0 and lc < h.local_components(rl));
    assert(ml >= 0 and ml < h.mglevels());
    assert(tl >= 0 and tl < timelevels(ml, rl));
    return (data<T> const *)storage.AT(ml).AT(rl).AT(lc).AT(tl);
  }
  data<T> *typed_data_pointer(int tl, int rl, int lc, int ml) {
    assert(rl >= 0 and rl < h.reflevels());
    assert(lc >= 0 and lc < h.local_components(rl));
    assert(ml >= 0 and ml < h.mglevels());
    assert(tl >= 0 and tl < timelevels(ml, rl));
    return (data<T> *)storage.AT(ml).AT(rl).AT(lc).AT(tl);
  }

  // Output
  virtual size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
  virtual ostream &output(ostream &os) const;

private:
  gf();                      // canonical default construtor
  gf(const gf &);            // canonical copy construtor
  gf &operator=(const gf &); // canonical copy
};

#endif // GF_HH
