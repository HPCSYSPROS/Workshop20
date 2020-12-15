#ifndef BBOXSET1_HH
#define BBOXSET1_HH

#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"

using namespace std;

// Choose the implementation of bboxset by #defining exactly one of
// these
#undef BBOXSET_SET     // outdated
#undef BBOXSET_LIST    // well tested
#define BBOXSET_VECTOR // brand new

namespace bboxset1 {

// Forward declaration
template <typename T, int D> class bboxset;

// template<typename T,int D>
// bboxset<T,D> operator+ (const bbox<T,D>& b1, const bbox<T,D>& b2);
// template<typename T,int D>
// bboxset<T,D> operator+ (const bbox<T,D>& b, const bboxset<T,D>& s);

// template<typename T,int D>
// bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2);
// template<typename T,int D>
// bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s);

// Input
template <typename T, int D> istream &operator>>(istream &is, bboxset<T, D> &s);

// Output
template <typename T, int D>
ostream &operator<<(ostream &os, const bboxset<T, D> &s);

// Bounding box set class
template <typename T, int D> class bboxset {

  // Cost annotations depend on the number of bset elements n, and
  // assume that normalization is skipped.

  struct skip_normalize_t {
    bboxset<T, D> &s;
    bool const saved_skip_normalize;
    skip_normalize_t(bboxset<T, D> &s_)
        : s(s_), saved_skip_normalize(s.skip_normalize) {
      s.skip_normalize = true;
    }
    ~skip_normalize_t() {
      s.skip_normalize = saved_skip_normalize;
      s.normalize();
    }
  };

  // Types
  typedef bbox<T, D> box;
#ifdef BBOXSET_SET
// S typedef set<box> bset;
#endif
#ifdef BBOXSET_LIST
  typedef list<box> bset;
#endif
#ifdef BBOXSET_VECTOR
  typedef vector<box> bset;
#endif

  // Fields
  bset bs;
  // Invariant:
  // All bboxes have the same stride.
  // No bbox is empty.
  // The bboxes don't overlap.

  bool skip_normalize;

public:
  // Constructors
  bboxset();                 // cost: O(1)
  bboxset(const box &b);     // cost: O(1)
  bboxset(const bboxset &s); // cost: O(n)

  bboxset(const list<box> &lb);
  bboxset(const set<box> &sb);
  bboxset(const vector<box> &vb);
  bboxset(const vector<list<box> > &vlb);
  template <typename U>
  bboxset(const vector<U> &vb, const bbox<T, D> U::*const v);
  template <typename U> bboxset(const vector<U> &vb, const bboxset U::*const v);

  static bboxset poison();

  // Invariant
  bool invariant() const CCTK_MEMBER_ATTRIBUTE_PURE;

private:
  // Normalisation
  void normalize();

public:
  // Accessors
  bool empty() const { return bs.empty(); } // cost: O(1)
  // T size () const;
  size_type size() const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  int setsize() const { return bs.size(); }          // cost: O(1)

  // Find out whether the bbox contains the point x
  bool
  contains(const vect<T, D> &x) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)

  // Find out whether this bboxset intersects the bbox b
  bool intersects(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)

  // Add (bboxes that don't overlap)
  bboxset &operator+=(const box &b);     // cost: O(1)
  bboxset &operator+=(const bboxset &s); // cost: O(n)
  bboxset &add_transfer(bboxset &s);     // cost: O(1)
  bboxset
  operator+(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset
  operator+(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)

  // Union
  bboxset &operator|=(const box &b);     // cost: O(n)
  bboxset &operator|=(const bboxset &s); // cost: O(n^2)
  bboxset
  operator|(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset
  operator|(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n^2)

  // Intersection
  bboxset
  operator&(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset
  operator&(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  bboxset &operator&=(const box &b);                            // cost: O(n)
  bboxset &operator&=(const bboxset &s);                        // cost: O(n)

  // Difference
  // friend bboxset operator- <T,D>(const box& b1, const box& b2);
  static bboxset minus(const box &b1,
                       const box &b2) CCTK_ATTRIBUTE_PURE; // cost: O(1)
  bboxset
  operator-(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  // cost: O(n)
  bboxset &operator-=(const box &b) {
    *this = *this - b;
    assert(invariant());
    return *this;
  }
  bboxset &operator-=(const bboxset &s); // cost: O(n^2)
  bboxset
  operator-(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n^2)
  // friend bboxset operator- <T,D>(const box& b, const bboxset& s);
  static bboxset minus(const box &b,
                       const bboxset &s) CCTK_ATTRIBUTE_PURE; // cost: O(n^2)

  /** Find a bbox containing the whole set.  */
  box container() const CCTK_MEMBER_ATTRIBUTE_PURE; // cost: O(n)
  /** Find the pseudo-inverse.  */
  bboxset pseudo_inverse(const int n) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Expand (enlarge) the bboxset by multiples of the stride.  */
  bboxset expand(const vect<T, D> &lo,
                 const vect<T, D> &hi) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bboxset expand(const vect<vect<T, D>, 2> &lohi) const {
    return expand(lohi[0], lohi[1]);
  }

  /** Shift the bboxset by multiples of the stride.  */
  bboxset shift(const vect<T, D> &v) const { return expand(-v, v); }

  /** Expand (enlarge) the bboxset by multiples of a fraction of the
      stride.  */
  // cost: O(n^2) in general, but only O(n) for shifting
  bboxset expand(const vect<T, D> &lo, const vect<T, D> &hi,
                 const vect<T, D> &denom) const CCTK_MEMBER_ATTRIBUTE_PURE;
  // cost: O(n^2) in general, but only O(n) for shifting
  bboxset expand(const vect<vect<T, D>, 2> &lohi,
                 const vect<T, D> &denom) const {
    return expand(lohi[0], lohi[1], denom);
  }

  /** Shift the bboxset by multiples of a fraction of the stride.  */
  // cost: O(n)
  bboxset shift(const vect<T, D> &v, const vect<T, D> &denom) const {
    return expand(-v, v, denom);
  }

  /** Find the smallest b-compatible box around this bbox.
      ("compatible" means having the same stride.)  */
  bboxset expanded_for(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

// TODO: this is incorrect
#if 1
  /** Find the largest b-compatible box inside this bbox.  */
  bboxset contracted_for(const box &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
#endif

  // Equality
  bool operator==(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator!=(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator<(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator<=(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator>(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator>=(const bboxset &s) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Iterators
  typedef typename bset::const_iterator const_iterator;
  typedef typename bset::iterator iterator;

  const_iterator begin() const { return bs.begin(); }
  const_iterator end() const { return bs.end(); }
  //   iterator begin () const { return bs.begin(); }
  //   iterator end () const   { return bs.end(); }

  template <typename C> void serialise(C &out) const;

  // Memory usage
  size_t memory() const { return memoryof(bs); }

  // Input
  istream &input(istream &is);

  // Output
  ostream &output(ostream &os) const;
};

template <typename T, int D>
inline bboxset<T, D> operator+(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  return bboxset<T, D>(b1) + bboxset<T, D>(b2);
}

template <typename T, int D>
inline bboxset<T, D> operator+(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) + s;
}

// cost: O(1)
template <typename T, int D>
inline bboxset<T, D> operator-(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  return bboxset<T, D>::minus(b1, b2);
}

// cost: O(n^2)
template <typename T, int D>
inline bboxset<T, D> operator-(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>::minus(b, s);
}

template <typename T, int D>
inline bboxset<T, D> operator|(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return s | b;
}

template <typename T, int D>
inline bboxset<T, D> operator&(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return s & b;
}

template <typename T, int D>
inline bool operator==(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) == s;
}

template <typename T, int D>
inline bool operator!=(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) != s;
}

template <typename T, int D>
inline bool operator<(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) < s;
}

template <typename T, int D>
inline bool operator<=(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) <= s;
}

template <typename T, int D>
inline bool operator>(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) > s;
}

template <typename T, int D>
inline bool operator>=(const bbox<T, D> &b, const bboxset<T, D> &s) {
  return bboxset<T, D>(b) >= s;
}

template <typename T, int D>
inline bool operator==(const bboxset<T, D> &s, const bbox<T, D> &b) {
  return s == bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator!=(const bboxset<T, D> &s, const bbox<T, D> &b) {
  return s != bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator<(const bboxset<T, D> &s, const bbox<T, D> &b) {
  return s < bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator<=(const bboxset<T, D> &s, const bbox<T, D> &b) {
  return s <= bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator>(const bboxset<T, D> &s, const bbox<T, D> &b) {
  return s > bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator>=(const bboxset<T, D> &s, const bbox<T, D> &b) {
  return s >= bboxset<T, D>(b);
}

// Memory usage
template <typename T, int D> inline size_t memoryof(bboxset<T, D> const &s) {
  return s.memory();
}

// Input
template <typename T, int D>
inline istream &operator>>(istream &is, bboxset<T, D> &s) {
  return s.input(is);
}

// Output
template <typename T, int D>
inline ostream &operator<<(ostream &os, const bboxset<T, D> &s) {
  return s.output(os);
}

////////////////////////////////////////////////////////////////////////////////
//  Implementation
////////////////////////////////////////////////////////////////////////////////

#define SKIP_NORMALIZE(s) skip_normalize_t skip_normalize_v(s)
// #define SKIP_NORMALIZE(s) do { } while(0)

// Constructors
template <typename T, int D> bboxset<T, D>::bboxset() : skip_normalize(false) {
  assert(invariant());
}

template <typename T, int D>
bboxset<T, D>::bboxset(const box &b)
    : skip_normalize(false) {
  // S if (not b.empty()) bs.insert(b);
  if (not b.empty())
    bs.push_back(b);
  assert(invariant());
}

template <typename T, int D>
bboxset<T, D>::bboxset(const bboxset &s)
    : bs(s.bs), skip_normalize(false) {
  assert(invariant());
}

template <typename T, int D>
bboxset<T, D>::bboxset(const list<box> &lb)
    : skip_normalize(false) {
  SKIP_NORMALIZE(*this);
  for (typename list<box>::const_iterator li = lb.begin(), le = lb.end();
       li != le; ++li) {
    *this |= *li;
  }
}

template <typename T, int D>
bboxset<T, D>::bboxset(const set<box> &sb)
    : skip_normalize(false) {
  SKIP_NORMALIZE(*this);
  for (typename set<box>::const_iterator vi = sb.begin(), ve = sb.end();
       vi != ve; ++vi) {
    *this |= *vi;
  }
}

template <typename T, int D>
bboxset<T, D>::bboxset(const vector<box> &vb)
    : skip_normalize(false) {
  SKIP_NORMALIZE(*this);
  for (typename vector<box>::const_iterator vi = vb.begin(), ve = vb.end();
       vi != ve; ++vi) {
    *this |= *vi;
  }
}

template <typename T, int D>
bboxset<T, D>::bboxset(const vector<list<box> > &vlb)
    : skip_normalize(false) {
  SKIP_NORMALIZE(*this);
  for (typename vector<list<box> >::const_iterator vli = vlb.begin(),
                                                   vle = vlb.end();
       vli != vle; ++vli) {
    *this |= bboxset(*vli);
  }
}

template <typename T, int D>
template <typename U>
bboxset<T, D>::bboxset(const vector<U> &vb, const bbox<T, D> U::*const v)
    : skip_normalize(false) {
  SKIP_NORMALIZE(*this);
  for (typename vector<U>::const_iterator vi = vb.begin(), ve = vb.end();
       vi != ve; ++vi) {
    *this |= (*vi).*v;
  }
}

template <typename T, int D>
template <typename U>
bboxset<T, D>::bboxset(const vector<U> &vb, const bboxset U::*const v)
    : skip_normalize(false) {
  SKIP_NORMALIZE(*this);
  for (typename vector<U>::const_iterator vi = vb.begin(), ve = vb.end();
       vi != ve; ++vi) {
    *this |= (*vi).*v;
  }
}

template <typename T, int D> bboxset<T, D> bboxset<T, D>::poison() {
  return bboxset(bbox<T, D>::poison());
}

// Invariant
template <typename T, int D> bool bboxset<T, D>::invariant() const {
// This is very slow when there are many bboxes
#if 0 && defined(CARPET_DEBUG)
  for (const_iterator bi=begin(), be=end(); bi!=be; ++bi) {
    if ((*bi).empty()) return false;
    if (not (*bi).is_aligned_with(*bs.begin())) return false;
    // check for overlap (quadratic -- expensive)
    for (const_iterator bi2=begin(); bi2!=bi; ++bi2) {
      if ((*bi2).intersects(*bi)) return false;
    }
  }
#endif
  return true;
}

// Normalisation
template <typename T, int D> void bboxset<T, D>::normalize() {
  if (skip_normalize)
    return;
  assert(invariant());

  bboxset const oldbs = *this;
  size_type const oldsize = this->size();

  // Split all bboxes into small pieces which have all their
  // boundaries aligned.
  for (int d = 0; d < D; ++d) {
    // Find all boundaries
    // S typedef set<T> buf;
    typedef vector<T> buf;
    buf sbnds;
    sbnds.reserve(2 * bs.size());
    for (typename bset::const_iterator si = bs.begin(), se = bs.end(); si != se;
         ++si) {
      box const &b = *si;
      int const bstr = b.stride()[d];
      int const blo = b.lower()[d];
      int const bhi = b.upper()[d] + bstr;
      // S sbnds.insert (blo);
      // S sbnds.insert (bhi);
      sbnds.push_back(blo);
      sbnds.push_back(bhi);
    }
    sort(sbnds.begin(), sbnds.end());
    typename buf::iterator const last = unique(sbnds.begin(), sbnds.end());
    sbnds.resize(last - sbnds.begin());
    // Split bboxes
    bset nbs;
    for (typename bset::const_iterator si = bs.begin(), se = bs.end(); si != se;
         ++si) {
      box const &b = *si;
      int const bstr = b.stride()[d];
      int const blo = b.lower()[d];
      int const bhi = b.upper()[d] + bstr;
      typename buf::const_iterator const ilo =
          lower_bound(sbnds.begin(), sbnds.end(), blo);
      typename buf::const_iterator const ihi =
          lower_bound(sbnds.begin(), sbnds.end(), bhi);
      assert(ilo != sbnds.end());
      assert(ihi != sbnds.end());
      assert(*ilo == blo);
      assert(*ihi == bhi);
      // Split one bbox
      for (typename buf::const_iterator curr = ilo; curr != ihi; ++curr) {
        typename buf::const_iterator next = curr;
        advance(next, 1);
        int const nblo = *curr;
        int const nbhi = *next;
        assert(nbhi > nblo); // ensure that the set remains sorted
        box const nb(b.lower().replace(d, nblo),
                     b.upper().replace(d, nbhi - bstr), b.stride());
        // S nbs.insert (nb);
        nbs.push_back(nb);
      }
    }
    // Replace old set
    bs.clear();
    bs.swap(nbs);
    assert(invariant());
  }

  // Combine bboxes if possible
  for (int d = 0; d < D; ++d) {
    bset nbs;
    while (not bs.empty()) {
      typename bset::iterator si = bs.begin();
      assert(si != bs.end());

      box const b = *si;
      int const bstr = b.stride()[d];
      int const blo = b.lower()[d];
      int const bhi = b.upper()[d] + bstr;

      for (typename bset::iterator nsi = nbs.begin(), nse = nbs.end();
           nsi != nse; ++nsi) {
        box const nb = *nsi;
        int const nblo = nb.lower()[d];
        int const nbhi = nb.upper()[d] + bstr;

        box const mb(nb.lower().replace(d, blo),
                     nb.upper().replace(d, bhi - bstr), nb.stride());

        // Check whether the other dimensions match
        if (b == mb) {
          // Check whether the bboxes are adjacent in this dimension
          if (nbhi == blo) {
            // Combine boxes, nb < b
            box const cb(b.lower().replace(d, nblo), b.upper(), b.stride());
            bs.erase(si);
            nbs.erase(nsi);
            // S bs.insert (cb);
            bs.push_back(cb);
            goto done;
          } else if (bhi == nblo) {
            // Combine boxes, b < nb
            box const cb(b.lower(), b.upper().replace(d, nbhi - bstr),
                         b.stride());
            bs.erase(si);
            nbs.erase(nsi);
            // S bs.insert (cb);
            bs.push_back(cb);
            goto done;
          }
        }
      }
      bs.erase(si);
      // S nbs.insert (b);
      nbs.push_back(b);
    done:;
    }
    bs.swap(nbs);
    assert(invariant());
  }

  // Sort bboxes lexicographically
  sort(bs.begin(), bs.end(), less<bbox<T, D> >());

  size_type const newsize = this->size();

  // Can't use operators on *this since these would call normalize again
  // assert (*this == oldbs);
  assert(newsize == oldsize);
}

// Accessors
// cost: O(n)
template <typename T, int D> size_type bboxset<T, D>::size() const {
  size_type s = 0;
  for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
    const size_type bsz = (*bi).size();
    assert(numeric_limits<size_type>::max() - bsz >= s);
    s += bsz;
  }
  return s;
}

// Queries

// Containment
// cost: O(n)
template <typename T, int D>
bool bboxset<T, D>::contains(const vect<T, D> &x) const {
  for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
    if ((*bi).contains(x))
      return true;
  }
  return false;
}

// Intersection
// cost: O(n)
template <typename T, int D>
bool bboxset<T, D>::intersects(const box &b) const {
  for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
    if ((*bi).intersects(b))
      return true;
  }
  return false;
}

// Add (bboxes that don't overlap)
// cost: O(1)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator+=(const box &b) {
  if (b.empty())
    return *this;
// This is very slow when there are many bboxes
#if 0 && defined(CARPET_DEBUG)
  // check for overlap
  for (const_iterator bi=begin(), be=end(); bi!=be; ++bi) {
    assert (not (*bi).intersects(b));
  }
#endif
  // S bs.insert(b);
  bs.push_back(b);
  assert(invariant());
  normalize();
  return *this;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator+=(const bboxset &s) {
  SKIP_NORMALIZE(*this);
  for (const_iterator bi = s.begin(), be = s.end(); bi != be; ++bi) {
    *this += *bi;
  }
  assert(invariant());
  return *this;
}

// L cost: O(1)
// cost O(n)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::add_transfer(bboxset &s) {
#ifdef BBOXSET_LIST
  bs.splice(bs.end(), s.bs);
#else
  bs.insert(bs.end(), s.bs.begin(), s.bs.end());
  s.bs.clear();
#endif
  assert(invariant());
  normalize();
  return *this;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator+(const box &b) const {
  bboxset r(*this);
  r += b;
  assert(r.invariant());
  return r;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator+(const bboxset &s) const {
  bboxset r(*this);
  r += s;
  assert(r.invariant());
  return r;
}

// Union
// cost: O(n)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator|=(const box &b) {
  if (b.empty())
    return *this;
#if 0
  // this has a cost of O(n^2)
  bboxset tmp = b - *this;
  add_transfer (tmp);
#else
  // this has a cost of O(n)
  bset oldbs;
  oldbs.swap(bs);
  bs.push_back(b);
  SKIP_NORMALIZE(*this);
  for (const_iterator bi = oldbs.begin(), be = oldbs.end(); bi != be; ++bi) {
    bboxset tmp = *bi - b;
    add_transfer(tmp);
  }
#endif
  assert(invariant());
  return *this;
}

// cost: O(n^2)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator|=(const bboxset &s) {
  bboxset tmp = s - *this;
  add_transfer(tmp);
  assert(invariant());
  return *this;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator|(const box &b) const {
  bboxset r(*this);
  r |= b;
  assert(r.invariant());
  return r;
}

// cost: O(n^2)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator|(const bboxset &s) const {
  bboxset r(*this);
  r |= s;
  assert(r.invariant());
  return r;
}

// Intersection
// cost: O(n)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator&(const box &b) const {
  // start with an empty set
  bboxset r;
  {
    SKIP_NORMALIZE(r);
    // walk all my elements
    for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
      // insert the intersection with the bbox
      r += *bi & b;
    }
  }
  assert(r.invariant());
  return r;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator&(const bboxset &s) const {
  // start with an empty set
  bboxset r;
  {
    SKIP_NORMALIZE(r);
    // walk all the bboxes
    for (const_iterator bi = s.begin(), be = s.end(); bi != be; ++bi) {
      // insert the intersection with this bbox
      bboxset tmp = *this & *bi;
      r.add_transfer(tmp);
    }
  }
  assert(r.invariant());
  return r;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator&=(const box &b) {
  *this = *this & b;
  assert(invariant());
  return *this;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator&=(const bboxset &s) {
  *this = *this & s;
  assert(invariant());
  return *this;
}

// Difference
// cost: O(1)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::minus(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  assert(b1.is_aligned_with(b2));
  if (b1.empty())
    return bboxset<T, D>();
  if (not b1.intersects(b2))
    return bboxset<T, D>(b1);
  vect<T, D> const b2lo = max(b1.lower(), b2.lower());
  vect<T, D> const b2up = min(b1.upper(), b2.upper());
  vect<T, D> const &b1lo = b1.lower();
  vect<T, D> const &b1up = b1.upper();
  vect<T, D> const &str = b1.stride();
  bboxset<T, D> r;
  SKIP_NORMALIZE(r);
  {
    for (int d = 0; d < D; ++d) {
      // make resulting bboxes as large as possible in x-direction
      // (for better consumption by Fortranly ordered arrays)
      vect<T, D> lb, ub;
      bbox<T, D> b;
      for (int dd = 0; dd < D; ++dd) {
        if (dd < d) {
          lb[dd] = b2lo[dd];
          ub[dd] = b2up[dd];
        } else if (dd > d) {
          lb[dd] = b1lo[dd];
          ub[dd] = b1up[dd];
        }
      }
      lb[d] = b1lo[d];
      ub[d] = b2lo[d] - str[d];
      b = bbox<T, D>(lb, ub, str);
      r += b;
      lb[d] = b2up[d] + str[d];
      ub[d] = b1up[d];
      b = bbox<T, D>(lb, ub, str);
      r += b;
    }
  }
  assert(r.invariant());
  return r;
}

// cost: O(n)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator-(const box &b) const {
  // start with an empty set
  bboxset r;
  {
    SKIP_NORMALIZE(r);
    // walk all my elements
    for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
      // insert the difference with the bbox
      bboxset tmp = *bi - b;
      r.add_transfer(tmp);
    }
  }
  assert(r.invariant());
  return r;
}

// cost: O(n^2)
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator-=(const bboxset &s) {
  SKIP_NORMALIZE(*this);
  for (const_iterator bi = s.begin(), be = s.end(); bi != be; ++bi) {
    *this -= *bi;
  }
  assert(invariant());
  return *this;
}

// cost: O(n^2)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator-(const bboxset &s) const {
  bboxset r(*this);
  r -= s;
  assert(r.invariant());
  return r;
}

// cost: O(n^2)
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::minus(const bbox<T, D> &b,
                                   const bboxset<T, D> &s) {
  bboxset<T, D> r = bboxset<T, D>(b) - s;
  assert(r.invariant());
  return r;
}

// cost: O(n)
template <typename T, int D>
typename bboxset<T, D>::box bboxset<T, D>::container() const {
  box b;
  for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
    b = b.expanded_containing(*bi);
  }
  return b;
}

template <typename T, int D>
bboxset<T, D> bboxset<T, D>::pseudo_inverse(const int n) const {
  assert(not empty());
  box const c = container().expand(n, n);
  return c - *this;
}

// cost: O(n^2) in general, but only O(n) for shifting
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::expand(const vect<T, D> &lo,
                                    const vect<T, D> &hi) const {
  bboxset res;
  {
    SKIP_NORMALIZE(res);
    if (all(lo == -hi)) {
      // Special case for shifting, since this is faster
      for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
        res += (*bi).expand(lo, hi);
      }
    } else {
      // We don't know (yet?) how to shrink a set
      assert(all(lo >= 0 and hi >= 0));
      for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
        res |= (*bi).expand(lo, hi);
      }
    }
  }
  return res;
}

// cost: O(n^2) in general, but only O(n) for shifting
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::expand(const vect<T, D> &lo, const vect<T, D> &hi,
                                    const vect<T, D> &denom) const {
  assert(all(denom > vect<T, D>(0)));
  bboxset res;
  {
    SKIP_NORMALIZE(res);
    if (all(lo == -hi)) {
      // Special case for shifting, since this is faster
      for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
        res += (*bi).expand(lo, hi, denom);
      }
    } else {
      // We don't know (yet?) how to shrink a set
      assert(all((lo >= 0 and hi >= 0) or (lo == hi)));
      for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
        res |= (*bi).expand(lo, hi, denom);
      }
    }
  }
  return res;
}

template <typename T, int D>
bboxset<T, D> bboxset<T, D>::expanded_for(const box &b) const {
  bboxset res;
  {
    SKIP_NORMALIZE(res);
    for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
      res |= (*bi).expanded_for(b);
    }
  }
  return res;
}

// TODO: this is incorrect
#if 1
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::contracted_for(const box &b) const {
  bboxset res;
  {
    SKIP_NORMALIZE(res);
    for (const_iterator bi = begin(), be = end(); bi != be; ++bi) {
      res |= (*bi).contracted_for(b);
    }
  }
  return res;
}
#endif

// Equality
template <typename T, int D>
bool bboxset<T, D>::operator<=(const bboxset<T, D> &s) const {
  return (*this - s).empty();
}

template <typename T, int D>
bool bboxset<T, D>::operator<(const bboxset<T, D> &s) const {
  return (*this - s).empty() && not(s - *this).empty();
}

template <typename T, int D>
bool bboxset<T, D>::operator>=(const bboxset<T, D> &s) const {
  return s <= *this;
}

template <typename T, int D>
bool bboxset<T, D>::operator>(const bboxset<T, D> &s) const {
  return s < *this;
}

template <typename T, int D>
bool bboxset<T, D>::operator==(const bboxset<T, D> &s) const {
  return (*this <= s) && (*this >= s);
}

template <typename T, int D>
bool bboxset<T, D>::operator!=(const bboxset<T, D> &s) const {
  return not(*this == s);
}

// Serialise
template <typename T, int D>
template <typename C>
void bboxset<T, D>::serialise(C &out) const {
  for (const_iterator it = begin(); it != end(); ++it) {
    out.insert(out.end(), *it);
  }
}

// Input
template <typename T, int D> istream &bboxset<T, D>::input(istream &is) {
  T Tdummy;
  try {
    skipws(is);
    consume(is, "bboxset<");
    consume(is, typestring(Tdummy));
    consume(is, ",");
    int D_;
    is >> D_;
    if (D_ != D) {
      cout << "Input error: Wrong bboxset dimension " << D_ << ", expected "
           << D << endl;
      throw input_error();
    }
    consume(is, ">:{");
    consume(is, "size=");
    size_type size_;
    is >> size_;
    consume(is, ",");
    consume(is, "setsize=");
    int setsize_;
    is >> setsize_;
    consume(is, ",");
    consume(is, "set=");
    is >> bs;
    consume(is, "}");
  } catch (input_error &err) {
    cout << "Input error while reading a bboxset<" << typestring(Tdummy) << ","
         << D << ">" << endl;
    throw err;
  }
  normalize();
  return is;
}

// Output
template <typename T, int D> ostream &bboxset<T, D>::output(ostream &os) const {
  T Tdummy;
  os << "bboxset<" << typestring(Tdummy) << "," << D << ">:{"
     << "size=" << size() << ","
     << "setsize=" << setsize() << ","
     << "set=" << bs << "}";
  return os;
}

#undef SKIP_NORMALIZE

} // namespace bboxset1

#endif // BBOXSET1_HH
