#ifndef FULLTREE_HH
#define FULLTREE_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <vect.hh>

using namespace std;

// This is a "full tree" data structure, i.e. a tree data structure
// which decomposes a cuboid domain into a set of non-overlapping
// cuboid subdomains.  It is an n-ary tree, i.e. each tree node can
// have arbitrarily many subtrees.  Each node splits a domain in
// exactly one direction.  Subdomains cannot be empty.

// All intervals are closed at the lower and open at the upper
// boundary.  This makes it possible to combine adjacent such
// intervals, obtaining again an interval with this property.  (In
// particular, if bboxes are stored, the upper bound of bboxes must be
// increased by the stride to arrive at such intervals.)

// Generic arithmetic search
template <typename T>
static int asearch(T t, vector<T> const &ts) CCTK_ATTRIBUTE_PURE;

// T: index space (usually integer, or maybe real)
// D: number of dimensions (rank)
// P: payload (e.g. region_t)
template <typename T, int D, typename P> class fulltree {

public:
  // Short name for a small vector
  typedef vect<T, D> tvect;

private:
  enum tree_type_t { type_empty, type_branch, type_leaf } type;

  // Direction in which the node is split (0 <= dir < D)
  int dir;

  // If this is a branch:
  // n+1 bounds, splitting the domain
  vector<T> bounds; // [n+1]
  // n pointers to subtrees
  vector<fulltree *> subtrees; // [n]

  // If this is a leaf:
  // just the payload
  P p;

public:
  // Create an empty tree
  fulltree();

  // Create a tree branch from a list of bounds and subtrees
  fulltree(int dir_, vector<T> const &bounds_,
           vector<fulltree *> const &subtrees_);

  // Create a tree leaf from a payload
  fulltree(P const &p_);

  // Create a tree as copy from another tree
  fulltree(fulltree const &t);

  // Copy a tree
  fulltree &operator=(fulltree const &t);

  // Delete a tree (and its subtrees)
  ~fulltree();

  // Inquire tree properties
  bool empty() const { return type == type_empty; }
  bool is_branch() const { return type == type_branch; }
  bool is_leaf() const { return type == type_leaf; }

  // Compare trees
  bool operator==(fulltree const &t) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator!=(fulltree const &t) const { return not(*this == t); }

  // Invariant
  bool invariant() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Access the payload
  P const &payload() const {
    assert(is_leaf());
    return p;
  }
  P &payload() {
    assert(is_leaf());
    return p;
  }

  // Find the leaf payload corresponding to a position
  P const *search(tvect const &ipos) const CCTK_MEMBER_ATTRIBUTE_PURE;
  P *search(tvect const &ipos) CCTK_MEMBER_ATTRIBUTE_PURE;

  class const_iterator {
    fulltree const &f;
    size_t i;
    const_iterator *it;

  public:
    const_iterator(fulltree const &f_);
    const_iterator(fulltree const &f_, int);
    ~const_iterator();

  private:
    const_iterator(const_iterator const &);
    const_iterator &operator=(const_iterator const &);

  public:
    fulltree const &operator*() const CCTK_MEMBER_ATTRIBUTE_PURE;
    bool operator==(const_iterator const &it2) const CCTK_MEMBER_ATTRIBUTE_PURE;
    bool operator!=(const_iterator const &it2) const {
      return not(*this == it2);
    }
    const_iterator &operator++();
    bool done() const CCTK_MEMBER_ATTRIBUTE_PURE;
  };

  class iterator {
    fulltree &f;
    size_t i;
    iterator *it;

  public:
    iterator(fulltree &f_);
    iterator(fulltree &f_, int);
    ~iterator();

  private:
    iterator(iterator const &);
    iterator &operator=(iterator const &);

  public:
    fulltree &operator*() const CCTK_MEMBER_ATTRIBUTE_PURE;
    bool operator==(iterator const &it2) const CCTK_MEMBER_ATTRIBUTE_PURE;
    bool operator!=(iterator const &it2) const { return not(*this == it2); }
    iterator &operator++();
    bool done() const CCTK_MEMBER_ATTRIBUTE_PURE;
  };

#if 0
  // ES 2008-09-04:
  // It seems that PGI does not handle comparisons to end() correctly.
  // We therefore disable these; please use the iterators'
  // constructors and their done() methods instead.
  
  const_iterator begin() const
  { return const_iterator (*this); }
  
  const_iterator end() const
  { return const_iterator (*this, 0); }
  
  iterator begin()
  { return iterator (*this); }
  
  iterator end()
  { return iterator (*this, 0); }
#endif

  // Memory usage
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Output helper
  void output(ostream &os) const;
};

// Memory usage
template <typename T, int D, typename P>
inline size_t memoryof(fulltree<T, D, P> const &f) {
  return f.memory();
}

// Output
template <typename T, int D, typename P>
ostream &operator<<(ostream &os, fulltree<T, D, P> const &f) {
  f.output(os);
  return os;
}

#endif // #ifndef FULLTREE_HH
