#ifndef BINTREE_HH
#define BINTREE_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <vect.hh>

using namespace std;

// This is a binary tree data structure, i.e. a tree data structure
// which decomposes a cuboid domain into two non-overlapping cuboid
// subdomains.  Each node splits a domain in exactly one direction.
// Subdomains cannot be empty.

// All intervals are closed at the lower and open at the upper
// boundary.  This makes it possible to combine adjacent such
// intervals, obtaining again an interval with this property.  (In
// particular, if bboxes are stored, the upper bound of bboxes must be
// increased by the stride to arrive at such intervals.)

// Generic arithmetic search
template <typename T, int D> static int asearch(T t, vect<T, D> const &ts);

// T: index space (usually integer, or maybe real)
// D: number of dimensions (rank)
// P: payload (e.g. region_t)
template <typename T, int D, typename P> class bintree {

public:
  // Short name for a small vector
  typedef vect<T, D> tvect;

private:
  enum tree_type_t { type_empty, type_branch, type_leaf } type;

  // Direction in which the node is split (0 <= dir < D)
  int dir;

  // If this is a branch:
  // 3 bounds, splitting the domain
  vect<T, 3> bounds;
  // 2 pointers to subtrees
  vect<bintree *, 2> subtrees;

  // If this is a leaf:
  // just the payload
  P p;

public:
  // Create an empty tree
  bintree();

  // Create a tree branch from a list of bounds and subtrees
  bintree(int dir_, vect<T, 3> const &bounds_,
          vect<bintree *, 2> const &subtrees_);

  // Create a tree leaf from a payload
  bintree(P const &p_);

  // Create a tree as copy from another tree
  bintree(bintree const &t);

  // Copy a tree
  bintree &operator=(bintree const &t);

  // Delete a tree (and its subtrees)
  ~bintree();

  // Inquire tree properties
  bool empty() const { return type == type_empty; }
  bool is_branch() const { return type == type_branch; }
  bool is_leaf() const { return type == type_leaf; }

  // Compare trees
  bool operator==(bintree const &t) const;
  bool operator!=(bintree const &t) const { return not(*this == t); }

  // Invariant
  bool invariant() const;

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
  P const *search(tvect const &ipos) const;
  P *search(tvect const &ipos);

  class const_iterator {
    bintree const &f;
    int i;
    const_iterator *it;

  public:
    const_iterator(bintree const &f_);
    const_iterator(bintree const &f_, int);
    ~const_iterator();

  private:
    const_iterator(const_iterator const &);
    const_iterator &operator=(const_iterator const &);

  public:
    bintree const &operator*() const;
    bool operator==(const_iterator const &it2) const;
    bool operator!=(const_iterator const &it2) const {
      return not(*this == it2);
    }
    const_iterator &operator++();
    bool done() const;
  };

  class iterator {
    bintree &f;
    int i;
    iterator *it;

  public:
    iterator(bintree &f_);
    iterator(bintree &f_, int);
    ~iterator();

  private:
    iterator(iterator const &);
    iterator &operator=(iterator const &);

  public:
    bintree &operator*() const;
    bool operator==(iterator const &it2) const;
    bool operator!=(iterator const &it2) const { return not(*this == it2); }
    iterator &operator++();
    bool done() const;
  };

  // Memory usage
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Output helper
  void output(ostream &os) const;
};

// Memory usage
template <typename T, int D, typename P>
inline size_t memoryof(bintree<T, D, P> const &f) {
  return f.memory();
}

// Output
template <typename T, int D, typename P>
ostream &operator<<(ostream &os, bintree<T, D, P> const &f) {
  f.output(os);
  return os;
}

#endif // #ifndef BINTREE_HH
