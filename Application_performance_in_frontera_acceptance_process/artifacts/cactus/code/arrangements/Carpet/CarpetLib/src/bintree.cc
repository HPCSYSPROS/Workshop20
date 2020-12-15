#include <algorithm>
#include <cmath>
#include <iostream>

#include "defs.hh"
#include "region.hh"

#include "bintree.hh"

// Create an empty tree
template <typename T, int D, typename P>
bintree<T, D, P>::bintree()
    : type(type_empty) {
  assert(invariant());
  // This is unused
  assert(0);
}

// Create a tree branch from a list of bounds and subtrees
template <typename T, int D, typename P>
bintree<T, D, P>::bintree(int const dir_, vect<T, 3> const &bounds_,
                          vect<bintree *, 2> const &subtrees_)
    : type(type_branch), dir(dir_), bounds(bounds_), subtrees(subtrees_) {
  assert(dir >= 0 and dir < D);
  assert(invariant());
}

// Create a tree leaf from a payload
template <typename T, int D, typename P>
bintree<T, D, P>::bintree(P const &p_)
    : type(type_leaf), p(p_) {
  assert(invariant());
}

// Create a tree as copy from another tree
template <typename T, int D, typename P>
bintree<T, D, P>::bintree(bintree const &t)
    : type(t.type) {
  switch (type) {
  case type_empty:
    // do nothing
    break;
  case type_branch:
    dir = t.dir;
    bounds = t.bounds;
    for (int i = 0; i < 2; ++i) {
      subtrees[i] = new bintree(*t.subtrees[i]);
    }
    break;
  case type_leaf:
    p = t.p;
    break;
  default:
    assert(0);
  }
  assert(invariant());
}

// Copy a tree
template <typename T, int D, typename P>
bintree<T, D, P> &bintree<T, D, P>::operator=(bintree const &t) {
  assert(&t !=
         this); // subtree delet handling is currently incorrect in this case
  if (&t == this)
    return *this; // nothing to do

  assert(invariant());
  if (is_branch()) {
    for (int i = 0; i < 2; ++i) {
      delete subtrees[i];
    }
  }
  type = t.type;
  switch (type) {
  case type_empty:
    // do nothing
    break;
  case type_branch:
    dir = t.dir;
    bounds = t.bounds;
    for (int i = 0; i < 2; ++i) {
      subtrees[i] = new bintree(*t.subtrees[i]);
    }
    break;
  case type_leaf:
    p = t.p;
    break;
  default:
    assert(0);
  }
  assert(invariant());
  return *this;
}

// Delete a tree (and its subtrees)
template <typename T, int D, typename P> bintree<T, D, P>::~bintree() {
  assert(invariant());
  if (is_branch()) {
    for (int i = 0; i < 2; ++i) {
      delete subtrees[i];
    }
  }
}

// Compare trees
template <typename T, int D, typename P>
bool bintree<T, D, P>::operator==(bintree const &t) const {
  assert(invariant());
  assert(t.invariant());
  if (type != t.type)
    return false;
  switch (type) {
  case type_empty:
    break;
  case type_branch:
    if (dir != t.dir)
      return false;
    if (any(bounds != t.bounds))
      return false;
    if (any(subtrees != t.subtrees))
      return false;
    break;
  case type_leaf:
    return p == t.p;
  default:
    assert(0);
  }
  return true;
}

// Invariant
template <typename T, int D, typename P>
bool bintree<T, D, P>::invariant() const {
  return empty() + is_branch() + is_leaf() == 1;
}

// Find the leaf payload corresponding to a position
template <typename T, int D, typename P>
P const *bintree<T, D, P>::search(tvect const &ipos) const {
  assert(not empty());
  if (is_leaf())
    return &p;
  int const i = ::asearch(ipos[dir], bounds);
  if (i < 0 or i >= 2)
    return NULL; // not found
  return subtrees[i]->search(ipos);
}

template <typename T, int D, typename P>
P *bintree<T, D, P>::search(tvect const &ipos) {
  assert(not empty());
  if (is_leaf())
    return &p;
  int const i = ::asearch(ipos[dir], bounds);
  if (i < 0 or i >= 2)
    return NULL; // not found
  return subtrees[i]->search(ipos);
}

// Const iterator
template <typename T, int D, typename P>
bintree<T, D, P>::const_iterator::const_iterator(bintree const &f_)
    : f(f_), i(0), it(0) {
  if (f.is_branch()) {
    it = new const_iterator(*f.subtrees[i]);
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new const_iterator(*f.subtrees[i]);
    }
    assert(done() or not(*it).done());
  }
}

template <typename T, int D, typename P>
bintree<T, D, P>::const_iterator::const_iterator(bintree const &f_, int)
    : f(f_), it(0) {
  if (f.empty()) {
    i = 0;
  } else if (f.is_leaf()) {
    i = 1;
  } else {
    i = 2;
  }
}

template <typename T, int D, typename P>
bintree<T, D, P>::const_iterator::~const_iterator() {
  if (it) {
    delete it;
  }
}

template <typename T, int D, typename P>
bintree<T, D, P> const &bintree<T, D, P>::const_iterator::operator*() const {
  assert(not done());
  assert(not f.empty());
  if (f.is_leaf()) {
    return f;
  } else {
    assert(it);
    return **it;
  }
}

template <typename T, int D, typename P>
bool bintree<T, D, P>::const_iterator::
operator==(const_iterator const &it2) const {
  // assert (f == it2.f);
  assert(&f == &it2.f);
  if (i != it2.i)
    return false;
  if (it == 0 and it2.it == 0)
    return true;
  if (it == 0 or it2.it == 0)
    return false;
  return *it == *it2.it;
}

template <typename T, int D, typename P>
typename bintree<T, D, P>::const_iterator &bintree<T, D, P>::const_iterator::
operator++() {
  assert(not done());
  assert(not f.empty());
  if (f.is_leaf()) {
    ++i;
  } else {
    ++*it;
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new const_iterator(*f.subtrees[i]);
    }
    assert(done() or not(*it).done());
  }
  return *this;
}

template <typename T, int D, typename P>
bool bintree<T, D, P>::const_iterator::done() const {
  if (f.empty()) {
    return true;
  } else if (f.is_leaf()) {
    return i > 0;
  } else {
    return i == 2;
  }
}

// Non-const iterator
template <typename T, int D, typename P>
bintree<T, D, P>::iterator::iterator(bintree &f_)
    : f(f_), i(0), it(0) {
  if (f.is_branch()) {
    it = new iterator(*f.subtrees[i]);
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new iterator(*f.subtrees[i]);
    }
    assert(done() or not(*it).done());
  }
}

template <typename T, int D, typename P>
bintree<T, D, P>::iterator::iterator(bintree &f_, int)
    : f(f_), it(0) {
  if (f.empty()) {
    i = 0;
  } else if (f.is_leaf()) {
    i = 1;
  } else {
    i = 2;
  }
}

template <typename T, int D, typename P>
bintree<T, D, P>::iterator::~iterator() {
  if (it) {
    delete it;
  }
}

template <typename T, int D, typename P>
bintree<T, D, P> &bintree<T, D, P>::iterator::operator*() const {
  assert(not done());
  assert(not f.empty());
  if (f.is_leaf()) {
    return f;
  } else {
    assert(it);
    return **it;
  }
}

template <typename T, int D, typename P>
bool bintree<T, D, P>::iterator::operator==(iterator const &it2) const {
  // assert (f == it2.f);
  assert(&f == &it2.f);
  if (i != it2.i)
    return false;
  if (it == 0 and it2.it == 0)
    return true;
  if (it == 0 or it2.it == 0)
    return false;
  return *it == *it2.it;
}

template <typename T, int D, typename P>
typename bintree<T, D, P>::iterator &bintree<T, D, P>::iterator::operator++() {
  assert(not done());
  assert(not f.empty());
  if (f.is_leaf()) {
    ++i;
  } else {
    ++*it;
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new iterator(*f.subtrees[i]);
    }
    assert(done() or not(*it).done());
  }
  return *this;
}

template <typename T, int D, typename P>
bool bintree<T, D, P>::iterator::done() const {
  if (f.empty()) {
    return true;
  } else if (f.is_leaf()) {
    return i > 0;
  } else {
    return i == 2;
  }
}

// Memory usage
template <typename T, int D, typename P>
size_t bintree<T, D, P>::memory() const {
  size_t size = sizeof *this;
  if (is_branch()) {
    for (int i = 0; i < 2; ++i) {
      size += memoryof(*subtrees[i]);
    }
  }
  return size;
}

// Output helper
template <typename T, int D, typename P>
void bintree<T, D, P>::output(ostream &os) const {
  os << "bintree{";
  if (empty()) {
    os << "empty";
  } else if (is_branch()) {
    os << "branch:"
       << "dir=" << dir << ","
       << "subtrees=[";
    for (int i = 0; i < 2; ++i) {
      os << bounds[i] << ":[" << i << "]=" << *subtrees[i] << ":";
    }
    os << bounds[2] << "]";
  } else {
    os << "leaf:"
       << "payload=" << p;
  }
  os << "}";
}

// Generic arithmetic search
template <typename T, int D>
static int asearch(T const t, vect<T, D> const &ts) {
  int imin = 0;
  int imax = D - 1;
  if (imax <= imin) {
    return imin;
  }
  T tmin = ts[imin];
  if (t < tmin) {
    return -1;
  }
  T tmax = ts[imax];
  if (t >= tmax) {
    return imax;
  }
  int isize = imax - imin;
  for (;;) {
    // check loop invariants
    assert(imin < imax);
    assert(t >= tmin);
    assert(t < tmax);
    if (imax == imin + 1) {
      return imin;
    }
    assert(tmax > tmin); // require that ts is strictly ordered
    CCTK_REAL const rguess =
        (imax - imin) * CCTK_REAL(t - tmin) / (tmax - tmin);
    int const iguess = imin + max(1, int(floor(rguess)));
    // handle round-off errors
    if (iguess == imax) {
      return imax - 1;
    }
    assert(iguess > imin and iguess < imax);
    T const tguess = ts[iguess];
    if (t < tguess) {
      imax = iguess;
      tmax = tguess;
    } else {
      imin = iguess;
      tmin = tguess;
    }
    // check loop variant
    int const newisize = imax - imin;
    assert(newisize < isize);
    isize = newisize;
  }
}

template class bintree<int, dim, pseudoregion_t>;

template size_t memoryof(bintree<int, dim, pseudoregion_t> const &f);

template ostream &operator<<(ostream &os,
                             bintree<int, dim, pseudoregion_t> const &f);
