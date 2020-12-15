#include <algorithm>
#include <cmath>
#include <iostream>

#include "defs.hh"
#include "region.hh"

#include "fulltree.hh"

// Create an empty tree
template <typename T, int D, typename P>
fulltree<T, D, P>::fulltree()
    : type(type_empty) {
  assert(invariant());
  // This is unused
  assert(0);
}

// Create a tree branch from a list of bounds and subtrees
template <typename T, int D, typename P>
fulltree<T, D, P>::fulltree(int const dir_, vector<T> const &bounds_,
                            vector<fulltree *> const &subtrees_)
    : type(type_branch), dir(dir_), bounds(bounds_), subtrees(subtrees_) {
  assert(dir >= 0 and dir < D);
  assert(bounds.size() == subtrees.size() + 1);
  assert(invariant());
}

// Create a tree leaf from a payload
template <typename T, int D, typename P>
fulltree<T, D, P>::fulltree(P const &p_)
    : type(type_leaf), p(p_) {
  assert(invariant());
}

// Create a tree as copy from another tree
template <typename T, int D, typename P>
fulltree<T, D, P>::fulltree(fulltree const &t)
    : type(t.type) {
  switch (type) {
  case type_empty:
    // do nothing
    break;
  case type_branch:
    dir = t.dir;
    bounds = t.bounds;
    subtrees.resize(t.subtrees.size());
    for (size_t i = 0; i < subtrees.size(); ++i) {
      subtrees.AT(i) = new fulltree(*t.subtrees.AT(i));
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
fulltree<T, D, P> &fulltree<T, D, P>::operator=(fulltree const &t) {
  assert(&t !=
         this); // subtree delet handling is currently incorrect in this case
  if (&t == this)
    return *this; // nothing to do

  assert(invariant());
  if (is_branch()) {
    for (size_t i = 0; i < subtrees.size(); ++i) {
      delete subtrees.AT(i);
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
    subtrees.resize(t.subtrees.size());
    for (size_t i = 0; i < subtrees.size(); ++i) {
      subtrees.AT(i) = new fulltree(*t.subtrees.AT(i));
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
template <typename T, int D, typename P> fulltree<T, D, P>::~fulltree() {
  assert(invariant());
  if (is_branch()) {
    for (size_t i = 0; i < subtrees.size(); ++i) {
      delete subtrees.AT(i);
    }
  }
}

// Compare trees
template <typename T, int D, typename P>
bool fulltree<T, D, P>::operator==(fulltree const &t) const {
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
    if (bounds != t.bounds)
      return false;
    if (subtrees != t.subtrees)
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
bool fulltree<T, D, P>::invariant() const {
  return empty() + is_branch() + is_leaf() == 1;
}

// Find the leaf payload corresponding to a position
template <typename T, int D, typename P>
P const *fulltree<T, D, P>::search(tvect const &ipos) const {
  assert(not empty());
  // cout << "fulltree::search ipos=" << ipos << endl;
  if (is_leaf())
    return &p;
  int const i = ::asearch(ipos[dir], bounds);
  // cout << "fulltree::search i=" << i << " size=" << subtrees.size() << endl;
  if (i < 0 or i >= int(subtrees.size()))
    return NULL; // not found
  return subtrees.AT(i)->search(ipos);
}

template <typename T, int D, typename P>
P *fulltree<T, D, P>::search(tvect const &ipos) {
  assert(not empty());
  // cout << "fulltree::search ipos=" << ipos << endl;
  if (is_leaf())
    return &p;
  int const i = ::asearch(ipos[dir], bounds);
  // cout << "fulltree::search i=" << i << " size=" << subtrees.size() << endl;
  if (i < 0 or i >= int(subtrees.size()))
    return NULL; // not found
  return subtrees.AT(i)->search(ipos);
}

// Const iterator
template <typename T, int D, typename P>
fulltree<T, D, P>::const_iterator::const_iterator(fulltree const &f_)
    : f(f_), i(0), it(0) {
  if (f.is_branch()) {
    assert(f.subtrees.size() > 0);
    it = new const_iterator(*f.subtrees.AT(i));
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new const_iterator(*f.subtrees.AT(i));
    }
    assert(done() or not(*it).done());
  }
}

template <typename T, int D, typename P>
fulltree<T, D, P>::const_iterator::const_iterator(fulltree const &f_, int)
    : f(f_), it(0) {
  if (f.empty()) {
    i = 0;
  } else if (f.is_leaf()) {
    i = 1;
  } else {
    i = f.subtrees.size();
  }
}

template <typename T, int D, typename P>
fulltree<T, D, P>::const_iterator::~const_iterator() {
  if (it) {
    delete it;
  }
}

template <typename T, int D, typename P>
fulltree<T, D, P> const &fulltree<T, D, P>::const_iterator::operator*() const {
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
bool fulltree<T, D, P>::const_iterator::
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
typename fulltree<T, D, P>::const_iterator &fulltree<T, D, P>::const_iterator::
operator++() {
  assert(not done());
  assert(not f.empty());
  if (f.is_leaf()) {
    ++i;
  } else {
    ++*it;
#if 0
    if ((*it).done()) {
      delete it;
      it = 0;
      ++ i;
      if (not done()) {
        // to do: use a new function "reset iterator" instead
        it = new const_iterator (* f.subtrees.AT(i));
        assert (not (*it).done());
      }
    }
#endif
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new const_iterator(*f.subtrees.AT(i));
    }
    assert(done() or not(*it).done());
  }
  return *this;
}

template <typename T, int D, typename P>
bool fulltree<T, D, P>::const_iterator::done() const {
  if (f.empty()) {
    return true;
  } else if (f.is_leaf()) {
    return i > 0;
  } else {
    return i == f.subtrees.size();
  }
}

// Non-const iterator
template <typename T, int D, typename P>
fulltree<T, D, P>::iterator::iterator(fulltree &f_)
    : f(f_), i(0), it(0) {
  if (f.is_branch()) {
    assert(f.subtrees.size() > 0);
    it = new iterator(*f.subtrees.AT(i));
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new iterator(*f.subtrees.AT(i));
    }
    assert(done() or not(*it).done());
  }
}

template <typename T, int D, typename P>
fulltree<T, D, P>::iterator::iterator(fulltree &f_, int)
    : f(f_), it(0) {
  if (f.empty()) {
    i = 0;
  } else if (f.is_leaf()) {
    i = 1;
  } else {
    i = f.subtrees.size();
  }
}

template <typename T, int D, typename P>
fulltree<T, D, P>::iterator::~iterator() {
  if (it) {
    delete it;
  }
}

template <typename T, int D, typename P>
fulltree<T, D, P> &fulltree<T, D, P>::iterator::operator*() const {
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
bool fulltree<T, D, P>::iterator::operator==(iterator const &it2) const {
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
typename fulltree<T, D, P>::iterator &fulltree<T, D, P>::iterator::
operator++() {
  assert(not done());
  assert(not f.empty());
  if (f.is_leaf()) {
    ++i;
  } else {
    ++*it;
#if 0
    if ((*it).done()) {
      delete it;
      it = 0;
      ++ i;
      if (not done()) {
        // to do: use a new function "reset iterator" instead
        it = new iterator (* f.subtrees.AT(i));
        assert (not (*it).done());
      }
    }
#endif
    while ((*it).done()) {
      delete it;
      it = 0;
      ++i;
      if (done())
        break;
      // to do: use a new function "reset iterator" instead
      it = new iterator(*f.subtrees.AT(i));
    }
    assert(done() or not(*it).done());
  }
  return *this;
}

template <typename T, int D, typename P>
bool fulltree<T, D, P>::iterator::done() const {
  if (f.empty()) {
    return true;
  } else if (f.is_leaf()) {
    return i > 0;
  } else {
    return i == f.subtrees.size();
  }
}

// Memory usage
template <typename T, int D, typename P>
size_t fulltree<T, D, P>::memory() const {
  size_t size = sizeof *this;
  if (is_branch()) {
    size += memoryof(bounds) + memoryof(subtrees);
    for (typename vector<fulltree *>::const_iterator i = subtrees.begin();
         i != subtrees.end(); ++i) {
      size += memoryof(**i);
    }
  }
  return size;
}

// Output helper
template <typename T, int D, typename P>
void fulltree<T, D, P>::output(ostream &os) const {
  os << "fulltree{";
  if (empty()) {
    os << "empty";
  } else if (is_branch()) {
    os << "branch:"
       << "dir=" << dir << ","
       << "subtrees=[";
    for (size_t i = 0; i < subtrees.size(); ++i) {
      os << bounds.AT(i) << ":[" << i << "]=" << *subtrees.AT(i) << ":";
    }
    os << bounds.AT(subtrees.size()) << "]";
  } else {
    os << "leaf:"
       << "payload=" << p;
  }
  os << "}";
}

// Generic arithmetic search
template <typename T> static int asearch(T const t, vector<T> const &ts) {
  // cout << "fulltree::asearch t=" << t << " ts=" << ts << endl;
  int imin = 0;
  int imax = int(ts.size()) - 1;
  if (imax <= imin) {
    // cout << "fulltree::asearch: no values" << endl;
    return imin;
  }
  T tmin = ts.AT(imin);
  // cout << "fulltree::asearch: imin=" << imin << " tmin=" << tmin << endl;
  if (t < tmin) {
    // cout << "fulltree::asearch: below minimum" << endl;
    return -1;
  }
  T tmax = ts.AT(imax);
  // cout << "fulltree::asearch: imax=" << imax << " tmax=" << tmax << endl;
  if (t >= tmax) {
    // cout << "fulltree::asearch: above maximum" << endl;
    return imax;
  }
  int isize = imax - imin;
  for (;;) {
    // check loop invariants
    assert(imin < imax);
    assert(t >= tmin);
    assert(t < tmax);
    // cout << "fulltree::asearch t=" << t << " imin=" << imin << " imax=" <<
    // imax << endl;
    if (imax == imin + 1) {
      // cout << "fulltree::asearch: found value" << endl;
      return imin;
    }
    assert(tmax > tmin); // require that ts is strictly ordered
    CCTK_REAL const rguess =
        (imax - imin) * CCTK_REAL(t - tmin) / (tmax - tmin);
    int const iguess = imin + max(1, int(floor(rguess)));
    // handle round-off errors
    if (iguess == imax) {
      // cout << "fulltree::asearch: accidental hit after roundoff error" <<
      // endl;
      return imax - 1;
    }
    assert(iguess > imin and iguess < imax);
    T const tguess = ts.AT(iguess);
    // cout << "fulltree::asearch: intersecting at iguess=" << iguess << "
    // tguess=" << tguess << endl;
    if (t < tguess) {
      imax = iguess;
      tmax = tguess;
      // cout << "fulltree::asearch: new imax=" << imax << " tmax=" << tmax <<
      // endl;
    } else {
      imin = iguess;
      tmin = tguess;
      // cout << "fulltree::asearch: new imin=" << imin << " tmin=" << tmin <<
      // endl;
    }
    // check loop variant
    int const newisize = imax - imin;
    assert(newisize < isize);
    isize = newisize;
  }
}

template class fulltree<int, dim, pseudoregion_t>;

template size_t memoryof(fulltree<int, dim, pseudoregion_t> const &f);

template ostream &operator<<(ostream &os,
                             fulltree<int, dim, pseudoregion_t> const &f);
