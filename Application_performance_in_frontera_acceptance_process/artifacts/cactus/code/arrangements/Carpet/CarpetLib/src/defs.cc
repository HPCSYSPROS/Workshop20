#include <cctk.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cctype>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <vector>
#include <utility>

#include "bbox.hh"
#include "defs.hh"
#include "dh.hh"
#include "region.hh"
#include "vect.hh"

using namespace std;

template <typename T>
inline T ipow_helper(T x, unsigned int y) CCTK_ATTRIBUTE_CONST;
template <typename T> inline T ipow_helper(T x, unsigned int y) {
  T z = y & 1 ? x : 1;
  while (y >>= 1) {
    x *= x;
    if (y & 1)
      z *= x;
  }
  return z;
}

template <class T> T ipow(T const x, int const y) {
  if (y < 0)
    return T(1) / ipow_helper(x, -y);
  else
    return ipow_helper(x, y);
}

template <class T> int ilog(T const b, T const x) {
  assert(b > T(1));
  int r = 0;
  T y(1);
  // TODO: This algorithm is slow; use a more clever one
  while (y < x) {
    y *= b;
    ++r;
  }
  assert((x < b && r == 0) || (x >= b && ipow(b, r) <= x));
  assert(ipow(b, r + 1) > x);
  return r;
}

// Access to CarpetLib parameters
CCTK_INT get_poison_value() {
  DECLARE_CCTK_PARAMETERS;
  return poison_value;
}

CCTK_INT get_deadbeef() {
  DECLARE_CCTK_PARAMETERS;
  return deadbeef;
}

void skipws(istream &is) {
  while (is.good() and isspace(is.peek())) {
    is.get();
  }
}

void expect(istream &is, const char c) {
  if (is.peek() == c)
    return;
  cout << "While reading characters from a stream:" << endl
       << "   Character '" << c << "' expected, but not found." << endl
       << "   The next up to 100 available characters are \"";
  for (int i = 0; i < 100; ++i) {
    const int uc = is.get();
    if (uc < 0)
      break;
    cout << (unsigned char)uc;
  }
  cout << "\"." << endl;
  throw input_error();
}

void consume(istream &is, const char c) {
  expect(is, c);
  is.get();
}

void consume(istream &is, const char *const s) {
  for (const char *p = s; *p; ++p) {
    expect(is, *p);
    is.get();
  }
}

// Container memory usage
template <class T> size_t memoryof(list<T> const &c) {
  size_t s = sizeof c;
  for (typename list<T>::const_iterator i = c.begin(); i != c.end(); ++i) {
    // Assume that there are two pointers per list element, pointing
    // to the previous and next element, respectively
    s += 2 * sizeof(void *) + memoryof(*i);
  }
  return s;
}

template <class T> size_t memoryof(set<T> const &c) {
  size_t s = sizeof c;
  for (typename set<T>::const_iterator i = c.begin(); i != c.end(); ++i) {
    // Assume that there are three pointers per list element, forming
    // a tree structure
    s += 3 * sizeof(void *) + memoryof(*i);
  }
  return s;
}

template <class S, class T> size_t memoryof(map<S, T> const &c) {
  size_t s = sizeof c;
  for (typename map<S, T>::const_iterator i = c.begin(); i != c.end(); ++i) {
    // Assume that there are three pointers per list element, forming
    // a tree structure
    s += 3 * sizeof(void *) + memoryof(i->second);
  }
  return s;
}

template <class T> size_t memoryof(stack<T> const &c) {
  size_t s = sizeof c;
#if 0
  for (typename stack<T>::const_iterator i=c.begin(); i!=c.end(); ++i) {
    // Assume that a stack is stored in an efficient manner
    s += memoryof(*i);
  }
#endif
  // Cannot access elements...  a stack is LIFO!
  s += c.size() * sizeof(T);
  return s;
}

template <class T> size_t memoryof(vector<T> const &c) {
  size_t s = sizeof c;
  for (typename vector<T>::const_iterator i = c.begin(); i != c.end(); ++i) {
    // Vectors are stored contiguously
    s += memoryof(*i);
  }
  return s;
}

// List input
template <class T> istream &input(istream &is, list<T> &l) {
  l.clear();
  try {
    skipws(is);
    consume(is, '[');
    skipws(is);
    while (is.good() and is.peek() != ']') {
      T elem;
      is >> elem;
      l.push_back(elem);
      skipws(is);
      if (is.peek() != ',')
        break;
      is.get();
      skipws(is);
    }
    skipws(is);
    consume(is, ']');
  } catch (input_error &err) {
    cout << "Input error while reading a list<>" << endl
         << "   The following elements have been read so far: " << l << endl;
    throw err;
  }
  return is;
}

// Set input
template <class T> istream &input(istream &is, set<T> &s) {
  s.clear();
  try {
    skipws(is);
    consume(is, '{');
    skipws(is);
    while (is.good() and is.peek() != '}') {
      T elem;
      is >> elem;
      s.insert(elem);
      skipws(is);
      if (is.peek() != ',')
        break;
      is.get();
      skipws(is);
    }
    skipws(is);
    consume(is, ']');
  } catch (input_error &err) {
    cout << "Input error while reading a set<>" << endl
         << "   The following elements have been read so far: " << s << endl;
    throw err;
  }
  return is;
}

// Vector input
template <class T> istream &input(istream &is, vector<T> &v) {
  v.clear();
  try {
    skipws(is);
    consume(is, '[');
    skipws(is);
    while (is.good() and is.peek() != ']') {
      T elem;
      is >> elem;
      v.push_back(elem);
      skipws(is);
      if (is.peek() != ',')
        break;
      is.get();
      skipws(is);
    }
    skipws(is);
    consume(is, ']');
  } catch (input_error &err) {
    cout << "Input error while reading a vector<>" << endl
         << "   The following elements have been read so far: " << v << endl;
    throw err;
  }
  return is;
}

// List output
template <class T> ostream &output(ostream &os, const list<T> &l) {
  os << "[";
  for (typename list<T>::const_iterator ti = l.begin(); ti != l.end(); ++ti) {
    if (ti != l.begin())
      os << ",";
    os << *ti;
  }
  os << "]";
  return os;
}

// Map output
template <class S, class T> ostream &output(ostream &os, const map<S, T> &m) {
  os << "{";
  for (typename map<S, T>::const_iterator ti = m.begin(); ti != m.end(); ++ti) {
    if (ti != m.begin())
      os << ",";
    os << ti->first << ":" << ti->second;
  }
  os << "}";
  return os;
}

// Pair output
template <class S, class T> ostream &output(ostream &os, const pair<S, T> &p) {
  os << "(" << p.first << "," << p.second << ")";
  return os;
}

// Set output
template <class T> ostream &output(ostream &os, const set<T> &s) {
  os << "{";
  for (typename set<T>::const_iterator ti = s.begin(); ti != s.end(); ++ti) {
    if (ti != s.begin())
      os << ",";
    os << *ti;
  }
  os << "}";
  return os;
}

#ifdef CARPET_ENABLE_BBOXSET2
// Shared pointer output
template <class T> ostream &output(ostream &os, const shared_ptr<T> &s) {
  return os << "(&" << *s.get() << ")";
}
#endif

// Stack output
template <class T> ostream &output(ostream &os, const stack<T> &s) {
  stack<T> s2(s);
  list<T> l;
  while (not s2.empty()) {
    l.insert(l.begin(), s2.top());
    s2.pop();
  }
  return output(os, l);
}

// Vector output
template <class T> ostream &output(ostream &os, const vector<T> &v) {
  os << "[";
  // Do not number the elements, as this would lead to a format that
  // cannot be read back in.
  //   int cnt=0;
  for (typename vector<T>::const_iterator ti = v.begin(); ti != v.end(); ++ti) {
    if (ti != v.begin())
      os << ",";
    //     os << cnt++ << ":";
    os << *ti;
  }
  os << "]";
  return os;
}

#include "bbox.hh"
#include "bboxset.hh"
#include "dh.hh"
#include "fulltree.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "region.hh"
#include "th.hh"
#include "vect.hh"

template int ipow(int x, int y);
template CCTK_REAL ipow(CCTK_REAL x, int y);
// template vect<int,dim> ipow (vect<int,dim> x, int y);
template vect<CCTK_REAL, dim> ipow(vect<CCTK_REAL, dim> x, int y);
template int ilog(int x, int y);
template int ilog(CCTK_REAL x, CCTK_REAL y);

template size_t memoryof(rvect const &v);
// template size_t memoryof (list<ibbox> const & l);
// template size_t memoryof (list<ivect> const & l);
template size_t memoryof(set<dh *> const &l);
template size_t memoryof(set<gh *> const &l);
template size_t memoryof(set<gdata *> const &l);
template size_t memoryof(set<ggf *> const &l);
template size_t memoryof(map<int, ggf *> const &l);
template size_t memoryof(set<th *> const &l);
template size_t memoryof(stack<void *> const &s);
template size_t memoryof(vector<bool> const &v);
template size_t memoryof(vector<int> const &v);
// template size_t memoryof (vector<CCTK_REAL> const & v);
template size_t memoryof(vector<dh *> const &v);
template size_t memoryof(vector<gh *> const &v);
template size_t memoryof(vector<bbox<int, 1> > const &v);
template size_t memoryof(vector<bbox<int, 2> > const &v);
template size_t memoryof(vector<bbox<int, 3> > const &v);
template size_t memoryof(vector<rbbox> const &v);
template size_t memoryof(vector<bboxset1::bboxset<int, dim> > const &v);
#ifdef CARPET_ENABLE_BBOXSET2
template size_t memoryof(vector<bboxset2::bboxset<int, dim> > const &v);
#endif
template size_t memoryof(vector<ivect> const &v);
template size_t memoryof(vector<i2vect> const &v);
template size_t memoryof(vector<rvect> const &v);
template size_t memoryof(vector<fulltree<int, dim, pseudoregion_t> *> const &f);
// template size_t memoryof (vector<pseudoregion_t> const & v);
// template size_t memoryof (vector<region_t> const & v);
template size_t memoryof(vector<sendrecv_pseudoregion_t> const &v);
template size_t memoryof(vector<vector<int> > const &v);
template size_t memoryof(vector<vector<CCTK_REAL> > const &v);
template size_t memoryof(vector<vector<ibbox> > const &v);
template size_t memoryof(vector<vector<rvect> > const &v);
template size_t memoryof(vector<vector<dh::fast_dboxes> > const &v);
// template size_t memoryof (vector<vector<dh::full_dboxes> > const & v);
template size_t memoryof(vector<vector<dh::level_dboxes> > const &v);
// template size_t memoryof (vector<vector<dh::light_dboxes> > const & v);
// template size_t memoryof (vector<vector<dh::local_dboxes> > const & v);
// template size_t memoryof (vector<vector<region_t> > const & v);
template size_t memoryof(vector<vector<vector<CCTK_REAL> > > const &v);
// template size_t memoryof (vector<vector<vector<dh::fast_dboxes> > > const &
// v);
// template size_t memoryof (vector<vector<vector<dh::full_dboxes> > > const &
// v);
// template size_t memoryof (vector<vector<vector<dh::level_dboxes> > > const &
// v);
template size_t memoryof(vector<vector<vector<dh::light_dboxes> > > const &v);
template size_t memoryof(vector<vector<vector<dh::local_dboxes> > > const &v);
template size_t memoryof(vector<vector<vector<region_t> > > const &v);
// template size_t memoryof (vector<vector<vector<gdata*> > > const & v);
template size_t memoryof(vector<vector<vector<vector<gdata *> > > > const &v);

// template istream& input (istream& os, list<ibbox>& l);
template istream &input(istream &os, set<bbox<int, 1> > &s);
template istream &input(istream &os, set<bbox<int, 2> > &s);
template istream &input(istream &os, set<bbox<int, 3> > &s);
template istream &input(istream &os, vector<int> &v);
template istream &input(istream &os, vector<CCTK_REAL> &v);
template istream &input(istream &os, vector<bbox<int, 1> > &v);
template istream &input(istream &os, vector<bbox<int, 2> > &v);
template istream &input(istream &os, vector<bbox<int, 3> > &v);
// template istream& input (istream& os, vector<rbbox>& v);
template istream &input(istream &os, vector<bboxset1::bboxset<int, dim> > &v);
#ifdef CARPET_ENABLE_BBOXSET2
template istream &input(istream &os, vector<bboxset2::bboxset<int, dim> > &v);
#endif
template istream &input(istream &os, vector<ivect> &v);
template istream &input(istream &os, vector<bbvect> &v);
template istream &input(istream &os, vector<i2vect> &v);
// template istream& input (istream& os, vector<region_t>& v);
// template istream& input (istream& os, vector<pseudoregion_t>& v);
template istream &input(istream &os, vector<sendrecv_pseudoregion_t> &v);
template istream &input(istream &os, vector<vector<int> > &v);
template istream &input(istream &os, vector<vector<CCTK_REAL> > &v);
template istream &input(istream &os, vector<vector<ibbox> > &v);
template istream &input(istream &os, vector<vector<rbbox> > &v);
template istream &input(istream &os, vector<vector<bbvect> > &v);
template istream &input(istream &os, vector<vector<i2vect> > &v);
template istream &input(istream &os, vector<vector<region_t> > &v);
template istream &input(istream &os, vector<vector<vector<CCTK_REAL> > > &v);
template istream &input(istream &os, vector<vector<vector<region_t> > > &v);

// template ostream& output (ostream& os, const list<ibbox>& l);
// template ostream& output (ostream& os, const list<region_t>& l);
template ostream &output(ostream &os, const pair<int const, ggf *> &p);
#ifdef CARPET_ENABLE_BBOXSET2
// template ostream& output (ostream& os, const
// map<int,shared_ptr<bboxset2::bboxset<int,0> > >& m);
// template ostream& output (ostream& os, const
// map<int,shared_ptr<bboxset2::bboxset<int,1> > >& m);
// template ostream& output (ostream& os, const
// map<int,shared_ptr<bboxset2::bboxset<int,2> > >& m);
// template ostream& output (ostream& os, const
// map<int,shared_ptr<bboxset2::bboxset<int,3> > >& m);
#endif
// template ostream& output (ostream& os, const map<string,Carpet::Timer*>& m);
// template ostream& output (ostream& os, const pair<bbox<int,0>,int>& p);
// template ostream& output (ostream& os, const set<bbox<int,0> >& s);
template ostream &output(ostream &os, const set<bbox<int, 1> > &s);
template ostream &output(ostream &os, const set<bbox<int, 2> > &s);
template ostream &output(ostream &os, const set<bbox<int, 3> > &s);
// template ostream& output (ostream& os, const set<bboxset1::bboxset<int,dim>
// >& s);
#ifdef CARPET_ENABLE_BBOXSET2
// template ostream& output (ostream& os, const
// shared_ptr<bboxset2::bboxset<int,0> >& s);
#endif
// template ostream& output (ostream& os, const stack<ibbox>& s);
template ostream &output(ostream &os, const vector<bool> &v);
template ostream &output(ostream &os, const vector<int> &v);
template ostream &output(ostream &os, const vector<CCTK_REAL> &v);
template ostream &output(ostream &os, const vector<bbox<int, 1> > &v);
template ostream &output(ostream &os, const vector<bbox<int, 2> > &v);
template ostream &output(ostream &os, const vector<bbox<int, 3> > &v);
template ostream &output(ostream &os, const vector<rbbox> &v);
template ostream &output(ostream &os,
                         const vector<bboxset1::bboxset<int, dim> > &v);
#ifdef CARPET_ENABLE_BBOXSET2
template ostream &output(ostream &os,
                         const vector<bboxset2::bboxset<int, dim> > &v);
#endif
template ostream &output(ostream &os, const vector<ivect> &v);
template ostream &output(ostream &os, const vector<rvect> &v);
// template ostream& output (ostream& os, const vector<bbvect>& v);
template ostream &output(ostream &os, const vector<i2vect> &v);
// template ostream& output (ostream& os, const vector<dh::fast_dboxes> & v);
// template ostream& output (ostream& os, const vector<dh::full_dboxes> & v);
// template ostream& output (ostream& os, const vector<dh::level_dboxes> & v);
// template ostream& output (ostream& os, const vector<dh::light_dboxes> & v);
// template ostream& output (ostream& os, const vector<dh::local_dboxes> & v);
template ostream &output(ostream &os, const vector<region_t> &v);
// template ostream& output (ostream& os, const vector<pseudoregion_t>& v);
template ostream &output(ostream &os, const vector<sendrecv_pseudoregion_t> &v);
// template ostream& output (ostream& os, const vector<list<ibbox> >& v);
// template ostream& output (ostream& os, const vector<pair<bbox<int,0>,int> >&
// v);
// template ostream& output (ostream& os, const vector<pair<bbox<int,1>,int> >&
// v);
// template ostream& output (ostream& os, const vector<pair<bbox<int,2>,int> >&
// v);
// template ostream& output (ostream& os, const vector<pair<bbox<int,3>,int> >&
// v);
template ostream &output(ostream &os, const vector<vector<int> > &v);
template ostream &output(ostream &os, const vector<vector<CCTK_REAL> > &v);
template ostream &output(ostream &os, const vector<vector<ibbox> > &v);
template ostream &output(ostream &os, const vector<vector<rbbox> > &v);
template ostream &output(ostream &os, const vector<vector<bbvect> > &v);
template ostream &output(ostream &os, const vector<vector<i2vect> > &v);
template ostream &output(ostream &os,
                         const vector<vector<dh::fast_dboxes> > &b);
// template ostream& output (ostream& os, const vector<vector<dh::full_dboxes> >
// & b);
template ostream &output(ostream &os,
                         const vector<vector<dh::level_dboxes> > &b);
// template ostream& output (ostream& os, const vector<vector<dh::light_dboxes>
// > & b);
// template ostream& output (ostream& os, const vector<vector<dh::local_dboxes>
// > & b);
template ostream &output(ostream &os, const vector<vector<region_t> > &v);
template ostream &output(ostream &os,
                         const vector<vector<vector<CCTK_REAL> > > &v);
// template ostream& output (ostream& os, const vector<vector<vector<ibbox> > >&
// v);
template ostream &output(ostream &os, const vector<vector<vector<ibset> > > &v);
// template ostream& output (ostream& os, const
// vector<vector<vector<dh::fast_dboxes> > > & b);
// template ostream& output (ostream& os, const
// vector<vector<vector<dh::full_dboxes> > > & b);
// template ostream& output (ostream& os, const
// vector<vector<vector<dh::level_dboxes> > > & b);
template ostream &output(ostream &os,
                         const vector<vector<vector<dh::light_dboxes> > > &b);
template ostream &output(ostream &os,
                         const vector<vector<vector<dh::local_dboxes> > > &b);
template ostream &output(ostream &os,
                         const vector<vector<vector<region_t> > > &v);
