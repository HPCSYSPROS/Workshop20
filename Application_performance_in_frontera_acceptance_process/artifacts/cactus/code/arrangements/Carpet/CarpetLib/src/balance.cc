#include <cctk.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <limits>
#include <list>
#include <vector>

#include "region.hh"

#include "balance.hh"

using namespace std;

namespace CarpetLib {

// Interface for one item
struct item_ifc {
  CCTK_REAL load() const;
  item_ifc split(CCTK_REAL ratio_new_over_old);
};

// Declaration of template functions implementing an equivalent
// interface
template <typename item_t> CCTK_REAL item_load(item_t const &item);

template <typename item_t>
item_t item_split(item_t &item, CCTK_REAL ratio_new_over_old);

// Default definitions of these template functions using the
// interface above
template <typename item_t> CCTK_REAL item_load(item_t const &item) {
  return item.load();
}

template <typename item_t>
item_t item_split(item_t &item, CCTK_REAL ratio_new_over_old) {
  return item.split(ratio_new_over_old);
}

// A collection of items
template <typename item_t> class items_t {
  typedef list<item_t> coll_t;
  coll_t items;
  typename coll_t::iterator find_largest_item();

public:
  items_t();
  items_t(vector<item_t> const &items_);
  void add_item(item_t const &item);
  bool empty() const;
  size_t size() const;
  CCTK_REAL load() const;
  item_t &get_one_item();
  item_t &get_largest_item();
  item_t get_and_remove_largest_item();
  void copy_out(vector<item_t> &dst) const;
};

// A worker
template <typename item_t> class worker_t : public items_t<item_t> {
public:
  CCTK_REAL strength() const;
  CCTK_REAL business() const;
};

// A collection of workers
template <typename item_t> class workers_t {
  typedef vector<worker_t<item_t> > coll_t;
  coll_t workers;
  typename coll_t::iterator find_least_busy_worker();
  typename coll_t::iterator find_most_busy_worker();

public:
  workers_t(int nworkers);
  bool empty() const;
  size_t size() const;
  CCTK_REAL load() const;
  CCTK_REAL strength() const;
  CCTK_REAL ideal_business() const;
  CCTK_REAL imbalance() const;
  worker_t<item_t> &get_least_busy_worker();
  worker_t<item_t> &get_most_busy_worker();
  void ensure_same_size();
  void copy_out(vector<vector<item_t> > &dst) const;
};

template <typename item_t> items_t<item_t>::items_t() {}

template <typename item_t>
items_t<item_t>::items_t(vector<item_t> const &items_) {
  for (typename vector<item_t>::const_iterator p = items_.begin();
       p != items_.end(); ++p) {
    add_item(*p);
  }
}

template <typename item_t> void items_t<item_t>::add_item(item_t const &item) {
  items.push_back(item);
}

template <typename item_t> bool items_t<item_t>::empty() const {
  return items.empty();
}

template <typename item_t> size_t items_t<item_t>::size() const {
  return items.size();
}

template <typename item_t> CCTK_REAL items_t<item_t>::load() const {
  CCTK_REAL total_load = 0.0;
  for (typename coll_t::const_iterator p = items.begin(); p != items.end();
       ++p) {
    total_load += item_load(*p);
  }
  return total_load;
}

template <typename item_t>
typename items_t<item_t>::coll_t::iterator
items_t<item_t>::find_largest_item() {
  typename coll_t::iterator max_item = items.end();
  CCTK_REAL max_load = -1.0;
  for (typename coll_t::iterator p = items.begin(); p != items.end(); ++p) {
    if (item_load(*p) > max_load) {
      max_item = p;
      max_load = item_load(*max_item);
    }
  }
  return max_item;
}

template <typename item_t> item_t &items_t<item_t>::get_one_item() {
  assert(not empty());
  return items.front();
}

template <typename item_t> item_t &items_t<item_t>::get_largest_item() {
  typename coll_t::iterator const max_item = find_largest_item();
  assert(max_item != items.end());
  return *max_item;
}

template <typename item_t>
item_t items_t<item_t>::get_and_remove_largest_item() {
  typename coll_t::iterator const max_item = find_largest_item();
  assert(max_item != items.end());
  item_t const item = *max_item;
  items.erase(max_item);
  return item;
}

template <typename item_t>
void items_t<item_t>::copy_out(vector<item_t> &dst) const {
  dst.resize(items.size());
  copy(items.begin(), items.end(), dst.begin());
}

template <typename item_t> CCTK_REAL worker_t<item_t>::strength() const {
  return 1.0; // All workers have the same strength
}

template <typename item_t> CCTK_REAL worker_t<item_t>::business() const {
  return this->load() / strength();
}

template <typename item_t>
workers_t<item_t>::workers_t(int const nworkers)
    : workers(nworkers) {}

template <typename item_t> bool workers_t<item_t>::empty() const {
  return workers.empty();
}

template <typename item_t> size_t workers_t<item_t>::size() const {
  return workers.size();
}

template <typename item_t> CCTK_REAL workers_t<item_t>::load() const {
  CCTK_REAL total_load = 0.0;
  for (typename coll_t::const_iterator w = workers.begin(); w != workers.end();
       ++w) {
    total_load += w->load();
  }
  return total_load;
}

template <typename item_t> CCTK_REAL workers_t<item_t>::strength() const {
  CCTK_REAL total_strength = 0.0;
  for (typename coll_t::const_iterator w = workers.begin(); w != workers.end();
       ++w) {
    total_strength += w->strength();
  }
  return total_strength;
}

template <typename item_t> CCTK_REAL workers_t<item_t>::ideal_business() const {
  return load() / strength();
}

template <typename item_t> CCTK_REAL workers_t<item_t>::imbalance() const {
  assert(not empty());
  CCTK_REAL max_load = 0.0;
  CCTK_REAL avg_load = 0.0;
  for (typename coll_t::const_iterator w = workers.begin(); w != workers.end();
       ++w) {
    max_load = max(max_load, w->business());
    avg_load += w->business();
  }
  avg_load /= size();
  return max_load - avg_load;
}

template <typename item_t>
typename workers_t<item_t>::coll_t::iterator
workers_t<item_t>::find_least_busy_worker() {
  typename coll_t::iterator min_worker = workers.end();
  CCTK_REAL min_business = numeric_limits<CCTK_REAL>::max();
  for (typename coll_t::iterator w = workers.begin(); w != workers.end(); ++w) {
    if (w->business() < min_business) {
      min_worker = w;
      min_business = min_worker->business();
    }
  }
  return min_worker;
}

template <typename item_t>
typename workers_t<item_t>::coll_t::iterator
workers_t<item_t>::find_most_busy_worker() {
  typename coll_t::iterator max_worker = workers.end();
  CCTK_REAL max_business = 0.0;
  for (typename coll_t::iterator w = workers.begin(); w != workers.end(); ++w) {
    if (w->business() > max_business) {
      max_worker = w;
      max_business = max_worker->business();
    }
  }
  return max_worker;
}

template <typename item_t>
worker_t<item_t> &workers_t<item_t>::get_least_busy_worker() {
  typename coll_t::iterator const min_worker = find_least_busy_worker();
  assert(min_worker != workers.end());
  return *min_worker;
}

template <typename item_t>
worker_t<item_t> &workers_t<item_t>::get_most_busy_worker() {
  typename coll_t::iterator const max_worker = find_most_busy_worker();
  assert(max_worker != workers.end());
  return *max_worker;
}

template <typename item_t> void workers_t<item_t>::ensure_same_size() {
  if (empty())
    return; // nothing to do

  size_t max_items = 0;
  typename coll_t::iterator nonempty_worker = workers.end();
  for (typename coll_t::iterator w = workers.begin(); w != workers.end(); ++w) {
    if (nonempty_worker == workers.end() and not w->empty()) {
      nonempty_worker = w;
    }
    max_items = max(max_items, w->size());
  }
  if (max_items == 0)
    return; // load is already equal

  for (typename coll_t::iterator w = workers.begin(); w != workers.end(); ++w) {
    // find a worker who has an item
    typename coll_t::iterator const worker = w->empty() ? nonempty_worker : w;
    while (w->size() < max_items) {
      CCTK_REAL const ratio = 0.0; // create empty fill-up items
      w->add_item(item_split(worker->get_one_item(), ratio));
    }
  }

  for (typename coll_t::const_iterator w = workers.begin(); w != workers.end();
       ++w) {
    assert(w->size() == max_items);
  }
}

template <typename item_t>
void workers_t<item_t>::copy_out(vector<vector<item_t> > &dst) const {
  dst.resize(workers.size());
  for (size_t w = 0; w < workers.size(); ++w) {
    workers.at(w).copy_out(dst.at(w));
  }
}

//////////////////////////////////////////////////////////////////////////////

template <typename item_t>
void assign_item(items_t<item_t> &items, workers_t<item_t> &workers) {
  // Assign the largest item to the least busy worker
  item_t const item = items.get_and_remove_largest_item();
  worker_t<item_t> &worker = workers.get_least_busy_worker();
  worker.add_item(item);
}

template <typename item_t>
void split_and_distribute(workers_t<item_t> &workers) {
  // Split the largest item of the most busy worker and give the
  // remainder to another worker
  worker_t<item_t> &worker = workers.get_most_busy_worker();
  item_t &item = worker.get_largest_item();

  // Determine how to split the item
  CCTK_REAL const imbalance = worker.business() - workers.ideal_business();
  // Should we even be here?
  assert(imbalance > 0.0);
  CCTK_REAL const item_business = item_load(item) / worker.strength();
  // This should be the largest item!
  assert(item_business > 0.0);
  // Determine how much of the item to give away
  CCTK_REAL const ratio = imbalance / item_business;
  // A ratio of one or more indicates that the item should be given
  // away in its entirety -- which would mean that something went
  // wrong before
  assert(ratio < 1.0);

  // Split the item...
  item_t const new_item = item_split(item, ratio);
  // ...and give the remainder to someone else
  worker_t<item_t> &new_worker = workers.get_least_busy_worker();
  // This should be someone else!
  assert(&new_worker != &worker);
  new_worker.add_item(new_item);
}

template <typename item_t>
void balance(vector<item_t> const &items_, vector<vector<item_t> > &split_items,
             int const nworkers, CCTK_REAL const max_imbalance,
             bool const ensure_same_size) {
  assert(max_imbalance >= 0.0);

  if (items_.empty())
    return; // nothing to do
  items_t<item_t> items(items_);
  assert(nworkers > 0);
  workers_t<item_t> workers(nworkers);

  // Assign items
  while (not items.empty()) {
    assign_item(items, workers);
  }

  // TODO: To parallelise this: Group workers into groups, and
  // assign work to these groups.  Then balance the load
  // recursively.

  // Balance items
  for (;;) {
    // Measure imbalance
    CCTK_REAL const imbalance = workers.imbalance();
    // Are we done?
    if (imbalance <= max_imbalance)
      break;
    // Should we even be here?
    assert(workers.size() > 1);

    // Do something
    split_and_distribute(workers);

    // Ensure progress
    assert(workers.imbalance() < imbalance);
  }

  if (ensure_same_size) {
    workers.ensure_same_size();
  }

  workers.copy_out(split_items);
}

template void balance(vector<region_t> const &items_,
                      vector<vector<region_t> > &split_items,
                      int const nworkers, CCTK_REAL const max_imbalance,
                      bool const ensure_same_size);

} // namespace CarpetLib
