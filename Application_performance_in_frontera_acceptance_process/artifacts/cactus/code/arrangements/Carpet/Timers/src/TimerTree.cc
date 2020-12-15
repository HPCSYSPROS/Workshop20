/*
 *
 * The MIT License
 *
 * Copyright (c) 1997-2010 Center for the Simulation of Accidental Fires and
 * Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
 * University of Utah.
 *
 * License for the specific language governing rights and limitations under
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <utility>

#include <dist.hh>
#include <mpi_string.hh>

#include "TimerTree.hh"

namespace Timers {

using namespace std;

TimerNode::TimerNode(TimerTree *tree, const string &name)
    : d_name(name), d_parent(0), d_tree(tree), d_running(false), d_timer(0) {}

TimerNode::~TimerNode() {
  for (map<string, TimerNode *>::iterator iter = d_children.begin();
       iter != d_children.end(); ++iter) {
    delete iter->second;
  }
  delete d_timer;
}

string TimerNode::pathName() const {
  assert(d_parent != this);
  if (d_parent)
    return d_parent->pathName() + string("/") + getName();
  else
    return getName();
}

void TimerNode::instantiate() {
  assert(!d_running);
  if (d_parent)
    assert(d_parent == d_tree->current);
  else
    d_parent = d_tree->current;
  d_tree->current = this;
  if (!d_timer)
    d_timer = new CactusTimer(pathName());
  d_tree->current = d_parent;
}

void TimerNode::start() {
  assert(!d_running);

  d_running = true;
  if (d_parent)
    assert(d_parent == d_tree->current);
  else
    d_parent = d_tree->current;
  d_tree->current = this;
  if (!d_timer)
    d_timer = new CactusTimer(pathName());
  assert(d_timer);
  d_timer->start();
}

void TimerNode::stop() {
  assert(d_running);

  // A timer can only be stopped if it is the current timer
  if (this != d_tree->current)
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Tried to stop non-current timer '%s'", getName().c_str());

  d_timer->stop();

  d_running = false;
  d_tree->current = d_parent;
}

/// Get the name of the timer
string TimerNode::getName() const {
  assert(not d_name.empty());
  return d_name;
}

/// Determine if the timer is running
bool TimerNode::isRunning() const { return d_running; }

/// Find the child timer that matches the name provided.  If it is
/// not found then a new timer with that name is allocated.
TimerNode *TimerNode::getChildTimer(const string &name) {
  // Find child
  TimerNode *child = d_children[name];

  // If the pointer is null then allocate it
  if (not child)
    d_children[name] = child = new TimerNode(d_tree, name);

  return child;
}

/// Get the time measured by this timer
double TimerNode::getTime() { return d_timer->getTime(); }

/// Get the global time measured by this timer
void TimerNode::getGlobalTime(double &avg, double &max) {
  return d_timer->getGlobalTime(avg, max);
}

/// Get the names of all clocks of this timer
vector<pair<string, string> > TimerNode::getAllTimerNames() const {
  return d_timer->getAllTimerNames();
}

/// Get the values of all clocks of this timer
vector<double> TimerNode::getAllTimerValues() {
  return d_timer->getAllTimerValues();
}

/// Print this node and its children as an ASCII tree
void TimerNode::print(ostream &out, double total, int level, double threshold,
                      int precision) {
  string space;

  // Compute the level of indentation for this depth
  for (int i = 0; i < level - 1; i++)
    space += "| ";

  if (level != 0)
    space += "|_";

  const int pcw = 6;
  const int tw = 8;
  const int tnw = 50; // timer name
  const int vw = 9;   // clock values
  const streamsize oldprecision = out.precision();
  const ios_base::fmtflags oldflags = out.flags();

  // const double t = getTime();
  double tavg, tmax;
  getGlobalTime(tavg, tmax);
  const vector<double> values = getAllTimerValues();
  const string hyphens = string(precision - 1, '-');
  const string spaces = string(precision - 1, ' ');

  if (level == 0) {
    const vector<pair<string, string> > names = getAllTimerNames();

    out << "--------" << hyphens << "--------" << hyphens << "--------"
        << hyphens << "--" << string(tnw, '-');
    for (size_t i = 0; i < values.size(); ++i) {
      out << "--" << string(vw, '-');
    }
    out << "\n";

    // timer names
    out << "Time    " << spaces << "  Time  " << spaces << " Imblnc " << spaces
        << "  " << setw(tnw) << left << "Timer" << right;
    for (size_t i = 0; i < names.size(); ++i) {
      out << "  " << setw(vw) << names[i].first.substr(0, vw);
    }
    out << "\n";

    // timer units
    out << "percent " << spaces << "  secs  " << spaces << " percent" << spaces
        << "  " << setw(tnw) << "     ";
    for (size_t i = 0; i < names.size(); ++i) {
      out << "  " << setw(vw) << names[i].second.substr(0, vw);
    }
    out << "\n";

    out << "--------" << hyphens << "--------" << hyphens << "--------"
        << hyphens << "--" << string(tnw, '-');
    for (size_t i = 0; i < values.size(); ++i) {
      out << "--" << string(vw, '-');
    }
    out << "\n";
  }

  // Print this timer value
  out << fixed << setw(pcw) << setprecision(precision) << 100.0 * tavg / total
      << "%"
      << " " << fixed << setw(tw) << setprecision(precision) << tavg << " "
      << fixed << setw(pcw) << setprecision(precision)
      << 100.0 * (1.0 - tavg / tmax) << "%"
      << "  " << space << setw(max(size_t(0), tnw - space.length())) << left
      << d_name.substr(0, max(size_t(10), tnw - space.length())) << right;
  for (size_t i = 0; i < values.size(); ++i) {
    out.unsetf(ios_base::floatfield);
    out << "  " << setw(vw) << setprecision(vw - 5) << values[i];
  }
  out << "\n";

  // TODO: Don't call getGlobalTime for all timers separately.
  // Instead, call a single function that takes a snapshot of all
  // timers, and reduces these snapshots to average and maximum.
  // While printing, access these snapshots.

  // double children_time = 0;
  double children_tavg = 0.0;
  bool printed_children = false;

  // Recursively print the children
  for (map<string, TimerNode *>::iterator iter = d_children.begin();
       iter != d_children.end(); iter++) {
    const string timername = iter->first;
    const string root_timername =
        CarpetLib::broadcast_string(dist::comm(), 0, timername);
    if (timername != root_timername) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Timers are inconsistent across processes: root process "
                  "expects timer %s, this process has timer %s instead",
                  root_timername.c_str(), timername.c_str());
    }
    double child_avg, child_max;
    iter->second->getGlobalTime(child_avg, child_max);
    if (child_max * 100.0 / total > threshold) {
      iter->second->print(out, total, level + 1, threshold, precision);
      printed_children = true;
    }
    // children_time += iter->second->getTime();
    children_tavg += child_avg;
  }
  {
    const string timername = "[done]";
    const string root_timername =
        CarpetLib::broadcast_string(dist::comm(), 0, timername);
    if (timername != root_timername) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Timers are inconsistent across processes: root process "
                  "expects timer %s, this process has timer %s instead",
                  root_timername.c_str(), timername.c_str());
    }
  }

  if (d_children.size() > 0 && printed_children) {
    // const double untimed = t - children_time;
    const double untimed = tavg - children_tavg;

    if (100.0 * untimed / total > threshold) {
      // Print the untimed portion
      out << fixed << setw(pcw) << setprecision(1) << 100.0 * untimed / total
          << "%"
          << " " << fixed << setw(tw) << setprecision(1) << untimed
          << "        "
          << "  | " << space << "untimed"
          << "\n";
    }
  }
  out.precision(oldprecision);
  out.setf(oldflags);

  if (level == 0) {
    out << "--------" << hyphens << "--------" << hyphens << "--------"
        << hyphens << "--" << string(tnw, '-');
    for (size_t i = 0; i < values.size(); ++i) {
      out << "--" << string(vw, '-');
    }
    out << "\n";
  }
}

void TimerNode::outputXML(const string &out_dir, int proc) {
  ostringstream filenamebuf;
  filenamebuf << out_dir << "/timertree." << proc << ".xml";
  string filenamestr = filenamebuf.str();
  const char *filename = filenamestr.c_str();
  ofstream file;
  file.open(filename, ios::out | ios::trunc);

  printXML(file, 0);

  file.close();
  assert(file.good());
}

/// Print this node and its children as an XML file
void TimerNode::printXML(ostream &out, int level) {
  string space;

  // Compute the level of indentation for this node
  for (int i = 0; i < level; i++)
    space += "  ";

  out << space << "<timer name = "
      << "\"" << escapeForXML(d_name) << "\"> ";
  out << getTime() << " ";
  const vector<double> values = getAllTimerValues();
  const vector<pair<string, string> > names = getAllTimerNames();

  out << endl;

  for (size_t i = 0; i < names.size(); ++i) {
    out << space << "  "
        << "<value clock=\"" << names[i].first << "\" unit=\""
        << names[i].second << "\">" << values[i] << "</value>" << endl;
  }

  // For compactness, only use multiple lines if there are children
  if (d_children.size() != 0) {
    out << "\n";

    // Recursively print the children
    for (map<string, TimerNode *>::iterator iter = d_children.begin();
         iter != d_children.end(); ++iter)
      iter->second->printXML(out, level + 1);
    out << space;
  }

  out << "</timer>"
      << "\n";
}

/// Make a string suitable for inclusion in an XML file
string TimerNode::escapeForXML(const string &s) const {
  ostringstream res;
  for (string::const_iterator si = s.begin(); si != s.end(); ++si) {
    switch (*si) {
    case '\'':
      res << "&apos;";
      break;
    case '"':
      res << "&quot;";
      break;
    case '&':
      res << "&amp;";
      break;
    case '<':
      res << "&lt;";
      break;
    case '>':
      res << "&gt;";
      break;
    default:
      res << *si;
      break;
    }
  }

  return res.str();
}

} // namespace Timers
