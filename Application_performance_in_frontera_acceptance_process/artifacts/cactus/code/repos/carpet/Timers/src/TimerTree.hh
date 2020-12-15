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

/* This class was originally written by Justin Luitjens and
   subsequently integrated with Cactus/Carpet by Ian Hinder and
   heavily modified. */

#ifndef TIMERTREE_HH
#define TIMERTREE_HH

#include <cassert>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <utility>

#include "CactusTimer.hh"

namespace Timers {

class TimerNode;

class TimerTree {
public:
  TimerTree() : root(0), current(0) {}
  TimerNode *root;
  TimerNode *current;
};

/**
   The TimerNode class implements a tree structure where each node
   represents a timer, implemented as a CactusTimer. Each node of
   the tree can have zero of more children, where the names of the
   child nodes are unique within a single parent, but not
   necessarily unique within the entire tree. A tree formed of
   TimerNode objects represents the execution profile of the
   program, inasmuch as it is instrumented by timers.

   Child nodes of a given name are accessed using the getChildTimer
   method, where a node is created with the given name if none
   exists already. This ensures that the names of the child timers
   are unique.
*/

class TimerNode {
public:
  TimerNode(TimerTree *root, const std::string &name);
  ~TimerNode();

  void instantiate();
  void start();
  void stop();

  std::string getName() const;
  std::string pathName() const;

  // Find the child timer that matches the name provided. If it is
  // not found then that timer is allocated.
  TimerNode *getChildTimer(const std::string &name);

  double getTime();
  void getGlobalTime(double &avg, double &max);
  std::vector<std::pair<std::string, std::string> > getAllTimerNames() const;
  std::vector<double> getAllTimerValues();
  bool isRunning() const;

  void print(std::ostream &out, double total, int level = 0,
             double threshold = 0.0, int precision = 1);
  void printXML(std::ostream &out, int level = 0);
  void outputXML(const std::string &out_dir, int proc);

private:
  std::string escapeForXML(const std::string &s) const;

  std::string d_name;
  std::map<std::string, TimerNode *> d_children;
  TimerNode *d_parent;
  TimerTree *d_tree;
  bool d_running;
  CactusTimer *d_timer;
};

} // namespace Timers

#endif // TIMERTREE_HH
