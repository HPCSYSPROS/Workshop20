/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Set.h"
#include "elements.h"

//#define DEBUG_IRSET

#ifdef DEBUG_IRSET
void NAMD_bug(const char *);
#endif

IRSet::IRSet() 
{
  head = (listNode *) 0;
  nElements = 0;
}

IRSet::~IRSet() {
  // delete all links; don't delete the InfoRecord objects
  listNode *tmp;
  for (listNode *link = head; link; link = tmp) {
    tmp = link->next;
    delete link;
  }
}

void IRSet::unchecked_insert(InfoRecord *info) 
{
#ifdef DEBUG_IRSET
  if (find(info)) NAMD_bug("IRSet::unchecked_insert duplicate");
#endif
    ++nElements;
    listNode *node = new listNode(info);
    node->next = head;
    head = node;
#ifdef DEBUG_IRSET
  int n = 0;
  while (node) { ++n; node = node->next; }
  if ( n != nElements ) NAMD_bug("IRSet::unchecked_insert count");
#endif
}


void IRSet::insert(InfoRecord *info) 
{
  if (!find(info))
  {
    ++nElements;
    listNode *node = new listNode(info);
    node->next = head;
    head = node;
#ifdef DEBUG_IRSET
    int n = 0;
    while (node) { ++n; node = node->next; }
    if ( n != nElements ) NAMD_bug("IRSet::insert count");
#endif
  }
   
}


void IRSet::myRemove(listNode **n, InfoRecord *r)
{
  if ((*n)->info == r)
    *n = (*n)->next;
  else 
    myRemove(&((*n)->next), r);
}

int IRSet::remove(InfoRecord * r) 
{
#ifdef DEBUG_IRSET
    listNode *node = head;
    int n = 0;
    while (node) { ++n; node = node->next; }
    if ( n != nElements ) NAMD_bug("IRSet::remove count");
#endif

  if (!head)
    return 0;

  listNode *p = head;
  listNode *q = p->next;

  if (p->info == r){
    head = q;
    delete p;
    --nElements;
    return 1;
  }

  while (q){
    if (q->info == r){
      p->next = q->next;
      delete q;
      --nElements;
      return 1;
    }
    else {
      p = q;
      q = q->next;
    }
  }
  return 0;
}

int IRSet::find(InfoRecord * r) 
{
  listNode *p = head;
  while (p) {
    if (p->info == r) return 1;
    else p = p->next;
  }
  return 0;
}

InfoRecord * IRSet::iterator(Iterator *iter)
{
  if (head){
    iter->next = head->next;
    return head->info;
  }
  return 0;
}

InfoRecord * IRSet::next(Iterator *iter)
{
  //  std::cout << "set::next: " << iter->next << "\n";
  if (!iter->next)
    { return 0;
    }
  //  std::cout << "set::next: iter->next->info=" << iter->next->info << "\n";
  InfoRecord *temp = iter->next->info;
  iter->next = iter->next->next;
  return temp;
}


int IRSet::numElements()
{
  return nElements;
}

int IRSet::hasElements()
{
  return ! ! head;
}

void IRSet::print() 
{
  listNode *p = head;
  while (p){
    if ( p->info ) iout << p->info->Id << " ";
    else iout << "NULL ";
    p = p->next;
  }
}
