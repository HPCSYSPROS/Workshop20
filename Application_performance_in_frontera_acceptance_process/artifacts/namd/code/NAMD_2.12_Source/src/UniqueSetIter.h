/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   UniqueSet template - (hashtable)
*/

#ifndef USITER_H
#define USITER_H

#include "InfoStream.h"
#include "UniqueSetRaw.h"
#include "UniqueSet.h"


template <class T> class UniqueSetIter {

  private:

    UniqueSet<T> *us;
    EntryGlob<T> *eg;
    int index;
  
  public:
  
    UniqueSetIter(void) { us = NULL; eg = NULL; index = 0; }

    UniqueSetIter(UniqueSet<T>& us_param) { 
       us = &us_param; eg = us_param.rep->globHead; index = 0;
    }

    UniqueSetIter(const UniqueSetIter<T>& iter) {
       us = iter.us; eg = iter.eg; index = iter.index;
    }

    UniqueSetIter<T>& operator= (const UniqueSetIter<T>& iter) {
       us = iter.us; eg = iter.eg; index = iter.index; return (*this);
    }

    ~UniqueSetIter(void) {}
  
    T *operator->(void) { 
      gotoUsed();
      if (eg)
        return (T *)&(eg->glob[index].obj);
      else { 
        index = 0;
        return(NULL);
      }
    }

    UniqueSetIter<T> begin(void) const {
        UniqueSetIter<T> iter;
        iter.us = us;
        iter.index = 0;
        iter.eg = us->rep->globHead;
        iter.gotoUsed();
        return(iter);
    }

    UniqueSetIter<T> end(void) const {
        UniqueSetIter<T> iter;
        iter.us = us;
        iter.index = 0;
        iter.eg = NULL;
        return(iter);
    }
        
    int operator!= (const UniqueSetIter<T> &iter) const {
        return (iter.index != index || iter.eg != eg);
    }

    int operator== (const UniqueSetIter<T> &iter) const {
        return (!operator!=(iter));
    }

    UniqueSetIter<T> operator++(void) {
      index++;
      gotoUsed();
      return (*this);
    }

    UniqueSetIter<T> operator++(int) {
       UniqueSetIter<T> tmp(*this);
       index++;
       gotoUsed();
       return (tmp);
    }
  
    T& operator* (void) {
       gotoUsed();
       return *((T *)&(eg->glob[index].obj));
    }
  
    void status(void) {
      std::cout << "Index is " << index << " addr is " << eg << std::endl;
    }
  
    void gotoUsed(void) {
      while (eg) {
        for(;index < us->rep->globSize; index++) {
	  if (eg->glob[index].isUsed()) break;
        }
        if (index < us->rep->globSize) break;
        index = 0;
        eg = eg->next();
      }
      if (!eg) index = 0;
    }
};

#endif
