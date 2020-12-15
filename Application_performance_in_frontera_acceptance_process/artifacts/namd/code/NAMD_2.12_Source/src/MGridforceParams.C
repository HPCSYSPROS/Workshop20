/*
 *  MGridforceParams.C
 *  
 *
 *  Created by Robert Brunner on 12/5/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 

#include "MGridforceParams.h"

MGridforceParams* MGridforceParamsList::find_key(const char* key)
{
    MGFElem* cur = head;
    MGFElem* found = NULL;
    MGridforceParams* result = NULL;
    
    while (found == NULL && cur != NULL) {
       if (!strcasecmp((cur->elem).gridforceKey,key)) {
        found = cur;
      } else {
        cur = cur->nxt;
      }
    }
    if (found != NULL) {
      result = &(found->elem);
    }
    return result;
}
  
int MGridforceParamsList::index_for_key(const char* key)
{
    MGFElem* cur = head;
    MGFElem* found = NULL;
    int result = -1;
    
    int idx = 0;
    while (found == NULL && cur != NULL) {
       if (!strcasecmp((cur->elem).gridforceKey,key)) {
        found = cur;
      } else {
        cur = cur->nxt;
	idx++;
      }
    }
    if (found != NULL) {
	result = idx;
    }
    return result;
}

MGridforceParams* MGridforceParamsList::at_index(int idx)
{
    MGFElem* cur = head;
    MGFElem* found = NULL;
    MGridforceParams* result = NULL;
    
    int counter = 0;
    while (found == NULL && cur != NULL) {
      if (counter == idx) {
	found = cur;
      } else {
        cur = cur->nxt;
	counter++;
      }
    }
    if (found != NULL) {
	result = &(found->elem);
    }
    return result;
}

 
MGridforceParams* MGridforceParamsList::add(const char* key) 
{
    // If the key is already in the list, we can't add it
    if (find_key(key)!=NULL) {
      return NULL;
    }
    
    MGFElem* new_elem = new MGFElem();
    int len = strlen(key);
    MGridforceParams* elem = &(new_elem->elem);
    elem->gridforceKey = new char[len+1];
    strncpy(elem->gridforceKey,key,len+1);
    elem->gridforceVfile = NULL;
    elem->gridforceFile = NULL;
    elem->gridforceCol = NULL;
    elem->gridforceQcol = NULL;
    elem->next = NULL;
    new_elem->nxt = NULL;
    if (head == NULL) {
      head = new_elem;
    }
    if (tail != NULL) {
      tail->nxt = new_elem;
      tail->elem.next = elem;
    }
    tail = new_elem;
    n_elements++;
    
    return elem;
}
  
void MGridforceParamsList::pack_data(MOStream *msg) 
{
    int i = n_elements;
    msg->put(n_elements);
    MGridforceParams *elem = get_first();
    while (elem != NULL) {
      int len;
      len = strlen(elem->gridforceKey) + 1;
      msg->put(len);
      msg->put(len,elem->gridforceKey);

      len = strlen(elem->gridforceVfile) + 1;
      msg->put(len);
      msg->put(len,elem->gridforceVfile);

      Vector v = elem->gridforceScale;
      msg->put(&v);
      
      len = strlen(elem->gridforceFile) + 1;
      msg->put(len);
      msg->put(len,elem->gridforceFile);
      
      len = strlen(elem->gridforceCol) + 1;
      msg->put(len);
      msg->put(len,elem->gridforceCol);
      
      if (elem->gridforceQcol == NULL) 
        msg->put(1); // Qcol_is_null = true
      else {
        msg->put(0); // Qcol_is_null = false
        len = strlen(elem->gridforceQcol) + 1;
        msg->put(len);
        msg->put(len,elem->gridforceQcol);
      }
      
      v = elem->gridforceVOffset;
      msg->put(&v);
      
      short boolvals[6];
      boolvals[0] = (elem->gridforceCont[0] ? 1 : 0);
      boolvals[1] = (elem->gridforceCont[1] ? 1 : 0);
      boolvals[2] = (elem->gridforceCont[2] ? 1 : 0);
      boolvals[3] = (elem->gridforceVolts ? 1 : 0);
      boolvals[4] = (elem->gridforceLite ? 1 : 0);
      boolvals[5] = (elem->gridforceCheckSize ? 1 : 0);
      msg->put(6,boolvals);
      
      i--;
      elem = elem->next;
    }
    if (i != 0) {
      NAMD_die("MGridforceParams message packing error\n");
    }
    return;
}
  
void MGridforceParamsList::unpack_data(MIStream *msg)
{
    int elements;
    msg->get(elements);
    
    for(int i=0; i < elements; i++) {
      // Get key
      int len;
      msg->get(len);
      char *key = new char[len];
      msg->get(len,key);
      MGridforceParams *elem = add(key);
      
      msg->get(len);
      char *str = new char[len];
      msg->get(len,str);
      elem->gridforceVfile = str;

      Vector v;
      msg->get(&v);
      elem->gridforceScale = v;
      
      msg->get(len);
      str = new char[len];
      msg->get(len,str);
      elem->gridforceFile = str;
      
      msg->get(len);
      str = new char[len];
      msg->get(len,str);
      elem->gridforceCol = str;
      
      int qcol_is_null;
      msg->get(qcol_is_null);
      if (qcol_is_null)
        elem->gridforceQcol = NULL;
      else {
        msg->get(len);
        str = new char[len];
        msg->get(len,str);
        elem->gridforceQcol = str;
      }
      
      msg->get(&v);
      elem->gridforceVOffset = v;
      
      short boolvals[6];
      msg->get(6,boolvals);
      elem->gridforceCont[0] = ( boolvals[0] != 0 );
      elem->gridforceCont[1] = ( boolvals[1] != 0 );
      elem->gridforceCont[2] = ( boolvals[2] != 0 );
      elem->gridforceVolts = ( boolvals[3] != 0 );
      elem->gridforceLite = ( boolvals[4] != 0 );
      elem->gridforceCheckSize = ( boolvals[5] != 0 );
      
      delete [] key;
    }
}
  
