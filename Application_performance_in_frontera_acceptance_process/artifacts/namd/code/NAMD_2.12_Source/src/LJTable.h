/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef LJTABLE_H
#define LJTABLE_H

#include "common.h"
#include "ProcessorPrivate.h"

class LJTable
{
public:
  struct TableEntry
  {
    // BigReal exclcut2;
#ifdef NAMD_KNL
    float A,B;
#else
    BigReal A;
    BigReal B;
#endif
  };

  LJTable(void);

  ~LJTable(void);

  const TableEntry *table_row(unsigned int i) const {
    return table + 2 * (i * table_dim);
  }

  const TableEntry *table_val(unsigned int i, unsigned int j) const {
    return table + 2 * (i * table_dim + j);
  }

  const TableEntry *table_val_scale14(unsigned int i, unsigned int j) const {
    return table + 2 * (i * table_dim + j) + 1;
  }

  const TableEntry *get_table() const { return table; }
  int get_table_dim() const { return table_dim; }

private:

  void compute_vdw_params(int i, int j, 
			  TableEntry *cur, TableEntry *cur_scaled);

  char *table_alloc;
  TableEntry *table;
  unsigned int table_dim;

};

#endif

