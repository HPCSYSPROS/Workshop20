#include <stdio.h>
#include <stdlib.h>
// #include <bgl_perfctr.h>
//#include <mpi.h>
#include <bglpersonality.h>
#include <rts.h>

extern "C" int init_personality();
extern "C" int get_rank_pset();
extern "C" int get_pset_num();
extern "C" int numnodesinpset();

extern "C" int get_xcoord();
extern "C" int get_ycoord();
extern "C" int get_zcoord();
extern "C" int get_xsize();
extern "C" int get_ysize();
extern "C" int get_zsize();

BGLPersonality p;

int get_xsize()
{
  int r;
  r = (int) p.getXsize();
  //printf("x Size= %d\n",r);
  return r;
}
int get_ysize()
{
  int r;
  r = (int) p.getYsize();
  //printf("y zise= %d\n",r);
  return r;
}
int get_zsize()
{
  int r;
  r = (int) p.getZsize();
  //printf("z size= %d\n",r);
  return r;
}

int get_xcoord()
  //int GET_XCOORD()
{
  int r;
  //r = (int) p.xCoord;
  r = (int) p.getXcoord();
  //printf("x Coord= %d\n",r);
  return r;
}

int get_ycoord()
  //int GET_YCOORD()
{
  int r;
  //r = (int) p.yCoord;
  r = (int) p.getYcoord();
  //printf("y Coord= %d\n",r);
  return r;
}

int get_zcoord()
  //int GET_ZCOORD()
{
  int r;
  //r = (int) p.zCoord;
  r = (int) p.getZcoord();
  //printf("z Coord= %d\n",r);
  return r;
}

int init_personality()
{
  int err;
  err = rts_get_personality(&p,sizeof(p));
  if( err > 0) printf("Error initializing personality\n"); 
}

int get_rank_pset()
  //int GET_RANK_PSET()
{
  int r;

  r = p.rankInPset();
  //  printf("Rank %d in pset\n",r);
  return r;
}

int get_pset_num()
  //int GET_PSET_NUM()
{
  int r;

  r = p.getPsetNum();
  //printf("Pset num. %d\n",r);
  return r;
}

int numnodesinpset()
  //int NUMNODESINPSET()
{
  int r;

  r = p.numNodesInPset();
  // printf("Nodes in Pset: %d\n",r);
  return r;
}

  


