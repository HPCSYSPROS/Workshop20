/* Test bed for Mortin ordering */

#include <stdio.h>
#include <stdlib.h>

static int squaresize[4];	   /* dimensions of sublattice for each MPI rank */
static int nsquares[4];	           /* number of hypercubes in each direction */
static int *m2hc;                  /* hash converting hypercube index to morton index */
static int *hc2m;                  /* hash converting morton index to hypercube index */

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )
# define MAXFACTOR 16

/* Would be in lattice.h */
static int nx, ny, nz, nt;
int sites_on_node;

/* External defines */
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3

/* Normally external functions */
static int ranks;
int numnodes(){
  return ranks;
}

void terminate(int status){
  exit(status);
}

/*------------------------------------------------------------------*/
/* Convert rank to coordinates */
static void lex_coords(int coords[], const int dim, const int size[], 
	   const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Parity of the coordinate */
static int coord_parity(int r[]){
  return (r[0] + r[1] + r[2] + r[3]) % 2;
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Map node number and hypercube index to coordinates  */
/* (The inverse of node_number and node_hc_index) */
/* Assumes even sites come first */
static void get_coords_hc(int coords[], int node, int hc_index){
  int mc[4];
  int ir;
  int meo, neven, xeo;
  int k = node;

  /* mc = the machine coordinates for node k */
  lex_coords(mc, 4, nsquares, k);

  /* meo = the parity of the machine coordinate */
  meo = coord_parity(mc);

  /* neven = the number of even sites on node k */
  neven = (sites_on_node + 1 - meo)/2;
  
  /* ir = the even part of the lexicographic index within the
     sublattice on node k */
  if(hc_index < neven){
    ir = 2*hc_index;
    xeo = 0;
  } else {
    ir = 2*(hc_index - neven);
    xeo = 1;
  }

  /* coords = the sublattice coordinate */
  lex_coords(coords, 4, squaresize, ir);

  /* Add offset to get full lattice coordinate (still a 2-fold ambiguity) */
  coords[XUP] += mc[XUP]*squaresize[XUP];
  coords[YUP] += mc[YUP]*squaresize[YUP];
  coords[ZUP] += mc[ZUP]*squaresize[ZUP];
  coords[TUP] += mc[TUP]*squaresize[TUP];

  /* Adjust coordinate according to parity */
  if( coord_parity(coords) != xeo ){
    coords[XUP]++;
    if(coords[XUP] >= squaresize[XUP]*(mc[XUP]+1)){
      coords[XUP] -= squaresize[XUP];
      coords[YUP]++;
      if(coords[YUP] >= squaresize[YUP]*(mc[YUP]+1)){
	coords[YUP] -= squaresize[YUP];
	coords[ZUP]++;
	if(coords[ZUP] >= squaresize[ZUP]*(mc[ZUP]+1)){
	  coords[ZUP] -= squaresize[ZUP];
	  coords[TUP]++;
	}
      }
    }
  }

}

/*------------------------------------------------------------------*/
/* Map coordinate to the node number (the MPI rank that has the coordinate */
int node_number(int x, int y, int z, int t) {
register int i;
    x /= squaresize[XUP]; y /= squaresize[YUP];
    z /= squaresize[ZUP]; t /= squaresize[TUP];
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
    return( i );
}

/*------------------------------------------------------------------*/
/* Map a global lattice coordinate to the coordinate on this node */
static void local_coords(int *xr, int *yr, int *zr, int *tr, int x, int y, int z, int t){
    *xr = x%squaresize[XUP]; *yr = y%squaresize[YUP];
    *zr = z%squaresize[ZUP]; *tr = t%squaresize[TUP];
}
/*------------------------------------------------------------------*/
/* Map a global lattice coordinate to the hypercubic node index */
int node_hc_index(int x, int y, int z, int t) {
  int i,xr,yr,zr,tr;
  local_coords(&xr, &yr, &zr, &tr, x, y, z, t);
  xr = x%squaresize[XUP]; yr = y%squaresize[YUP];
  zr = z%squaresize[ZUP]; tr = t%squaresize[TUP];
  i = xr + squaresize[XUP]*( yr + squaresize[YUP]*( zr + squaresize[ZUP]*tr));
  if( (x+y+z+t)%2==0 ){	/* even site */
    return( i/2 );
  }
  else {
    return( (i + sites_on_node)/2 );
  }
}

/*------------------------------------------------------------------*/
/* Map a coordinate to the node index (the rank of the site on this node) */
int node_index(int x, int y, int z, int t) {
  return hc2m[node_hc_index(x, y, z, t)];
}

/*------------------------------------------------------------------*/
/* Map node number and Morton index to coordinates  */
/* (The inverse of node_number and node_hc_index) */
/* Assumes even sites come first */
void get_coords(int coords[], int node, int m_index){
  int k;

  get_coords_hc(coords, node, m2hc[m_index]);

  /* Consistency checks for debugging */
  if((k = node_number(coords[0], coords[1], coords[2], coords[3])) 
     != node){
    printf("get_coords: coords %d %d %d %d for node %d index %d map to wrong node %d\n",
	   coords[0], coords[1], coords[2], coords[3], node, m_index, k);
    terminate(1);
  }
  if((k = node_index(coords[0], coords[1], coords[2], coords[3]))
      != m_index){
    printf("get_coords: coords %d %d %d %d for node %d index %d map to wrong index %d\n",
	   coords[0], coords[1], coords[2], coords[3], node, m_index, k);
    terminate(1);
  }
}

/*--------------------------------------------------------------------*/
/* Find the least prime factor of n from the prime list */

static int least_prime(int n){

  int k;
  for(k = 0; k < MAXPRIMES; k++){
    if(n % prime[k] == 0) return prime[k];
  }
  printf("Can't find prime factor in %d\n", n);
  terminate(1);
  return 0;
}

/*--------------------------------------------------------------------*/
/* Print the prime factors                                            */

static void print_factors(int factors[], int nfactors){

  int i;
  for(i = 0; i < nfactors; i++)
    printf("%d ",factors[i]);
  printf("\n");
}

/*--------------------------------------------------------------------*/
/* Get a list of all prime factors of n                              */

static void prime_factors(int n, int factors[], int *nfactors, int max){

  int i;
  int m = n;

  if(m <= 1){
    printf("prime_factors: improper factor request: %d\n", m);
    terminate(1);
  }

  for(i = 0; i < max; i++){
    if(m == 1){
      *nfactors = i;
      break;
    }
    factors[i] = least_prime(m);
    m = m/factors[i];
  }
  if(m != 1){
    printf("Can't factor %d: got ", n);
    print_factors(factors, *nfactors);
    terminate(1);
  }
}

/*--------------------------------------------------------------------------*/
/* Convert the Cartesian coordinate to the Morton basis expressed in digits */

static void decimal2morton(int digits[], int factors[], int nfactors, int coord){
  int k;
  int c = coord;

  for(k = 0; k < nfactors; k++){
    digits[k] = c % factors[k];
    c /= factors[k];
  }
}

/*--------------------------------------------------------------------*/
/* Convert the coordinate in the Morton basis to the Cartesian basis */

static void morton2decimal(int *coord, int digits[], int factors[], int nfactors){
  int k;
  int c = digits[nfactors-1];

  for(k = nfactors-1; k >= 1; k--)
    c = c*factors[k] + digits[k-1];

  *coord = c;
}
/*--------------------------------------------------------------------*/
/* Count up to the specified limit                                    */

static void morton_hash(int m2hc[], int hc2m[], 
			int factors[4][MAXFACTOR], int nfactors[], 
			int lastfactor, int limit){

  int i, k, dir, parity;
  int digits[4][MAXFACTOR];
  int coords[4];
  int hypercube_index, morton_index;

  for(dir = 0; dir < 4; dir++)
    for(k = 0; k < nfactors[dir]; k++)
      digits[dir][k] = 0;

  morton_index = 0;
  for(parity = 0; parity < 2; parity++){
    for(i = 0; i < sites_on_node; i++){
      for(dir = 0; dir < 4; dir++){
	morton2decimal(&coords[dir], digits[dir], factors[dir], nfactors[dir]);
	hypercube_index = node_hc_index(coords[XUP], coords[YUP], coords[ZUP], coords[TUP]);
      }
      if(coord_parity(coords) == parity){
	m2hc[hypercube_index] = morton_index;
	hc2m[morton_index] = hypercube_index;
	printf("Morton %d hypercube %d\n", morton_index, hypercube_index);
	for(dir = 0; dir < 4; dir++){
	  printf("%d: ",coords[dir]);
	  print_factors(digits[dir], nfactors[dir]);
	}
	morton_index++;
      }
      for(k = 0; k < lastfactor; k++){
	for(dir = 0; dir < 4; dir++){
	  if(k < nfactors[dir]){
	    digits[dir][k]++;
	    if(digits[dir][k] < factors[dir][k])break;
	    digits[dir][k] = 0;
	  }
	}
	if(dir < 4)break;
      }
    }
  }
}

/*--------------------------------------------------------------------*/
/* Determine which sites reside on which nodes (MPI ranks)            */
/* Use a hypercubic division                                          */
static void setup_hyper_prime(){
  int i,j,k,dir,p;

  /* Figure out dimensions of rectangle */
  squaresize[XUP] = nx; squaresize[YUP] = ny;
  squaresize[ZUP] = nz; squaresize[TUP] = nt;
  nsquares[XUP] = nsquares[YUP] = nsquares[ZUP] = nsquares[TUP] = 1;
  
  i = 1;	/* current number of hypercubes */
  while(i<numnodes()){
    /* figure out which prime to divide by starting with largest */
    while( (numnodes()/i)%prime[k] != 0 && k>0 ) --k;
    p = prime[k];

    /* figure out which direction to divide */
    
    /* find largest dimension of h-cubes divisible by p */
    for(j=0,dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]>j && squaresize[dir]%p == 0 )
	j=squaresize[dir];
    
    /* if one direction with largest dimension has already been
       divided, divide it again.  Otherwise divide first direction
       with largest dimension. */
    for(dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]==j && nsquares[dir]>1 )break;
    if( dir > TUP)for(dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]==j )break;
    /* This can fail if I run out of prime factors in the dimensions */
    if(dir > TUP){
      printf("LAYOUT: Can't lay out this lattice, not enough factors of %d\n" ,p);
      terminate(1);
    }
    
    /* do the surgery */
    i*=p; squaresize[dir] /= prime[k]; nsquares[dir] *= prime[k];
  }
}

/*------------------------------------------------------------------*/
/* Set up the Morton order layout.  This is defined in a lazy way
   by setting up a hash table for converting the hypercubic order
   to Morton order */
static void setup_morton(void){
  int dir;
  int factors[4][MAXFACTOR];
  int nfactors[4];
  int lastfactor;

  /* Determine the Morton basis for local coordinates */
  for(dir = 0; dir < 4; dir++){
    prime_factors(squaresize[dir], factors[dir], &nfactors[dir], MAXFACTOR);
    printf("%d: ", squaresize[dir]);
    print_factors(factors[dir], nfactors[dir]);
  }

  /* The maximum number of Morton digits */
  lastfactor = nfactors[0];
  for(dir = 1; dir < 4; dir++){
    if(lastfactor < nfactors[dir])
      lastfactor = nfactors[dir];
  }

  /* Set up hash tables */
  hc2m = (int *)malloc(sizeof(int)*sites_on_node);
  m2hc = (int *)malloc(sizeof(int)*sites_on_node);

  morton_hash(m2hc, hc2m, factors, nfactors, lastfactor, sites_on_node);
}

/*------------------------------------------------------------------*/
/* Initialization entry point.  Determine the distribution of sites
   across the nodes and the order of sites on each node            */

void setup_layout(void){


  /* Set the sublattice dimensions according to the hyper_prime algorithm */
  setup_hyper_prime();
  printf("automatic hyper_prime layout\n");

  /* Number of sites on node */
  sites_on_node =
    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];

  printf("ON EACH NODE %d x %d x %d x %d\n",squaresize[XUP],squaresize[YUP],
	 squaresize[ZUP],squaresize[TUP]);

  /* Set up the Morton order */
  setup_morton();

  /* The Morton-order hash tables must be the same for all nodes,
     so we must require an even number of sites per node.  Otherwise
     number of even sites on a node varies with the node number and
     therefore the hash table varies.  We could fix this by maintaining
     two hash tables, one for nodes with even machine coordinates and
     one for nodes with odd machine coordinates. */
  if(sites_on_node % 2 != 0){
    printf("ODD SITES ON NODE NOT SUPPORTED WITH THIS MORTON ORDERING\n");
    terminate(1);
  }
}

/*------------------------------------------------------------------*/
/* Test bed for the layout routines                                 */

int main(int argc, char *argv[]){
  int flag;
  int coords[4];
  int x,y,z,t;
  int node, morton_index;
  int status;

  printf("Enter the lattice dimensions nx, ny, nz, nt ");
  scanf("%d%d%d%d", &nx, &ny, &nz, &nt);

  printf("Enter the number of nodes (MPI ranks) ");
  scanf("%d", &ranks);

  /* Set up the lattice layout */
  setup_layout();

  while(1){
    printf("Enter 0 for mapping coord->index; 1 for mapping index->coord ");
    status = scanf("%d", &flag);
    if(status != 1)return 0;

    if(flag == 0){
      while(1){
	printf("Enter coordinates x, y, z, t ");
	status = scanf("%d%d%d%d", &x, &y, &z, &t);
	if(status != 4)break;
	printf("node %d index %d\n", node_number(x, y, z, t),
	       node_index(x, y, z, t));
      }
      printf("\n");
    } else {
      while(1){
	printf("Enter node number and Morton node index ");
	status = scanf("%d%d", &node, &morton_index);
	if(status != 2)break;
	get_coords(coords, node, morton_index);
	printf("coords %d %d %d %d\n", coords[0], coords[1], coords[2], coords[3]);
      }
    }
  }

  return 0;
}
