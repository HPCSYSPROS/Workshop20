/******** layout_hyper_morton_prime.c *********/
/* MIMD version 7 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This is a version of layout_hyper_prime that uses the traditional
   hypercubic partitioning of the lattice to distribute sublattices
   on the nodes, but then Morton ordering for sites within the node.

   For creating sublsttices, this code divides the lattice by factors
   of prime numbers in any of the four directions.  It prefers to
   divide the longest dimensions, which mimimizes the area of the
   surfaces.  Similarly, it prefers to divide dimensions which have
   already been divided, thus not introducing more off-node
   directions.

	S. Gottlieb, May 18, 1999
	The code will start trying to divide with the largest prime factor
	and then work its way down to 2.  The current maximum prime is 53.
	The array of primes on line 46 may be extended if necessary.

   This requires that the lattice volume be divisible by the number
   of nodes.  Each dimension must be divisible by a suitable factor
   such that the product of the four factors is the number of nodes.

   3/29/00 EVENFIRST is the rule now. CD.
   12/10/00 Fixed so k = MAXPRIMES-1 DT
*/

// $Log: layout_hyper_prime.c,v $
// Revision 1.18  2012/11/24 04:43:47  detar
// Fix nsquares print out.
//
// Revision 1.17  2012/01/21 21:28:12  detar
// Support new QMP
//
// Revision 1.16  2011/11/29 20:11:30  detar
// Cosmetic fix to initialization
//
// Revision 1.15  2008/04/18 15:36:46  detar
// Permit odd number of lattice sites per node
//
// Revision 1.14  2008/04/11 15:36:00  detar
// Allow an odd number of sites per node
//

/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
     This routine sets the global variables "sites_on_node",
     "even_sites_on_node" and "odd_sites_on_node".
   num_sites(node) returns the number of sites on a node
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   get_logical_dimensions() returns the machine dimensions
   get_logical_coordinates() returns the mesh coordinates of this node
   get_coords() returns the coordinates for a given node and index
       (the inverse of node_number + node_index)
   io_node(node) maps nodes to their I/O node (for I/O partitions)
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
#ifdef HAVE_QMP
#include <qmp.h>
#endif

static int squaresize[4];	   /* dimensions of hypercubes */
static int nsquares[4];	           /* number of hypercubes in each direction */
static int machine_coordinates[4]; /* logical machine coordinates */ 
static int *m2hc;                  /* hash converting hypercube index to morton index */
static int *hc2m;                  /* hash converting morton index to hypercube index */


static int nodes_per_ionode[4];    /* dimensions of ionode partition */
static int *ionodegeomvals = NULL; /* ionode partitions */

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )
# define MAXFACTOR 16

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

#ifdef HAVE_QMP

/*--------------------------------------------------------------------*/
/* Sets the QMP logical topology if we need one */
static void set_qmp_logical_topology(const int *geom, int n){

  /* Has a geometry already been specified by the -geom command-line
     argument or on the input parameter line "node_geometry"? */
  /* If not, don't set the grid geometry here */
  if(geom == NULL)return;
  /* If so, then pass the grid dimensions to QMP now */
  if(QMP_declare_logical_topology(geom, n) != QMP_SUCCESS){
    node0_printf("setup_layout: QMP_declare_logical_topology failed on %d %d %d %d \n",
		 geom[0], geom[1], geom[2], geom[3] );
    terminate(1);
  }
}

/*--------------------------------------------------------------------*/
static void setup_qmp_grid(const int *nsquares2, int ndim2){
  int ndim = 4;
  int len[4] = {nx, ny, nz, nt};
  int i;

  if(mynode()==0){
    printf("qmp_grid,");
    printf("\n");
  }

  for(i=0; i<ndim; i++) {
    if(i<ndim2) nsquares[i] = nsquares2[i];
    else nsquares[i] = 1;
  }

  if(mynode()==0){
    printf("Using machine geometry: ");
    for(i=0; i<ndim; i++){
      printf("%d ",nsquares[i]);
      if(i < ndim-1)printf("X ");
    }
    printf("\n");
  }

  for(i=0; i<ndim; i++) {
    if(len[i]%nsquares[i] != 0) {
      node0_printf("LATTICE SIZE DOESN'T FIT GRID\n");
      QMP_abort(0);
    }
    squaresize[i] = len[i]/nsquares[i];
  }
}
#endif

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
/* Map a global lattice coordinate to the hypercubic node index */
int node_hc_index(int x, int y, int z, int t) {
register int i,xr,yr,zr,tr;
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
    c = c*factors[k-1] + digits[k-1];

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
    for(i = 0; i < limit; i++){
      FORALLUPDIR(dir){
	morton2decimal(&coords[dir], digits[dir], factors[dir], nfactors[dir]);
	hypercube_index = node_hc_index(coords[XUP], coords[YUP], coords[ZUP], coords[TUP]);
      }
      if(coord_parity(coords) == parity){
	m2hc[hypercube_index] = morton_index;
	hc2m[morton_index] = hypercube_index;
#if 0
	printf("Morton %d hypercube %d\n", morton_index, hypercube_index);
	FORALLUPDIR(dir){
	  printf("%d: ",coords[dir]);
	  print_factors(digits[dir], nfactors[dir]);
	}
#endif
	morton_index++;
      }
      for(k = 0; k < lastfactor; k++){
	FORALLUPDIR(dir){
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
  int i,j,k,dir;

  if(mynode()==0){
    printf("hyper_prime,");
    printf("\n");
  }

  /* Figure out dimensions of rectangle */
  squaresize[XUP] = nx; squaresize[YUP] = ny;
  squaresize[ZUP] = nz; squaresize[TUP] = nt;
  nsquares[XUP] = nsquares[YUP] = nsquares[ZUP] = nsquares[TUP] = 1;
  
  i = 1;	/* current number of hypercubes */
  while(i<numnodes()){
    /* figure out which prime to divide by starting with largest */
    k = MAXPRIMES-1;
    while( (numnodes()/i)%prime[k] != 0 && k>0 ) --k;
    /* figure out which direction to divide */
    
    /* find largest dimension of h-cubes divisible by prime[k] */
    for(j=0,dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]>j && squaresize[dir]%prime[k]==0 )
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
      if(mynode()==0)
	printf("LAYOUT: Can't lay out this lattice, not enough factors of %d\n"
	       ,prime[k]);
      terminate(1);
    }
    
    /* do the surgery */
    i*=prime[k]; squaresize[dir] /= prime[k]; nsquares[dir] *= prime[k];
  }
}

/*--------------------------------------------------------------------*/

void setup_fixed_geom(int const *geom, int n){
  int i;
  int node_count;
  int len[4];
  int status;

#ifdef FIX_NODE_GEOM
  if(geom != NULL){
      node0_printf("setup_layout: Preallocated machine geometry overrides request\n");
  }
#endif

  len[0] = nx; len[1] = ny; len[2] = nz; len[3] = nt;

  node_count = 1;
  status = 0;
  for(i = 0; i < 4; i++){
    nsquares[i] = geom[i];
    node_count *= geom[i];
    if(len[i] % nsquares[i] != 0)status++;
    squaresize[i] = len[i]/nsquares[i];
  }

  if(node_count != numnodes()){
    node0_printf("/nsetup_fixed_geom: Requested geometry %d %d %d %d ",
		 geom[0], geom[1], geom[2], geom[3]);
    node0_printf("does not match number of nodes %d\n",numnodes());
    terminate(1);
  }

  if(status){
    node0_printf("setup_fixed_geom: Requested geometry %d %d %d %d ",
		 geom[0], geom[1], geom[2], geom[3]);
    node0_printf("is not commensurate with the lattice dims %d %d %d %d\n",
		 nx, ny, nz, nt);
    terminate(1);
  }
}

/*------------------------------------------------------------------*/
/* Initialize io_node function */


#ifdef FIX_IONODE_GEOM

static void init_io_node(){
  int i;
  int status = 0;

  if(ionodegeom() == NULL){
    ionodegeomvals = ionode_geometry;
  } else {
    node0_printf("init_io_node: Command line ionode geometry overrides request\n");
    ionodegeomvals = ionodegeom();
  }

  if(ionodegeomvals == NULL)return;

  /* Compute the number of nodes per I/O node along each direction */
  for(i = 0; i < 4; i++){
    if(nsquares[i] % ionodegeomvals[i] != 0)status++;
    nodes_per_ionode[i] = nsquares[i]/ionodegeomvals[i];
  }
  
  if(status){
    node0_printf("init_io_node: ionode geometry %d %d %d %d \n",
		 ionodegeomvals[0], ionodegeomvals[1],
		 ionodegeomvals[2], ionodegeomvals[3]);
    node0_printf("is incommensurate with node geometry %d %d %d %d\n",
		 nsquares[0], nsquares[1], nsquares[2], nsquares[3]);
    terminate(1);
  }
}
#endif

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
  FORALLUPDIR(dir){
    prime_factors(squaresize[dir], factors[dir], &nfactors[dir], MAXFACTOR);
    printf("%d: ", squaresize[dir]);
    print_factors(factors[dir], nfactors[dir]);
  }
  
  /* The maximum number of Morton digits */
  lastfactor = nfactors[XUP];
  FORALLUPDIRBUT(XUP,dir){
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

void setup_layout(){
  int k = mynode();
  int nd = 0;
  int const *geom;

  if(k == 0) printf("LAYOUT = Hypercubes+Morton, options = ");

#ifdef HAVE_QMP

  /* QMP treatment */

  /* The layout dimensions (geometry) are set as follows:
     1. If the command line has both -qmp-geom and -job-geom we use
        the job geometry. 
     2. Otherwise if -qmp-geom is specified use the allocated geometry
     3. Otherwise if FIX_NODE_GEOM is in force and node_geometry is defined
        use node_geometry
     4. Otherwise use the layout_hyper_prime algorithm to set the geometry
  */
  
  nd = QMP_get_number_of_job_geometry_dimensions();
  if(nd > 0){
    /* Use job geometry */
    geom = QMP_get_job_geometry();
    setup_qmp_grid(geom, nd);
    node0_printf("QMP using job_geometry_dimensions\n");
  } else {
    nd = QMP_get_allocated_number_of_dimensions();
    if(nd > 0) {
      geom = QMP_get_allocated_dimensions();
      /* use allocated geometry */
      setup_qmp_grid(geom, nd);
      node0_printf("QMP using allocated_dimension\n");
    } else {
#ifdef FIX_NODE_GEOM
      if(node_geometry != NULL){
	nd = 4;
	geom = node_geometry;
	/* take geometry from input parameter node_geometry line */
	setup_fixed_geom(geom, nd);
	node0_printf("QMP with specified node_geometry\n");
      } else {
#endif
	setup_hyper_prime();
	nd = 4;
	geom = nsquares;
	node0_printf("QMP with automatic hyper_prime layout\n");
#ifdef FIX_NODE_GEOM
      }
#endif
    }
  }
  
  set_qmp_logical_topology(geom, nd);
  
#else
  
  /* Non QMP treatment */
  
  /* The layout dimensions (geometry) are set as follows:
     1. If the command line has -geom use it
     2. Otherwise, if FIX_NODE_GEOM is in force and the
     node_geometry parameters are specified, use them
     3. Otherwise set the geometry with the layout_hyper_prime 
     algorithm
  */
  
  nd = 4;
  geom = nodegeom();  /* Command line values */
  
#ifdef FIX_NODE_GEOM
  if(geom == NULL){
    geom = node_geometry; /* Input parameter values */
  }
#endif
  
  if(geom != NULL){
    /* Set the sublattice dimensions according to the specified geometry */
    node0_printf("with fixed node_geometry\n");
    setup_fixed_geom(geom, nd);
  } else {
    /* Set the sublattice dimensions according to the hyper_prime algorithm */
    setup_hyper_prime();
    node0_printf("automatic hyper_prime layout\n");
  }
  
#endif
  

  /* Initialize I/O node function */
#ifdef FIX_IONODE_GEOM
  init_io_node();
#endif
  
  /* Compute machine coordinates for this node */
  lex_coords(machine_coordinates, 4, nsquares, k);

  /* Number of sites on node */
  sites_on_node =
    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];

  if( mynode()==0)
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

  even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

/*------------------------------------------------------------------*/
size_t num_sites(int node) {
    return( sites_on_node );
}

/*------------------------------------------------------------------*/
const int *get_logical_dimensions(){
  return nsquares;
}

/*------------------------------------------------------------------*/
/* Coordinates simulate a mesh architecture and must correspond
   to the node_number result */
const int *get_logical_coordinate(){
  return machine_coordinates;
}


/*------------------------------------------------------------------*/
/* io_node(node) maps a node to its I/O node.  The nodes are placed on
   a node lattice with dimensions nsquares.  The I/O partitions are
   hypercubes of the node lattice.  The dimensions of the hypercube are
   given by nodes_per_ionode.  The I/O node is at the origin of that
   hypercube. */


/* Map any node to its I/O node */
int io_node(const int node){
  int i; 
  int io_node_coords[4];

  /* If we don't have I/O partitions, each node does its own I/O */
  if(ionodegeomvals == NULL)
    return node;

  /* Get the machine coordinates for the specified node */
  lex_coords(io_node_coords, 4, nsquares, node);

  /* Round the node coordinates down to get the io_node coordinate */
  for(i = 0; i < 4; i++)
    io_node_coords[i] = nodes_per_ionode[i] * 
      (io_node_coords[i]/nodes_per_ionode[i]);
  
  /* Return the linearized machine coordinates of the I/O node */
  return (int)lex_rank(io_node_coords, 4, nsquares);
}


