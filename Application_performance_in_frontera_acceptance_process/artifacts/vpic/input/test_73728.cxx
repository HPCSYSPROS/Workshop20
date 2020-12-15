////////////////////////////////////////////////////////////////////////
//
//  Reconnection Problem --> single Force-Free Current Sheet with conductive BC
//
///////////////////////////////////////////////////////////////////////

#define NUM_TURNSTILES 8192

// structure to hold the data for energy diagnostics
struct edata {
  species_id sp_id;         /* species id */
  double     vth;          /* thermal energy */
  char fname[256];        /* file to save data */
};

// naming convention for the hydro dump files
#define HYDRO_FILE_FORMAT "hydro/T.%d/%s.%d.%d"

begin_globals {

  int restart_interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int quota_check_interval;  //  How frequently to check if quote exceeded

  int rtoggle;             // enables save of last two restart dumps for safety
  double quota_sec;        // Run quota in seconds
  double b0;               // B0
  double bg;               // Guide field
  double topology_x;       // domain topology
  double topology_y;
  double topology_z;

//  Variables for new output format

  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  // Variables for the energy diagnostics

  edata ede;                        // parameters for electron species
  edata edi;                        // parameters for ion species
  double emax;                       // maximum energy (in units of vth*2/2)
  int nex;                           // number of energy bins

};

begin_initialization {

 // use natural PIC units

 double ec   = 1;         // Charge normalization
 double me   = 1;         // Mass normalization
 double c    = 1;         // Speed of light
 double de   = 1;         // Length normalization (electron inertial length)
 double eps0 = 1;         // Permittivity of space

  double cfl_req   = 0.99;  // How close to Courant should we try to run
  double wpedt_max = 0.36;  // How big a timestep is allowed if Courant is not too restrictive
  double damp      = 0.0;   // Level of radiation damping
  int rng_seed     = 1;     // Random number seed increment

  // Physics parameters

  double mi_me   = 100.0;    // Ion mass / electron mass
  double L_di    = 0.5;      // Sheet thickness / ion inertial length
  double Ti_Te   = 1.0;      // Ion temperature / electron temperature
  double vthe    = 0.15;     // Electron thermal speed over c
  double wpe_wce = 2.0;      // electron plasma freq / electron cyclotron freq
  double bg = 4.0;           // electron plasma freq / electron cyclotron freq
  double theta   = 0;        // B0 = Bx
  double taui    = 2;      // simulation wci's to run

  double quota   = 22.8;   // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds

  double pi = 3.1415927;
  double cs   = cos(theta/180.0*pi);
  double sn   = sin(theta/180.0*pi);

  //derived qunatities
  double mi = me*mi_me;                                   // Ion mass
  double vthi = vthe*sqrt(Ti_Te/mi_me);                   // Ion thermal velocity
  double wci  = 1.0/(mi_me*wpe_wce);                      // Ion cyclotron frequency
  double wce  = wci*mi_me;                                // Electron cyclotron freqeuncy
  double wpe  = wce*wpe_wce;                              // electron plasma frequency
  double wpi  = wpe/sqrt(mi_me);                          // ion plasma frequency
  double di   = c/wpi;                                    // ion inertial length
  double L    = L_di*di;                                  // Sheet thickness in c/wpe

  double ion_sort_interval = 50;        //  Injector moments are also updated at this internal
  double electron_sort_interval=25;    //  Injector moments are also updated at this internal

  // Numerical parameters

  double nppc          = 160; // Average number of macro particle per cell per species

  //double Lx            = 24.0*pi*L_di*di; // size of box in x dimension
  //double Ly            = 24.0*pi*L;     // size of box in y dimension
  //double Lz            = 9.0*pi*L_di*di; // size of box in z dimension

  double Lx            = 40.0*di; // size of box in x dimension
  double Ly            = 40.0*di;     // size of box in y dimension
  double Lz            = 15.0*di; // size of box in z dimension

  double topology_x = 48;  // Number of domains in x, y, and z
  double topology_y = 48;
  double topology_z = 32;

  double nx = 1536;
  double ny = 1536;
  double nz = 1536;

  double hx = Lx/nx;
  double hy = Ly/ny;
  double hz = Lz/nz;

  double b0  = me*c*wce/ec; // Asymptotic magnetic field strength
  double n0  = me*eps0*wpe*wpe/(ec*ec);  // Peak electron (ion) density
  double Ne  = nppc*nx*ny*nz;  // total macro electrons in box
  double Np  = n0*Lx*Ly*Lz;  //  total number of physical electrons
  Ne  = trunc_granular(Ne,nproc()); // Make it divisible by number of processors
  double qe = -ec*Np/Ne;  // Charge per macro electron
  double qi =  ec*Np/Ne;  // Charge per macro ion

  //L.O. switch to perturbation for single layer
  //double eps = 0.0; // Perturbation
  double Lpert = 1.0*Lx;   // wavelength of perturbation
  double dbz =  0.03*b0; //  Perturbation in Bz relative to Bo (Only change here)
  double dbx = -dbz*Lpert/(2*Lz); // Set Bx perturbation so that div(B) = 0

  // Determine the time step

  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);        // courant length
  double dt = cfl_req*dg/c;                             // courant limited time step
  if( wpe*dt>wpedt_max) dt=wpedt_max/wpe;               // override timestep if plasma frequency limited

  //  int restart_interval = int(300.0/(wci*dt));
  int restart_interval = 10000;
  int energies_interval = 500;
  int interval = 100*int(taui/60.0/(wci*dt));
  int fields_interval = interval;
  int ehydro_interval = interval;
  int Hhydro_interval = interval;
  int eparticle_interval = 8*interval;
  int Hparticle_interval = 8*interval;
  int quota_check_interval     = 100;

  //  Determine which domains area along the boundaries - Use macro from grid/partition.c

# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                   \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(topology_x);   /* iy = iy+gpy*iz */                    \
    _ix -= _iy*int(topology_x);   /* ix = ix */                           \
    _iz  = _iy/int(topology_y);   /* iz = iz */                           \
    _iy -= _iz*int(topology_y);   /* iy = iy */ 	        	  \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \
  } END_PRIMITIVE

  int ix, iy, iz ;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );


  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  num_step             = int(taui/(wci*dt));
  status_interval      = 200;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = status_interval/2;
  clean_div_b_interval = status_interval/2;

  global->restart_interval   = restart_interval;
  global->energies_interval  = energies_interval;
  global->fields_interval    = fields_interval;
  global->ehydro_interval    = ehydro_interval;
  global->Hhydro_interval    = Hhydro_interval;
  global->eparticle_interval = eparticle_interval;
  global->Hparticle_interval = Hparticle_interval;
  global->quota_check_interval     = quota_check_interval;
  global->quota_sec          = quota_sec;

  global->rtoggle            = 0;
  global->b0  = b0;
  global->bg  = bg;

  global->topology_x  = topology_x;
  global->topology_y  = topology_y;
  global->topology_z  = topology_z;


  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;
  grid->damp = damp;

  // Define the grid

//L.O.
// define_periodic_grid(  0,   0, -0.5*Lz,                      // Low corner
//			  Lx, Ly, 0.5*Lz,                      // High corner
//			  nx, ny, nz,                          // Resolution
//			  topology_x, topology_y, topology_z); // Topology

  define_periodic_grid(  0, -0.5*Ly, -0.5*Lz,    // Low corner
                          Lx, 0.5*Ly, 0.5*Lz,     // High corner
                          nx, ny, nz,             // Resolution
                          topology_x, topology_y, topology_z); // Topology

 // ***** Set Field Boundary Conditions *****

  // sim_log("Conducting fields on X & Z-boundaries");
  if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY( 0,0,1), pec_fields );

  if ( iz==0 )            set_domain_particle_bc( BOUNDARY(0,0,-1), reflect_particles );
  if ( iz==topology_z-1 ) set_domain_particle_bc( BOUNDARY(0,0,1), reflect_particles );
//LO


  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the species

//L.O.? slight different with Bill's example, should I take care of it?
  sim_log("Setting up species. ");
  species_t *electron = define_species("electron",-ec/me,2.0*Ne/nproc(),-1,electron_sort_interval,0);
  species_t *ion = define_species("ion",ec/mi,2.0*Ne/nproc(),-1,ion_sort_interval,0);

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. ");

  define_material( "vacuum", 1 );

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

////////////////////////////////////////////////////////////////////////////////////////////
//  Finalize Field Advance

  sim_log("Finalizing Field Advance");

  finalize_field_advance(standard_field_advance);

  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z);
  sim_log ( "L_di   = " << L_di );
  sim_log ( "Ti/Te = " << Ti_Te ) ;
  sim_log ( "wpe/wce = " << wpe_wce );
  sim_log ( "mi/me = " << mi_me );
  sim_log ( "theta = " << theta );
  sim_log ( "taui = " << taui );
  sim_log ( "num_step = " << num_step );
  sim_log ( "Lx/di = " << Lx/di );
  sim_log ( "Lx/de = " << Lx/de );
  sim_log ( "Ly/di = " << Ly/di );
  sim_log ( "Ly/de = " << Ly/de );
  sim_log ( "Lz/di = " << Lz/di );
  sim_log ( "Lz/de = " << Lz/de );
  sim_log ( "nx = " << nx );
  sim_log ( "ny = " << ny );
  sim_log ( "nz = " << nz );
  sim_log ( "damp = " << damp );
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << nproc ()  );
  sim_log ( "nppc = " << nppc );
  sim_log ( " b0 = " << b0 );
  sim_log ( " di = " << di );
  sim_log ( " Ne = " << Ne );
  sim_log ( "total # of particles = " << 2*Ne );
  sim_log ( " qi = " << qi );
  sim_log ( " qe = " << qe );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( "dt*wce = " << wce*dt );
  sim_log ( "dt*wci = " << wci*dt );
  sim_log ( " energies_interval: " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/rhoe = " << (Lx/nx)/(vthe/wce)  );
  sim_log ( "L/debye = " << L/(vthe/wpe)  );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "n0 = " << n0 );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );


  // Dump simulation information to file "info"
  if (rank() == 0 ) {

    FileIO fp_info;

    if ( ! (fp_info.open("info", io_write)==ok) ) ERROR(("Cannot open file."));

    fp_info.print( "           ***** Simulation parameters ***** \n");
    fp_info.print( " L/di=%e\n", L_di);
    fp_info.print( " L/de=%e\n", L/de);
    fp_info.print( " Ti/Te=%e\n", Ti_Te );
    fp_info.print( " wpe/wce = %e\n", wpe_wce );
    fp_info.print( " mi/me =%e\n", mi_me );
    fp_info.print( " theta =%e\n", theta );
    fp_info.print( " taui =%e\n", taui );
    fp_info.print( " num_step = %i\n", num_step );
    fp_info.print( " Lx/de = %e\n", Lx/de );
    fp_info.print( " Ly/de = %e\n", Ly/de );
    fp_info.print( " Lz/de =%e\n", Lz/de );
    fp_info.print( " Lx/di = %e\n", Lx/di );
    fp_info.print( " Ly/di = %e\n", Ly/di );
    fp_info.print( " Lz/di =%e\n", Lz/di );
    fp_info.print( " nx = %e\n", nx );
    fp_info.print( " ny = %e\n", ny );
    fp_info.print( " nz =%e\n", nz );
    fp_info.print( " damp =%e\n", damp );
    fp_info.print( " courant = %e\n", c*dt/dg );
    fp_info.print( " nproc = %e\n", nproc() );
    fp_info.print( " nppc = %e\n", nppc );
    fp_info.print( " b0 =%e\n", b0 );
    fp_info.print( " di = %e\n", di );
    fp_info.print( " Ne = %e\n", Ne );
    fp_info.print( " total # of particles = %e\n", 2*Ne );
    fp_info.print( " dt*wpe = %e\n", wpe*dt );
    fp_info.print( " dt*wce = %e\n", wce*dt );
    fp_info.print( " dt*wci = %e\n", wci*dt );
    fp_info.print( " energies_interval: %i\n", energies_interval);
    fp_info.print( " dx/de =%e\n", Lx/(de*nx) );
    fp_info.print( " dy/de =%e\n", Ly/(de*ny) );
    fp_info.print( " dz/de =%e\n", Lz/(de*nz) );
    fp_info.print( " L/debye =%e\n", L/(vthe/wpe) );
    fp_info.print( " dx/rhoi =%e\n", (Lx/nx)/(vthi/wci) );
    fp_info.print( " dx/rhoe = %e\n", (Lx/nx)/(vthe/wce) );
    fp_info.print( " dx/debye = %e\n", (Lx/nx)/(vthe/wpe) );
    fp_info.print( " n0 =            %e\n", n0 );
    fp_info.print( " vthi/c =%e\n", vthi/c );
    fp_info.print( " vthe/c =%e\n", vthe/c );
    fp_info.print( " ***************************\n");
    fp_info.close();


    // for the parallized translate.f90 written by Vadim
    // write binary info file

    if ( ! (fp_info.open("info.bin", io_write)==ok) ) ERROR(("Cannot open file."));

        fp_info.write(&topology_x, 1 );
        fp_info.write(&topology_y, 1 );
        fp_info.write(&topology_z, 1 );

        fp_info.write(&Lx, 1 );
        fp_info.write(&Ly, 1 );
        fp_info.write(&Lz, 1 );

        fp_info.write(&nx, 1 );
        fp_info.write(&ny, 1 );
        fp_info.write(&nz, 1 );

        fp_info.write(&dt, 1 );

        fp_info.write(&mi_me, 1 );
        fp_info.write(&wpe_wce, 1 );
        fp_info.write(&vthe, 1 );
        fp_info.write(&vthi, 1 );

        fp_info.close();

}

  ////////////////////////////
  // Load fields


  // Define some function to load profiles

//L.O.
//# define BX b0*(tanh((z+Lz/4)/L)-tanh((z-Lz/4)/L)-1.0)
//# define BY sqrt(b0*b0 + bg*bg*b0*b0 - BX*BX)
//# define UDY -(b0/L)*(1/(cosh((z+Lz/4)/L)*cosh((z+Lz/4)/L)) - 1/(cosh((z-Lz/4)/L) *cosh((z-Lz/4)/L)))
//# define UDX UDY*BX/BY
//# define PI 3.1415927
//# define DBX -(eps*b0*2*PI/Lz)*sin(2*PI*x/Lx)*sin(4*PI*z/Lz)
//# define DBZ  (eps*b0*2*PI/Lx)*cos(2*PI*x/Lx)*sin(2*PI*z/Lz)*sin(2*PI*z/Lz)

# define BX b0*tanh(z/L)
# define BY sqrt(b0*b0 + bg*bg*b0*b0 - BX*BX)
# define UDY -(b0/L)/(cosh(z/L)*cosh(z/L))
# define UDX UDY*BX/BY
# define PI 3.1415927
# define DBX dbx*cos(2.0*PI*(x-0.5*Lx)/Lpert)*sin(PI*z/Lz)
# define DBZ dbz*cos(PI*z/Lz)*sin(2.0*PI*(x-0.5*Lx)/Lpert)
//LO

  sim_log( "Loading fields" );
  //L.O. set_region_field( everywhere, 0, 0, 0, BX+DBX, BY, DBZ);
   set_region_field( everywhere, 0, 0, 0, (BX+DBX)*cs+BY*sn, -(BX+DBX)*sn+BY*cs, DBZ);

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

  // LOAD PARTICLES

  sim_log( "Loading particles" );

  // Do a fast load of the particles

  seed_rand( rng_seed*nproc() + rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  // Load Harris population

  sim_log( "-> Force Free Sheet" );

  repeat ( Ne/nproc() ) {
    double x, y, z, ux, uy, uz ;

   x = uniform_rand(xmin,xmax);
   y = uniform_rand(ymin,ymax);
   z = uniform_rand(zmin,zmax);

   // inject_particles() will return an error for particles no on this
   // node and will not inject particle locally

   //  Load electrons as drifting Maxwellian with velocity specified to be consistent with B field

    ux = maxwellian_rand(vthe) + UDX;
    uy = maxwellian_rand(vthe) + UDY;
    uz = maxwellian_rand(vthe);

    //L.O. inject_particle(electron, x, y, z, ux, uy, uz, qe, 0, 0 );
    inject_particle(electron, x, y, z, ux*cs+uy*sn,-ux*sn+uy*cs, uz, qe, 0, 0 );

    //  Ions are spatially uniform Maxwellian with no drifts

    ux = maxwellian_rand(vthi);
    uy = maxwellian_rand(vthi);
    uz = maxwellian_rand(vthi);

    //L.O. inject_particle(ion, x, y, z, ux, uy, uz, qi, 0, 0 );
    inject_particle(ion, x, y, z, ux*cs+uy*sn,-ux*sn+uy*cs, uz, qi, 0, 0 );

  }

  sim_log( "Finished loading particles" );

   /*--------------------------------------------------------------------------
     * New dump definition
     *------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------
   * Set data output format
   *
   * This option allows the user to specify the data format for an output
   * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
   * format is the native storage format for data in VPIC.  For field data,
   * this looks something like:
   *
   *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
   *
   * Banded data format stores all data of a particular state variable as a
   * contiguous array, and is easier for ParaView to process efficiently.
   * Banded data looks like:
   *
   *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
   *
     *------------------------------------------------------------------------*/

  global->fdParams.format = band;

  sim_log ( "Fields output format = band" );

  global->hedParams.format = band;

  sim_log ( "Electron species output format = band" );

  global->hHdParams.format = band;

  sim_log ( "Ion species output format = band" );

    /*--------------------------------------------------------------------------
   * Set stride
   *
   * This option allows data down-sampling at output.  Data are down-sampled
   * in each dimension by the stride specified for that dimension.  For
   * example, to down-sample the x-dimension of the field data by a factor
   * of 2, i.e., half as many data will be output, select:
   *
   *   global->fdParams.stride_x = 2;
   *
   * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
   * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
   * Setting the strides in x and y to equal 2 results in an output grid of
   * nx = 4, ny = 4, with actual extents 6x6.
   *
     * G G G G G G G G G
     * G X X X X X X X G
     * G X X X X X X X G         G G G G G G
     * G X X X X X X X G         G X X X X G
     * G X X X X X X X G   ==>   G X X X X G
     * G X X X X X X X G         G X X X X G
     * G X X X X X X X G         G X X X X G
     * G X X X X X X X G         G G G G G G
     * G G G G G G G G G
     *
   * Note that grid extents in each dimension must be evenly divisible by
   * the stride for that dimension:
   *
   *   nx = 150;
   *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
   *
   *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
     *------------------------------------------------------------------------*/

  // relative path to fields data from global header
  sprintf(global->fdParams.baseDir, "fields");

  // base file name for fields output
  sprintf(global->fdParams.baseFileName, "fields");

  global->fdParams.stride_x = 1;
  global->fdParams.stride_y = 1;
  global->fdParams.stride_z = 1;

  // add field parameters to list
  global->outputParams.push_back(&global->fdParams);

  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  // relative path to electron species data from global header
  sprintf(global->hedParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(global->hedParams.baseFileName, "ehydro");

  global->hedParams.stride_x = 1;
  global->hedParams.stride_y = 1;
  global->hedParams.stride_z = 1;

  // add electron species parameters to list
  global->outputParams.push_back(&global->hedParams);

  sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  // relative path to electron species data from global header
  sprintf(global->hHdParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(global->hHdParams.baseFileName, "Hhydro");

  global->hHdParams.stride_x = 1;
  global->hHdParams.stride_y = 1;
  global->hHdParams.stride_z = 1;

  sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

  // add electron species parameters to list
  global->outputParams.push_back(&global->hHdParams);

    /*--------------------------------------------------------------------------
   * Set output fields
   *
   * It is now possible to select which state-variables are output on a
   * per-dump basis.  Variables are selected by passing an or-list of
   * state-variables by name.  For example, to only output the x-component
   * of the electric field and the y-component of the magnetic field, the
   * user would call output_variables like:
   *
   *   global->fdParams.output_variables( ex | cby );
   *
   * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
   * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
   * 'output_variables' WILL HAVE NO EFFECT.
   *
   * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
   * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
   *
   * For convenience, the output variable 'all' is defined:
   *
   *   global->fdParams.output_variables( all );
     *------------------------------------------------------------------------*/
  /* CUT AND PASTE AS A STARTING POINT
   * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

    output_variables( all );

    output_variables( electric | div_e_err | magnetic | div_b_err |
                      tca      | rhob      | current  | rhof |
                      emat     | nmat      | fmat     | cmat );

    output_variables( current_density  | charge_density |
                      momentum_density | ke_density     | stress_tensor );
  */

  global->fdParams.output_variables( electric | magnetic );
  global->hedParams.output_variables( current_density | charge_density | stress_tensor );
  global->hHdParams.output_variables( current_density | charge_density | stress_tensor );

  //global->fdParams.output_variables( all );
  //global->hedParams.output_variables( all );
  //global->hHdParams.output_variables( all );

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  create_field_list(varlist, global->fdParams);

  sim_log ( "Fields variable list: " << varlist );

  create_hydro_list(varlist, global->hedParams);

  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hHdParams);

  sim_log ( "Ion species variable list: " << varlist );


  /* ---------------------------------------------

     now add parameters for the energy diagnostics

   --------------------------------------------- */

  global->ede.sp_id = electron->id;
  global->ede.vth = sqrt(2.0)*vthe;
  sprintf(global->ede.fname,global->hedParams.baseFileName);

  global->edi.sp_id = ion->id;
  global->edi.vth = sqrt(2.0)*vthi;
  sprintf(global->edi.fname, global->hHdParams.baseFileName);

  global->nex  = 5;
  global->emax = 25;

  sim_log("*** Finished with user-specified initialization ***");


  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message

} //begin_initialization

#define should_dump(x) \
  (global->x##_interval>0 && remainder(step, global->x##_interval) == 0)

#include <FileIO.hxx>

begin_diagnostics {


  /*--------------------------------------------------------------------------
   * NOTE: YOU CANNOT DIRECTLY USE C FILE DESCRIPTORS OR SYSTEM CALLS ANYMORE
   *
   * To create a new directory, use:
   *
   *   dump_mkdir("full-path-to-directory/directoryname")
   *
   * To open a file, use: FileIO class
   *
   * Example for file creation and use:
   *
   *   // declare file and open for writing
   *   // possible modes are: io_write, io_read, io_append,
   *   // io_read_write, io_write_read, io_append_read
   *   FileIO fileIO;
   *   FileIOStatus status;
   *   status= fileIO.open("full-path-to-file/filename", io_write);
   *
   *   // formatted ASCII  output
   *   fileIO.print("format string", varg1, varg2, ...);
   *
   *   // binary output
   *   // Write n elements from array data to file.
   *   // T is the type, e.g., if T=double
   *   // fileIO.write(double * data, size_t n);
   *   // All basic types are supported.
   *   fileIO.write(T * data, size_t n);
   *
   *   // close file
   *   fileIO.close();
     *------------------------------------------------------------------------*/

     /*--------------------------------------------------------------------------
   * Data output directories
   * WARNING: The directory list passed to "global_header" must be
   * consistent with the actual directories where fields and species are
   * output using "field_dump" and "hydro_dump".
   *
   * DIRECTORY PATHES SHOULD BE RELATIVE TO
   * THE LOCATION OF THE GLOBAL HEADER!!!
     *------------------------------------------------------------------------*/


  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
  if(step==0) {
    dump_mkdir("fields");
    dump_mkdir("hydro");
    dump_mkdir("rundata");
    dump_mkdir("restart0");
    dump_mkdir("restart1");  // 1st backup
    dump_mkdir("restart2");  // 2nd backup
    dump_mkdir("particle");

    dump_grid("rundata/grid");
    dump_materials("rundata/materials");
    dump_species("rundata/species");
    global_header("global", global->outputParams);
  } // if

  /*--------------------------------------------------------------------------
   * Normal rundata energies dump
   *------------------------------------------------------------------------*/
  if(should_dump(energies)) {
    dump_energies("rundata/energies", step == 0 ? 0 : 1);
  } // if

  /*--------------------------------------------------------------------------
   * Field data output
   *------------------------------------------------------------------------*/

  if(step == 1 || should_dump(fields)) field_dump(global->fdParams);

  /*--------------------------------------------------------------------------
   * Electron species output
   *------------------------------------------------------------------------*/

  if(should_dump(ehydro)) hydro_dump("electron", global->hedParams);

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

  if(should_dump(Hhydro)) hydro_dump("ion", global->hHdParams);

  /*--------------------------------------------------------------------------
   * Energy Spectrum Output
   *------------------------------------------------------------------------*/

  //         #include "energy.cxx"   //  Subroutine to compute energy spectrum diagnostic

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  global->restart_interval=8000;
  if(step && !(step%global->restart_interval)) {
    double dumpstart = mp_elapsed(grid->mp);
    begin_turnstile(NUM_TURNSTILES);
    if(!global->rtoggle) {
      global->rtoggle = 1;
      dump_restart("restart1/restart", 0);
    }
    else {
      global->rtoggle = 0;
      dump_restart("restart2/restart", 0);
    } // if
    end_turnstile;
    double dumpelapsed = mp_elapsed(grid->mp) - dumpstart;
    sim_log("Restart duration "<<dumpelapsed);
  } // if


  // Dump particle data

  char subdir[36];
  if ( should_dump(eparticle) && step !=0 && step > 20*(global->fields_interval)  ) {
      //	if ( should_dump(eparticle) && step !=0 ) {
    sprintf(subdir,"particle/T.%d",step);
    dump_mkdir(subdir);
    sprintf(subdir,"particle/T.%d/eparticle",step);
    dump_particles("electron",subdir);
  }

//   if ( should_dump(Hparticle) ) {
//     dump_particles("ion",  "Hparticle");
//   }

  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the mp_elapsed() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive mp_elapsed call from every
  // timestep. mp_elapsed has an ALL_REDUCE in it!

  //   	global->quota_sec = 7.5*3600;
  if  ( step>0 && global->quota_check_interval>0 && (step&global->quota_check_interval)==0 ) {
    if( mp_elapsed( grid->mp ) > global->quota_sec ) {
      sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
      double dumpstart = mp_elapsed(grid->mp);
      begin_turnstile(NUM_TURNSTILES);
      dump_restart("restart0/restart",0);
      end_turnstile;
      mp_barrier( grid->mp ); // Just to be safe
      sim_log( "Restart dump restart completed." );
      double dumpelapsed = mp_elapsed(grid->mp) - dumpstart;
      sim_log("Restart duration "<<dumpelapsed);
      mp_finalize( grid->mp );
      exit(0); // Exit or abort?
    }
  }

} // end diagnostics

// *******************  PARTICLE INJECTION  - OPEN BOUNDARY ***************************

begin_particle_injection {


}  // end particle injection


//   *******************  CURRENT INJECTION ***************************

begin_current_injection {

  // No current injection for this simulation

}

//   *******************  FIELD INJECTION ***************************

begin_field_injection {


}  // end field injection


begin_particle_collisions {

  // No particle collisions in this simulation


}
