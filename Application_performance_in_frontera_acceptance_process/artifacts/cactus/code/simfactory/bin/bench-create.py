#!/usr/bin/env python

import optparse
import os
import re

def flatten(l):
    for el in l:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def ifthen(cond, iftrue, iffalse):
    if cond: return iftrue
    return iffalse

def unique(elts):
    # return sorted(list(set(elts)))
    keys = {}
    for elt in elts:
        keys[elt] = 1
    uelts = keys.keys()
    uelts.sort()
    return uelts



# User options

version = '3'                   # Increase this after every release

parser = optparse.OptionParser()
parser.add_option('--path', type='string', default='bench',
                  help='output directory')
parser.add_option('--prefix', type='string', default='bench',
                  help='benchmark name prefix (cannot be empty)')
parser.add_option('--machine', type='string',
                  help='target machine name prefix')
parser.add_option('--configuration', type='string',
                  help='configuration name')
parser.add_option('--optionlist', type='string',
                  help='optionlist', default='')
parser.add_option('--submitscript', type='string',
                  help='submitscript', default='')
parser.add_option('--runscript', type='string',
                  help='runscript', default='')
parser.add_option('--memory', type='float', default=500.0e+6,
                  help='memory to use per node in bytes')
parser.add_option('--walltime', type='float', default=300.0,
                  help='wall time to use in seconds')
parser.add_option('--levels', action='append', type='int', default=[],
                  help='list of number of refinement levels')
parser.add_option('--fdorders', action='append', type='int', default=[],
                  help='list of finite differencing orders')
parser.add_option('--multipatch', action='append', type='int', default=[],
                  help='list of 0 or 1')
parser.add_option('--hydro', action='append', type='int', default=[],
                  help='list of 0 or 1')
parser.add_option('--flop-per-second', type='float', default=(2.0e+9*4),
                  help='theoretical maximum flop per second for each core')
parser.add_option('--efficiency', type='float', default=(0.1),
                  help='estimated floating point efficiency')
parser.add_option('--min-cores', type='int', default=1,
                  help='minimum number of cores to use; should probably be the number of cores per node')
parser.add_option('--max-cores', type='int', default=1,
                  help='maximum number of cores to use')
parser.add_option('--min-nodes', type='int', default=1,
                  help='minimum number of nodes to use')
parser.add_option('--max-nodes', type='int', default=1,
                  help='maximum number of nodes to use')
parser.add_option('--threads', action='append', type='int', default=[],
                  help='list of numbers of threads to use')
parser.add_option('--smt', action='append', type='int', default=[],
                  help='list of numbers of SMT threads to use')
parser.add_option('--ppn', type='int',
                  help='numbers of cores per node')
parser.add_option('--ppn-used', action='append', type='int', default=[],
                  help='list of numbers of cores per node to use, or -1 for unspecified')
parser.add_option('--run', action='append', type='int', default=[],
                  help='list of integers to identify benchmark runs')
(options, args) = parser.parse_args()
if args:
    raise 'Unused command line arguments'
print options

path            = options.path
prefix          = options.prefix
machine         = options.machine
configuration   = options.configuration
optionlist      = options.optionlist
submitscript    = options.submitscript
runscript       = options.runscript
memory          = options.memory
walltime        = options.walltime
levels_list     = unique(options.levels)
fdorder_list    = unique(options.fdorders)
multipatch_list = unique(options.multipatch)
hydro_list      = unique(options.hydro)
flop_per_second = options.flop_per_second
efficiency      = options.efficiency
min_cores       = options.min_cores
max_cores       = options.max_cores
min_nodes       = options.min_nodes
max_nodes       = options.max_nodes
threads_list    = unique(options.threads)
smt_list        = unique(options.smt)
ppn             = options.ppn
ppn_used_list   = unique(options.ppn_used)
run_list        = unique(options.run)

# We can't have default values for appending options, since the
# default is then always present as first element
if not levels_list:
    levels_list = [1]
if not fdorder_list:
    fdorder_list = [4]
if not multipatch_list:
    multipatch_list = [0]
if not hydro_list:
    hydro_list = [0]
if not threads_list:
    threads_list = [1]
if not smt_list:
    smt_list = [1]
if not ppn_used_list:
    ppn_used_list = [0]
if not run_list:
    run_list = [0]

print 'Benchmark suite parameters:'
print
print '   Path:               %s' % path
print '   Prefix:             %s' % prefix
print '   Memory [bytes]:     %g' % memory
print '   Walltime [seconds]: %g' % walltime
print '   Levels:             %s' % levels_list
print '   FD orders:          %s' % fdorder_list
print '   Multipatch:         %s' % multipatch_list
print '   Hydro:              %s' % hydro_list
print '   Min cores:          %d' % min_cores
print '   Max cores:          %d' % max_cores
print '   Min nodes:          %d' % min_nodes
print '   Max nodes:          %d' % max_nodes
print '   Threads:            %s' % threads_list
print '   SMT:                %s' % smt_list
print '   PPN:                %s' % ppn
print '   PPN used:           %s' % ppn_used_list
print '   Runs:               %s' % run_list



# Write out one benchmark
def create_benchmark(hydro, multipatch, levels, fdorder,
                     nodes, cores, ppn_used, procs, threads, smt):
    
    # Basic settings
    maxlevels = 10
    
    minsize = ifthen(hydro, 10.0, 1.0) # interior of finest grid
    
    numvars = 25 + ifthen(hydro, 10, 0) # BSSN / GRHydro
    timelevels      = 5                 # RK4
    bytes_per_var   = 8                 # double
    ghosts          = fdorder/2+1
    memory_overhead = 1.5               # estimate
    
    amr_overhead   = 1.0 + 0.5 * (levels-1) # estimate
    ops_per_point  = 5000 + ifthen(hydro, 5000, 0) # BSSN / GRHydro
    steps_per_iter = 4                             # RK4
    
    
    
    buffers = ghosts * (steps_per_iter - 1)
    
    # Determine grid size
    
    # Bytes per core per grid point
    bytes_per_point = (
        1.0 * timelevels * numvars * bytes_per_var * memory_overhead)
    
    # Memory for a given grid size
    def grid_memory(gridsize):
        # TODO: take additional overlap points into account
        npoints = gridsize**3
        return bytes_per_point * npoints * (levels + 6*multipatch)
    def grid_memory_with_ghosts(gridsize):
        # TODO: take additional overlap points into account
        npoints = (
            (gridsize / procs**(1.0/3.0) + 2*ghosts)**3 * procs)
        return bytes_per_point * npoints * (levels + 6*multipatch)
    
    min_gridsize = ifthen(levels==1, 1, buffers)
    gridsize = min_gridsize - 1
    while grid_memory(gridsize) < nodes * memory:
        gridsize += 1
    if gridsize < min_gridsize:
        print (
            "Not enough memory for levels=%d nodes=%d cores=%d hydro=%d multipatch=%d with memory=%g" %
            (levels, nodes, cores, hydro, multipatch, memory))
        return
    real_memory = grid_memory_with_ghosts(gridsize)
    
    # Determine iterations
    
    seconds_per_point = (
        1.0 * ops_per_point * steps_per_iter * amr_overhead /
        (cores * flop_per_second * efficiency))
    
    # Wall time for a single iteration
    def iteration_walltime(iteration):
        # TODO: take additional overlap points into account
        npoints = gridsize**3
        
        count = 0
        for level in range(0, levels):
            if iteration % (2 ** (maxlevels - level - 1)) != 0:
                continue
            count += npoints
            if multipatch and level==0:
                count += 6*npoints
        return seconds_per_point * count
    
    iterations = 0
    real_walltime = 0.0
    while real_walltime < walltime:
        real_walltime += iteration_walltime(iterations)
        iterations += 1
    
    
    
    # Dependent parameters
    
    # Add buffers explicitly to coarse grid extent, so
    # that all levels have the same number of active grid
    # points
    coarsesize = (
        1.0 * minsize / (gridsize - buffers) * gridsize * (2 ** (levels - 1)))
    
    # bssn = ifthen(multipatch, 'ML_BSSN_MP', 'ML_BSSN')
    bssn = 'ML_BSSN'
    
    filename = ('%s-%s-%s-lev%02d-fdo%01d-grid%06d-iter%06d' %
                (prefix,
                 ifthen(hydro, 'tov', 'minkowski'),
                 ifthen(multipatch, 'mp', 'amr'),
                 levels, fdorder, gridsize, iterations))
    exename = ('%s-%s-%s-lev%02d-fdo%01d' %
               (prefix,
                ifthen(hydro, 'tov', 'minkowski'),
                ifthen(multipatch, 'mp', 'amr'),
                levels, fdorder))
    
    
    
    print
    print
    print
    print 'Benchmark parameters:'
    print
    print '   Physics:           %s' % ifthen(hydro, 'TOV', 'Minkowski')
    print '   Grid:              %s' % ifthen(multipatch, 'multi-block', 'AMR')
    print '   Levels:            %d' % levels
    print '   FD order:          %d' % fdorder
    print '   Cores:             %d' % cores
    print '   Total grid points: %d' % gridsize**3
    print '   Iterations:        %d' % iterations
    print
    print '   Memory requested [bytes]:     %g' % memory
    print ('   Memory used      [bytes]:     %g (%.1f%%)' %
           (real_memory, 100.0*real_memory/memory))
    print '   Walltime requested [seconds]: %g' % walltime
    print ('   Walltime real      [seconds]: %g (%.1f%%)' %
           (walltime, 100.0*real_walltime/walltime))
    print
    print '   Parameter file name: %s.par' % filename
    
    
    
    lines = [
        'Cactus::cctk_run_title = "Einstein Toolkit Benchmark version %s"' % version,
        '',
        '# Target requirements:',
        '#    Physics:           %s' % ifthen(hydro, 'TOV', 'Minkowski'),
        '#    Refinement levels: %d' % levels,
        '#    Memory:            %g bytes' % memory,
        '#    Wall time:         %g seconds' % walltime,
        '#    Number of cores:   %d' % cores,
        '# Estimated actual requirements:',
        '#    Memory:            %g bytes' % real_memory,
        '#    Wall time:         %g seconds' % real_walltime,
        '# Specifications:',
        '#    Grid size:         %d' % gridsize**3,
        '#    Iterations:        %d' % iterations,
        '',
        '',
        '',
        'Cactus::cctk_full_warnings         = yes',
        'Cactus::highlight_warning_messages = no',
        '',
        'Cactus::cctk_itlast = %d' % iterations,
        '',
        '',
        '',
        'ActiveThorns = "IOUtil"',
        '',
        'IO::out_dir = $parfile',
        '',
        '',
        ifthen(
            hydro,
            [
                '',
                'ActiveThorns = "Constants"'],
            []),
        '',
        'ActiveThorns = "CycleClock"',
        '',
        'ActiveThorns = "GenericFD"',
        '',
        'ActiveThorns = "hwloc"',
        '',
        'ActiveThorns = "InitBase"',
        # '',
        # 'ActiveThorns = "PAPI"',
        '',
        'ActiveThorns = "SystemTopology"',
        '',
        'ActiveThorns = "Timers"',
        '',
        'ActiveThorns = "Vectors"',
        '',
        '',
        '',
        'ActiveThorns = "LoopControl"',
        '',
        'LoopControl::use_smt_threads = yes',
        '',
        'LoopControl::tilesize_i = 4',
        'LoopControl::tilesize_j = 4',
        'LoopControl::tilesize_k = 4',
        '',
        'LoopControl::statistics_filename = ""',
        '',
        '',
        '',
        'ActiveThorns = "Carpet CarpetLib CarpetInterp CarpetInterp2 CarpetReduce"',
        '',
        'Carpet::verbose           = no',
        'Carpet::veryverbose       = no',
        'Carpet::schedule_barriers = no',
        'Carpet::storage_verbose   = no',
        'CarpetLib::output_bboxes  = no',
        'Timers::verbose           = no',
        '',
        ifthen(
            not multipatch,
            [
                'Carpet::domain_from_coordbase  = yes'],
            [
                'Carpet::domain_from_multipatch = yes']),
        'Carpet::max_refinement_levels  = %d' % maxlevels,
        '',
        'driver::ghost_size       = %d' % ghosts,
        'Carpet::use_buffer_zones = yes',
        '',
        'Carpet::prolongation_order_space = %d' % min(5, fdorder+1),
        'Carpet::prolongation_order_time  = 2',
        '',
        'Carpet::init_fill_timelevels = yes',
        '',
        'Carpet::poison_new_timelevels = no',
        'CarpetLib::poison_new_memory  = no',
        '',
        'CarpetLib::pad_to_cachelines = no # TODO yes',
        '',
        'Carpet::output_timers_every      = 0 # %d' % iterations,
        'CarpetLib::print_timestats_every = 0 # %d' % iterations,
        'CarpetLib::print_memstats_every  = %d' % (iterations/2),
        '',
        '',
        '',
        ifthen(
            not multipatch,
            [
                'ActiveThorns = "Boundary CartGrid3D CoordBase ReflectionSymmetry SymBase"',
                '',
                'CoordBase::domainsize = "minmax"',
                'CoordBase::spacing    = "numcells"',
                '',
                'CoordBase::xmin     = 0.0',
                'CoordBase::ymin     = 0.0',
                'CoordBase::zmin     = 0.0',
                'CoordBase::xmax     = %d' % coarsesize,
                'CoordBase::ymax     = %d' % coarsesize,
                'CoordBase::zmax     = %d' % coarsesize,
                'CoordBase::ncells_x = %d' % gridsize,
                'CoordBase::ncells_y = %d' % gridsize,
                'CoordBase::ncells_z = %d' % gridsize,
                '',
                'CoordBase::boundary_shiftout_x_lower = 1',
                'CoordBase::boundary_shiftout_y_lower = 1',
                'CoordBase::boundary_shiftout_z_lower = 1',
                'CoordBase::boundary_size_x_lower     = %d' % ghosts,
                'CoordBase::boundary_size_y_lower     = %d' % ghosts,
                'CoordBase::boundary_size_z_lower     = %d' % ghosts,
                'CoordBase::boundary_size_x_upper     = %d' % ghosts,
                'CoordBase::boundary_size_y_upper     = %d' % ghosts,
                'CoordBase::boundary_size_z_upper     = %d' % ghosts,
                '',
                'CartGrid3D::type = "coordbase"',
                '',
                'ReflectionSymmetry::reflection_x   = yes',
                'ReflectionSymmetry::reflection_y   = yes',
                'ReflectionSymmetry::reflection_z   = yes',
                'ReflectionSymmetry::avoid_origin_x = no',
                'ReflectionSymmetry::avoid_origin_y = no',
                'ReflectionSymmetry::avoid_origin_z = no'],
            [
                'ActiveThorns = "Boundary CartGrid3D CoordBase Coordinates Interpolate2 SymBase TensorTypes"',
                '',
                'Coordinates::coordinate_system = "Thornburg04"',
                'Coordinates::h_cartesian       = %.17g' % (1.0*coarsesize/(gridsize-1)),
                'Coordinates::h_radial          = %.17g' % (2.0*coarsesize/(gridsize-1)),
                '',
                'Coordinates::sphere_inner_radius = %.17g' % (0.5*coarsesize),
                'Coordinates::sphere_outer_radius = %.17g' % (2.5*coarsesize),
                'Coordinates::n_angular           = %d' % (gridsize-1),
                '',
                'Coordinates::patch_boundary_size     = %d' % ghosts,
                'Coordinates::additional_overlap_size = %d' % ghosts,
                'Coordinates::outer_boundary_size     = %d' % ghosts,
                '',
                'Interpolate::interpolator_order = %d' % min(4, fdorder),
                'CarpetInterp::tree_search       = yes',
                'CarpetInterp::check_tree_search = no',
                '',
                'CartGrid3D::type                     = "multipatch"',
                'CartGrid3D::set_coordinate_ranges_on = "all maps"']),
        ifthen(
            levels>1,
            [
                '',
                '',
                '',
                'ActiveThorns = "CarpetRegrid2"',
                '',
                'CarpetRegrid2::regrid_every = 0',
                'CarpetRegrid2::num_centres  = 1',
                '',
                'CarpetRegrid2::num_levels_1 = %d' % levels,
                [
                    ('CarpetRegrid2::radius_1[%d] = %.17g' %
                     (level, 1.0 * minsize * (2 ** (levels - level - 1))))
                    for level in range(1, levels)]],
            []),
        '',
        '',
        '',
        'ActiveThorns = "SphericalSurface"',
        '',
        '',
        '',
        'ActiveThorns = "MoL Time"',
        '',
        'MoL::ODE_Method             = "RK4"',
        'MoL::MoL_Intermediate_Steps = 4',
        'MoL::MoL_Num_Scratch_Levels = 1',
        '',
        'MoL::skip_initial_copy = yes',
        'MoL::init_RHS_zero     = no',
        '',
        '# Evolve for a short time only to avoid any instabilities',
        'Time::dtfac = 0.25e-6',
        '',
        '',
        '',
        'ActiveThorns = "ADMBase ADMCoupling ADMMacros CoordGauge SpaceMask StaticConformal TmunuBase"',
        ifthen(
            hydro,
            [
                '',
                '',
                '',
                'Activethorns = "HydroBase"',
                '',
                'HydroBase::timelevels = 3',
                '',
                'TmunuBase::stress_energy_storage           = yes',
                'TmunuBase::stress_energy_at_RHS            = yes',
                'TmunuBase::timelevels                      = 1',
                'TmunuBase::prolongation_type               = "none"',
                'TmunuBase::support_old_CalcTmunu_mechanism = no',
                '',
                'SpaceMask::use_mask = yes'],
            []),
        '',
        '',
        '',
        'ActiveThorns = "%s %s_Helper NewRad"' % (bssn, bssn),
        # Dissipation
        ifthen(
            hydro,
            [],
            [
                '',
                'ADMBase::initial_data    = "Cartesian Minkowski"',
                'ADMBase::initial_lapse   = "one"',
                'ADMBase::initial_shift   = "zero"',
                'ADMBase::initial_dtlapse = "zero"',
                'ADMBase::initial_dtshift = "zero"']),
        '',
        'ADMBase::evolution_method         = "%s"' % bssn,
        'ADMBase::lapse_evolution_method   = "%s"' % bssn,
        'ADMBase::shift_evolution_method   = "%s"' % bssn,
        'ADMBase::dtlapse_evolution_method = "%s"' % bssn,
        'ADMBase::dtshift_evolution_method = "%s"' % bssn,
        '',
        'GenericFD::assume_stress_energy_state = %d' % ifthen(hydro, 1, 0),
        'GenericFD::assume_use_jacobian        = %d' % ifthen(multipatch, 1, 0),
        ifthen(
            multipatch,
            [
                'GenericFD::jacobian_group            = "Coordinates::jacobian"',
                'GenericFD::jacobian_derivative_group = "Coordinates::jacobian2"'],
            []),
        '',
        '%s::timelevels = 3' % bssn,
        '',
        '%s::fdOrder = %d' % (bssn, fdorder),
        ''
        '%s::conformalMethod  = 0' % bssn,
        '%s::shiftFormulation = 0' % bssn,
        '%s::shiftAlphaPower  = 0' % bssn,
        '',
        '# Disable BSSN constraints',
        '%s::other_timelevels                          =  0' % bssn,
        '%s::%s_ConstraintsEverywhere_calc_offset = -1' % (bssn, bssn),
        '%s::%s_ConstraintsInterior_calc_offset   = -1' % (bssn, bssn),
        '',
        '%s::harmonicF                    = 2.0    # 1+log' % bssn,
        '%s::harmonicN                    = 1      # 1+log' % bssn,
        '%s::shiftGammaCoeff              = 0.75'           % bssn,
        '%s::evolveA                      = 0'              % bssn,
        '%s::evolveB                      = 1'              % bssn,
        '%s::alphaDriver                  = 1.0'            % bssn,
        '%s::betaDriver                   = 1.0'            % bssn,
        '%s::useSpatialBetaDriver         = 0'              % bssn,
        ifthen(
            hydro,
            [
                '%s::useSpatialShiftGammaCoeff    = 1'              % bssn,
                '%s::spatialShiftGammaCoeffRadius = 20.0'           % bssn],
            [
                '%s::useSpatialShiftGammaCoeff    = 0'              % bssn,
                '%s::spatialShiftGammaCoeffRadius = 1.0e+12'        % bssn]),
        '%s::advectLapse                  = 1'              % bssn,
        '%s::advectShift                  = 1'              % bssn,
        '',
        '%s::epsDiss = 0.2' % bssn,
        '',
        '%s::initial_boundary_condition = "extrapolate-gammas"' % bssn,
        '%s::rhs_boundary_condition     = "NewRad"'             % bssn,
        'Boundary::radpower                 = 2',
        '',
        '%s::minimumLapse = 1.0e-8' % bssn,
        '',
        '%s::rhs_evaluation = "splitBy"' % bssn,
        # '',
        # 'Dissipation::epsdis = 0.2',
        # 'Dissipation::order  = %d' % (fdorder+1),
        # 'Dissipation::vars   = "',
        # '        %s::ML_confac'     % bssn,
        # '        %s::ML_metric'     % bssn,
        # '        %s::ML_Gamma'      % bssn,
        # '        %s::ML_trace_curv' % bssn,
        # '        %s::ML_curv'       % bssn,
        # '        %s::ML_lapse'      % bssn,
        # '        %s::ML_dtlapse'    % bssn,
        # '        %s::ML_shift'      % bssn,
        # '        %s::ML_dtshift'    % bssn,
        # '"',
        ifthen(
            hydro,
            [
                '',
                '',
                '',
                'ActiveThorns = "EOS_Omni GRHydro"',
                '',
                'HydroBase::evolution_method = "GRHydro"',
                '',
                'EOS_Omni::poly_gamma =   2.0',
                'EOS_Omni::poly_k     = 100.0',
                '',
                'GRHydro::Riemann_solver    = "Marquina"',
                'GRHydro::GRHydro_EOS_type  = "General"',
                'GRHydro::GRHydro_EOS_table = "Ideal_Fluid"',
                'GRHydro::recon_method      = "ppm"',
                'GRHydro::GRHydro_stencil   = 3',
                'GRHydro::bound             = "none"',
                'GRHydro::rho_abs_min       = 1.e-10',
                '',
                '',
                '',
                'ActiveThorns = "TOVSolver"',
                '',
                'ADMBase::initial_data    = "TOV"',
                'ADMBase::initial_lapse   = "TOV"',
                'ADMBase::initial_shift   = "zero"',
                'ADMBase::initial_dtlapse = "zero"',
                'ADMBase::initial_dtshift = "zero"',
                '',
                'TOVSolver::TOV_Rho_Central[0] = 1.28e-3',
                'TOVSolver::TOV_Gamma[0]       = 2.0',
                'TOVSolver::TOV_K[0]           = 100.0'],
            []),
        '',
        '',
        '',
        'ActiveThorns = "CarpetIOBasic"',
        '',
        'IOBasic::outInfo_every      = 1',
        'IOBasic::outInfo_reductions = "norm2"',
        'IOBasic::outInfo_vars       = "',
        '        Carpet::physical_time_per_hour',
        '        Carpet::local_grid_points_per_second',
        '"',
        '',
        '',
        '',
        'ActiveThorns = "CarpetIOScalar"',
        '',
        'IOScalar::one_file_per_group = yes',
        '',
        'IOScalar::outScalar_every = %d' % iterations,
        'IOScalar::outScalar_vars  = "',
        '        ADMBase::lapse',
        ifthen(
            hydro,
            [
                '        HydroBase::rho'],
            []),
        '"',
        '',
        '',
        '',
        'ActiveThorns = "CarpetIOASCII"',
        '',
        'IOASCII::one_file_per_group = yes',
        '',
        'IOASCII::out3D_ghosts = no',
        '',
        'IOASCII::out0D_every = %d' % iterations,
        'IOASCII::out0D_vars  = "',
        '        Carpet::timing',
        '"',
        # '',
        # 'IOASCII::out1D_every = %d' % iterations,
        # 'IOASCII::out1D_d     = no',
        # 'IOASCII::out1D_vars  = "',
        # ifthen(
        #     multipatch,
        #     [
        #     '        grid::coordinates'],
        #     []),
        # '        ADMBase::lapse',
        # '        %s::ML_lapse' % bssn,
        # ifthen(
        #     hydro,
        #     [
        #         '        HydroBase::rho',
        #         '        GRHydro::dens'],
        #     []),
        # '"',
        '',
        '',
        '',
        'ActiveThorns = "NaNChecker"',
        '',
        'NaNChecker::check_every     = 1',
        'NaNChecker::action_if_found = "terminate"',
        'NaNChecker::check_vars      = "',
        '        ADMBase::lapse',
        ifthen(
            hydro,
            [
                '        HydroBase::rho'],
            []),
        '"',
        '',
        '',
        '',
        'ActiveThorns = "Formaline"',
        '',
        '',
        '',
        'ActiveThorns = "TimerReport"',
        '',
        'CycleClock::register_clock = no',
        '',
        'Carpet::output_initialise_timer_tree = yes',
        'Carpet::output_timer_tree_every      = %d' % iterations,
        '',
        'TimerReport::out_every                  = %d' % iterations,
        '# TimerReport::out_filename               = "TimerReport"',
        'TimerReport::output_all_timers_together = yes',
        'TimerReport::output_all_timers_readable = yes',
        'TimerReport::n_top_timers               = 100']
    
    # TODO: use defaults if not found in parameter file
    fixed_parameters = [
        '%s::fdOrder' % bssn,
        '%s::conformalMethod' % bssn,
        '%s::shiftFormulation' % bssn,
        '%s::shiftAlphaPower' % bssn,
        '%s::%s_ConstraintsEverywhere_calc_offset' % (bssn, bssn),
        '%s::%s_ConstraintsInterior_calc_offset' % (bssn, bssn),
        '%s::harmonicF' % bssn,
        '%s::harmonicN' % bssn,
        '%s::shiftGammaCoeff' % bssn,
        '%s::evolveA' % bssn,
        '%s::evolveB' % bssn,
        '%s::alphaDriver' % bssn,
        '%s::betaDriver' % bssn,
        '%s::useSpatialShiftGammaCoeff' % bssn,
        '%s::spatialShiftGammaCoeffRadius' % bssn,
        '%s::advectLapse' % bssn,
        '%s::advectShift' % bssn,
        '%s::epsDiss' % bssn,
        '%s::minimumLapse' % bssn,
        'GenericFD::assume_stress_energy_state',
        'GenericFD::assume_use_jacobian']
    
    # TODO: determine thorns automatically
    thorns = [
        ifthen(multipatch, '', '# ') + 'AEIThorns/TensorTypes',
        'CactusBase/Boundary',
        'CactusBase/CartGrid3D',
        'CactusBase/CoordBase',
        'CactusBase/IOUtil',
        'CactusBase/InitBase',
        'CactusBase/SymBase',
        'CactusBase/Time',
        # 'CactusNumerical/Dissipation',
        'CactusNumerical/MoL',
        ifthen(not multipatch, '', '# ') + 'CactusNumerical/ReflectionSymmetry',
        'CactusNumerical/SpaceMask',
        'CactusNumerical/SphericalSurface',
        'CactusUtils/Formaline',
        'CactusUtils/NaNChecker',
        'CactusUtils/SystemTopology',
        'CactusUtils/TimerReport',
        'CactusUtils/Vectors',
        'Carpet/Carpet',
        'Carpet/CarpetIOASCII',
        'Carpet/CarpetIOBasic',
        'Carpet/CarpetIOScalar',
        'Carpet/CarpetInterp',
        'Carpet/CarpetInterp2',
        'Carpet/CarpetLib',
        'Carpet/CarpetReduce',
        ifthen(levels>1, '', '# ') + 'Carpet/CarpetRegrid2',
        'Carpet/CycleClock',
        'Carpet/LoopControl',
        'Carpet/Timers',
        'EinsteinBase/ADMBase',
        'EinsteinBase/ADMCoupling',
        'EinsteinBase/ADMMacros',
        ifthen(hydro, '', '# ') + 'EinsteinBase/Constants',
        'EinsteinBase/CoordGauge',
        ifthen(hydro, '', '# ') + 'EinsteinBase/HydroBase',
        'EinsteinBase/StaticConformal',
        'EinsteinBase/TmunuBase',
        ifthen(hydro, '', '# ') + 'EinsteinEOS/EOS_Omni',
        ifthen(hydro, '', '# ') + 'EinsteinEvolve/GRHydro',
        'EinsteinEvolve/NewRad',
        ifthen(hydro, '', '# ') + 'EinsteinInitialData/TOVSolver',
        'EinsteinUtils/TGRtensor',
        # 'ExternalLibraries/BLAS',
        # 'ExternalLibraries/LAPACK',
        'ExternalLibraries/MPI',
        'ExternalLibraries/hwloc',
        'ExternalLibraries/OpenBLAS',
        'ExternalLibraries/zlib',
        'KrancNumericalTools/GenericFD',
        ifthen(multipatch, '', '# ') + 'Llama/Coordinates',
        ifthen(multipatch, '', '# ') + 'Llama/Interpolate2',
        'McLachlan/%s' % bssn,
        'McLachlan/%s_Helper' % bssn]
    
    parfile = open('%s/%s.par' % (path, filename), 'w')
    parlines = []
    for line in flatten(lines):
        parfile.write('%s\n' % line)
        parlines.append(line)
    parfile.close()
    
    fixedparlines = []
    for fp in fixed_parameters:
        (thorn, parameter) = fp.split('::')
        thorn = thorn.upper()
        for parline in parlines:
            m = re.match(r'^\s*(?P<thorn>\w+)'
                         r'\s*::\s*(?P<parameter>\w+)'
                         r'\s*=\s*(?P<value>[^#]+?)\s*(?:#.*)?$', parline)
            if (m and
                m.group('thorn').upper() == thorn and
                m.group('parameter') == parameter):
                break
        else:
            raise RuntimeError("Couldn't find parameter to fix (%s::%s)" % (thorn,parameter))
        value = m.group('value')
        fixedparline = ('-DCCTK_PARAMETER__%s__%s=%s' %
                        (thorn, parameter, value))
        fixedparlines.append(fixedparline)
    
    thornfile = open('%s/%s.th' % (path, exename), 'w')
    thornfile.write('# Custom thorn list for executable %s\n' % exename)
    thornfile.write('\n')
    for thorn in thorns:
        thornfile.write('%s\n' % thorn)
    thornfile.close()
    
    cmd = ['./simfactory/bin/sim',
           'build',
           exename,
           ifthen(machine, '--machine=%s' % machine, ''),
           '--clean',
           ifthen(optionlist, '--optionlist=%s' % optionlist, ''),
           ifthen(submitscript, '--submitscript=%s' % submitscript, ''),
           ifthen(runscript, '--runscript=%s' % runscript, ''),
           '--thornlist=%s/%s.th' % (path, exename),
           ifthen('-aligned' in exename,
                  '--replace=VECTORISE_ALIGNED_ARRAYS=yes', ''),
           # TODO: This assumes that fixedparlines contain neither
           # single quotes nor spaces
           ('--append=CPPFLAGS=\'-DCARPET_OPTIMISE -DNDEBUG %s\'' %
            ' '.join(fixedparlines)),
           ('--append=FPPFLAGS=\'-DCARPET_OPTIMISE -DNDEBUG %s\'' % 
            ' '.join(fixedparlines))]
    buildlines.append(' '.join(cmd))
    
    for run in run_list:
        sim = ('%s-nodes%04d-cores%06d-procs%06d-threads%04d-smt%02d-run%04d' %
               (filename, nodes, cores, procs, threads, smt, run))
        cmd = ['./simfactory/bin/sim',
               'purge',
               sim,
               ifthen(machine, '--machine=%s' % machine, '')]
        submit.write(' '.join(cmd) + '\n')
        cmd = ['./simfactory/bin/sim',
               'submit',
               sim,
               ifthen(machine, '--machine=%s' % machine, ''),
               '--configuration=%s' % exename,
               '--parfile=%s/%s.par' % (path, filename),
               '--procs=%d' % (procs * threads),
               '--num-threads=%d' % threads,
               '--num-smt=%d' % smt,
               ifthen(ppn, '--ppn=%d' % ppn, ''),
               ifthen(ppn_used, '--ppn-used=%d' % ppn_used, ''),
               '--walltime=0:59:0']
        submit.write(' '.join(cmd) + '\n')
        simulations.write(sim + '\n')



# All power-of-two multiples of min_cores, and all divisors of max_cores
cores_list = ([min_cores * 2**p for p in range(0, 20)] +
              [n for n in range(min_cores, max_cores+1) if max_cores % n == 0])
cores_list = [n for n in cores_list if n <= max_cores]
cores_list = unique(cores_list)

try:
    os.makedirs(path)
except OSError:
    if not os.path.exists(path):
        raise
    pass
global buildlines
buildlines = []
submit = open('%s/submit.sh' % path, 'w')
# os.fchmod(submit.fileno(), 0755)
os.chmod('%s/submit.sh' % path, 0755)
submit.write('''\
#! /bin/bash

# This script contains Simulation Factory commands to submit all benchmarks

''')
simulations = open('%s/simulations.txt' % path, 'w')
simulations.write('BEGIN_SIMULATIONS\n')

# Requested from the queuing system:
#    nodes:    nodes (as many as are used)
#    ppn:      cores per node
# Actually used:
#    nodes:    nodes
#    cores:    cores
#    ppn_used: cores per node
# Processes and threads:
#    procs:    OpenMP threads
#    threads:  OpenMP threads per process
#    smt:      threads per core
# Relation:
#    cores * smt = procs

count = 0

# Loop over physics parameters
for hydro in hydro_list:
    for multipatch in multipatch_list:
        if hydro and multipatch: continue
        for levels in levels_list:
            for fdorder in fdorder_list:
                
                # Loop over system configurations
                for cores in cores_list:
                    for ppn_used in ppn_used_list:
                        if ppn_used:
                            if cores % ppn_used != 0: continue
                            nodes = cores / ppn_used
                            if nodes < min_nodes or nodes > max_nodes: continue
                            else:
                                nodes = 1
                                for threads in threads_list:
                                    for smt in smt_list:
                                        if threads % smt != 0: continue
                                        if cores * smt % threads != 0: continue
                                        procs = cores * smt / threads
                                        create_benchmark(hydro, multipatch,
                                                         levels, fdorder,
                                                         nodes, cores, ppn_used,
                                                         procs, threads, smt)
                                        count += 1

submit.close()
simulations.write('END_SIMULATIONS\n')
simulations.close()

build = open('%s/build.sh' % path, 'w')
# os.fchmod(build.fileno(), 0755)
os.chmod('%s/build.sh' % path, 0755)
build.write('''\
#! /bin/bash

# This script contains Simulation Factory commands to build executables
# for all benchmarks

''')
for buildline in unique(buildlines):
    build.write(buildline + '\n')
build.close()

collect = open('%s/collect.sh' % path, 'w')
# os.fchmod(collect.fileno(), 0755)
os.chmod('%s/collect.sh' % path, 0755)
collect.write('''\
#! /bin/bash

# This script collects results from all benchmarks

echo '# column 1: total grid point updates/second'
echo '# column 2: seconds/grid point update (smaller is better)'
echo '# column 3: flop/grid point update (smaller is better)'
echo '# column 4: simulation name'
path=$(dirname $0)
simcmd="$path/../simfactory/bin/sim"
machine=$($simcmd whoami 2>/dev/null | grep -v 'Warning: Unable to determine CACTUS_PATH' | sed -e 's/Current machine: *//')
user=$($simcmd print-mdb-entry $machine 2>/dev/null | grep '^ *user[ =]' | sed -e 's/^ *user *= *//;s/ *#.*$//')
basedir=$($simcmd print-mdb-entry $machine 2>/dev/null | grep '^ *basedir[ =]' | sed -e 's/^ *basedir *= *//;s/ *#.*$//')
basedir=$(echo " $basedir" | sed -e 's/^ //;s/@USER@/'"$user"'/g')
echo BEGIN_TIMINGS
flop_sec='%g'
for sim in $(grep -v _SIMULATIONS $path/simulations.txt); do
    cores=$(echo $sim | sed -e 's/^.*-cores0*\\([0-9]*\\).*$/\\1/')
    awk '/^[^#]/ {
        flop_sec = '$flop_sec';
        cores = '$cores';
        total_gpu_sec = $26;
        sec_gpu = total_gpu_sec==0 ? 0 : cores / total_gpu_sec;
        flop_gpu = flop_sec * sec_gpu;
        print total_gpu_sec, sec_gpu, flop_gpu, "'$sim'";
    }' "$basedir"/$sim/output-0000/*/carpet-timing..asc 2>/dev/null |
    tail -n 1
done
echo END_TIMINGS
''' % flop_per_second)
collect.close()

print
print
print
print 'Created %d parameter files.' % count
