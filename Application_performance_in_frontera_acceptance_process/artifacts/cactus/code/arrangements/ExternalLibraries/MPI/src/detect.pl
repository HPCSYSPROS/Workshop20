#!/usr/bin/perl
use Carp;
use strict;
use FileHandle;
use Cwd;
$/ = undef;

################################################################################
# Prepare
################################################################################

# Set up shell
my $verbose = 0;
$verbose = 1 if $ENV{VERBOSE} =~ /^yes$/i;

################################################################################
# Check for old mechanism
################################################################################

if (defined($ENV{MPI})) {
    error("Setting the option \"MPI\" is incompatible with the MPI thorn. Please remove the option MPI = $ENV{MPI}.", 2);
}

################################################################################
# Determine whether to build and/or search
################################################################################

my $info = undef;
my $mpi_info_set = 0;
my $mpi_dir = undef;
my $mpi_cmd = undef;
my $mpi_search = 1;
my $mpi_build = 0;
my $mpi_manual = 0;

# Note:
# - look for "mpiCC" last, since this is equivalent to mpicc on
#   case-insensitive file systems
# - don't look for "mpicc", since linking with a C compiler won't work
my @mpicxx_names = ("mpic++", "mpicxx", "mpicxx-openmpi-mp", "mpiCC");

if (!is_set("MPI_DIR")) {
    message("MPI selected, but MPI_DIR is not set. Computing settings...");
    $mpi_build = 1;
    $mpi_search = 1;
} elsif ($ENV{MPI_DIR} eq "NO_BUILD") {
    $mpi_dir = $ENV{MPI_DIR};
    $mpi_build = 0;
    $mpi_search = 1;
} elsif ($ENV{MPI_DIR} eq "BUILD") {
    $mpi_build = 1;
    $mpi_search = 0;
} elsif ($ENV{MPI_DIR} eq "NONE") {
    $mpi_build = 0;
    $mpi_search = 0;
    $mpi_info_set = 1;
    $mpi_dir = '';
    $info = '';
} else {
    if (!-d $ENV{MPI_DIR}) {
        message("MPI_DIR is set to a directory that does not exist (MPI_DIR = $ENV{MPI_DIR}); continuing anyway");
    }
    $mpi_dir = $ENV{MPI_DIR};
    $mpi_build = 0;
    $mpi_search = 0;
    if (is_set("MPI_INC_DIRS") or is_set("MPI_LIB_DIRS") or is_set("MPI_LIBS"))
    {
        # If some of the MPI variables are set, this is a completely
        # manual configuration.
        $mpi_manual = 1;
    } else {
        # If none of the MPI variables are set, check for the compiler
        # wrapper under MPI_DIR
        $mpi_manual = 0;
        for my $name (@mpicxx_names) {
            my $full_name = "$ENV{MPI_DIR}/bin/$name";
            if (-x $full_name) {
                $mpi_cmd = $full_name;
                last;
            }
        }
        if (defined($mpi_cmd)) {
            message("Found MPI compiler wrapper at $mpi_cmd!");
            mpi_get_info();
        } else {
            message("No MPI compiler wrapper found beneath MPI_DIR (MPI_DIR = $ENV{MPI_DIR})");
        }
    }
}

################################################################################
# Search
################################################################################
if ($mpi_search and !defined($mpi_cmd)) {
    for my $name (@mpicxx_names) {
        last if defined($mpi_cmd);
        $mpi_cmd = which($name);
    }
    if (defined($mpi_cmd)) {
        $mpi_dir = $mpi_cmd;
        $mpi_dir =~ s{/mpi(c\+\+|CC|cc|cxx)[^/]*$}{};
        $mpi_dir =~ s{/bin$}{};
        message("Found MPI compiler wrapper at $mpi_cmd!");
        mpi_get_info();
    }
}

my $THORN = "MPI";

################################################################################
# Build
################################################################################

if ($mpi_build and !$mpi_info_set) {
    # Check for required tools. Do this here so that we don't require
    # them when using the system library.
    unless(defined($ENV{TAR}) and $ENV{TAR} =~ /\S/ and -x which($ENV{TAR})) {
        error(
            "ENV{TAR} = $ENV{TAR}\n" .
            "Could not find tar command. Please make sure that (GNU) tar is present,\n" .
            "and that the TAR variable is set to its location.\n", 3);
    }
    unless(defined($ENV{PATCH}) and $ENV{PATCH} =~ /\S/ and
           -x which($ENV{PATCH}))
    {
        error(
            "Could not find patch command. Please make sure that (GNU) patch is present,\n" .
            "and that the PATCH variable is set to its location.\n", 4);
    }

    # Set locations
    my $NAME = "openmpi-1.10.1";
    my $INSTALL_DIR = undef;
    my $BUILD_DIR = undef;
    my $SRCDIR = $0;
    $SRCDIR =~ s{(.*)/.*}{$1};
    ${BUILD_DIR} = "$ENV{SCRATCH_BUILD}/build/${THORN}";
    if (defined($ENV{MPI_INSTALL_DIR}) and $ENV{MPI_INSTALL_DIR} =~ /\S/) {
        $INSTALL_DIR = "$ENV{MPI_INSTALL_DIR}/${THORN}";
    } else {
        $INSTALL_DIR = "$ENV{SCRATCH_BUILD}/external/${THORN}";
    }
    message("Installing MPI into ${INSTALL_DIR}");
    $mpi_dir = ${INSTALL_DIR};

    $mpi_manual = 1;
    $ENV{MPI_BUILD} = '1';
    $ENV{MPI_DIR} = $mpi_dir;
    $ENV{MPI_INC_DIRS} = "$mpi_dir/include";
    $ENV{MPI_LIB_DIRS} = "$mpi_dir/lib";
    my $mpi_fortranlibs = '';
    if ($ENV{F90} ne 'none') {
        $mpi_fortranlibs = "mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh";
    }
    $ENV{MPI_LIBS} = "$mpi_fortranlibs mpi_cxx mpi open-rte open-pal";
} else {
    $ENV{MPI_BUILD} = '';
    my $DONE_FILE = "$ENV{SCRATCH_BUILD}/done/${THORN}";
    if (! -e $DONE_FILE) {
        mkdir("$ENV{SCRATCH_BUILD}/done");
        system("date > ${DONE_FILE}") == 0 or die;
    }
}

################################################################################
# Configure MPI options
################################################################################

if ($mpi_info_set) {
    my @incdirs = ();
    my @libdirs = ();
    my @libs = ();
    while($info =~ /\s-I\s*(\S+)/g) {
        push @incdirs, $1;
    }
    while($info =~ /\s-L\s*(\S+)/g) {
        push @libdirs, $1;
    }
    while($info =~ /\s-l(\S+)/g) {
        push @libs, $1;
    }

    $ENV{MPI_DIR} = $mpi_dir;
    $ENV{MPI_INC_DIRS} = join(" ",@incdirs);
    $ENV{MPI_LIB_DIRS} = join(" ",@libdirs);
    $ENV{MPI_LIBS} = join(" ",@libs);

    message("Successfully configured MPI.");
} elsif ($mpi_manual) {
    my @incdirs = ();
    my @libdirs = ();
    my @libs = ();
    if (is_set("MPI_INC_DIRS")) {
        push @incdirs, $ENV{MPI_INC_DIRS};
    } else {
        push @incdirs, "$ENV{MPI_DIR}/include";
    }
    if (is_set("MPI_LIB_DIRS")) {
        push @libdirs, $ENV{MPI_LIB_DIRS};
    } else {
        push @libdirs, "$ENV{MPI_DIR}/lib64";
        push @libdirs, "$ENV{MPI_DIR}/lib";
    }
    if (is_set("MPI_LIBS")) {
        push @libs, $ENV{MPI_LIBS};
    } else {
        # do nothing
    }

    $ENV{MPI_INC_DIRS} = join(" ",@incdirs);
    $ENV{MPI_LIB_DIRS} = join(" ",@libdirs);
    $ENV{MPI_LIBS} = join(" ",@libs);

    message("MPI was manually configured.");
} else {
    error("MPI could not be configured: neither automatic nor manual configuration succeeded",5);
}

################################################################################
# Configure Cactus
################################################################################

# Strip standard paths
$ENV{MPI_INC_DIRS} = strip_inc_dirs($ENV{MPI_INC_DIRS});
$ENV{MPI_LIB_DIRS} = strip_lib_dirs($ENV{MPI_LIB_DIRS});

# Pass configuration options to build script
print "BEGIN MAKE_DEFINITION\n";
print "MPI_BUILD       = $ENV{MPI_BUILD}\n";
print "MPI_INSTALL_DIR = $ENV{MPI_INSTALL_DIR}\n";
print "HWLOC_DIR       = $ENV{HWLOC_DIR}\n";
print "END MAKE_DEFINITION\n";

# Pass options to Cactus

print "BEGIN DEFINE\n";
print "CCTK_MPI 1\n";
print "END DEFINE\n";

print "BEGIN MAKE_DEFINITION\n";
print "CCTK_MPI     = 1\n";
print "MPI_DIR      = $ENV{MPI_DIR}\n";
print "MPI_INC_DIRS = $ENV{MPI_INC_DIRS}\n";
print "MPI_LIB_DIRS = $ENV{MPI_LIB_DIRS}\n";
print "MPI_LIBS     = $ENV{MPI_LIBS}\n";
print "END MAKE_DEFINITION\n";

print "INCLUDE_DIRECTORY \$(MPI_INC_DIRS)\n";
print "INCLUDE_DIRECTORY_FORTRAN \$(MPI_LIB_DIRS)\n";
print "LIBRARY_DIRECTORY \$(MPI_LIB_DIRS)\n";
print "LIBRARY           \$(MPI_LIBS)\n";

################################################################################
# Functions
################################################################################

sub which {
    my $cmd = shift;
    for my $path (split(/:/,$ENV{PATH})) {
        my $full_cmd = "$path/$cmd";
        if (-x $full_cmd) {
            return $full_cmd;
        }
    }
    return undef;
}

sub mpi_get_info {
    my $fd = new FileHandle;
    open($fd,"$mpi_cmd -compile_info 2>/dev/null|");
    $info = <$fd>;
    close($fd);
    if ($info eq "") {
        open($fd,"$mpi_cmd --showme 2>/dev/null|");
        $info = <$fd>;
        close($fd);
    }
    if ($info eq "") {
        # The command, mpicc, is quite often a shell script.
        # Run it with -x to trace, and find the compile command.
        open($fd,"sh -x $mpi_cmd /dev/null 2>/dev/null|");
        my $contents = <$fd>;
        if ($contents =~ /\w+cc.*-I\/.*-lmpi.*/) {
            $info = $&;
        }
    }
    if ($info =~ /^\s*$/) {
        $mpi_info_set = 0;
    } else {
        $mpi_info_set = 1;
    }
}

sub message {
    my $msg = shift;
    $msg =~ s/\n$//;
    print "BEGIN MESSAGE\n";
    print "$msg\n";
    print "END MESSAGE\n";
}
sub error {
    my ($msg,$errno) = @_;
    $msg =~ s/\n$//;
    print "BEGIN ERROR\n";
    print "$msg\n";
    print "END ERROR\n";
    exit $errno;
}

sub is_set {
    my $var = shift;
    return (defined($ENV{$var}) and !($ENV{$var} =~ /^\s*$/));
}

sub strip_inc_dirs {
    my $dirlist = shift;
    my @dirs = split / /, $dirlist;
    map { s{//}{/}g } @dirs;
    @dirs = grep { !m{^/(usr/(local/)?)?include/?$} } @dirs;
    return join ' ', @dirs;
}
sub strip_lib_dirs {
    my $dirlist = shift;
    my @dirs = split / /, $dirlist;
    map { s{//}{/}g } @dirs;
    @dirs = grep { !m{^/(usr/(local/)?)?lib(64?)/?$} } @dirs;
    return join ' ', @dirs;
}
