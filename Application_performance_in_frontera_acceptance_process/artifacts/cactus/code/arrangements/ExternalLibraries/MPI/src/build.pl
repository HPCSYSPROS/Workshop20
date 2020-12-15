#!/usr/bin/perl
use Carp;
use strict;
use FileHandle;
use Cwd;
$/ = undef;

# Set locations
my $THORN = "MPI";
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
print "Installing MPI into ${INSTALL_DIR}\n";
my $DONE_FILE = "$ENV{SCRATCH_BUILD}/done/${THORN}";
my $mpi_dir = ${INSTALL_DIR};

chdir($ENV{SCRATCH_BUILD});

# Set up environment
# Disable ccache: remove "ccache" and all options that follow
$ENV{CC} =~ s/^\s*ccache\s+(-\S*\s+)*//;
$ENV{CXX} =~ s/^\s*ccache\s+(-\S*\s+)*//;
if ($ENV{F90} eq 'none') {
    print "No Fortran 90 compiler available. Building MPI library without Fortran support.\n";
    $ENV{FC} = undef;
    $ENV{FCFLAGS} = undef;
} else {
    $ENV{FC} = $ENV{F90};
    $ENV{FCFLAGS} = $ENV{F90FLAGS};
}
$ENV{LIBS} = undef;
$ENV{RPATH} = undef;
if ($ENV{ARFLAGS} =~ /64/) {
    $ENV{OBJECT_MODE} = "64";
}

print "MPI: Preparing directory structure...\n";
mkdir("build");
mkdir("external");
mkdir("done");
system("rm -rf ${BUILD_DIR} ${INSTALL_DIR}") == 0 or die;
mkdir(${BUILD_DIR});
mkdir(${INSTALL_DIR});
error("${INSTALL_DIR} does not exist.",6) unless(-e ${INSTALL_DIR});
error("${INSTALL_DIR} is not a directory.",6) unless(-d ${INSTALL_DIR});
error("${INSTALL_DIR} is not readable.",7) unless(-r ${INSTALL_DIR});
error("${INSTALL_DIR} is not writeable.",8) unless(-w ${INSTALL_DIR});
error("${INSTALL_DIR} is not executable.",8) unless(-x ${INSTALL_DIR});

print "MPI: Unpacking archive...\n";
chdir(${BUILD_DIR});
system("$ENV{TAR} xzf ${SRCDIR}/../dist/${NAME}.tar.gz") == 0 or die;

print "MPI: Configuring...\n";
chdir(${NAME});
my $hwloc_opts = '';
if ($ENV{HWLOC_DIR} ne '') {
    $hwloc_opts = "--with-hwloc='$ENV{HWLOC_DIR}'";
}
# Cannot have a memory manager with a static library on some systems
# (e.g. Linux); see
# <http://www.open-mpi.org/faq/?category=mpi-apps#static-mpi-apps>
system("./configure --prefix='$mpi_dir' $hwloc_opts --without-memory-manager --enable-shared=no --enable-static=yes") == 0 or die;

print "MPI: Building...\n";
system("$ENV{MAKE}") == 0 or die;

print "MPI: Installing...\n";
system("$ENV{MAKE} install") == 0 or die;
chdir($ENV{SCRATCH_BUILD});

print "MPI: Cleaning up...\n";
system("rm -rf ${BUILD_DIR}") == 0 or die;

system("date > ${DONE_FILE}") == 0 or die;
print "MPI: Done.\n";
