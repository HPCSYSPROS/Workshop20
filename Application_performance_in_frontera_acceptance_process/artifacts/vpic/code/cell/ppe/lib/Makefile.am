################################################################################
## VPIC library Makefile.am
## THIS FILE IS GENERATED BY config/bootstrap DO NOT EDIT
################################################################################

lib_LTLIBRARIES = libvpic.la

INCLUDES = -I@HOST_ACDIR@/src/include \
   @EXTRA_CPPFLAGS@

libvpic_includedir = ${includedir}

libvpic_include_HEADERS = \
	@HOST_ACDIR@/src/boundary/boundary.h \
	@HOST_ACDIR@/src/emitter/emitter.h \
	@HOST_ACDIR@/src/field_advance/field_advance.h \
	@HOST_ACDIR@/src/field_advance/standard/sfa_private.h \
	@HOST_ACDIR@/src/field_advance/standard/v4/sfa_v4_private.h \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/v4/vfa_v4_private.h \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_private.h \
	@HOST_ACDIR@/src/grid/grid.h \
	@HOST_ACDIR@/src/material/material.h \
	@HOST_ACDIR@/src/sf_interface/sf_interface.h \
	@HOST_ACDIR@/src/sf_interface/sf_interface_private.h \
	@HOST_ACDIR@/src/species_advance/species_advance.h \
	@HOST_ACDIR@/src/species_advance/standard/spa.h \
	@HOST_ACDIR@/src/species_advance/standard/spa_private.h \
	@HOST_ACDIR@/src/util/io/FileIO.hxx \
	@HOST_ACDIR@/src/util/io/FileIOData.hxx \
	@HOST_ACDIR@/src/util/io/FileUtils.hxx \
	@HOST_ACDIR@/src/util/io/P2PIOPolicy.hxx \
	@HOST_ACDIR@/src/util/io/P2PUtilsPolicy.hxx \
	@HOST_ACDIR@/src/util/io/StandardIOPolicy.hxx \
	@HOST_ACDIR@/src/util/io/StandardUtilsPolicy.hxx \
	@HOST_ACDIR@/src/util/mp/dmp/DMPPolicy.hxx \
	@HOST_ACDIR@/src/util/mp/dmp/mp_dmp.h \
	@HOST_ACDIR@/src/util/mp/dmp/mp_t.h \
	@HOST_ACDIR@/src/util/mp/mp.h \
	@HOST_ACDIR@/src/util/mp/mp.hxx \
	@HOST_ACDIR@/src/util/mp/mp_handle.h \
	@HOST_ACDIR@/src/util/mp/MPWrapper.hxx \
	@HOST_ACDIR@/src/util/mp/relay/RelayPolicy.hxx \
	@HOST_ACDIR@/src/util/mtrand/drandn_table.h \
	@HOST_ACDIR@/src/util/mtrand/frandn_table.h \
	@HOST_ACDIR@/src/util/mtrand/mtrand.h \
	@HOST_ACDIR@/src/util/mtrand/mtrand_conv.h \
	@HOST_ACDIR@/src/util/pipelines/pipelines.h \
	@HOST_ACDIR@/src/util/relay/CMPolicyDaCS.hxx \
	@HOST_ACDIR@/src/util/relay/CMPolicyMPIDaCS.hxx \
	@HOST_ACDIR@/src/util/relay/CMPolicyMultipleContextMPI.hxx \
	@HOST_ACDIR@/src/util/relay/ConnectionManager.hxx \
	@HOST_ACDIR@/src/util/relay/dmp/DMPConnection.hxx \
	@HOST_ACDIR@/src/util/relay/dmp/DMPPolicyMPI.hxx \
	@HOST_ACDIR@/src/util/relay/MPData.hxx \
	@HOST_ACDIR@/src/util/relay/MPRelay.hxx \
	@HOST_ACDIR@/src/util/relay/p2p/DaCSUtils.h \
	@HOST_ACDIR@/src/util/relay/p2p/P2PConnection.hxx \
	@HOST_ACDIR@/src/util/relay/p2p/P2PPolicyDaCS.hxx \
	@HOST_ACDIR@/src/util/relay/p2p/P2PPolicyMPI.hxx \
	@HOST_ACDIR@/src/util/relay/p2p/P2PTag.hxx \
	@HOST_ACDIR@/src/util/relay/Type2DaCSSwapType.hxx \
	@HOST_ACDIR@/src/util/relay/Type2MPIType.hxx \
	@HOST_ACDIR@/src/util/swap.h \
	@HOST_ACDIR@/src/util/util.h \
	@HOST_ACDIR@/src/util/util_base.h \
	@HOST_ACDIR@/src/util/BitField.hxx \
	@HOST_ACDIR@/src/util/v4/v4.h \
	@HOST_ACDIR@/src/util/v4/v4_altivec.hxx \
	@HOST_ACDIR@/src/util/v4/v4_portable.hxx \
	@HOST_ACDIR@/src/util/v4/v4_spu.hxx \
	@HOST_ACDIR@/src/util/v4/v4_sse.hxx \
	@HOST_ACDIR@/src/util/v4/v4c_spu.h \
	@HOST_ACDIR@/src/vpic/vpic.hxx \
	@HOST_ACDIR@/src/vpic/dumpmacros.h \
	@HOST_ACDIR@/src/util/CheckSum.hxx

libvpic_la_SOURCES = \
	@HOST_ACDIR@/src/boundary/absorb_tally.c \
	@HOST_ACDIR@/src/boundary/link.c \
	@HOST_ACDIR@/src/boundary/maxwellian_reflux.c \
	@HOST_ACDIR@/src/emitter/ccube.c \
	@HOST_ACDIR@/src/emitter/child-langmuir.c \
	@HOST_ACDIR@/src/emitter/emitter.c \
	@HOST_ACDIR@/src/emitter/ivory.c \
	@HOST_ACDIR@/src/field_advance/field_advance.c \
	@HOST_ACDIR@/src/field_advance/standard/advance_b.c \
	@HOST_ACDIR@/src/field_advance/standard/advance_e.c \
	@HOST_ACDIR@/src/field_advance/standard/clean_div_b.c \
	@HOST_ACDIR@/src/field_advance/standard/clean_div_e.c \
	@HOST_ACDIR@/src/field_advance/standard/compute_curl_b.c \
	@HOST_ACDIR@/src/field_advance/standard/compute_div_b_err.c \
	@HOST_ACDIR@/src/field_advance/standard/compute_div_e_err.c \
	@HOST_ACDIR@/src/field_advance/standard/compute_rhob.c \
	@HOST_ACDIR@/src/field_advance/standard/compute_rms_div_b_err.c \
	@HOST_ACDIR@/src/field_advance/standard/compute_rms_div_e_err.c \
	@HOST_ACDIR@/src/field_advance/standard/energy_f.c \
	@HOST_ACDIR@/src/field_advance/standard/local.c \
	@HOST_ACDIR@/src/field_advance/standard/remote.c \
	@HOST_ACDIR@/src/field_advance/standard/sfa.c \
	@HOST_ACDIR@/src/field_advance/standard/v4/sfa_v4.c \
	@HOST_ACDIR@/src/field_advance/standard/v4/v4_advance_b.cxx \
	@HOST_ACDIR@/src/field_advance/standard/v4/v4_advance_e.cxx \
	@HOST_ACDIR@/src/field_advance/standard/v4/v4_clean_div_b.cxx \
	@HOST_ACDIR@/src/field_advance/standard/v4/v4_compute_curl_b.cxx \
	@HOST_ACDIR@/src/field_advance/standard/v4/v4_compute_div_b_err.cxx \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/v4/vfa_v4.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/v4/vfa_v4_advance_e.cxx \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_advance_e.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_clean_div_e.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_compute_curl_b.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_compute_div_e_err.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_compute_rhob.c \
	@HOST_ACDIR@/src/field_advance/standard/vacuum/vfa_energy_f.c \
	@HOST_ACDIR@/src/grid/add_boundary.c \
	@HOST_ACDIR@/src/grid/distribute_voxels.c \
	@HOST_ACDIR@/src/grid/grid_comm.c \
	@HOST_ACDIR@/src/grid/grid_structors.c \
	@HOST_ACDIR@/src/grid/ops.c \
	@HOST_ACDIR@/src/grid/partition.c \
	@HOST_ACDIR@/src/material/material.c \
	@HOST_ACDIR@/src/sf_interface/clear_accumulators.c \
	@HOST_ACDIR@/src/sf_interface/hydro.c \
	@HOST_ACDIR@/src/sf_interface/load_interpolator.cxx \
	@HOST_ACDIR@/src/sf_interface/reduce_accumulators.cxx \
	@HOST_ACDIR@/src/sf_interface/sf_interface.c \
	@HOST_ACDIR@/src/sf_interface/unload_accumulator.cxx \
	@HOST_ACDIR@/src/species_advance/species_advance.c \
	@HOST_ACDIR@/src/species_advance/standard/advance_p.cxx \
	@HOST_ACDIR@/src/species_advance/standard/boundary_p.c \
	@HOST_ACDIR@/src/species_advance/standard/center_p.cxx \
	@HOST_ACDIR@/src/species_advance/standard/energy_p.cxx \
	@HOST_ACDIR@/src/species_advance/standard/hydro_p.c \
	@HOST_ACDIR@/src/species_advance/standard/move_p.c \
	@HOST_ACDIR@/src/species_advance/standard/rho_p.c \
	@HOST_ACDIR@/src/species_advance/standard/sort_p.c \
	@HOST_ACDIR@/src/species_advance/standard/uncenter_p.cxx \
	@HOST_ACDIR@/src/util/mp/dmp/mp_dmp.c \
	@HOST_ACDIR@/src/util/mp/mp.cxx \
	@HOST_ACDIR@/src/util/mtrand/mtrand.c \
	@HOST_ACDIR@/src/util/pipelines/pipelines_serial.c \
	@HOST_ACDIR@/src/util/pipelines/pipelines_thread.c \
	@HOST_ACDIR@/src/util/util.c \
	@HOST_ACDIR@/src/vpic/advance.cxx \
	@HOST_ACDIR@/src/vpic/dump.cxx \
	@HOST_ACDIR@/src/vpic/finalize.cxx \
	@HOST_ACDIR@/src/vpic/initialize.cxx \
	@HOST_ACDIR@/src/vpic/misc.cxx \
	@HOST_ACDIR@/src/vpic/vpic.cxx

if USE_CELL_SPUS
libvpic_la_SOURCES += \
	@HOST_ACDIR@/src/util/pipelines/pipelines_spu.c
endif

libvpic_la_LDFLAGS = @VERSION_INFORMATION@ \
	@EXTRA_LDFLAGS@
libvpic_la_DEPENDENCIES = @EXTRA_O_FILES@
libvpic_la_LIBADD = @EXTRA_LO_FILES@ @EXTRA_LIBS@

include @top_builddir@/config/rules.mk
