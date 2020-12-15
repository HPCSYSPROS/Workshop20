# makefile for DNS2d.x on multiple platforms (DD07)
# Note: on DS, use "gmake" instead of "make"
# Note: on non-IBM machines use FFTW instead of ESSL
SHELL=/bin/bash


#MACHINE = hopper

#MACHINE = craycce
# suggested by D. Whitaker (Cray)
# run 'module swap PrgEnv-pgi PrgEnv-cray' first
#MACHINE = intel
#MACHINE = intel.edison
#MACHINE = intel.comet
#MACHINE = BW

#MACHINE = mira
#MACHINE = mira.essl
#MACHINE = fermi

MACHINE = intel

# Allow command line to add FLAGS
CPPFLAGS=$(ARG)

#====================================================================
# hybrid MPI and OpenMP
# CAF replaces MPI_AlltoAll and MPI_AlltoAllV
# CAF_TEST compares CAF call results to MPI_AlltoALL(V) call results
# For CAF executables, type 'module load craype-hugepages8M' first
#-------------
CPPFLAGS += -P -DOPENMP
# CPPFLAGS += -P -DCAF
#CPPFLAGS += -P -DMPI -DCAF
#CPPFLAGS += -P -DCAF -DOPENMP
#CPPFLAGS += -P -DMPI -DCAF_TEST #-DOPENMP 
#CPPFLAGS += -P -DMPI -DCAFI #-DOPENMP
#CPPFLAGS += -P -DMPI -DCAFI_TEST #-DOPENMP 

#CPPFLAGS += -DAUG1912
CPPFLAGS += -DKXCOMM2

CPPFLAGS += -DCLEAN
#====================================================================
# Fresh start 
#-------------
CPPFLAGS += -DFRESH 
#===================================================================
#Option for getting full convective terms (for extracting the pressure)
#CPPFLAGS += -P -DCVC_PRESS
#===================================================================
# Option for getting vorticity in physicsl space (for Track 1 benchmark)
#CPPFLAGS += -P -DVORPHY

CPPFLAGS += -P -DBOXPHY

#====================================================================
# 4th order R-K (for homogeneous code): RKFOUR and
CPPFLAGS += -DRKFOUR #-DJUN1411 #-DMAY28 -DMAY30
# CPPFLAGS += -DRKFOUR
# MAY28: pretend to do an RK2 under RKFOUR

#====================================================================
# Random number generator
CPPFLAGS += -DRANDOM_SEED             
#====================================================================
# Forcing (choose one)
#---------------------
CPPFLAGS += -DEPFOR #                # Eswaran & Pope 1988 forcing
#CPPFLAGS += -DEPFOR -DSBP_2010 #                # Eswaran & Pope 1988 forcing
#CPPFLAGS += -DFEK_FORC       # freezes energy in forced wavenumber shell
#                                    #(use kforce=0 for EPFOR in input)
CPPFLAGS += -DDOUBLE_PREC
CPPFLAGS += -DREADR
#====================================================================
# FFT (choose one)
#-----------------
#CPPFLAGS += -DESSL #      # uses ESSL
CPPFLAGS += -DFFTW #       # uses FFTWs
#
#====================================================================
# optimiz. for y transforms
#---------------------------
#CPPFLAGS += -DONEYTRAN # creates plans and executes transforms in y in one shot (FFTW)
#
#====================================================================
# Homogeneous / Inhomogeneous Flow
#---------
CPPFLAGS += -DHOMOGENEOUS     # homogeneous, isotropic flow
#CPPFLAGS += -DCHANNEL       # channel flow
#CPPFLAGS += -DBL            # boundary layer
#
# _====================================================================
# Scalars?
#---------
#CPPFLAGS += -DNOSCALAR 
#
#====================================================================
# optionally set prefered blocking factors for loop blocking in memory transposes
#---------------------------
#CPPFLAGS += -DNBL1_X=8 -DNBL1_Z=64 -DNBL2_Z=4 -DNBL2_Y=4
#
#====================================================================
# Optionally set period for task remapping routine (to reduce load imbalance)
#
CPPFLAGS += -DBAL_PERIOD=4
#CPPFLAGS += -DBAL_PERIOD=8
#
#====================================================================
# Debugging
#CPPFLAGS += -DDEBUG
#
#====================================================================
# Profiling
#----------
#CPPFLAGS += -DTIMERS #    # gets some timings on alltoall's and output in file "timers"
#CPPFLAGS += -DTIMERS_IO #  # adds I/O timings to the file "timers"
#
#====================================================================
# Lagrangian stuff
#-----------------
#CPPFLAGS += -DLAG -DLAGSC2 -DLGRAD -DSPXYZM -DMOL # -DCF_LAG
#Use next line on BW with CAF
#CPPFLAGS += -DLAG_DOUBLE #

#CPPFLAGS += -DLAG -DLAGSC2 -DLGRAD -DSPXYZM -DMOL -DCF_LAG
#CPPFLAGS += -DMOL2
#CPPFLAGS += -DPARTIN_AR

#
#====================================================================
# Some optimizations regarding communications
#--------------------------------------------
CPPFLAGS += -DREVERSE_DIMS #
CPPFLAGS += -DUSE_EVEN #   # uses ALLTOALL instead of ALLTOALLV in *comm*.f subroutine 
# (it shows better performance on some machines, e.g. Cray XT3/XT4)
#CPPFLAGS += -DMULTIVAR     # transpose multiple variables together to decrease number of messages 
#				  
#====================================================================
# outputs 1d spectra and correlations
#------------------------------------
CPPFLAGS += -DOUT1D
#
#====================================================================
# each proc. generates a individual stdout: stdout.taskid
#------------------------------------------------------
#CPPFLAGS += -DIND_STDOUT  #
#
#====================================================================
# Buoyancy
#---------
#CPPFLAGS += -DLINDBORG -DBOUSS -DBUOY #-DBIF -DBUOY2
#CPPFLAGS += -DBOUSS #-DBUOY #-DBIF -DBUOY2
#
#====================================================================
# Checkpointing in old (pencils) format
CPPFLAGS += -DINPUT_PEN -DOUTPUT_PEN
# Used to turn on/off new HDF5 features
#CPPFLAGS += -DNEW_HDF5
#--------------------------------------------------------------------

# Output spectral quantities in shells of thickness "dk" instead
# of unity (for runs with say 4pi^3 or 8pi^3 domains)
CPPFLAGS += -DSHELL_DK

# Allow for restart data files distributed among multiple sub-directories
CPPFLAGS += -DSUBDIR

# Update velocity field only once every several steps
# (this may be useful for simulations of scalars of very low Sc)
CPPFLAGS += -DLVSTEP

#CPPFLAGS += -DFROZEN

#CPPFLAGS += -DROTD
#CPPFLAGS += -DMHD

#====================================================================
# The following directives were implemented by Matthew Clay
# Code to set the initial spectrum to the Pope model spectrum.
CPPFLAGS += -DMODEL_SPECTRUM
# #CPPFLAGS += -DSPECTRAL_TRUNCATION
# #CPPFLAGS += -DBETA_UPDATE
# #CPPFLAGS += -DROGALLO_EXACT
# #CPPFLAGS += -DSPATIAL_STRAIN -DERF_VELOCITY
# #CPPFLAGS += -DTIME_STRAIN
# #CPPFLAGS += -DNODE_DENSITY
# #--------------------------------------------------------------------



OTHERLIBS = 
CFLAGS =
LFLAGS =

#====================================================================
#--------------------------------------------------------------------
#            * setting up hook structure *
#
# VPATH specifies a list of directories that make should search
#
#--------------------------------------------------------------------
#====================================================================
#save paths 
HOOK=hooks
NHOOK=nullhooks

# inhomogeneous capability -- can only choose one
ifneq (,$(findstring CHANNEL,$(CPPFLAGS)))
	VPATH +=inhomogeneous/channel 
else
ifneq (,$(findstring BL,$(CPPFLAGS)))
	VPATH +=inhomogeneous/bl
endif
endif

# no inhomogeneous chosen -- make homogeneous
ifneq (,$(findstring HOMOGENEOUS,$(CPPFLAGS)))
	VPATH +=homogeneous
endif

# Lagrangian Particle Tracking Hook
ifneq (,$(findstring LAG,$(CPPFLAGS)))
#    	VPATH +=$(HOOK)/LAG 
#	VPATH +=$(HOOK)/LAG_PIPE
 	VPATH +=$(HOOK)/LAG_LOCAL
else
	VPATH +=$(NHOOK)/LAG
endif

# BUOY Hook (must come before BOUSS)
ifneq (,$(findstring BUOY,$(CPPFLAGS)))
	VPATH +=$(HOOK)/stratification/BUOY
endif
# Bouss Hook 
ifneq (,$(findstring BOUSS,$(CPPFLAGS)))
	VPATH +=$(HOOK)/stratification/BOUSS 
endif

# Add generic paths
VPATH += $(HOOK)/initial_conditions $(NHOOK)/diagnostics
#VPATH += $(NHOOK)/initial_conditions $(NHOOK)/diagnostics $(NHOOK)/stratification
#VPATH += $(HOOK)/stratification

#######################################################################
ifeq ($(MACHINE),intel) ###############################################
  FF = mpif90
#  FFLAGS = -132 -no-vec
#  FFLAGS = -132 -O3 #-xW #-fast
#FFLAGS = -132 -O3 -convert big_endian -assume buffered_io -xCOMMON-AVX512 #-xhost
#FFLAGS90 = -O3 -convert big_endian -assume buffered_io -xCOMMON-AVX512 #-xhost
FFLAGS = -132 -O3 -convert big_endian -assume buffered_io -xCORE-AVX2 #-xhost
FFLAGS90 = -O3 -convert big_endian -assume buffered_io -xCORE-AVX2 #-xhost
#FFLAGS = -132 -O3 -convert big_endian -assume buffered_io -xCORE-AVX2 -no-vec -no-simd #-xhost
#FFLAGS90 = -O3 -convert big_endian -assume buffered_io -xCORE-AVX2 -no-vec -no-simd #-xhost
#FFLAGS = -132 -O3 -convert big_endian -assume buffered_io -xMIC-AVX512 #-xhost
#FFLAGS90 = -O3 -convert big_endian -assume buffered_io -xMIC-AVX512 #-xhost
#FFLAGS = -132 -O3 -convert big_endian -assume buffered_io -xCORE-AVX512 #-xhost
#FFLAGS90 = -O3 -convert big_endian -assume buffered_io -xCORE-AVX512 #-xhost
#-O3 -xW -fast
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
#    FFLAGS += -openmp
    FFLAGS += -qopenmp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -r8
  endif
  CPP = /lib/cpp -traditional
  # run "module load fftw3" before make
  FFT_LIB = -L$(TACC_FFTW3_LIB) -lfftw3 -lfftw3f
  FFT_INC = -I$(TACC_FFTW3_INC)
# On Stampede: must do a 
# swap intel intel/14.0.1.106
# module load fftw3 phdf5
  HDF_LIB = -Wl,-rpath,$(TACC_HDF5_LIB) -L$(TACC_HDF5_LIB) -lhdf5_fortran -lz
  HDF_INC = -I$(TACC_HDF5_INC)
  # run "module load mkl" before make
#  MKL_LIB = -L$(TACC_MKL_LIB) -lmkl -lguide -lpthread
#  MKL_INC = -I$(TACC_MKL_INC)
endif
#######################################################################
ifeq ($(MACHINE),intel.edison) ###############################################
  FF = ftn
FFLAGS = -132 -O3 -convert big_endian -assume buffered_io
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -openmp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -r8
  endif
  CPP = /lib/cpp -traditional
  # run "module load fftw3" before make
  FFT_LIB = -L$(FFTW_DIR) -lfftw3 -lfftw3f
  FFT_INC = -I$(FFTW_INC)
  # run "module load mkl" before make
endif
#######################################################################
ifeq ($(MACHINE),intel.comet) ###############################################
CPPFLAGS += -DALLTOALLV
  FF = mpif90
  FFLAGS = -132 -no-vec
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -openmp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -r8
  endif
  CPP = /lib/cpp -traditional
  #run "module load fftw" before make
  FFT_LIB = -L/opt/fftw/3.3.4/intel/mvapich2_ib/lib -lfftw3 -lfftw3f
  FFT_INC = -I/opt/fftw/3.3.4/intel/mvapich2_ib/include
  endif
#######################################################################
ifeq ($(MACHINE),craycce) ############################################
# to use Cray Fortran module,
# module swap PrgEnv-pgi PrgEnv-cray
  # run module load fftw/ before make
  FF = ftn #-C
   FFLAGS = -I. -Ofp3,scalar3,thread0 -N132 -h byteswapio -rm  #-Mframe
#   FFLAGS = -I. -Ofp3,scalar3,thread0 -N132 #-Rbc #-h byteswapio -rm  #-Mframe
#   FFLAGS = -I. -g -N132 -h byteswapio -rm -Rbc  #-Mframe
#-Mipa=fast #-Munroll
  ifneq (,$(findstring CAF,$(CPPFLAGS)))
    FFLAGS += -hcaf
  endif
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -Oomp,thread2
  else
    FFLAGS += -hnoomp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -sreal64
  endif
  CPP = /lib/cpp -traditional
  FFT_LIB = -L$(FFTW_DIR) -lfftw3 -lfftw3f
  FFT_INC = -I$(FFTW_INC)
  HDF_LIB =
  HDF_INC =
#  OTHERLIBS = -L ~/lib -lprofiler
endif
#######################################################################
ifeq ($(MACHINE),BW) ############################################
# compile with
# module load fftw craype-hugepages2M cray-hdf5-parallel ; module unload darshan
# to use Cray Fortran module,
# module swap PrgEnv-pgi PrgEnv-cray
  # run module load fftw/ before make
  FF = ftn #-C
  FFLAGS = -I. -Ofp3,scalar3,thread0 -N132 #-h byteswapio -rm  #-Mframe
  FFLAGS90 = -I. -Ofp3,scalar3,thread0 -N132 #-h byteswapio -rm  #-Mframe
#    FFLAGS = -I. -Ofp3,scalar3,thread0 -N132 -h byteswapio -rm  #-Mframe
#   FFLAGS = -I. -Ofp3,scalar3,thread0 -N132 #-Rbc #-h byteswapio -rm  #-Mframe
#-Mipa=fast #-Munroll
  ifneq (,$(findstring CAF,$(CPPFLAGS)))
    FFLAGS += -hcaf
    FFLAGS90 += -hcaf
  endif
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -Oomp,thread2
    FFLAGS90 += -Oomp,thread2
  else
    FFLAGS += -hnoomp
    FFLAGS90 += -hnoomp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -sreal64
    FFLAGS90 += -sreal64
  endif
  CPP = /lib/cpp -traditional
  FFT_LIB = -L$(FFTW_DIR) -lfftw3 -lfftw3f
  FFT_INC = -I$(FFTW_INC)
  MKL_LIB =
  MKL_INC =
  HDF_LIB =
  HDF_INC =
#  OTHERLIBS = -L ~/lib -lprofiler
endif
#######################################################################
ifeq ($(MACHINE),mira) ############################################
CPPFLAGS += -DIBM
  FF = mpixlf95_r
  CC = mpixlc_r
  CFLAGS += -DNUS_XCOMP
  FFLAGS = -O3 -qstrict -qfixed=132 -qarch=qp -qtune=qp -qfloat=hsflt:fltint
FFLAGS90 = -O3 -qstrict -qfree=F90 -qarch=qp -qtune=qp -qfloat=hsflt:fltint
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -qsmp=omp
  else
    FFLAGS += -hnoomp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -qrealsize=8
  endif
  CPP = /lib/cpp -traditional
  FFT_LIB = -L /soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3 -lfftw3f
  FFT_INC = -I /soft/libraries/alcf/current/xl/FFTW3/include
  HDF_LIB = -L /soft/libraries/hdf5/current/cnk-xl/current/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L /soft/libraries/alcf/current/xl/ZLIB/lib -lz
  HDF_INC = -I /soft/libraries/hdf5/current/cnk-xl/current/include
LDFLAGS = -L /soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3 -lfftw3f
    endif
#  #######################################################################
ifeq ($(MACHINE),mira.essl) ############################################
CPPFLAGS += -DESSL
CPPFLAGS += -DIBM
  FF = mpixlf95_r
  CC = mpixlc_r
  CFLAGS += -DNUS_XCOMP
  FFLAGS = -O3 -qstrict -qfixed=132 -qarch=qp -qtune=qp -qfloat=hsflt:fltint #for FFT kernel
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -qsmp=omp
  else
    FFLAGS += -hnoomp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -qrealsize=8
  endif
  CPP = /lib/cpp -traditional
  #FFT_LIB = -L$(DEPSPATH)/fftw-3.2.2/lib -lfftw3
# FFT_LIB = -L /soft/libraries/essl/current/lib64 -lesslbg
  FFT_LIB = -L /soft/libraries/essl/current/lib64 -lesslsmpbg
  #FFT_INC = -I$(DEPSPATH)/fftw-3.2.2/include -I$(DRIVER)/comm/xl.legacy/include
  FFT_INC = -I /soft/libraries/essl/current/include
   endif
#######################################################################
ifeq ($(MACHINE),fermi) ############################################
CPPFLAGS += -DIBM
  FF = mpixlf95_r
  CC = mpixlc_r
  CFLAGS += -DNUS_XCOMP
  FFLAGS = -O3 -qstrict -qfixed=132 -qarch=qp -qtune=qp -qfloat=hsflt:fltint
FFLAGS90 = -O3 -qstrict -qfree=F90 -qarch=qp -qtune=qp -qfloat=hsflt:fltint
  ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
    FFLAGS += -qsmp=omp
  else
    FFLAGS += -hnoomp
  endif
  ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
    FFLAGS += -qrealsize=8
  endif
  CPP = /lib/cpp -traditional
  FFT_LIB = -L /soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3 -lfftw3f
  FFT_INC = -I /soft/libraries/alcf/current/xl/FFTW3/include
  HDF_LIB = -L /soft/libraries/hdf5/current/cnk-xl/current/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L /soft/libraries/alcf/current/xl/ZLIB/lib -lz
  HDF_INC = -I /soft/libraries/hdf5/current/cnk-xl/current/include
LDFLAGS = -L /soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3 -lfftw3f
    endif
########################################################################




SRCDIR = .srcmake

# choose initial condition file to compile against
# these files should all be located in hooks/initial_conditions
INITCOND=pfield.o
#INITCOND=taylorgreen.o

LDFLAGS= $(FFLAGS) $(LFLAGS) $(FLAGS_PROF_DEBUG)
ifeq ($(MACHINE),BW)
LDFLAGS += -dynamic
endif


# these are core routines associated with the timestepping
OBJ2d_a = rungekutta_module.o proc4a.o proc4b.o wavespace_module.o wavespace.o yitransform.o itransform.o timestep_module.o realspace.o transform.o advanc_1.o advanc_2.o advanc_3.o advanc_4.o phshift.o phase1.o  itransform_vel.o realspace_vel.o transform_vel.o wavespace_vel.o phase1_om.o itransform_om.o om_write.o helical.o

# initialization routines
OBJ2d_b = $(OBJ2d_a) param_set.o com_set.o comp_set.o comsp_set.o input.o mpisetup.o abrt.o 

# modules
#OBJ2d = module.o lagmodule.o stratification_module.o  $(OBJ2d_b)
OBJ2d = module.o lagmodule.o IO.o $(OBJ2d_b)

OBJEUL = main.o chkstf.o clreop.o covcof.o init.o masks.o pcheck.o field.o hist1.o inpen.o inpen_overlap.o inpen_read.o inpen_read_r.o inpvf.o incond.o phymom.o scdlpdf.o scgout.o scgpdf.o scgrad.o scjout.o shifts.o sptr.o step.o trfsp.o vgmomt.o vorsp1.o vorsp2.o vorout.o waveno.o escout.o eulout.o fopen1.o fclose1.o clapf.o blanks.o ranseq.o seqput.o ranu.o ranu2.o vrand.o raniuo.o eModel.o enfsym.o openf.o store.o rstow.o istow.o irest.o rrest.o force.o forcdt.o vfuo.o vfiuo.o rann2.o scnpdf.o prni.o time_stamp.o out1d.o corr1d.o out_timers.o ek_init.o ek_force.o ek_shell.o epfor_param.o outpen.o outpen_write.o sptvar.o set_ia_xz.o write_timings.o write_iotimings.o outpen_write_r.o press_stat.o press_write.o consym.o  inpen_read_bw.o outpen_write_bw.o incomm.o incomm_row.o outcomm.o outcomm_row.o dissenst.o phyave.o logspace.o phypdfkrsp.o hist1vip.o phymom4.o check_consym2.o check_consym.o OUTPEN_HDF5.o OUTROW_HDF5.o RESTART_HDF5_DRIVER.o RESTART_HDF5.o AXI_SPTR.o bijdecomp.o

# these are files that link either to a hook or a null hook
OBJHOOK= partic.o compart_set.o #bouss.o gauss3d.o hsptr.o noise3d.o #all_fields_wave.o all_fields_real.o stratification_module.o noise3d.o

# core FFT pieces
OBJFFT = compi_alltoall.o escfft.o esrfft.o escpft.o esrpft.o kxtran.o xktran.o xkcomm1.o xkcomm2.o kxcomm1.o kxcomm2.o kxcomm2t.o  fft_init.o divide_thr.o kxcomm1_clean.o xkcomm2_clean.o

# linking channel flow files
ifneq (,$(findstring CHANNEL,$(CPPFLAGS)))
  OBJCHANNEL= laps.o cfdiff.o
endif

# linking lagrangian particle tracking files
ifneq (,$(findstring LAG,$(CPPFLAGS)))
OBJLAG = bs_setup.o inputlag.o intbf_ars.o lagout1.o lagout2.o order.o partin.o part2sc.o partsp.o  partvg.o solsp1.o spcal_ars.o splx.o sply.o splz.o splz_m.o spxyz_m.o tetra.o udxphy.o wrtpop1.o wrtpos.o valid1.o pop_file4.o pop_write4.o  pop2_file5.o pop2_write5.o rdate.o rk_part.o rk_advpos.o part_interp.o check_pvrms.o lagpart.o partic_rk4.o partsp_1.o rk4_advpos.o partsp_2.o partsp_3.o partsp_4.o gltime2.o lag_iotimings.o lag_timings.o init_mol_random.o molin.o rk_advmpos.o mol_file3.o mol_write3.o ppminmax.o ppminmaxv.o part2in.o bsminmax.o bsminmaxv.o splx_minmax.o splz_minmax.o  spcal_caf.o setrandlist.o compi_a2ar.o spcal_tree.o part2in_ar.o molrann.o boxmul.o spcal_pipe.o wrtpos_pipe.o spcal_local.o  pptoid.o update_part.o pop1_write_local.o pop2_write_local.o tetra3.o partin_ar.o tetra3_local.o wrtpos_h5.o readpos_h5.o pop1_file_local.o pop1_write_local2.o wrtpos_comm_h5.o readpos_comm_h5.o lagcomm_setup.o lagout_comm_h5.o lagout_local.o
endif


# include stratification files
ifneq (,$(findstring BOUSS,$(CPPFLAGS)))
  OBJSTRAT = # force_lindborg.o force_eirate.o #buoy.o bouss.o #noise.o
endif

ifneq (,$(findstring BGW,$(CPPFLAGS)))
  OBJ2d += pset_func.o
endif

# add objects for FFTWs
ifneq (,$(findstring FFTW,$(CPPFLAGS)))
  OBJFFT += epfftw.o 
ifeq (,$(findstring NOTWISE,$(CPPFLAGS)))
  OBJFFT += wisdom.o
endif
endif

ifneq (,$(findstring MHD,$(CPPFLAGS)))
  OBJEUL += sptr_mhd.o
endif

#----------------------------------------------------
# ----- set name of executable based on choice of major directives:
# ----- single vs double prec, Pure MPI vs OpenMP, CAF vs non-CAF

ifneq (,$(findstring OUT1D,$(CPPFLAGS)))
  OBJFFT += rfft99.o fourt.o
endif

ifneq (,$(findstring OPENMP,$(CPPFLAGS)))
NAME1 = DNS2d
else
NAME1 = DNS2d_mpi
endif

NAME2 = p4
ifneq (,$(findstring DOUBLE_PREC,$(CPPFLAGS)))
 NAME2 = p8
endif

#NAME3 = _
ifneq (,$(findstring CAF,$(CPPFLAGS)))
 NAME3 = _CAF
endif

ifneq (,$(findstring LOCAL,$(VPATH)))
NAME4 = _LOCAL
else
NAME4 = 
endif

EXEC = ${NAME1}_${NAME2}${NAME3}${NAME4}.x

#----------------------------------------------------


OBJ =  $(OBJ2d) $(OBJHOOK) $(INITCOND) $(OBJEUL) $(OBJSTRAT) $(OBJFFT) $(OBJLAG) $(OBJCHANNEL) check_velmax.o

$(EXEC): $(OBJ)
	$(FF) $(LDFLAGS) -o $(EXEC) $(MKL_INC) $(FFT_INC) $(HDF_INC) $(OBJ) $(MKL_LIB) $(FFT_LIB) $(HDF_LIB) $(OTHERLIBS)
#	$(FF) $(LDFLAGS) -dynamic -o $(EXEC) $(MKL_INC) $(FFT_INC) $(OBJ) $(MKL_LIB) $(FFT_LIB) $(OTHERLIBS)

ifdef FLAGS_PROF_DEBUG 
	@echo  "-------------------------------------------------------"
	@echo  "NOTE : code compiled with profiling/debugging options."
	@echo  "       It may slow down execution."
	@echo  "FLAGS_PROF_DEBUG = " $(FLAGS_PROF_DEBUG)
endif

$(OBJ2d) $(OBJEUL) $(OBJFFT) : module.f IO.F90
lagmodule.o: lagmodule.f module.f IO.F90
$(OBJLAG) : lagmodule.o module.f IO.F90
$(INITCOND): module.f IO.F90

ifneq (,$(findstring BGW,$(CPPFLAGS)))
pset_func.o: pset_func.C
	$(CC) -c pset_func.C
endif	

%.o: %.c
		$(CC) -c $(CFLAGS) $(MKL_INC) $(FFT_INC) $(HDF_INC) *.c
%.o: %.f
ifndef MACHINE
		@echo "Error: MACHINE is not defined"; exit 1;
endif

ifneq (exists, $(shell if [ -d $(SRCDIR) ]; then echo exists; fi))
		@echo "Error: $(SRCDIR) folders have not been created"
		@echo "Use the command 'make srcmake'"; exit 1
endif
		$(CPP) $(CPPFLAGS) $< $(SRCDIR)/$<
		$(FF) -c $(FFLAGS) $(FLAGS_PROF_DEBUG) $(MKL_INC) $(FFT_INC) $(HDF_INC) $(SRCDIR)/$<

%.o: %.f90
		$(CPP) $(CPPFLAGS) $< $(SRCDIR)/$<
		$(FF) -c $(FFLAGS) $(FLAGS_PROF_DEBUG) $(MKL_INC) $(FFT_INC) $(HDF_INC) $(SRCDIR)/$<

%.o: %.F90
		$(CPP) $(CPPFLAGS) $< $(SRCDIR)/$<
		$(FF) -c $(FFLAGS90) $(FLAGS_PROF_DEBUG) $(MKL_INC) $(FFT_INC) $(HDF_INC) $(SRCDIR)/$<
#this rule constructs the folder structure needed for the source directories
srcmake:
	@echo deleting $(SRCDIR)
	rm -rf $(SRCDIR)
	@echo constructing $(SRCDIR)

        # building null hook structure
	mkdir -v -p $(SRCDIR)/$(NHOOK)/diagnostics
	mkdir -v -p $(SRCDIR)/$(NHOOK)/initial_conditions
	mkdir -v -p $(SRCDIR)/$(NHOOK)/LAG
	mkdir -v -p $(SRCDIR)/$(NHOOK)/stratification

        # building hook directory
	mkdir -v -p $(SRCDIR)/$(HOOK)/diagnostics
	mkdir -v -p $(SRCDIR)/$(HOOK)/initial_conditions
	mkdir -v -p $(SRCDIR)/$(HOOK)/LAG
	mkdir -v -p $(SRCDIR)/$(HOOK)/LAG_PIPE
	mkdir -v -p $(SRCDIR)/$(HOOK)/LAG_LOCAL
	mkdir -v -p $(SRCDIR)/$(HOOK)/stratification
	mkdir -v -p $(SRCDIR)/$(HOOK)/stratification/BOUSS

        # building homogeneous / inhomogeneous directory
	mkdir -v -p $(SRCDIR)/homogeneous
	mkdir -v -p $(SRCDIR)/inhomogeneous/channel
	mkdir -v -p $(SRCDIR)/inhomogeneous/bl
	@echo "srcmake folders generated: code can now be compiled"

#remove object files
clean:
	- /bin/rm *.o *.mod *~ *.lst
	@echo files cleaned

# this removes the srcmake directory in addition to all object files
cleanest:
	- /bin/rm *.o *.mod *.lst
	@echo files cleaned
	rm -rf $(SRCDIR)
	@echo removed $(SRCDIR)

# this removes object files for routines involving Co-Array Fortran
cleancaf:
	rm compi_alltoall.o kxcomm1_clean.o kxcomm1.o kxcomm2.o kxcomm2t.o transform_vel.o xkcomm1.o xkcomm2_clean.o xkcomm2.o
