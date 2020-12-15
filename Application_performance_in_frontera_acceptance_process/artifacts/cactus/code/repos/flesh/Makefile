# /*@@
#   @file      Makefile
#   @date      Sun Jan 17 22:26:05 1999
#   @author    Tom Goodale
#   @desc
#   gnu Makefile for the CCTK.
#
# WARNING: This makefile may not function with "make".  Errors like
# make: file `Makefile' line 36: Must be a separator (: or ::) for rules (bu39)
# mean you should have used gmake.  gmake is available free from
#
#   http://www.gnu.org/software/make/
#
# and should be installed on all production systems.
#
# For information on how to use this makefile, type
#
#   gmake help
#
# Use VERBOSE=yes on the make command line to get verbose output of
# the compiler commands that make is running
#
#   @enddesc
# @@*/

################################################################################
# Version number
################################################################################
CCTK_VERSION_MAJOR = 4
CCTK_VERSION_MINOR = 2
CCTK_VERSION_OTHER = 3
################################################################################
CCTK_VERSION=$(CCTK_VERSION_MAJOR).$(CCTK_VERSION_MINOR).$(CCTK_VERSION_OTHER)
################################################################################
export CCTK_VERSION_MAJOR CCTK_VERSION_MINOR CCTK_VERSION_OTHER CCTK_VERSION
################################################################################

# Stop with prompts unless told not to
ifeq ($(strip $(PROMPT)), )
PROMPT = "yes"
endif

export PROMPT

# backwards compatibility for VERBOSE=yes (SILENT=no)
ifeq ($(shell echo $(strip $(SILENT)) | tr '[:upper:]' '[:lower:]'),no)
VERBOSE = yes
endif
# Make quietly unless told not to
ifneq ($(shell echo $(strip $(VERBOSE)) | tr '[:upper:]' '[:lower:]'),yes)
.SILENT:
endif

# Stuff for parallel makes
# TJOBS is the number of thorns to compile in parallel
ifeq ($(strip $(TJOBS)), )
TPARFLAGS =
else
TPARFLAGS = -j $(TJOBS)
endif

# FJOBS is the number of files within a thorn to compile in parallel
ifeq ($(strip $(FJOBS)), )
FPARFLAGS =
else
FPARFLAGS = -j $(FJOBS)
endif

export TPARFLAGS FPARFLAGS

# Directory for configuration options
# Set THORNLIST_DIR to "." if it is not set already, and if THORNLIST
# does not contain an absolute patch (starting with a slash)
ifeq ($(strip $(THORNLIST_DIR)), )
ifeq ($(shell echo $(THORNLIST) | grep '^/'), )
THORNLIST_DIR = "."
endif
endif

# End of parallel make stuff


# Set the options to pass to the setup script
ifneq ($(strip $(options)),)
SETUP_OPTIONS = -config_file=$(options)
else
SETUP_OPTIONS =
endif


# Allow various options to be passed to the configure script

SETUP_ENV =

ifdef CC
ifneq ($(strip $(origin CC)), default)
SETUP_ENV += CC="$(CC)" ; export CC ;
endif
endif

ifdef F90
ifneq ($(strip $(origin F90)), default)
SETUP_ENV += F90="$(F90)" ; export F90 ;
endif
endif

ifdef F77
ifneq ($(strip $(origin F77)), default)
SETUP_ENV += F77="$(F77)" ; export F77 ;
endif
endif

ifdef LD
ifneq ($(strip $(origin LD)), default)
SETUP_ENV += LD="$(LD)" ; export LD ;
endif
endif

ifdef CFLAGS
ifneq ($(strip $(origin CFLAGS)), default)
SETUP_ENV += CFLAGS="$(CFLAGS)" ; export CFLAGS;
endif
endif

ifdef F90FLAGS
ifneq ($(strip $(origin F90FLAGS)), default)
SETUP_ENV += F90FLAGS="$(F90FLAGS)" ; export F90FLAGS ;
endif
endif

ifdef F77FLAGS
ifneq ($(strip $(origin F77FLAGS)), default)
SETUP_ENV += F77FLAGS="$(F77FLAGS)" ; export F77FLAGS ;
endif
endif

ifdef LDFLAGS
ifneq ($(strip $(origin LDFLAGS)), default)
SETUP_ENV += LDFLAGS="$(LDFLAGS)" ; export LDFLAGS ;
endif
endif

ifdef REAL_PRECISION
ifneq ($(strip $(origin REAL_PRECISION)), default)
SETUP_ENV += REAL_PRECISION=$(REAL_PRECISION) ; export REAL_PRECISION ;
endif
endif

ifdef INTEGER_PRECISION
ifneq ($(strip $(origin INTEGER_PRECISION)), default)
SETUP_ENV += INTEGER_PRECISION=$(INTEGER_PRECISION) ; export INTEGER_PRECISION ;
endif
endif

# Arrangement options

ifdef MPI
ifneq ($(strip $(origin MPI)), default)
SETUP_ENV += MPI="$(MPI)" ; export MPI ;
endif
endif

# Debug options
ifdef DEBUG
ifneq ($(strip $(origin DEBUG)), default)
SETUP_ENV += DEBUG="$(DEBUG)" ; export DEBUG ;
endif
endif

# Optimisation options
ifdef OPTIMISE
ifneq ($(strip $(origin OPTIMISE)), default)
SETUP_ENV += OPTIMISE="$(OPTIMISE)" ; export OPTIMISE ;
endif
endif

# Compile-time options

# Warning options
ifdef WARN
ifneq (, $(filter-out yes no, $(WARN)))
$(error Didn't recognize setting of WARN="$(WARN)" (should be either "yes" or "no"))
endif
CCTK_WARN_MODE=$(WARN)
else
CCTK_WARN_MODE=no
endif

export CCTK_WARN_MODE

# Various auxilary programs
PERL = perl
SETUP    = lib/make/setup_configuration.pl
NEWTHORN = lib/make/new_thorn.pl
BUILD_ACTIVETHORNS = lib/sbin/BuildActiveThorns.pl

# Dividers to make the screen output slightly nicer
DIVEL   =  __________________
DIVIDER =  $(DIVEL)$(DIVEL)$(DIVEL)$(DIVEL)
ifeq ($(shell echo $(VERBOSE) | tr '[:upper:]' '[:lower:]'),yes)
define NOTIFY_DIVIDER
	echo $(DIVIDER)
endef
else
define NOTIFY_DIVIDER
	:
endef
endif

# Work out where we are
export CCTK_HOME := $(shell pwd)


# Work out where the configuration directory is
ifdef CACTUS_CONFIGS_DIR
CONFIGS_DIR = $(CACTUS_CONFIGS_DIR)
else
CONFIGS_DIR = $(CCTK_HOME)/configs
endif

export CONFIGS_DIR

# Work out which configurations are available
CONFIGURATIONS = $(patsubst $(CONFIGS_DIR)/%,%,$(wildcard $(CONFIGS_DIR)/*))
CONFIGINFOS = $(wildcard $(CONFIGS_DIR)/*/config-info)

# Default target does nothing.
# Used to set up a default based upon uname or something.
.PHONY:default-target

default-target:
ifeq ($(strip $(CONFIGURATIONS)),)
	$(NOTIFY_DIVIDER)
	@echo No configurations defined.
	@echo Please use \'$(MAKE) \<name\>\' to setup a configuration called \<name\>.
	$(NOTIFY_DIVIDER)
	@echo \'$(MAKE) help\' lists all $(MAKE) options.
else
ifeq ($(words $(CONFIGURATIONS)), 1)
	@echo Using the only available configuration: $(CONFIGURATIONS).
	$(NOTIFY_DIVIDER)
	$(MAKE) $(CONFIGURATIONS)
else
	@echo Known configurations are: $(CONFIGURATIONS)
	@echo Please use $(MAKE) \<configuration\>
	$(NOTIFY_DIVIDER)
	@echo \'$(MAKE) help\' lists all $(MAKE) options.
endif
endif
	$(NOTIFY_DIVIDER)

# Target to build a configuration
.PHONY: $(CONFIGURATIONS)

$(CONFIGURATIONS):
	if test ! -f "$(CONFIGS_DIR)/$@/config-data/cctk_Config.h" ; then \
	  $(NOTIFY_DIVIDER);\
	  echo "Cactus - version: $(CCTK_VERSION)";\
          if test "x$(PROMPT)" = 'xno'; then\
            if (! $(SETUP_ENV) $(PERL) -s $(SETUP) $(SETUP_OPTIONS) $@) ; then \
              echo "" ;                                                        \
              echo "Error reconfiguring configuration $@" ;                    \
              rm -f "$(CONFIGS_DIR)/$@/config-data/cctk_Config.h" ;            \
              exit 2 ;                                                         \
            fi                                                                 \
          else                                                                 \
	    echo "Error: Configuration $@ is incomplete.";\
	    echo "Please check the files in $(CONFIGS_DIR)/$@/config-data for error messages.";\
	    echo "You can try again to configure using $(MAKE) $@-config";\
	    echo "or delete this configuration with $(MAKE) $@-delete.";\
	    $(NOTIFY_DIVIDER);\
	    exit 1; \
          fi                                                                 \
	fi
	if ($(PERL) -e 'exit ((stat shift)[9] > (stat shift)[9])' $(CONFIGS_DIR)/$@/config-info $(CCTK_HOME)/lib/make/force-reconfigure); then    \
	  $(NOTIFY_DIVIDER);\
	  echo "Cactus - version: $(CCTK_VERSION)";\
	  echo "Error: Configuration $@ is out of date.";\
	  echo "  Please reconfigure your configuration by running the command"; \
          echo ;\
          echo "    $(MAKE) $@-reconfig"; \
          echo ;\
	  echo "  (It is likely that recent changes to the flesh require this.)";\
	  $(NOTIFY_DIVIDER);\
	  exit 1;\
	fi
	if test "x${MAKELEVEL}" = "x0" ; then \
	  $(NOTIFY_DIVIDER);\
	  echo "Cactus - version: $(CCTK_VERSION)"; \
	  echo "Building configuration $@"; \
	  $(NOTIFY_DIVIDER);\
	fi
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$@ CCTK_HOME=$(CCTK_HOME) $(TPARFLAGS) rebuild
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$@ CCTK_HOME=$(CCTK_HOME) $(TPARFLAGS)

# Clean target
.PHONY: distclean

distclean:
	$(NOTIFY_DIVIDER)
	@echo Deleting all your configurations !
	rm -rf $(CONFIGS_DIR)
	$(NOTIFY_DIVIDER)

# Targets to make tags files

# Mark these targets phony to force an update when gmake TAGS is done.
.PHONY: TAGS tags

TAGS:
	$(NOTIFY_DIVIDER)
	@echo Updating the Emacs TAGS file
	: > TAGS
	find src arrangements -follow \( \
		-name '*.[chCfF]' -o -name '*.[fF]77' -o -name '*.[fF]90' -o \
		-name '*.cc' -o -name '*.cxx' -o -name '*.hh' -o \
		-name '*.[ch]pp' -o -name '*.inc' \) -print | \
	xargs etags -a
	$(PERL) -pi.bak -e 's/(subroutine\s*)([a-zA-Z0-9_]+)/\1\L\2/g;' TAGS
	rm -f TAGS.bak
	$(NOTIFY_DIVIDER)

tags:
	$(NOTIFY_DIVIDER)
	@echo Updating the vi tags file
	: > tags
	find src arrangements -follow \( \
		-name '*.[chCfF]' -o -name '*.[fF]77' -o -name '*.[fF]90' -o \
		-name '*.cc' -o -name '*.cxx' -o -name '*.hh' -o \
		-name '*.[ch]pp' -o -name '*.inc' \) -print | \
	xargs ctags -a
	$(PERL) -pi.bak -e 's/(subroutine\s*)([a-zA-Z0-9_]+)/\1\L\2/g;' tags
	rm -f tags.bak
	sort tags > sortedtags
	mv sortedtags tags
	$(NOTIFY_DIVIDER)

# Make a new configuration with a default name
.PHONY: default

default:
	$(NOTIFY_DIVIDER)
	@echo Running the configuration program
	$(SETUP_ENV) $(PERL) -s $(SETUP) $(SETUP_OPTIONS)
	$(NOTIFY_DIVIDER)
	@echo You are now ready to build the CCTK.
	@echo This is done by $(MAKE) \<configuration\>
	$(NOTIFY_DIVIDER)

# The help system.
.PHONY: help

help:
	$(NOTIFY_DIVIDER)
	@echo "****************************** "
	@echo "* Welcome to the Cactus Code *"
	@echo "******************************"
ifeq ($(strip $(CONFIGURATIONS)),)
	@echo There are no configurations currently specified.
	@echo \'$(MAKE) \<name\>\' will run a setup script to set up a configuration called \'\<name\>\'.
else
	@echo The following configurations are currently specified
	@echo
	@echo "  $(CONFIGURATIONS)"
	$(NOTIFY_DIVIDER)
	@echo "To build a configuration: "
	@echo "  run $(MAKE) followed by the name of a configuration."
	$(NOTIFY_DIVIDER)
	@echo There is a range of options available to act on a configuration.
	@echo These are activated by '$(MAKE) \<conf-name\>-\<option\>'
	@echo Valid options are
	@echo "  -build         : build individual thorns of a configuration."
	@echo "  -clean         : clean a configuration"
	@echo "                  (deletes all object and dependency files"
	@echo "                   in the configuration)."
	@echo "  -cleandeps     : clean a configuration's dependency files."
	@echo "  -cleanobjs     : clean a configuration's object files."
	@echo "  -config        : create a new configuration, or reconfigure an existing one"
	@echo "                  (overwrites previous configuration options)."
	@echo "  -configinfo    : display the configuration options for a configuration."
	@echo "  -delete        : delete a configuration."
	@echo "  -editthorns    : edit the ThornList file."
	@echo "  -realclean     : restore a configuration to almost a new state."
	@echo "                  (deletes all but the config-data directory"
	@echo "                   and the ThornList file)."
	@echo "  -rebuild       : rebuild a configuration (forces the CST to be rerun)."
	@echo "  -reconfig      : reconfigure an existing configuration"
	@echo "                   using its previous configuration options."
	@echo "  -utils         : build a configuration's utility programs."
	@echo "  -testsuite     : run the test suites."
	@echo "  -thornlist     : regenerate the ThornList file."
	@echo "  -ThornGuide    : create the thorn manual for a specific configuration."
	@echo "  -update        : update the files for a specific configuration from CVS and/or SVN."
	@echo "  -examples      : copy thorn parameter files to examples directory."
endif
	$(NOTIFY_DIVIDER)
	@echo There are options available to act on thorns or arrangements.
	@echo These are activated by \'$(MAKE) \<thorn-name\>-\<option\>\' or
	@echo   \'$(MAKE) \<arrangement-name\>-\<option\>\' respectively.
	@echo Valid options are
	@echo "  -ThornDoc           - produce the documentation for the thorn"
	@echo "                        in doc/ThornDoc/<arrangement>/<thorn-name>/."
	@echo "  -ThornDocHTML       - produce the documentation for the thorn in HTML format"
	@echo "                        in doc/ThornDoc/<arrangement>/<thorn-name>/."
	@echo "  -ArrangementDoc     - produce documentation for the arrangement"
	@echo "                        in doc/ArrangementDoc/<arrangement-name>/."
	@echo "  -ArrangementDocHTML - produce documentation for the arrangement in HTML format"
	@echo "                        in doc/ArrangementDoc/<arrangement-name>/."
	$(NOTIFY_DIVIDER)
	@echo $(MAKE) also knows the following targets
	@echo
	@echo "  checkout            - checkout public arrangements/thorns."
	@echo "  cvsdiff             - show differences between installed Cactus and"
	@echo "                        version in CVS repository."
	@echo "  cvsstatus           - report on status of Cactus (when installed from CVS)."
	@echo "  update              - update flesh and arrangements from CVS and/or SVN."
	@echo "  default             - create a new configuration with a default name."
	@echo "  distclean           - delete all existing configurations."
	@echo "  downsize            - remove non-essential files."
	@echo "  newthorn            - create a new thorn."
	@echo "  TAGS                - create an emacs TAGS file."
	@echo "  tags                - create a vi TAGS file."
	@echo "  thorninfo           - give information about all available thorns."
	@echo "  MaintGuide          - create maintainers manual."
	@echo "  MaintGuideHTML      - create maintainers manual in HTML format in"
	@echo "                        doc/HTML/MaintGuide."
	@echo "  UsersGuide          - create users manual doc/UsersGuide.pdf."
	@echo "  UsersGuideHTML      - create users manual in HTML format in"
	@echo "                        doc/HTML/UsersGuide/."
	@echo "  ReferenceManual     - create reference manual doc/ReferenceManual.pdf."
	@echo "  ReferenceManualHTML - create reference manual in HTML format in"
	@echo "                        doc/HTML/ReferenceManual/."
	@echo "                        doc/HTML/ThornGuide/."
	@echo "  ThornDoc            - create documentation for all thorns in doc/ThornDoc."
	@echo "  ThornDocHTML        - create documentation for all thorns in HTML format"
	@echo "                        in doc/ThornHTML."
	@echo "  ArrangementDoc      - create documentation for all arrangements"
	@echo "                        in doc/ArrangementDoc."
	@echo "  ArrangementDocHTML  - create documentation for all arrangements in HTML format"
	@echo "                        format in doc/ArrangementHTML."
	@echo "  AllDoc              - build all documentation."
	@echo "  AllDocHTML          - build all documentation in HTML format."
	@echo "  <anything else>     - prompt to create a configuration with that name."
	$(NOTIFY_DIVIDER)


# Version information
.PHONY: version

version:
	@echo "Cactus - version: $(CCTK_VERSION)"


# Version information for internal use
.PHONY: int_version

int_version:
	$(NOTIFY_DIVIDER)
	@echo "Cactus - version: $(CCTK_VERSION)"


############################################
# Build individual thorns of a configuration
############################################
.PHONY: build

build: int_version
	@echo Please specify a configuration to build.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -build, $(CONFIGURATIONS))

$(addsuffix -build, $(CONFIGURATIONS)): int_version
	if test "x$(BUILDLIST)" = "x"; then \
	  $(NOTIFY_DIVIDER); \
	  echo "Please specify the thorns to build with \"BUILDLIST=<list of thorns>\""; \
	  $(NOTIFY_DIVIDER); \
	else \
	  echo Building thorns \'$(BUILDLIST)\' of configuration $(@:%-build=%); \
	  cd $(CONFIGS_DIR)/$(@:%-build=%); \
	  $(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-build=%) CCTK_HOME=$(CCTK_HOME) build; \
	fi;
endif

%-build:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-build=%) does not exist.
	@echo Build aborted.


#######################
# Clean a configuration
#######################
.PHONY: clean

clean:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to clean.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -clean,$(CONFIGURATIONS))

$(addsuffix -clean,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Cleaning configuration $(@:%-clean=%)
	cd $(CONFIGS_DIR)/$(@:%-clean=%)
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-clean=%) CCTK_HOME=$(CCTK_HOME) clean
	$(NOTIFY_DIVIDER)

endif

%-clean:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-clean=%) does not exist.
	@echo Cleaning aborted.


#############################
# Clean just dependency files
#############################
.PHONY: cleandeps

cleandeps:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to clean the dependencies of.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -cleandeps,$(CONFIGURATIONS))

$(addsuffix -cleandeps,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Cleaning dependencies in configuration $(@:%-cleandeps=%)
	cd $(CONFIGS_DIR)/$(@:%-cleandeps=%)
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-cleandeps=%) CCTK_HOME=$(CCTK_HOME) cleandeps
	$(NOTIFY_DIVIDER)

endif

%-cleandeps:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-cleandeps=%) does not exist.
	@echo Cleaning dependencies aborted.


#########################
# Clean just object files
#########################
.PHONY: cleanobjs

cleanobjs:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to clean the object files of.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -cleanobjs,$(CONFIGURATIONS))

$(addsuffix -cleanobjs,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Cleaning configuration $(@:%-cleanobjs=%)
	cd $(CONFIGS_DIR)/$(@:%-cleanobjs=%)
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-cleanobjs=%) CCTK_HOME=$(CCTK_HOME) cleanobjs
	$(NOTIFY_DIVIDER)

endif

%-cleanobjs:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-cleanobjs=%) does not exist.
	@echo Cleaning object files aborted.


##########################################################
# Clean away all produced files (doesn't delete ThornList)
##########################################################
.PHONY: realclean

realclean:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to really clean.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -realclean,$(CONFIGURATIONS))

$(addsuffix -realclean,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Cleaning configuration $(@:%-realclean=%)
	cd $(CONFIGS_DIR)/$(@:%-realclean=%)
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-realclean=%) CCTK_HOME=$(CCTK_HOME) realclean
	$(NOTIFY_DIVIDER)

endif

%-realclean:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-realclean=%) does not exist.
	@echo Cleaning aborted.


########################
# Delete a configuration
########################
.PHONY: delete

delete:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to delete.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -delete,$(CONFIGURATIONS))

$(addsuffix -delete,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	if test "x$(DELETE_CONFIRMATION)" = "xyes" ; then          \
	  echo "Really delete configuration $(@:%-delete=%) (no)?";\
	  read confirm rest;                                       \
	  if test $$? -ne 0 ; then                                 \
	    confirm='no';                                          \
	  fi;                                                      \
	else                                                       \
	  confirm=yes ;                                            \
	fi ;                                                       \
	if test "x$$confirm" = "xyes" ; then                       \
	  echo Deleting configuration $(@:%-delete=%);             \
	  cd $(CONFIGS_DIR) ; rm -rf $(@:%-delete=%) ;             \
	fi
	$(NOTIFY_DIVIDER)
endif

%-delete:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-delete=%) does not exist.
	@echo Deletion aborted.


#########################
# Rebuild a configuration
#########################
.PHONY: rebuild

rebuild: int_version
	@echo Please specify a configuration to rebuild.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -rebuild,$(CONFIGURATIONS))

$(addsuffix -rebuild,$(CONFIGURATIONS)): int_version
	@echo Rebuilding configuration $(@:%-rebuild=%)
	if [ -r $(CONFIGS_DIR)/$(@:%-rebuild=%)/config-data/make.thornlist ] ; then rm  $(CONFIGS_DIR)/$(@:%-rebuild=%)/config-data/make.thornlist ; fi
	$(MAKE) $(@:%-rebuild=%)
endif

%-rebuild:
	@echo Configuration $(@:%-rebuild=%) does not exist.
	@echo Rebuild aborted.


#####################################
# Regenerate the compiled thorns list
#####################################
.PHONY: thornlist

thornlist:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to regenerate the thornlist of.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -thornlist,$(CONFIGURATIONS))

$(addsuffix -thornlist,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Regenerating compiled ThornList $(@:%-thornlist=%)
	if [ -r $(CONFIGS_DIR)/$(@:%-thornlist=%)/ThornList ] ; then rm $(CONFIGS_DIR)/$(@:%-thornlist=%)/ThornList ; fi
	$(MAKE) $(@:%-thornlist=%)
endif

%-thornlist:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-thornlist=%) does not exist.
	@echo Regeneration of compiled ThornList aborted.


####################
# Edit the thornlist
####################
.PHONY: editthorn

editthorns:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to edit the thornlist of.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -editthorns,$(CONFIGURATIONS))

$(addsuffix -editthorns,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Editing compiled ThornList $(@:%-editthorn=%)
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-editthorns=%) CCTK_HOME=$(CCTK_HOME) editthorns
endif

%-editthorns:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-editthorns=%) does not exist.
	@echo Editing of compiled ThornList aborted.


################################
# Rerun the configuration script
################################
.PHONY: config

config: int_version
	@echo Please specify a configuration to configure.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -config,$(CONFIGURATIONS))

$(addsuffix -config,$(CONFIGURATIONS)): int_version
	if test -z "$(THORNLIST)" || (test -n "$(THORNLIST)" && test -r "$(THORNLIST_DIR)/$(THORNLIST)") ; \
	then \
	  if ($(SETUP_ENV) $(PERL) -s $(SETUP) $(SETUP_OPTIONS) $(@:%-config=%)) ; then : ; else \
            echo "" ;                                                      \
            echo "Error reconfiguring $@" ;                                \
            rm -f "$(CONFIGS_DIR)/$(@:%-config=%)/config-data/cctk_Config.h";\
            exit 2                                 ;                       \
          fi ;                                                             \
	  if test -n "$(THORNLIST)" ; \
	  then \
	    cp $(THORNLIST_DIR)/$(THORNLIST) $(CONFIGS_DIR)/$(@:%-config=%)/ThornList;\
	  fi ; \
	  if test -n "$(THORNS)" ; \
          then \
            echo $(THORNS) >> $(CONFIGS_DIR)/$(@:%-config=%)/ThornList ; \
          fi ; \
	  $(NOTIFY_DIVIDER) ; \
	  if test "x$(PROMPT)" = "xno" ; then \
	    $(MAKE) $(@:%-config=%) WARN=$(WARN); \
	  else \
	    echo Use $(MAKE) $(@:%-config=%) to build the configuration. ; \
	  fi; \
	else \
	  echo "ThornList $(THORNLIST_DIR)/$(THORNLIST) does not exist" ; \
	  exit 2; \
	fi
	$(NOTIFY_DIVIDER)
endif

%-config:
	@echo Configuration $(@:%-config=%) does not exist.;
	if test "x$(PROMPT)" = "xyes" ; then \
	  echo Setup configuration $(@:%-config=%) \(yes\)?; \
          yesno='yes';\
	  if test $$? -ne 0 ; then \
	    yesno='no'; \
	  fi; \
	fi; \
	if [ "x$$yesno" = "xno" -o "x$$yesno" = "xn" -o "x$$yesno" = "xNO" -o "x$$yesno" = "xN" ] ;\
	then \
	  echo Setup of configuration $(@:%-config=%) cancelled ;     \
	else \
	  echo Setting up new configuration $(@:%-config=%); \
	  if test -z "$(THORNLIST)" || (test -n "$(THORNLIST)" && test -r "$(THORNLIST_DIR)/$(THORNLIST)") ; \
	  then \
	    if ($(SETUP_ENV) $(PERL) -s $(SETUP) $(SETUP_OPTIONS) $(@:%-config=%)) ; then : ; else \
              echo "" ; \
              echo "Error creating configuration $(@:%-config=%)" ; \
              exit 2; \
            fi ; \
	    if test -n "$(THORNLIST)"; \
	    then \
	      cp $(THORNLIST_DIR)/$(THORNLIST) $(CONFIGS_DIR)/$(@:%-config=%)/ThornList;\
	    fi ;\
	    if test -n "$(THORNS)" ; then \
	      echo $(THORNS) >> $(CONFIGS_DIR)/$(@:%-config=%)/ThornList ; \
	    fi ; \
	    $(NOTIFY_DIVIDER)   ;  \
	    if test "x$(PROMPT)" = "xno" ; then \
	      $(MAKE) $(@:%-config=%) WARN=$(WARN); \
	    else \
	      echo Use $(MAKE) $(@:%-config=%) to build the configuration.; \
	    fi; \
	  else \
	    echo "ThornList $(THORNLIST_DIR)/$(THORNLIST) does not exist" ; \
	    exit 2; \
	  fi ; \
	fi
	$(NOTIFY_DIVIDER)


################################################################
# Reconfigure an existing configuration (using the same options)
################################################################
.PHONY: reconfig

reconfig: int_version
	@echo Please specify a configuration to reconfigure.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -reconfig,$(CONFIGURATIONS))

$(addsuffix -reconfig,$(CONFIGURATIONS)): int_version
	if test ! -r "$(CONFIGS_DIR)/$(@:%-reconfig=%)/config-info"; then \
	  echo ""; \
	  echo "Error reconfiguring '$(@:%-reconfig=%)': configuration is incomplete."; \
	  echo "Use '$(MAKE) $(@:%-config=%)' to configure the configuration."; \
	  exit 2; \
	elif ! head -n 1 $(CONFIGS_DIR)/$(@:%-reconfig=%)/config-info | grep -q '# CONFIGURATION'; then \
	  echo "Error reconfiguring '$(@:%-reconfig=%)': unrecognized config-info file format" ; \
	  echo "This probably means that the Cactus config-info file format has changed"; \
	  echo "since this configuration was last configured.  You will have to freshly"; \
	  echo "configure the configuration with"; \
	  echo "  $(MAKE) $(@:%-reconfig=%-config) OPTION1=VALUE1 OPTION2=VALUE2 ..."; \
	  echo "Note that this will overwrite all of this configuration's current"; \
	  echo "configuration options.  See the config-info file"; \
	  echo "  $(CONFIGS_DIR)/$(@:%-reconfig=%)/config-info"; \
	  echo "to see what these are."; \
	  exit 2; \
	fi; \
	if ($(SETUP_ENV) $(PERL) -s $(SETUP) -config_file=$(CONFIGS_DIR)/$(@:%-reconfig=%)/config-info $(@:%-reconfig=%)) ; then \
	  $(MAKE) $(@:%-reconfig=%)-cleandeps ; \
        else \
          echo "" ;                                                      \
          echo "Error reconfiguring $@" ;                                \
          rm -f "$(CONFIGS_DIR)/$(@:%-reconfig=%)/config-data/cctk_Config.h";\
          exit 2                                 ;                       \
        fi ;                                                             \
	$(NOTIFY_DIVIDER) ; \
	if test "x$(PROMPT)" = "xno" ; then \
	  $(MAKE) $(@:%-reconfig=%) WARN=$(WARN); \
	else \
	  echo Use '$(MAKE) $(@:%-reconfig=%)' to build the configuration. ; \
	fi
	$(NOTIFY_DIVIDER)
endif

%-reconfig:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-reconfig=%) does not exist.;
	@echo Reconfiguration aborted.


#####################
# Build the utilities
#####################
.PHONY: utils

utils:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to build the utilities of.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -utils,$(CONFIGURATIONS))

$(addsuffix -utils,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Building utilities for $(@:%-utils=%)
	cd $(CONFIGS_DIR)/$(@:%-utils=%)
	$(MAKE) -f $(CCTK_HOME)/lib/make/make.configuration TOP=$(CONFIGS_DIR)/$(@:%-utils=%) CCTK_HOME=$(CCTK_HOME) utils UTILS=$(UTILS) CONFIG_NAME=$(@:%-utils=%)
	$(NOTIFY_DIVIDER)

endif

%-utils:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-utils=%) does not exist.
	@echo Building of utilities aborted.


##################
# Make a new thorn
##################
.PHONY: newthorn
newthorn:
	$(NOTIFY_DIVIDER)
	@echo Creating a new thorn
	$(PERL) -s $(NEWTHORN);
	$(NOTIFY_DIVIDER)


###################
# Run the testsuite
###################
.PHONY: testsuite

testsuite: int_version
	@echo Please specify a configuration to test.
	$(NOTIFY_DIVIDER)


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -testsuite,$(CONFIGURATIONS))

$(addsuffix -testsuite,$(CONFIGURATIONS)):
	@echo Running test suite $(@:%-thornlist=%)
	if [ -r $(CONFIGS_DIR)/$(@:%-testsuite=%)/ThornList ] ; then $(PERL) -s lib/sbin/RunTest.pl $(PROMPT) $(CCTK_HOME) $(@:%-testsuite=%); fi
endif

%-testsuite:
	@echo Configuration $(@:%-testsuite=%) does not exist.
	@echo Test suite aborted.



##########################################

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -examples,$(CONFIGURATIONS))

# Copy thorn parameter files
.PHONY: examples

examples:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration.
	$(NOTIFY_DIVIDER)


$(addsuffix -examples,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Copying parameter files $(@:%-examples=%)
	if [ -r $(CONFIGS_DIR)/$(@:%-examples=%)/ThornList ] ; then $(PERL) lib/sbin/CopyParFiles.pl $(@:%-examples=%) ; fi
endif

%-examples:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-examples=%) does not exist.
	@echo Parameter file copying aborted.


# Show configuration information

.PHONY: configinfo

configinfo:
ifeq ($(strip $(CONFIGURATIONS)),)
	$(NOTIFY_DIVIDER)
	@echo No configurations defined.
	$(NOTIFY_DIVIDER)
else
	cat $(CONFIGINFOS)
endif
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -configinfo,$(CONFIGURATIONS))

$(addsuffix -configinfo,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Displaying configuration information
	cat configs/$(@:%-configinfo=%)/config-info
endif

%-configinfo:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-configinfo=%) does not exist.
	@echo Displaying configuration information aborted.


# Create sysinfo file

.PHONY: sysinfo

sysinfo:
	$(NOTIFY_DIVIDER)
	@echo Please specify a configuration to run sysinfo with.
	$(NOTIFY_DIVIDER)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -sysinfo,$(CONFIGURATIONS))

$(addsuffix -sysinfo,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Running SystemInfo
	$(PERL) ./lib/sbin/SystemInfo.pl $(@:%-sysinfo=%)
endif

%-sysinfo:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-sysinfo=%) does not exist.
	@echo Getting system info aborted.


###############################################################################
#                      Documentation targets
###############################################################################

# Make all documentation
.PHONY: AllDoc
AllDoc: UsersGuide ReferenceManual MaintGuide ThornDoc ArrangementDoc

# Make all HTML documentation
.PHONY: AllDocHTML
AllDocHTML: UsersGuideHTML ReferenceManualHTML MaintGuideHTML ThornDocHTML ArrangementDocHTML

# Make the Users Guide

.PHONY: UsersGuide.pdf
UsersGuide.pdf: UsersGuide

.PHONY: UsersGuide
UsersGuide:
	$(NOTIFY_DIVIDER)
	@echo Creating user documentation UsersGuide.pdf
	cd doc/UsersGuide;                          \
	echo "  Running pdflatex....";              \
	pdflatex  -interaction=nonstopmode UsersGuide.tex > LATEX_MESSAGES 2>&1; \
	pdflatex  -interaction=nonstopmode UsersGuide.tex > LATEX_MESSAGES 2>&1; \
	pdflatex  -interaction=nonstopmode UsersGuide.tex > LATEX_MESSAGES 2>&1; \
	if grep "^\! " "LATEX_MESSAGES"; then                            \
	  echo "  Problem in $<.  See doc/UsersGuide/LATEX_MESSAGES.";      \
	  exit 1;                                                           \
	elif grep "^LaTeX Warning:" "LATEX_MESSAGES"; then                  \
	  echo "  For more information see doc/UsersGuide/LATEX_MESSAGES."; \
	fi;                                                                 \
	mv UsersGuide.pdf $(CCTK_HOME)/doc/UsersGuide.pdf 
	@echo "  UsersGuide.pdf created in doc directory."
	@echo "  Done."
	$(NOTIFY_DIVIDER)

.PHONY: UsersGuideHTML

doc/UsersGuide/bincactus2.ps: doc/UsersGuide/bincactus2.pdf
	pdf2ps doc/UsersGuide/bincactus2.pdf doc/UsersGuide/bincactus2.ps

UsersGuideHTML: doc/UsersGuide/bincactus2.ps
	$(NOTIFY_DIVIDER)
	@echo "Creating user documentation (HTML)"
	cd doc/UsersGuide;                         \
	$(CCTK_HOME)/lib/sbin/ConvertFigures;      \
	echo "  Running htlatex....";              \
	htlatex UsersGuide.tex "html,2,fn-in" "" "" -interaction=nonstopmode > LATEX_MESSAGES 2>&1; \
	if [ $$? -ne 0 ]; then                                                \
	  echo "  Problem executing htlatex. See doc/UsersGuide/LATEX_MESSAGES."; \
		exit 1;                                                           \
	fi;                                                                 \
	if grep "^\! " "LATEX_MESSAGES"; then                               \
	  echo "  Problem in $<.  See doc/UsersGuide/LATEX_MESSAGES.";      \
	  exit 1;                                                           \
	elif grep "^LaTeX Warning:" "LATEX_MESSAGES"; then                  \
	  echo "  For more information see doc/UsersGuide/LATEX_MESSAGES."; \
	fi;                                                                 \
  mkdir -p $(CCTK_HOME)/doc/HTML/UsersGuide;            \
	cp UsersGuide*.png $(CCTK_HOME)/doc/HTML/UsersGuide/; \
	cp UsersGuide*.html $(CCTK_HOME)/doc/HTML/UsersGuide/; \
	cp UsersGuide.css $(CCTK_HOME)/doc/HTML/UsersGuide/
	@echo "  Users Guide (HTML) created in doc/HTML/UsersGuide directory."
	@echo "  Done."
	$(NOTIFY_DIVIDER)

# Make the Reference Manual

.PHONY: ReferenceManual.pdf
ReferenceManual.pdf: ReferenceManual

.PHONY: ReferenceManual
ReferenceManual:
	$(NOTIFY_DIVIDER)
	@echo Creating user reference manual ReferenceManual.pdf
	cd doc/ReferenceManual;                     \
	echo "  Running pdflatex....";              \
	pdflatex  -interaction=nonstopmode ReferenceManual.tex > LATEX_MESSAGES 2>&1; \
	pdflatex  -interaction=nonstopmode ReferenceManual.tex > LATEX_MESSAGES 2>&1; \
	pdflatex  -interaction=nonstopmode ReferenceManual.tex > LATEX_MESSAGES 2>&1; \
	if grep "^\! " "LATEX_MESSAGES"; then                            \
	  echo "  Problem in $<.  See doc/ReferenceManual/LATEX_MESSAGES."; \
	  exit 1;                                                           \
	elif grep "^LaTeX Warning:" "LATEX_MESSAGES"; then    \
	  echo "  For more information see doc/ReferenceManual/LATEX_MESSAGES."; \
	fi;                                         \
	mv ReferenceManual.pdf $(CCTK_HOME)/doc/ReferenceManual.pdf
	@echo "  ReferenceManual.pdf created in doc directory."
	@echo "  Done."
	$(NOTIFY_DIVIDER)

.PHONY: ReferenceManualHTML
ReferenceManualHTML:
	$(NOTIFY_DIVIDER)
	@echo "Creating user reference manual (HTML)"
	cd doc/ReferenceManual;                    \
	$(CCTK_HOME)/lib/sbin/ConvertFigures;      \
	echo "  Running htlatex....";              \
	htlatex ReferenceManual.tex "html,2,fn-in" "" "" -interaction=nonstopmode  > LATEX_MESSAGES 2>&1; \
	if [ $$? -ne 0 ]; then                                                \
	  echo "  Problem executing htlatex. See doc/UsersGuide/LATEX_MESSAGES."; \
		exit 1;                                                           \
	fi;                                                                 \
	if grep "^\! " "LATEX_MESSAGES"; then                            \
	  echo "  Problem in $<.  See doc/ReferenceManual/LATEX_MESSAGES."; \
	  exit 1;                                                           \
	elif grep "^LaTeX Warning:" "LATEX_MESSAGES"; then    \
	  echo "  For more information see doc/ReferenceManual/LATEX_MESSAGES."; \
	fi;                                         \
	mkdir -p $(CCTK_HOME)/doc/HTML/ReferenceManual;                \
	cp ReferenceManual*.png $(CCTK_HOME)/doc/HTML/ReferenceManual; \
	cp ReferenceManual*.html $(CCTK_HOME)/doc/HTML/ReferenceManual; \
	cp ReferenceManual.css  $(CCTK_HOME)/doc/HTML/ReferenceManual
	@echo "  HTML ReferenceManual created in doc/HTML/ReferenceManual directory."
	@echo "  Done."
	$(NOTIFY_DIVIDER)

# Make the Maintainers' Guide

.PHONY: MaintGuide.pdf
MaintGuide.pdf: MaintGuide

.PHONY: MaintGuide
MaintGuide:
	$(NOTIFY_DIVIDER)
	@echo Creating maintainers documentation MaintGuide.pdf
	cd doc/MaintGuide;                          \
	echo "  Running pdflatex....";              \
	pdflatex  -interaction=nonstopmode MaintGuide.tex > LATEX_MESSAGES 2>&1; \
	pdflatex  -interaction=nonstopmode MaintGuide.tex > LATEX_MESSAGES 2>&1; \
	pdflatex  -interaction=nonstopmode MaintGuide.tex > LATEX_MESSAGES 2>&1; \
	if grep "^\! " "LATEX_MESSAGES"; then                               \
	  echo "  Problem in $<.  See doc/MaintGuide/LATEX_MESSAGES.";      \
	  exit 1;                                                           \
	elif grep "^LaTeX Warning:" "LATEX_MESSAGES"; then                  \
	  echo "  For more information see doc/MaintGuide/LATEX_MESSAGES."; \
	fi;                                                                 \
	mv MaintGuide.pdf $(CCTK_HOME)/doc/MaintGuide.pdf
	@echo "  MaintGuide.pdf created in doc directory."
	@echo "  Done."
	$(NOTIFY_DIVIDER)

.PHONY: MaintGuideHTML
MaintGuideHTML:
	$(NOTIFY_DIVIDER)
	@echo "Creating maintainers documentation (HTML)"
	cd doc/MaintGuide;                          \
	$(CCTK_HOME)/lib/sbin/ConvertFigures;             \
	echo "  Running htlatex....";              \
	htlatex MaintGuide.tex "html,2,fn-in" "" "" -interaction=nonstopmode  > LATEX_MESSAGES 2>&1; \
	if [ $$? -ne 0 ]; then                                                \
	  echo "  Problem executing htlatex. See doc/UsersGuide/LATEX_MESSAGES."; \
		exit 1;                                                           \
	fi;                                                                 \
	if grep "^\! " "LATEX_MESSAGES"; then                               \
	  echo "  Problem in $<.  See doc/MaintGuide/LATEX_MESSAGES.";      \
	  exit 1;                                                           \
	elif grep "^LaTeX Warning:" "LATEX_MESSAGES"; then                  \
	  echo "  For more information see doc/MaintGuide/LATEX_MESSAGES."; \
	fi;                                                                 \
	mkdir -p $(CCTK_HOME)/doc/HTML/MaintGuide;                          \
	cp MaintGuide*.png $(CCTK_HOME)/doc/HTML/MaintGuide;                \
	cp MaintGuide*.html $(CCTK_HOME)/doc/HTML/MaintGuide;                \
	cp MaintGuide.css  $(CCTK_HOME)/doc/HTML/MaintGuide
	@echo "  Maintainers Guide (HTML) created in doc/HTML/MaintGuide directory."
	@echo "  Done."
	$(NOTIFY_DIVIDER)

# Run ThornGuide on a configuration

DOCDIR		= $(CCTK_HOME)/doc
CONFIGNAME	= $(@:%-ThornGuide=%)
CONFIGDIR	= $(CONFIGS_DIR)/$(CONFIGNAME)
CONFIGOCDIR	= $(CONFIGDIR)/doc
CONFIGBUILDDIR	= $(CONFIGDIR)/doc/build
GUIDENAME	= ThornGuide-$(CONFIGNAME)

ifneq ($strip($(CONFIGURATIONS)),)
.PHONY: $(addsuffix -ThornGuide,$(CONFIGURATIONS))

$(addsuffix -ThornGuide,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Creating ThornGuide for configuration $(CONFIGNAME)
	cd $(CONFIGDIR); \
	mkdir -p doc
	rm -rf $(CONFIGBUILDDIR)
	mkdir $(CONFIGBUILDDIR)
	if test -r $(CONFIGDIR)/ThornList ; then \
	  cd $(CONFIGBUILDDIR); \
	  $(MAKE) -f $(DOCDIR)/ThornGuide/Makefile THORNLIST=$(CONFIGDIR)/ThornList MASTER_FILE=$(GUIDENAME) DOCBUILDDIR=$(CONFIGBUILDDIR); \
	  if test -e "$(CONFIGBUILDDIR)/$(GUIDENAME).pdf"; then \
	    mv "$(CONFIGBUILDDIR)/$(GUIDENAME).pdf" $(DOCDIR)/$(GUIDENAME).pdf; \
	    echo "  $(GUIDENAME).pdf created in doc directory."; \
	    echo "  Done."; \
	  fi \
        else \
          echo "  Error: $(CONFIGDIR)/ThornList not found."; \
	fi
endif

%-ThornGuide:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-ThornGuide=%) does not exist.
	@echo Thorn Guide creation aborted.

# Make the ThornGuide

.PHONY: ThornDoc
%-ThornDoc:
	$(NOTIFY_DEVIDER)
	@lib/sbin/ThornDoc $(@:%-ThornDoc=%)
	$(NOTIFY_DIVIDER)
ThornDoc:
	$(NOTIFY_DEVIDER)
	@lib/sbin/ThornDoc
	$(NOTIFY_DIVIDER)
%-ArrangementDoc:
	$(NOTIFY_DEVIDER)
	@lib/sbin/ArrangementDoc $(@:%-ArrangementDoc=%)
	$(NOTIFY_DIVIDER)
ArrangementDoc:
	$(NOTIFY_DEVIDER)
	@lib/sbin/ArrangementDoc
	$(NOTIFY_DIVIDER)

.PHONY: ThornDocHTML
%-ThornDocHTML: doc/UsersGuide/bincactus2.ps
	$(NOTIFY_DEVIDER)
	@lib/sbin/ThornDocHTML $(@:%-ThornDocHTML=%)
	$(NOTIFY_DIVIDER)
ThornDocHTML: doc/UsersGuide/bincactus2.ps
	$(NOTIFY_DEVIDER)
	@lib/sbin/ThornDocHTML
	$(NOTIFY_DIVIDER)
%-ArrangementDocHTML: doc/UsersGuide/bincactus2.ps
	$(NOTIFY_DEVIDER)
	@lib/sbin/ArrangementDocHTML $(@:%-ArrangementDocHTML=%)
	$(NOTIFY_DIVIDER)
ArrangementDocHTML: doc/UsersGuide/bincactus2.ps
	$(NOTIFY_DEVIDER)
	@lib/sbin/ArrangementDocHTML
	$(NOTIFY_DIVIDER)

###############################################################################
#                      End of documentation targets
###############################################################################

# Rule to show thorn information

.PHONY: thorninfo
thorninfo:
	$(NOTIFY_DIVIDER)
	@echo Displaying info for all thorns in the arrangements directory
	$(PERL) -s $(BUILD_ACTIVETHORNS) $(CCTK_HOME)/arrangements | cat;
	$(NOTIFY_DIVIDER)

# Processed CVS information

.PHONY: cvsstatus

cvsstatus:
	$(PERL) -s $(CCTK_HOME)/lib/sbin/CVSStatus.pl

# run cvsudpate on a configuration

.PHONY: update

cvsupdate:
	$(PERL) -s $(CCTK_HOME)/lib/sbin/CVSUpdate.pl arrangements
update:
	$(PERL) -s $(CCTK_HOME)/lib/sbin/CVSUpdate.pl arrangements


ifneq ($strip($(CONFIGURATIONS)),)
.PHONY $(addsuffix -cvsupdate,$(CONFIGURATIONS)):
.PHONY $(addsuffix -update,$(CONFIGURATIONS)):

$(addsuffix -cvsupdate,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Updating files for configuration $(@:%-cvsupdate=%)
	if test -r $(CONFIGS_DIR)/$(@:%-cvsupdate=%)/ThornList ; then \
          $(PERL) -s lib/sbin/CVSUpdate.pl arrangements $(CONFIGS_DIR)/$(@:%-cvsupdate=%)/ThornList; \
        fi
	@echo Done.

$(addsuffix -update,$(CONFIGURATIONS)):
	$(NOTIFY_DIVIDER)
	@echo Updating files for configuration $(@:%-update=%)
	if test -r $(CONFIGS_DIR)/$(@:%-update=%)/ThornList ; then \
          $(PERL) -s lib/sbin/CVSUpdate.pl arrangements $(CONFIGS_DIR)/$(@:%-update=%)/ThornList; \
        fi
	@echo Done.
endif

%-cvsupdate:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-cvsupdate=%) does not exist.
	@echo CVS Update aborted.

%-update:
	$(NOTIFY_DIVIDER)
	@echo Configuration $(@:%-update=%) does not exist.
	@echo Update aborted.

.PHONY: cvsdiff

cvsdiff:
	$(PERL) -s $(CCTK_HOME)/lib/sbin/CVSStatus.pl -case=diff

# Remove non-essential files

.PHONY: downsize
downsize:
	$(NOTIFY_DIVIDER)
	@echo Remove flesh and thorn documentation \(\no\)?
	read yesno rest; \
	if test $$? -ne 0 ; then \
	  yesno='no'; \
	fi; \
	if [ "x$$yesno" = "xyes" -o "x$$yesno" = "xy" -o "x$$yesno" = "xYES" -o "x$$yesno" = "xY" ] ;\
	then  \
	rm -rf doc; rm -rf arrangements/*/*/doc; \
	$(NOTIFY_DIVIDER)   ;  \
	fi
	@echo Remove thorn testsuites \(\no\)?
	read yesno rest; \
	if test $$? -ne 0 ; then \
	  yesno='no'; \
	fi; \
	if [ "x$$yesno" = "xyes" -o "x$$yesno" = "xy" -o "x$$yesno" = "xYES" -o "x$$yesno" = "xY" ] ;\
	then  \
	rm -rf arrangements/*/*/test; \
	$(NOTIFY_DIVIDER)   ;  \
	fi
	@echo Remove all configurations \(\no\)?
	read yesno rest; \
	if test $$? -ne 0 ; then \
	  yesno='no'; \
	fi; \
	if [ "x$$yesno" = "xyes" -o "x$$yesno" = "xy" -o "x$$yesno" = "xYES" -o "x$$yesno" = "xY" ] ;\
	then  \
	$(MAKE) distclean; \
	$(NOTIFY_DIVIDER)   ;  \
	fi

# Last resort rule.  Assume it is the name of a configuration

%::
	$(NOTIFY_DIVIDER)
	@echo "Cactus - version: $(CCTK_VERSION)"
	if test "x$(PROMPT)" = "xyes" ; then \
	  echo Setup configuration $@ \(yes\)?; \
	  read yesno rest; \
	  if test $$? -ne 0 ; then \
	    yesno='no'; \
	  fi; \
	fi; \
	if [ "x$$yesno" = "xno" -o "x$$yesno" = "xn" -o "x$$yesno" = "xNO" -o "x$$yesno" = "xN" ] ; \
	then  \
	  echo Setup of configuration $@ cancelled ; \
	else \
	  echo Setting up new configuration $@ ; \
	  if test -z "$(THORNLIST)" || (test -n "$(THORNLIST)" && test -r "$(THORNLIST_DIR)/$(THORNLIST)") ; \
	  then \
	  if ($(SETUP_ENV) $(PERL) -s $(SETUP) $(SETUP_OPTIONS) $@) ; then : ; else\
            echo "" ;                                                      \
            echo "Error creating configuration $@" ;                       \
            rm -f "$(CONFIGS_DIR)/$@/config-data/cctk_Config.h";           \
            exit 2                                 ;                       \
          fi ;                                                             \
	  if test -n "$(THORNLIST)" ; \
	  then \
	    echo Using ThornList $(THORNLIST_DIR)/$(THORNLIST) ; \
	    cp $(THORNLIST_DIR)/$(THORNLIST) $(CONFIGS_DIR)/$@/ThornList ; \
	  fi ; \
	  if test -n "$(THORNS)" ; \
          then \
            echo $(THORNS) >> $(CONFIGS_DIR)/$@/ThornList ; \
          fi ; \
	  $(NOTIFY_DIVIDER) ;  \
	  if test "x$(PROMPT)" = "xno" ; then \
	    $(MAKE) $(@:%-config=%) WARN=$(WARN); \
	  else \
	    echo Use $(MAKE) $@ to build the configuration. ; \
	  fi; \
	  else \
	    echo "ThornList $(THORNLIST_DIR)/$(THORNLIST) does not exist" ; \
	    exit 2; \
	  fi ; \
	fi
	$(NOTIFY_DIVIDER)
