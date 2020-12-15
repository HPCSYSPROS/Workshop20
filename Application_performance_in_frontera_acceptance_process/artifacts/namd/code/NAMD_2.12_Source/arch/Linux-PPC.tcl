
TCLDIR=/Projects/namd2/tcl/linux-ppc
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -lnamdtcl8.3
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

