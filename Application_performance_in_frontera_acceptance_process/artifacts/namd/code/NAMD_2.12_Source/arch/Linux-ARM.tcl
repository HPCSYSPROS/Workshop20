
#TCLDIR=/Projects/namd2/tcl/tcl8.5.9-linux-arm
TCLDIR=/Projects/namd2/tcl/tcl8.5.9-linux-arm-threaded
TCLINCL=-I$(TCLDIR)/include
#TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl
TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl -lpthread
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

