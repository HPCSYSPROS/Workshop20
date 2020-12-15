
#TCLDIR=$(HOME)/tcl8.5.9-linux-arm64
TCLDIR=$(HOME)/tcl8.5.9-linux-arm64-threaded
TCLINCL=-I$(TCLDIR)/include
#TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl
TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl -lpthread
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

