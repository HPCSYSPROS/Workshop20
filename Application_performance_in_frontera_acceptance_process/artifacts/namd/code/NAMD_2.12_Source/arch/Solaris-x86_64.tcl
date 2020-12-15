
TCLDIR=/Projects/namd2/tcl/tcl8.5.9-solaris-x86_64-threaded
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -lnsl -lsocket
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

