//
// The extended adaptive biasing force method has been contributed to NAMD by the following authors:
//
//    Haohao Fu and Christophe Chipot
//    Laboratoire International Associ\'e
//    Centre National de la Recherche Scientifique et University of Illinois at Urbana--Champaign
//    Unit\'e Mixte de Recherche No. 7565, Universit\'e de Lorraine
//    B.P. 70239, 54506 Vand\oe uvre-lès-Nancy cedex, France
//
// Copyright 2016, Centre National de la Recherche Scientifique
//

#ifndef LINK_CPP
#define LINK_CPP

#include "eabf.h"
#include "eabf1D.h"
#include "eabf2D.h"
#include "eabffunc.h"

#ifdef NAMD_TCL
#include <tcl.h>

/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
////////////// Here is NAMD interface ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// the point to the runtime parseABF object
static eABF* parse = 0;

// the interface to Tcl, startrun
static int Tcl_startrun(ClientData useless, Tcl_Interp* stupidTcl, int argc, Tcl_Obj *CONST objv[])
{
    std::string* argv = new std::string[argc];
    // the last objv is the dimension of free-energy surface
    for(int i=0; i<argc; i++)
    {
        argv[i] = Tcl_GetString(objv[i]);
    }

	// 1D eABF
	if (eabffunc::chartoint(argv[argc - 1]) == 1)
    {
		parse = new eABF1D(eabffunc::chartodouble(argv[1]),      // lowerboundary
			eabffunc::chartodouble(argv[2]),      // upperboundary
			eabffunc::chartodouble(argv[3]),      // width
			eabffunc::chartodouble(argv[4]),    //krestr
			std::string(argv[5]),             // outputfile
			eabffunc::chartoint(argv[6]),         // outputfreq
			bool(eabffunc::chartoint(argv[7])),   // restart
			std::string(argv[8]),             // inputfile
			bool(eabffunc::chartoint(argv[9])),  // outputgrad
			eabffunc::chartoint(argv[10]),        // gradfreq
			eabffunc::chartodouble(argv[11]));     // temperature
    }

	// 2D eABF
	if (eabffunc::chartoint(argv[argc - 1]) == 2)
    {
		parse = new eABF2D(eabffunc::chartodouble(argv[1]),      // lowerboundary
			eabffunc::chartodouble(argv[2]),      // upperboundary
			eabffunc::chartodouble(argv[3]),      // width
			eabffunc::chartodouble(argv[4]),      // krestr
			eabffunc::chartodouble(argv[5]),      // lowerboundary2
			eabffunc::chartodouble(argv[6]),      // upperboundary2
			eabffunc::chartodouble(argv[7]),      // width2
			eabffunc::chartodouble(argv[8]),       // krestr2
			std::string(argv[9]),             // outputfile
			eabffunc::chartoint(argv[10]),         // outputfreq
			bool(eabffunc::chartoint(argv[11])),   // restart
			std::string(argv[12]),             // inputfile
			bool(eabffunc::chartoint(argv[13])),  // outputgrad
			eabffunc::chartoint(argv[14]),        // gradfreq
			eabffunc::chartodouble(argv[15]));     // temperature
    }

	delete[] argv;
    return TCL_OK;
}

// the interface to tcl, updaterun
static int Tcl_updaterun(ClientData useless, Tcl_Interp* stupidTcl, int argc, Tcl_Obj *CONST argv[])
{
	std::string colvar = Tcl_GetString(argv[1]);
	parse->update(colvar);
    return TCL_OK;
}

// set col
static int Tcl_setcol(ClientData useless, Tcl_Interp* stupidTcl, int argc, Tcl_Obj *CONST argv[])
{
	std::vector<int> cols;
	for (int i = 1; i<argc; i++)
	{
		cols.push_back(eabffunc::chartoint(Tcl_GetString(argv[i])));
	}
	parse->setcolumn(cols);
	return TCL_OK;
}

// input: dimension lowerboundary upperboundary width outputfile file1 file2 ...
// or: dimension lowerboundary upperboundary width lowerboundary2 upperboundary2 width2 outputfile file1 file2...
static int Tcl_mergefile(ClientData useless, Tcl_Interp* stupidTcl, int argc, Tcl_Obj *CONST objv[])
{
	std::string* argv = new std::string[argc];
	// the last objv is the dimension of free-energy surface
	for (int i = 0; i < argc; i++)
	{
		argv[i] = Tcl_GetString(objv[i]);
	}

	// 1D eABF
	if (eabffunc::chartoint(argv[1]) == 1)
	{
		parse = new eABF1D(eabffunc::chartodouble(argv[2]),      // lowerboundary
			eabffunc::chartodouble(argv[3]),      // upperboundary
			eabffunc::chartodouble(argv[4]),      // width
			0,                                 //krestr
			std::string(argv[5]),             // outputfile
			0,         // outputfreq
			0,   // restart
			"",             // inputfile
			0,  // outputgrad
			0,        // gradfreq
			0);    // temperature

		// merge file
		for (int i = 6; i < argc; i++)
		{
			parse->mergefile(argv[i], std::string(argv[5]));
		}
	}

	if (eabffunc::chartoint(argv[1]) == 2)
	{
		parse = new eABF2D(eabffunc::chartodouble(argv[2]),      // lowerboundary
			eabffunc::chartodouble(argv[3]),      // upperboundary
			eabffunc::chartodouble(argv[4]),      // width
			0,                                  //krestr
			eabffunc::chartodouble(argv[5]),      // lowerboundary2
			eabffunc::chartodouble(argv[6]),      // upperboundary2
			eabffunc::chartodouble(argv[7]),      // width2
			0,                                 //krestr2
			std::string(argv[8]),             // outputfile
			0,         // outputfreq
			0,   // restart
			"",             // inputfile
			0,  // outputgrad
			0,   // gradfreq
			0);        //temperature

			// merge file 
		for (int i = 9; i < argc; i++)
		{
			parse->mergefile(argv[i], std::string(argv[8]));
		}
	}

	delete[] argv;
	return TCL_OK;
}

extern "C" int Eabf_Init(Tcl_Interp* stupidTcl)
{
	// create Tcl process "startrun" and "updaterun"
	Tcl_CreateObjCommand(stupidTcl, "startrun", Tcl_startrun, (ClientData)0, (Tcl_CmdDeleteProc*)0);
	Tcl_CreateObjCommand(stupidTcl, "updaterun", Tcl_updaterun, (ClientData)0, (Tcl_CmdDeleteProc*)0);
	Tcl_CreateObjCommand(stupidTcl, "mergefile", Tcl_mergefile, (ClientData)0, (Tcl_CmdDeleteProc*)0);
	Tcl_CreateObjCommand(stupidTcl, "setcol", Tcl_setcol, (ClientData)0, (Tcl_CmdDeleteProc*)0);
	Tcl_PkgProvide(stupidTcl, "Eabf", "1.0.0");

	Tcl_Eval(stupidTcl, "puts \"eABF [package require Eabf]\"");

	return TCL_OK;
}

int eabf_static_init(Tcl_Interp *interp) {
  Tcl_StaticPackage(0,"Eabf",Eabf_Init,0);
  return Tcl_Eval(interp,"package ifneeded Eabf 1.0.0 {load {} Eabf}");
}

#endif // NAMD_TCL
#endif // LINK_CPP
