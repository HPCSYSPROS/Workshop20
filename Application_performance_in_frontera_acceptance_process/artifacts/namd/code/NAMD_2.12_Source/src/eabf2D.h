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

#ifndef EABF2D_H
#define EABF2D_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "eabf.h"

class eABF2D : public eABF
{
public:
    // called when (re)start a namd run
	eABF2D(const double _lowerboundary,
		const double _upperboundary,
		const double _width,
		const double _krestr,
		const double _lowerboundary2,
		const double _upperboundary2,
		const double _width2,
		const double _krestr2,
		const std::string& _outputfile,
		const int _outputfreq,
		const bool _restart,
		const std::string& _inputfile,
		const bool _outputgrad,
		const int _gradfreq,
		const double _temperature);

    virtual ~eABF2D()
    {
		for (int i = 0; i < bins[0] * bins[0]; i++)
			delete[] countall[i];
		delete[] countall;
    }

    // called by colvar each step, update the sampling point
    bool update(const std::string&);


private:
    // for 2D eABF, much memories are need, I have to use C-style pointer and array
    // BE VERY CAREFUL ABOUT IT!
    int** countall; // the distribution of samples, needed for eABF
    std::vector<std::vector<double> > sumx1; // the sum of x in each bin, needed for eABF
    std::vector<std::vector<double> > sumx21; // the sum of x^2, needed for eABF
    std::vector<std::vector<double> > sumx2;
	std::vector<std::vector<double> > sumx22;
    std::vector<std::vector <int> > county; // the distribution in each bin, needed for eABF
    //////////////////////////////////////////////////////////

    // initialize the variables
    bool initialize();

    // read variables from a restart file
    bool readfile();

    // write restart file
    bool writefile() const;

    // write the head of grad file
    bool writehead(ofstream_namd&) const;

    // calculate grad and pmf
	bool calpmf() const;

	// convert the scale in file to the scale in this eabf
	virtual int convertscale(double lowerboundary, int window, int dimension) const
	{
		if (dimension == 1)
			return int((lowerboundary - this->lowerboundary[0]) / width[0] + window);
		if (dimension == 2)
			return int((lowerboundary - this->lowerboundary[1]) / width[1] + window);
		return 0;
	}
};
#endif // EABF2D_H
