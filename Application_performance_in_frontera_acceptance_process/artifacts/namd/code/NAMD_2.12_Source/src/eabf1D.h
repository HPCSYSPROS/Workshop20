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

#ifndef EABF1D_H
#define EABF1D_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "eabf.h"

class eABF1D: public eABF
{
public:
	eABF1D() {}
	// called when (re)start a namd run
	eABF1D(const double _lowerboundary,
		const double _upperboundary,
		const double _width,
		const double _krestr,
		const std::string& _outputfile,
		const int _outputfreq,
		const bool _restart,
		const std::string& _inputfile,
		const bool _outputgrad,
		const int _gradfreq,
		const double _temperature);
	virtual ~eABF1D() {}

	bool update(const std::string&);

protected:

	std::vector<std::vector<int> > countall; // the distribution of samples, needed for eABF
	std::vector<double> sumx; // the sum of x in each bin, needed for eABF
	std::vector<double> sumx2; // the sum of x^2, needed for eABF
	std::vector<int> county; // the distribution in each bin, needed for eABF

	bool initialize();

	bool readfile();

	bool writefile() const;

	bool writehead(ofstream_namd&) const;

	bool calpmf() const;

	int convertscale(double lowerboundary, int window) const
	{
		return int((lowerboundary - this->lowerboundary[0]) / width[0] + window);
	}

};

#endif // EABF1D_H
