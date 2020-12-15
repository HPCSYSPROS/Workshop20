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

#ifndef EABFFUNC_H
#define EABFFUNC_H

#include <string>
#include <sstream>
#include <vector>

namespace eabffunc
{

    const double BOLTZMANN=0.00198721;

    // the trim method of string
    // accept a string, remove the space in the left and right of the string
    // return the reference of the same string
	std::string& trim(std::string &s);


    // the split of string
    // accept a string, return a vector<string>
	void split(const std::string &s, std::vector<std::string>& ret);

    // convert string to int
	int chartoint(const std::string& c);

    // convert string to double
	double chartodouble(const std::string& c);

	// convert double to int
	int doubletoint(const double);
}

#endif // EABFFUNC_H
