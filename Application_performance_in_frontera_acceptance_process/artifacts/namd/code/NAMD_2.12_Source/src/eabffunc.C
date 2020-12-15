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

#ifndef FUNC_CPP
#define FUNC_CPP

#include "eabffunc.h"


// the trim method of string
// accept a string, remove the space in the left and right of the string
// return the reference of the same string
std::string& eabffunc::trim(std::string &s)
{
	if (s.empty())
	{
		return s;
	}

	s.erase(0, s.find_first_not_of(" "));
	s.erase(s.find_last_not_of(" ") + 1);
	return s;
}


// the split of string
// accept a string, return a vector<string>
void eabffunc::split(const std::string &s, std::vector<std::string>& ret)
{
	std::stringstream temp(s);
	std::string token;
	while (temp >> token)
	{
		ret.push_back(token);
	}
}

// convert string to int
int eabffunc::chartoint(const std::string& c)
{
	std::stringstream temp(c);
	int token;
	temp >> token;
	return token;
}

// convert string to double
double eabffunc::chartodouble(const std::string& c)
{
	std::stringstream temp(c);
	double token;
	temp >> token;
	return token;
}

// convert double to int
int eabffunc::doubletoint(const double a)
{
	if (a > 0)
		return int(a + 0.000001);
	else
		return int(a - 0.000001);
}

#endif
