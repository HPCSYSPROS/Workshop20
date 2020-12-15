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

#ifndef EABF1D_CPP
#define EABF1D_CPP

#include <cmath>
#include <fstream>
#include <string>
#include <memory>
#include <iomanip>
#include "fstream_namd.h"

#include "eabf1D.h"
#include "eabffunc.h"


// constructor
eABF1D::eABF1D(const double _lowerboundary, const double _upperboundary, const double _width,
	const double _krestr,
	const std::string& _outputfile,
	const int _outputfreq,
	const bool _restart, const std::string& _inputfile,
	const bool _outputgrad, const int _gradfreq, const double _temperature) :
	eABF(_outputfile, _outputfreq, _restart, _inputfile, _outputgrad, _gradfreq, _temperature)
{
	lowerboundary.push_back(_lowerboundary);
	upperboundary.push_back(_upperboundary);
	width.push_back(_width);
	krestr.push_back(_krestr);
	this->initialize();
	if (restart == true)
		this->readfile();
}

bool eABF1D::initialize()
{
	// calculate the number of bins in the calculation
	// the bins is equal to "real bins + 20" to prevent boundary effect
    min.push_back(eabffunc::doubletoint(lowerboundary[0] / width[0]) - 10);
    max.push_back(eabffunc::doubletoint(upperboundary[0] / width[0]) + 10);
	bins.push_back(max[0] - min[0]);

    // initialize the distribution of "countall"
	for (int i = 0; i < bins[0]; i++)
    {
		countall.push_back(std::vector<int>(bins[0], 0));
    }

    // initialize other variables
	sumx = std::vector<double>(bins[0], 0);
	sumx2 = std::vector<double>(bins[0], 0);
	county = std::vector<int>(bins[0], 0);
    line = 0;
    return true;
}

bool eABF1D::update(const std::string& colvarstring)
{
    // this is the string obtained from colvar each step
    // the format of this string is the same as those in colvar.traj file

	std::vector<std::string> data;
	eabffunc::split(colvarstring, data);  // split the string to a vector<string>

    int step;

    double abf_force;       // the value of the two columns characterizing the value of CV and extended CV
    double eabf_force;

	step = eabffunc::chartoint(data[0]);    // read the variables fron the "colvarstring"
	abf_force = eabffunc::chartodouble(data[col[0]]);
	eabf_force = eabffunc::chartodouble(data[col[1]]);


	// for dihedral RC, it is possible that abf_force = 179 and eabf_force = -179, should correct it
    if(abf_force > 150 && eabf_force < -150)
        eabf_force += 360;
    else if(abf_force < -150 && eabf_force > 150)
        eabf_force -= 360;

	// normalize the index of the value of CV and extended CV
	int binx = int(floor(abf_force / width[0])) - min[0];
	int biny = int(floor(eabf_force / width[0])) - min[0];

    // update the 
	if (binx < 0 || binx >= bins[0] || biny < 0 || biny >= bins[0])
        return false;
    else
    {
		countall[binx][biny]++;
        // if this is the first time add "countall", then the line in outputfile will increase
        if(countall[binx][biny]==1)
            line++;
        county[biny]++;
    }

    sumx[biny] += abf_force;
    sumx2[biny] += abf_force * abf_force;

	if (step%outputfreq == 0)
    {
        writefile();
    }

    if(outputgrad==1&&(step%gradfreq)==0)
    {
        calpmf();
    }

    return true;
}

bool eABF1D::readfile()
{
    std::ifstream file(inputfile.c_str(), std::ios::in);

    // output file format: (also see writefile function)
    // (number of lines) lowerboundary width (number of bins) krestr temperature
	// (normalized CV) (normalized extended CV) numbers
	// ......  "only non-zero value is recorded here"
	// (normalized CV) sumx sumx2 county
	// ......
    int temp_line, temp_bins;
    double temp_lowerboundary, temp_width;
	file >> temp_line >> temp_lowerboundary >> temp_width >> temp_bins >> krestr[0] >> temperature;

    temp_bins = temp_bins + 20;

	int x = 0, y = 0, temp_countall;
	for (int i = 0; i < temp_line; i++)
    {
		file >> x >> y;
		file >> temp_countall;
		//if (x > temp_bins - 10 || x < 10)
		//	continue;
        ;
		if (countall[convertscale(temp_lowerboundary, x)][convertscale(temp_lowerboundary, y)] == 0)
			line++;
        countall[convertscale(temp_lowerboundary,x)][convertscale(temp_lowerboundary,y)] += temp_countall;
    }

	double temp_sumx, temp_sumx2;
	int temp_county;
	for (int i = 0; i < temp_bins; i++)
    {
		file >> x >> temp_sumx >> temp_sumx2 >> temp_county;
		//if (x > temp_bins - 10 || x < 10)
		//	continue;
        ;
		sumx[convertscale(temp_lowerboundary, x)] += temp_sumx;
		sumx2[convertscale(temp_lowerboundary, x)] += temp_sumx2;
		county[convertscale(temp_lowerboundary, x)] += temp_county;
    }

    file.close();
    return true;
}

bool eABF1D::writefile() const
{
    ofstream_namd file(outputfile.c_str(),".old");
    file<<std::setprecision(15);

	file << line << " " << lowerboundary[0] << " " << width[0] << " " << bins[0] - 20 << " " << krestr[0] << " " << temperature << std::endl;

	for (int i = 0; i < bins[0]; i++)
		for (int j = 0; j < bins[0]; j++)
			if (countall[i][j] != 0)
				// the number of line should be "this->line", I guess :)
				file << i << " " << j << " " << countall[i][j] << std::endl;

	for (int i = 0; i < bins[0]; i++)
		file << i << " " << sumx[i] << " " << sumx2[i] << " " << county[i] << " " << std::endl;

    file.close();
    return true;
}

bool eABF1D::writehead(ofstream_namd& os) const
{
	os << "# 1" << std::endl;
	os << "# " << lowerboundary[0] << " " << width[0] << " " << bins[0] - 20 << " " << 0 << std::endl;
	os << std::endl;
    return true;
}

bool eABF1D::calpmf() const
{
	// implementation of YangWei's estimator
	// mimic Jerome's script
	int norm;
	std::vector<double> x_av(bins[0]);
	std::vector<double> sigma2(bins[0]);
	for (int biny = 0; biny < bins[0]; biny++)
	{
		norm = county[biny] > 0 ? county[biny] : 1;
		x_av[biny] = sumx[biny] / norm;
		sigma2[biny] = sumx2[biny] / norm - x_av[biny] * x_av[biny];
	}

	static std::string gridfilename = outputfile + ".grad";
	static std::string histfilename = outputfile + ".hist.grad";
	static std::string pmffilename = outputfile + ".pmf";
	static std::string countfilename = outputfile + ".count";

	ofstream_namd ofile(gridfilename.c_str(),".old");
	ofstream_namd ofile_hist(histfilename.c_str(), std::ios::app);
	ofstream_namd ofile_pmf(pmffilename.c_str(),".old");
	ofstream_namd ofile_count(countfilename.c_str(),".old");

	writehead(ofile);
	writehead(ofile_hist);
	writehead(ofile_count);

	ofile_pmf << lowerboundary[0] << " 0.0" << std::endl;

	double pmf = 0, pmf_x = 0;
	double x, y, av, diff_av, grad;
	for (int binx = 0; binx < bins[0]; binx++)
	{
		x = (binx + min[0] + 0.5) * width[0];
		norm = 0;
		av = 0;
		diff_av = 0;

		for (int biny = 0; biny < bins[0]; biny++)
		{
			y = (biny + min[0] + 0.5) * width[0];
			norm += countall[binx][biny];

			if (sigma2[biny] > 0.00001 || sigma2[biny] < -0.00001)
			{
				av += countall[binx][biny] * (x - x_av[biny]) / sigma2[biny];
			}
			else if (countall[binx][biny] > 1)
			{
				// TODO, print error
			}
			diff_av += countall[binx][biny] * (x - y);
		}

		diff_av /= (norm > 0 ? norm : 1);
		av = eabffunc::BOLTZMANN * temperature * av / (norm > 0 ? norm : 1);
		grad = av - krestr[0] * diff_av;


		if (x >= lowerboundary[0] && x < upperboundary[0])
		{
			ofile << x << " " << grad << std::endl;
			ofile_hist << x << " " << grad << std::endl;
			pmf += grad * width[0];
			pmf_x = x + 0.5 * width[0];
			ofile_pmf << pmf_x << " " << pmf << std::endl;
			ofile_count << x << " " << norm << std::endl;
			// TODO, rescale the pmf file to prevent negetive free energy
		}
	}
	ofile.close();
	ofile_pmf.close();
	ofile_hist.close();
	ofile_count.close();
	return true;
}

#endif // EABF
