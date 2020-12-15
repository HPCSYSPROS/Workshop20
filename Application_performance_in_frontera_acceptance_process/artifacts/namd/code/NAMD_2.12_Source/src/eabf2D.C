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

#ifndef EABF2D_CPP
#define EABF2D_CPP

#include <cmath>
#include <fstream>
#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include "fstream_namd.h"

#include "eabf2D.h"
#include "eabffunc.h"

// constructor
eABF2D::eABF2D(const double _lowerboundary,
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
	const double _temperature) :
	eABF(_outputfile, _outputfreq, _restart, _inputfile, _outputgrad, _gradfreq, _temperature)
{
	lowerboundary.push_back(_lowerboundary);
	lowerboundary.push_back(_lowerboundary2);
	upperboundary.push_back(_upperboundary);
	upperboundary.push_back(_upperboundary2);
	width.push_back(_width);
	width.push_back(_width2);
	krestr.push_back(_krestr);
	krestr.push_back(_krestr2);
	this->initialize();
	this->restart = _restart;
    if(restart == true)
        this->readfile();
}

bool eABF2D::initialize()
{
	for (int i = 0; i < 2; i++)
	{
		min.push_back(eabffunc::doubletoint(lowerboundary[i] / width[i]) - 10);
		max.push_back(eabffunc::doubletoint(upperboundary[i] / width[i]) + 10);
		bins.push_back(max[i] - min[i]);
	}

    // initialize the distribution of "countall"
    // be very careful about it
	countall = new int*[bins[0] * bins[0]];
	for (int i = 0; i < bins[0] * bins[0]; i++)
	{
		countall[i] = new int[bins[1] * bins[1]];
		for (int j = 0; j < bins[1] * bins[1]; j++)
			countall[i][j] = 0;
    }

    // initialize other variables
	for (int i = 0; i < bins[0]; i++)
	{
		sumx1.push_back(std::vector<double>(bins[1], 0));
		sumx21.push_back(std::vector<double>(bins[1], 0));
		sumx2.push_back(std::vector<double>(bins[1], 0));
		sumx22.push_back(std::vector<double>(bins[1], 0));
		county.push_back(std::vector<int>(bins[1], 0));
	}
    line = 0;

    return true;
}

bool eABF2D::update(const std::string& colvarstring)
{
    // this is the string obtained from colvar each step
    // the format of this string is the same as those in colvar.traj file

	std::vector<std::string> data;
	eabffunc::split(colvarstring, data);

    int step;

    double abf_force;
    double eabf_force;
    double abf_force2;
    double eabf_force2;

	step = eabffunc::chartoint(data[0]);
	abf_force = eabffunc::chartodouble(data[col[0]]);
	eabf_force = eabffunc::chartodouble(data[col[1]]);
	abf_force2 = eabffunc::chartodouble(data[col[2]]);
	eabf_force2 = eabffunc::chartodouble(data[col[3]]);


	if (abf_force > 150 && eabf_force < -150)
		eabf_force += 360;
	else if (abf_force < -150 && eabf_force > 150)
		eabf_force -= 360;

	if (abf_force2 > 150 && eabf_force2 < -150)
		eabf_force2 += 360;
	else if (abf_force2 < -150 && eabf_force2 > 150)
        eabf_force2 -= 360;

	int binx = int(floor(abf_force / width[0])) - min[0];
	int biny = int(floor(eabf_force / width[0])) - min[0];
	int binx2 = int(floor(abf_force2 / width[1])) - min[1];
	int biny2 = int(floor(eabf_force2 / width[1])) - min[1];

	if (step%outputfreq == 0)
    {
        writefile();
    }

	if (outputgrad == true && (step%gradfreq) == 0)
    {
        calpmf();
    }


	if (binx < 0 || binx >= bins[0] || biny < 0 || biny >= bins[0] || binx2 < 0 || binx2 >= bins[1] || biny2 < 0 || biny2 >= bins[1])
		return false;
    else
    {
        // distribution
		countall[binx*bins[0] + biny][binx2*bins[1] + biny2]++;

        // if this is the first time add "countall", then the line in outputfile will increase
		if (countall[binx*bins[0] + biny][binx2*bins[1] + biny2] == 1)
			line++;
        county[biny][biny2]++;
    }

    sumx1[biny][biny2] += abf_force;
    sumx21[biny][biny2] += abf_force * abf_force;
    sumx2[biny][biny2] += abf_force2;
	sumx22[biny][biny2] += abf_force2 * abf_force2;


    // TODO: add convergence judgement
    return true;
}


// WARNGING: do not test whether this function works !!!
bool eABF2D::readfile()
{
    std::ifstream file(inputfile.c_str(), std::ios::in);

    // output file format, the first line is "this->line"
	int temp_line, temp_bins, temp_bins2;
	double temp_lowerboundary, temp_width, temp_lowerboundary2, temp_width2;
	file >> temp_line >> temp_lowerboundary >> temp_width >> temp_bins >> krestr[0] >> temp_lowerboundary2 >> temp_width2 >> temp_bins2 >> krestr[1] >> temperature;


	temp_bins = temp_bins + 20;
	temp_bins2 = temp_bins2 + 20;

	int x = 0, y = 0, m = 0, n = 0;
	int pos0 = 0, pos1 = 0, temp_countall = 0;
	for (int i = 0; i < temp_line; i++)
    {
		file >> x >> y >> m >> n;
		//if (x > temp_bins - 10 || x < 10 || y > temp_bins2 - 10 || y < 10)
		//	continue;
		file >> temp_countall;
		;
		pos0 = convertscale(temp_lowerboundary, x, 1) * bins[0] + convertscale(temp_lowerboundary, y, 1);
		pos1 = convertscale(temp_lowerboundary2, m, 2) * bins[1] + convertscale(temp_lowerboundary2, n, 2);
		// find the correct place to increase countall
		if (countall[pos0][pos1] == 0)
			line++;
		countall[pos0][pos1] += temp_countall;
    }

	double temp_sumx1, temp_sumx21, temp_sumx2, temp_sumx22;
	int temp_county;
	for (int i = 0; i < temp_bins; i++)
    {
		for (int j = 0; j < temp_bins2; j++)
        {
			file >> x >> y >> temp_sumx1 >> temp_sumx21 >> temp_sumx2 >> temp_sumx22 >> temp_county;
			//if (x > temp_bins - 10 || x < 10 || y > temp_bins2 - 10 || y < 10)
			//	continue;
			;
			sumx1[convertscale(temp_lowerboundary, x, 1)][convertscale(temp_lowerboundary2, y, 2)] += temp_sumx1;
			sumx21[convertscale(temp_lowerboundary, x, 1)][convertscale(temp_lowerboundary2, y, 2)] += temp_sumx21;
			sumx2[convertscale(temp_lowerboundary, x, 1)][convertscale(temp_lowerboundary2, y, 2)] += temp_sumx2;
			sumx22[convertscale(temp_lowerboundary, x, 1)][convertscale(temp_lowerboundary2, y, 2)] += temp_sumx22;
			county[convertscale(temp_lowerboundary, x, 1)][convertscale(temp_lowerboundary2, y, 2)] += temp_county;
        }
    }

    file.close();
    return true;
}

bool eABF2D::writefile() const
{
    ofstream_namd file(outputfile.c_str(),".old");
    file<<std::setprecision(15);

	file << line << " " << lowerboundary[0] << " " << width[0] << " " << bins[0] - 20 << " " << krestr[0] << " ";
	file << lowerboundary[1] << " " << width[1] << " " << bins[1] - 20 << " " << krestr[1] << " " << temperature << std::endl;

	for (int i = 0; i < bins[0]; i++)
		for (int j = 0; j < bins[0]; j++)
			for (int m = 0; m < bins[1]; m++)
				for (int n = 0; n < bins[1]; n++)
					if (countall[i*bins[0] + j][m*bins[1] + n] != 0)
						// the number of line should be "this->line", I guess :)
						file << i << " " << j << " " << m << " " << n << " " << countall[i*bins[0] + j][m*bins[1] + n] << std::endl;

	for (int i = 0; i < bins[0]; i++)
		for (int j = 0; j < bins[1]; j++)
			file << i << " " << j << " " << sumx1[i][j] << " " << sumx21[i][j] << " " << sumx2[i][j] << " " << sumx22[i][j] << " " << county[i][j] << " " << std::endl;

    file.close();
    return true;
}


//bool eABF2D::writecount(ofstream_namd& os) const
//{
//	double x1,x2;
//    writehead(os);
//    for(int binx=0; binx<bins; binx++)
//    {
//       x1 = (binx + _min + 0.5) * width;
//        for(int binx2=0; binx2<bins2; binx2++)
//       {
//            x2 = (binx2 + _min2 + 0.5) * width2;
//			if(x1>=lowerboundary && x1<upperboundary && x2>=lowerboundary2 && x2<upperboundary2)
//				os<<x1<<" "<<x2<<" "<<county[binx][binx2]<<std::endl;
//       }
//	}
//    return true;
//}

bool eABF2D::writehead(ofstream_namd& os) const
{
    os<<"# 2"<<std::endl;
	// the "real" bin numbers is "bins - 20"
	os << "# " << lowerboundary[0] << " " << width[0] << " " << (bins[0] - 20) << " " << 0 << std::endl;
	os << "# " << lowerboundary[1] << " " << width[1] << " " << (bins[1] - 20) << " " << 0 << std::endl;
    os<<std::endl;
    return true;
}

bool eABF2D::calpmf() const
{
	// double integration
	// this is just a conversion from Chipot's script
	int norm;
	std::vector<std::vector<double> > x_av(bins[0]);
	std::vector<std::vector<double> > x_av2(bins[0]);
	std::vector<std::vector<double> > sigma21(bins[0]);
	std::vector<std::vector<double> > sigma22(bins[0]);
	for (int i = 0; i < bins[0]; i++)
	{
		x_av[i] = (std::vector<double>(bins[1], 0));
		x_av2[i] = (std::vector<double>(bins[1], 0));
		sigma21[i] = (std::vector<double>(bins[1], 0));
		sigma22[i] = (std::vector<double>(bins[1], 0));
	}
	for (int biny = 0; biny < bins[0]; biny++)
	{
		for (int biny2 = 0; biny2 < bins[1]; biny2++)
		{
			norm = county[biny][biny2] > 0 ? county[biny][biny2] : 1;
			x_av[biny][biny2] = sumx1[biny][biny2] / norm;
			x_av2[biny][biny2] = sumx2[biny][biny2] / norm;
			sigma21[biny][biny2] = sumx21[biny][biny2] / norm - x_av[biny][biny2] * x_av[biny][biny2];
			sigma22[biny][biny2] = sumx22[biny][biny2] / norm - x_av2[biny][biny2] * x_av2[biny][biny2];
		}

	}

	static std::string gridfilename = outputfile + ".grad";
	static std::string histfilename = outputfile + ".hist.grad";
	static std::string countfilename = outputfile + ".count";
	//static std::string pmffilename = outputfile + ".pmf";

	ofstream_namd ofile(gridfilename.c_str(),".old");
	ofstream_namd ofile_hist(histfilename.c_str(), std::ios::app);
	ofstream_namd ofile_count(countfilename.c_str(),".old");

	writehead(ofile);
	writehead(ofile_hist);
	writehead(ofile_count);
	//writecount(ofile_count);
	//ofstream_namd ofile_pmf(pmffilename.c_str(),".old");

	//ofile_pmf<<lowerboundary<<" 0.0"<<std::endl;

	//double pmf = 0, pmf_x = 0;
	double x1, x2, y1, y2, av1, av2, diff_av1, diff_av2, grad1, grad2;
	for (int binx = 0; binx < bins[0]; binx++)
	{
		x1 = (binx + min[0] + 0.5) * width[0];
		for (int binx2 = 0; binx2 < bins[1]; binx2++)
		{
			x2 = (binx2 + min[1] + 0.5) * width[1];
			norm = 0;
			av1 = 0;
			av2 = 0;
			diff_av1 = 0;
			diff_av2 = 0;

			for (int biny = 0; biny < bins[0]; biny++)
			{
				y1 = (biny + min[0] + 0.5) * width[0];
				for (int biny2 = 0; biny2 < bins[1]; biny2++)
				{
					y2 = (biny2 + min[1] + 0.5) * width[1];
					norm += countall[binx * bins[0] + biny][binx2 * bins[1] + biny2];

					if (sigma21[biny][biny2]>0.00001 || sigma21[biny][biny2] < -0.00001)
					{
						av1 += countall[binx * bins[0] + biny][binx2 * bins[1] + biny2] * (x1 - x_av[biny][biny2]) / sigma21[biny][biny2];
					}
					else if (countall[binx * bins[0] + biny][binx2 * bins[1] + biny2] > 1)
					{
						// TODOC:\Users\ChinaFu\eabf\eabf2d.cpp|165|error: passing 'const eABF2D' as 'this' argument of 'int eABF2D::position(int, int, int, int)' discards qualifiers [-fpermissive]|, print error
					}
					diff_av1 += countall[binx * bins[0] + biny][binx2 * bins[1] + biny2] * (x1 - y1);

					if (sigma22[biny][biny2]>0.00001 || sigma22[biny][biny2] < -0.00001)
					{
						av2 += countall[binx * bins[0] + biny][binx2 * bins[1] + biny2] * (x2 - x_av2[biny][biny2]) / sigma22[biny][biny2];
					}
					else if (countall[binx * bins[0] + biny][binx2 * bins[1] + biny2] > 1)
					{
						// TODO, print error
					}

					diff_av2 += countall[binx * bins[0] + biny][binx2 * bins[1] + biny2] * (x2 - y2);
				}
			}

			diff_av1 /= (norm > 0 ? norm : 1);
			diff_av2 /= (norm > 0 ? norm : 1);

			av1 = eabffunc::BOLTZMANN * temperature * av1 / (norm > 0 ? norm : 1);
			av2 = eabffunc::BOLTZMANN * temperature * av2 / (norm > 0 ? norm : 1);

			grad1 = av1 - krestr[0] * diff_av1;
			grad2 = av2 - krestr[1] * diff_av2;

			if (x1 >= lowerboundary[0] && x1 < upperboundary[0] && x2 >= lowerboundary[1] && x2 < upperboundary[1])
			{
				ofile << x1 << " " << x2 << " " << grad1 << " " << grad2 << std::endl;
				ofile_hist << x1 << " " << x2 << " " << grad1 << " " << grad2 << std::endl;
				ofile_count << x1 << " " << x2 << " " << norm << std::endl;
				//pmf += grad * width;
				//pmf_x = x + 0.5 * width;
				//ofile_pmf<<pmf_x<<" "<<pmf<<std::endl;
				// TODO, rescale the pmf file to prevent negetive free energy
			}
		}
		if (x1 > lowerboundary[0] && x1 < upperboundary[0])
		{
			ofile << std::endl;
			ofile_count << std::endl;
		}
	}
	ofile.close();
	ofile_hist.close();
	ofile_count.close();
	//ofile_pmf.close();
	return true;
}

#endif // EABF2D
