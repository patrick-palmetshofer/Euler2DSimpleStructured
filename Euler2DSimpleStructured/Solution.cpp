#include "Solution.h"
#include <fstream>
#include <iostream>

Solution::Solution()
{
}


Solution::~Solution()
{
}

//Residuals calculation. Either returning the maximum residual or residual in each variable.
double Solution::calcResidualL2(StateMatrix2D &o, StateMatrix2D &n)
{
	double sumres = 0;
	for (int i = 0; i < n.size(); i++)
	{
		for (int j = 0; j < n[0].size(); j++)
		{
			for (int k = 0; k < 4; k++)
			{
				double normk = o[0][0][k];
				if (normk == 0)
					normk = o[0][0][1];
				double temp = (n[i][j][k] - o[i][j][k]) / normk;
				sumres += temp * temp;
			}
		}
	}
	if (sumres < 0)
		throw;
	sumres = std::sqrt(sumres);
	return sumres;
}

StateVector2D Solution::calcResidualsL2(StateMatrix2D &o, StateMatrix2D &n)
{
	StateVector2D res;
	for (int k = 0; k < n[0][0].size(); k++)
	{
		double sumres = 0;
		double normk = o[0][0][k];
		if (normk == 0)
			normk = o[0][0][1];

		for (int i = 0; i < n.size(); i++)
		{
			for (int j = 0; j < n[0].size(); j++)
			{
				double temp = (n[i][j][k] - o[i][j][k]) / normk;
				sumres += temp * temp;
			}
		}

		if (sumres < 0)
			throw;

		sumres = std::sqrt(sumres);
		res[k] = sumres;
	}
	return res;
}

double Solution::calcResidualLinfty(StateMatrix2D &o, StateMatrix2D &n)
{
	double res = 0;
	int imax, jmax, kmax;
	for (int i = 0; i < n.size(); i++)
	{
		for (int j = 0; j < n[0].size(); j++)
		{
			for (int k = 0; k < n[0][0].size(); k++)
			{
				double normk = o[0][0][k];
				if (normk == 0)
					normk = o[0][0][1];
				double temp = std::abs((n[i][j][k] - o[i][j][k]) / normk);
				if (temp > res)
				{
					imax = i;
					jmax = j;
					kmax = k;
					res = temp;
				}
			}
		}
	}
	return res;
}

StateVector2D Solution::calcResidualsLinfty(StateMatrix2D &o, StateMatrix2D &n)
{
	StateVector2D res;
	int imax, jmax, kmax;
	for (int k = 0; k < 4; k++)
	{
		res[k] = 0;
		double normk = o[0][0][k];
		if (normk == 0)
			normk = o[0][0][1];
		for (int i = 0; i < n.size(); i++)
		{
			for (int j = 0; j < n[0].size(); j++)
			{
				double temp = std::abs((n[i][j][k] - o[i][j][k]) / normk);
				if (temp > res[k])
				{
					imax = i;
					jmax = j;
					kmax = k;
					res[k] = temp;
				}
			}
		}
	}
	return res;
}


StateVector2D Solution::getResidualsL2()
{
	return residuals_L2;
}

StateVector2D Solution::getResidualsLinfty()
{
	return residuals_Linfty;
}

void Solution::writeSolution(std::string filename, StateMatrix2D p)
{
	//Write a set of primitive variables, with Temperature in place of e_t
	for (int i = 0; i < p.size(); i++)
	{
		for (int j = 0; j < p[0].size(); j++)
		{
			p[i][j] = fluid->cons2prim(p[i][j]);
			p[i][j][3] = (p[i][j][3] - 0.5*(p[i][j][2] * p[i][j][2] + p[i][j][1] * p[i][j][1]))*fluid->getGamma() / fluid->getCp();
		}
	}

	//Write to file
	std::ofstream stream;
	try
	{
		stream.open(filename, std::ofstream::out);
		stream << p.size() << ",\t" << p[0].size() << "\n";
		stream << "rho" << ",\t" << "u" << ",\t" << "v" << ",\t" << "T" << "\n";
		//Reserve add here
		for (int i = 0; i < p.size() ; i++)
		{
			for (int j = 0; j < p[0].size(); j++)
			{
				for (int k = 0; k < 3; k++)
				{
					stream << p[i][j][k] << ",\t";
				}
				stream << p[i][j][3] << "\n";
			}
		}
		stream.flush();
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}

void writeResidual(std::vector<StateVector2D> &residuals, std::string &filename)
{
	std::ofstream stream;
	try
	{
		stream.open(filename, std::ofstream::out);
		stream << "rho" << ",\t" << "rhou" << ",\t" << "rhov" << ",\t" << "rhoet" << "\n";
		//Reserve add here
		for (int i = 0; i < residuals.size(); i++)
		{
			for (int k = 0; k < 3; k++)
			{
				stream << residuals[i][k] << ",\t";
			}
			stream << residuals[i][3] << "\n";
		}
		stream.flush();
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}