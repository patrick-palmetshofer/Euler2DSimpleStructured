#include "Grid.h"

#include <fstream>
#include <stdexcept>
#include <iostream>

Grid::Grid()
{
}


Grid::~Grid()
{
}

double Grid::getnComponent(int i, int j, int compdim, int physdim)
{
	if (compdim == 0)
	{
		if (physdim == 0)
			return getnXiXs(i, j);
		else if (physdim == 1)
			return getnXiYs(i, j);
	}
	else if (compdim == 1)
	{
		if (physdim == 0)
			return getnEtaXs(i, j);
		else if (physdim == 1)
			return getnEtaYs(i, j);
	}
	throw;
	return std::numeric_limits<double>::infinity();
}

//Reads a GridPro grid file. Structured grid.
void Grid::readGridPro(std::string filename)
{
	std::ifstream stream;

	//Number of points in both computational directions
	int nxi_points;
	int neta_points;
	//std::vector<std::vector<std::array<double,2>>> points;

	try
	{
		stream.open(filename, std::ifstream::in);
		if (stream.fail())
			throw;
		stream >> nxi_points;
		stream >> neta_points;
		int zpoints;
		stream >> zpoints;

		nxi_cells = nxi_points - 1;
		neta_cells = neta_points - 1;

		points.resize(nxi_points);
		for (int i = 0; i < nxi_points; i++)
		{
			points[i].resize(neta_points);
		}
		double zcord;


		for (int i = 0; i < nxi_points; i++)
		{
			for (int j = 0; j < neta_points; j++)
			{
				stream >> points[i][j][0];
				stream >> points[i][j][1];
				stream >> zcord;
			}
		}
		stream.close();
	}
	catch (std::ifstream::failure e) {
		std::cerr << "Exception reading file\n";
	}


	//Allocate all matrices
	nxi_xs.resize(nxi_cells + 1, std::vector<double>(neta_cells));
	nxi_ys.resize(nxi_cells + 1, std::vector<double>(neta_cells));
	Sxis.resize(nxi_cells + 1, std::vector<double>(neta_cells));

	neta_xs.resize(nxi_cells, std::vector<double>(neta_cells+1));
	neta_ys.resize(nxi_cells, std::vector<double>(neta_cells + 1));
	Setas.resize(nxi_cells, std::vector<double>(neta_cells + 1));

	volumes.resize(nxi_cells, std::vector<double>(neta_cells));


	//Fill matrices
	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells + 1; j++)
		{
			double Setax = -(points[i + 1][j][1] - points[i][j][1]);
			double Setay = (points[i + 1][j][0] - points[i][j][0]);

			double Seta = std::sqrt(Setax*Setax + Setay * Setay);

			neta_xs[i][j] = Setax / Seta;
			neta_ys[i][j] = Setay / Seta;

			Setas[i][j] = Seta;
		}
	}
	for (int i = 0; i < nxi_cells + 1; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			double Sxix = points[i][j + 1][1] - points[i][j][1];
			double Sxiy = -(points[i][j + 1][0] - points[i][j][0]);

			double Sxi = std::sqrt(Sxix*Sxix + Sxiy * Sxiy);

			nxi_xs[i][j] = Sxix / Sxi;
			nxi_ys[i][j] = Sxiy / Sxi;

			Sxis[i][j] = Sxi;
		}
	}

	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			volumes[i][j] = std::abs(0.5*((points[i + 1][j + 1][0] - points[i][j][0])*(points[i][j + 1][1] - points[i + 1][j][1]) - (points[i][j + 1][0] - points[i + 1][j][0])*(points[i + 1][j + 1][1] - points[i][j][1])));
		}
	}

	double xmax = 0, xmin = 0, ymax = 0, ymin = 0;
	for (auto& row : points)
	{
		for (auto& elem : row)
		{
			xmax = std::max(xmax, elem[0]);
			xmin = std::min(xmin, elem[0]);
			ymax = std::max(ymax, elem[1]);
			ymin = std::min(ymin, elem[1]);

		}
	}
	maxDistance = { xmax-xmin,ymax-ymin };
}

