#include "Grid.h"
#include <cmath>
#include <iostream>
#include <fstream>


Grid::Grid()
{
}


Grid::~Grid()
{
}

void Grid::readGrid(std::string filename)
{
	std::ifstream stream;
	try
	{
		stream.open(filename, std::ifstream::in);
		stream >> nxi_points;
		stream >> neta_points;


		std::string check;
		stream >> check;
		if (check != "NDIME")
		{
			std::cerr << "File Format Invalid: NDIME not found\n";
			throw std::exception();
		}
		int dim;
		stream >> dim;
		if (dim != 2)
		{
			std::cerr << "File Format Invalid: Dimension not 2\n";
			throw std::exception();
		}

		


		double i;
		double j;

		while (!stream.eof())
		{
			stream >> i;
			stream >> j;
			stream >> points[i][j][0];
			stream >> points[i][j][1];

		}

		////Reserve add here
		//for (int i = 0; i < nxi_points; i++)
		//{
		//	for (int j = 0; j < neta_points; j++)
		//	{
		//		for (int d = 0; d < dim; d++)
		//		{
		//			stream >> points[i][j][d];
		//		}
		//	}
		//}
		stream.close();
	}
	catch (std::ifstream::failure e) {
		std::cerr << "Exception reading file\n";
	}
}

void Grid::constructCells()
{

	double i;
	double j;
	//Reserve one more cell than points

	for (int i = 0; i < nxi_points-1; i++)
	{
		for (int j = 0; j < neta_points-1; j++)
		{
			double Setay = points[i + 1][j][0] - points[i][j][0];
			double Sxiy = -(points[i][j + 1][0] - points[i][j][0]);
			double Setax = -(points[i + 1][j][1] - points[i][j][1]);
			double Sxix = points[i][j + 1][1] - points[i][j][1];

			double Sxi = std::sqrt(Sxix*Sxix + Sxiy * Sxiy);
			double Seta = std::sqrt(Setax*Setax + Setay * Setay);

			xi_faces[i][j] = std::make_shared<Face>(Sxi, Sxix, Sxiy, *cfg);
			eta_faces[i][j] = std::make_shared<Face>(Seta, Setax, Setay, *cfg);
		}
	}

	for (int i = 0; i < nxi_points - 1; i++)
	{
		for (int j = 0; j < neta_points - 1; j++)
		{
			double volume = std::abs(0.5*((points[i + 1][j + 1][0] - points[i][j][0])*(points[i][j + 1][1] - points[i + 1][j][1]) - (points[i][j + 1][0] - points[i + 1][j][0])*(points[i + 1][j + 1][1] - points[i][j][1])));
			cells[i][j] = std::make_shared<Cell>(volume, &xi_faces[i][j], &xi_faces[i + 1][j], &eta_faces[i][j], &eta_faces[i][j + 1], cfg);
		}
	}

	for (int i = 0; i < nxi_points - 1; i++)
	{
		for (int j = 0; j < neta_points - 1; j++)
		{
			cells[i][j]->setNeighbors(cells[i-1][j]);
		}
	}
}