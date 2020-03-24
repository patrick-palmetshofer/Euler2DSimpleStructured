#include "Grid.h"
#include <cmath>


Grid::Grid()
{
}


Grid::~Grid()
{
}

void Grid::readGrid(std::string filename)
{
	double i;
	double j;
	double volume = std::abs(0.5*((points[i + 1][j + 1][0] - points[i][j][0])*(points[i][j + 1][1] - points[i+1][j][1]) - (points[i][j + 1][0] - points[i+1][j][0])*(points[i + 1][j + 1][1] - points[i][j][1])));

}
