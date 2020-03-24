#pragma once
#include "Cell.h"
#include <string>
#include <vector>
class Grid
{
protected:
	std::vector<std::vector<std::array<double,2>>> points;
	std::vector<std::vector<Cell>> cells;
public:
	Grid();
	~Grid();

	void readGrid(std::string filename);
};

