#pragma once
#include <vector>
#include <memory>
#include "Grid.h"
#include "Cell.h"

class Solution
{
private:
	std::vector<std::shared_ptr<Cell>> cells;
public:
	Solution(Grid &grid);
	~Solution();

	void write(std::string filename);
};

