#pragma once
#include <vector>
#include "Grid.h"
#include "Cell.h"

class TimeStepper
{
protected:
	Grid grid;
	std::vector<std::vector<Cell>> structuredCells;
public:
	TimeStepper();
	virtual ~TimeStepper();
};

