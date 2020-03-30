#pragma once
#include <vector>
#include "Grid.h"
#include "Cell.h"

class TimeStepper
{
protected:
	std::vector<std::shared_ptr<Cell>> cells;
	std::vector<std::shared_ptr<Face>> faces;
	std::vector<std::shared_ptr<GhostCell>> boundaryCells;
public:
	TimeStepper();
	virtual ~TimeStepper();
};

