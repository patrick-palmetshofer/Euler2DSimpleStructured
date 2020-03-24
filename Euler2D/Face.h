#pragma once
#include "Cell.h"
#include <vector>
class Face
{
private:
	double Area;
	std::array<std::array<double, 2>,2> physcoord;
	PrimitiveVariables leftstate;
	PrimitiveVariables rightstate;

public:
	Face();
	Face(double x1, double x2, double y1, double y2);
	~Face();
};

