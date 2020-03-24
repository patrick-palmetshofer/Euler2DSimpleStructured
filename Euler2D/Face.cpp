#include "Face.h"



Face::Face()
{
}

Face::Face(double x1, double x2, double y1, double y2)
{
	physcoord[0][0] = x1;
	physcoord[0][1] = x2;
	physcoord[1][0] = y1;
	physcoord[1][1] = y2;
}


Face::~Face()
{
}
