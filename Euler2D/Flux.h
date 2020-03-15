#pragma once
#include "Cell.h"
class Flux
{
protected:
	Cell *leftCell;
	Cell *rightCell;
public:
	Flux();
	double calcFlux();
	virtual ~Flux();
};

