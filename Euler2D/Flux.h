#pragma once
#include "Cell.h"

class Flux
{
protected:
public:
	Flux();
	virtual ~Flux();

	virtual double calcFlux(Cell *leftCell, Cell *rightCell);
};