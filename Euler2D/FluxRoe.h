#pragma once
#include "Flux.h"

class RoeAverages
{
private:
	PrimitiveVariables prim;
public:
	double u;
	double v;
	double H;

	RoeAverages();
	RoeAverages(Cell *cell);
};

class FluxRoe :
	public Flux
{
protected:
	Fluid * fluid;
public:
	FluxRoe();
	~FluxRoe();

	void calcRoeAverages(Cell *cell);

	double calcFlux(Cell *leftCell, Cell *rightCell);
};

