#pragma once
#include "Flux.h"
#include "Fluid.h"

class RoeAverages
{
private:
	PrimitiveVariables prim;
public:
	double u;
	double v;
	double H;

	RoeAverages();
};

class FluxRoe :
	public Flux
{
protected:
	Fluid * fluid;
	ConservativeVariables calcPhysicalFlux(PrimitiveVariables &prim);
public:
	FluxRoe();
	FluxRoe(Fluid * newfluid);
	~FluxRoe();

	ConservativeVariables calcFlux(PrimitiveVariables &prim_left, PrimitiveVariables &prim_right);
};

