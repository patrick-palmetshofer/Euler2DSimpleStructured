#pragma once
#include "Fluid.h"
#include "Cell.h"
class PerfectGas :
	public Fluid
{
private:
	const double Rm = 8.31446261815324;
	double W;
	double cp;
	double R;
	double cv;
	double gamma;
public:
	PerfectGas(double W, double cp);
	~PerfectGas();

	double getSoundSpeed(PrimitiveVariables &prim);
	double getPressure(PrimitiveVariables &prim);
	double getTemperature(PrimitiveVariables &prim);
	double getEnthalpy(PrimitiveVariables &prim);
	double getGamma(PrimitiveVariables &prim);
};

