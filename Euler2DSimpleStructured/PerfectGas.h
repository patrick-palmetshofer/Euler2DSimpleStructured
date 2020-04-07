#pragma once
#include "GlobalTypes.h"
class PerfectGas
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

	double getSoundSpeed(StateVector2D &prim);
	double getPressure(StateVector2D &prim);
	double getTemperature(StateVector2D &prim);
	double getEnthalpy(StateVector2D &prim);
	double getGamma();
};

