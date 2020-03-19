#pragma once
#include "Fluid.h"
class PerfectGas :
	public Fluid
{
public:
	PerfectGas();
	~PerfectGas();

	double getSoundSpeed();
	double getPressure();
	double getTemperature();
	double getEnthalpy();

	double getGamma();
};

