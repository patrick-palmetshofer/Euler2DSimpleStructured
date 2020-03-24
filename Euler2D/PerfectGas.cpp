#include "PerfectGas.h"
#include <cmath>


//W in g/mol; cp in J/kgK
PerfectGas::PerfectGas(double W,double cp)
{
	this->W = W;
	this->cp = cp;

	R = 1e3*Rm / W;
	cv = cp - R;
	gamma = cp / cv;
}


PerfectGas::~PerfectGas()
{
}

double PerfectGas::getSoundSpeed(PrimitiveVariables & prim)
{
	double T = getTemperature(prim);
	return std::sqrt(gamma*R*T);
}

double PerfectGas::getPressure(PrimitiveVariables & prim)
{
	return prim.rho*R*getTemperature(prim);
}

double PerfectGas::getTemperature(PrimitiveVariables & prim)
{
	return prim.e/cv;
}

double PerfectGas::getEnthalpy(PrimitiveVariables & prim)
{
	return gamma*prim.e;
}

double PerfectGas::getGamma(PrimitiveVariables & prim)
{
	return gamma;
}
