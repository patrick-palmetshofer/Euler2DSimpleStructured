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

double PerfectGas::getSoundSpeed(StateVector2D &prim)
{
	double T = getTemperature(prim);
	return std::sqrt(gamma*R*T);
}

double PerfectGas::getPressure(StateVector2D &prim)
{
	return prim[0]*R*getTemperature(prim);
}

double PerfectGas::getTemperature(StateVector2D &prim)
{
	return prim[3]/cv;
}

double PerfectGas::getEnthalpy(StateVector2D &prim)
{
	return gamma*prim[3];
}

double PerfectGas::getGamma()
{
	return gamma;
}
