#pragma once
class Fluid
{
protected:
	double density;
	double pressure;
	double temperature;
	double soundspeed;
	double internalenergy;
	double enthalpy;

public:
	Fluid();
	virtual ~Fluid();
	//virtual double pEOS(rho,)

	virtual double getSoundSpeed();
	virtual double getPressure();
	virtual double getTemperature();
	virtual double getEnthalpy();

	virtual double getGamma();
};

