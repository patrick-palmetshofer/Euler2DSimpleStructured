#pragma once
class Fluid
{
public:
	Fluid();
	virtual ~Fluid();
	//virtual double pEOS(rho,)

	virtual double getSoundSpeed(PrimitiveVariables &prim);
	virtual double getPressure(PrimitiveVariables &prim);
	virtual double getTemperature(PrimitiveVariables &prim);
	virtual double getEnthalpy(PrimitiveVariables &prim);
	virtual double getGamma(PrimitiveVariables &prim);
};

