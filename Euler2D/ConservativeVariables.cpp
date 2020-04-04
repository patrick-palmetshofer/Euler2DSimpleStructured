#include "ConservativeVariables.h"

ConservativeVariables ConservativeVariables::operator=(ConservativeVariables & p)
{
	ConservativeVariables v = p;
	return v;
}

ConservativeVariables ConservativeVariables::operator=(PrimitiveVariables & p)
{
	ConservativeVariables c;
	c.rho = p.rho;
	c.rhou = p.rho*p.u;
	c.rhov = p.rho*p.v;
	double kinetic = 0.5*p.rho*(p.u*p.u + p.v*p.v);
	c.rhoet = p.rho*p.e + kinetic;
	return c;
}

ConservativeVariables ConservativeVariables::operator*(double scalar)
{
	ConservativeVariables v = *(this);
	v *= scalar;
	return v;
}
