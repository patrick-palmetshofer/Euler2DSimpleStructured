#include "PrimitiveVariables.h"

PrimitiveVariables PrimitiveVariables::operator=(PrimitiveVariables & p)
{
	PrimitiveVariables v = p;
	return v;
}

PrimitiveVariables PrimitiveVariables::operator*(double scalar)
{
	PrimitiveVariables v = *(this);
	v *= scalar;
	return v;
}

PrimitiveVariables PrimitiveVariables::operator=(ConservativeVariables & c)
{
	PrimitiveVariables p;
	p.rho = c.rho;
	p.u = c.rhou / p.rho;
	p.v = c.rhov / p.rho;
	double kinetic = 0.5*p.rho*(p.u*p.u + p.v*p.v);
	p.e = (c.rhoet - kinetic) / p.rho;
	return p;
}
