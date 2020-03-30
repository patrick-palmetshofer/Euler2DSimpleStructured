#pragma once
#include <valarray>
#include "ConservativeVariables.h"

class PrimitiveVariables :
	public std::valarray<double>
{
public:
	using std::valarray<double>::valarray;
	using std::valarray<double>::operator*=;

	double &rho = (*this)[0];
	double &u = (*this)[1];
	double &v = (*this)[2];
	double &e = (*this)[3];

	PrimitiveVariables operator= (PrimitiveVariables &p)
	{
		PrimitiveVariables v = p;
		return v;
	}

	PrimitiveVariables operator* (double scalar)
	{
		PrimitiveVariables v = *(this);
		v *= scalar;
		return v;
	}

	PrimitiveVariables operator= (ConservativeVariables &c)
	{
		PrimitiveVariables p;
		p.rho = c.rho;
		p.u = c.rhou / p.rho;
		p.v = c.rhov / p.rho;
		double kinetic = 0.5*p.rho*(p.u*p.u + p.v*p.v);
		p.e = (c.rhoet - kinetic) / p.rho;
		return p;
	}
};