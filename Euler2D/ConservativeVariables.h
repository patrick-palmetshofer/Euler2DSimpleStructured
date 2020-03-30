#pragma once
#include <valarray>
#include "PrimitiveVariables.h"

class ConservativeVariables :
	public std::valarray<double>
{
public:
	using std::valarray<double>::valarray;
	using std::valarray<double>::operator*=;

	double &rho = (*this)[0];
	double &rhou = (*this)[1];
	double &rhov = (*this)[2];
	double &rhoet = (*this)[3];

	ConservativeVariables operator= (ConservativeVariables &p)
	{
		ConservativeVariables v = p;
		return v;
	}

	ConservativeVariables operator= (PrimitiveVariables &p)
	{
		ConservativeVariables c;
		c.rho = p.rho;
		c.rhou = p.rho*p.u;
		c.rhov = p.rho*p.v;
		double kinetic = 0.5*p.rho*(p.u*p.u + p.v*p.v);
		c.rhoet = p.rho*p.e + kinetic;
		return c;
	}

	//ConservativeVariables()
	//{
	//}

	//ConservativeVariables(PrimitiveVariables &p)
	//{
	//	*this = p;
	//}

	ConservativeVariables operator* (double scalar)
	{
		ConservativeVariables v = *(this);
		v *= scalar;
		return v;
	}
};
