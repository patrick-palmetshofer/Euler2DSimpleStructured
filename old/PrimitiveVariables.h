#pragma once
#include <valarray>
#include "DeclareConsPrim.h"
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

	PrimitiveVariables operator= (PrimitiveVariables &p);


	PrimitiveVariables operator* (double scalar);


	PrimitiveVariables operator= (ConservativeVariables &c);

};