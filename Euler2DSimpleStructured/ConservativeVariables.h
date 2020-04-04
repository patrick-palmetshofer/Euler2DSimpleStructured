#pragma once
#include <valarray>
#include "DeclareConsPrim.h"
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

	ConservativeVariables operator= (ConservativeVariables &p);
	ConservativeVariables operator= (PrimitiveVariables &p);


	//ConservativeVariables()
	//{
	//}

	//ConservativeVariables(PrimitiveVariables &p)
	//{
	//	*this = p;
	//}

	ConservativeVariables operator* (double scalar);

};
