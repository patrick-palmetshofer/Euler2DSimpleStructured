#pragma once
#include <valarray>
#include "Fluid.h"
#include "Face.h"

class Euler2DVector :
	public std::valarray<double>
{};

class PrimitiveVariables :
	public Euler2DVector
{
public:
	double &rho = (*this)[0];
	double &u = (*this)[1];
	double &v = (*this)[2];
	double &e = (*this)[3];
};

class ConservativeVariables :
	public Euler2DVector
{
public:
	double &rho = (*this)[0];
	double &rhou = (*this)[1];
	double &rhov = (*this)[2];
	double &rhoet = (*this)[3];
};

class Cell
{
protected:
	double volume;

	Fluid *fluid;

	std::vector<Face*> faces;

	Cell *left;
	Cell *right;
	Cell *up;
	Cell *down;


	// Primitive variables consisting of rho, u, v, et
	PrimitiveVariables primitive;
	// Conservative variables consisting of rho, rho*u, rho*v, rho*et
	ConservativeVariables conservative;

	void Prim2Cons();
	void Cons2Prim();

public:
	Cell();
	~Cell();

	void setPrimitive(PrimitiveVariables new_primitive);
	void setConservative(ConservativeVariables new_conservative);

	PrimitiveVariables getPrimitive();
	ConservativeVariables getConservative();
	Fluid * getFluid();
};

