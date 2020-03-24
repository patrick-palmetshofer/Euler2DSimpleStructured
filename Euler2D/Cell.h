#pragma once
#include <valarray>
#include "Fluid.h"
#include "Face.h"
#include "Flux.h"
#include "Limiter.h"
#include "Config.h"

class Euler2DVector :
	public std::valarray<double>
{
public:
	using std::valarray<double>::valarray;
	using std::valarray<double>::operator*=;
	using std::valarray<double>::operator=;

	Euler2DVector operator* (double scalar)
	{
		Euler2DVector v = *(this);
		v *= scalar;
		return v;
	}	
};

class PrimitiveVariables :
	public Euler2DVector
{
public:
	using std::valarray<double>::valarray;
	using std::valarray<double>::operator*=;
	using Euler2DVector::Euler2DVector;

	double &rho = (*this)[0];
	double &u = (*this)[1];
	double &v = (*this)[2];
	double &e = (*this)[3];

	PrimitiveVariables operator= (PrimitiveVariables p)
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
	}
};

class ConservativeVariables :
	public Euler2DVector
{
public:
	using std::valarray<double>::valarray;
	using std::valarray<double>::operator*=;
	using Euler2DVector::Euler2DVector;
	using Euler2DVector::operator*;

	double &rho = (*this)[0];
	double &rhou = (*this)[1];
	double &rhov = (*this)[2];
	double &rhoet = (*this)[3];

	ConservativeVariables operator= (ConservativeVariables p)
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

	ConservativeVariables operator* (double scalar)
	{
		ConservativeVariables v = *(this);
		v *= scalar;
		return v;
	}
};

class Cell
{
protected:
	double volume;

	Fluid *fluid;
	Flux *flux;

	Limiter *limiter;

	Face* leftface;
	Face* rightface;
	Face* upface;
	Face* downface;

	Cell *left;
	Cell *right;
	Cell *up;
	Cell *down;

	//NOTE: Cell variables are in computational space!
	// Primitive variables consisting of rho, u, v, et
	PrimitiveVariables primitive;
	// Conservative variables consisting of rho, rho*u, rho*v, rho*et
	ConservativeVariables conservative;

	//PrimitiveVariables xslope;
	//PrimitiveVariables yslope;

	double reconstructepsilon;
	double reconstructkappa;

	PrimitiveVariables leftstate;
	PrimitiveVariables upstate;
	PrimitiveVariables rightstate;
	PrimitiveVariables downstate;

	ConservativeVariables rightflux;
	ConservativeVariables downflux;

	template<class T> T swapUV(T oldT)
	{
		T newT;
		newT[1] = oldT[2];
		newT[2] = oldT[1];
		return newT;
	}

public:
	Cell(Config &cfg);
	Cell(double volume,
	double Sxi,
	double Seta,
	double Sxix,
	double Sxiy,
	double Setax,
	double Setay, Config &cfg);
	~Cell();

	void setNeighbors(Cell* leftcell, Cell* rightcell, Cell* upcell, Cell* downcell);
	void setPrimitive(PrimitiveVariables new_primitive);
	void setConservative(ConservativeVariables new_conservative);

	PrimitiveVariables *getPrimitive();
	ConservativeVariables *getConservative();

	void reconstructStates();

	PrimitiveVariables *getLeftState();
	PrimitiveVariables *getRightState();
	PrimitiveVariables *getUpState();
	PrimitiveVariables *getDownState();

	ConservativeVariables *getLeftFlux();
	ConservativeVariables *getUpFlux();
	ConservativeVariables *getRightFlux();
	ConservativeVariables *getDownFlux();

	Cell *getLeftCell();
	Cell *getRightCell();
	Cell *getUpCell();
	Cell *getDownCell();

	Fluid * getFluid();
	double getCdeltat();

	void calcRightFlux();
	void calcDownFlux();
};
