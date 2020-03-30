#pragma once
#include <valarray>
#include <memory>
#include "Fluid.h"
#include "Face.h"
#include "Flux.h"
#include "Limiter.h"
#include "Config.h"

#include "PrimitiveVariables.h"
#include "ConservativeVariables.h"
#include "Jacobian.h"

class Cell
{
protected:
	double volume;

	Fluid *fluid;
	Flux *flux;

	Limiter *limiter;

	Face* xi_pos;
	Face* xi_neg;
	Face* eta_pos;
	Face* eta_neg;

	std::shared_ptr<Cell> left;
	std::shared_ptr<Cell> right;
	std::shared_ptr<Cell> up;
	std::shared_ptr<Cell> down;

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

	template<class T> T swapUV(T oldT)
	{
		T newT;
		newT[1] = oldT[2];
		newT[2] = oldT[1];
		return newT;
	}

public:
	Cell(double volume, Face * xi_neg, Face * xi_pos, Face * eta_neg, Face * eta_pos, Config & cfg);
	~Cell();

	void setNeighbors(std::shared_ptr<Cell> left, std::shared_ptr<Cell> right, std::shared_ptr<Cell> up, std::shared_ptr<Cell> down);
	void setPrimitive(PrimitiveVariables new_primitive);
	void setConservative(ConservativeVariables new_conservative);

	PrimitiveVariables *getPrimitive();
	ConservativeVariables *getConservative();

	double getVolume();

	void reconstructStates();

	PrimitiveVariables *getLeftState();
	PrimitiveVariables *getRightState();
	PrimitiveVariables *getUpState();
	PrimitiveVariables *getDownState();

	std::shared_ptr<Cell> getLeftCell();
	std::shared_ptr<Cell> getRightCell();
	std::shared_ptr<Cell> getUpCell();
	std::shared_ptr<Cell> getDownCell();

	ConservativeVariables getSumFaceFluxes();

	Fluid * getFluid();
};
