#include "Cell.h"

Cell::Cell(Config &cfg)
{
	Cell(1, 1, 1, 1, 1, 1, 1, cfg);
}

Cell::Cell(double volume, double Sxi, double Seta, double Sxix, double Sxiy, double Setax, double Setay, Config &cfg)
{
	this->volume= volume;


	limiter = cfg.getLimiter();
	flux = cfg.getFlux();
	fluid = cfg.getFluid();
}


Cell::~Cell()
{
}

void Cell::setNeighbors(Cell * leftcell, Cell * rightcell, Cell * upcell, Cell * downcell)
{
	left = leftcell;
	right = rightcell;
	up = upcell;
	down = downcell;
}

void Cell::reconstructStates()
{
	PrimitiveVariables rshoriz = (*right->getPrimitive() - primitive) / (primitive - *left->getPrimitive());
	PrimitiveVariables xphis = limiter->calcPhis(rshoriz);

	PrimitiveVariables rsvert = (*down->getPrimitive() - primitive) / (primitive - *up->getPrimitive());
	PrimitiveVariables yphis = limiter->calcPhis(rsvert);

	rightstate = primitive + 0.25*reconstructepsilon*(primitive - *(left->getPrimitive()))*((1 - reconstructkappa)*xphis + (1 + reconstructkappa)*rshoriz / xphis);
	leftstate = primitive - 0.25*reconstructepsilon*(primitive - *(left->getPrimitive()))*((1 - reconstructkappa)*xphis + (1 + reconstructkappa)*rshoriz / xphis);

	upstate = primitive + 0.25*reconstructepsilon*(primitive - *(up->getPrimitive()))*((1 - reconstructkappa)*xphis + (1 + reconstructkappa)*rsvert / xphis);
	downstate = primitive - 0.25*reconstructepsilon*(primitive - *(up->getPrimitive()))*((1 - reconstructkappa)*xphis + (1 + reconstructkappa)*rsvert / xphis);
}

void Cell::setPrimitive(PrimitiveVariables new_primitive)
{
	primitive = new_primitive;
	conservative = primitive;
}

void Cell::setConservative(ConservativeVariables new_conservative)
{
	conservative = new_conservative;
	primitive = conservative;
}

PrimitiveVariables* Cell::getLeftState()
{
	return &leftstate;
}

PrimitiveVariables *Cell::getDownState()
{
	return &downstate;
}

ConservativeVariables * Cell::getLeftFlux()
{
	return left->getRightFlux();
}

ConservativeVariables * Cell::getUpFlux()
{
	return up->getDownFlux();
}

ConservativeVariables * Cell::getRightFlux()
{
	return &rightflux;
}

ConservativeVariables * Cell::getDownFlux()
{
	return &downflux;
}

Cell * Cell::getLeftCell()
{
	return left;
}

Cell * Cell::getRightCell()
{
	return right;
}

Cell * Cell::getUpCell()
{
	return up;
}

Cell * Cell::getDownCell()
{
	return down;
}

PrimitiveVariables *Cell::getRightState()
{
	return &rightstate;
}

PrimitiveVariables *Cell::getUpState()
{
	return &upstate;
}

PrimitiveVariables* Cell::getPrimitive()
{
	return &primitive;
}

ConservativeVariables* Cell::getConservative()
{
	return &conservative;
}

Fluid * Cell::getFluid()
{
	return fluid;
}

double Cell::getCdeltat()
{
	return primitive.u / xsize + primitive.v / ysize;
}

void Cell::calcRightFlux()
{
	//rightstate is left state of interface
	rightflux = flux->calcFlux(rightstate, *(right->getLeftState()));
}

void Cell::calcDownFlux()
{
	PrimitiveVariables rotatedright = swapUV(primitive);
	PrimitiveVariables rotatedleft = swapUV(*(right->getLeftState()));
	ConservativeVariables rflux = flux->calcFlux(rotatedright, rotatedleft);
	downflux = swapUV(rflux);
}
