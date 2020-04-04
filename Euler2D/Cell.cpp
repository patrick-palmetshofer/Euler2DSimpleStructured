#include "Cell.h"

Cell::Cell(double volume, Face * xi_neg, Face * xi_pos, Face * eta_neg, Face * eta_pos, Config & cfg)
{
	this->volume = volume;
	this->xi_neg = xi_neg;
	this->xi_pos = xi_pos;
	this->eta_neg = eta_neg;
	this->eta_pos = eta_pos;

	limiter = cfg.getLimiter();
	fluid = cfg.getFluid();
}

Cell::~Cell()
{
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

void Cell::setNeighbors(std::shared_ptr<Cell> left, std::shared_ptr<Cell> right, std::shared_ptr<Cell> up, std::shared_ptr<Cell> down)
{
	this->left = left;
	this->right = right;
	this->up = up;
	this->down = down;
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

std::shared_ptr<Cell> Cell::getLeftCell()
{
	return left;
}

std::shared_ptr<Cell> Cell::getRightCell()
{
	return right;
}

std::shared_ptr<Cell> Cell::getUpCell()
{
	return up;
}

std::shared_ptr<Cell> Cell::getDownCell()
{
	return down;
}

ConservativeVariables Cell::getSumFaceFluxes()
{
	return xi_pos->getIntFlux() - xi_neg->getIntFlux() + eta_pos->getIntFlux() - eta_neg->getIntFlux();
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

double Cell::getVolume()
{
	return volume;
}

Fluid * Cell::getFluid()
{
	return fluid;
}