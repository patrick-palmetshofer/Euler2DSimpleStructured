#include "Cell.h"

void Cell::Prim2Cons()
{
	conservative.rho = primitive.rho;
	conservative.rhou = primitive.rho*primitive.u;
	conservative.rhov = primitive.rho*primitive.v;
	double kinetic = 0.5*primitive.rho*(conservative.rhou*primitive.u + conservative.rhov*primitive.v);
	conservative.rhoet = primitive.rho*primitive.e+kinetic;
}

void Cell::Cons2Prim()
{
	primitive.rho = conservative.rho;
	primitive.u = conservative.rhou / primitive.rho;
	primitive.v = conservative.rhov / primitive.rho;
	double kinetic = 0.5*primitive.rho*(conservative.rhou*primitive.u + conservative.rhov*primitive.v);
	primitive.e = (conservative.rhoet - kinetic) / primitive.rho;
}

Cell::Cell()
{
}


Cell::~Cell()
{
}

void Cell::setPrimitive(PrimitiveVariables new_primitive)
{
	primitive = new_primitive;
	Prim2Cons();
}

void Cell::setConservative(ConservativeVariables new_conservative)
{
	conservative = new_conservative;
	Cons2Prim();
}

PrimitiveVariables Cell::getPrimitive()
{
	return primitive;
}

ConservativeVariables Cell::getConservative()
{
	return conservative;
}

Fluid * Cell::getFluid()
{
	return fluid;
}
