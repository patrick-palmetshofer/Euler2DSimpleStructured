#pragma once
#include "Cell.h"

class Flux
{
public:
	Flux();
	virtual ~Flux();

	virtual ConservativeVariables calcFlux(PrimitiveVariables &prim_left, PrimitiveVariables &prim_right);
	
};