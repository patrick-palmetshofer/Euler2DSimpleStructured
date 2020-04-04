#pragma once
#include "PrimitiveVariables.h"

class Flux
{
public:
	Flux();
	virtual ~Flux();

	virtual ConservativeVariables calcFlux(PrimitiveVariables &prim_left, PrimitiveVariables &prim_right);
	
};