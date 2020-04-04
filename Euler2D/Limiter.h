#pragma once
#include "PrimitiveVariables.h"
#include <array>

class Limiter
{
public:
	Limiter();
	virtual ~Limiter();
	virtual PrimitiveVariables calcPhis(PrimitiveVariables &rs) const;
};