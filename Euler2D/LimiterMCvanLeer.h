#pragma once
#include "Limiter.h"
#include <array>
class LimiterMCvanLeer :
	public Limiter
{
protected:
public:
	LimiterMCvanLeer();
	~LimiterMCvanLeer();
	PrimitiveVariables calcPhis(PrimitiveVariables &rs) const;
};

