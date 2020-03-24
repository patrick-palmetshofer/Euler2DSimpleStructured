#pragma once
#include "Limiter.h"
class LimiterMinmod :
	public Limiter
{
public:
	LimiterMinmod();
	~LimiterMinmod();

	PrimitiveVariables calcPhis(PrimitiveVariables &rs) const;
};

