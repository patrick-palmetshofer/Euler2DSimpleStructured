#pragma once
#include "Limiter.h"
#include "Flux.h"
#include "Fluid.h"

class Config
{
public:
	Config();
	~Config();

	Limiter * getLimiter();
	Flux * getFlux();
	Fluid * getFluid();
};