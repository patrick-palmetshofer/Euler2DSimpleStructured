#pragma once
#include "TimeStepper.h"
class TimeStepperExplicit :
	public TimeStepper
{
public:
	TimeStepperExplicit();
	~TimeStepperExplicit();

	void executeTimeStep();
};

