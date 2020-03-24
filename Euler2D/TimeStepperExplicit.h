#pragma once
#include "TimeStepper.h"
class TimeStepperExplicit :
	public TimeStepper
{
private:
	const double maxCFL = 0.8;

	ConservativeVariables spaceDisc(Cell &cons);
	ConservativeVariables execute4thOrderRungeKutta(Cell &cell, double dt);
	ConservativeVariables executeFirstOrderBackward(Cell &cell, double dt);
	double reconstructAndMaxCFLGlobal();
	void calcFluxes();
public:
	TimeStepperExplicit();
	~TimeStepperExplicit();

	void executeTimeStepGlobal();
};

