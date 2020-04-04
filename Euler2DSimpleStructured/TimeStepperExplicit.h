#pragma once
#include "TimeStepper.h"
#include <memory>
#include "Grid.h"
class TimeStepperExplicit :
	public TimeStepper
{
private:
	const double maxCFL = 0.8;
	double time;

	//ConservativeVariables spaceDisc(Cell &cons);
	ConservativeVariables execute4thOrderRungeKutta(Cell &cell, double dt);
	//ConservativeVariables executeFirstOrderBackward(Cell &cell, double dt);
	ConservativeVariables executeFirstOrderBackward(std::shared_ptr<Cell>& cell, double dt);
	void reconstruct();
	double maxTimeStep();
	void calcFluxes();

public:
	TimeStepperExplicit(std::vector<std::shared_ptr<Cell>> &cells,std::vector<std::shared_ptr<Face>> &faces);
	TimeStepperExplicit(Grid &grid);
	~TimeStepperExplicit();

	ConservativeVariables spaceDisc(std::shared_ptr<Cell>& cell);

	void executeTimeStepGlobal();
	double getTime();
};

