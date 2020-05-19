#pragma once
#include "GlobalTypes.h"
#include "Grid.h"
#include "Flux.h"
class TimeStepper
{
protected:
	Grid * grid;
	Fluid * fluid;

	//Calculate time step for global time stepping
	double calcTimeStep(StateMatrix2D * conservative);
	//Method calls:
	//Calculate timestep for local timestepping
	double calcTimeStep(int i, int j, StateMatrix2D * conservative);

	double maxCFL = 0.5;

public:
	TimeStepper();
	virtual ~TimeStepper();

	void setFluid(Fluid * new_fluid) { fluid = new_fluid; };
	void setGrid(Grid * new_grid) { grid = new_grid; };

	virtual void execute(StateMatrix2D * conservative, StateMatrix2D * xi_fluxes, StateMatrix2D * eta_fluxes) = 0;
};

