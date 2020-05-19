#include "TimeStepperGlobal.h"

TimeStepperGlobal::TimeStepperGlobal()
{
}


TimeStepperGlobal::~TimeStepperGlobal()
{
}


void TimeStepperGlobal::execute(StateMatrix2D * conservative, StateMatrix2D *xi_fluxes, StateMatrix2D *eta_fluxes)
{
	//Calculate time step with CFL
	dt = maxCFL * calcTimeStep(conservative);

	//checkNaN(xi_fluxes);
	//checkNaN(eta_fluxes);

	//Execute time step
	StateVector2D Dxi, Deta;
	for (int i = 0; i < conservative->size(); i++)
	{
		for (int j = 0; j < conservative->at(0).size(); j++)
		{
			Dxi = xi_fluxes->at(i + 1).at(j) * grid->getSxi(i + 1, j) - xi_fluxes->at(i).at(j) * grid->getSxi(i, j);
			Deta = eta_fluxes->at(i).at(j + 1) * grid->getSeta(i, j + 1) - eta_fluxes->at(i).at(j) * grid->getSeta(i, j);
			conservative->at(i).at(j) = conservative->at(i).at(j) - dt / grid->getVolume(i, j) * (Dxi + Deta);
		}
	}

	time += dt;
}