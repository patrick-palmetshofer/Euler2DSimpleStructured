#include "TimeStepperExplicit.h"



TimeStepperExplicit::TimeStepperExplicit()
{
}


TimeStepperExplicit::~TimeStepperExplicit()
{
}

ConservativeVariables TimeStepperExplicit::spaceDisc(Cell &cell)
{
	ConservativeVariables spacedisc = *cell.getRightFlux() - *cell.getLeftFlux() + *cell.getDownFlux() - *cell.getUpFlux();
	return spacedisc;
}

//ConservativeVariables TimeStepperExplicit::execute4thOrderRungeKutta(Cell &cell, double dt)
//{
//	ConservativeVariables consold = *(cell.getConservative());
//	ConservativeVariables k1 = spaceDisc(consold) * dt;
//	ConservativeVariables k2 = spaceDisc(consold+k1 * 0.5)*dt;
//	ConservativeVariables k3 = spaceDisc(consold+k2 * 0.5)*dt;
//	ConservativeVariables k4 = spaceDisc(consold+k3)*dt;
//	ConservativeVariables consnew = consold - (k1 + k2*2.0 + k3*2.0 + k4) / 6.0;
//	return consnew;
//}

ConservativeVariables TimeStepperExplicit::executeFirstOrderBackward(Cell &cell, double dt)
{
	ConservativeVariables consnew = *(cell.getConservative()) - spaceDisc(cell)*dt;
	return consnew;
}

//TODO:: Move all from loops to recursion

//Returns global time step and calculates fluxes in cells
double TimeStepperExplicit::reconstructAndMaxCFLGlobal()
{
	double dt = 1e19;
	for (auto row : *structuredCells)
	{
		for (auto cell : row)
		{
			cell.reconstructStates();
			double dtlocal = maxCFL/cell.getCdeltat();
			if (dtlocal < dt)
				dt = dtlocal;
		}
	}
	return dt;
}

void TimeStepperExplicit::calcFluxes()
{
	for (auto row : *structuredCells)
	{
		for (auto cell : row)
		{
			cell.calcDownFlux();
			cell.calcRightFlux();
		}
	}
}

void TimeStepperExplicit::executeTimeStepGlobal()
{
	double timestep = reconstructAndMaxCFLGlobal();
	calcFluxes();

	for (auto rowit = structuredCells->begin(); rowit != structuredCells->end(); rowit++)
	{
		auto &row = *rowit;
		for (auto cellit = row.begin(); cellit != row.end(); cellit++)
		{
			auto &cell = *cellit;
			ConservativeVariables consnew = executeFirstOrderBackward(cell, timestep);
			cell.setConservative(consnew);
		}
	}
}
