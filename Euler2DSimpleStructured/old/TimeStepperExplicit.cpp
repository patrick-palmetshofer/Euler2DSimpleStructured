#include "TimeStepperExplicit.h"

TimeStepperExplicit::TimeStepperExplicit(std::vector<std::shared_ptr<Cell>>& cells, std::vector<std::shared_ptr<Face>>& faces)
{
	this->cells = cells;
	this->faces = faces;
	time = 0;
}

TimeStepperExplicit::TimeStepperExplicit(Grid & grid)
{
	cells = grid.getCells();
	faces = grid.getFaces();
	time = 0;
}

TimeStepperExplicit::~TimeStepperExplicit()
{
}

ConservativeVariables TimeStepperExplicit::spaceDisc(std::shared_ptr<Cell> &cell)
{
	ConservativeVariables spacedisc = cell->getSumFaceFluxes();
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

ConservativeVariables TimeStepperExplicit::executeFirstOrderBackward(std::shared_ptr<Cell> &cell, double dt)
{
	ConservativeVariables consnew = *(cell->getConservative()) - spaceDisc(cell)*dt/cell->getVolume();
	return consnew;
}

//TODO:: Move all from loops to recursion

//Returns global time step and calculates fluxes in cells
void TimeStepperExplicit::reconstruct()
{
	
	for (auto cell : cells)
	{
		cell->reconstructStates();
	}
	
}
double TimeStepperExplicit::maxTimeStep()
{
	double dt = 1e19;
	for (auto face : faces)
	{
		double dtlocal = maxCFL / face->getVelNorm();
		if (dtlocal < dt)
			dt = dtlocal;
	}
	return dt;
}

void TimeStepperExplicit::executeTimeStepGlobal()
{
	reconstruct();
	double timestep = maxTimeStep();
	
	for (auto cell : boundaryCells)
	{
		cell->applyBC();
	}
	for (auto face : faces)
	{
		calcFluxes();
	}
	for (auto cell : cells)
	{
		ConservativeVariables consnew = executeFirstOrderBackward(cell, timestep);
		cell->setConservative(consnew);
	}
	time += timestep;
}

double TimeStepperExplicit::getTime()
{
	return time;
}