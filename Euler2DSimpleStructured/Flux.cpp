#include "Flux.h"

using namespace Euler;

Flux::Flux(int new_dim) : dim(new_dim)
{
}


Flux::~Flux()
{
}

void Flux::setConservative(StateMatrix2D * cons)
{
	conservative = cons;
	int xdim = 0;
	int ydim = 1;
	if (dim)
	{
		xdim = 1;
		ydim = 0;
	}

	fluxes.resize(grid->getnComponentCells(xdim) + ydim);
	for (int i = 0; i < fluxes.size(); i++)
	{
		fluxes[i].resize(grid->getnComponentCells(ydim) + xdim);
	}
}

void Flux::calcFluxes()
{
	int xdim = 0;
	int ydim = 1;
	if (dim)
	{
		xdim = 1;
		ydim = 0;
	}
	for (int i = 0; i < grid->getnxiCells()+xdim; i++)
	{
		for (int j = 0; j < grid->getnetaCells()+ydim; j++)
		{
			std::pair<StateVector2D, StateVector2D> leftrightstates = reconstruct->reconstructStates(i, j, dim);
			double nx = grid->getnComponent(i, j, xdim, xdim);
			double ny = grid->getnComponent(i, j, xdim, ydim);

			std::pair<StateVector2D, StateVector2D> swapped = { swap(leftrightstates.first), swap(leftrightstates.second) };
			StateVector2D flux = calcFlux(swapped,nx,ny);

			fluxes.at(i).at(j) = swap(flux);
		}
	}
}

StateVector2D Flux::calcFlux(std::pair<StateVector2D, StateVector2D> leftrightstates, double nx, double ny)
{
	return 0.5*(fluid->calcPhysFlux(leftrightstates.first,nx,ny)+ fluid->calcPhysFlux(leftrightstates.second,nx,ny)-calcDissip(leftrightstates,nx,ny));
}
