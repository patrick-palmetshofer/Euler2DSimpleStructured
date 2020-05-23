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

	fluxes.resize(grid->getnComponentCells(0) + ydim,std::vector<StateVector2D>(grid->getnComponentCells(1) + xdim));
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

	for (int i = ydim; i < grid->getnxiCells(); i++)
	{
		for (int j = xdim; j < grid->getnetaCells(); j++)
		{
			std::pair<StateVector2D, StateVector2D> leftrightstates = reconstruct->reconstructStates(i, j, dim);
			double nx = grid->getnComponent(i, j, xdim, xdim);
			double ny = grid->getnComponent(i, j, xdim, ydim);

			std::pair<StateVector2D, StateVector2D> swapped = { swap(leftrightstates.first), swap(leftrightstates.second) };
			StateVector2D flux = calcFlux(swapped,nx,ny);

			fluxes[i][j] = swap(flux);
		}
	}
}

StateVector2D Flux::calcFlux(std::pair<StateVector2D, StateVector2D> leftrightstates, double nx, double ny)
{
	return 0.5*(fluid->calcPhysFlux(leftrightstates.first,nx,ny)+ fluid->calcPhysFlux(leftrightstates.second,nx,ny)-calcDissip(leftrightstates,nx,ny));
}
