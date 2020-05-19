#include "Boundary.h"



Boundary::Boundary()
{
}


Boundary::~Boundary()
{
}

void Boundary::setBottom()
{
	dim = 1;
	dir = 1;
}

void Boundary::setTop()
{
	dim = 1;
	dir = -1;
}

void Boundary::setRight()
{
	dim = 0;
	dir = -1;
}

void Boundary::setLeft()
{
	dim = 0;
	dir = 1;
}

void Boundary::init()
{
	//int nmax = 0, nconst = 0;
	//if (dim)
	//{
	//}
	//if (dir < 0)
	//{

	//}
	//cell_inds.reserve(nmax);
	//for (int i = 0; i < nxi; ++i)
	//{
	//	cell_inds[i][0] = i;
	//	cell_inds[i][1] = nconst;
	//}
}

StateVector2D SupersonicInlet::getGhostState(std::array<int, 2> ind)
{
	return in_state;
}

StateVector2D SlipWall::getGhostState(std::array<int, 2> ind)
{
	int i = ind[0];
	int j = ind[1];
	int ineg = i, jneg = j;
	double nx = 0, ny = 0;
	switch (dim)
	{
	case 0:
		ineg += dir;
		nx = grid->getnXiXs(i, j);
		ny = grid->getnXiYs(i, j);
		break;
	case 1:
		jneg += dir;
		nx = grid->getnEtaXs(i, j);
		ny = grid->getnEtaYs(i, j);
		break;
	}

	double unorm = conservative->at(i).at(j)[1] * nx + conservative->at(i).at(j)[2]* ny;
	ghost_state[0] = conservative->at(i).at(j)[0];
	ghost_state[1] = conservative->at(i).at(j)[1] - 2 * unorm*nx;
	ghost_state[2] = conservative->at(i).at(j)[2] - 2 * unorm*ny;
	ghost_state[3] = conservative->at(i).at(j)[3];

	return ghost_state;
}

StateVector2D SupersonicOutlet::getGhostState(std::array<int, 2> ind)
{
	int i = ind[0];
	int j = ind[1];
	int ineg = i, jneg = j;
	double nx = 0, ny = 0;
	switch (dim)
	{
	case 0:
		ineg += dir;
		nx = grid->getnXiXs(i, j);
		ny = grid->getnXiYs(i, j);
		break;
	case 1:
		jneg += dir;
		nx = grid->getnEtaXs(i, j);
		ny = grid->getnEtaYs(i, j);
		break;
	}
	ghost_state = conservative->at(i).at(j);
	
	return ghost_state;
}

StateVector2D BoundaryLODI::getGhostState(std::array<int, 2> ind)
{
	int i = ind[0];
	int j = ind[1];
	int ineg = i, jneg = j;
	double nx = 0, ny = 0;
	switch (dim)
	{
	case 0:
		ineg += dir;
		nx = grid->getnXiXs(i, j);
		ny = grid->getnXiYs(i, j);
		break;
	case 1:
		jneg += dir;
		nx = grid->getnEtaXs(i, j);
		ny = grid->getnEtaYs(i, j);
		break;
	}

	StateVector2D &c = conservative->at(i).at(j);
	StateVector2D &cneg = conservative->at(ineg).at(jneg);

	/*	cpos = c + (c - cneg);*/
	double u = c[2] / c[0];
	double sound = fluid->calcSoundSpeedcons(c);
	double rho = c[0];

	double dx = (grid->getPoint(i,j)[dim] - grid->getPoint(ineg, jneg)[dim]);

	double drho = c[0] - cneg[0];
	double dp = fluid->calcPcons(c) - fluid->calcPcons(cneg);
	double du = (c[1] / c[0] - cneg[1] / cneg[0]);

	calcWaveStrengths(u,sound,rho,dx,drho,dp,du,c,cneg);

	double p = fluid->calcPcons(c);

	double ppos = p + dx * 0.5*(L5 / (u + sound) + L1 / (u - sound));
	StateVector2D prim = fluid->cons2prim(c);
	StateVector2D primpos;
	primpos[0] = prim[0] + dx * (L2 / u + 0.5*(L5 / (u + sound) + L1 / (u - sound))) / (sound*sound);
	primpos[1] = prim[1] + dx * (L5 / (u + sound) - L1 / (u - sound)) / (2 * sound*rho);
	primpos[2] = prim[2] + dx * L3 / u;
	primpos[3] = ppos / ((fluid->getGamma() - 1)*primpos[0]) + 0.5*(primpos[1] * primpos[1] + primpos[2] * primpos[2]);

	ghost_state = fluid->prim2cons(primpos);

	if (u == 0 || u + sound == 0 || u - sound == 0) {
		ghost_state = c + (c - cneg);
	}

	return ghost_state;
}

void Outlet::calcWaveStrengths(double u, double sound, double rho, double dx, double drho, double dp, double du, StateVector2D& c, StateVector2D& cneg)
{
	L1 = (u - sound)*(dp - rho * sound*du) / dx;
	L2 = u * (sound*sound*drho - dp) / dx;
	L3 = u * (c[1] / c[0] - cneg[1] / cneg[0]) / dx;
	L5 = (u + sound)*(dp + rho * sound*du) / dx;

	double p = fluid->calcPcons(c);
	double Ma = u / sound;
	if (Ma*Ma < 1)
	{
		double K = 0.58*(1 - Ma * Ma)*sound / grid->getMaxDistance(dim);
		L1 = K * (p - p_infty);
	}
}
