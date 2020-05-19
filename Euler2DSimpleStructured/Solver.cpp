#include "Solver.h"
#include <fstream>
#include <iostream>

//Ugly, hard-coded initial conditions for Sod problem

//Constructor for solver class. Initializes time and MUSCL parameters
Solver::Solver() : Solver(1,1.0/3.0)
{
}

Solver::Solver(double eps, double kappa)
{
	time = 0;
	limit = true;
	//reconstruct_eps = eps;
	//reconstruct_kappa = kappa;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1 / 3;
}

Solver::Solver(std::string filename, double eps, double kappa) : Solver(eps, kappa)
{
	grid->readGridPro(filename);
	initSizeFromGrid();
}

//Consolidated constructor and file reader
Solver::Solver(std::string filename) : Solver()
{
	grid->readGridPro(filename);
	initSizeFromGrid();
}

void Solver::initSizeFromGrid()
{
	conservative.resize(grid->getnxiCells());
	for (int i = 0; i < conservative.size(); i++)
	{
		conservative[i].resize(grid->getnetaCells());

	}
}

//Destructor. As everything is implemented using STL containers, data deallocation is handled by STL
Solver::~Solver()
{
}

std::vector<std::unique_ptr<Boundary>>* Solver::setBoundaries(std::vector<std::unique_ptr<Boundary>>& new_boundaries)
{
	boundaries.clear(); 
	boundaries.reserve(new_boundaries.size());
	for (auto &b : new_boundaries) 
	{ 
		boundaries.push_back(std::move(b));
	}
	return &boundaries;
}

void Solver::crossPopulatePointers()
{
	stepper->setGrid(grid.get());
	stepper->setFluid(fluid.get());

	reconstruct->setGrid(grid.get());

	solution->setFluid(fluid.get());

	xi_fluxes->setFluid(fluid.get());
	xi_fluxes->setGrid(grid.get());
	xi_fluxes->setReconstruct(reconstruct.get());

	eta_fluxes->setFluid(fluid.get());
	eta_fluxes->setGrid(grid.get());
	eta_fluxes->setReconstruct(reconstruct.get());

	for (auto &b : boundaries)
		b->setGrid(grid.get());
}

//Sets conservative variables at inlet, takes user primitive variables
void Solver::setConsInlet(double p, double u, double v, double T)
{
	cons_inlet = fluid->user2cons(p, u, v, T);
}

//Sets conservative variables at initial time, takes user primitive variables
//Sets all cells in domain to the specified value
void Solver::setConsInitial(double p, double u, double v, double T)
{
	p_infty = p;
	cons_initial = fluid->user2cons(p, u, v, T);
	setInitialCondition();
}

//Takes cons-initial Solver parameter and sets all cells in domain to these conditions
void Solver::setInitialCondition()
{
	for (int i = 0; i < conservative.size(); i++)
	{
		for (int j = 0; j < conservative[0].size(); j++)
		{
			conservative[i][j] = cons_initial;
		}
	}
}

//void Solver::setCharacteristicBoundaryRightOutlet()
//{
//	for (int j = 0; j < neta_cells + 1; j++)
//	{
//		const StateVector2D &cpos = conservative[nxi_cells + 1][j];
//		const StateVector2D &c = conservative[nxi_cells][j];
//		const StateVector2D &cneg = conservative[nxi_cells-1][j];
//		double u = c[1] / c[0];
//		double sound = calcSoundSpeedcons(c);
//		double rho = c[0];
//
//		double dx = (points[nxi_cells][j][0] - points[nxi_cells-1][j][0]);
//		double drho = (c[0] - cneg[0]);
//		double dp = calcPcons(c) - calcPcons(cneg);
//		double du = (c[1] / c[0] - cneg[1] / cneg[0]);
//
//		double L1 = (u - sound)*(dp - c[0] * sound*du)/dx;
//		double L2 = u * (sound*sound*drho - dp) / dx;
//		double L3 = u * (c[2]/c[0] - cneg[2]/cneg[0]) / dx;
//		double L5 = (u + sound)*(dp + c[0] * sound*du) / dx;
//
//		double p = calcPcons(cpos);
//		p = p - dt * 0.5*(L5 + L1);
//
//		cpos[2] = cpos[2] - dt * cpos[0] * L3;
//		cpos[0] = cpos[0] - dt * (L2 + 0.5*(L5 + L1)) / (sound*sound);
//		cpos[1] = cpos[1] - dt*(L5 - L1) / (2 * sound);
//		cpos[3] = p / (gamma - 1) + 0.5*(cpos[1] * cpos[1] + cpos[2] * cpos[2]) / cpos[0];
//	}
//}
//
//void Solver::setCharacteristicBoundaryUpperOutlet()
//{
//	for (int i = 0; i < nxi_cells + 1; i++)
//	{
//		const StateVector2D &cpos = conservative[i][neta_cells + 1];
//		const StateVector2D &c = conservative[i][neta_cells];
//		const StateVector2D &cneg = conservative[i][neta_cells-1];
//		double u = c[2] / c[0];
//		double sound = calcSoundSpeedcons(c);
//		double rho = c[0];
//
//		double dx = (points[i][neta_cells][1] - points[i][neta_cells-1][1]);
//
//		double drho = c[0] - cneg[0];
//		double dp = calcPcons(c) - calcPcons(cneg);
//		double du = (c[2]/c[0] - cneg[2]/cneg[0]);
//
//		//double L1 = (u - sound)*(dp - c[0] * sound*du)/dx;
//		double L2 = u * (sound*sound*drho - dp)/dx;
//		double L3 = u * (c[1]/c[0] - cneg[1]/cneg[0]) / dx;
//		double L5 = (u + sound)*(dp + rho * sound*du)/dx;
//		//double MaMax = calcMacons(conservative[0][1]);
//		double K = 0.58*(1 - 0.5*0.5)*sound / (points[0][neta_cells][1] - points[0][0][1]);
//		double p = calcPcons(cpos);
//		double L1 = K * (p - p_infty);
//		//double L1 = 0;
//
//		
//		p = p - dt * 0.5*(L5 + L1);
//
//		cpos[1] = cpos[1] - dt * cpos[0] * L3;
//		cpos[0] = cpos[0] - dt * (L2 + 0.5*(L5 + L1)) / (sound*sound);
//		cpos[2] = cpos[2] - dt*(L5 - L1) / (2 * sound);
//		cpos[3] = p / (gamma - 1) + 0.5*(cpos[1] * cpos[1] + cpos[2] * cpos[2]) / cpos[0];
//
//		//cpos[0] = cpos[0] + dx / (sound*sound)*(L2 / (u + sound) + 0.5*(L5 / (u + sound) + L1 / (u - sound));
//		//cpos[1] =
//		double x = 0;
//		x++;
//	/*	StateVector2D prim = cons2prim(cpos);
//		prim[0] = prim[0] - dt * (L2 + +0.5*(L5 + L1) / (sound*sound));
//		prim[2] = prim[2] - (L5 - L1) / (2 * rho*sound)*rho;
//		prim[1] = prim[1] - dt * rho*L3;
//		prim[3] = p / (gamma - 1)*rho + 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0];*/
//	}
//}

double Solver::getTime()
{
	return time;
}

void Solver::allocateConservative()
{
	conservative.resize(grid->getnxiCells());
	for (auto &c : conservative)
		c.resize(grid->getnetaCells());
}

//void Solver::solve(double vel, int angle, double eps, double kappa, bool limit, double convcrit)
//{
//	std::string meshpath = "../../../mesh/";
//	std::string respath = "../../../solution/";
//
//	if (limit)
//		respath += "withlimiter/";
//	else
//		respath += "nolimiter/";
//
//	std::string gridname = "Grid" + std::to_string(angle) + "deg";
//	std::string velstr = std::to_string(vel);
//	std::string infile = meshpath + gridname + ".grd";
//	std::string outfile = respath + gridname + velstr + ".res";
//	std::string residualfile = respath + gridname + velstr + "_residualsL2.csv";
//	std::string residualinffile = respath + gridname + velstr + "_residualsLinf.csv";
//
//	Solver s(infile, eps, kappa);
//
//	if (!limit)
//		s.disableLimiter();
//
//	s.setConsInlet(1.01325e5, vel, 0, 300);
//	s.setConsInitial(1.01325e5, vel, 0, 300);
//
//	int maxsteps = 20000;
//	std::vector<StateVector2D> residuals(maxsteps + 1);
//	std::vector<StateVector2D> residualsinf(maxsteps + 1);
//
//	int i = 0;
//	for (i = 0; i <= maxsteps; i++)
//	{
//		s.executeTimeStepLocal();
//		residuals[i] = s.getResidualsL2();
//		residualsinf[i] = s.getResidualsLinfty();
//		if (residuals[i][0] < convcrit)
//			break;
//		if (i % 1000 == 0)
//		{
//			std::cout << residuals[i][0] << "\t" << residuals[i][1] << "\t" << residuals[i][2] << "\t" << residuals[i][3] << "\n";
//			//if (s.getResidual() < 1e-9)
//			//	break;
//		}
//	}
//
//
//	s.writeSolution(outfile);
//	residuals.resize(i + 1);
//	residualsinf.resize(i + 1);
//
//	writeResidual(residuals, residualfile);
//	writeResidual(residualsinf, residualinffile);
//}

void Solver::solve()
{
	if (!grid)
		throw;
	if (!stepper)
		throw;
	if (!reconstruct)
		throw;
	if (!fluid)
		throw;
	if (!solution)
		throw;
	if (!xi_fluxes)
		throw;
	if (!eta_fluxes)
		throw;

	allocateConservative();

	for (auto &b : boundaries)
		b->init();
	reconstruct->setConservative(&conservative);
	xi_fluxes->setConservative(&conservative);
	eta_fluxes->setConservative(&conservative);
	for (int i = 0; i <= maxIter; i++)
	{
		//StateVector2D calcFluxMUSCL(const StateVector2D & c_left_left, const StateVector2D & c_left, const StateVector2D & c_right, const StateVector2D & c_right_right, double nx, double ny);
		//StateVector2D calcFluxMUSCL(const StateVector2D & c_left_left, const StateVector2D & c_left, const StateVector2D & c_right, const StateVector2D & c_right_right, double Sx_left_left, double Sx_left, double Sx_right, double Sx_right_right, double nx, double ny);
		//for (auto &b : boundaries)
		//	b->apply();

		//Calculate fluxes
		xi_fluxes->calcFluxes();
		eta_fluxes->calcFluxes();

		old_conservative = conservative;

		stepper->execute(&conservative, xi_fluxes->get(), xi_fluxes->get());

		solution->calcResidualsL2(old_conservative,conservative);
		solution->calcResidualsLinfty(old_conservative,conservative);

		if (i % iter_write_interval == 0)
			solution->writeSolution("sol.sol", conservative);
	}
}

//Apply Boundary condtitions
//setBoundaryInletLeft();
//setBoundaryLowerWall();
////setBoundaryOutlet();
////setBoundaryUpperOutlet();
////setWalls();
//setCharacteristicBoundaryRightOutlet();
//setCharacteristicBoundaryUpperOutlet();