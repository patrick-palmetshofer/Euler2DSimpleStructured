#include "Solver.h"
#include <fstream>
#include <iostream>

//Ugly, hard-coded initial conditions for Sod problem
void Solver::setSodXInitial()
{
	StateVector2D v1 = user2cons(5e5, 0, 0, 300);
	StateVector2D v2 = user2cons(1e5, 0, 0, 300);
	for (int i = 0; i < (nxi_cells + 2)/2; i++)
	{
		for (int j = 0; j < neta_cells + 2; j++)
		{
			conservative[i][j] = v1;
		}
	}
	for (int i = (nxi_cells + 2) / 2; i < nxi_cells + 2; i++)
	{
		for (int j = 0; j < neta_cells + 2; j++)
		{
			conservative[i][j] = v2;
		}
	}
}
void Solver::setSodYInitial()
{
	StateVector2D v1 = user2cons(5e5, 0, 0, 300);
	StateVector2D v2 = user2cons(1e5, 0, 0, 300);
	for (int i = 0; i < nxi_cells + 2; i++)
	{
		for (int j = 0; j < (neta_cells + 2) / 2; j++)
		{
			conservative[i][j] = v1;
		}
	}
	for (int i = 0; i < nxi_cells + 2; i++)
	{
		for (int j = (neta_cells + 2) / 2; j < neta_cells + 2; j++)
		{
			conservative[i][j] = v2;
		}
	}
}

//Checks for errors in matrices, used for debugging
bool checkNaN(StateMatrix2D &m)
{
	bool ret = false;
	for (int i = 0; i < m.size(); i++)
	{
		for (int j = 0; j < m[i].size(); j++)
		{
			for (int k = 0; k < m[i][j].size(); k++)
			{
				if (!std::isfinite(m[i][j][k]))
				{
					ret = true;
					std::cout << "Calculation aborted due to NaN value at i=" << i << " j=" << j << "\n";
					//throw;
				}
			}
		}
	}
	return ret;
}
bool checkNaN(StateVector2D &m)
{
	bool ret = false;
	for (int k = 0; k < m.size(); k++)
	{
		if (!std::isfinite(m[k]))
		{
			ret = true;
			//throw;
		}
	}
	return ret;
}


//Constructor for solver class. Initializes time and MUSCL parameters
Solver::Solver() : Solver(1,1.0/3.0)
{
}

Solver::Solver(double eps, double kappa)
{
	time = 0;
	limit = true;
	reconstruct_eps = eps;
	reconstruct_kappa = kappa;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1.0 / 3.0;// 1 / 3;
}

Solver::Solver(std::string filename, double eps, double kappa) : Solver(eps, kappa)
{
	readGridGridPro(filename);
}

//Consolidated constructor and file reader
Solver::Solver(std::string filename) : Solver()
{
	readGridGridPro(filename);
}

//Destructor. As everything is implemented using STL containers, data deallocation is handled by STL
Solver::~Solver()
{
}

//Reads a GridPro grid file. Structured grid.
void Solver::readGridGridPro(std::string filename)
{
	std::ifstream stream;

	//Number of points in both computational directions
	int nxi_points;
	int neta_points;
	//std::vector<std::vector<std::array<double,2>>> points;

	try
	{
		stream.open(filename, std::ifstream::in);
		if (stream.fail())
			throw;
		stream >> nxi_points;
		stream >> neta_points;
		int zpoints;
		stream >> zpoints;

		nxi_cells = nxi_points-1;
		neta_cells = neta_points-1;

		points.resize(nxi_points);
		for (int i = 0; i < nxi_points; i++)
		{
			points[i].resize(neta_points);
		}
		double zcord;


		for (int i = 0; i < nxi_points; i++)
		{
			for (int j = 0; j < neta_points; j++)
			{
				stream >> points[i][j][0];
				stream >> points[i][j][1];
				stream >> zcord;
			}
		}
		stream.close();
	}
	catch (std::ifstream::failure e) {
		std::cerr << "Exception reading file\n";
	}


	//Allocate all matrices
	xi_fluxes.resize(nxi_cells + 1);
	nxi_xs.resize(nxi_cells + 1);
	nxi_ys.resize(nxi_cells + 1);
	Sxis.resize(nxi_cells + 1);
	
	eta_fluxes.resize(nxi_cells);
	neta_xs.resize(nxi_cells);
	neta_ys.resize(nxi_cells);
	Setas.resize(nxi_cells);

	for (int i = 0; i < xi_fluxes.size(); i++)
	{
		nxi_xs[i].resize(neta_cells);
		nxi_ys[i].resize(neta_cells);
		Sxis[i].resize(neta_cells);
		xi_fluxes[i].resize(neta_cells);
	}
	for (int i = 0; i < eta_fluxes.size(); i++)
	{
		neta_xs[i].resize(neta_cells + 1);
		neta_ys[i].resize(neta_cells + 1);
		Setas[i].resize(neta_cells + 1);
		eta_fluxes[i].resize(neta_cells + 1);
	}

	volumes.resize(nxi_cells);
	for (int i = 0; i < volumes.size(); i++)
	{
		volumes[i].resize(neta_cells);
	}
	conservative.resize(nxi_cells + 2);
	for (int i = 0; i < conservative.size(); i++)
	{
		conservative[i].resize(neta_cells + 2);
		
	}


	//Fill matrices
	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells+1; j++)
		{
			double Setax = -(points[i + 1][j][1] - points[i][j][1]);
			double Setay = (points[i + 1][j][0] - points[i][j][0]);

			double Seta = std::sqrt(Setax*Setax + Setay * Setay);

			neta_xs[i][j] = Setax / Seta;
			neta_ys[i][j] = Setay / Seta;

			Setas[i][j] = Seta;
		}
	}
	for (int i = 0; i < nxi_cells+1; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			double Sxix = points[i][j + 1][1] - points[i][j][1];
			double Sxiy = -(points[i][j + 1][0] - points[i][j][0]);

			double Sxi = std::sqrt(Sxix*Sxix + Sxiy * Sxiy);

			nxi_xs[i][j] = Sxix / Sxi;
			nxi_ys[i][j] = Sxiy / Sxi;

			Sxis[i][j] = Sxi;
		}
	}

	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			volumes[i][j] = std::abs(0.5*((points[i+1][j + 1][0] - points[i][j][0])*(points[i][j + 1][1] - points[i + 1][j][1]) - (points[i][j + 1][0] - points[i + 1][j][0])*(points[i + 1][j + 1][1] - points[i][j][1])));
		}
	}
}

//Sets conservative variables at inlet, takes user primitive variables
void Solver::setConsInlet(double p, double u, double v, double T)
{
	cons_inlet = user2cons(p, u, v, T);
}

//Sets conservative variables at initial time, takes user primitive variables
//Sets all cells in domain to the specified value
void Solver::setConsInitial(double p, double u, double v, double T)
{
	p_infty = p;
	cons_initial = user2cons(p, u, v, T);
	setInitialCondition();
}

//Takes cons-initial Solver parameter and sets all cells in domain to these conditions
void Solver::setInitialCondition()
{
	for (int i = 0; i < nxi_cells + 2; i++)
	{
		for (int j = 0; j < neta_cells + 2; j++)
		{
			conservative[i][j] = cons_initial;
		}
	}
}

//Functions to convert between conservative and primitive variables
StateVector2D Solver::prim2cons(StateVector2D & p)
{
	StateVector2D c;
	c[0] = p[0];
	c[1] = p[0]*p[1];
	c[2] = p[0]*p[2];
	c[3] = p[0]*p[3];
	return c;
}

StateVector2D Solver::cons2prim(StateVector2D & c)
{
	StateVector2D p;
	p[0] = c[0];
	p[1] = c[1] / c[0];
	p[2] = c[2] / c[0];
	p[3] = c[3] / c[0];
	return p;
}

//Functions to calculate physical properties (and charactersitic variables) from conservative/primite variables
double Solver::calcPcons(StateVector2D & c)
{
	double p = (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2])/c[0]);
	return p;
}

double Solver::calcMacons(StateVector2D & c)
{
	double Ma = std::sqrt(c[1] * c[1] + c[2] * c[2]) / (c[0] * calcSoundSpeedcons(c));
	return Ma;
}

double Solver::calcSoundSpeedcons(StateVector2D &c)
{
	double sound = std::sqrt((gamma - 1)*gamma*(c[3] - 0.5*(c[1]*c[1] + c[2] * c[2])/c[0])/c[0]);
	return sound;
}

double Solver::calcPprim(StateVector2D & prim)
{
	double p = prim[0] * (gamma - 1)*(prim[3] - 0.5*(prim[1] * prim[1] + prim[2] * prim[2]));
	return p;
}

//Limiter function in terms of r. can lead to problems if r->infty
StateVector2D Solver::limiterMinmod(StateVector2D &rs)
{
	StateVector2D phis;
	for (int i = 0; i < rs.size(); i++)
	{
		double r = rs[i];

		double phi = std::min(1.0, r);
		phi = std::max(phi, 0.0);
		phis[i] = phi;// *2 / (1 + r);

		if (!std::isfinite(r))
			phis[i] = 1.0;
	}
	return phis;
}

//Safer form of the minmod limiter
StateVector2D Solver::limiterMinmod(StateVector2D &x, StateVector2D &y)
{
	StateVector2D phis;
	for (int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		double yi = y[i];
		double phi = 0;

		if (xi*yi > 0 && std::isfinite(xi) && std::isfinite(yi))
		{
			if (std::abs(xi) < std::abs(yi))
				phi = xi/yi;
			else
				phi = 1;
		}
		phis[i] = phi;// *2 / (1 + r);
	}
	return phis;
}

//Safer form of the minmod limiter
StateVector2D Solver::limitervanAlbada(StateVector2D &x, StateVector2D &y)
{
	StateVector2D phis;
	for (int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		double yi = y[i];
		double phi = 0;

		if (xi*yi > 0 && std::isfinite(xi) && std::isfinite(yi))
		{
			phi = 2 / (xi / yi + yi / xi);
		}
		phis[i] = phi;// *2 / (1 + r);
	}
	return phis;
}

//Monotone centered limiter
StateVector2D Solver::limiterMC(StateVector2D &rs)
{
	StateVector2D phis;
	for (int i = 0; i < rs.size(); i++)
	{
		double r = rs[i];
		double m = std::min(2 * r, 2.0);
		double phi = std::min(m, 0.5*(1+r));
		phi = std::max(phi, 0.0);
		phis[i] = phi;
		if (!std::isfinite(r))
			phis[i] = 0;
	}
	return phis;
}

StateVector2D Solver::limiter(StateVector2D &r)
{
	return limiterMinmod(r);
}

StateVector2D Solver::limiter(StateVector2D &x, StateVector2D &y)
{
	return limiterMinmod(x,y);
}

//Hard-coded boundary conditions. Not pretty
void Solver::setBoundaryInletLeft()
{
	for (int j = 0; j < neta_cells+1; j++)
	{
		conservative[0][j] = cons_inlet;
	}
}

void Solver::setBoundaryLowerWall()
{
	for (int i = 0; i < nxi_cells; i++)
	{
		double unorm = conservative[i + 1][1][1] * neta_xs[i][0] + conservative[i + 1][1][2] * neta_ys[i][0];

		conservative[i + 1][0][0] = conservative[i + 1][1][0];
		conservative[i + 1][0][1] = conservative[i + 1][1][1] - 2 * unorm*neta_xs[i][0];
		conservative[i + 1][0][2] = conservative[i + 1][1][2] - 2 * unorm*neta_ys[i][0];
		conservative[i + 1][0][3] = conservative[i + 1][1][3];
	}
}

void Solver::setBoundaryUpperLowerWalls()
{
	for (int i = 0; i < nxi_cells; i++)
	{
		double unorm = conservative[i+1][1][1] * neta_xs[i][0] + conservative[i + 1][1][2] * neta_ys[i][0];
		//double upar = conservative[i][1][1] * neta_ys[i][0] + conservative[i][1][2] * neta_xs[i][0];
		//conservative[i][0][1] = (unorm * neta_xs[i][0] - upar * neta_ys[i][0]) / (neta_xs[i][0] * neta_xs[i][0] - neta_ys[i][0] * neta_ys[i][0]);
		//conservative[i][0][2] = (-unorm * neta_ys[i][0] - upar * neta_xs[i][0]) / (neta_ys[i][0] * neta_ys[i][0] - neta_xs[i][0] * neta_xs[i][0]);

		conservative[i + 1][0][0] = conservative[i + 1][1][0];
		conservative[i + 1][0][1] = conservative[i + 1][1][1] - 2 * unorm*neta_xs[i][0];
		conservative[i + 1][0][2] = conservative[i + 1][1][2] - 2 * unorm*neta_ys[i][0];
		conservative[i + 1][0][3] = conservative[i + 1][1][3];

		unorm = conservative[i + 1][neta_cells][1] * neta_xs[i][neta_cells] + conservative[i + 1][neta_cells][2] * neta_ys[i][neta_cells];

		conservative[i + 1][neta_cells + 1][0] = conservative[i + 1][neta_cells][0];
		conservative[i + 1][neta_cells + 1][1] = conservative[i + 1][neta_cells][1] - 2 * unorm*neta_xs[i][neta_cells];
		conservative[i + 1][neta_cells + 1][2] = conservative[i + 1][neta_cells][2] - 2 * unorm*neta_ys[i][neta_cells];
		conservative[i + 1][neta_cells + 1][3] = conservative[i + 1][neta_cells][3];
		//upar = conservative[i][neta_cells][1] * neta_ys[i][neta_cells + 1] + conservative[i][neta_cells][2] * neta_xs[i][neta_cells + 1];

		
		//conservative[i][neta_cells + 1][1] = (unorm * neta_xs[i][neta_cells + 1] + upar * neta_ys[i][neta_cells + 1]) / (neta_xs[i][neta_cells + 1] * neta_xs[i][neta_cells + 1] - neta_ys[i][neta_cells + 1] * neta_ys[i][neta_cells + 1]);
		//conservative[i][neta_cells + 1][2] = (-unorm * neta_ys[i][neta_cells + 1] - upar * neta_xs[i][neta_cells + 1]) / (neta_ys[i][neta_cells + 1] * neta_ys[i][neta_cells + 1] - neta_xs[i][neta_cells + 1] * neta_xs[i][neta_cells + 1]);
	}
}

void Solver::setWalls()
{
	for (int j = 0; j < neta_cells; j++)
	{
		double unorm = conservative[1][j+1][1] * nxi_xs[0][j] + conservative[1][j+1][2] * nxi_ys[0][j];
		//double upar = conservative[i][1][1] * neta_ys[i][0] + conservative[i][1][2] * neta_xs[i][0];
		//conservative[i][0][1] = (unorm * neta_xs[i][0] - upar * neta_ys[i][0]) / (neta_xs[i][0] * neta_xs[i][0] - neta_ys[i][0] * neta_ys[i][0]);
		//conservative[i][0][2] = (-unorm * neta_ys[i][0] - upar * neta_xs[i][0]) / (neta_ys[i][0] * neta_ys[i][0] - neta_xs[i][0] * neta_xs[i][0]);

		conservative[0][j + 1][0] = conservative[1][j + 1][0];
		conservative[0][j + 1][1] = conservative[1][j + 1][1] - 2 * unorm*nxi_xs[0][j];
		conservative[0][j + 1][2] = conservative[1][j + 1][2] - 2 * unorm*nxi_ys[0][j];
		conservative[0][j + 1][3] = conservative[1][j + 1][3];

		unorm = conservative[nxi_cells][j + 1][1] * nxi_xs[nxi_cells][j] + conservative[nxi_cells][j + 1][2] * nxi_ys[nxi_cells][j];

		conservative[nxi_cells + 1][j + 1][0] = conservative[nxi_cells][j + 1][0];
		conservative[nxi_cells + 1][j + 1][1] = conservative[nxi_cells][j + 1][1] - 2 * unorm*nxi_xs[nxi_cells][j];
		conservative[nxi_cells + 1][j + 1][2] = conservative[nxi_cells][j + 1][2] - 2 * unorm*nxi_ys[nxi_cells][j];
		conservative[nxi_cells + 1][j + 1][3] = conservative[nxi_cells][j + 1][3];
		//upar = conservative[i][neta_cells][1] * neta_ys[i][neta_cells + 1] + conservative[i][neta_cells][2] * neta_xs[i][neta_cells + 1];


		//conservative[i][neta_cells + 1][1] = (unorm * neta_xs[i][neta_cells + 1] + upar * neta_ys[i][neta_cells + 1]) / (neta_xs[i][neta_cells + 1] * neta_xs[i][neta_cells + 1] - neta_ys[i][neta_cells + 1] * neta_ys[i][neta_cells + 1]);
		//conservative[i][neta_cells + 1][2] = (-unorm * neta_ys[i][neta_cells + 1] - upar * neta_xs[i][neta_cells + 1]) / (neta_ys[i][neta_cells + 1] * neta_ys[i][neta_cells + 1] - neta_xs[i][neta_cells + 1] * neta_xs[i][neta_cells + 1]);
	}
	for (int i = 0; i < nxi_cells; i++)
	{
		double unorm = conservative[i + 1][1][1] * neta_xs[i][0] + conservative[i + 1][1][2] * neta_ys[i][0];
		//double upar = conservative[i][1][1] * neta_ys[i][0] + conservative[i][1][2] * neta_xs[i][0];
		//conservative[i][0][1] = (unorm * neta_xs[i][0] - upar * neta_ys[i][0]) / (neta_xs[i][0] * neta_xs[i][0] - neta_ys[i][0] * neta_ys[i][0]);
		//conservative[i][0][2] = (-unorm * neta_ys[i][0] - upar * neta_xs[i][0]) / (neta_ys[i][0] * neta_ys[i][0] - neta_xs[i][0] * neta_xs[i][0]);

		conservative[i + 1][0][0] = conservative[i + 1][1][0];
		conservative[i + 1][0][1] = conservative[i + 1][1][1] - 2 * unorm*neta_xs[i][0];
		conservative[i + 1][0][2] = conservative[i + 1][1][2] - 2 * unorm*neta_ys[i][0];
		conservative[i + 1][0][3] = conservative[i + 1][1][3];

		unorm = conservative[i + 1][neta_cells][1] * neta_xs[i][neta_cells] + conservative[i + 1][neta_cells][2] * neta_ys[i][neta_cells];


		conservative[i + 1][neta_cells + 1][0] = conservative[i + 1][neta_cells][0];
		conservative[i + 1][neta_cells + 1][1] = conservative[i + 1][neta_cells][1] - 2 * unorm*neta_xs[i][neta_cells];
		conservative[i + 1][neta_cells + 1][2] = conservative[i + 1][neta_cells][2] - 2 * unorm*neta_ys[i][neta_cells];
		conservative[i + 1][neta_cells + 1][3] = conservative[i + 1][neta_cells][3];
		//upar = conservative[i][neta_cells][1] * neta_ys[i][neta_cells + 1] + conservative[i][neta_cells][2] * neta_xs[i][neta_cells + 1];
		//if (std::abs(unorm) > 0)
		//	std::cout << "Boundary reached" << std::endl;

		//conservative[i][neta_cells + 1][1] = (unorm * neta_xs[i][neta_cells + 1] + upar * neta_ys[i][neta_cells + 1]) / (neta_xs[i][neta_cells + 1] * neta_xs[i][neta_cells + 1] - neta_ys[i][neta_cells + 1] * neta_ys[i][neta_cells + 1]);
		//conservative[i][neta_cells + 1][2] = (-unorm * neta_ys[i][neta_cells + 1] - upar * neta_xs[i][neta_cells + 1]) / (neta_ys[i][neta_cells + 1] * neta_ys[i][neta_cells + 1] - neta_xs[i][neta_cells + 1] * neta_xs[i][neta_cells + 1]);
	}
}

void Solver::setBoundaryOutlet()
{
	for (int j = 0; j < neta_cells + 2; j++)
	{
		conservative[nxi_cells+1][j] = conservative[nxi_cells][j];
	}
}

void Solver::setBoundaryUpperOutlet()
{
	for (int i = 0; i < nxi_cells + 2; i++)
	{
		conservative[i][neta_cells +1] = conservative[i][neta_cells];
	}
}

//void Solver::setCharacteristicBoundaryRightOutlet()
//{
//	for (int j = 0; j < neta_cells + 1; j++)
//	{
//		StateVector2D &cpos = conservative[nxi_cells + 1][j];
//		StateVector2D &c = conservative[nxi_cells][j];
//		StateVector2D &cneg = conservative[nxi_cells-1][j];
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
//		StateVector2D &cpos = conservative[i][neta_cells + 1];
//		StateVector2D &c = conservative[i][neta_cells];
//		StateVector2D &cneg = conservative[i][neta_cells-1];
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

//Hard coded LODI boundary conditions. Still not pretty, but better for subsonic.
//Really only necessary for Upper Outlet
void Solver::setCharacteristicBoundaryRightOutlet()
{
	for (int j = 1; j < neta_cells + 1; j++)
	{
		StateVector2D &cpos = conservative[nxi_cells + 1][j];
		StateVector2D &c = conservative[nxi_cells][j];
		StateVector2D &cneg = conservative[nxi_cells - 1][j];
		double u = c[1] / c[0];
		double sound = calcSoundSpeedcons(c);
		double rho = c[0];

		double dx = (points[nxi_cells][j][0] - points[nxi_cells - 1][j][0]);
		double drho = (c[0] - cneg[0]);
		double dp = calcPcons(c) - calcPcons(cneg);
		double du = (c[1] / c[0] - cneg[1] / cneg[0]);

		//cpos = c + (c - cneg);
		double p = calcPcons(c);

		double Ma = u / sound;
		double L1 = (u - sound)*(dp - rho * sound*du) / dx;
		if (Ma*Ma < 1)
		{
			double K = 0.58*(1 - Ma * Ma)*sound / (points[0][neta_cells][1] - points[0][0][1]);

			
			L1 = K * (p - p_infty);
		}

		double L2 = u * (sound*sound*drho - dp) / dx;
		double L3 = u * (c[2] / c[0] - cneg[2] / cneg[0]) / dx;
		double L5 = (u + sound)*(dp + rho * sound*du) / dx;

		double ppos = p + dx * 0.5*(L5 / (u + sound) + L1 / (u - sound));

		StateVector2D prim = cons2prim(c);
		StateVector2D primpos = cons2prim(cpos);
		primpos[0] = prim[0] + dx * (L2 / u + 0.5*(L5 / (u + sound) + L1 / (u - sound))) / (sound*sound);
		primpos[1] = prim[1] + dx * (L5 / (u + sound) - L1 / (u - sound)) / (2 * sound*rho);
		primpos[2] = prim[2] + dx * L3/u;
		primpos[3] = ppos / ((gamma - 1)*primpos[0]) + 0.5*(primpos[1] * primpos[1] + primpos[2] * primpos[2]);

		cpos = prim2cons(primpos);
		//StateVector2D cpos2 = c + (c - cneg);
		//for (int k = 0; k < 4; k++)
		//{
		//	if (std::abs((cpos2[k] - cpos[k])/cpos[k]) > 0.1 &&)
		//	{
		//		double palt = calcPcons(cpos2);
		//		throw;
		//	}
		//}
		if (u == 0 || u + sound == 0 || u - sound == 0) {
			cpos = c + (c - cneg);
		}
	}
}

void Solver::setCharacteristicBoundaryUpperOutlet()
{
	for (int i = 1; i < nxi_cells + 1; i++)
	{
		StateVector2D &cpos = conservative[i][neta_cells + 1];
		StateVector2D &c = conservative[i][neta_cells];
		StateVector2D &cneg = conservative[i][neta_cells - 1];

	/*	cpos = c + (c - cneg);*/
		double u = c[2] / c[0];
		double sound = calcSoundSpeedcons(c);
		double rho = c[0];

		double dx = (points[i][neta_cells][1] - points[i][neta_cells - 1][1]);

		double drho = c[0] - cneg[0];
		double dp = calcPcons(c) - calcPcons(cneg);
		double du = (c[2] / c[0] - cneg[2] / cneg[0]);

		//double L1 = (u - sound)*(dp - rho * sound*du) / dx;
		double L2 = u * (sound*sound*drho - dp) / dx;
		double L3 = u * (c[1] / c[0] - cneg[1] / cneg[0]) / dx;
		double L5 = (u + sound)*(dp + rho * sound*du) / dx;

		double p = calcPcons(c);
		double Ma = u/sound;
		double L1 = (u - sound)*(dp - rho * sound*du) / dx;
		if (Ma*Ma < 1)
		{
			double K = 0.58*(1 - Ma * Ma)*sound / (points[0][neta_cells][1] - points[0][0][1]);
			L1 = K * (p - p_infty);
		}
		//double L1 = 0;

		double ppos = p + dx * 0.5*(L5 / (u + sound) + L1 / (u - sound));

		StateVector2D prim = cons2prim(c);
		StateVector2D primpos = cons2prim(cpos);
		primpos[0] = prim[0] + dx * (L2/u +0.5*(L5/(u+sound) + L1/(u-sound))) / (sound*sound);
		primpos[2] = prim[2] + dx * (L5 / (u + sound) - L1 / (u - sound)) / (2 * sound*primpos[0]);
		primpos[1] = prim[1] + dx * L3 / u;
		primpos[3] = ppos / ((gamma - 1)*primpos[0]) + 0.5*(primpos[1] * primpos[1] + primpos[2] * primpos[2]);
		cpos = prim2cons(primpos);
		//StateVector2D cpos2 = c + (c - cneg);
		//for (int k = 0; k < 4; k++)
		//{
		//	if (std::abs((cpos2[k] - cpos[k]) / cpos[k]) > 0.1)
		//		throw;
		//}

		if (u == 0 || u + sound == 0 || u - sound == 0) {
			cpos = c + (c - cneg);
		}
	}
}

//Easy input variables to conservative
StateVector2D Solver::user2cons(double p, double u, double v, double T)
{
	double e = cp / gamma * T;
	double et = e + 0.5*(u*u + v * v);
	double rho = p / ((gamma - 1)*e);
	StateVector2D prim = { rho,u,v,et };
	StateVector2D cons = prim2cons(prim);
	return cons;
}

//Horizontal fluxes
void Solver::calcFluxesXi()
{
	StateVector2D c_left_left;
	StateVector2D c_left;
	StateVector2D c_right;
	StateVector2D c_right_right;
	StateVector2D flux;

	//First order treatment of boundaries
	for (int j = 0; j < neta_cells; j++)
	{
		c_left = conservative[0][j + 1];
		c_right = conservative[1][j + 1];

		flux = calcFlux(c_left, c_right, nxi_xs[0][j], nxi_ys[0][j]);
		xi_fluxes[0][j] = flux;

		c_left = conservative[nxi_cells][j + 1];
		c_right = conservative[nxi_cells+1][j + 1];

		flux = calcFlux(c_left, c_right, nxi_xs[nxi_cells][j], nxi_ys[nxi_cells][j]);
		xi_fluxes[nxi_cells][j] = flux;

	}

	//MUSCL treatment of inner cells
	for (int i = 1; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			c_left_left = conservative[i-1][j+1];
			c_left = conservative[i][j + 1];
			c_right = conservative[i + 1][j + 1];
			c_right_right = conservative[i + 2][j + 1];

			double Seta_left_left = Setas[i - 1][j];
			if (i != 1)
				Seta_left_left = Setas[i - 2][j];

			double Seta_right_right = Setas[i][j];
			if (i != nxi_cells-1)
				Seta_right_right = Setas[i + 1][j];

			flux = calcFluxMUSCL(c_left_left, c_left, c_right, c_right_right, Seta_left_left, Setas[i-1][j], Setas[i][j], Seta_right_right, nxi_xs[i][j], nxi_ys[i][j]);
			//flux = calcFluxMUSCL(c_left_left, c_left, c_right, c_right_right, nxi_xs[i][j], nxi_ys[i][j]);
			xi_fluxes[i][j] = flux;
		}
	}
}

//Utilities to swap coordinates in Statevectors
template<typename T>
T swap(T &data, double ind1, double ind2)
{
	T res = data;
	res[ind1] = data[ind2];
	res[ind2] = data[ind1];
	return res;
}

template<typename T>
T swap(T &data)
{
	T res = data;
	res[1] = data[2];
	res[2] = data[1];
	return res;
}

//Vertical fluxes
void Solver::calcFluxesEta()
{
	StateVector2D c_left_left;
	StateVector2D c_left;
	StateVector2D c_right;
	StateVector2D c_right_right;
	StateVector2D flux;

	for (int i = 0; i < nxi_cells; i++)
	{
		c_left = swap(conservative[i + 1][0]);
		c_right = swap(conservative[i + 1][1]);

		//c_left[2] *= -1;
		//c_right[2] *= -1;
		flux = calcFlux(c_left, c_right, neta_ys[i][0], neta_xs[i][0]);
		//flux[2] *= -1;
		eta_fluxes[i][0] = swap(flux);


		c_left = swap(conservative[i + 1][neta_cells]);
		c_right = swap(conservative[i + 1][neta_cells + 1]);

		//c_left[2] *= -1;
		//c_right[2] *= -1;
		flux = calcFlux(c_left, c_right, neta_ys[i][neta_cells], -neta_xs[i][neta_cells]);
		//flux[2] *= -1;
		eta_fluxes[i][neta_cells] = swap(flux);
	}

	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 1; j < neta_cells; j++)
		{
			c_left_left = swap(conservative[i + 1][j-1]);
			c_left = swap(conservative[i + 1][j]);
			c_right = swap(conservative[i + 1][j+1]);
			c_right_right = swap(conservative[i + 1][j+2]);

			double Sxi_left_left = Sxis[i][j - 1];
			if (j != 1)
				Sxi_left_left = Sxis[i][j - 2];

			double Sxi_right_right = Sxis[i][j];
			if (j != neta_cells-1)
				Sxi_right_right = Sxis[i][j + 1];

			//c_left_left[2] *= -1;
			//c_left[2] *= -1;
			//c_right[2] *= -1;
			//c_right_right[2] *= -1;
			flux = calcFluxMUSCL(c_left_left, c_left, c_right, c_right_right, Sxi_left_left, Sxis[i][j-1], Sxis[i][j], Sxi_right_right, neta_ys[i][j], neta_xs[i][j]);
			//flux = calcFluxMUSCL(c_left_left, c_left, c_right, c_right_right, neta_ys[i][j], neta_xs[i][j]);
			//flux[3] *= -1;
			eta_fluxes[i][j] = swap(flux);
		}
	}
}

//Higher-order flux calculation. Reconstructs states and then calls first order flux
StateVector2D Solver::calcFluxMUSCL(StateVector2D & c_left_left, StateVector2D & c_left, StateVector2D & c_right, StateVector2D & c_right_right, double nx, double ny)
{
	return calcFluxMUSCL(c_left_left, c_left, c_right, c_right_right, 1.0, 1.0, 1.0, 1.0, nx, ny);
}
StateVector2D Solver::calcFluxMUSCL(StateVector2D & c_left_left, StateVector2D & c_left, StateVector2D & c_right, StateVector2D & c_right_right, double Sx_left_left, double Sx_left, double Sx_right, double Sx_right_right, double nx, double ny)
{
	StateVector2D leftdiff = (c_left - c_left_left) * 2.0 / (Sx_left_left + Sx_left);
	StateVector2D diff = (c_right - c_left) * 2.0 / (Sx_left + Sx_right);
	StateVector2D rightdiff = (c_right_right - c_right) * 2.0 / (Sx_right_right + Sx_right);

	StateVector2D slope_left = 0.25*reconstruct_eps*((1 + reconstruct_kappa)*diff + (1 - reconstruct_kappa)*leftdiff);
	StateVector2D slope_right = 0.25*reconstruct_eps*((1 + reconstruct_kappa)*diff + (1 - reconstruct_kappa)*rightdiff);
	if (limit)
	{
		//StateVector2D r_left = (c_left - c_left_left) / (c_right - c_left);
		//StateVector2D r_right = (c_right - c_left) / (c_right_right - c_right);

		//slope_left = slope_left * limiter(r_left);
		//slope_right = slope_right * limiter(r_right);

		slope_left = slope_left * limiter(leftdiff,diff);
		slope_right = slope_right * limiter(diff,rightdiff);
	}
	StateVector2D reconstruct_left = c_left + slope_left*Sx_left;
	StateVector2D reconstruct_right = c_right - slope_right*Sx_right;

	return Solver::calcFlux(reconstruct_left, reconstruct_right, nx, ny);
}

//First order flux calculation.
StateVector2D Solver::calcFlux(StateVector2D & c_left, StateVector2D & c_right, double nx, double ny)
{
	StateVector2D flux_left = calcPhysFlux(c_left, nx, ny);
	StateVector2D flux_right = calcPhysFlux(c_right, nx, ny);
	StateVector2D flux_dissip = calcDissip(c_left, c_right, nx, ny);

	StateVector2D flux = 0.5*(flux_left + flux_right - flux_dissip);

	return flux;
}

//Calculate physical flux without numerical dissipation scheme
StateVector2D Solver::calcPhysFlux(StateVector2D & c, double nx, double ny)
{
	StateVector2D flux;
	flux[0] = c[1] * nx + c[2] * ny;
	flux[1] = c[1] * c[1] / c[0] * nx + c[1] * c[2] / c[0] *ny + (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0])*nx;
	flux[2] = c[2] * c[1] / c[0] * nx + c[2] * c[2] / c[0] *ny + (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0])*ny;
	flux[3] = (gamma*c[3] - 0.5*(gamma - 1)*(c[1] * c[1] + c[2] * c[2]) / c[0])*(c[1] * nx + c[2] * ny) / c[0];
	return flux;
}

//Call dissipation scheme. Hard-coded to Roe
StateVector2D Solver::calcDissip(StateVector2D & cons_left, StateVector2D & cons_right, double nx, double ny)
{
	StateVector2D prim_left = cons2prim(cons_left);
	StateVector2D prim_right = cons2prim(cons_right);
	StateVector2D cons_roe = calcRoeDissip(prim_left, prim_right, nx, ny);
	return cons_roe;
}

//Roe solver
StateVector2D Solver::calcRoeDissip(StateVector2D &prim_left, StateVector2D &prim_right, double nx, double ny)
{
	double sqrtrho_right = std::sqrt(prim_right[0]);
	double sqrtrho_left = std::sqrt(prim_left[0]);

	double rho = sqrtrho_right * sqrtrho_left;

	// Calculate Roe averaged variables
	double u = (sqrtrho_right*prim_right[1] + sqrtrho_left * prim_left[1]) / (sqrtrho_right + sqrtrho_left);
	double v = (sqrtrho_right*prim_right[2] + sqrtrho_left * prim_left[2]) / (sqrtrho_right + sqrtrho_left);

	double pleft = calcPprim(prim_left);
	double pright = calcPprim(prim_right);

	double Hleft = prim_left[3] + pleft / prim_left[0];
	double Hright = prim_right[3] + pright / prim_right[0];
	//double Hright = prim_right[3] + (gamma - 1)*(prim_right[3] - 0.5*(prim_right[1] * prim_right[1] + prim_right[2] * prim_right[2]));
	//double Hleft = prim_left[3] + (gamma - 1)*(prim_left[3] - 0.5*(prim_left[1] * prim_left[1] + prim_left[2] * prim_left[2]));

	double h = (sqrtrho_right*Hright + sqrtrho_left * Hleft) / (sqrtrho_right + sqrtrho_left);
	
	double unorm = u * nx + v * ny;
	double upar = -u * ny + v * nx; //Check this

	double c = std::sqrt((gamma - 1)*(h - 0.5*(u*u + v * v)));

	//Calculate transformation matrices
	std::vector<StateVector2D> roevectors;
	roevectors.resize(4);
	roevectors[0] = {
		1,
		u - c * nx,
		v - c * ny,
		h - unorm * c
	};
	roevectors[1] = {
		0,
		-ny,
		nx,
		upar
	};
	roevectors[2] = {
		1,
		u,
		v,
		0.5*(u*u + v * v)
	};
	roevectors[3] = {
		1,
		u + c * nx,
		v + c * ny,
		h + unorm * c
	};

	double drho = prim_right[0] - prim_left[0];
	double dp = pright - pleft;
	//Check this
	double dunorm = (prim_right[1] - prim_left[1])* nx + (prim_right[2] - prim_left[2])*ny;
	double dupar = -(prim_right[1] - prim_left[1])* ny + (prim_right[2] - prim_left[2])*nx;

	StateVector2D roefactors;
	roefactors[0] = (dp - rho * c*dunorm) / (2 * c*c);
	roefactors[1] = rho * dupar;
	roefactors[2] = -(dp - c * c * drho) / (c*c);
	roefactors[3] = (dp + rho * c*dunorm) / (2 * c*c);

	//Entropy fix
	double unorm_L = prim_left[1] * nx + prim_left[2] * ny;
	double unorm_R = prim_right[1] * nx + prim_right[2] * ny;

	StateVector2D eigenvals;
	eigenvals[0] = unorm - c;// std::min(unorm - c, unorm_L - c);
	eigenvals[1] = unorm;
	eigenvals[2] = unorm;
	eigenvals[3] = unorm + c; //std::max(unorm + c, unorm_R + c);

	StateVector2D dissip;
	dissip.fill(0.0);
	for (int i = 0; i < 4; i++)
	{
		dissip = dissip + roefactors[i] * std::abs(eigenvals[i])*roevectors[i];
	}
	//if (std::abs(dissip[0]) > 5)
	//	checkNaN(dissip)di;
	return dissip;
}

StateVector2D Solver::getResidualsL2()
{
	return residuals_L2;
}

StateVector2D Solver::getResidualsLinfty()
{
	return residuals_Linfty;
}

//Calculates maximum time step for a single cell (does not apply maxCFL)
double Solver::calcTimeStep(int i, int j)
{
	StateVector2D p = cons2prim(conservative[i + 1][j + 1]);
	double c = calcSoundSpeedcons(conservative[i + 1][j + 1]);// std::sqrt(gamma*(gamma - 1)*(p[3] - 0.5*(p[1] * p[1] + p[2] * p[2])));

	double Uxi_velnorm = nxi_xs[i][j] * p[1] + nxi_ys[i][j] * p[2];
	double Ueta_velnorm = neta_xs[i][j] * p[1] + neta_ys[i][j] * p[2];

	double r_xi = std::abs(Uxi_velnorm) + c;
	double r_eta = std::abs(Ueta_velnorm) + c;

//	double dtlocal = std::min((points[i + 1][j][0] - points[i][j][0]) / r_xi, (points[i][j + 1][1] - points[i][j][1]) / r_eta);
	double dtlocal = std::min(Sxis[i][j] / r_xi, Setas[i][j]  / r_eta);

	return dtlocal;
}

//Calculates maximum time step for a all cells
double Solver::calcTimeStep()
{
	double dt = 1e19;
	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			//double dtlocal = std::min(1 / r_xi, 1 / r_eta);
			double dtlocal = calcTimeStep(i, j);
			if (dtlocal < dt)
				dt = dtlocal;
		}
	}
	return dt;
}

//Residuals calculation. Either returning the maximum residual or residual in each variable.
double Solver::calcResidualL2(StateMatrix2D &o, StateMatrix2D &n)
{
	double sumres = 0;
	for (int i = 1; i < nxi_cells; i++)
	{
		for (int j = 1; j < neta_cells; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				double normk = cons_initial[k];
				if (normk == 0)
					normk = cons_initial[1];
				double temp = (n[i][j][k] - o[i][j][k])/normk;
				sumres += temp * temp;
			}
		}
	}
	if (sumres < 0)
		throw;
	sumres = std::sqrt(sumres);
	return sumres;
}

StateVector2D Solver::calcResidualsL2(StateMatrix2D &o, StateMatrix2D &n)
{
	StateVector2D res;
	for (int k = 0; k < 4; k++)
	{
		double sumres = 0;
		double normk = cons_initial[k];
		if (normk == 0)
			normk = cons_initial[1];

		for (int i = 1; i < nxi_cells; i++)
		{
			for (int j = 1; j < neta_cells; j++)
			{
				double temp = (n[i][j][k] - o[i][j][k]) / normk;
				sumres += temp * temp;
			}
		}

		if (sumres < 0)
			throw;

		sumres = std::sqrt(sumres);
		res[k] = sumres;
	}	
	return res;
}

double Solver::calcResidualLinfty(StateMatrix2D &o, StateMatrix2D &n)
{
	double res = 0;
	int imax, jmax, kmax;
	for (int i = 1; i < nxi_cells; i++)
	{
		for (int j = 1; j < neta_cells; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				double normk = cons_initial[k];
				if (normk == 0)
					normk = cons_initial[1];
				double temp = std::abs((n[i][j][k] - o[i][j][k]) / normk);
				if (temp > res)
				{
					imax = i;
					jmax = j;
					kmax = k;
					res = temp;
				}
			}
		}
	}
	return res;
}

StateVector2D Solver::calcResidualsLinfty(StateMatrix2D &o, StateMatrix2D &n)
{
	StateVector2D res;
	int imax, jmax, kmax;
	for (int k = 0; k < 4; k++)
	{
		res[k] = 0;
		double normk = cons_initial[k];
		if (normk == 0)
			normk = cons_initial[1];
		for (int i = 1; i < nxi_cells; i++)
		{
			for (int j = 1; j < neta_cells; j++)
			{
				double temp = std::abs((n[i][j][k] - o[i][j][k]) / normk);
				if (temp > res[k])
				{
					imax = i;
					jmax = j;
					kmax = k;
					res[k] = temp;
				}
			}
		}
	}
	return res;
}

void Solver::executeTimeStepGlobal()
{
	//Calculate time step with CFL
	dt = maxCFL * calcTimeStep();

	//Apply Boundary condtitions
	setBoundaryInletLeft();
	setBoundaryLowerWall();
	//setBoundaryOutlet();
	//setBoundaryUpperOutlet();
	//setWalls();
	setCharacteristicBoundaryRightOutlet();
	setCharacteristicBoundaryUpperOutlet();

	//Calculate fluxes
	calcFluxesXi();
	calcFluxesEta();

	//checkNaN(xi_fluxes);
	//checkNaN(eta_fluxes);

	//Execute time step
	StateVector2D Dxi, Deta;
	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			Dxi = xi_fluxes[i+1][j] * Sxis[i+1][j] - xi_fluxes[i][j] * Sxis[i][j];
			Deta = eta_fluxes[i][j+1] * Setas[i][j+1] - eta_fluxes[i][j] * Setas[i][j];
			conservative[i+1][j+1] = conservative[i+1][j+1] - dt / volumes[i][j] * (Dxi + Deta);
			//if (std::abs(dt / volumes[i][j]*Dxi[0]) > 0.9*conservative[i][j][0] || std::abs(dt / volumes[i][j]*Deta[0]) > 0.9*conservative[i][j][0])
			//	std::cout << "Large Flux at " << i << " , " << j << std::endl;
			//conservative[i][j] = conservative[i][j] - dt / volumes[i][j] * (xi_fluxes[i][j] * Sxis[i][j] - xi_fluxes[i - 1][j] * Sxis[i - 1][j] + eta_fluxes[i][j] * Setas[i][j] - eta_fluxes[i][j - 1] * Setas[i][j - 1]);
		}
	}

	checkNaN(conservative);
	time += dt;
}

void Solver::executeTimeStepLocal()
{
	//Apply Boundary condtitions
	setBoundaryInletLeft();
	setBoundaryLowerWall();
	//setBoundaryOutlet();
	//setBoundaryUpperOutlet();
	//setWalls();
	setCharacteristicBoundaryRightOutlet();
	setCharacteristicBoundaryUpperOutlet();

	//Calculate fluxes
	calcFluxesXi();
	calcFluxesEta();

	//Store old values for calculation of residuals
	StateMatrix2D old = conservative;

	StateVector2D Dxi, Deta;
	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells; j++)
		{
			dt = maxCFL * calcTimeStep(i, j); //local time step
			Dxi = xi_fluxes[i + 1][j] * Sxis[i + 1][j] - xi_fluxes[i][j] * Sxis[i][j];
			Deta = eta_fluxes[i][j + 1] * Setas[i][j + 1] - eta_fluxes[i][j] * Setas[i][j];
			conservative[i + 1][j + 1] = conservative[i + 1][j + 1] - dt / volumes[i][j] * (Dxi + Deta);

			//if (std::abs(dt / volumes[i][j]*Dxi[0]) > 0.9*conservative[i][j][0] || std::abs(dt / volumes[i][j]*Deta[0]) > 0.9*conservative[i][j][0])
			//	std::cout << "Large Flux at " << i << " , " << j << std::endl;
			//conservative[i][j] = conservative[i][j] - dt / volumes[i][j] * (xi_fluxes[i][j] * Sxis[i][j] - xi_fluxes[i - 1][j] * Sxis[i - 1][j] + eta_fluxes[i][j] * Setas[i][j] - eta_fluxes[i][j - 1] * Setas[i][j - 1]);
		}
	}
	if (checkNaN(conservative))
	{
		//maxCFL *= 0.5;
		std::cout << "Aborting calculation with CFL number " << maxCFL << "\n";
		//setInitialCondition();
		conservative = old;
	}

	time = 0;
	residuals_L2 = calcResidualsL2(old, conservative);
	residuals_Linfty = calcResidualsLinfty(old, conservative);
}

double Solver::getTime()
{
	return time;
}

void Solver::writeSolution(std::string filename)
{
	StateMatrix2D p = conservative;
	//Write a set of primitive variables, with Temperature in place of e_t
	for (int i = 0; i < nxi_cells+2; i++)
	{
		for (int j = 0; j < neta_cells+2; j++)
		{
			p[i][j] = cons2prim(conservative[i][j]);
			p[i][j][3] = (p[i][j][3] - 0.5*(p[i][j][2] * p[i][j][2] + p[i][j][1] * p[i][j][1]))*gamma / cp;
		}
	}

	//Write to file
	std::ofstream stream;
	try
	{
		stream.open(filename, std::ofstream::out);
		stream << nxi_cells << ",\t" << neta_cells << "\n";
		stream << "rho" << ",\t" << "u" << ",\t" << "v" << ",\t" << "T" << "\n";
		//Reserve add here
		for (int i = 0; i < nxi_cells+2; i++)
		{
			for (int j = 0; j < neta_cells+2; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					stream << p[i][j][k] << ",\t";
				}
				stream << p[i][j][3] << "\n";
			}
		}
		stream.flush();
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}