#include "Solver.h"
#include <fstream>
#include <iostream>


bool checkNaN(StateMatrix2D &m)
{
	bool ret = false;
	for (int i = 0; i < m.size(); i++)
	{
		for (int j = 0; j < m[i].size(); j++)
		{
			for (int k = 0; k < m[i][j].size(); k++)
			{
				if (std::isnan(m[i][j][k]))
				{
					ret = true;
					throw;
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
		if (std::isnan(m[k]))
		{
			ret = true;
			throw;
		}
	}
	return ret;
}


Solver::Solver()
{
	time = 0;
	reconstruct_eps = 0;
	reconstruct_kappa = -1;
}


Solver::~Solver()
{
}


void Solver::readGridGridPro(std::string filename)
{
	std::ifstream stream;
	int nxi_points;
	int neta_points;
	std::vector<std::vector<std::array<double,2>>> points;

	try
	{
		stream.open(filename, std::ifstream::in);
		if (stream.fail())
			throw;
		stream >> nxi_cells;
		stream >> neta_cells;
		int zpoints;
		stream >> zpoints;

		nxi_points = nxi_cells;
		neta_points = neta_cells;

		nxi_cells--;
		neta_cells--;

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

	xi_fluxes.resize(nxi_cells + 1);
	eta_fluxes.resize(nxi_cells + 1);
	nxi_xs.resize(nxi_cells + 1);
	nxi_ys.resize(nxi_cells + 1);
	neta_xs.resize(nxi_cells + 1);
	neta_ys.resize(nxi_cells + 1);
	Sxis.resize(nxi_cells + 1);
	Setas.resize(nxi_cells + 1);

	for (int i = 0; i < xi_fluxes.size(); i++)
	{
		nxi_xs[i].resize(neta_cells + 1);
		nxi_ys[i].resize(neta_cells + 1);
		neta_xs[i].resize(neta_cells + 1);
		neta_ys[i].resize(neta_cells + 1);
		Sxis[i].resize(neta_cells + 1);
		Setas[i].resize(neta_cells + 1);
		xi_fluxes[i].resize(neta_cells + 1);
		eta_fluxes[i].resize(neta_cells + 1);
	}

	volumes.resize(nxi_cells + 2);
	conservative.resize(nxi_cells + 2);
	for (int i = 0; i < conservative.size(); i++)
	{
		conservative[i].resize(neta_cells + 2);
		volumes[i].resize(neta_cells + 2);
	}

	for (int i = 0; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells+1; j++)
		{
			double Setay = (points[i + 1][j][0] - points[i][j][0]);
			double Setax = -(points[i + 1][j][1] - points[i][j][1]);

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

	for (int i = 1; i < nxi_cells+1; i++)
	{
		for (int j = 1; j < neta_cells+1; j++)
		{
			volumes[i][j] = std::abs(0.5*((points[i][j][0] - points[i-1][j-1][0])*(points[i-1][j][1] - points[i][j-1][1]) - (points[i-1][j][0] - points[i][j-1][0])*(points[i][j][1] - points[i-1][j-1][1])));
		}
	}
	
	for (int i = 0; i < nxi_cells+2; i++)
	{
		volumes[i][0] = volumes[i][1];
		volumes[i][neta_cells+1] = volumes[i][neta_cells];
	}
	for (int j = 0; j < neta_cells + 2; j++)
	{
		volumes[0][j] = volumes[1][j];
		volumes[nxi_cells+1][j] = volumes[nxi_cells][j];
	}
}

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

double Solver::calcPprim(StateVector2D & prim)
{
	double p = prim[0] * (gamma - 1)*(prim[3] - 0.5*(prim[1] * prim[1] + prim[2] * prim[2]));
	return p;
}

StateVector2D Solver::limiter(StateVector2D &r)
{
	return limiterMinmod(r);
}

void Solver::setBoundaryInletLeft()
{
	for (int j = 0; j < neta_cells+1; j++)
	{
		conservative[0][j] = cons_inlet;
	}
}

void Solver::setBoundaryUpperLowerWalls()
{
	for (int i = 1; i < nxi_cells+1; i++)
	{
		double unorm = conservative[i][1][1] * neta_xs[i][0] + conservative[i][1][2] * neta_ys[i][0];
		//double upar = conservative[i][1][1] * neta_ys[i][0] + conservative[i][1][2] * neta_xs[i][0];
		//conservative[i][0][1] = (unorm * neta_xs[i][0] - upar * neta_ys[i][0]) / (neta_xs[i][0] * neta_xs[i][0] - neta_ys[i][0] * neta_ys[i][0]);
		//conservative[i][0][2] = (-unorm * neta_ys[i][0] - upar * neta_xs[i][0]) / (neta_ys[i][0] * neta_ys[i][0] - neta_xs[i][0] * neta_xs[i][0]);

		conservative[i][0] = conservative[i][1];
		conservative[i][0][1] = conservative[i][1][1] - 2 * unorm*neta_xs[i][0];
		conservative[i][0][2] = conservative[i][1][2] - 2 * unorm*neta_ys[i][0];

		unorm = conservative[i][neta_cells][1] * neta_xs[i][neta_cells + 1] + conservative[i][neta_cells][2] * neta_ys[i][neta_cells + 1];

		conservative[i][neta_cells + 1] = conservative[i][neta_cells];
		conservative[i][neta_cells + 1][1] = conservative[i][neta_cells][1] - 2 * unorm*neta_xs[i][neta_cells];
		conservative[i][neta_cells + 1][2] = conservative[i][neta_cells][2] - 2 * unorm*neta_ys[i][neta_cells];
		//upar = conservative[i][neta_cells][1] * neta_ys[i][neta_cells + 1] + conservative[i][neta_cells][2] * neta_xs[i][neta_cells + 1];

		
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

StateVector2D Solver::limiterMinmod(StateVector2D &rs)
{
	StateVector2D phis;
	for (int i = 0; i < rs.size(); i++)
	{
		double r = rs[i];
		double phi = std::min(1.0, r);
		phi = std::max(phi, 0.0);
		phis[i] = phi;
	}
	return phis;
}

void Solver::setConsInlet(double p, double u, double v, double T)
{
	double e = cp / gamma * T;
	double et = e + 0.5*(u*u + v * v);
	double rho = p / ((gamma - 1)*e);
	StateVector2D prim = { rho,u,v,et };
	cons_inlet = prim2cons(prim);
}

void Solver::calcFluxesXi()
{
	StateVector2D c_left_left;
	StateVector2D c_left;
	StateVector2D c_right;
	StateVector2D c_right_right;
	StateVector2D flux;

	for (int j = 0; j < neta_cells+1; j++)
	{
		c_left_left = conservative[0][j];
		c_left = conservative[0][j];
		c_right = conservative[1][j];
		c_right_right = conservative[2][j];

		flux = calcFlux(c_left_left, c_left, c_right, c_right_right, nxi_xs[0][j], nxi_ys[0][j]);
		xi_fluxes[0][j] = flux;
	}

	for (int j = 0; j < neta_cells+1; j++)
	{
		c_left_left = conservative[nxi_cells-1][j];
		c_left = conservative[nxi_cells][j];
		c_right = conservative[nxi_cells+1][j];
		c_right_right = conservative[nxi_cells+1][j];

		flux = calcFlux(c_left_left, c_left, c_right, c_right_right, nxi_xs[nxi_cells][j], nxi_ys[nxi_cells][j]);
		xi_fluxes[nxi_cells][j] = flux;
	}

	for (int i = 1; i < nxi_cells; i++)
	{
		for (int j = 0; j < neta_cells+1; j++)
		{
			c_left_left = conservative[i-1][j];
			c_left = conservative[i][j];
			c_right = conservative[i + 1][j];
			c_right_right = conservative[i + 2][j];

			flux = calcFlux(c_left_left, c_left, c_right, c_right_right, nxi_xs[i][j], nxi_ys[i][j]);
			xi_fluxes[i][j] = flux;
		}
	}
}

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

void Solver::calcFluxesEta()
{
	StateVector2D c_left_left;
	StateVector2D c_left;
	StateVector2D c_right;
	StateVector2D c_right_right;
	StateVector2D flux;

	for (int i = 0; i < nxi_cells+1; i++)
	{
		c_left_left = swap(conservative[i][0]);
		c_left = swap(conservative[i][0]);
		c_right = swap(conservative[i][1]);
		c_right_right = swap(conservative[i][2]);

		flux = calcFlux(c_left_left, c_left, c_right, c_right_right, neta_ys[i][0], neta_xs[i][0]);
		eta_fluxes[i][0] = swap(flux);
	}

	for (int i = 0; i < nxi_cells+1; i++)
	{
		c_left_left = swap(conservative[i][neta_cells - 1]);
		c_left = swap(conservative[i][neta_cells]);
		c_right = swap(conservative[i][neta_cells + 1]);
		c_right_right = swap(conservative[i][neta_cells + 1]);

		flux = calcFlux(c_left_left, c_left, c_right, c_right_right, neta_ys[i][neta_cells], neta_xs[i][neta_cells]);
		eta_fluxes[i][neta_cells] = swap(flux);
	}

	for (int i = 0; i < nxi_cells+1; i++)
	{
		for (int j = 1; j < neta_cells; j++)
		{
			c_left_left = swap(conservative[i][j-1]);
			c_left = swap(conservative[i][j]);
			c_right = swap(conservative[i][j+1]);
			c_right_right = swap(conservative[i][j+2]);

			flux = calcFlux(c_left_left, c_left, c_right, c_right_right, neta_ys[i][j], neta_xs[i][j]);

			eta_fluxes[i][j] = swap(flux);
		}
	}
}

StateVector2D Solver::calcFlux(StateVector2D & c_left_left, StateVector2D & c_left, StateVector2D & c_right, StateVector2D & c_right_right, double nx, double ny)
{
	StateVector2D div_left = (c_left - c_left_left);
	StateVector2D div_right = (c_right_right - c_right);

	StateVector2D r_left = (c_right - c_left) / div_left;
	StateVector2D r_right = (c_right - c_left) / div_right;

	for (int i = 0; i < 4; i++)
	{
		if (div_left[i] == 0)
			r_left[i] = 1;
		if (div_right[i] == 0)
			r_right[i] = 1;
	}
	
	StateVector2D limit_left = limiter(r_left);
	StateVector2D limit_right = limiter(r_right);

	StateVector2D rinv_left = 1 / r_left;
	StateVector2D rinv_right = 1 / r_right;

	StateVector2D limitinv_left = limiter(rinv_left);
	StateVector2D limitinv_right = limiter(rinv_right);

	for (int i = 0; i < 4; i++)
	{
		if (limit_left[i] == 0)
			limit_left[i] = 0;
		if (limit_right[i] == 0)
			limit_right[i] = 0;
	}

	StateVector2D reconstruct_left = c_left + 0.25*reconstruct_eps*(c_right - c_left)*((1 - reconstruct_kappa)*limit_left + (1 + reconstruct_kappa)*r_left * limitinv_left);
	StateVector2D reconstruct_right = c_right - 0.25*reconstruct_eps*(c_right_right - c_right)*((1 - reconstruct_kappa)*r_right * limitinv_right + (1 + reconstruct_kappa)*limit_right);

	StateVector2D flux_left = calcPhysFlux(reconstruct_left, nx, ny);
	StateVector2D flux_right = calcPhysFlux(reconstruct_right, nx, ny);
	StateVector2D flux_dissip = calcDissip(reconstruct_left, reconstruct_right, nx, ny);

	return 0.5*(flux_left + flux_right - flux_dissip);
}

StateVector2D Solver::calcPhysFlux(StateVector2D & c, double nx, double ny)
{
	StateVector2D flux;
	flux[0] = c[1] * nx + c[2] * ny;
	flux[1] = c[1] * c[1] / c[0] * nx + c[1] * c[2] / c[0] *ny + (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] + c[2] * c[2]) / c[0])*nx;
	flux[2] = c[2] * c[1] / c[0] * nx + c[2] * c[2] / c[0] *ny + (gamma - 1)*(c[3] - 0.5*(c[1] * c[1] +c[0] + c[2]) / c[0])*ny;
	flux[3] = (gamma*c[3] - 0.5*(gamma - 1)*(c[1] * c[1] + c[2] * c[2]) / c[0])*(c[1] * nx + c[2] * ny) / c[0];
	return flux;
}

StateVector2D Solver::calcDissip(StateVector2D & cons_left, StateVector2D & cons_right, double nx, double ny)
{
	StateVector2D prim_left = cons2prim(cons_left);
	StateVector2D prim_right = cons2prim(cons_right);
	StateVector2D cons_roe = calcRoeDissip(prim_left, prim_right, nx, ny);
	return cons_roe;
}

StateVector2D Solver::calcRoeDissip(StateVector2D &prim_left, StateVector2D &prim_right, double nx, double ny)
{
	double sqrtrho_right = std::sqrt(prim_right[0]);
	double sqrtrho_left = std::sqrt(prim_left[0]);

	double rho = sqrtrho_right * sqrtrho_left;

	// Calculate Roe averaged variables
	double u = (sqrtrho_right*prim_right[1] + sqrtrho_left * prim_left[1]) / (sqrtrho_right + sqrtrho_left);
	double v = (sqrtrho_right*prim_right[2] + sqrtrho_left * prim_left[2]) / (sqrtrho_right + sqrtrho_left);

	double Hright = prim_right[3] + (gamma - 1)*(prim_right[3] - 0.5*(prim_right[1] * prim_right[1] + prim_right[2] * prim_right[2]));
	double Hleft = prim_left[3] + (gamma - 1)*(prim_left[3] - 0.5*(prim_left[1] * prim_left[1] + prim_left[2] * prim_left[2]));

	double h = (sqrtrho_right*Hright + sqrtrho_left * Hleft) / (sqrtrho_right + sqrtrho_left);

	double pleft = calcPprim(prim_left);
	double pright = calcPprim(prim_right);
	
	double unorm = u * nx + v * ny;
	double upar = u * ny + v * nx;

	double c = std::sqrt((gamma - 1)*(h - 0.5*(u*u + v * v)));

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
	double dupar = (prim_right[1] - prim_left[1])* ny + (prim_right[2] - prim_left[2])*nx;
	double dunorm = (prim_right[1] - prim_left[1])* nx + (prim_right[2] - prim_left[2])*ny;

	StateVector2D roefactors;
	roefactors[0] = (dp - rho * c*dunorm) / (2 * c*c);
	roefactors[1] = rho * dupar;
	roefactors[2] = -(dp - c * drho) / (c*c);
	roefactors[3] = (dp + rho * c*dunorm) / (2 * c*c);

	StateVector2D eigenvals;
	eigenvals[0] = unorm - c;
	eigenvals[1] = unorm;
	eigenvals[2] = unorm;
	eigenvals[3] = unorm + c;

	StateVector2D dissip;
	dissip.fill(0.0);
	for (int i = 0; i < 4; i++)
	{
		dissip = dissip + roefactors[i] * std::abs(eigenvals[i])*roevectors[i];
	}
	checkNaN(dissip);
	return dissip;
}

void Solver::setInitialCondition()
{
	for (int i = 0; i < nxi_cells + 2; i++)
	{
		for (int j = 0; j < neta_cells + 2; j++)
		{
			conservative[i][j] = cons_inlet;
		}
	}
}

double Solver::getResidual()
{
	return residual;
}


double Solver::calcTimeStep()
{
	double dt = 1e19;
	for (int i = 1; i < nxi_cells - 1; i++)
	{
		for (int j = 1; j < neta_cells - 1; j++)
		{
			StateVector2D p = cons2prim(conservative[i][j]);
			double c = std::sqrt(gamma*(gamma - 1)*(p[3] - 0.5*(p[1] * p[1] + p[2] * p[2])));

			double Uxi_velnorm = nxi_xs[i][j] * p[1] + nxi_ys[i][j] * p[2];
			double Ueta_velnorm = neta_ys[i][j] * p[1] + neta_xs[i][j] * p[2];

			double r_xi = std::abs(Uxi_velnorm) + c;
			double r_eta = std::abs(Ueta_velnorm) + c;

			double dtlocal = std::min(1 / r_xi, 1 / r_eta);

			if (dtlocal < dt)
				dt = dtlocal;
		}
	}
	return dt;
}

double Solver::calcResidualL2(StateMatrix2D &o, StateMatrix2D &n)
{
	double sumres = 0;
	for (int i = 1; i < nxi_cells - 1; i++)
	{
		for (int j = 1; j < neta_cells - 1; j++)
		{
			for (int k = 0; k < 1; k++)
			{
				double temp = n[i][j][k] - o[i][j][k];
				sumres += temp * temp;
			}
		}
	}
	sumres = std::sqrt(sumres);
	return sumres;
}

void Solver::executeTimeStepGlobal()
{
	setBoundaryInletLeft();
	setBoundaryUpperLowerWalls();
	setBoundaryOutlet();

	calcFluxesXi();
	calcFluxesEta();

	checkNaN(xi_fluxes);
	checkNaN(eta_fluxes);

	dt = maxCFL*calcTimeStep();

	StateMatrix2D old = conservative;

	StateVector2D Dxi, Deta;

	for (int i = 1; i < nxi_cells+1; i++)
	{
		for (int j = 1; j < neta_cells+1; j++)
		{
			Dxi = xi_fluxes[i][j] * Sxis[i][j] - xi_fluxes[i - 1][j] * Sxis[i - 1][j];
			Deta = eta_fluxes[i][j] * Setas[i][j] - eta_fluxes[i][j - 1] * Setas[i][j - 1];
			conservative[i][j] = conservative[i][j] - dt / volumes[i][j] * (Dxi + Deta);
			if (std::abs(dt / volumes[i][j]*Dxi[0]) > 0.9*conservative[i][j][0] || std::abs(dt / volumes[i][j]*Deta[0]) > 0.9*conservative[i][j][0])
				std::cout << "Large Flux at " << i << " , " << "j" << std::endl;
			//conservative[i][j] = conservative[i][j] - dt / volumes[i][j] * (xi_fluxes[i][j] * Sxis[i][j] - xi_fluxes[i - 1][j] * Sxis[i - 1][j] + eta_fluxes[i][j] * Setas[i][j] - eta_fluxes[i][j - 1] * Setas[i][j - 1]);
		}
	}
	checkNaN(conservative);
	time += dt;
	residual = calcResidualL2(old, conservative);
}

void Solver::writeSolution(std::string filename)
{
	StateMatrix2D p = conservative;
	for (int i = 1; i < nxi_cells+1; i++)
	{
		for (int j = 1; j < neta_cells+1; j++)
		{
			p[i][j] = cons2prim(conservative[i][j]);
			p[i][j][3] = (p[i][j][3] - 0.5*(p[i][j][2] * p[i][j][2] + p[i][j][1] * p[i][j][1]))*gamma / cp;
		}
	}

	std::ofstream stream;
	try
	{
		stream.open(filename, std::ofstream::out);
		stream << nxi_cells << ",\t" << neta_cells << "\n";
		stream << "rho" << ",\t" << "u" << ",\t" << "v" << ",\t" << "T" << "\n";
		//Reserve add here
		for (int i = 1; i < nxi_cells+1; i++)
		{
			for (int j = 1; j < neta_cells+1; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					stream << p[i][j][k] << ",\t";
				}
				stream << p[i][j][3] << "\n";
			}
		}
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}