#include "FluxRoe.h"



FluxRoe::FluxRoe(int new_dim) : Flux(new_dim)
{
}


FluxRoe::~FluxRoe()
{
}

//Roe solver
StateVector2D FluxRoe::calcDissip(std::pair<StateVector2D, StateVector2D> leftrightstates, double nx, double ny)
{
	StateVector2D prim_left = fluid->cons2prim(leftrightstates.first);
	StateVector2D prim_right = fluid->cons2prim(leftrightstates.second);

	double sqrtrho_right = std::sqrt(prim_right[0]);
	double sqrtrho_left = std::sqrt(prim_left[0]);

	double rho = sqrtrho_right * sqrtrho_left;

	// Calculate Roe averaged variables
	double u = (sqrtrho_right*prim_right[1] + sqrtrho_left * prim_left[1]) / (sqrtrho_right + sqrtrho_left);
	double v = (sqrtrho_right*prim_right[2] + sqrtrho_left * prim_left[2]) / (sqrtrho_right + sqrtrho_left);

	double pleft = fluid->calcPprim(prim_left);
	double pright = fluid->calcPprim(prim_right);

	double Hleft = prim_left[3] + pleft / prim_left[0];
	double Hright = prim_right[3] + pright / prim_right[0];
	//double Hright = prim_right[3] + (gamma - 1)*(prim_right[3] - 0.5*(prim_right[1] * prim_right[1] + prim_right[2] * prim_right[2]));
	//double Hleft = prim_left[3] + (gamma - 1)*(prim_left[3] - 0.5*(prim_left[1] * prim_left[1] + prim_left[2] * prim_left[2]));

	double h = (sqrtrho_right*Hright + sqrtrho_left * Hleft) / (sqrtrho_right + sqrtrho_left);

	double unorm = u * nx + v * ny;
	double upar = -u * ny + v * nx; //Check this

	double c = std::sqrt((fluid->getGamma() - 1)*(h - 0.5*(u*u + v * v)));

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
