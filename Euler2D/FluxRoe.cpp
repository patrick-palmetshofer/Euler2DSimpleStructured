#include "FluxRoe.h"
#include <cmath>
#include <valarray>

RoeAverages::RoeAverages()
{
}

RoeAverages::RoeAverages(Cell * cell)
{
	prim = cell->getPrimitive();
}


FluxRoe::FluxRoe()
{
}


FluxRoe::~FluxRoe()
{
}



//Calculates the exact solution of the linear equation
// u_t + A_lr u_x = 0
// Where A_lr(u,u) = A(u)
double FluxRoe::calcFlux(Cell * leftCell, Cell * rightCell)
{
	PrimitiveVariables prim_right = rightCell->getPrimitive();
	PrimitiveVariables prim_left = leftCell->getPrimitive();

	double gamma =fluid->getGamma();

	double sqrtrho_right = std::sqrt(prim_right.rho);
	double sqrtrho_left = std::sqrt(prim_left.rho);


	// Calculate Roe averaged variables
	double u = (sqrtrho_right*prim_right.u + sqrtrho_left*prim_left.u)/(sqrtrho_right + sqrtrho_left);
	double v = (sqrtrho_right*prim_right.v + sqrtrho_left * prim_left.v) / (sqrtrho_right + sqrtrho_left);

	double Hright = rightCell->getFluid()->getEnthalpy(); //TODO 
	double Hleft = rightCell->getFluid()->getEnthalpy(); // TODO

	double h = (sqrtrho_right*Hright + sqrtrho_left * Hleft) / (sqrtrho_right + sqrtrho_left);

	double c = fluid->getSoundSpeed();
	//double vel;

	//Calculate Eigenvalues
	Euler2DVector eigenvals;
	eigenvals[0] = u - c;
	eigenvals[1] = u;
	eigenvals[2] = u;
	eigenvals[3] = u + c;

	//Calculate Eigenvectors
	std::valarray<std::valarray<double>> roevectors;
	roevectors[0] = {
		1,
		u - c,
		v,
		h - u * c
	};
	roevectors[1] = {
		1,
		u,
		v,
		0.5*(u^2 + v ^ 2)
	};
	roevectors[2] = {
		0,
		0,
		1,
		v
	};
	roevectors[3] = {
		1,
		u + c,
		v,
		h + u * c
	};

	
	double Drho = prim_right.rho / prim_left.rho;
	double Dmu	= prim_right.rho*prim_right.u - prim_left.rho*prim_left.u;
	double Dmv	= prim_right.rho*prim_right.v - prim_left.rho*prim_left.v;
	double Deavg = prim_right.e - prim_left.e - (Dmv-v*Drho)*v;

	//Calculate Roe Factors
	Euler2DVector roefactors;
	roefactors[1] = -(gamma - 1) / c ^ 2 * (Drho*(u^2-h)-u*Dmu+Deavg);
	roefactors[0] = -1 /(2*c) * (Dmu-Drho*(u+c)-0.5*roefactors[1]);
	roefactors[3] = Drho - roefactors[0] - roefactors[1];
	roefactors[2] = Dmv - v * Drho;

	std::valarray<double> flux = 0.5*(rightCell->getFlux(), leftCell->getFlux());
	for (int i = 0; i < 4; i++)
	{
		flux += roefactors[i] * std::abs(eigenvals[i])*roevectors[i];
	}
	return &flux;
}