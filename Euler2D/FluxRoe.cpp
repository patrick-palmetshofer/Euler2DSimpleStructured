#include "FluxRoe.h"
#include <cmath>
#include <valarray>

RoeAverages::RoeAverages()
{
}

//Calculate phsysical flux. Note: Coordinate system must have been rotated perpendicular to face
ConservativeVariables FluxRoe::calcPhysicalFlux(PrimitiveVariables &prim)
{
	ConservativeVariables flux;
	flux[0] = prim.rho*prim.u;
	flux[1] = prim.rho*prim.u*prim.u + fluid->getPressure(prim);
	flux[2] = prim.rho*prim.u*prim.v;
	flux[3] = prim.rho*(fluid->getEnthalpy(prim)+ prim.u*prim.u + prim.v*prim.v);
	return flux;
}

FluxRoe::FluxRoe()
{
}

FluxRoe::FluxRoe(Fluid * newfluid)
{
	fluid = newfluid;
}


FluxRoe::~FluxRoe()
{
}



//Calculates the exact solution of the linear equation
// u_t + A_lr u_x = 0
// Where A_lr(u,u) = A(u)
ConservativeVariables FluxRoe::calcFlux(PrimitiveVariables &prim_left, PrimitiveVariables &prim_right)
{

	double gamma = fluid->getGamma(prim_left);

	double sqrtrho_right = std::sqrt(prim_right.rho);
	double sqrtrho_left = std::sqrt(prim_left.rho);


	// Calculate Roe averaged variables
	double u = (sqrtrho_right*prim_right.u + sqrtrho_left*prim_left.u)/(sqrtrho_right + sqrtrho_left);
	double v = (sqrtrho_right*prim_right.v + sqrtrho_left * prim_left.v) / (sqrtrho_right + sqrtrho_left);

	double Hright = fluid->getEnthalpy(prim_right); 
	double Hleft = fluid->getEnthalpy(prim_left);

	double h = (sqrtrho_right*Hright + sqrtrho_left * Hleft) / (sqrtrho_right + sqrtrho_left);

	//UGLYYY!!!
	PrimitiveVariables prim_avg;
	prim_avg.rho = 0.5 * (prim_right.rho + prim_left.rho);
	prim_avg.u = u;
	prim_avg.u = v;
	prim_avg.e = (sqrtrho_right*prim_right.e + sqrtrho_left * prim_left.e) / (sqrtrho_right + sqrtrho_left);

	double c = fluid->getSoundSpeed(prim_avg);
	//double vel;

	//Calculate Eigenvalues
	PrimitiveVariables eigenvals;
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
		0.5*(u*u + v*v)
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
	PrimitiveVariables roefactors;
	roefactors[1] = -(gamma - 1) / (c*c) * (Drho*(u*u-h)-u*Dmu+Deavg);
	roefactors[0] = -1 /(2*c) * (Dmu-Drho*(u+c)-0.5*roefactors[1]);
	roefactors[3] = Drho - roefactors[0] - roefactors[1];
	roefactors[2] = Dmv - v * Drho;

	ConservativeVariables flux = 0.5*(calcPhysicalFlux(prim_right) + calcPhysicalFlux(prim_left));
	for (int i = 0; i < 4; i++)
	{
		flux += roefactors[i] * std::abs(eigenvals[i])*roevectors[i];
	}
	return flux;
}