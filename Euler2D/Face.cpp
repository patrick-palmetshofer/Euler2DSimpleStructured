#include "Face.h"


void Face::calcNormalStates()
{
	leftstate_normal = *leftstate;
	rightstate_normal = *rightstate;

	leftstate_normal[1] = (*leftstate)[1] * SxS + (*leftstate)[2] * SyS;
	leftstate_normal[2] = (*leftstate)[2] * SxS - (*leftstate)[1] * SyS;

	rightstate_normal[1] = (*rightstate)[1] * SxS + (*rightstate)[2] * SyS;
	rightstate_normal[2] = (*rightstate)[2] * SxS - (*rightstate)[1] * SyS;
}

void Face::calcPhysFlux()
{
	PrimitiveVariables p;
	p = numflux_normal;
	numflux = numflux_normal;
	numflux[1] = -(numflux_normal[1] * SxS + numflux_normal[2] * SyS);
	numflux[2] = -(numflux_normal[1] * SxS - numflux_normal[2] * SyS);
}

Face::Face(double S, double Sx, double Sy, Config &cfg)
{
	this->cfg = cfg;
	this->S = S;
	this->Sy = Sy;
	this->Sx = Sx;

	SxS = Sx / S;
	SyS = Sy / S;

	fluid = cfg.getFluid();
}

void Face::setLeftRight(PrimitiveVariables *leftstate, PrimitiveVariables *rightstate)
{
	this->leftstate = leftstate;
	this->rightstate = rightstate;
}

Face::~Face()
{
}

//void Face::calcJacobians()
//{
//	ConservativeVariables rightcons;
//	rightcons = rightstate;
//	calcJacobian(right->getVolume(), rightcons, rightJacobian);
//	ConservativeVariables leftcons;
//	leftcons = leftstate;
//	calcJacobian(left->getVolume(), leftcons, leftJacobian);
//}



void Face::calcFlux()
{
	calcNormalStates();
	ConservativeVariables c = flux->calcFlux(leftstate_normal, rightstate_normal);;
	numflux_normal = c;
	calcPhysFlux();
}

ConservativeVariables Face::getFlux()
{
	return numflux;
}

ConservativeVariables Face::getIntFlux()
{
	return numflux*S;
}

double Face::getVelNorm()
{
	double c = std::max(fluid->getSoundSpeed(leftstate_normal), fluid->getSoundSpeed(rightstate_normal));
	return std::abs(numflux_normal[1])+c;
}
