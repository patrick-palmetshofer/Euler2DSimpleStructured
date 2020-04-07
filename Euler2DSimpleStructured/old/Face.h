#pragma once
#include "Config.h"
#include "Flux.h"
#include "Fluid.h"
#include <vector>
class Face
{
private:
	Config &cfg;
	Flux *flux;
	Fluid *fluid;

	PrimitiveVariables *leftstate;
	PrimitiveVariables *rightstate;

	PrimitiveVariables leftstate_normal;
	PrimitiveVariables rightstate_normal;

	ConservativeVariables numflux;
	PrimitiveVariables numflux_normal;

	double S; //Area

	double Sx; //x projected Area
	double Sy; //y projected Area

	double SxS;
	double SyS;

	void calcNormalStates();

	void calcPhysFlux();

public:
	Face(double S, double Sx, double Sy, Config &cfg);
	void setLeftRight(PrimitiveVariables * leftstate, PrimitiveVariables * rightstate);
	~Face();

	void calcFlux();
	ConservativeVariables getFlux();

	ConservativeVariables getIntFlux();

	double getVelNorm();

};

