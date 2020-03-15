#pragma once

class PrimVars
{
	double rho;
	double u;
	double v;
	double et;
};

class ConsVars
{
	double rho;
	double rhou;
	double rhov;
	double rhoet;
};

class Cell
{
protected:
	PrimVars prim;
	ConsVars cons;
public:
	Cell();
	~Cell();

	void PrimCons();
	void ConsPrim();
};

