#pragma once
#include "Cell.h"
class GhostCell :
	protected Cell
{
protected:
	Face * boundary;
	Cell * neighbor;
public:
	GhostCell(double volume, Face * boundary, Config & cfg);
	~GhostCell();

	PrimitiveVariables *getPrimitive();
	ConservativeVariables *getConservative();

	double getVolume();

	void reconstructStates();

	PrimitiveVariables *getBoundaryState();

	void setNeighbor(Cell *neighbor);
	Cell *getNeighborCell();

	Fluid * getFluid();

	virtual void applyBC();
};

class SupersonicInflowCell :
	public GhostCell
{
private:
	PrimitiveVariables conststate;
	SupersonicInflowCell(PrimitiveVariables state, double volume, Face * boundary, Config & cfg) : GhostCell(volume, boundary, cfg)
	{
		conststate = state;
	}
	void applyBC() 
	{
		setPrimitive(conststate);
	}
};

class SupersonicOutflowCell :
	public GhostCell
{
	void applyBC()
	{
		setPrimitive(*neighbor->getPrimitive());
	}
};

//class SubsonicInflowCell :
//	public GhostCell
//{
//	void applyBC();
//};

//class SubsonicOutflowCell :
//	public GhostCell
//{
//	void applyBC();
//};

class SlipWallCell :
	public GhostCell
{
private:
	double Sx, Sy, S;
public:
	SlipWallCell(double S, double Sx, double Sy, double volume, Face * boundary,  Config &cfg) : GhostCell(volume, boundary, cfg)
	{
		this->S = S;
		this->Sy = Sy;
		this->Sx = Sx;
	}
	void applyBC()
	{
		PrimitiveVariables p;
		p = *neighbor->getPrimitive();
		PrimitiveVariables normal;
		normal[1] = p[1] * Sx/S + p[2] * Sy/S;
		normal[2] = p[2] * Sx/S - p[1] * Sy/S;
		normal[1] *= -1;

		p[1] = -(normal[1] * Sx/S + normal[2] * Sy/S);
		p[2] = -(normal[1] * Sx/S - normal[2] * Sy/S);

	}
};