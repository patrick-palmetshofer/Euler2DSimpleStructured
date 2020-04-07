#pragma once
#include "GlobalTypes.h"

class Solver
{
private:
	int nxi_cells;
	int neta_cells;

	const double maxCFL = 0.01;

	double dt;
	double time;

	ValueMatrix2D volumes;

	StateVector2D cons_inlet;
	//StateMatrix2D primitive;

	double reconstruct_eps;
	double reconstruct_kappa;
	double residual;

	ValueMatrix2D Sxis;
	ValueMatrix2D Setas;

	ValueMatrix2D nxi_xs;
	ValueMatrix2D nxi_ys;

	ValueMatrix2D neta_xs;
	ValueMatrix2D neta_ys;

	//StateMatrix2D xi_face_prim;
	//StateMatrix2D eta_face_prim;

	StateMatrix2D xi_fluxes;
	StateMatrix2D eta_fluxes;
	StateMatrix2D conservative;

	double gamma = 1.4;
	double cp = 1005;

public:
	Solver();
	Solver(std::string filename);
	~Solver();

	void readGridGridPro(std::string filename);

	StateVector2D prim2cons(StateVector2D &p);
	StateVector2D cons2prim(StateVector2D &c);
	double calcPprim(StateVector2D & prim);

	StateVector2D limiter(StateVector2D &r);

	void setBoundaryInletLeft();

	void setBoundaryUpperLowerWalls();

	void setBoundaryOutlet();

	StateVector2D limiterMinmod(StateVector2D &r);

	StateVector2D calcFlux();

	StateVector2D spaceDisc(int i, int j);

	StateVector2D spaceDisc(StateVector2D & c, double nx, double ny);

	StateVector2D calcFlux(StateVector2D & c, double nx, double ny);

	StateVector2D calcPhysFluxesXi();

	void setConsInlet(double p, double u, double v, double T);

	void calcFluxesXi();

	StateVector2D calcPhysFlux(StateVector2D & c, double nx, double ny);

	void calcFluxesEta();

	StateVector2D calcFlux(StateVector2D & c_left_left, StateVector2D & c_left, StateVector2D & c_right, StateVector2D & c_right_right, double nx, double ny);

	StateVector2D calcPhysFlux(StateVector2D & prim_left, StateVector2D & prim_right, double nx, double ny);

	StateVector2D calcDissip(StateVector2D & prim_left, StateVector2D & prim_right, double nx, double ny);

	void physicalFluxRight(StateVector2D &c);
	void physicalFluxUp(StateVector2D &c);

	void RoeDissipRight(int i, int j);

	double calcTimeStep();

	double calcResidualL2(StateMatrix2D & o, StateMatrix2D & n);

	void executeTimeStepGlobal();

	void writeSolution(std::string filename);

	StateVector2D calcRoeDissip(StateVector2D & prim_left, StateVector2D & prim_right, double nx, double ny);
	StateVector2D mainloop();
	StateVector2D calcDissip();

	void setInitialCondition();

	double getResidual();
};

