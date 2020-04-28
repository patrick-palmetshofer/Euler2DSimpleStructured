#pragma once
#include "GlobalTypes.h"

class Solver
{
private:
	//Numbers of cells in both computational directions, excluding ghost cells
	int nxi_cells;
	int neta_cells;

	//Use limiter?
	bool limit;

	//Maximum CFL number, use 0.5 for sure stability
	const double maxCFL = 0.5;
	double p_infty;

	//2D Vector of points. 
	std::vector<std::vector<std::array<double, 2>>> points;

	//Time step and current time in simulation
	double dt;
	double time;

	//Volumes of cells in domain
	ValueMatrix2D volumes;

	//conservative variables for Inlet and initial condition
	StateVector2D cons_inlet;
	StateVector2D cons_initial;
	//StateMatrix2D primitive;

	//MUSCL parameters
	double reconstruct_eps;
	double reconstruct_kappa;

	//Total residual. Can be calculated through L2 or Linfty. See corresponding methods
	StateVector2D residuals_L2;
	StateVector2D residuals_Linfty;

	//Face areas in both computational directions
	ValueMatrix2D Sxis;
	ValueMatrix2D Setas;

	//Face normal vector components for each computational direction (2 physical directions each)
	ValueMatrix2D nxi_xs;
	ValueMatrix2D nxi_ys;
	ValueMatrix2D neta_xs;
	ValueMatrix2D neta_ys;

	//StateMatrix2D xi_face_prim;
	//StateMatrix2D eta_face_prim;

	//Fluxes at Faces and Conservative variables of all cells
	StateMatrix2D xi_fluxes;
	StateMatrix2D eta_fluxes;
	StateMatrix2D conservative;

	//Perfect gas constants
	double gamma = 1.4;
	double cp = 1005;

public:
	//Utilities for Sod problem tests
	void setSodXInitial();
	void setSodYInitial();

	//Constructors/destructors
	Solver();
	Solver(double eps, double kappa);
	Solver(std::string filename, double eps, double kappa);
	Solver(std::string filename);
	~Solver();

	//Reads GridPro grid file (Warning, no treatment of Gridpro Boundary Conditions)
	void readGridGridPro(std::string filename);

	//Utilities for conversion between variable types
	StateVector2D prim2cons(StateVector2D &p);
	StateVector2D cons2prim(StateVector2D &c);
	StateVector2D user2cons(double p, double u, double v, double T);
	double calcPcons(StateVector2D & c);
	double calcMacons(StateVector2D & c);
	double calcSoundSpeedcons(StateVector2D & c);
	double calcPprim(StateVector2D & prim);

	//Limiter function using the ratio. Can lead to problems with 0/0 inputs
	StateVector2D limiter(StateVector2D &r);

	//Limiter function using explicit values. Returns Limiter function in terms of r
	StateVector2D limiter(StateVector2D & x, StateVector2D & y);

	//Sets of limiter function called by previous function
	StateVector2D limiterMinmod(StateVector2D &r);
	StateVector2D limiterMinmod(StateVector2D & x, StateVector2D & y);
	StateVector2D limiterMC(StateVector2D & rs);

	//Treatment of Boundary conditions. Hard coded :(
	//First-order (constant boundary) conditions
	void setBoundaryInletLeft();
	void setBoundaryLowerWall();
	void setBoundaryUpperLowerWalls();
	void setWalls();
	void setBoundaryOutlet();
	void setBoundaryUpperOutlet();
	//LODI Boundary conditions
	void setCharacteristicBoundaryRightOutlet();
	void setCharacteristicBoundaryUpperOutlet();

	//Set initial and boundary condition values for initialization
	void setConsInlet(double p, double u, double v, double T);
	void setConsInitial(double p, double u, double v, double T);

	//Fill all cells with initial condition
	void setInitialCondition();

	//Calculate fluxes of all faces.
	void calcFluxesXi(); //Xi=const faces
	void calcFluxesEta(); //Eta=const faces

	//Flux calculation
	StateVector2D calcFluxMUSCL(StateVector2D & c_left_left, StateVector2D & c_left, StateVector2D & c_right, StateVector2D & c_right_right, double nx, double ny);
	StateVector2D calcFlux(StateVector2D & c_left, StateVector2D & c_right, double nx, double ny);

	//Previous method calls:
	//Calculate physical fluxes without numerical dissipation
	StateVector2D calcPhysFlux(StateVector2D & c, double nx, double ny);
	//Calculate numerical dissipation
	StateVector2D calcDissip(StateVector2D & prim_left, StateVector2D & prim_right, double nx, double ny);
	//Options for numerical dissipation calculation: Roe flux (no entropy correction)
	StateVector2D calcRoeDissip(StateVector2D & prim_left, StateVector2D & prim_right, double nx, double ny);





	//Calculate time step for global time stepping
	double calcTimeStep();
	//Method calls:
	//Calculate timestep for local timestepping
	double calcTimeStep(int i, int j);

	//Calculate residuals, total, not per equation in terms of old conservatives and new conservatives
	double calcResidualL2(StateMatrix2D & o, StateMatrix2D & n);
	double calcResidualLinfty(StateMatrix2D & o, StateMatrix2D & n);

	//Calculate residuals per value
	StateVector2D calcResidualsL2(StateMatrix2D & o, StateMatrix2D & n);
	StateVector2D calcResidualsLinfty(StateMatrix2D & o, StateMatrix2D & n);

	//Actual time stepping routines. To be called on each iteration.
	//Calls BCs, timestep calculator, Flux calculators and performs time integration. 
	//Saves residual
	void executeTimeStepGlobal();
	void executeTimeStepLocal();

	//Get current time for global time stepping. Returns 0 for local time steps
	double getTime();

	StateVector2D getResidualsL2();
	StateVector2D getResidualsLinfty();

	//Write solution to file. Same format as input grid
	void writeSolution(std::string filename);

	//limiter control
	void enableLimiter() {
		limit = true;
	};
	void disableLimiter() {
		limit = false;
	};;

	// Deprecated prototypes
	//StateVector2D calcFlux();
	//StateVector2D spaceDisc(int i, int j);
	//StateVector2D spaceDisc(StateVector2D & c, double nx, double ny);
	//StateVector2D calcFlux(StateVector2D & c, double nx, double ny);
	//StateVector2D calcPhysFluxesXi();
	//StateVector2D mainloop();
	//StateVector2D calcDissip();
	//void physicalFluxRight(StateVector2D &c);
	//void physicalFluxUp(StateVector2D &c);
	//void RoeDissipRight(int i, int j);
	//StateVector2D calcPhysFlux(StateVector2D & prim_left, StateVector2D & prim_right, double nx, double ny);
};

