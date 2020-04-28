#include "Solver.h"
#include <ctime>
#include <iostream>
#include <fstream>

const double eps = 1;
const double kappa = 1.0/3.0;// 1.0 / 3.0;
const bool limit = true;

const double convcrit = 1e-10;


void writeResidual(std::vector<StateVector2D> &residuals, std::string &filename)
{
	std::ofstream stream;
	try
	{
		stream.open(filename, std::ofstream::out);
		stream << "rho" << ",\t" << "rhou" << ",\t" << "rhov" << ",\t" << "rhoet" << "\n";
		//Reserve add here
		for (int i = 0; i < residuals.size(); i++)
		{
			for (int k = 0; k < 3; k++)
			{
				stream << residuals[i][k] << ",\t";
			}
			stream << residuals[i][3] << "\n";
		}
		stream.flush();
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}

void executeSolver(double vel, int angle, double eps, double kappa, bool limit, double convcrit)
{
	std::string meshpath = "../../../mesh/";
	std::string respath = "../../../solution/";

	if (limit)
		respath += "withlimiter/";
	else
		respath += "nolimiter/";

	std::string gridname = "Grid" + std::to_string(angle) + "deg";
	std::string velstr = std::to_string(vel);
	std::string infile = meshpath + gridname + ".grd";
	std::string outfile = respath + gridname + velstr + ".res";
	std::string residualfile = respath + gridname + velstr + "_residualsL2.csv";
	std::string residualinffile = respath + gridname + velstr + "_residualsLinf.csv";

	Solver s(infile,eps,kappa);

	if (!limit)
		s.disableLimiter();

	s.setConsInlet(1.01325e5, vel, 0, 300);
	s.setConsInitial(1.01325e5, vel, 0, 300);

	int maxsteps = 10000;
	std::vector<StateVector2D> residuals(maxsteps+1);
	std::vector<StateVector2D> residualsinf(maxsteps + 1);

	int i = 0;
	for (i = 0; i <= maxsteps; i++)
	{
		s.executeTimeStepLocal();
		residuals[i] = s.getResidualsL2();
		residualsinf[i] = s.getResidualsLinfty();
		if (residuals[i][0] < convcrit)
			break;
		if (i % 1000 == 0)
		{
			std::cout << residuals[i][0] << "\t" << residuals[i][1] << "\t" << residuals[i][2] << "\t" << residuals[i][3] << "\n";
			//if (s.getResidual() < 1e-9)
			//	break;
		}
	}


	s.writeSolution(outfile);
	residuals.resize(i+1);
	residualsinf.resize(i + 1);

	writeResidual(residuals, residualfile);
	writeResidual(residualsinf, residualinffile);
}

int main()
{
	// 555.50,867.97,1041.566
	std::vector<double> vels = { 1041.566, 867.97, 555.50 };
	std::vector<int> angles = { 10,20,30 };

	for (auto vel : vels)
	{
		for (auto angle : angles)
		{
			std::clock_t c_start = std::clock();
			executeSolver(vel, angle,eps,kappa,limit,convcrit);
			std::clock_t c_end = std::clock();
			double cpu_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
			std::cout << "Run complete after " << cpu_time << "ms on CPU for v=" << std::to_string(vel) << " and alpha=" << std::to_string(angle) << "\n";
		}
	}

	return 0;
}




//int main()
//{
//	Solver s;
//	std::string meshpath = "../../../mesh/";
//	std::string respath = "../../../solution/";
//	std::string gridname = "sodX";
//	std::string infile = meshpath + gridname + ".grd";
//	std::string outfile = respath + gridname + ".res";
//
//	s.readGridGridPro(infile);
//	s.setSodXInitial();
//	int maxsteps = 1000;
//	for (int i = 1; i<=maxsteps;i++)
//	{
//		s.executeTimeStepGlobal();
//		if (i % 10 == 0)
//		{
//			outfile = respath + gridname + ".res";
//			outfile = respath + gridname + std::to_string(i) + ".res";
//			s.writeSolution(outfile);
//			std::cout << s.getResidual() << "\n";
//			//if (s.getResidual() < 1e-9)
//			//	break;
//		}
//	}
//	return 0;
//}
//
//int main()
//{
//	Solver s;
//	std::string meshpath = "../../../mesh/";
//	std::string respath = "../../../solution/";
//	std::string gridname = "Grid20deg";
//	std::string infile = meshpath + gridname + ".grd";
//	std::string outfile = respath + gridname + ".res";
//
//	s.readGridGridPro(infile);
//	// 555.50,867.97,1041.566
//	double vel = 555.50;
//	s.setConsInlet(1.01325e5, vel, 0, 300);
//	s.setConsInitial(1.01325e5, vel, 0, 300);
//	int maxsteps = 10000;
//	for (int i = 1; i <= maxsteps; i++)
//	{
//		s.executeTimeStepLocal();
//		if (i % (maxsteps / 1000) == 0)
//		{
//			outfile = respath + gridname + ".res";
//			outfile = respath + gridname + std::to_string(i) + ".res";
//			s.writeSolution(outfile);
//			std::cout << i << "\t" << s.getTime() << "\t" << s.getResiduals()[0] << "\n";
//			if (s.getResiduals()[0] < 1e-6)
//				break;
//		}
//	}
//	return 0;
//}