#include "Solver.h"
#include <iostream>

int main()
{
	Solver s;
	std::string meshpath = "../../../mesh/";
	std::string respath = "../../../solution/";
	std::string gridname = "sodX";
	std::string infile = meshpath + gridname + ".grd";
	std::string outfile = respath + gridname + ".res";

	s.readGridGridPro(infile);
	s.setSodXInitial();
	int maxsteps = 1000;
	for (int i = 1; i<=maxsteps;i++)
	{
		s.executeTimeStepGlobal();
		if (i % 1 == 0)
		{
			outfile = respath + gridname + ".res";
			outfile = respath + gridname + std::to_string(i) + ".res";
			s.writeSolution(outfile);
			std::cout << s.getResidual() << "\n";
			//if (s.getResidual() < 1e-9)
			//	break;
		}
	}
	return 0;
}

//int main()
//{
//	Solver s;
//	std::string meshpath = "../../../mesh/";
//	std::string respath = "../../../solution/";
//	std::string gridname = "Grid10deg";
//	std::string infile = meshpath + gridname + ".grd";
//	std::string outfile = respath + gridname + ".res";
//
//	s.readGridGridPro(infile);
//	s.setConsInlet(1e5, 4e2, 0, 300);
//	s.setConsInitial(1e5, 4e2, 0, 300);
//	int maxsteps = 10000;
//	for (int i = 1; i<=maxsteps;i++)
//	{
//		s.executeTimeStepGlobal();
//		if (i % 10 == 0)
//		{
//			outfile = respath + gridname + ".res";
//			outfile = respath + gridname + std::to_string(i) + ".res";
//			s.writeSolution(outfile);
//			std::cout << s.getResidual() << "\n";
//			if (s.getResidual() < 1e-6)
//				break;
//		}
//	}
//	return 0;
//}