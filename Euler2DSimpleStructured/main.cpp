#include "Solver.h"
#include <iostream>

int main()
{
	Solver s;
	std::string meshpath = "../../../mesh/";
	std::string respath = "../../../solution/";
	std::string gridname = "Grid20deg";
	std::string infile = meshpath + gridname + ".grd";
	std::string outfile = respath + gridname + ".res";

	s.readGridGridPro(infile);
	s.setConsInlet(1e5, 6e2, 0, 300);
	s.setInitialCondition();
	int maxsteps = 100000;
	for (int i = 0; i<maxsteps;i++)
	{
		if ((i+1) % 10000 == 0)
		{
			outfile = respath + gridname + ".res";
			outfile = respath + gridname + std::to_string(i) + ".res";
			s.writeSolution(outfile);
			std::cout << s.getResidual() << std::endl;
			if (s.getResidual() < 1e-3)
				break;
		}
		s.executeTimeStepGlobal();
	}
	return 0;
}