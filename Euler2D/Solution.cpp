#include "Solution.h"
#include <fstream>
#include <iostream>



Solution::Solution(Grid & grid)
{
	cells = grid.getCells();
}


Solution::~Solution()
{
}

void Solution::write(std::string filename)
{
	std::ofstream stream;
	try
	{
		stream.open(filename, std::ifstream::out);
		PrimitiveVariables p;
		//Reserve add here
		for (auto cell : cells)
		{
			p = *(cell->getPrimitive());
			stream << p.rho << "\t" << p.u << "\t" << p.v << "\t" << p.e << "\n";
		}
		stream.close();
	}
	catch (std::ofstream::failure e) {
		std::cerr << "Exception writing file\n";
	}
}
