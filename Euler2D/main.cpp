#include <iostream>
#include "Grid.h"
#include "TimeStepperExplicit.h"
#include "Solution.h"

int main()
{
	Grid grid;
	grid.readGrid("grid.su2");
	TimeStepperExplicit timestepper(grid);
	Solution solution(grid);

	while (timestepper.getTime() < 10)
	{
		timestepper.executeTimeStepGlobal();
	}

	solution.write("Result.txt");
}