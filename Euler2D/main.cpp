#include <iostream>

#include "Grid.h"
#include "TimeStepperExplicit.h"
#include "Solution.h"
#include "StructuredBoundary.h"

int main()
{
	Grid grid;
	grid.readGrid("grid.su2");
	SlipWallStructured lowerwall();
	SlipWallStructured upperwall();
	SupersonicInletStructured inlet();
	SupersonicOutletStructured outlet();

	grid.setLeftBoundaryCondition(&inlet);
	grid.setRightBoundaryCondition(&outlet);
	grid.setUpBoundaryCondition(&upperwall);
	grid.setDownBoundaryCondition(&lowerwall);

	TimeStepperExplicit timestepper(grid);
	Solution solution(grid);

	while (timestepper.getTime() < 10)
	{
		timestepper.executeTimeStepGlobal();
	}

	solution.write("Result.txt");
}