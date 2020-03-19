#include "TimeStepperExplicit.h"



TimeStepperExplicit::TimeStepperExplicit()
{
}


TimeStepperExplicit::~TimeStepperExplicit()
{
}

void TimeStepperExplicit::executeTimeStep()
{
	for (auto rowit = structuredCells.begin(); rowit != structuredCells.end(); rowit++)
	{
		auto &row = *rowit;
		for (auto cellit = row.begin(); cellit != row.end(); cellit++)
		{
			auto &cell = *cellit;
			ConservativeVariables consold = cell.getConservative();
			ConservativeVariables consnew = ;
			cell.setConservative(consnew);
		}
	}

}
