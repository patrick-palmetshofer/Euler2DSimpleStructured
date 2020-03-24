#include "LimiterMCvanLeer.h"
#include <cmath>
#include <algorithm>



LimiterMCvanLeer::LimiterMCvanLeer()
{
}


LimiterMCvanLeer::~LimiterMCvanLeer()
{
}

PrimitiveVariables LimiterMCvanLeer::calcPhis(PrimitiveVariables &rs) const
{
	PrimitiveVariables phis;
	for (auto i = 0; i < rs.size(); i++)
	{
		double r = rs[i];
		double phi = std::min(2 * r, 0.5*(1 + r));
		phi = std::min(phi, 2.0);
		phi = std::max(phi, 0.0);
		phis[i] = phi;
	}
	return phis;
}
