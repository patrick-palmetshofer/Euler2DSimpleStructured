#include "LimiterMinmod.h"



PrimitiveVariables LimiterMinmod::calcPhis(PrimitiveVariables &rs) const
{
	PrimitiveVariables phis;
	for (auto i = 0; i < rs.size(); i++)
	{
		double r = rs[i];
		double phi = std::min(1.0, r);
		phi = std::max(phi, 0.0);
		phis[i] = phi;
	}
	return phis;
}

LimiterMinmod::LimiterMinmod()
{
}


LimiterMinmod::~LimiterMinmod()
{
}
