#pragma once
#include "BoundaryDirichlet.h"
#include "BoundaryNeumann.h"
class BoundaryMixed :
	public BoundaryNeumann, 
	public BoundaryDirichlet
{
public:
	BoundaryMixed();
	~BoundaryMixed();
};

