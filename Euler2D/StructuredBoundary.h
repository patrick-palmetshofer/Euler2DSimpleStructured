#pragma once
#include <memory>
#include "PrimitiveVariables.h"
#include "Grid.h"

class StructuredBoundary
{
public:
	StructuredBoundary();
	~StructuredBoundary();
};

class SupersonicInletStructured : public StructuredBoundary
{
private:
	PrimitiveVariables &p;

public:
	SupersonicInletStructured();
	~SupersonicInletStructured();
};

class SupersonicOutletStructured : public StructuredBoundary
{
private:
};

class SlipWallStructured : public StructuredBoundary
{
private:

};
