#pragma once
class TimeStep
{
public:
	TimeStep();
	virtual ~TimeStep();

	virtual void calculateTimeStep();
};

