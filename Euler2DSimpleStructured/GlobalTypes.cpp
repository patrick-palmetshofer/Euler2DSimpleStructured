#include "GlobalTypes.h"

StateVector2D operator* (double scalar, StateVector2D v)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = scalar * v[i];
	return res;
}

StateVector2D operator* (StateVector2D v, double scalar)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = scalar * v[i];
	return res;
}

StateVector2D operator/ (StateVector2D v, double scalar)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = v[i]/scalar;
	return res;
}

StateVector2D operator/ (StateVector2D v1, StateVector2D v2)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = v1[i] / v2[i];
	return res;
}

StateVector2D operator+ (StateVector2D v1, StateVector2D v2)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = v1[i] + v2[i];
	return res;
}

StateVector2D operator* (StateVector2D v1, StateVector2D v2)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = v1[i] * v2[i];
	return res;
}

StateVector2D operator- (StateVector2D v1, StateVector2D v2)
{
	StateVector2D res;
	for (int i = 0; i < 2 + 2; i++)
		res[i] = v1[i] - v2[i];
	return res;
}