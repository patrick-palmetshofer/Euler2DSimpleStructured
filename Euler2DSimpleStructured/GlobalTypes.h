#pragma once
#include <memory>
#include <vector>
#include <valarray>


//Generalized State Vector specified to 2 dimensions
template <uint ndim>
using StateVector = std::array<double,ndim+2>;
using StateVector2D = StateVector<2>;

//Matrices of StateVectors or Values defined on the grid
using ValueMatrix2D = std::vector<std::vector<double>>;
using StateMatrix2D = std::vector<std::vector<StateVector2D>>;
//using Matrix3D = std::vector<std::vector<std::vector<stateVector3D>>>;

//template <uint ndim>
//using Jacobian = std::array<std::array<double,ndim+2>,ndim+2>;

//using Jacobian2D = Jacobian<2>;

//template <uint ndim>
//StateVector<ndim> operator* (StateVector<ndim> &v1, StateVector<ndim> &v2)
//{
//	StateVector<ndim> res;
//	for (int i = 0; i < ndim+2; i++)
//		res[i] = v1[i] * v2[i];
//	return res;
//}
//
//template <uint ndim>
//StateVector<ndim> operator+ (StateVector<ndim> &v1, StateVector<ndim> &v2)
//{
//	StateVector<ndim> res;
//	for (int i = 0; i < ndim+2; i++)
//		res[i] = v1[i] + v2[i];
//	return res;
//}
//
//template <uint ndim>
//StateVector<ndim> operator+= (StateVector<ndim> &v1, StateVector<ndim> &v2)
//{
//	v1 = v1 + v2;
//}
//
//template <uint ndim>
//StateVector<ndim> operator- (StateVector<ndim> &v1, StateVector<ndim> &v2)
//{
//	StateVector<ndim> res;
//	for (int i = 0; i < ndim+2; i++)
//		res[i] = v1[i] - v2[i];
//	return res;
//}
//
//template <uint ndim>
//StateVector<ndim> operator* (double scalar, StateVector<ndim> &v)
//{
//	StateVector<ndim> res;
//	for (int i = 0; i < ndim+2; i++)
//		res[i] = v[i]*scalar;
//	return res;
//}
//
//template <uint ndim>
//StateVector<ndim> operator* (StateVector<ndim> &v, double scalar)
//{
//	StateVector<ndim> res;
//	for (int i = 0; i < ndim+2; i++)
//		res[i] = v[i] * scalar;
//	return res;
//}
//
//template <uint ndim>
//StateVector<ndim> operator* (Jacobian<ndim> J, StateVector<ndim> &v)
//{
//	StateVector<ndim> res;
//	for (int i = 0; i < ndim + 2; i++)
//	{
//		res[i] = 0;
//		for (int j = 0; j < ndim + 2; j++)
//		{
//			res[i] += J[i, j] * v[j];
//		}
//	}
//	return res;
//}
//
//template <uint ndim>
//double sum(StateVector<ndim> &v)
//{
//	double res = 0;
//	for (int i = 0; i < ndim + 2; i++)
//	{
//		res += v[i];
//	}
//	return res;
//}



//Operator overloads for all StateVectors. Implements easy arithmetics
StateVector2D operator* (double scalar, StateVector2D v);
StateVector2D operator* (StateVector2D v, double scalar);
StateVector2D operator/ (StateVector2D v1, StateVector2D v2);
StateVector2D operator/ ( StateVector2D v, double scalar);
StateVector2D operator+ (StateVector2D v1, StateVector2D v2);
StateVector2D operator* (StateVector2D v1, StateVector2D v2);
StateVector2D operator- (StateVector2D v1, StateVector2D v2);