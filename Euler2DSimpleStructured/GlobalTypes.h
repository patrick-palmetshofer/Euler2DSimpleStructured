#pragma once
#include <memory>
#include <vector>
#include <valarray>
#include <array>

#include <iostream>

//Generalized State Vector specified to 2 dimensions
template <unsigned int ndim>
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

namespace Euler
{
	//Utilities to swap coordinates in Statevectors
	template<typename T>
	T swap(T &data, int ind1, int ind2)
	{
		if (ind1 == ind2)
			return data;
		T res = data;
		res[ind1] = data[ind2];
		res[ind2] = data[ind1];
		return res;
	}

	template<typename T>
	T swap(T &data, int i)
	{
		return Euler::swap(data, 1, i);
	}


	//Checks for errors in matrices, used for debugging
	inline bool checkNaN(StateMatrix2D &m)
	{
		bool ret = false;
		for (int i = 0; i < m.size(); i++)
		{
			for (int j = 0; j < m[i].size(); j++)
			{
				for (int k = 0; k < m[i][j].size(); k++)
				{
					if (!std::isfinite(m[i][j][k]))
					{
						ret = true;
						std::cout << "Calculation aborted due to NaN value at i=" << i << " j=" << j << "\n";
						//throw;
					}
				}
			}
		}
		return ret;
	}

	inline bool checkNaN(StateVector2D &m)
	{
		bool ret = false;
		for (int k = 0; k < m.size(); k++)
		{
			if (!std::isfinite(m[k]))
			{
				ret = true;
				//throw;
			}
		}
		return ret;
	}
}


//namespace std {
//	template<class T> struct _Unique_if {
//		typedef unique_ptr<T> _Single_object;
//	};
//
//	template<class T> struct _Unique_if<T[]> {
//		typedef unique_ptr<T[]> _Unknown_bound;
//	};
//
//	template<class T, size_t N> struct _Unique_if<T[N]> {
//		typedef void _Known_bound;
//	};
//
//	template<class T, class... Args>
//	typename _Unique_if<T>::_Single_object
//		make_unique(Args&&... args) {
//		return unique_ptr<T>(new T(std::forward<Args>(args)...));
//	}
//
//	template<class T>
//	typename _Unique_if<T>::_Unknown_bound
//		make_unique(size_t n) {
//		typedef typename remove_extent<T>::type U;
//		return unique_ptr<T>(new U[n]());
//	}
//
//	template<class T, class... Args>
//	typename _Unique_if<T>::_Known_bound
//		make_unique(Args&&...) = delete;
//}