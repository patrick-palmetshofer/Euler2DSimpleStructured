#pragma once
#include <memory>
#include <vector>

using stateVector2D = std::array<double,4>;

//using stateVector3D = std::array<double,5>;

using Matrix2D = std::vector<std::vector<stateVector2D>>;
//using Matrix3D = std::vector<std::vector<std::vector<stateVector3D>>>;

Matrix2D conservative;
Matrix2D primitive;
Matrix2D rightstate;
Matrix2D upstate;
Matrix2D xi_fluxes;
Matrix2D eta_fluxes;

