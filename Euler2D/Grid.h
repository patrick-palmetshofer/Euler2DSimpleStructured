#pragma once
#include "Cell.h"
#include "GhostCell.h"
#include <string>
#include <vector>
#include <memory>

class Grid
{
protected:
	std::vector<std::vector<std::valarray<double>>> points;
	std::vector<std::vector<std::shared_ptr<Cell>>> cells;
	std::vector<std::vector<std::shared_ptr<Face>>> xi_faces;
	std::vector<std::vector<std::shared_ptr<Face>>> eta_faces;
	std::vector<std::shared_ptr<GhostCell>> boundaryCells;

	std::vector<std::shared_ptr<Cell>> all_cells;
	std::vector<std::shared_ptr<Cell>> all_faces;

	Config *cfg;

	double nxi_points;
	double neta_points;
	const int dim = 2;
public:
	Grid();
	~Grid();
	 
	void readGrid(std::string filename);
	void constructCells();
	std::vector<std::shared_ptr<Cell>> getCells();
	std::vector<std::shared_ptr<Face>> getFaces();

};

