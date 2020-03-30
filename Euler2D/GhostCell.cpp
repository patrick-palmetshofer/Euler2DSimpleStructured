#include "GhostCell.h"



GhostCell::GhostCell()
{
}


GhostCell::GhostCell(double volume, Face * boundary, Cell * neighbor, Config & cfg) : Cell(volume,boundary,boundary,boundary,boundary,cfg)
{
	this->boundary = boundary;
	this->boundary = boundary;
}

GhostCell::~GhostCell()
{
}

void GhostCell::setNeighbor(Cell * neighbor)
{
	this->neighbor = neighbor;
}
