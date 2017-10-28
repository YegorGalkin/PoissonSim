#pragma once

enum BoundaryConditions
{
	PERIODIC,
	KILLING,
	REFLECTIVE
};
struct Cell
{

	vector<double> coords_x;
	vector<double> coords_y;
	vector<double> death_rates;

	Cell() {}
};
