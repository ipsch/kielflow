#ifndef GRID_BASE_HPP_
#define GRID_BASE_HPP_


#include<cmath>
#include<iostream>

#include "axis.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif


class grid
{
public :
	axis * x_axis;
	axis * y_axis;
	axis * z_axis;

	grid(const axis &Ax, const axis &Ay, const axis &Az);
	grid(const grid &that);

	~grid();

	const int &Nx;
	const int &Ny;
	const int &Nz;

	void resize(const int &Nx, const int &Ny, const int &Nz)
	{
		x_axis->resize(Nx);
		y_axis->resize(Ny);
		z_axis->resize(Nz);
	}

};








inline double Norm_Euklid(const double &x, const double &y, const double &z)
{
	return sqrt(x*x+y*y+z*z);
}




#endif
