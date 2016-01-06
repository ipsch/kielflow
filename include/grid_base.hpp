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



	int N;
	int Nx;
	int Ny;
	int Nz;

	axis * x_axis;
	axis * y_axis;
	axis * z_axis;

	grid(const int &Nx_, const int &Ny_, const int &Nz_,
			const axis &Ax, const axis &Ay, const axis &Az,
			const ptr_axis_factory_fkt &factory_fkt);
	grid(const grid &that);


	virtual ~grid();
	virtual  int index_at(const int &i, const int &j, const int &k) const = 0 ;
	virtual void ijk_at(const int &index, int &i, int &j, int &k) const = 0;
	virtual grid * clone() const = 0;


};








inline double Norm_Euklid(const double &x, const double &y, const double &z)
{
	return sqrt(x*x+y*y+z*z);
}




#endif
