#ifndef FIELD_BASE_HPP_
#define FIELD_BASE_HPP_

#include <iostream>

#include "fftw3.h"
#include "grid.hpp"
#include "masks.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif

class field
{
public :
	int Nx;
	int Ny;
	int Nz;
	int N;
	grid my_grid;
	field(const grid &domain, const double &Nx_, const double &Ny_, const double &Nz_) :
		my_grid(domain), Nx(Nx_), Ny(Ny_), Nz(Nz_), N(Nx*Ny*Nz)
		{	}
	~field()
		{	};
	int index(const int &i, const int &j, const int &k) const
		{return (k + Nz*(j + i*Ny));}
private :
	//grid my_grid_;
protected :
	field();
};

/*
fftw_complex& operator=(const fftw_complex & rhs)
{
	fftw_complex lhs;
	lhs[0] = rhs[0];
	lhs[1] = rhs[1];
	return *lhs;
}
*/

/*
double& operator=(const double & rhs)
{
	fftw_complex lhs;
	lhs[0] = rhs[0];
	lhs[1] = rhs[1];

	double z =0;

	return z;
}
*/

/*
fftw_complex& operator =(fftw_complex &rhs)
{

}
*/











#endif
