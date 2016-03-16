#ifndef FIELD_REAL_HPP_
#define FIELD_REAL_HPP_

#include "field_base.hpp"
#include "interface_3d_fkt.hpp"

class field_real : public field
{
public :
	field_real(const grid_Co &Omega);
	field_real(const grid_Fo &FOmega);
	field_real(const field_real &that);
	~field_real();

	void fill(double (*fill_fkt)(const double &,const double &, const double &)) const;
	void fill(interface_3d_fkt const &rhs);
	void add(double (*fill_fkt)(const double &,const double &, const double &)) const;
	void mulitply(const double &lambda);

	double& operator() (const int &ix, const int &iy, const int &iz); // Subscript operators often come in pairs
	double  operator() (const int &ix, const int &iy, const int &iz) const; // Subscript operators often come in pairs
	double  operator() (const double &x, const double &y, const double &z) const;
	double& operator() (int N); // Subscript operators often come in pairs
	//double  operator() (int N) const; // Subscript operators often come in pairs

	grid_Co * my_grid;
	double * val;
	const int &N;
	const int &Nx;
	const int &Ny;
	const int &Nz;

	//void resize(const int &NX, const int &NY, const int &NZ);
	void resize(int NX, int NY, int NZ);
	// ToDO : functoren implementieren

	void swap(field_real &that);
	field_real& operator= (const field_real &rhs);
	field_real& operator+= (const field_real &rhs);

	//field_real& operator/= (double const& lambda);
};





#endif
