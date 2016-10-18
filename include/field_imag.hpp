#ifndef FIELD_IMAG_HPP_
#define FIELD_IMAG_HPP_

#include "field_base.hpp"
#include "fftw3.h"

class field_imag : public field
{
public :
	field_imag(const grid &Omega);
	field_imag(const field_imag &that);
	~field_imag();

	void fill_RE(double (*fill_fkt)(const double &, const double &, const double &)) const;
	void fill_IM(double (*fill_fkt)(const double &, const double &, const double &)) const;
	bool IsNan();

	fftw_complex * val;

	field_imag& operator+=(const field_imag& rhs);
	field_imag& operator-=(const field_imag& rhs);


	// ToDo : field_imag functoren implementieren
	fftw_complex& operator() (int ix, int iy, int iz); // Subscript operators often come in pairs
	//fftw_complex& operator() (int ix, int iy, int iz) const; // Subscript operators often come in pairs
	fftw_complex& operator() (const int &N); // Subscript operators often come in pairs
	//fftw_complex& operator() (const int &N) const; // Subscript operators often come in pairs
	double& operator() (int ix, int iy, int iz, int re_im); // Subscript operators often come in pairs
	//fftw_complex& operator() (int ix, int iy, int iz) const; // Subscript operators often come in pairs
	double& operator() (const int &N, int re_im); // Subscript operators often come in pairs
	//fftw_complex& operator() (const int &N) const; // Subscript operators often come in pairs



	/*
	ToDo : field_imag Operatoren implementieren
	field_imag& operator= (field_imag const& rhs);
	field_imag& operator/= (double const& lambda);
	// #### Operatoren (real) ##########################
	real_field operator+(const real_field &c1, const real_field &c2);
	real_field operator-(const real_field &c1, const real_field &c2);
	real_field operator*(const real_field &c1, const real_field &c2);
	real_field operator*(const double &lambda, const real_field &in);
	real_field operator*( const real_field &in, const double &lambda);

	// #### Operatoren (imag) ##########################
	field_imag operator+(const field_imag &c1, const field_imag &c2);
	field_imag operator-(const field_imag &c1, const field_imag &c2);
	field_imag operator*(const field_imag &c1, const field_imag &c2);
	field_imag operator*(const double &lambda, const field_imag &in);
	field_imag operator*(const field_imag &in, const double &lambda);
*/
};









#endif
