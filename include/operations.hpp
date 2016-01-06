#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "fftw3.h"   // google it (something about discret fourier Transformations)
#include "field.hpp"

#include <iostream>
#include <cmath>

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif

void handle_dealiasing(field_imag &DEST);




class OP_partial_derivative
{
public :
	OP_partial_derivative(const field_real &field_caller, const direction &e_i);
	OP_partial_derivative(const field_imag &field_caller, const direction &e_i);
	~OP_partial_derivative();
	void execute(const field_imag &in, field_imag &out);
private :
	bool IsLinear;
	int my_i, my_j, my_k;
	int * my_X;
	const axis * my_axis; // which is really a "axis_Co"
};



void FFT (const field_real &in, field_imag &out);
void iFFT(const field_imag &in, field_real &out);

double supremum(const field_real &in);
double supremum(const field_imag &in);

#endif // ende der funktionssammlung für Feldoperationen


