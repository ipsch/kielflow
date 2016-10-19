#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "fftw3.h"   // google it (something about discret fourier Transformations)
#include "field.hpp"

#include <iostream>
#include <cmath>
#include "o_math.hpp"


#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

void handle_dealiasing(field_imag &DEST);
void dealiasing_36er(field_imag &DEST);
void dealiasing_23rd(field_imag &DEST);
void dealiasing_theta(field_imag &DEST, const double &x = 0.66666);

void FFT (const field_real &in, field_imag &out);
void iFFT(const field_imag &in, field_real &out);

double supremum(const field_real &in);
double supremum(const field_imag &in);

#endif // ende der funktionssammlung für Feldoperationen


