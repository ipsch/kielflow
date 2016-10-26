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
template<typename Func>
void dealiasing_undesignated(field_imag &DEST, Func filter)
{
	double kx_max = DEST.my_grid.x_axis->k_val_at(DEST.Nx/2);
	double ky_max = DEST.my_grid.y_axis->k_val_at(DEST.Ny/2);
	double kz_max = DEST.my_grid.z_axis->k_val_at(DEST.Nz/2);

	for(int i=0;i<DEST.Nx; ++i)
	{
		double kx = DEST.my_grid.x_axis->k_val_at(i);
		for(int j=0;j<DEST.Ny; ++j)
		{
			double ky = DEST.my_grid.y_axis->k_val_at(j);
			for (int k=0; k< DEST.Nz; ++k)
			{
				double kz = DEST.my_grid.z_axis->k_val_at(k);

				int ijk = DEST.index(i,j,k);
				double kabs = sqrt(pow(kx/kx_max,2.)+pow(ky/ky_max,2.)+pow(kz/kz_max,2.));
				DEST.val[ijk][0] *= filter(kabs);
				DEST.val[ijk][1] *= filter(kabs);

			}
		}
	}
	return;
}

void FFT (const field_real &in, field_imag &out);
void iFFT(const field_imag &in, field_real &out);

double supremum(const field_real &in);
double supremum(const field_imag &in);

#endif // ende der funktionssammlung für Feldoperationen


