#ifndef EFFECT_SPECTRAL_VISCOSITY_HPP
#define EFFECT_SPECTRAL_VISCOSITY_HPP

#include <iostream> // ToDO : remove
#include <cmath>
#include "effect_base.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "operations.hpp"



#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

class effect_spectral_viscosity : public effect
{
public :
	effect_spectral_viscosity(const double &fraction, const grid_Fo &Omega);
	void execute(const field_imag &in, field_imag &out) const;
private :
	effect_spectral_viscosity() : eps(0.), cutoff(0.), kx_max(0.), ky_max(0.), kz_max(0.) { }
	double filter(const double &kx, const double &ky, const double &kz) const;
	double eps;
	double cutoff;
	double kx_max;
	double ky_max;
	double kz_max;
};








#endif
