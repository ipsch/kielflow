#ifndef EFFECT_SPECTRAL_VISCOSITY_HPP
#define EFFECT_SPECTRAL_VISCOSITY_HPP

#include <cmath>
#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"



#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

class effect_spectral_viscosity : public effect
{
public :
	effect_spectral_viscosity(const double &fraction, const double &L, const int &N);
	void execute(const field_imag &in, field_imag &out) const;
private :
	effect_spectral_viscosity() : eps(0.), cutoff(0.) { }
	double filter(const double &k2) const;
	double eps;
	int cutoff;
};








#endif
