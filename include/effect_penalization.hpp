#ifndef EFFECT_PENALIZATION_HPP
#define EFFECT_PENALIZATION_HPP


#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif

class effect_penalization : public effect
{
public :
	effect_penalization(const double &eta, field_real &H, const double &boundary_value);
	void execute(const field_imag &in, field_imag &out) const;
	void penalize_density(const field_imag &in, field_imag &out) const;

private :

	double damping;
	double fix_value;
	field_real &mask;
};


#endif
