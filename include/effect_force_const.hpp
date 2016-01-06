#ifndef EFFECT_FORCE_CONST_HPP_
#define EFFECT_FORCE_CONST_HPP_

#include "effect_base.hpp"
#include "field.hpp"

class effect_force_const : public effect
{
public :
	effect_force_const(const double &k);
	void execute(const field_imag &in, field_imag &out) const;
private :
	double a0;
};

#endif
