#ifndef EFFECT_FORCE_LINEAR_HPP_
#define EFFECT_FORCE_LINEAR_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"


class effect_force_linear : public effect
{


public :
	effect_force_linear(const double &k);
	void execute(const field_imag &in, field_imag &out);
private :
	double hook_constant;

};





#endif



