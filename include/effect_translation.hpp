#ifndef EFFECT_TRANSLATION_HPP_
#define EFFECT_TRANSLATION_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"


class effect_translation : public effect
// (d/dt) u = v * (d/d e_i) u
// translates a field "u" with constant velocity "v"
{
public :
	effect_translation(const double &v, OP_partial_derivative &Dx);
	void execute(const field_imag &in, field_imag &out);
private :
	double my_velocity;
	OP_partial_derivative &d_dx;

};



#endif
