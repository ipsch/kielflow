#ifndef EFFECT_FORCE_GRAD_HPP_
#define EFFECT_FORCE_GRAD_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

class effect_force_grad : public effect
{

public :
	effect_force_grad(double coupling);
	void execute(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
			field_imag &Buffer_Fni);
private :
	const double k;

};




#endif
