#ifndef EFFECT_FORCE_E_HPP_
#define EFFECT_FORCE_E_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

class effect_force_E : public effect
{

public :
	effect_force_E();
	void execute(const field_imag &FPhi,
			field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz);

};

#endif
