#ifndef EFFECT_CONTINUITYEQ_HPP_
#define EFFECT_CONTINUITYEQ_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

class effect_continuityEQ : public effect
{

public :
	effect_continuityEQ() { };
	void execute(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni,
			field_imag &Buffer_Fni);

};




#endif
