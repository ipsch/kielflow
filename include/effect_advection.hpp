#ifndef EFFECT_ADVECTION_HPP_
#define EFFECT_ADVECTION_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

class effect_advection : public effect
{

public :
	effect_advection();
	void execute(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz, const field_imag FXX,
					field_imag &return_Buffer);

};




#endif
