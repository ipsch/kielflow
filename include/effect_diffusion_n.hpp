#ifndef EFFECT_DIFFUSION_N_HPP_
#define EFFECT_DIFFUSION_N_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

class effect_diffusion_n : public effect
{

public :
	effect_diffusion_n();
	void execute(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
			field_imag &Buffer_Fni);

};




#endif
