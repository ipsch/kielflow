#ifndef EFFECT_DIFFUSION_U_HPP_
#define EFFECT_DIFFUSION_U_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"
#include "OP_partial_derivative.hpp"
#include "IO.hpp" // ToDo : remove this

class effect_diffusion_u : public effect
{

public :
	effect_diffusion_u(double temperature_ratio);
	void execute(const field_imag &Fni,
			field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz);

private :
	const double theta;
};




#endif
