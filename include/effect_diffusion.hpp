#ifndef EFFECT_DIFFUSION_HPP_
#define EFFECT_DIFFUSION_HPP_

#include "effect_base.hpp"
#include "field.hpp"
#include "operations.hpp"

class effect_diffusion : public effect
{

public :
	effect_diffusion(const double &eta);
	void execute(const field_imag &in, field_imag &out) const;
private :
	double viscosity;
};









#endif
