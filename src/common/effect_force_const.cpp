#include "effect_force_const.hpp"


effect_force_const::effect_force_const(const double &k) :
a0(k)
{

}


void effect_force_const::execute(const field_imag &in, field_imag &out) const
{

	static double N = sqrt(in.my_grid->x_axis->N*in.my_grid->y_axis->N*in.my_grid->z_axis->N);
	out.val[0][0] += a0*N;
	return;
}
