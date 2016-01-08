#include "effect_translation.hpp"


effect_translation::effect_translation(const double &v, const direction &e_i) :
	my_velocity(v), my_direction(e_i)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("effect_translation::effect_translation(const double &v, const direction &e_i)");
	my_log << v;
	my_log << e_i;
   #endif
}

void effect_translation::execute(const field_imag &in, field_imag &out)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("effect_translation");
	my_log << "start";
   #endif

	field_imag dX_in(*in.my_grid);
	static OP_partial_derivative d_dX(in,my_direction);
	d_dX.execute(in,dX_in);

	for(int i=0; i<in.N; ++i)
	{
		out.val[i][0] -= my_velocity*dX_in.val[i][0];
		out.val[i][1] -= my_velocity*dX_in.val[i][1];
	}


   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}
