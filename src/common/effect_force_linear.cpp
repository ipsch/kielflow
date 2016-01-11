#include "effect_force_linear.hpp"

effect_force_linear::effect_force_linear(const double &k) :
	hook_constant(k)
{

}


void effect_force_linear::execute(const field_imag &in, field_imag &out)
{
   #if defined(_MY_VERBOSE)
	logger my_log("effect_force_linear");
	my_log <<"execute(...)";
   #endif

	for(int i=0; i<in.N; ++i)
	{
		out.val[i][0] += hook_constant*in.val[i][0];
		out.val[i][1] += hook_constant*in.val[i][1];
	}

  #if defined(_MY_VERBOSE_TEDIOUS)
	my_log <<"done";
   #endif

	return;
}




