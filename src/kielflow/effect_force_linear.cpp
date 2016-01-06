#include "effect_force_linear.hpp"

effect_force_linear::effect_force_linear(const double &k) :
	hook_constant(k)
{

}


void effect_force_linear::execute(const field_imag &in, field_imag &out)
{
   #ifdef _MY_VERBOSE
	logger log("effect_force_linear");
	log <<"execute(...)";
   #endif

	for(int i=0; i<in.N; ++i)
	{
		out.val[i][0] += hook_constant*in.val[i][0];
		out.val[i][1] += hook_constant*in.val[i][1];
	}

  #ifdef _MY_VERBOSE
	log <<"done";
   #endif

	return;
}




