#include "effect_penalization.hpp"


effect_penalization::effect_penalization(const double &eta, field_real &H, const double &boundary_value) :
damping(eta), mask(H), fix_value(boundary_value)
{
	if(damping>0)
	{
		std::cout << "damping in penalization method is positiv;\n";
		std::cout << "method would blow up - stopping\n";
		throw("error in penalization method due to damping");
	}

   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("effect_penalization");
	my_log <<"effect_penalization(...)";
   #endif
}


void effect_penalization::execute(const field_imag &in, field_imag &out) const
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("effect_penalization");
	my_log << "execute(const field_imag &in, field_imag &out) const";
   #endif

	field_real f(*in.my_grid);
	field_real u(*in.my_grid);

	iFFT(in,u);
	iFFT(out,f);

	for(int i=0; i<f.N; ++i)
		f.val[i] = (1.-mask.val[i])*f.val[i] + mask.val[i]*damping*(u.val[i]-fix_value);

	FFT(f,out);

   #if defined(_MY_VERBOSE_MORE)
	my_log << "done";
   #endif

	return;
}

void effect_penalization::penalize_density(const field_imag &in, field_imag &out) const
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("effect_penalization");
	my_log << "execute(const field_imag &in, field_imag &out) const";
   #endif

	field_real  tmp(*in.my_grid);

	iFFT(in,tmp);

	for(int i=0; i<tmp.N; ++i)
		tmp.val[i] = tmp.val[i]*(1.- mask.val[i]/damping);

	FFT(tmp,out);

   #ifdef _MY_VERBOSE
	my_log << "done";
   #endif
}
