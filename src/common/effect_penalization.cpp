#include "effect_penalization.hpp"


effect_penalization::effect_penalization(const double &eta, field_real &H, const double &boundary_value) :
damping(eta), mask(H), fix_value(boundary_value)
{
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

	field_real  tmp(*in.my_grid);
	field_imag Ftmp(*in.my_grid);

	iFFT(in,tmp);

	for(int i=0; i<tmp.N; ++i)
		tmp.val[i] = (tmp.val[i]-fix_value)*mask.val[i];

	FFT(tmp,Ftmp);

	for(int i=0; i<out.N; ++i)
	{
		out.val[i][0] += -damping*Ftmp.val[i][0];
		out.val[i][1] += -damping*Ftmp.val[i][1];
	}

   #ifdef _MY_VERBOSE
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
