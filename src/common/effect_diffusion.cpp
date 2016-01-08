#include "effect_diffusion.hpp"




effect_diffusion::effect_diffusion(const double &eta) :
	viscosity(eta)
{

}





void effect_diffusion::execute(const field_imag &in, field_imag &out) const
{

	field_imag d2dx(*in.my_grid);
	static OP_partial_derivative d_dx(in,e_x);
	d_dx.execute(in,d2dx);
	d_dx.execute(d2dx,d2dx);


	field_imag d2dy(*in.my_grid);
	static OP_partial_derivative d_dy(in,e_y);
	d_dy.execute(in,d2dy);
	d_dy.execute(d2dy,d2dy);


	field_imag d2dz(*in.my_grid);
	static OP_partial_derivative d_dz(in,e_z);
	d_dz.execute(in,d2dz);
	d_dz.execute(d2dz,d2dz);


	for(int i=0; i<in.N; ++i)
	{
		out.val[i][0] -= viscosity*( d2dx.val[i][1] + d2dy.val[i][1] + d2dz.val[i][1] );
		out.val[i][1] -= viscosity*( d2dx.val[i][1] + d2dy.val[i][1] + d2dz.val[i][1] );
	}

	return;
}
