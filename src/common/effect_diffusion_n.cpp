#include "effect_diffusion_n.hpp"

effect_diffusion_n::effect_diffusion_n()
{

}

void effect_diffusion_n::execute(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
		field_imag &Buffer_Fni)
	// - grad(u)
{
	int N = FUx.N;

	static OP_partial_derivative d_dx(FUx,e_x);
	static OP_partial_derivative d_dy(FUx,e_y);
	static OP_partial_derivative d_dz(FUx,e_z);

	field_imag FBuffer_x(FUx.my_grid);
	field_imag FBuffer_y(FUy.my_grid);
	field_imag FBuffer_z(FUz.my_grid);

	d_dx(FUx,FBuffer_x);
	d_dy(FUy,FBuffer_y);
	d_dz(FUz,FBuffer_z);

	for(int i=0; i<N;++i)
	{
		Buffer_Fni.val[i][0] += -(FBuffer_x.val[i][0] + FBuffer_y.val[i][0] + FBuffer_z.val[i][0]);
		Buffer_Fni.val[i][1] += -(FBuffer_x.val[i][1] + FBuffer_y.val[i][1] + FBuffer_z.val[i][1]);
	}

	return;
}
