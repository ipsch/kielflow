#include "effect_diffusion_u.hpp"



effect_diffusion_u::effect_diffusion_u(double temperature_ratio) :
theta(temperature_ratio)

{


}


void effect_diffusion_u::execute(const field_imag &Fni,
		field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz)
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_diffusion_u(...)";
   #endif

	static OP_partial_derivative d_dx(Fni,e_x);
	static OP_partial_derivative d_dy(Fni,e_y);
	static OP_partial_derivative d_dz(Fni,e_z);

	field_imag FBuffer(*Fni.my_grid);
	field_real Buffer(*Fni.my_grid);


	field_real ni(*Fni.my_grid);
	iFFT(Fni,ni);
	/*


	subdim my_dim;
	my_dim.default_xpos = ni.Nx/2;
	my_dim.default_ypos = ni.Ny/2;
	my_dim.default_zpos = ni.Nz/2;
	my_dim.default_direction = 0;
	my_dim.default_plane = 2;

	save_2d("ni.txt" , ni, ni, ni, ni, ni, my_dim);
*/

	// x-Komponente
	d_dx.execute(Fni,FBuffer);
	iFFT(FBuffer, Buffer);


	for(int i=0; i<ni.N; ++i)
	{
		// ToDo : Add secure flag here
		if( (ni.val[i]<=-1.) && (Buffer.val[i]<=1.e-6) )
			Buffer.val[i] = 0.;
		Buffer.val[i] = Buffer.val[i]/(theta*(1.+ni.val[i]));
	}

	FFT(Buffer, FBuffer);

	for(int i=0; i<Buffer_FUx.N; i++)
	{
		Buffer_FUx.val[i][0] += -FBuffer.val[i][0];
		Buffer_FUx.val[i][1] += -FBuffer.val[i][1];
	}


	// y-Komponente
	d_dy.execute(Fni,FBuffer);
	iFFT(FBuffer, Buffer);


	for(int i=0; i<ni.N; ++i)
	{
		if( (ni.val[i]<=-1.) && (Buffer.val[i]<=1.e-6) )
			Buffer.val[i] = 0.;
		Buffer.val[i] = Buffer.val[i]/(theta*(1.+ni.val[i]));
	}

	FFT(Buffer, FBuffer);

	for(int i=0; i<Buffer_FUy.N; i++)
	{
		Buffer_FUy.val[i][0] += -FBuffer.val[i][0];
		Buffer_FUy.val[i][1] += -FBuffer.val[i][1];
	}


	// z-Komponente
	d_dz.execute(Fni,FBuffer);
	iFFT(FBuffer, Buffer);


	for(int i=0; i<ni.N; ++i)
	{
		if( (ni.val[i]<=-1.) && (Buffer.val[i]<=1.e-6) )
			Buffer.val[i] = 0.;
		Buffer.val[i] = Buffer.val[i]/(theta*(1.+ni.val[i]));
	}

	FFT(Buffer, FBuffer);

	for(int i=0; i<Buffer_FUz.N; i++)
	{
		Buffer_FUz.val[i][0] += -FBuffer.val[i][0];
		Buffer_FUz.val[i][1] += -FBuffer.val[i][1];
	}


   #ifdef _MY_VERBOSE
	log <<"done";
   #endif
	return;
}
