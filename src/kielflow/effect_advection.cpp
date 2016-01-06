#include "effect_advection.hpp"


effect_advection::effect_advection()
{

}


void effect_advection::execute(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz, const field_imag FXX,
		field_imag &return_Buffer)
{
   #ifdef _MY_VERBOSE
	logger log("effect_advection");
	log << "execute(...)";
   #endif

	// input  F[Ux], F[Uy], F[Uz] and F[XX]  (XX denotes a wildcard that can either be Ux, Uy or Uz)
	// output F[ Ux*(d/dx)XX + Uy*(d/dy)XX + Uz*(d/dz)XX ]

	// step1: berechne hier F^{-1}[F[Ux]], F^{-1}[F[Uy]] und F^{-1}[F[Uz]]
	// dies ist identisch mit
	// Ux, Uy und Uz
	field_real Ux(*FUx.my_grid);
	field_real Uy(*FUy.my_grid);
	field_real Uz(*FUz.my_grid);


	static OP_partial_derivative d_dx(FUx,e_x);
	static OP_partial_derivative d_dy(FUx,e_y);
	static OP_partial_derivative d_dz(FUx,e_z);

	iFFT(FUx,Ux);
	iFFT(FUy,Uy);
	iFFT(FUz,Uz);


	// step2: berechne hier F^{-1}[îkx*F[Ux]], F^{-1}[îky*F[Ux]] und F^{-1}[îkz*F[Ux]]
	// dies ist identisch mit
	// (d/dx) Ux, (d/dy) Ux und (d/dz) Ux

	field_imag FBuffer(*FXX.my_grid);


	d_dx.execute(FXX,FBuffer);
	field_real dx_XX(*FBuffer.my_grid);
	iFFT(FBuffer,dx_XX);


	d_dy.execute(FXX,FBuffer);
	field_real dy_XX(*FBuffer.my_grid);
	iFFT(FBuffer,dy_XX);

	d_dz.execute(FXX,FBuffer);
	field_real dz_XX(*FBuffer.my_grid);
	iFFT(FBuffer,dz_XX);

	// step3: berechne hier F[ Ux*(d/dx)Ux + Uy*(d/dy)Ux + Uz*(d/dz)Ux ]

	field_real Buffer(*FXX.my_grid);


	for(int i=0; i<Ux.N; ++i)
		Buffer.val[i] = Ux.val[i]*dx_XX.val[i] + Uy.val[i]*dy_XX.val[i] + Uz.val[i]*dz_XX.val[i];


	FFT(Buffer, FBuffer);


	// Add result to RHS-Buffer
	for(int i=0; i<FUx.N; ++i)
	{
		return_Buffer.val[i][0] += -FBuffer.val[i][0];
		return_Buffer.val[i][1] += -FBuffer.val[i][1];
	}

   #ifdef _MY_VERBOSE
	log << "done";
   #endif

	return;
}
