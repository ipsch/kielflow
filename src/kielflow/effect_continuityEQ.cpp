#include "effect_continuityEQ.hpp"

void effect_continuityEQ::execute(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni,
		field_imag &Buffer_Fni)
     // - grad(n*u)

// This function evaluates the continuity equation. First a copy of Ux,Uy,Uz and ni is created in
// real-space. In real-space the products, Ux*ni, Uy*ni, Uz*ni are evaluated.
// There after the products are transformed back to fourier-space and are multiplied by ik_x btw. ik_y, ik_z (which
// accounts for the derivative in space). At last the products are added to the output-buffer.
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_continuityEQ(...)";
   #endif

	static OP_partial_derivative d_dx(FUx,e_x);
	static OP_partial_derivative d_dy(FUx,e_y);
	static OP_partial_derivative d_dz(FUx,e_z);

	// Temporäre größen im Ortsraum (iFFT impliziert)
	field_real ni(*Fni.my_grid);
	iFFT(Fni,ni);
	field_real Ux(*FUx.my_grid);
	iFFT(FUx,Ux);
	field_real Uy(*FUy.my_grid);
	iFFT(FUy,Uy);
	field_real Uz(*FUz.my_grid);
	iFFT(FUz,Uz);

	// evaluate non-linearity in realspace FUx \cdot ni (and so on)
	for(int i=0; i<Ux.N; ++i)
	{
		Ux.val[i] = Ux.val[i]*ni.val[i];
		Uy.val[i] = Uy.val[i]*ni.val[i];
		Uz.val[i] = Uz.val[i]*ni.val[i];
	}

	// backtransformation
	field_imag FBuffer_x(*Ux.my_grid);
	FFT(Ux,FBuffer_x);
	field_imag FBuffer_y(*Uy.my_grid);
	FFT(Uy,FBuffer_y);
	field_imag FBuffer_z(*Uz.my_grid);
	FFT(Uz,FBuffer_z);

	// Ableitung bilden
	d_dx.execute(FBuffer_x,FBuffer_x);
	d_dy.execute(FBuffer_y,FBuffer_y);
	d_dz.execute(FBuffer_z,FBuffer_z);


	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] += -(FBuffer_x.val[i][0] + FBuffer_y.val[i][0] + FBuffer_z.val[i][0]);
		Buffer_Fni.val[i][1] += -(FBuffer_x.val[i][1] + FBuffer_y.val[i][1] + FBuffer_z.val[i][1]);
	}


   #ifdef _MY_VERBOSE
	log << "done";
   #endif
	return;
}
