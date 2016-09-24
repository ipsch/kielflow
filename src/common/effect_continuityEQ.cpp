#include "effect_continuityEQ.hpp"

void effect_continuityEQ::execute(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni,
		field_imag &Buffer_Fni)
     // - div(n*\mathbf{u})

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


	static field_imag FBuffer(FUx.my_grid);
	static field_real Buffer_n(FUx.my_grid);
	static field_real Buffer_u(FUx.my_grid);


	d_dx(FUx,FBuffer);
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] = -FBuffer.val[i][0];
		Buffer_Fni.val[i][1] = -FBuffer.val[i][1];
	}

	d_dy(FUy,FBuffer);
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] -= FBuffer.val[i][0];
		Buffer_Fni.val[i][1] -= FBuffer.val[i][1];
	}

	d_dz(FUz,FBuffer);
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] -= FBuffer.val[i][0];
		Buffer_Fni.val[i][1] -= FBuffer.val[i][1];
	}


	d_dx(Fni,FBuffer);
	iFFT(FBuffer,Buffer_n);
	iFFT(FUx,Buffer_u);
	for(int i=0; i<Buffer_u.N; ++i)
		Buffer_n.val[i] = Buffer_n.val[i]*Buffer_u.val[i];
	FFT(Buffer_n, FBuffer);
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] -= FBuffer.val[i][0];
		Buffer_Fni.val[i][1] -= FBuffer.val[i][1];
	}

	d_dy(Fni,FBuffer);
	iFFT(FBuffer,Buffer_n);
	iFFT(FUy,Buffer_u);
	for(int i=0; i<Buffer_u.N; ++i)
		Buffer_n.val[i] = Buffer_n.val[i]*Buffer_u.val[i];
	FFT(Buffer_n, FBuffer);
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] -= FBuffer.val[i][0];
		Buffer_Fni.val[i][1] -= FBuffer.val[i][1];
	}

	d_dz(Fni,FBuffer);
	iFFT(FBuffer,Buffer_n);
	iFFT(FUz,Buffer_u);
	for(int i=0; i<Buffer_u.N; ++i)
		Buffer_n.val[i] = Buffer_n.val[i]*Buffer_u.val[i];
	FFT(Buffer_n, FBuffer);
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] -= FBuffer.val[i][0];
		Buffer_Fni.val[i][1] -= FBuffer.val[i][1];
	}






/*
	// n in logscale alternative formulation
	static OP_partial_derivative d_dx(FUx,e_x);
	static OP_partial_derivative d_dy(FUy,e_y);
	static OP_partial_derivative d_dz(FUz,e_z);

	field_imag FxBuffer(*FUx.my_grid);
	field_imag FyBuffer(*FUx.my_grid);
	field_imag FBuffer(*FUx.my_grid);
	field_real Buffer_n(*FUx.my_grid);
	field_real Buffer_u(*FUx.my_grid);


	iFFT(Fni, Buffer_n);

	// drei parallele stränge
	iFFT(FUx, Buffer_u);
	// Ableitung bilden
	for(int i=0; i<Buffer_u.N; ++i)
		Buffer_u.val[i] = Buffer_u.val[i]*exp(Buffer_n.val[i]);
	FFT(Buffer_u, FBuffer);
	d_dx.execute(FBuffer,FxBuffer);

	iFFT(FUy, Buffer_u);
	// Ableitung bilden
	for(int i=0; i<Buffer_u.N; ++i)
		Buffer_u.val[i] = Buffer_u.val[i]*exp(Buffer_n.val[i]);
	FFT(Buffer_u, FBuffer);
	d_dy.execute(FBuffer,FyBuffer);

	iFFT(FUz, Buffer_u);
	// Ableitung bilden
	for(int i=0; i<Buffer_u.N; ++i)
		Buffer_u.val[i] = Buffer_u.val[i]*exp(Buffer_n.val[i]);
	FFT(Buffer_u, FBuffer);
	d_dz.execute(FBuffer,FBuffer);

	// vereinigung
	for(int i=0; i<FBuffer.N; ++i)
	{
		FBuffer.val[i][0] += FxBuffer.val[i][0] + FyBuffer.val[i][0];
		FBuffer.val[i][1] += FxBuffer.val[i][1] + FyBuffer.val[i][1];
	}

	// letzter schritt teilen durch exp fkt
	iFFT(FBuffer,Buffer_u);
    for(int i=0; i<Buffer_u.N; ++i)
    	Buffer_u.val[i] = exp(-Buffer_n.val[i])*Buffer_u.val[i];
    FFT(Buffer_u,FBuffer);

    // Übertrag
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] = -FBuffer.val[i][0];
		Buffer_Fni.val[i][1] = -FBuffer.val[i][1];
	}
*/



   #if defined (_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	log << "done";
   #endif
	return;
}
