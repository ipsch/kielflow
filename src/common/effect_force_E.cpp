#include "effect_force_E.hpp"

effect_force_E::effect_force_E()
{

}

void effect_force_E::execute(const field_imag &FPhi,
		field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz)
// Evaluate the stato-electrical force upon the velocity field of the ions.
// F = -(q/m)*grad(Phi)
// The factor q/m is already included by equations due to dimensionless quantities
// evaluation of grad(Phi) is done in fourier-space by multiplying fourier coeffcients by (ik) (pay attention to sign an re<->im beeing crossed)
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_force_E(...)";
   #endif


	static OP_partial_derivative d_dx(FPhi,e_x);
	static OP_partial_derivative d_dy(FPhi,e_y);
	static OP_partial_derivative d_dz(FPhi,e_z);

	field_imag FBuffer(FPhi.my_grid);



	d_dx(FPhi,FBuffer);
	for(int i=0; i<FPhi.N; ++i)
	{
		Buffer_FUx.val[i][0] -= FBuffer.val[i][0];
		Buffer_FUx.val[i][1] -= FBuffer.val[i][1];
	}

	d_dy(FPhi,FBuffer);
	for(int i=0; i<FPhi.N; ++i)
	{
		Buffer_FUy.val[i][0] -= FBuffer.val[i][0];
		Buffer_FUy.val[i][1] -= FBuffer.val[i][1];
	}

	d_dz(FPhi,FBuffer);
	for(int i=0; i<FPhi.N; ++i)
	{
		Buffer_FUz.val[i][0] -= FBuffer.val[i][0];
		Buffer_FUz.val[i][1] -= FBuffer.val[i][1];
	}


   #ifdef _MY_VERBOSE
	log <<"done";
   #endif
	return;
}
