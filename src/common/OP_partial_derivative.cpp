#include "OP_partial_derivative.hpp"


OP_partial_derivative::OP_partial_derivative(const field_real &field_caller, const direction &e_i) :
		OP_partial_derivative(field_caller.my_grid, e_i)
{

}

OP_partial_derivative::OP_partial_derivative(const grid &domain, const direction &e_i) :
		my_FFT(domain), my_iFFT(domain), Buffer(domain),
		my_dealiasing(domain, [] (const double &k) {return k<(2./3.) ? 1. : 0.;})
{
	switch(e_i)
	{
	case e_x :
		my_axis = domain.x_axis;
		my_X = &my_i;
		break;

	case e_y :
		my_axis = domain.y_axis;
		my_X = &my_j;
		break;

	case e_z :
		my_axis = domain.z_axis;
		my_X = &my_k;
		break;
	} // END OF SWITCH;


	switch(my_axis->type_id)
	{
	case 1:
		IsLinear = true;
		break;
	case 2:
		IsLinear = false;
		break;
	case 3:
		IsLinear = false;
		break;
	} // END OF SWITCH

	double N_re = domain.Nx*domain.Ny*(domain.Nz);
	double N_im = domain.Nx*domain.Ny*(domain.Nz/2+1);
	k_val = (double*) fftw_malloc(sizeof(double) * N_im);
	dS_val = (double*) fftw_malloc(sizeof(double) * N_re);

	// set up wavevector (dependent of direction of derivative)
	for(my_i=0; my_i<domain.Nx; ++my_i)
		for(my_j=0; my_j<domain.Ny; ++my_j)
			for(my_k=0; my_k<domain.Nz/2+1; ++my_k)
			{
				int ijk = (my_k + (domain.Nz/2+1)*(my_j + my_i*domain.Ny));
				k_val[ijk] = my_axis->k_val_at(*my_X);
			}

	// set up derivative of mapping function (needed in evaluation of
	// partial derivatives on non uniform grids
	for(my_i=0; my_i<domain.Nx; ++my_i)
		for(my_j=0; my_j<domain.Ny; ++my_j)
			for(my_k=0; my_k<domain.Nz; ++my_k)
			{
				int ijk = (my_k + (domain.Nz)*(my_j + my_i*domain.Ny));
				dS_val[ijk] = 1./my_axis->dS(*my_X);
			}

	// dS_val needs to be dealiased since partial derivatives are
	// nonlinear operations on non uniform grids :
	for(int ijk=0; ijk<Buffer.N; ++ijk)
		Buffer.val[ijk] = dS_val[ijk];
	field_imag FBuffer(domain);
	my_FFT(Buffer,FBuffer);
	my_dealiasing(FBuffer);
	my_iFFT(FBuffer,Buffer);

	for(int ijk=0; ijk<Buffer.N; ++ijk)
		dS_val[ijk] = Buffer.val[ijk];

	my_axis = 0L;

}

OP_partial_derivative::OP_partial_derivative(const field_imag &field_caller, const direction &e_i) :
		OP_partial_derivative(field_caller.my_grid,e_i)
{

}


OP_partial_derivative::~OP_partial_derivative()
{
	fftw_free(k_val);
	fftw_free(dS_val);
}






void OP_partial_derivative::operator()(const field_imag &in, field_imag &out)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("OP_partial_derivative");
	my_log << "nonlinear_call(field_imag &in, field_imag &out)";
   #endif

	for(int ijk=0; ijk<in.N; ++ijk)
	{
		double tmp = in.val[ijk][0]; // needed in case of in==out
		out.val[ijk][0] = -k_val[ijk]*in.val[ijk][1];
		out.val[ijk][1] =  k_val[ijk]*tmp;
	}

	if(IsLinear)
		return; // done if grid is uniform

	my_dealiasing(out);
	field_real tmp(in.my_grid);
	my_iFFT(out,tmp);

	// multiply by derivative of mapping function in real-space
	for(int ijk=0; ijk<tmp.N; ++ijk)
		tmp.val[ijk] = tmp.val[ijk]*dS_val[ijk];

	my_FFT(tmp,out);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}



