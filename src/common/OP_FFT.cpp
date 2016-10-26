#include "OP_FFT.hpp"

#define MY_FFTW_MODE FFTW_PATIENT
#ifndef MY_FFTW_MODE
#define MY_FFTW_MODE FFTW_ESTIMATE
#endif

OP_FFT::OP_FFT(const grid &domain) :
	N(domain.Nx*domain.Ny*domain.Nz)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
    logger my_log("OP_FFT::OP_FFT(const grid &domain)");
	my_log << "start";
   #endif
	//fftw_import_wisdom_from_filename("./data/wisdom.fftw");

	input = (double*)  fftw_malloc(sizeof(double) * domain.Nx*domain.Ny*domain.Nz);
	output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * domain.Nx*domain.Ny*(domain.Nz/2+1));
	//my_plan = fftw_plan_dft_r2c_3d(domain.Nx, domain.Ny, domain.Nz, input, output,
	//		FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT); // ToDo : change if program works
	my_plan = fftw_plan_dft_r2c_3d(domain.Nx, domain.Ny, domain.Nz, input, output, MY_FFTW_MODE);

	// optimization level
	// FFTW_ESTIMATE
	// FFTW_MEASURE
	// FFTW_PATIENT
	// FFTW_EXHAUSTIVE

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
}


OP_FFT::~OP_FFT()
{
	//fftw_export_wisdom_to_filename("./data/wisdom.fftw");
	fftw_destroy_plan(my_plan);
	fftw_free(input);
	fftw_free(output);
}


void OP_FFT::operator()(const field_real &in, field_imag &out)
{
	for(int i=0; i<in.N; ++i)
		input[i] = in.val[i];

	fftw_execute(my_plan);

	for (int i=0; i < out.N; ++i)
	{
		out.val[i][0] = output[i][0];
		out.val[i][1] = output[i][1];
	}

	return;
}
