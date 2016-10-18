#include "OP_iFFT.hpp"


OP_iFFT::OP_iFFT(const grid &domain) :
	N(domain.Nx*domain.Ny*domain.Nz)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("OP_iFFT::OP_iFFT(const grid &domain)");
	my_log << "start";
   #endif
	//fftw_import_wisdom_from_filename("./data/wisdom.fftw");

	input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * domain.Nx*domain.Ny*(domain.Nz/2+1));
	output = (double*)  fftw_malloc(sizeof(double) * domain.Nx*domain.Ny*domain.Nz);
	//my_plan = fftw_plan_dft_c2r_3d(domain.Nx, domain.Ny, domain.Nz, input, output,
	//		FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT); // ToDo : change if program works
	my_plan = fftw_plan_dft_c2r_3d(domain.Nx, domain.Ny, domain.Nz, input, output, FFTW_ESTIMATE);

	// optimization level
	// FFTW_ESTIMATE
	// FFTW_MEASURE
	// FFTW_PATIENT
	// FFTW_EXHAUSTIVE

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
}


OP_iFFT::~OP_iFFT()
{
	//fftw_export_wisdom_to_filename("./data/wisdom.fftw");
	fftw_destroy_plan(my_plan);
	fftw_free(input);
	fftw_free(output);
}


void OP_iFFT::operator()(const field_imag &in, field_real &out)
{
	for(int i=0; i<in.N; ++i)
	{
		input[i][0] = in.val[i][0];
		input[i][1] = in.val[i][1];
	}

	fftw_execute(my_plan);

	double Norm = 1./((double) N);
	for (int i=0; i<out.N; ++i)
		out.val[i] = output[i]*Norm;

	return;
}
