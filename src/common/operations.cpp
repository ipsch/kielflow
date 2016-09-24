#include "operations.hpp"

//#define NORMALIZATION_SYMMETRIC_

/*
void handle_dealiasing(field_imag &DEST)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE)
	logger log("operation");
	log << "handle_dealiasing(field_imag &DEST)";
   #endif

	double kx_max = DEST.my_grid.x_axis->val_at(DEST.Nx/2);
	double ky_max = DEST.my_grid.y_axis->val_at(DEST.Ny/2);
	double kz_max = DEST.my_grid.z_axis->val_at(DEST.Nz/2);

	for(int i=0;i<DEST.Nx; ++i)
	{
		double kx = DEST.my_grid.x_axis->val_at(i);
		for(int j=0;j<DEST.Ny; ++j)
		{
			double ky = DEST.my_grid.y_axis->val_at(j);
			for (int k=0; k< DEST.Nz; ++k)
			{
				double kz = DEST.my_grid.z_axis->val_at(k);

				// ToDo : Filter ist zur hart eingestellt.
				if(( kx*kx/(kx_max*kx_max) + ky*ky/(ky_max*ky_max) + kz*kz/(kz_max*kz_max) ) > 2./9.)
				{
					int .index = DEST.index(i,j,k);
					DEST.val[index][0] = 0.;
					DEST.val[index][1] = 0.;

				}
			}
		}
	}

	return;
}
*/

/*
double filter_fkt(const double &kx, const double &kx_max,
		          const double &ky, const double &ky_max,
		          const double &kz, const double &kz_max)
{
	const double C = 36.; // default = 36.

	return exp(-C*pow(sqrt( (kx*kx)/(kx_max*kx_max) + (ky*ky)/(ky_max*ky_max) + (kz*kz)/(kz_max*kz_max) ),C));

	return 0;
}
*/
double filter_fkt(const int &i, const int &i_max,
		          const int &j, const int &j_max,
		          const int &k, const int &k_max)
{

	double kx_sqare = (double(i)*double(i))/(double(i_max)*double(i_max));
	double ky_sqare = (double(j)*double(j))/(double(j_max)*double(j_max));
	double kz_sqare = (double(k)*double(k))/(double(k_max)*double(k_max));

	double kr = kx_sqare + ky_sqare + kz_sqare;

	return (kr<(1./9.)) ? 1. : 0. ;
}


void dealaising_36er(field_imag &DEST)
{
	double kx_max = DEST.my_grid.x_axis->k_val_at(DEST.Nx/2);
	double ky_max = DEST.my_grid.y_axis->k_val_at(DEST.Ny/2);
	double kz_max = DEST.my_grid.z_axis->k_val_at(DEST.Nz/2);

	for(int i=0;i<DEST.Nx; ++i)
	{
		double kx = DEST.my_grid.x_axis->k_val_at(i);
		for(int j=0;j<DEST.Ny; ++j)
		{
			double ky = DEST.my_grid.y_axis->k_val_at(j);
			for (int k=0; k< DEST.Nz; ++k)
			{
				double kz = DEST.my_grid.z_axis->k_val_at(k);

				int ijk = DEST.index(i,j,k);
				double filter_val = exp(-36.*pow(sqrt( pow(kx/kx_max,2.) + pow(ky/ky_max,2.) + pow(kz/kz_max,2.) ),36.));
				DEST.val[ijk][0] *= filter_val;
				DEST.val[ijk][1] *= filter_val;

			}
		}
	}

	return;
}

void dealaising_23rd(field_imag &DEST)
{
	double kx_max = DEST.my_grid.x_axis->k_val_at(DEST.Nx/2);
	double ky_max = DEST.my_grid.y_axis->k_val_at(DEST.Ny/2);
	double kz_max = DEST.my_grid.z_axis->k_val_at(DEST.Nz/2);

	for(int i=0;i<DEST.Nx; ++i)
	{
		double kx = DEST.my_grid.x_axis->k_val_at(i);
		for(int j=0;j<DEST.Ny; ++j)
		{
			double ky = DEST.my_grid.y_axis->k_val_at(j);
			for (int k=0; k< DEST.Nz; ++k)
			{
				double kz = DEST.my_grid.z_axis->k_val_at(k);

				int ijk = DEST.index(i,j,k);
				double filter_val = sqrt( pow(kx/kx_max,2.) + pow(ky/ky_max,2.) + pow(kz/kz_max,2.) ) < 2./3. ? 1 : 0;
				DEST.val[ijk][0] *= filter_val;
				DEST.val[ijk][1] *= filter_val;

			}
		}
	}

	return;
}


void dealaising_theta(field_imag &DEST, const double &x)
{
	double kx_max = DEST.my_grid.x_axis->k_val_at(DEST.Nx/2);
	double ky_max = DEST.my_grid.y_axis->k_val_at(DEST.Ny/2);
	double kz_max = DEST.my_grid.z_axis->k_val_at(DEST.Nz/2);

	for(int i=0;i<DEST.Nx; ++i)
	{
		double kx = DEST.my_grid.x_axis->k_val_at(i);
		for(int j=0;j<DEST.Ny; ++j)
		{
			double ky = DEST.my_grid.y_axis->k_val_at(j);
			for (int k=0; k< DEST.Nz; ++k)
			{
				double kz = DEST.my_grid.z_axis->k_val_at(k);

				int ijk = DEST.index(i,j,k);
				double filter_val = sqrt( pow(kx/kx_max,2.) + pow(ky/ky_max,2.) + pow(kz/kz_max,2.) ) < x ? 1 : 0;
				DEST.val[ijk][0] *= filter_val;
				DEST.val[ijk][1] *= filter_val;

			}
		}
	}

	return;
}

void handle_dealiasing(field_imag &DEST)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	logger my_log("operation");
	my_log << "handle_dealiasing(field_imag &DEST)";
   #endif


	double kx_max = DEST.my_grid.x_axis->k_val_at(DEST.Nx/2);
	double ky_max = DEST.my_grid.y_axis->k_val_at(DEST.Ny/2);
	double kz_max = DEST.my_grid.z_axis->k_val_at(DEST.Nz/2);


	for(int i=0;i<DEST.Nx; ++i)
	{
		double kx = DEST.my_grid.x_axis->k_val_at(i);
		for(int j=0;j<DEST.Ny; ++j)
		{
			double ky = DEST.my_grid.y_axis->k_val_at(j);
			for (int k=0; k< DEST.Nz; ++k)
			{
				double kz = DEST.my_grid.z_axis->k_val_at(k);

				int ijk = DEST.index(i,j,k);
				DEST.val[ijk][0] *= filter_fkt(kx, kx_max, ky, ky_max, kz, kz_max);
				DEST.val[ijk][1] *= filter_fkt(kx, kx_max, ky, ky_max, kz, kz_max);

			}
		}
	}

	return;
}


// ToDO : Das hier funktioniert irgendwie noch nicht richtig :(
void handle_dealiasing_BACKUP(field_imag &DEST)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("handle_dealiasing_BACKUP(field_imag &DEST)");
	my_log << "start";
   #endif



	for(int i=0;i<DEST.Nx; ++i)
	{
		for(int j=0;j<DEST.Ny; ++j)
		{
			for (int k=0; k< DEST.Nz; ++k)
			{
				int ijk = DEST.index(i,j,k);
				DEST.val[ijk][0] *= filter_fkt(i, DEST.Nx, j, DEST.Ny, k, DEST.Nz);
				DEST.val[ijk][1] *= filter_fkt(i, DEST.Nx, j, DEST.Ny, k, DEST.Nz);
			}
		}
	}

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}


void FFT(const field_real &in, field_imag &out)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("FFT(const field_real &in, field_imag &out)");
	my_log << "start";
   #endif
	field_real tmp = field_real(in);

	fftw_plan FFT = fftw_plan_dft_r2c_3d(in.Nx, in.Ny, in.Nz, tmp.val, out.val, FFTW_ESTIMATE);
	fftw_execute(FFT);
	fftw_destroy_plan(FFT);

   #if defined (NORMALIZATION_SYMMETRIC_)
	for (int i = 0; i < out.N; ++i)
	{
		out.val[i][0] *= 1.0 / sqrt(in.N);
		out.val[i][1] *= 1.0 / sqrt(in.N);
	}
   #endif

   #if defined(_MY_VERBOSE_MORE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}


void iFFT(const field_imag &in, field_real &out)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("operation");
	my_log << "iFFT(const field_imag &in, field_real &out)";
   #endif

	field_imag tmp = field_imag(in);

	//handle_dealiasing(tmp);

	fftw_plan iFFT = fftw_plan_dft_c2r_3d(out.Nx, out.Ny, out.Nz, tmp.val, out.val, FFTW_ESTIMATE);
	fftw_execute(iFFT);
	fftw_destroy_plan(iFFT);

   #if defined (NORMALIZATION_SYMMETRIC_)
	for (int i=0; i < out.N; ++i)
		out.val[i] *= 1.0 / sqrt(out.N); // normierung
   #else
	for (int i=0; i < out.N; ++i)
		out.val[i] *= 1.0 / out.N; // normierung
   #endif

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}






OP_partial_derivative::OP_partial_derivative(const field_real &field_caller, const direction &e_i) :
		OP_partial_derivative(field_caller.my_grid, e_i)
{

}

OP_partial_derivative::OP_partial_derivative(const grid &domain, const direction &e_i) :
		my_FFT(domain), my_iFFT(domain), Buffer(domain)
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
	dealaising_36er(FBuffer);
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

	dealaising_36er(out);
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







double supremum(const field_real &in)
{
	double sup = 0.;
	for(int i=0; i<in.N; ++i)
		sup = max<double>(sup, fabs(in.val[i]));
	return sup;
}


double supremum(const field_imag &in)
{
	double sup = 0.;
	double val = 0.;
	for(int i=0; i<in.N; ++i)
	{
		val = in.val[i][0]*in.val[i][0] + in.val[i][1]*in.val[i][1];
		sup = max<double>(sup, val);
	}
	return sqrt(sup);
}

