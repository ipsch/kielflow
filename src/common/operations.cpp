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


void dealiasing_36er(field_imag &DEST)
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
				double kabs = sqrt(pow(kx/kx_max,2.)+pow(ky/ky_max,2.)+pow(kz/kz_max,2.));
				double filter_val = exp(-36.*pow(kabs,36.));
				DEST.val[ijk][0] *= filter_val;
				DEST.val[ijk][1] *= filter_val;

			}
		}
	}

	return;
}

void dealiasing_23rd(field_imag &DEST)
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
				double kabs = sqrt(pow(kx/kx_max,2.)+pow(ky/ky_max,2.)+pow(kz/kz_max,2.));
				double filter_val =  kabs < 2./3. ? 1 : 0;
				DEST.val[ijk][0] *= filter_val;
				DEST.val[ijk][1] *= filter_val;

			}
		}
	}

	return;
}


void dealiasing_theta(field_imag &DEST, const double &x)
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

