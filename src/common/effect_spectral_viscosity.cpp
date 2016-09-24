#include "effect_spectral_viscosity.hpp"


effect_spectral_viscosity::effect_spectral_viscosity(const double &fraction, const grid &Omega)
{
	const double pi = acos(-1.);
	kx_max = (pi*Omega.x_axis->N/Omega.x_axis->L);
	ky_max = (pi*Omega.y_axis->N/Omega.y_axis->L);
	kz_max = (pi*Omega.z_axis->N/Omega.z_axis->L);
	cutoff = fraction;
	//cutoff = cutoff*cutoff;
	//eps = 0.25/(fraction*double(N));
	eps = 0.01;

}


inline double effect_spectral_viscosity::filter(const double &kx, const double &ky, const double &kz) const
{
	double k2 = sqrt(pow(kx/kx_max,2.) + pow(ky/ky_max,2.) + pow(kz/kz_max,2.));
	//return (k2<cutoff) ? 0. : x;
	double x = 1. - exp(-36.*pow(sqrt( k2 ),15.));
	return x;
}


void effect_spectral_viscosity::execute(const field_imag &in, field_imag &out) const
{
	double kx, ky, kz;
	for(int i=0; i<in.Nx;++i)
	{
		kx = in.my_grid.x_axis->val_at(i);
		for(int j=0; j<in.Ny; ++j)
		{
			ky = in.my_grid.y_axis->val_at(j);
			for(int k=0; k<in.Nz; ++k)
			{
				kz = in.my_grid.z_axis->val_at(k);
				int index= in.index(i,j,k);
				double k2 = kx*kx + ky*ky + kz*kz;

				out.val[index][0] += -eps*k2*this->filter(kx,ky,kz)*in.val[index][0];
				out.val[index][1] += -eps*k2*this->filter(ky,ky,kz)*in.val[index][1];
			}
		}
	}
	return;
}
