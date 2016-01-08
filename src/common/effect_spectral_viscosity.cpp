#include "effect_spectral_viscosity.hpp"






effect_spectral_viscosity::effect_spectral_viscosity(const double &fraction, const double &L, const int &N)
{
	const double pi = acos(-1.);

	double k_max = (pi*N/L);
	cutoff = fraction*k_max;
	cutoff = pow(cutoff,2.);
	eps = 0.25/(fraction*double(N));
}

inline double effect_spectral_viscosity::filter(const double &k2) const
{
	return (k2<cutoff) ? 0 : 1;
}

void effect_spectral_viscosity::execute(const field_imag &in, field_imag &out) const
{

	field_imag FTemp(*in.my_grid);
	double kx, ky, kz;

	double retrn;

	for(int i=0; i<in.Nx;++i)
	{
		kx = in.my_grid->x_axis->val_at(i);
		for(int j=0; j<in.Ny; ++j)
		{
			ky = in.my_grid->y_axis->val_at(j);
			for(int k=0; k<in.Nz; ++k)
			{
				kz = in.my_grid->z_axis->val_at(k);
				int index= in.my_grid->index_at(i,j,k);
				double k2 = kx*kx + ky*ky * kz*kz;
				out.val[index][0] += -eps*k2*this->filter(i)*in.val[index][0];
				out.val[index][1] += -eps*k2*this->filter(i)*in.val[index][1];
			}
		}
	}


	return;

}
