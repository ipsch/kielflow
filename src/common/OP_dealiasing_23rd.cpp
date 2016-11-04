#include "OP_dealiasing_23rd.hpp"



OP_dealiasing_23rd::OP_dealiasing_23rd(const grid &domain) :
N(domain.Nx*domain.Ny*(domain.Nz/2+1))
{
	std::cout << "OP_dea: " << N << std::endl;
	filter = new double[N];

	double kx_max=domain.x_axis->k_val_at(domain.Nx/2);
	double ky_max=domain.y_axis->k_val_at(domain.Ny/2);
	double kz_max=domain.z_axis->k_val_at(domain.Nz/2);

	// setup dealiasing-filter 2/3-rule
	for(int i=0; i<domain.Nx; i++)
		for(int j=0; j<domain.Ny; j++)
			for (int k=0; k<(domain.Nz/2)+1; k++)
			{
				int ijk = k + (domain.Nz/2+1)*(j + i*domain.Ny);
				double kx = domain.x_axis->k_val_at(i);
				double ky = domain.y_axis->k_val_at(j);
				double kz = domain.z_axis->k_val_at(k);
				double kabs = sqrt(pow(kx/kx_max,2.)+pow(ky/ky_max,2.)+pow(kz/kz_max,2.));
				filter[ijk] = kabs<(2./3.) ? 1 : 0;
			}

	return;
}

void OP_dealiasing_23rd::operator() (field_imag &target) const
{
	for(int ijk=0; ijk<N; ijk++)
	{
		target.val[ijk][0]*=filter[ijk];
		target.val[ijk][1]*=filter[ijk];
	}
	return;
}

OP_dealiasing_23rd::~OP_dealiasing_23rd()
{
	delete[] filter;
}
