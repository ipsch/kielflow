#ifndef OP_DEALIASING_HPP_
#define OP_DEALIASING_HPP_

#include <iostream>
#include <cmath>
#include "grid.hpp"
#include "field_imag.hpp"

class OP_dealiasing
{
public :
	template<typename Func>
	OP_dealiasing(const grid &domain, Func fkt) :
	N(domain.Nx*domain.Ny*(domain.Nz/2+1))
	{
		filter = new double[N];
		double kx_max=domain.x_axis->k_val_at(domain.Nx/2);
		double ky_max=domain.y_axis->k_val_at(domain.Ny/2);
		double kz_max=domain.z_axis->k_val_at(domain.Nz/2);
		for(int i=0; i<domain.Nx; i++)
			for(int j=0; j<domain.Ny; j++)
				for (int k=0; k<(domain.Nz/2)+1; k++)
				{
					int ijk = k + (domain.Nz/2+1)*(j + i*domain.Ny);
					double kx = domain.x_axis->k_val_at(i);
					double ky = domain.y_axis->k_val_at(j);
					double kz = domain.z_axis->k_val_at(k);
					double kabs = sqrt(pow(kx/kx_max,2.)+pow(ky/ky_max,2.)+pow(kz/kz_max,2.));
					filter[ijk] = fkt(kabs);
				}
		return;
	}
	~OP_dealiasing();
	void operator() (field_imag &target) const;
private :
	const unsigned int N;
	double * filter;
	OP_dealiasing() : N(0) {filter = 0L;}
};


#endif
