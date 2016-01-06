#include "CFL.hpp"


// ToDo : I've copied this part (from solver_poisson_jacobi_nlin), make it reuseable
// from some where else
double get_delta(const axis * const A, const int &i)
{
	if( (i==A->N-1) && (i==-1) )
		return A->val_at(1) - A->val_at(0);

	return A->val_at(i+1) - A->val_at(i);
}


// ToDo : Implement if needed
/*
double cfl_single(const double &dt, const field_real &XX)
{
	for(int i=0; i<XX.Nx; ++i)
	{
		double delta_x;
	}
	return 0;
}
*/



double cfl_all(const double &dt, const field_imag &FUx, const field_imag &FUy, const field_imag &FUz)
{
	field_real Ux(*FUx.my_grid);
	field_real Uy(*FUy.my_grid);
	field_real Uz(*FUz.my_grid);

	iFFT(FUx,Ux);
	iFFT(FUy,Uy);
	iFFT(FUz,Uz);

	return cfl_all(dt,Ux,Uy,Uz);
}


double cfl_all(const double &dt, const field_real &Ux, const field_real &Uy, const field_real &Uz)
{
	double max_overall;

	double CFLx = 0.;
	double CFLy = 0.;
	double CFLz = 0.;
	double CFLxory = 0.;
	double CFLmax = 0.;

	for(int i=0; i<Ux.Nx; ++i)
	{
		double delta_x = min<double>(get_delta(Ux.my_grid->x_axis, i), get_delta(Ux.my_grid->x_axis, i+1));
		for(int j=0; j<Ux.Ny; ++j)
		{
			double delta_y = min<double>(get_delta(Uy.my_grid->y_axis,j), get_delta(Uy.my_grid->y_axis,j+1) );
			for(int k=0; k<Ux.Nz; ++k)
			{
				double delta_z = min<double>(get_delta(Uz.my_grid->z_axis,k), get_delta(Uz.my_grid->z_axis,k+1) );

				int index = Ux.my_grid->index_at(i,j,k);

				CFLx = max<double>(CFLx, dt*fabs(Ux.val[index])/delta_x );
				CFLy = max<double>(CFLy, dt*fabs(Uy.val[index])/delta_y );
				CFLz = max<double>(CFLz, dt*fabs(Uz.val[index])/delta_z );

				CFLmax = max<double>(CFLx, max<double>(CFLy,CFLz));
			}
		}
	}

	std::cout << "CFLmax= " << CFLmax << std::endl;
	std::cout << "dt= " << dt << std::endl;

	if(.4*dt/CFLmax>= 0.01)
		return 0.01;

	return .4*dt/CFLmax;

}
