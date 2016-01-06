#include "field_interpolation.hpp"


void OP_hto2h(const field_real &in, field_real &out)
// level 2 interpolation Operator
{

	out.resize(2*in.Nx, 2*in.Ny, 2*in.Nz);


	for(int i=0; i<out.Nx; i++)
	{
		for(int j=0; j<out.Ny; j++)
		{
			for(int k=0; k<out.Nz; k++)
			{
				int index_out = out.my_grid->index_at(i,j,k);
				out.val[index_out] = 0.;
			}
		}
	}


	for(int i=0; i<in.Nx; i++)
	{
		for(int j=0; j<in.Ny; j++)
		{
			for(int k=0; k<in.Nz; k++)
			{
				int index_in  = in.my_grid->index_at(  i,  j,  k);
				int index_out = out.my_grid->index_at(2*i,2*j,2*k);
				out.val[index_out] = in.val[index_in];
			}
		}
	}



	for(int i=0; i<out.Nx-1; i+=2)
	{
		for(int j=0; j<out.Ny; j+=2)
		{
			for(int k=1; k<out.Nz; k+=2)
			{
				double z1 = out.my_grid->z_axis->val_at(k-1);
				double z2 = out.my_grid->z_axis->val_at(k+1);
				double z  = out.my_grid->z_axis->val_at(k);
				double gamma_z = (z-z1)/(z2-z1);

				out(i, j, k) = (1.-gamma_z)*out(i,j,k-1) + gamma_z*out(i,j,k+1);
			}
		}
	}


	for(int i=0; i<out.Nx; i+=2)
	{
		for(int j=1; j<out.Ny; j+=2)
		{
			double y1 = out.my_grid->y_axis->val_at(j-1);
			double y2 = out.my_grid->y_axis->val_at(j+1);
			double y  = out.my_grid->y_axis->val_at(j);
			double gamma_y = (y-y1)/(y2-y1);

			for(int k=0; k<out.Nz; k+=1)
			{
				out(i, j, k) = (1.-gamma_y)*out(i,j-1,k) + gamma_y*out(i,j+1,k);

			}
		}
	}



	for(int i=1; i<out.Nx; i+=2)
	{
		double x1 = out.my_grid->x_axis->val_at(i-1);
		double x2 = out.my_grid->x_axis->val_at(i+1);
		double x  = out.my_grid->x_axis->val_at(i);
		double gamma_x = (x-x1)/(x2-x1);

		for(int j=0; j<out.Ny; j+=1)
		{
			for(int k=0; k<out.Nz; k+=1)
			{
				out(i, j, k) = (1.-gamma_x)*out(i-1,j,k) + gamma_x*out(i+1,j,k);
			}
		}
	}


	return;
}



void OP_2htoh_lvl0(const field_real &in, field_real &out)
{
	if( (in.Nx<4) || (in.Ny<4) || (in.Nz<4))
		throw("system too small");

	if( ((in.Nx % 2)!=0) || ((in.Ny % 2)!=0) || ((in.Nz % 2) != 0) )
		throw("Feld Groesse ungerade");

	out.resize(in.Nx/2, in.Ny/2, in.Nz/2);

	for(int i=0; i<out.Nx; i++)
	{
		for(int j=0; j<out.Ny; j++)
		{
			for(int k=0; k<out.Nz; k++)
			{
				int index_in  = in.my_grid->index_at( 2*i, 2*j, 2*k);
				int index_out = out.my_grid->index_at(   i,   j,   k);

				out.val[index_out] = in.val[index_in];
			}
		}
	}


	return;
}

void OP_2htoh_lvl1(const field_real &in, field_real &out)
{
	if( (in.Nx<4) || (in.Ny<4) || (in.Nz<4))
		throw("system too small");

	if( ((in.Nx % 2)!=0) || ((in.Ny % 2)!=0) || ((in.Nz % 2) != 0) )
		throw("Feld Groesse ungerade");

	out.resize(in.Nx/2, in.Ny/2, in.Nz/2);

	field_real tmp(*in.my_grid);





	for(int i=1; i<tmp.Nx-1; i++)
	{
		for(int j=0; j<tmp.Ny; j++)
		{
			for(int k=0; k<tmp.Nz; k++)
			{
				int index_zero = tmp.my_grid->index_at( i, j, k);
				int index_plus = tmp.my_grid->index_at( i, j, k+1);
				int index_minu = tmp.my_grid->index_at( i, j, k-1);
				tmp.val[index_zero] =         in.val[index_zero];
				tmp.val[index_zero] = (1./6.)*in.val[index_plus];
				tmp.val[index_zero] = (1./6.)*in.val[index_minu];
			}
		}
	}


	// Direkt
	for(int i=1; i<tmp.Nx-1; i++)
	{
		for(int j=0; j<tmp.Ny; j++)
		{
			for(int k=0; k<tmp.Nz; k++)
			{
				int index = tmp.my_grid->index_at( i, j, k);
				tmp.val[index] =         in.val[index];
			}
		}
	}


	// x-shift +1
	for(int i=0; i<tmp.Nx-1; i++)
	{
		for(int j=0; j<tmp.Ny; j++)
		{
			for(int k=0; k<tmp.Nz; k++)
			{
				int index       = tmp.my_grid->index_at( i  , j  , k );
				int index_shift = tmp.my_grid->index_at( i+1, j  , k );

				tmp.val[index] += (1./6.)*in.val[index_shift];
			}
		}
	}

	// x-shift -1
	for(int i=1; i<tmp.Nx; i++)
	{
		for(int j=0; j<tmp.Ny; j++)
		{
			for(int k=0; k<tmp.Nz; k++)
			{
				int index       = tmp.my_grid->index_at( i  , j  , k );
				int index_shift = tmp.my_grid->index_at( i-1, j  , k );

				tmp.val[index] += (1./6.)*in.val[index_shift];
			}
		}
	}


	// y-shift +1
	for(int i=0; i<tmp.Nx; i++)
	{
		for(int j=0; j<tmp.Ny-1; j++)
		{
			for(int k=0; k<tmp.Nz; k++)
			{
				int index       = tmp.my_grid->index_at( i  , j   , k );
				int index_shift = tmp.my_grid->index_at( i  , j+1 , k );

				tmp.val[index] += (1./6.)*in.val[index_shift];
			}
		}
	}

	// y-shift -1
	for(int i=0; i<tmp.Nx; i++)
	{
		for(int j=1; j<tmp.Ny; j++)
		{
			for(int k=0; k<tmp.Nz; k++)
			{
				int index       = tmp.my_grid->index_at( i  , j   , k );
				int index_shift = tmp.my_grid->index_at( i  , j-1 , k );

				tmp.val[index] += (1./6.)*in.val[index_shift];
			}
		}
	}


	// z-shift +1
	for(int i=0; i<tmp.Nx; i++)
	{
		for(int j=0; j<tmp.Ny; j++)
		{
			for(int k=0; k<tmp.Nz-1; k++)
			{
				int index       = tmp.my_grid->index_at( i  , j   , k   );
				int index_shift = tmp.my_grid->index_at( i  , j   , k+1 );

				tmp.val[index] += (1./6.)*in.val[index_shift];
			}
		}
	}


	// z-shift -1
	for(int i=0; i<tmp.Nx; i++)
	{
		for(int j=0; j<tmp.Ny; j++)
		{
			for(int k=1; k<tmp.Nz; k++)
			{
				int index       = tmp.my_grid->index_at( i  , j   , k   );
				int index_shift = tmp.my_grid->index_at( i  , j   , k-1 );

				tmp.val[index] += (1./6.)*in.val[index_shift];
			}
		}
	}

	for(int i=0; i<out.Nx; i++)
	{
		for(int j=0; j<out.Ny; j++)
		{
			for(int k=0; k<out.Nz; k++)
			{
				int index_in  = in.my_grid->index_at( 2*i, 2*j, 2*k);
				int index_out = out.my_grid->index_at(   i,   j,   k);

				out.val[index_out] = tmp.val[index_in];
			}
		}
	}


	return;
}


void OP_smoothing(const field_real &in, field_real &out)
{

	try
	{
		if( (in.Nx!=out.Nx) || (in.Ny!=out.Ny) || (in.Nz!=out.Nz))
			throw 2;
	}
	catch(int execption_number)
	{
		std::cout << "An exception occurred. Exception Nr. ";
		std::cout << execption_number << std::endl;

		std::cout << in.Nx << "\t << out.Nx" << std::endl;
		std::cout << in.Ny << "\t << out.Ny" << std::endl;
		std::cout << in.Nz << "\t << out.Nz" << std::endl;
		return;
	}

	field_real tmp(*in.my_grid);

	for(int i=1; i<in.Nx-1; i++)
	{
		for(int j=1; j<in.Ny-1; j++)
		{
			for(int k=1; k<in.Nz-1; k++)
			{
				int index  = tmp.my_grid->index_at( i, j, k);

				tmp.val[index] = 0.2*in(i,j,k);
				tmp.val[index]+= 0.1*in(i+1,j,k);
				tmp.val[index]+= 0.1*in(i-1,j,k);
				tmp.val[index]+= 0.1*in(i,j+1,k);
				tmp.val[index]+= 0.1*in(i,j-1,k);
				tmp.val[index]+= 0.1*in(i,j,k+1);
				tmp.val[index]+= 0.1*in(i,j,k-1);
			}
		}
	}

	for(int i=0; i<in.N; ++i)
		out.val[i] = tmp.val[i];

	return;
}
