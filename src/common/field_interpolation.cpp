#include "field_interpolation.hpp"

#define _MY_VERBOSE_TEDIOUS

void OP_xto2x(const field_real &in, field_real &out)
{
	out.resize(2*in.Nx, in.Ny, in.Nz);

	// fill every second plain in x-direction with already
	// existing data from "in"
	for(int i=0; i<in.Nx; i++)
	{
		for(int j=0; j<in.Ny; j++)
		{
			for(int k=0; k<in.Nz; k++)
			{
				int index_in  = in.my_grid->index_at(  i,  j,  k);
				int index_out = out.my_grid->index_at(2*i,j,k);
				out.val[index_out] = in.val[index_in];
			}
		}
	}

	// fill every other plain (the ones missing from the first step
	// with data interpolated from the two neighbouring plains
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
				out(i, j, k) = (1.-gamma_x)*out.val_at(i-1,j,k) + gamma_x*out.val_at(i+1,j,k);
			}
		}
	}
	return;
}

void OP_yto2y(const field_real &in, field_real &out)
{
	// ToDo : implement OP_yto2y
	throw("OP_yto2y not implemented yet");
	return;
}

void OP_zto2z(const field_real &in, field_real &out)
{
	// ToDo : implement OP_zto2z
	throw("OP_yto2y not implemented yet");
	return;
}

void OP_hto2h(const field_real &in, field_real &out)
// level 2 interpolation Operator
{

	out.resize(2*in.Nx, 2*in.Ny, 2*in.Nz);


	// ToDo : delete this after testing
	/*
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
	*/


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

				out(i, j, k) = (1.-gamma_z)*out.val_at(i,j,k-1) + gamma_z*out.val_at(i,j,k+1);
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
				out(i, j, k) = (1.-gamma_y)*out.val_at(i,j-1,k) + gamma_y*out.val_at(i,j+1,k);

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
				out(i, j, k) = (1.-gamma_x)*out.val_at(i-1,j,k) + gamma_x*out.val_at(i+1,j,k);
			}
		}
	}


	return;
}



void OP_2xtox(const field_real &in, field_real &out)
{
	if( (in.Nx<4) || (in.Ny<2) || (in.Nz<2))
		throw("system too small");

	if( (in.Nx % 2)!=0 )
		throw("Feld Groesse ungerade");

	out.resize(in.Nx/2, in.Ny, in.Nz);

	for(int i=0; i<out.Nx; i++)
	{
		for(int j=0; j<out.Ny; j++)
		{
			for(int k=0; k<out.Nz; k++)
			{
				int index_in  = in.my_grid->index_at(2*i, j, k);
				int index_out = out.my_grid->index_at(i, j, k);

				out.val[index_out] = in.val[index_in];
			}
		}
	}

	return;
}

void OP_2ytoy(const field_real &in, field_real &out)
{
	if( (in.Nx<2) || (in.Ny<4) || (in.Nz<2))
		throw("system too small");

	if( (in.Ny % 2)!=0 )
		throw("Feld Groesse ungerade");

	out.resize(in.Nx, in.Ny/2, in.Nz);

	for(int i=0; i<out.Nx; i++)
	{
		for(int j=0; j<out.Ny; j++)
		{
			for(int k=0; k<out.Nz; k++)
			{
				int index_in  = in.my_grid->index_at(i, 2*j, k);
				int index_out = out.my_grid->index_at(i, j, k);

				out.val[index_out] = in.val[index_in];
			}
		}
	}

	return;
}

void OP_2ztoz(const field_real &in, field_real &out)
{
	if( (in.Nx<2) || (in.Ny<2) || (in.Nz<4))
		throw("system too small");

	if( (in.Nz % 2)!=0 )
		throw("Feld Groesse ungerade");

	out.resize(in.Nx, in.Ny, in.Nz/2);

	for(int i=0; i<out.Nx; i++)
	{
		for(int j=0; j<out.Ny; j++)
		{
			for(int k=0; k<out.Nz; k++)
			{
				int index_in  = in.my_grid->index_at(i, j, 2*k);
				int index_out = out.my_grid->index_at(i, j, k);

				out.val[index_out] = in.val[index_in];
			}
		}
	}

	return;
}


void OP_2htoh_lvl0(const field_real &in, field_real &out)
/* truncation of a field to another field with half the resolution */
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
/* interpolation of a field onto another with double the resolution */
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

void OP_XhtoYh_lvl1(const field_real &in, field_real &out, const int type, const double &c1value, const double  &c2value)
/* interpolation of a field with arbitrary resolution onto another field
 * with a different arbitrary resolution */
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("OP_XhtoYh_lvl1(const field_real &in, field_real &out, const int type)");
	my_log << "start";
   #endif

	double xmin, xmax, ymin, ymax, zmin, zmax, xlim;

	switch(type)
	{


		case 0 :
			// if input is too small use periodicity
			std::cout << "case 0" << std::endl;
			for(int i=0; i<out.Nx; ++i)
			{
				double x = out.my_grid->x_axis->val_at(i);
				for(int j=0; j<out.Ny; ++j)
				{
					double y = out.my_grid->y_axis->val_at(j);
					for(int k=0; k<out.Nz; ++k)
					{
						double z = out.my_grid->z_axis->val_at(k);
						out(i,j,k) = in(x,y,z);
		               #if defined(_MY_VERBOSE_TEDIOUS)
				    	std::stringstream sstr;
				    	std::string msg;
				    	sstr << "(" << i << "," << j << "," << k << ") ";
				    	sstr << "(" << x << "," << y << "," << z << ") ";
		            	sstr << in(x,y,z);
		            	msg = sstr.str();
		            	my_log << msg;
	                   #endif
					}
				}
			}
			break;

		case 1 :
			// if input is too small continue with given value c1value
			std::cout << "case 1" << std::endl;
			xmin = in.my_grid->x_axis->l0;
			xmax = xmin + in.my_grid->x_axis->L;
			ymin = in.my_grid->y_axis->l0;
			ymax = ymin + in.my_grid->y_axis->L;
			zmin = in.my_grid->z_axis->l0;
			zmax = zmin + in.my_grid->z_axis->L;

			for(int i=0; i<out.Nx; ++i)
			{
				double x = out.my_grid->x_axis->val_at(i);
				for(int j=0; j<out.Ny; ++j)
				{
					double y = out.my_grid->y_axis->val_at(j);
					for(int k=0; k<out.Nz; ++k)
					{
						double z = out.my_grid->z_axis->val_at(k);
						out(i,j,k) = c1value;
						if( ((xmin <= x) && (x <xmax)) &&
							((ymin <= y) && (y <ymax)) &&
							((zmin <= z) && (z <zmax)) )
						{
							out(i,j,k) = in(x,y,z);
						}

			           #if defined(_MY_VERBOSE_TEDIOUS)
						std::stringstream sstr;
						std::string msg;
						sstr << "(" << i << "," << j << "," << k << ") ";
						sstr << "(" << x << "," << y << "," << z << ") ";
			            sstr << in(x,y,z);
			            msg = sstr.str();
			            my_log << msg;
		               #endif
					}
				}
			}
		break;


		case 2 :
			// mixed case of 0 and 1:
			// if input is too small use periodicity in positive x-direction for x .ge. c2value and
			// continue with given value c1value in negative xdirection for x .lt. c2value
			std::cout << "case 1" << std::endl;
			xmin = in.my_grid->x_axis->l0;
			xmax = xmin + in.my_grid->x_axis->L;
			xlim = c2value;
			ymin = in.my_grid->y_axis->l0;
			ymax = ymin + in.my_grid->y_axis->L;
			zmin = in.my_grid->z_axis->l0;
			zmax = zmin + in.my_grid->z_axis->L;

			for(int i=0; i<out.Nx; ++i)
			{
				double x = out.my_grid->x_axis->val_at(i);
				for(int j=0; j<out.Ny; ++j)
				{
					double y = out.my_grid->y_axis->val_at(j);
					for(int k=0; k<out.Nz; ++k)
					{
						double z = out.my_grid->z_axis->val_at(k);
						out(i,j,k) = c1value;
						if(  (xlim <= x) &&
							((ymin <= y) && (y <ymax)) &&
							((zmin <= z) && (z <zmax)) )
						{
							out(i,j,k) = in(x,y,z);
						}
						else
						{
							out(i,j,k) = c1value;
						}

			           #if defined(_MY_VERBOSE_TEDIOUS)
						std::stringstream sstr;
						std::string msg;
						sstr << "(" << i << "," << j << "," << k << ") ";
						sstr << "(" << x << "," << y << "," << z << ") ";
			            sstr << in(x,y,z);
			            msg = sstr.str();
			            my_log << msg;
		               #endif
					}
				}
			}
		break;

	}

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}




double get_PG(const field_real &in, int i, int j, int k)
// returns PHI or boundary value depending on position (i,j,k)
{
	while(i<0)
		i+=in.Nx;
	while(i>=in.Nx)
		i-=in.Nx;

	while(j<0)
		j+=in.Ny;
	while(j>=in.Ny)
		j-=in.Ny;

	while(k<0)
		k+=in.Nz;
	while(k>=in.Nz)
		k-=in.Nz;

	return in(i,j,k);
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

	for(int i=0; i<in.Nx; i++)
	{
		for(int j=0; j<in.Ny; j++)
		{
			for(int k=0; k<in.Nz; k++)
			{
				int index  = tmp.my_grid->index_at( i, j, k);

				tmp.val[index] = (2./8.)*get_PG(in,i,j,k);
				tmp.val[index]+= (1./8.)*get_PG(in,i+1,j,k);
				tmp.val[index]+= (1./8.)*get_PG(in,i-1,j,k);
				tmp.val[index]+= (1./8.)*get_PG(in,i,j+1,k);
				tmp.val[index]+= (1./8.)*get_PG(in,i,j-1,k);
				tmp.val[index]+= (1./8.)*get_PG(in,i,j,k+1);
				tmp.val[index]+= (1./8.)*get_PG(in,i,j,k-1);
			}
		}
	}

	for(int i=0; i<in.N; ++i)
		out.val[i] = tmp.val[i];

	return;
}


void OP_smoothing_lvl2(const field_real &in, field_real &out)
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

	for(int i=0; i<in.Nx; i++)
	{
		for(int j=0; j<in.Ny; j++)
		{
			for(int k=0; k<in.Nz; k++)
			{
				int index  = tmp.my_grid->index_at( i, j, k);

				tmp.val[index] = 3.*get_PG(in,i,j,k);
				tmp.val[index]+= 1.*get_PG(in,i+2,j,k);
				tmp.val[index]+= 1.*get_PG(in,i+1,j,k);
				tmp.val[index]+= 1.*get_PG(in,i-1,j,k);
				tmp.val[index]+= 1.*get_PG(in,i-2,j,k);
				tmp.val[index]+= 1.*get_PG(in,i,j+2,k);
				tmp.val[index]+= 1.*get_PG(in,i,j+1,k);
				tmp.val[index]+= 1.*get_PG(in,i,j-1,k);
				tmp.val[index]+= 1.*get_PG(in,i,j-2,k);
				tmp.val[index]+= 1.*get_PG(in,i,j,k+2);
				tmp.val[index]+= 1.*get_PG(in,i,j,k+1);
				tmp.val[index]+= 1.*get_PG(in,i,j,k-1);
				tmp.val[index]+= 1.*get_PG(in,i,j,k-2);
				tmp.val[index] = (1./15.)*tmp.val[index];
			}
		}
	}

	for(int i=0; i<in.N; ++i)
		out.val[i] = tmp.val[i];

	return;
}
