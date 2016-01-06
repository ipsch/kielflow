#include "field_real.hpp"


// ###################### coordinate space #############################
field_real::field_real(const grid_Co &Omega) :
my_grid(Omega.clone()),
N(my_grid->N), Nx(my_grid->Nx), Ny(my_grid->Ny), Nz(my_grid->Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "field_real(const grid_Co &Omega)";
   #endif

	val = (double*) fftw_malloc(sizeof(double) * my_grid->N);
}

field_real::field_real(const grid_Fo &FOmega) :
my_grid(FOmega.reziprocal()),
N(my_grid->N), Nx(my_grid->Nx), Ny(my_grid->Ny), Nz(my_grid->Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "field_real(const grid_Fo &FOmega)";
   #endif

	val = (double*) fftw_malloc(sizeof(double) * my_grid->N);
}

field_real::field_real(const field_real &that) :
my_grid(that.my_grid->clone()),
N(my_grid->N), Nx(my_grid->Nx), Ny(my_grid->Ny), Nz(my_grid->Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "field_real(const grid &Omega)";
   #endif

	val = (double*) fftw_malloc(sizeof(double) * my_grid->N);

	for( int i=0; i<that.N; ++i)
	{
		this->val[i] = that.val[i];
	}
}

field_real::~field_real()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "~field_real()";
   #endif
	fftw_free(val);
	delete my_grid;
}

void field_real::fill(double (*fill_fkt)(const double &,const double &, const double &)) const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "fill(double (*fill_fkt)(double,double,double))";
   #endif
	for( int i=0; i < Nx; ++i)
		for( int j=0; j < Ny; ++j)
			for( int k=0; k < Nz; ++k)
			{
				double tmp_x = my_grid->x_axis->val_at(i);
				double tmp_y = my_grid->y_axis->val_at(j);
				double tmp_z = my_grid->z_axis->val_at(k);
				val[my_grid->index_at(i,j,k)] = fill_fkt(tmp_x,tmp_y,tmp_z);
			}
	return;
}

void field_real::fill(interface_3d_fkt const &rhs)
{
	for(int i=0; i<Nx; ++i)
	{
		double x = my_grid->x_axis->val_at(i);
		for(int j=0; j<Ny; ++j)
		{
			double y = my_grid->y_axis->val_at(j);
			for(int k=0; k<Nz; ++k)
			{
				double z = my_grid->z_axis->val_at(k);
				val[my_grid->index_at(i,j,k)] = rhs(x,y,z);
			}
		}
	}

	return;
}

void field_real::add(double (*fill_fkt)(const double &,const double &, const double &)) const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "add(double (*fill_fkt)(double,double,double))";
   #endif
	for( int i=0; i < Nx; ++i)
		for( int j=0; j < Ny; ++j)
			for( int k=0; k < Nz; ++k)
			{
				double tmp_x = my_grid->x_axis->val_at(i);
				double tmp_y = my_grid->y_axis->val_at(j);
				double tmp_z = my_grid->z_axis->val_at(k);
				val[my_grid->index_at(i,j,k)] += fill_fkt(tmp_x,tmp_y,tmp_z);
			}
	return;
}

void field_real::mulitply(const double &lambda)
{
	for( int i=0; i<N; ++i)
	{
		val[i] *= lambda;
	}

}




//inline
double& field_real::operator() (const int &ix, const int &iy, const int &iz)
{
	// ToDo: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("Matrix subscript out of bounds");
	// ToDo : secure als Flag einbauen


	return val[my_grid->index_at(ix % Nx,iy % Ny,iz % Nz)];
	//return val[my_grid->index_at(ix,iy,iz)];
}


//inline
double field_real::operator() (const int &ix, const int &iy, const int &iz) const
{
	// ToDO: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("const Matrix subscript out of bounds");
	// ToDo : secure als Flag einbauen

	return val[my_grid->index_at(ix % Nx,iy % Ny,iz % Nz)];
	//return val[my_grid->index_at(ix,iy,iz)];
}


inline
double& field_real::operator() (int N)
{
	// ToDo: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("Matrix subscript out of bounds");
return val[N];
}






field_real& field_real::operator= (const field_real &rhs)
{
	if(this==&rhs)
		return *this;
	if(this->Nx!= rhs.Nx)
		throw("field size not matching");
	if(this->Ny!= rhs.Ny)
		throw("field size not matching");
	if(this->Nz!= rhs.Nz)
		throw("field size not matching");

	// ToDO : secure flag hier!
	for(int i=0; i<rhs.N; ++i)
		val[i] = rhs.val[i];
	return *this;
}


field_real& field_real::operator+= (const field_real &rhs)
{
	for(int i=0; i<rhs.N; ++i)
		val[i] += rhs.val[i];
	return *this;
}


//void field_real::resize(const int &NX, const int &NY, const int &NZ)
void field_real::resize(int NX, int NY, int NZ)
{
	fftw_free(val);
	my_grid->set_resolution(NX,NY,NZ);
	val = (double*) fftw_malloc(sizeof(double) * my_grid->N);

	return;
}


/*
inline
double field_real::operator() (int N) const
{
	// ToDO: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("const Matrix subscript out of bounds");
return val[N];
}
*/
