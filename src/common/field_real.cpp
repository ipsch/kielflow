#include "field_real.hpp"


bool field_real::IsNan() const
{
	for(int i=0; i<N; ++ i)
		if(val[i]!=val[i])
			return true;
	return false;
}

bool field_real::IsInf() const
{
	for(int i=0; i<N; ++ i)
		if(std::isinf(val[i]))
			return true;
	return false;
}


// ###################### coordinate space #############################
field_real::field_real(const grid &Omega) :
field(Omega, Omega.Nx, Omega.Ny, Omega.Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "field_real(const grid_Co &Omega)";
   #endif

	val = (double*) fftw_malloc(sizeof(double) * N);
}



field_real::field_real(const field_real &that) :
		field(that.my_grid, that.my_grid.Nx, that.my_grid.Ny, that.my_grid.Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_real");
	log << "field_real(const grid &Omega)";
   #endif

	val = (double*) fftw_malloc(sizeof(double) * N);

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
				double tmp_x = my_grid.x_axis->val_at(i);
				double tmp_y = my_grid.y_axis->val_at(j);
				double tmp_z = my_grid.z_axis->val_at(k);
				val[index(i,j,k)] = fill_fkt(tmp_x,tmp_y,tmp_z);
			}
	return;
}

void field_real::fill(interface_3d_fkt const &rhs)
{
	for(int i=0; i<Nx; ++i)
	{
		double x = my_grid.x_axis->val_at(i);
		for(int j=0; j<Ny; ++j)
		{
			double y = my_grid.y_axis->val_at(j);
			for(int k=0; k<Nz; ++k)
			{
				double z = my_grid.z_axis->val_at(k);
				val[index(i,j,k)] = rhs(x,y,z);
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
				double tmp_x = my_grid.x_axis->val_at(i);
				double tmp_y = my_grid.y_axis->val_at(j);
				double tmp_z = my_grid.z_axis->val_at(k);
				val[index(i,j,k)] += fill_fkt(tmp_x,tmp_y,tmp_z);
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



double field_real::val_at(int ix, int iy, int iz) const
{
	while(ix<0) ix+=Nx;
	while(ix>Nx-1) ix-=Nx;
	while(iy<0) iy+=Ny;
	while(iy>Ny-1) iy-=Ny;
	while(iz<0) iz+=Nz;
	while(iz>Nz-1) iz-=Nz;
	return val[index(ix, iy, iz)];
}

//inline
double& field_real::operator() (const int &ix, const int &iy, const int &iz)
{
	return val[index(ix, iy, iz)];
}


//inline
double field_real::operator() (const int &ix, const int &iy, const int &iz) const
{
	return val[index(ix, iy, iz)];
}


inline
double& field_real::operator() (int N)
{
	// ToDo: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("Matrix subscript out of bounds");
return val[N];
}


double field_real::operator() (const double &x, const double &y, const double &z) const
{
	int i = my_grid.x_axis->index_at(x); // index_at returns nearest int
	int j = my_grid.y_axis->index_at(y);
	int k = my_grid.z_axis->index_at(z);
	return val[index(i,j,k)];
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
	my_grid.resize(NX,NY,NZ);
	Nx=NX;
	Ny=NY;
	Nz=NZ;
	N=Nx*Ny*Nz;
	val = (double*) fftw_malloc(sizeof(double) * N);

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
