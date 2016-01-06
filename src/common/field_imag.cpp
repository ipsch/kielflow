#include "field_imag.hpp"




// ###################### fourier space #############################


field_imag::field_imag(const grid_Fo &FOmega) :
my_grid(FOmega.clone()),
N(my_grid->N), Nx(my_grid->Nx), Ny(my_grid->Ny), Nz(my_grid->Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_imag");
	log << "field_imag(const grid_Fo &FOmega)";
   #endif

	val = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * my_grid->N);


}

field_imag::field_imag(const grid_Co &Omega) :
my_grid(Omega.reziprocal()),
N(my_grid->N), Nx(my_grid->Nx), Ny(my_grid->Ny), Nz(my_grid->Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_imag");
	log << "field_imag(const grid_Co &Omega)";
   #endif

	val = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * my_grid->N);

}


field_imag::field_imag(const field_imag &that) :
my_grid(that.my_grid->clone()),
N(my_grid->N), Nx(my_grid->Nx), Ny(my_grid->Ny), Nz(my_grid->Nz)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_imag");
	log << "field_imag(const field_imag &that)";
   #endif

	val = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * my_grid->N);


	for( int i=0; i<that.N; ++i)
	{
		this->val[i][0] = that.val[i][0];
		this->val[i][1] = that.val[i][1];
	}
}

field_imag::~field_imag()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field_imag");
	log << "~field_imag()";
   #endif

	fftw_free(val);
	delete my_grid;

   #ifdef _MY_VERBOSE_MORE
	log << "done";
   #endif
}

void field_imag::fill_RE(double (*fill_fkt)(const double &,const double &, const double &)) const
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
				val[my_grid->index_at(i,j,k)][0] = fill_fkt(tmp_x,tmp_y,tmp_z);
			}
	return;
}

void field_imag::fill_IM(double (*fill_fkt)(const double &,const double &, const double &)) const
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
				val[my_grid->index_at(i,j,k)][1] = fill_fkt(tmp_x,tmp_y,tmp_z);
			}
	return;
}

inline
fftw_complex& field_imag::operator() (int ix, int iy, int iz)
{
	// ToDo: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("Matrix subscript out of bounds");
return val[my_grid->index_at(ix,iy,iz)];
}

/*
inline
fftw_complex& field_imag::operator() (int ix, int iy, int iz) const
{
	// ToDO: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("const Matrix subscript out of bounds");
return val[my_grid->index_at(ix,iy,iz)];
}
*/

inline
fftw_complex& field_imag::operator() (const int &N)
{
	// ToDo: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("Matrix subscript out of bounds");
return val[N];
}

/*
inline
fftw_complex& field_imag::operator() (const int &N) const
{
	// ToDO: Fehler abfangen
	//if (row >= rows_ || col >= cols_)
	//throw BadIndex("const Matrix subscript out of bounds");
return val[N];
}
*/

double& field_imag::operator() (int ix, int iy, int iz, int re_im)
{

	return val[my_grid->index_at(ix,iy,iz)][re_im];
}

double& field_imag::operator() (const int &N, int re_im)
{
	return val[N][re_im];
}
