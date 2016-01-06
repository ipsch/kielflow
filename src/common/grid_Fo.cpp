#include "grid_Fo.hpp"


// ###################### grid in Fourier Space ########################

grid_Fo::grid_Fo(const axis_Fo &x, const axis_Fo &y, const axis_Fo &z) :
		grid(x.N, y.N, z.N/2+1, x, y, z, &axis::clone)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Fo");
	log << "grid_Fo(const axis_Fo &x, const axis_Fo &y, const axis_Fo &z)";
   #endif
}


grid_Fo::grid_Fo(const grid_Fo &that) :
		grid(that.x_axis->N, that.y_axis->N, that.z_axis->N/2+1,
			*that.x_axis,   *that.y_axis,   *that.z_axis,
				&axis::clone)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Fo");
	log << "grid_Fo(const grid_Fo &that)";
   #endif
}


grid_Fo::grid_Fo(const grid_Co &that) :
		grid(that.x_axis->N, that.y_axis->N, that.z_axis->N/2+1,
			*that.x_axis,   *that.y_axis,   *that.z_axis,
				&axis_Co::create_reciprocal)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Fo");
	log << "grid_Fo(const grid_Co &that)";
   #endif
}


grid_Fo * grid_Fo::clone() const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Fo");
	log << "grid_Fo::clone() const";
   #endif
	return new grid_Fo(*this);
}


grid_Co * grid_Fo::reziprocal() const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Fo");
	log << "grid_Co::reziprocal() const";
   #endif
	return new grid_Co(*this);
}


grid_Fo::~grid_Fo()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Fo");
	log << "~grid_Fo()";
   #endif
}


 int grid_Fo::index_at(const int &i, const  int &j, const int &k) const
{
	return k + Nz*(j + i*Ny);
	//return k + (Nz/2 +1)*(j + i*Ny);
}


void grid_Fo::ijk_at(const  int &index,  int &i,  int &j,  int &k) const
{
	i =  index/(y_axis->N*z_axis->N);
	j = (index - i*(y_axis->N*z_axis->N)) / z_axis->N;
	k =  index - i*(y_axis->N*z_axis->N) - j*z_axis->N;
	return;
}





// ######### Operatoren #######################################################
bool operator==(const grid_Co &c1, const grid_Co &c2)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid");
	log << "operator==(const grid_Co &c1, const grid_Co &c2)";
   #endif
    if (c1.N!=c2.N) return false;
    if (c1.x_axis->N!=c2.x_axis->N) return false;
    if (c1.y_axis->N!=c2.y_axis->N) return false;

    return true;
}

bool operator==(const grid_Fo &c1, const grid_Fo &c2)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid");
	log << "operator==(const grid &c1, const grid &c2)";
   #endif
    if (c1.N!=c2.N) return false;
    if (c1.x_axis->N!=c2.x_axis->N) return false;
    if (c1.y_axis->N!=c2.y_axis->N) return false;

    return true;
}

/*
int domain::index_real(const int &i, const int &j, const int &k)
{

	return (k + z_axis->N *(j + i*y_axis->N));
}


int domain::index_imag(const int &i, const int &j, const int &k)
{
	return (k + (z_axis->N/2 +1)*(j + i*y_axis->N));
}
*/

/*
 int domain::Nx = 0;
 int domain::Ny = 0;
 int domain::Nz = 0;
 int domain::N_real_ges = Nx*Ny*Nz;
 int domain::N_imag_ges = Nx*Ny*(Nz/2 +1);
double domain::x0 = 0;
double domain::y0 = 0;
double domain::z0 = 0;
double domain::Lx = 0;
double domain::Ly = 0;
double domain::Lz = 0;
double domain::Lt = 0;
double domain::old_dt = 0;
double domain::old_t = 0;

double domain::grid_scale_factor = 3.;

double domain::k1(int i, int N, double L)
{
	double pi = acos(-1);
	return ( (i<=N/2) ? (2*pi*i/L) : (2*pi*(i-N)/L) );
}

double domain::k2(const int &i, const int &N, const double &L)
{
	double pi = acos(-1);
	return ( (i<=N/2) ? (-pow((2.*pi*i/L),2)) : (-pow((2*pi*(i-N)/L),2.)) );
}

double domain::k2(const int &i, const int &j, const int &k)
{
	return (k2(i, domain::Nx, domain::Lx) + k2(j, domain::Ny, domain::Ly) + k2(k, domain::Nz, domain::Lz));
}


#ifndef _NEQ_GRID
// see domain.hpp for inline funktions
#else

double domain::grid_scale_function(const double &x, const double &L)
{
	return (0.5*L)*sinh(2.*x*grid_scale_factor/L)/sinh(grid_scale_factor);
}

double domain::x_at(const int &i)
{
	double x_val;

	x_val = x0 + Lx*double(i)/(Nx);

	if     (x_val>= .5*Lx) x_val = Lx+grid_scale_function(x_val-Lx,Lx);
	else if(x_val< -.5*Lx) x_val =-Lx+grid_scale_function(x_val+Lx,Lx);
	else                   x_val =    grid_scale_function(x_val   ,Lx);

	return x_val;
}

double domain::y_at(const int &j)
{
	return y0 + j*(Ly/Ny);
}

double domain::z_at(const int &k)
{
	return z0 + k*(Lz/Nz);
}

int domain::i_at(double x)
{
	return Nx*(x-x0)/Lx;
}

int domain::j_at(double y)
{
	return Ny*(y-y0)/Ly;
}

int domain::k_at(double z)
{
	return Nz*(z-z0)/Lz;
}

#endif
*/


