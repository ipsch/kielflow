#include "grid_Co.hpp"

// ###################### grid in Coordinate Space ########################













grid_Co::grid_Co(const axis_Co &x, const axis_Co &y, const axis_Co &z) :
		grid(x.N, y.N, z.N, x, y, z, &axis::clone)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field");
	log << "grid_Co(const axis_Co &x, const axis_Co &y, const axis_Co &z)";
   #endif
}

grid_Co::grid_Co(const grid_Co &that) :
		grid(that.x_axis->N, that.y_axis->N, that.z_axis->N,
			*that.x_axis,   *that.y_axis,   *that.z_axis,
				&axis::clone)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field");
	log << "grid_Co(const grid_Co &that)";
   #endif
}

grid_Co::grid_Co(const grid_Fo &that) :
		grid(that.x_axis->N, that.y_axis->N, that.z_axis->N,
			*that.x_axis,   *that.y_axis,   *that.z_axis,
				&axis::create_reciprocal)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field");
	log << "grid_Co(const grid_Fo &that) ";
   #endif
}

grid_Co * grid_Co::clone() const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Co");
	log << "grid_Co::clone() const";
   #endif
	return new grid_Co(*this);
}

grid_Fo * grid_Co::reziprocal() const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Co");
	log << "grid_Co::reziprocal() const";
   #endif
	return new grid_Fo(*this);
}

grid_Co::~grid_Co()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid_Co");
	log << "~grid_Co()";
   #endif
}


 int grid_Co::index_at(const int &i, const int &j, const int &k) const
{
	return (k + z_axis->N *(j + i*y_axis->N));
}

void grid_Co::ijk_at(const  int &index,  int &i,  int &j,  int &k) const
{
	i =  index/(y_axis->N*z_axis->N);
	j = (index - i*(y_axis->N*z_axis->N)) / z_axis->N;
	k =  index - i*(y_axis->N*z_axis->N) - j*z_axis->N;
	return;
}

void grid_Co::set_resolution(const int &N_x, const int &N_y, const int &N_z)
{
	x_axis->N = Nx = N_x;
	y_axis->N = Ny = N_y;
	z_axis->N = Nz = N_z;
	N = (Nx*Ny*Nz);
}
