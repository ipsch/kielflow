#include "grid_base.hpp"



// ###################### parent grid ########################



grid::grid(const axis &Ax, const axis &Ay, const axis &Az) :
	x_axis(Ax.clone()), y_axis(Ay.clone()), z_axis(Az.clone()),
	Nx(x_axis->N), Ny(y_axis->N), Nz(z_axis->N)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("grid");
	my_log << "grid(const int &Nx_, const int &Ny_, const int &Nz_, axis &Ax, axis &Ay, axis &Az, ptr_axis_factory_fkt &factory_fkt) ";
   #endif
}


grid::grid(const grid &that) :
	grid(*that.x_axis, *that.y_axis, *that.z_axis)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("grid");
	my_log << "grid(const  grid &that)";
   #endif
}


grid::~grid()
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("grid");
	my_log << "~grid()";
   #endif
	delete x_axis;
	delete y_axis;
	delete z_axis;
};








