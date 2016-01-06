#include "grid_base.hpp"



// ###################### parent grid ########################



grid::grid(const int &Nx_, const int &Ny_, const int &Nz_,
		const axis &Ax, const axis &Ay, const axis &Az,
		const ptr_axis_factory_fkt &factory_fkt) :
		Nx(Nx_), Ny(Ny_), Nz(Nz_), N(Nx_*Ny_*Nz_),
		x_axis( (Ax.*factory_fkt)()), y_axis( (Ay.*factory_fkt)()), z_axis( (Az.*factory_fkt)() )
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid");
	log << "grid(const int &Nx_, const int &Ny_, const int &Nz_, axis &Ax, axis &Ay, axis &Az, ptr_axis_factory_fkt &factory_fkt) ";
   #endif
}


grid::grid(const grid &that) :
Nx(that.Nx), Ny(that.Ny), Nz(that.Nz), N(that.N),
x_axis(that.x_axis->clone()), y_axis(that.y_axis->clone()), z_axis(that.z_axis->clone())

{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid");
	log << "grid(const  grid &that)";
   #endif
}


grid::~grid()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("grid");
	log << "~grid()";
   #endif
	delete x_axis;
	delete y_axis;
	delete z_axis;
};








