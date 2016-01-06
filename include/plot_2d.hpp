#ifndef PLOT_2D_HPP_
#define PLOT_2D_HPP_


#include "plot_base.hpp"


class plot_2d : public plot
{
public :
	plot_2d(const field_real &Fdata, const parameters &Pdata);
	~plot_2d();

	field_real my_field;

	parameters my_parameters;

	static int default_plane;
	static int default_direction;
	static int default_axis;
	static int default_xpos;
	static int default_ypos;
	static int default_zpos;

	std::string path;
	std::string file;
	std::string title;
	std::string x_label;
	std::string y_label;
	std::string z_label;
	std::string cb_label;
	double cb_max;
	double cb_min;

	int dim2_plane_int;
	std::string gnuplot_plane;
	std::string gnuplot_v_vec;
	void save_(std::string script);
protected :

};



void save_density_pl4(const grid_Co &Omega, const parameters &Params, std::string file_name);
void save_velocity_pl4(const grid_Co &Omega, const parameters &Params, std::string file_name);




#endif
