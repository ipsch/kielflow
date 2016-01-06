#ifndef PLOT_2D_VELOCITY_HPP_
#define PLOT_2D_VELOCITY_HPP_

#include "plot_2d.hpp"



class velocity_plot: public plot_2d
{
public :
	velocity_plot(const field_real &xFdata, const field_real &yFdata, const field_real &zFdata, const parameters &Pdata);
	void auto_set (double *field);
	std::string script_velocity(void);
	void save(void) {save_(this->script_velocity());}
private :
};





#endif
