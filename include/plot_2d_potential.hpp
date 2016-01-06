#ifndef PLOT_2D_POTENTIAL_HPP_
#define PLOT_2D_POTENTIAL_HPP_

#include "plot_2d.hpp"




class potential_plot: public plot_2d
{
public :
	potential_plot(const field_real &Fdata, const parameters &Pdata);
	void auto_set(const field_real &field);
	std::string script_potential(void);
	void save(void) {save_(this->script_potential());}
private :
};





#endif
