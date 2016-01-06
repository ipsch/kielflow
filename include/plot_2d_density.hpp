#ifndef PLOT_2D_DENSITY_HPP_
#define PLOT_2D_DENSITY_HPP_

#include "plot_2d.hpp"


class density_plot: public plot_2d
{
public :
	density_plot(const field_real &Fdata, const parameters &Pdata);
	void auto_set(const field_real &field);
	std::string script_density(void);
	void save(void) {save_(this->script_density());}
private :
};






#endif
