#ifndef AXIS_BASE_HPP_
#define AXIS_BASE_HPP_

#include "axis_interface.hpp"
#include <iostream>
#include <cmath>
#include <tr1/memory>

#if defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif



class grid;

typedef enum
{
    e_x, e_y, e_z,
    e_r, e_phi, e_theta
} direction;




// parent class for all types of axis
class axis : public axis_interface
{
public :

	virtual ~axis() = 0;

	unsigned int type_id;
	direction e_i;

	int N;
	double l0;
	double L;
	double m;


	virtual axis * clone() const = 0;
	virtual void resize(const int &N_) = 0;

	virtual double dS(const int &index) const;
protected :
	axis();
	void set_type(const int &id) {type_id = id;}
	int  get_type(void) const {return type_id;}

};

typedef axis * (axis::*ptr_axis_factory_fkt)() const;



#endif
