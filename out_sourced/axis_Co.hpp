#ifndef AXIS_CO_HPP_
#define AXIS_CO_HPP_

#include "axis_base.hpp"
#include "axis_Fo.hpp"



// ############################ parent Class for all axis types of coordinate-space type
class axis_Co : public axis
{
friend class grid;
public :
	//virtual axis_Co * clone() const = 0;
protected :
	//virtual axis_Fo * create_reciprocal() const = 0;
};

//typedef axis * (axis_Co::*ptr_axis_Co_factory_fkt)() const;

#endif



