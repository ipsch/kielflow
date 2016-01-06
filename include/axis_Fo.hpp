#ifndef AXIS_FO_HPP_
#define AXIS_FO_HPP_

#include "axis_base.hpp"
#include "axis_Co.hpp"


class axis_Co;
// ############################ parent Class for all axis types of fourier-space type
class axis_Fo : public axis
{
friend class grid;
public :
	//virtual axis_Fo * clone() const = 0;
protected :
	//virtual axis_Co * create_reciprocal() const = 0;
};

//typedef axis_Fo * (axis::*ptr_axis_Fo_factory_fkt)() const;

#endif
