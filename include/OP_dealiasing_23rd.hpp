#ifndef OP_DEALIASING_23RD_HPP_
#define OP_DEALIASING_23RD_HPP_

#include <iostream>
#include <cmath>
#include "grid.hpp"
#include "field_imag.hpp"

class OP_dealiasing_23rd
{
public :
	OP_dealiasing_23rd(const grid &domain);
	~OP_dealiasing_23rd();
	void operator() (field_imag &target) const;
private :
	const unsigned int N;
	double * filter;
	OP_dealiasing_23rd() : N(0) {filter = 0L;}
};


#endif
