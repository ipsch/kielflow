#ifndef INTERFACE_RHS_HPP_
#define INTERFACE_RHS_HPP_

#include "parameters.hpp"
#include "field_imag.hpp"

class interface_rhs
{
public :
	virtual ~interface_rhs() { };
	virtual void solve(const double &t, field_imag &FUx, field_imag &FUy, field_imag &FUz,
			field_imag &Fni) = 0;
protected :
	interface_rhs() { };
};





#endif /* END INTERFACE_RHS_HPP_ */
