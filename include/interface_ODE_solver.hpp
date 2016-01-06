#ifndef INTERFACE_ODE_SOLVER_HPP_
#define INTERFACE_ODE_SOLVER_HPP_

#include "field_imag.hpp"

class interface_ode_solver
{

public :
	virtual ~interface_ode_solver() { };
	virtual void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz,
			field_imag &Fni) = 0;
protected :
	interface_ode_solver() { };
};






#endif /* END INTERFACE_ODE_SOLVER_HPP_ */
