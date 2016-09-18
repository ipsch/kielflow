#ifndef SOLVER_RHS_NACA_HPP_
#define SOLVER_RHS_NACA_HPP_

#include "rhs_base.hpp"
#include "effect.hpp"
#include <string>
#include <sstream>

class solver_rhs_naca : public solver_rhs
{
public :
	solver_rhs_naca(const parameters &Params, solver_poisson &Solver);
	void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni, field_imag &FPh);
	~solver_rhs_naca() {	}
};




#endif
