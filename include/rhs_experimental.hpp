#ifndef SOLVER_RHS_EXPERIMENTAL_HPP_
#define SOLVER_RHS_EXPERIMENTAL_HPP_

#include "rhs_base.hpp"
#include "effect.hpp"

class solver_rhs_experimental : public solver_rhs
{
public :
	solver_rhs_experimental(const parameters &Params, solver_poisson &Solver);
	void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni, field_imag &FPh);
	~solver_rhs_experimental() {	}

};




#endif
