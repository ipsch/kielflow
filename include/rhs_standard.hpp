#ifndef RHS_STANDARD_HPP_
#define RHS_STANDARD_HPP_


#include "interface_rhs.hpp"
#include "interface_poisson_solver.hpp"

#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <omp.h>



#include "parameters.hpp"
#include "field.hpp"
#include "operations.hpp"
#include "effect.hpp"
#include "masks.hpp"


#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif
// magic stuff




class rhs_standard : public interface_rhs
{
public :
	rhs_standard(parameters &Params, interface_poisson_solver &Solver,
			field_real &potential, field_real &density_dust,
			field_real &solid_mask);
	void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni);
	~rhs_standard() {	}

protected :
	parameters &my_params;
	interface_poisson_solver &my_poisson_solver;
	field_real &Phi;
	field_real &nd;

	field_real &my_solids;


};



#endif /* RHS_STANDARD_HPP_ */



