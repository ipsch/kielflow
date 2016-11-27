#ifndef SOLVER_RHS_BASE_HPP_
#define SOLVER_RHS_BASE_HPP_


//#define _OLD_FORCE
//#define _NEW_FORCE
//#define DIFFUSION_TRUE
//#define B_FIELD_TRUE

// gcc standard C++ libraries (common stuff)
#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)



// own includes

#include "parameters.hpp"
#include "field.hpp"
#include "operations.hpp"
#include "solver_poisson.hpp"
#include "masks.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif
// magic stuff

#include <omp.h>





class solver_rhs
{


public :

	solver_rhs(const parameters &Params, solver_poisson &Solver);
	virtual void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni, field_imag &FPh) = 0;
	virtual ~solver_rhs() {	}

protected :
	parameters my_params;
	solver_poisson * const my_poisson_solver;



	void evaluate_advection(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz, const field_imag FXX,
			field_imag &return_Buffer);
	//void evaluate_transport(const field_imag &FXX, field_imag &Buffer_FXX);
	void evaluate_dissipation(const field_imag &FXX,
			field_imag &Buffer_FXX);
	void evaluate_force_E(const field_imag &FPhi,
			field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz);
	void evaluate_force_B(const field_imag &FUy, const field_imag &FUz,
			field_imag &Buffer_FUy, field_imag &Buffer_FUz);
	void evaluate_diffusion_mu(const field_imag &FXX, field_imag &Buffer_FXX);
	void evaluate_diffusion_u(const field_imag &Fni,
			field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz);

	void evaluate_diffusion_n(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
			field_imag &Buffer_Fni);
	void evaluate_continuityEQ(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni,
			field_imag &Buffer_Fni);


};




#endif
