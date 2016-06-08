#ifndef EULER_METHOD_HPP_
#define EULER_METHOD_HPP_

/*
 *  Euler-method to solve ODEs
 *
 *
 */


#include <omp.h>

#include "interface_ODE_solver.hpp"
#include "interface_rhs.hpp"
#include "field_imag.hpp"
#include "CFL.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif





class euler_method : public interface_ode_solver
{

public :
	euler_method(interface_rhs &rhs, double dt = 0.1);
	~euler_method() {	}
	void solve(const double &t, field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni);

private :
	//euler_method() : rhs(NULL), dt_(0.), t_(0) { }; // can never be called!
	double t_;
	double dt_;

	interface_rhs &my_rhs;
};


#endif /* END EULER_METHOD_HPP_ */
