#ifndef EULER_METHOD_HPP_
#define EULER_METHOD_HPP_

/* Euler-method to solve ODEs */

#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <omp.h>
#include <fstream>

#include "interface_ODE_solver.hpp"
#include "interface_rhs.hpp"
#include "field.hpp"
#include "field_imag.hpp"
#include "CFL.hpp"
#include "operations.hpp"
#include "subdim.hpp"
#include "IO.hpp"

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif




class euler_method : public interface_ode_solver
{
public :
	euler_method(interface_rhs &rhs, double dt = 0.1);
	~euler_method() {	}
	void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni);
private :
	//euler_method() : rhs(NULL), dt_(0.), t_(0) { }; // can never be called!
	std::string my_logfile;
	double t_;
	double dt_;
	interface_rhs &my_rhs;
};


#endif /* END EULER_METHOD_HPP_ */
