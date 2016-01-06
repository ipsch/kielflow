#ifndef RKO4_HPP_
#define RKO4_HPP_



#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <omp.h>
#include <fstream>

#include "interface_ODE_solver.hpp"
#include "interface_rhs.hpp"
#include "field_imag.hpp"
#include "CFL.hpp"
#include "operations.hpp"

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif
// magic stuff





class Runge_kutta_O4 : public interface_ode_solver
{

public :
	Runge_kutta_O4(interface_rhs &rhs, double dt = 0.01);
	~Runge_kutta_O4() {	}
	void solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni);

private :
	std::string my_logfile;
	double t_;
	double dt_;
	interface_rhs &my_rhs;
};


#endif /* END RKO4_HPP_ */
