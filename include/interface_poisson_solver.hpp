#ifndef INTERFACE_POISSON_SOLVER_HPP_
#define INTERFACE_POISSON_SOLVER_HPP_

/*
 *
 *   Interface Class to a Poisson solver
 *
 *
 */

#include "field_real.hpp"


class interface_poisson_solver
{
public :
	virtual ~interface_poisson_solver() { };
	virtual void solve(field_real &Phi_IO, field_real &rho) = 0;
protected :
	interface_poisson_solver() {};
};


class interface_relaxation_solver : public interface_poisson_solver
{
public :
	virtual ~interface_relaxation_solver() { };
	virtual void solve(field_real &Phi_IO, field_real &rho) = 0;
	virtual void set_max_iterations(const int &iter) = 0;
	virtual bool converged(void) const = 0;
protected :
	interface_relaxation_solver() {};
};


#endif
