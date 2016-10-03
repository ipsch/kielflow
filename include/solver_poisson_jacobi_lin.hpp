#ifndef SOLVER_POISSON_JACOBI_LIN_HPP_
#define SOLVER_POISSON_JACOBI_LIN_HPP_

/*
 * nonlinear jacobi-method to solve
 * the nonlinear poission-equation
 * div(grad(Phi)) + rho_i + exp(-Phi) = 0 in Omega
 * Phi = Phi_delta on Omega_delta
 *
 */

#include "interface_poisson_solver.hpp"
#include "masks.hpp"

// ToDo : These won't be necessary in final version
//ToDo : write own diagnostics class
#include "o_math.hpp"
#include <sstream>
#include <iostream>

#ifdef _MY_VERBOSE
#include <fstream> // ToDo : Quick and dirty logging to file
#include "logger.hpp"
#endif

#if defined(__FRONTEND__)
extern double global_theta;
#endif

class solver_poisson_jacobi_lin : public interface_relaxation_solver
// relaxation solver for the nonlinear poission-equation
// on a nonuniform grid
// with direchlet boundaries (Phi=0 at x=x0, y=y0 and z=z0)
// with additional boundaries in form of a (metallic) submerged body
{


public :
	solver_poisson_jacobi_lin(interface_3d_fkt &boundary, interface_3d_fkt &val_boundary);
	~solver_poisson_jacobi_lin() {};
	void solve(field_real &Phi_IO, field_real &rho);
	void set_max_iterations(const int &iter) {max_iterations = iter;}
private :
	void main_loop(void);
	void iteration_loop(const field_real &in, field_real &out, const field_real &rho_i);

	double get_PG(const field_real &in, const int &i, const int &j, const int &k) const;
	double get_HXX(const axis * const A, const int &i, double &hp, double &hm) const;

	void check_convergence(const field_real &field_new, const field_real &field_old);

	double norm_maximum;
	double norm_sum_total;
	double norm_sum_i;
	double norm_sum_j;
	double norm_sum_k;

	interface_3d_fkt &H;
	interface_3d_fkt &val_H;

	int iteration;
	int max_iterations;
	int invocations;
	double hxp;
	double hxm;
	double hyp;
	double hym;
	double hzp;
	double hzm;
};



#endif /* SOLVER_POISSON_JACOBI_LIN_HPP_ */

