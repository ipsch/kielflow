#ifndef SOLVER_POISSON_JACOBI_NLIN_HPP_
#define SOLVER_POISSON_JACOBI_NLIN_HPP_

/*
 * nonlinear jacobi-method to solve
 * the nonlinear poission-equation
 * div(grad(Phi)) + rho_i + exp(-Phi) = 0 in Omega
 * Phi = Phi_delta on Omega_delta
 *
 */

#include "interface_poisson_solver.hpp"
#include "masks.hpp"
#include "newton_method.hpp"


// ToDo : These won't be necessary in final version
#include "o_math.hpp"
#include "o_string.hpp"
#include <fstream>
#include <thread>
#include "IO.hpp"
#include <omp.h>


#include <iostream>
#include <iomanip>
#include <string>
#include <functional>

#if defined(_MY_VERBOSE)  || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TIDEOUS)
#include "logger.hpp"
#endif



class solver_poisson_jacobi_nlin : public interface_relaxation_solver
// relaxation solver for the nonlinear poission-equation
// on a nonuniform grid
// with direchlet boundaries (Phi=0 at x=x0, y=y0 and z=z0)
// with additional boundaries in form of a (metallic) submerges body
{


public :
	solver_poisson_jacobi_nlin(interface_3d_fkt &boundary, interface_3d_fkt &val_boundary, const double &w = 1.);
	void solve(field_real &Phi_IO, field_real &rho);
	bool converged(void) const {return converged_;}
	void set_max_iterations(const int &iter) {max_iterations = iter;}
	void set_tolerance(const double &tolerance){limit_max = tolerance;};

	double limit_max;
	double limit_sum;

private :
	void main_loop(void);
	void iteration_loop(const field_real &in, field_real &out, const field_real &rho);
    void save_evolution(const field_real &Phi, const field_real &rho) const;



	double norm_max;
	double norm_sum;
	double supremum;
	double infinum;

	interface_3d_fkt &H;
	interface_3d_fkt &val_H;

	double * HX;
	double * HY;
	double * HZ;
	void H_create(const grid &Omega);
	void H_delete();

	std::string my_logfile;

	double omega_SOR; // Successive Over-Relaxation parameter
    double omega_NEWTON;

    double norm_sum_old;
    double norm_max_old;

	int iteration;
	int invocation;
	int max_iterations;
	static int iterations_total;
	bool converged_;
	bool use_boundary_;

	// ToDo : move newton-method to own header (make it stand-alone)
	double newton(const int &i, const int j, const int k,\
			const field_real &Phi, const double &rho_ijk) const;
	double f_df(const double &x, \
			const int &i, const int j, const int k,\
			const field_real &Phi, const double &rho_ijk) const;
	const double eps;

};






#endif
