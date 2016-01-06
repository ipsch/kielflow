#ifndef SOLVER_POISSON_MULTIGRID_HPP_
#define SOLVER_POISSON_MULTIGRID_HPP_

/*
 *
 * 	Multigrid solver scheme
 * 	utilizuing a relaxations solver
 *
 */

#include <string>
#include <sstream>
#include "interface_poisson_solver.hpp"
//#include "solver_poisson_jacobi_lin.hpp"
//#include "solver_poisson_jacobi_nlin.hpp"

#include "field.hpp"
#include "field_interpolation.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif


typedef enum
{
    lvl_down, lvl_keep, lvl_up
} MG_lvl_control;

class solver_poisson_multigrid : public interface_poisson_solver
{
public :
	// inherited from interface
	solver_poisson_multigrid(interface_relaxation_solver &method);
	~solver_poisson_multigrid();
	void solve(field_real &Phi_IO, field_real &rho);
	// my own methods
	void set_method(interface_relaxation_solver &method) {ptr_relaxation_method = &method;}
	void set_level_control(const int &N, int * lvl_steps, MG_lvl_control * lvl_C);
private :
	solver_poisson_multigrid();

	interface_relaxation_solver * ptr_relaxation_method; // Method used for soving
	int my_cascades;                                     // Number of cascades (big steps) in the MG-cycle
	int * my_steps;                                      // Array with number of sub-steps per cascade
	void (*I_hto2h)(const field_real&, field_real&);     // Interpolation Operator
	void (*I_2htoh)(const field_real&, field_real&);     // Averageing Operator
	MG_lvl_control * my_lvl;                             //


};


#endif
