#include "solver_poisson_base.hpp"


//void (*solver_poisson::solve)(const field_imag &, field_imag &) = &solver_poisson::method_debye_linear;

/*
void solver_poisson::set_method(int opt)
{
   #ifdef _MY_VERBOSE
	logger log("solver_poisson::set_method");
	log << "start";
   #endif
	switch (opt)
	{
	case 1:
       #ifdef _MY_VERBOSE
		log << "setting solver to linear debye-hückel-approximation (opt=1)";
       #endif
		solver_poisson::solve = &solver_poisson::method_debye_linear;
		break;

	default :
       #ifdef _MY_VERBOSE
		log << "WARNING! No solver choosen";
		log << "setting solver to linear debye-hückel-approximation (default)";
       #endif
		solver_poisson::solve = &solver_poisson::method_debye_linear;
	}

   #ifdef _MY_VERBOSE
	log << "done";
   #endif
	return;
}
*/


solver_poisson::solver_poisson()
{
   #ifdef _MY_VERBOSE
	logger log("solver_poisson");
	log << "solver_poisson(const grid_Fo &domain)";

	invocations = 0;
   #endif
}


solver_poisson::~solver_poisson()
{
   #ifdef _MY_VERBOSE
	logger log("solver_poisson");
	log << "~solver_poisson()";
   #endif
}
