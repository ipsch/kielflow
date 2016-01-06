#ifndef SOLVER_POISSON_BASE_HPP_
#define SOLVER_POISSON_BASE_HPP_




// gcc standard C++ libraries (common stuff)
#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <vector>

// uncommon libraries
#include "fftw3.h"

// own includes

#include "field.hpp"
#include "operations.hpp"
#include "IO.hpp"


#include "counter.hpp"
#include "particle.hpp"
#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif
// magic stuff

#include <sstream>
#include <omp.h>
#include <fstream>


class solver_poisson
{
public :
	solver_poisson();
	virtual ~solver_poisson();
	virtual void solve(field_real &X, const field_real &b)= 0;


private :

protected :
	int invocations;
};





#endif
