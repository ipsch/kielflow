#ifndef NEWTON_METHOD_HPP_
#define NEWTON_METHOD_HPP_

#include <limits>
#include <cmath>
#include "CFL.hpp"

// ToDo : built a interface
// ToDo : move CFL to another base

class newton_method
{
public :
	newton_method(double (&function)(double), double (&derivative)(double), const double &accuracy = 0.e-5);
	double solve(const double &x_i) const;
private :
	double (&f)(double);
	double (&df)(double);
	const double eps;
};




#endif
