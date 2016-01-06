#ifndef CFL_HPP_
#define CFL_HPP_

/*
 *
 * function to check if the cauchy-friedrich-levy criterion
 * is fullfilled
 *
 */

#include <cmath> // to check if new timestep isn't equal inf

#include "field.hpp"
#include "operations.hpp"
#include "o_math.hpp"

// ToDo : richtig loggen
#include <iostream>

// ToDo : implement if needed
//double cfl_single(const double &dt, const field_real  &XX);
//double cfl_single(const double &dt, const field_imag &FXX);
double cfl_all(const double &dt, const field_real  &Ux, const field_real  &Uy, const field_real  &Uz);
double cfl_all(const double &dt, const field_imag &FUx, const field_imag &FUy, const field_imag &FUz);


#endif /* END CFL_HPP_ */
