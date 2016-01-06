#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTIVS_HPP


//#include <complex.h> // simple arithmetics with complex numbers (defines symbol "I" for imaginary unit)
#include <stdio.h>      /* printf, scanf, puts, NULL */

#include <cmath>
#include "fftw3.h"
#include <iostream>

#include "grid.hpp"
#include "field.hpp"

#include "o_math.hpp"

class diagnostics
{

public :
	diagnostics();
	void analyze(const field_real &that);
	double get_supremum(void) {return supremum;}
	double get_infinum(void) {return infinum;}
	int get_i_supremum(void) {return i_supremum;}
	int get_j_supremum(void) {return j_supremum;}
	int get_k_supremum(void) {return k_supremum;}
	int get_i_infinum(void) {return i_infinum;}
	int get_j_infinum(void) {return j_infinum;}
	int get_k_infinum(void) {return k_infinum;}
	//void set_ijk(int index, int &i, int &j, int &k);
	double supremum;
	int i_supremum;
	int j_supremum;
	int k_supremum;
	double infinum;
	int i_infinum;
	int j_infinum;
	int k_infinum;
};


# endif // ende der Diagnose-tool Klasse
