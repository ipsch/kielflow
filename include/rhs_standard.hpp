#ifndef RHS_STANDARD_HPP_
#define RHS_STANDARD_HPP_


#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)

#include "interface_rhs.hpp"
#include "interface_poisson_solver.hpp"
#include "parameters.hpp"
#include "field.hpp"
#include "operations.hpp"
#include "OP_FFT.hpp"
#include "OP_iFFT.hpp"
#include "OP_partial_derivative.hpp"
#include "OP_dealiasing_23rd.hpp"


#if defined(_MY_VERBOSE)
#include "logger.hpp"
#endif

#define _RHS_POISSON
#define _RHS_DISSIPATION
#define _RHS_ADVECTION    // CHECKED
#define _RHS_DIFFUSION
//#define _RHS_E_INT
//#define _RHS_E_EXT
#define _RHS_CONTINUITY
//#define _RHS_PENALIZATION
//#define _RHS_SVISCOSITY

class rhs_standard : public interface_rhs
{
public :
	rhs_standard(parameters &Params, interface_poisson_solver &Solver,
			field_real &potential, field_real &density_dust,
			field_real &solid_mask, const grid &domain);
	void solve(const double &t, field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni);
	~rhs_standard() { fftw_free(field_SV);	}

protected :
	const int N;
	parameters &my_params;
	interface_poisson_solver &my_poisson_solver;
	field_real &Phi;
	field_real &nd;
	field_real &my_solids;

	field_real Ux;
	field_real Uy;
	field_real Uz;

	field_imag dFUx_dx;
	field_imag dFUy_dx;
	field_imag dFUz_dx;

	field_imag dFUx_dy;
	field_imag dFUy_dy;
	field_imag dFUz_dy;

	field_imag dFUx_dz;
	field_imag dFUy_dz;
	field_imag dFUz_dz;

	field_real Buffer_1st;
	field_real Buffer_2nd;

	field_imag FBuffer_Ux;
	field_imag FBuffer_Uy;
	field_imag FBuffer_Uz;
	field_imag FBuffer_ni;

	field_imag FBuffer_1st;
	field_imag FBuffer_2nd;

	OP_FFT my_FFT;
	OP_iFFT my_iFFT;
	OP_partial_derivative d_dx;
	OP_partial_derivative d_dy;
	OP_partial_derivative d_dz;
	OP_dealiasing my_dealiasing;

	double * field_SV;
	double * penalization_barrier;
	double * penalization_grain;
	double * penalization_ramp;
	int * penalization_mask;

};



#endif /* RHS_STANDARD_HPP_ */



