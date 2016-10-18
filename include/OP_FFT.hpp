#ifndef FFT_HPP_
#define FFT_HPP_

#include "grid.hpp"
#include "field_real.hpp"
#include "field_imag.hpp"
#include "fftw3.h"

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

class OP_FFT
{
	public :
		OP_FFT(const grid &domain);
		~OP_FFT();
		void operator() (const field_real &in, field_imag &out);

	private :
		OP_FFT() : N(0) {	};
		const int N;
		fftw_plan my_plan;
		double * input;
		fftw_complex * output;
};

#endif

