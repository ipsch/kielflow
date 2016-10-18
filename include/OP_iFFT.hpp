#ifndef OP_IFFT_HPP_
#define OP_IFFT_HPP_

#include "grid.hpp"
#include "field_real.hpp"
#include "field_imag.hpp"
#include "fftw3.h"

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

class OP_iFFT
{
	public :
		OP_iFFT(const grid &domain);
		~OP_iFFT();
		void operator() (const field_imag &in, field_real &out);

	private :
		OP_iFFT() : N(0) {	}
		const int N;
		fftw_plan my_plan;
		fftw_complex * input;
		double * output;
};

#endif
